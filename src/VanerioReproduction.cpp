#include <SprayThicknessPrediction/VanerioReproduction.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace spraythickness::published
{
    namespace
    {
        constexpr double kEpsilon = 1.0e-12;

        struct GridCell
        {
            std::int64_t x{ 0 };
            std::int64_t y{ 0 };

            bool operator==(const GridCell& other) const
            {
                return x == other.x && y == other.y;
            }
        };

        struct GridCellHash
        {
            std::size_t operator()(const GridCell& cell) const
            {
                return static_cast<std::size_t>(cell.x) * 73856093u
                    ^ static_cast<std::size_t>(cell.y) * 19349663u;
            }
        };

        struct ClosestVertex
        {
            std::size_t index{ 0 };
            double distanceSquared{ std::numeric_limits<double>::infinity() };
        };

        double elapsedMilliseconds(
            const std::chrono::steady_clock::time_point& started)
        {
            return std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - started).count();
        }

        void localBasis(
            const Eigen::Vector3d& axis,
            Eigen::Vector3d& x,
            Eigen::Vector3d& y)
        {
            const Eigen::Vector3d reference = std::abs(axis.z()) < 0.9
                ? Eigen::Vector3d::UnitZ()
                : Eigen::Vector3d::UnitY();
            x = reference.cross(axis).normalized();
            y = axis.cross(x).normalized();
        }

        double profileIntegral(double radius, double k2, double stretch)
        {
            const double limit = std::min(1.0, 1.0 / stretch);
            const double scaledLimit = stretch * limit;
            const double gaussianIntegral = -std::expm1(
                -k2 * scaledLimit * scaledLimit)
                / (2.0 * k2 * stretch * stretch);
            const double baselineIntegral =
                0.5 * std::exp(-k2) * limit * limit;
            return radius * radius
                * (gaussianIntegral - baselineIntegral)
                / (-std::expm1(-k2));
        }

        double particleDistribution(
            double radialDistance,
            double stretch,
            double radius,
            double k2)
        {
            const double transformedRadius = stretch * radialDistance;
            if(radialDistance >= radius || transformedRadius >= radius) {
                return 0.0;
            }
            const double base = (std::exp(-k2
                    * (transformedRadius / radius)
                    * (transformedRadius / radius))
                - std::exp(-k2)) / (1.0 - std::exp(-k2));
            const double initialIntegral = profileIntegral(radius, k2, 1.0);
            const double stretchedIntegral = profileIntegral(radius, k2, stretch);
            const double scale = initialIntegral / stretchedIntegral;
            return scale * base;
        }

        std::vector<std::vector<std::size_t>> adjacentFaces(
            const TriangleMesh& mesh)
        {
            std::vector<std::vector<std::size_t>> adjacency(mesh.vertices.size());
            for(std::size_t faceIndex = 0; faceIndex < mesh.faces.size(); ++faceIndex) {
                for(const std::uint32_t vertex : mesh.faces[faceIndex]) {
                    adjacency[vertex].push_back(faceIndex);
                }
            }
            return adjacency;
        }

        std::vector<bool> visibleFaces(
            const TriangleMesh& mesh,
            const SprayPose& pose,
            double gridStep,
            double radius)
        {
            const Eigen::Vector3d axis = pose.axis.normalized();
            Eigen::Vector3d basisX;
            Eigen::Vector3d basisY;
            localBasis(axis, basisX, basisY);
            std::unordered_map<GridCell, ClosestVertex, GridCellHash> closest;
            const double rayRadiusSquared = 0.5 * gridStep * gridStep;
            for(std::size_t vertexIndex = 0;
                vertexIndex < mesh.vertices.size(); ++vertexIndex) {
                const Eigen::Vector3d relative =
                    mesh.vertices[vertexIndex] - pose.position;
                const double axial = relative.dot(axis);
                if(axial <= 0.0) {
                    continue;
                }
                const double x = relative.dot(basisX);
                const double y = relative.dot(basisY);
                if(x * x + y * y > radius * radius) {
                    continue;
                }
                const double distanceSquared = relative.squaredNorm();
                const std::int64_t firstX =
                    static_cast<std::int64_t>(std::floor(x / gridStep));
                const std::int64_t firstY =
                    static_cast<std::int64_t>(std::floor(y / gridStep));
                // Grid points are nozzle-parallel ray axes; each covers a
                // cylindrical sub-volume with a circular cross-section.
                for(std::int64_t rayX = firstX; rayX <= firstX + 1; ++rayX) {
                    for(std::int64_t rayY = firstY; rayY <= firstY + 1; ++rayY) {
                        const double dx = x - static_cast<double>(rayX) * gridStep;
                        const double dy = y - static_cast<double>(rayY) * gridStep;
                        if(dx * dx + dy * dy > rayRadiusSquared) {
                            continue;
                        }
                        const GridCell cell{ rayX, rayY };
                        const auto found = closest.find(cell);
                        if(found == closest.end()
                            || distanceSquared < found->second.distanceSquared) {
                            closest[cell] = { vertexIndex, distanceSquared };
                        }
                    }
                }
            }

            const auto adjacency = adjacentFaces(mesh);
            std::vector<bool> visible(mesh.faces.size(), false);
            for(const auto& item : closest) {
                for(const std::size_t face : adjacency[item.second.index]) {
                    visible[face] = true;
                }
            }
            return visible;
        }

        double maximumEdgeLength(const TriangleMesh& mesh)
        {
            double maximum = 0.0;
            for(const auto& face : mesh.faces) {
                for(int edge = 0; edge < 3; ++edge) {
                    maximum = std::max(maximum,
                        (mesh.vertices[face[(edge + 1) % 3]]
                            - mesh.vertices[face[edge]]).norm());
                }
            }
            return maximum;
        }

        bool beginsNewPass(const std::vector<SprayPose>& trajectory,
            std::size_t poseIndex)
        {
            if(poseIndex == 0 || !trajectory[poseIndex].sprayEnabled) {
                return false;
            }
            if(!trajectory[poseIndex - 1].sprayEnabled) {
                return true;
            }
            if(poseIndex + 1 >= trajectory.size()) {
                return false;
            }
            const Eigen::Vector3d before = trajectory[poseIndex].position
                - trajectory[poseIndex - 1].position;
            const Eigen::Vector3d after = trajectory[poseIndex + 1].position
                - trajectory[poseIndex].position;
            return before.squaredNorm() > kEpsilon
                && after.squaredNorm() > kEpsilon
                && before.dot(after) < 0.0;
        }

        TriangleMesh bisectLongEdges(
            const TriangleMesh& input,
            double maximumEdge,
            std::vector<double>& vertexThickness)
        {
            TriangleMesh current = input;
            while(maximumEdgeLength(current) > maximumEdge) {
                TriangleMesh next;
                next.vertices = current.vertices;
                std::vector<double> nextThickness = vertexThickness;
                next.faces.reserve(current.faces.size() * 4);
                std::unordered_map<std::uint64_t, std::uint32_t> midpoints;
                const auto midpoint = [&](std::uint32_t first,
                                          std::uint32_t second) {
                    const std::uint64_t key =
                        (static_cast<std::uint64_t>(std::min(first, second)) << 32)
                        | std::max(first, second);
                    const auto found = midpoints.find(key);
                    if(found != midpoints.end()) {
                        return found->second;
                    }
                    const std::uint32_t index =
                        static_cast<std::uint32_t>(next.vertices.size());
                    next.vertices.push_back(0.5
                        * (current.vertices[first] + current.vertices[second]));
                    nextThickness.push_back(0.5
                        * (vertexThickness[first] + vertexThickness[second]));
                    midpoints.emplace(key, index);
                    return index;
                };
                for(const auto& face : current.faces) {
                    std::array<double, 3> lengths{
                        (current.vertices[face[1]] - current.vertices[face[0]]).norm(),
                        (current.vertices[face[2]] - current.vertices[face[1]]).norm(),
                        (current.vertices[face[0]] - current.vertices[face[2]]).norm()
                    };
                    const std::array<bool, 3> split{
                        lengths[0] > maximumEdge,
                        lengths[1] > maximumEdge,
                        lengths[2] > maximumEdge
                    };
                    const int count = static_cast<int>(split[0])
                        + static_cast<int>(split[1])
                        + static_cast<int>(split[2]);
                    if(count == 0) {
                        next.faces.push_back(face);
                        continue;
                    }
                    const std::uint32_t a = face[0];
                    const std::uint32_t b = face[1];
                    const std::uint32_t c = face[2];
                    const std::uint32_t ab = split[0] ? midpoint(a, b) : 0;
                    const std::uint32_t bc = split[1] ? midpoint(b, c) : 0;
                    const std::uint32_t ca = split[2] ? midpoint(c, a) : 0;
                    if(count == 1) {
                        if(split[0]) {
                            next.faces.push_back({ a, ab, c });
                            next.faces.push_back({ ab, b, c });
                        } else if(split[1]) {
                            next.faces.push_back({ b, bc, a });
                            next.faces.push_back({ bc, c, a });
                        } else {
                            next.faces.push_back({ c, ca, b });
                            next.faces.push_back({ ca, a, b });
                        }
                    } else if(count == 2) {
                        if(!split[0]) {
                            next.faces.push_back({ c, ca, bc });
                            next.faces.push_back({ ca, a, b });
                            next.faces.push_back({ ca, b, bc });
                        } else if(!split[1]) {
                            next.faces.push_back({ a, ab, ca });
                            next.faces.push_back({ ab, b, c });
                            next.faces.push_back({ ab, c, ca });
                        } else {
                            next.faces.push_back({ b, bc, ab });
                            next.faces.push_back({ bc, c, a });
                            next.faces.push_back({ bc, a, ab });
                        }
                    } else {
                        next.faces.push_back({ a, ab, ca });
                        next.faces.push_back({ ab, b, bc });
                        next.faces.push_back({ ca, bc, c });
                        next.faces.push_back({ ab, bc, ca });
                    }
                }
                current = std::move(next);
                vertexThickness = std::move(nextThickness);
            }
            return current;
        }
    }

    void VanerioParameters::validate() const
    {
        depositionEfficiencyByTangentAngle.validate(
            "Vanerio angle-dependent deposition efficiency");
        depositionEfficiencyByDistance.validate(
            "Vanerio distance-dependent deposition efficiency");
        profileStretchByDistance.validate(
            "Vanerio distance-dependent profile stretch");
        if(!std::isfinite(growthRateCoefficientMetersPerSecond)
            || growthRateCoefficientMetersPerSecond <= 0.0
            || !std::isfinite(jetRadiusMeters) || jetRadiusMeters <= 0.0
            || !std::isfinite(jetShapeCoefficientK2)
            || jetShapeCoefficientK2 <= 0.0
            || !std::isfinite(maximumMeshEdgeMeters)
            || maximumMeshEdgeMeters <= 0.0
            || !std::isfinite(shadowGridStepMeters)
            || shadowGridStepMeters <= 0.0) {
            throw std::invalid_argument(
                "Vanerio reproduction requires A, rn, k2, mesh and shadow "
                "resolutions, and all experimental calibration curves.");
        }
    }

    VanerioResult VanerioReproducer::run(
        const VanerioInputModel& input,
        const VanerioParameters& parameters,
        const ReproductionExecution& execution)
    {
        reportDiagnostic(execution, "Input validation",
            "initial faces="
                + std::to_string(input.initialStlSurface.faces.size())
                + ", nozzle poses="
                + std::to_string(input.nozzleTrajectory.size()));
        parameters.validate();
        input.initialStlSurface.validate();
        if(input.nozzleTrajectory.size() < 2) {
            throw std::invalid_argument(
                "Vanerio reproduction requires at least two timed nozzle poses.");
        }

        const auto started = std::chrono::steady_clock::now();
        VanerioResult result;
        result.vertexThicknessMeters.assign(
            input.initialStlSurface.vertices.size(), 0.0);
        reportDiagnostic(execution, "Mesh remeshing",
            "maximum edge m="
                + std::to_string(parameters.maximumMeshEdgeMeters));
        result.evolvedStlSurface = bisectLongEdges(
            input.initialStlSurface, parameters.maximumMeshEdgeMeters,
            result.vertexThicknessMeters);
        reportDiagnostic(execution, "Visibility, deposition and surface growth",
            "remeshed faces="
                + std::to_string(result.evolvedStlSurface.faces.size()));
        result.implementationNotes.push_back(
            "The paper specifies a maximum-edge remeshing constraint but does "
            "not publish the remesher; shared-edge conforming bisection is used.");
        result.implementationNotes.push_back(
            "This calibration uses the paper's Gaussian jet-profile branch; "
            "non-Gaussian material profiles require measured input not provided here.");
        result.implementationNotes.push_back(
            "Nozzle-parallel ray cylinders use grid step / sqrt(2) as their "
            "radius to cover the outlet grid; the paper does not specify it.");
        result.implementationNotes.push_back(
            "Pass boundaries are inferred from a spraying direction reversal "
            "or restart after a non-spraying pose; the trajectory has no pass IDs.");
        result.statistics.trajectorySampleCount = input.nozzleTrajectory.size();
        std::size_t visibleFaceCount = 0;
        std::size_t invalidGeometryCount = 0;
        std::size_t outsideProfileCount = 0;
        double minimumScaledRadius = std::numeric_limits<double>::infinity();

        for(std::size_t poseIndex = 0;
            poseIndex + 1 < input.nozzleTrajectory.size(); ++poseIndex) {
            if(canceled(execution)) {
                result.canceled = true;
                break;
            }
            const SprayPose& pose = input.nozzleTrajectory[poseIndex];
            const double dt = input.nozzleTrajectory[poseIndex + 1].timeSeconds
                - pose.timeSeconds;
            if(dt <= 0.0) {
                throw std::invalid_argument(
                    "Vanerio trajectory times must be strictly increasing.");
            }
            if(beginsNewPass(input.nozzleTrajectory, poseIndex)
                && maximumEdgeLength(result.evolvedStlSurface)
                    > parameters.maximumMeshEdgeMeters) {
                result.evolvedStlSurface = bisectLongEdges(
                    result.evolvedStlSurface,
                    parameters.maximumMeshEdgeMeters,
                    result.vertexThicknessMeters);
                reportDiagnostic(execution, "Inter-pass remeshing",
                    "before pass at pose " + std::to_string(poseIndex + 1)
                        + ", faces="
                        + std::to_string(result.evolvedStlSurface.faces.size()));
            }
            if(!pose.sprayEnabled) {
                continue;
            }
            if(pose.axis.squaredNorm() <= kEpsilon) {
                throw std::invalid_argument("Vanerio nozzle axis is zero.");
            }
            const Eigen::Vector3d particleAxis = pose.axis.normalized();
            const Eigen::Vector3d growthDirection = -particleAxis;
            const std::vector<bool> visible = visibleFaces(
                result.evolvedStlSurface,
                pose,
                parameters.shadowGridStepMeters,
                parameters.jetRadiusMeters);
            std::vector<double> faceIncrements(
                result.evolvedStlSurface.faces.size(), 0.0);
            for(std::size_t faceIndex = 0;
                faceIndex < result.evolvedStlSurface.faces.size(); ++faceIndex) {
                if(!visible[faceIndex]) {
                    ++result.statistics.hiddenElementCount;
                    continue;
                }
                ++visibleFaceCount;
                const Eigen::Vector3d centroid =
                    faceCentroid(result.evolvedStlSurface, faceIndex);
                const Eigen::Vector3d normal =
                    faceNormal(result.evolvedStlSurface, faceIndex);
                const Eigen::Vector3d relative = centroid - pose.position;
                const double distance = relative.norm();
                const double axialDistance = relative.dot(particleAxis);
                if(distance <= kEpsilon || axialDistance <= 0.0) {
                    ++invalidGeometryCount;
                    continue;
                }
                const double normalDotAxis =
                    std::abs(normal.dot(particleAxis));
                if(normalDotAxis <= kEpsilon) {
                    ++invalidGeometryCount;
                    continue;
                }
                const double tangentAngle =
                    normal.cross(particleAxis).norm() / normalDotAxis;
                const double radialDistance =
                    (relative - axialDistance * particleAxis).norm();
                const double stretch = parameters.profileStretchByDistance.interpolate(
                    distance,
                    "Vanerio distance-dependent profile stretch");
                minimumScaledRadius = std::min(minimumScaledRadius,
                    stretch * radialDistance / parameters.jetRadiusMeters);
                const double distribution = particleDistribution(
                    radialDistance,
                    stretch,
                    parameters.jetRadiusMeters,
                    parameters.jetShapeCoefficientK2);
                if(distribution <= 0.0) {
                    ++outsideProfileCount;
                    continue;
                }
                const double angleEfficiency =
                    parameters.depositionEfficiencyByTangentAngle.interpolate(
                        tangentAngle,
                        "Vanerio angle-dependent deposition efficiency");
                const double distanceEfficiency =
                    parameters.depositionEfficiencyByDistance.interpolate(
                        distance,
                        "Vanerio distance-dependent deposition efficiency");
                faceIncrements[faceIndex] =
                    parameters.growthRateCoefficientMetersPerSecond
                    * angleEfficiency * distanceEfficiency * distribution * dt;
                ++result.statistics.candidatePairCount;
                ++result.statistics.visibilityQueryCount;
            }

            const auto adjacency = adjacentFaces(result.evolvedStlSurface);
            for(std::size_t vertexIndex = 0;
                vertexIndex < result.evolvedStlSurface.vertices.size(); ++vertexIndex) {
                if(adjacency[vertexIndex].empty()) {
                    continue;
                }
                double increment = 0.0;
                for(const std::size_t face : adjacency[vertexIndex]) {
                    increment += faceIncrements[face];
                }
                increment /= static_cast<double>(adjacency[vertexIndex].size());
                result.vertexThicknessMeters[vertexIndex] += increment;
                result.evolvedStlSurface.vertices[vertexIndex] +=
                    increment * growthDirection;
            }
            reportProgress(execution,
                static_cast<double>(poseIndex + 1)
                    / static_cast<double>(input.nozzleTrajectory.size() - 1),
                "Vanerio 2021 reproduction");
        }

        reportDiagnostic(execution, "Result assembly",
            "visible faces=" + std::to_string(visibleFaceCount)
                + ", invalid geometry="
                + std::to_string(invalidGeometryCount)
                + ", outside profile="
                + std::to_string(outsideProfileCount)
                + ", minimum scaled radius="
                + (std::isfinite(minimumScaledRadius)
                    ? std::to_string(minimumScaledRadius) : "N/A")
                + ", accepted contributions="
                + std::to_string(result.statistics.candidatePairCount));
        result.statistics.evaluatedElementCount =
            result.evolvedStlSurface.faces.size();
        result.statistics.elapsedMilliseconds = elapsedMilliseconds(started);
        return result;
    }
}
