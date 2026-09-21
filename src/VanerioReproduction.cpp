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
            double axialDistance{ std::numeric_limits<double>::infinity() };
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

        double profileIntegral(double radius, double k2)
        {
            constexpr int intervals = 512;
            const double step = radius / static_cast<double>(intervals);
            double sum = 0.0;
            for(int index = 0; index <= intervals; ++index) {
                const double r = static_cast<double>(index) * step;
                const double profile = r < radius
                    ? (std::exp(-k2 * (r / radius) * (r / radius))
                        - std::exp(-k2)) / (1.0 - std::exp(-k2))
                    : 0.0;
                const double weight = index == 0 || index == intervals
                    ? 0.5 : 1.0;
                sum += weight * r * profile;
            }
            return sum * step;
        }

        double particleDistribution(
            double radialDistance,
            double stretch,
            double radius,
            double k2)
        {
            const double transformedRadius = stretch * radialDistance;
            if(transformedRadius >= radius) {
                return 0.0;
            }
            const double base = (std::exp(-k2
                    * (transformedRadius / radius)
                    * (transformedRadius / radius))
                - std::exp(-k2)) / (1.0 - std::exp(-k2));
            const double initialIntegral = profileIntegral(radius, k2);
            const double stretchedIntegral = initialIntegral / (stretch * stretch);
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
                const GridCell cell{
                    static_cast<std::int64_t>(std::floor(x / gridStep)),
                    static_cast<std::int64_t>(std::floor(y / gridStep))
                };
                auto found = closest.find(cell);
                if(found == closest.end() || axial < found->second.axialDistance) {
                    closest[cell] = { vertexIndex, axial };
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

        TriangleMesh bisectLongEdges(
            const TriangleMesh& input,
            double maximumEdge)
        {
            TriangleMesh current = input;
            while(maximumEdgeLength(current) > maximumEdge) {
                TriangleMesh next;
                next.vertices = current.vertices;
                next.faces.reserve(current.faces.size() * 2);
                for(const auto& face : current.faces) {
                    std::array<double, 3> lengths{
                        (current.vertices[face[1]] - current.vertices[face[0]]).norm(),
                        (current.vertices[face[2]] - current.vertices[face[1]]).norm(),
                        (current.vertices[face[0]] - current.vertices[face[2]]).norm()
                    };
                    const int longest = static_cast<int>(std::distance(
                        lengths.begin(), std::max_element(lengths.begin(), lengths.end())));
                    if(lengths[longest] <= maximumEdge) {
                        next.faces.push_back(face);
                        continue;
                    }
                    const std::uint32_t first = face[longest];
                    const std::uint32_t second = face[(longest + 1) % 3];
                    const std::uint32_t opposite = face[(longest + 2) % 3];
                    const std::uint32_t midpoint =
                        static_cast<std::uint32_t>(next.vertices.size());
                    next.vertices.push_back(0.5
                        * (current.vertices[first] + current.vertices[second]));
                    next.faces.push_back({ first, midpoint, opposite });
                    next.faces.push_back({ midpoint, second, opposite });
                }
                current = std::move(next);
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
        parameters.validate();
        input.initialStlSurface.validate();
        if(input.nozzleTrajectory.size() < 2) {
            throw std::invalid_argument(
                "Vanerio reproduction requires at least two timed nozzle poses.");
        }

        const auto started = std::chrono::steady_clock::now();
        VanerioResult result;
        result.evolvedStlSurface = bisectLongEdges(
            input.initialStlSurface, parameters.maximumMeshEdgeMeters);
        result.implementationNotes.push_back(
            "The paper specifies a maximum-edge remeshing constraint but does "
            "not publish the remesher; longest-edge bisection is used.");
        result.statistics.trajectorySampleCount = input.nozzleTrajectory.size();

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
                const Eigen::Vector3d centroid =
                    faceCentroid(result.evolvedStlSurface, faceIndex);
                const Eigen::Vector3d normal =
                    faceNormal(result.evolvedStlSurface, faceIndex);
                const Eigen::Vector3d relative = centroid - pose.position;
                const double distance = relative.norm();
                const double axialDistance = relative.dot(particleAxis);
                if(distance <= kEpsilon || axialDistance <= 0.0) {
                    continue;
                }
                const double normalDotAxis =
                    std::abs(normal.dot(particleAxis));
                if(normalDotAxis <= kEpsilon) {
                    continue;
                }
                const double tangentAngle =
                    normal.cross(particleAxis).norm() / normalDotAxis;
                const double radialDistance =
                    (relative - axialDistance * particleAxis).norm();
                const double stretch = parameters.profileStretchByDistance.interpolate(
                    distance,
                    "Vanerio distance-dependent profile stretch");
                const double distribution = particleDistribution(
                    radialDistance,
                    stretch,
                    parameters.jetRadiusMeters,
                    parameters.jetShapeCoefficientK2);
                if(distribution <= 0.0) {
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
                result.evolvedStlSurface.vertices[vertexIndex] +=
                    increment * growthDirection;
            }
            if(maximumEdgeLength(result.evolvedStlSurface)
                > parameters.maximumMeshEdgeMeters) {
                result.evolvedStlSurface = bisectLongEdges(
                    result.evolvedStlSurface,
                    parameters.maximumMeshEdgeMeters);
            }
            reportProgress(execution,
                static_cast<double>(poseIndex + 1)
                    / static_cast<double>(input.nozzleTrajectory.size() - 1),
                "Vanerio 2021 reproduction");
        }

        result.statistics.evaluatedElementCount =
            result.evolvedStlSurface.faces.size();
        result.statistics.elapsedMilliseconds = elapsedMilliseconds(started);
        return result;
    }
}
