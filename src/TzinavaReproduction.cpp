#include <SprayThicknessPrediction/TzinavaReproduction.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace spraythickness::published
{
    namespace
    {
        constexpr double kPi = 3.14159265358979323846;
        constexpr double kEpsilon = 1.0e-12;

        struct EdgeKey
        {
            std::uint32_t first{ 0 };
            std::uint32_t second{ 0 };

            bool operator==(const EdgeKey& other) const
            {
                return first == other.first && second == other.second;
            }
        };

        struct EdgeKeyHash
        {
            std::size_t operator()(const EdgeKey& key) const
            {
                return static_cast<std::size_t>(key.first) * 73856093u
                    ^ static_cast<std::size_t>(key.second) * 19349663u;
            }
        };

        struct GunState
        {
            Eigen::Vector3d position = Eigen::Vector3d::Zero();
            Eigen::Vector3d axis = -Eigen::Vector3d::UnitZ();
            Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
            bool sprayEnabled{ true };
        };

        struct Candidate
        {
            std::size_t faceIndex{ 0 };
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            Eigen::Vector3d normal = Eigen::Vector3d::Zero();
            std::array<Eigen::Vector2d, 3> projected;
            double standOffDistance{ 0.0 };
            double impactAngleDegrees{ 0.0 };
            double facingCosine{ 0.0 };
            double radialDistance{ 0.0 };
            double spotSpeedMetersPerSecond{ 0.0 };
            double projectedNormalDotAxis{ 1.0 };
            bool hidden{ false };
        };

        double elapsedMilliseconds(
            const std::chrono::steady_clock::time_point& started)
        {
            return std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - started).count();
        }

        std::size_t lowerInterval(
            const std::vector<double>& values,
            double value)
        {
            const auto upper = std::upper_bound(values.begin(), values.end(), value);
            if(upper == values.begin()) {
                return 0;
            }
            if(upper == values.end()) {
                return values.size() - 2;
            }
            return static_cast<std::size_t>(upper - values.begin() - 1);
        }

        std::uint32_t midpoint(
            TriangleMesh& output,
            std::uint32_t first,
            std::uint32_t second,
            std::unordered_map<EdgeKey, std::uint32_t, EdgeKeyHash>& cache)
        {
            const EdgeKey key{ std::min(first, second), std::max(first, second) };
            const auto found = cache.find(key);
            if(found != cache.end()) {
                return found->second;
            }
            const std::uint32_t index =
                static_cast<std::uint32_t>(output.vertices.size());
            output.vertices.push_back(
                0.5 * (output.vertices[first] + output.vertices[second]));
            cache.emplace(key, index);
            return index;
        }

        TriangleMesh subdivide(const TriangleMesh& input, double beamRadius)
        {
            TriangleMesh output;
            output.vertices = input.vertices;
            std::unordered_map<EdgeKey, std::uint32_t, EdgeKeyHash> midpointCache;
            const double threshold = 0.5 * beamRadius;
            for(const auto& face : input.faces) {
                const std::array<bool, 3> split{
                    (input.vertices[face[1]] - input.vertices[face[0]]).norm()
                        > threshold,
                    (input.vertices[face[2]] - input.vertices[face[1]]).norm()
                        > threshold,
                    (input.vertices[face[0]] - input.vertices[face[2]]).norm()
                        > threshold
                };
                const int count = static_cast<int>(split[0])
                    + static_cast<int>(split[1]) + static_cast<int>(split[2]);
                if(count == 0) {
                    output.faces.push_back(face);
                    continue;
                }
                const std::uint32_t a = face[0];
                const std::uint32_t b = face[1];
                const std::uint32_t c = face[2];
                const std::uint32_t ab = split[0]
                    ? midpoint(output, a, b, midpointCache) : 0;
                const std::uint32_t bc = split[1]
                    ? midpoint(output, b, c, midpointCache) : 0;
                const std::uint32_t ca = split[2]
                    ? midpoint(output, c, a, midpointCache) : 0;
                if(count == 1) {
                    if(split[0]) {
                        output.faces.push_back({ a, ab, c });
                        output.faces.push_back({ ab, b, c });
                    } else if(split[1]) {
                        output.faces.push_back({ b, bc, a });
                        output.faces.push_back({ bc, c, a });
                    } else {
                        output.faces.push_back({ c, ca, b });
                        output.faces.push_back({ ca, a, b });
                    }
                    continue;
                }
                if(count == 2) {
                    if(!split[0]) {
                        output.faces.push_back({ c, ca, bc });
                        output.faces.push_back({ ca, a, b });
                        output.faces.push_back({ ca, b, bc });
                    } else if(!split[1]) {
                        output.faces.push_back({ a, ab, ca });
                        output.faces.push_back({ ab, b, c });
                        output.faces.push_back({ ab, c, ca });
                    } else {
                        output.faces.push_back({ b, bc, ab });
                        output.faces.push_back({ bc, c, a });
                        output.faces.push_back({ bc, a, ab });
                    }
                    continue;
                }
                output.faces.push_back({ a, ab, ca });
                output.faces.push_back({ ab, b, bc });
                output.faces.push_back({ ca, bc, c });
                output.faces.push_back({ ab, bc, ca });
            }
            return output;
        }

        std::size_t trajectorySegmentAt(
            const std::vector<SprayPose>& trajectory,
            double time)
        {
            const auto next = std::upper_bound(trajectory.begin() + 1,
                trajectory.end(), time,
                [](double value, const SprayPose& pose) {
                    return value < pose.timeSeconds;
                });
            return std::min(static_cast<std::size_t>(
                next - trajectory.begin() - 1), trajectory.size() - 2);
        }

        GunState gunStateAt(
            const std::vector<SprayPose>& trajectory,
            double time)
        {
            const std::size_t segment = trajectorySegmentAt(trajectory, time);
            const SprayPose& first = trajectory[segment];
            const SprayPose& second = trajectory[segment + 1];
            const double duration = second.timeSeconds - first.timeSeconds;
            if(duration <= 0.0) {
                throw std::invalid_argument(
                    "Tzinava gun trajectory times must be strictly increasing.");
            }
            const double ratio = std::clamp(
                (time - first.timeSeconds) / duration, 0.0, 1.0);
            GunState state;
            state.position = (1.0 - ratio) * first.position + ratio * second.position;
            state.axis = ((1.0 - ratio) * first.axis + ratio * second.axis).normalized();
            state.velocity = (second.position - first.position) / duration;
            state.sprayEnabled = first.sprayEnabled;
            return state;
        }

        TriangleMesh rotatedMesh(
            const TriangleMesh& mesh,
            const Eigen::Vector3d& origin,
            const Eigen::Vector3d& axis,
            double angle)
        {
            TriangleMesh output = mesh;
            const Eigen::AngleAxisd rotation(angle, axis.normalized());
            for(Eigen::Vector3d& vertex : output.vertices) {
                vertex = origin + rotation * (vertex - origin);
            }
            return output;
        }

        void projectionBasis(
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

        double orientation(
            const Eigen::Vector2d& a,
            const Eigen::Vector2d& b,
            const Eigen::Vector2d& c)
        {
            return (b.x() - a.x()) * (c.y() - a.y())
                - (b.y() - a.y()) * (c.x() - a.x());
        }

        bool pointInsideTriangle(
            const Eigen::Vector2d& point,
            const std::array<Eigen::Vector2d, 3>& triangle)
        {
            if(std::abs(orientation(triangle[0], triangle[1], triangle[2]))
                <= kEpsilon) {
                return false;
            }
            const double first = orientation(triangle[0], triangle[1], point);
            const double second = orientation(triangle[1], triangle[2], point);
            const double third = orientation(triangle[2], triangle[0], point);
            return (first > kEpsilon && second > kEpsilon
                    && third > kEpsilon)
                || (first < -kEpsilon && second < -kEpsilon
                    && third < -kEpsilon);
        }

        bool edgesIntersect(const Eigen::Vector2d& first,
            const Eigen::Vector2d& second,
            const Eigen::Vector2d& third,
            const Eigen::Vector2d& fourth)
        {
            if(std::max(std::min(first.x(), second.x()),
                    std::min(third.x(), fourth.x()))
                    > std::min(std::max(first.x(), second.x()),
                        std::max(third.x(), fourth.x())) + kEpsilon
                || std::max(std::min(first.y(), second.y()),
                    std::min(third.y(), fourth.y()))
                    > std::min(std::max(first.y(), second.y()),
                        std::max(third.y(), fourth.y())) + kEpsilon) {
                return false;
            }
            const double a = orientation(first, second, third);
            const double b = orientation(first, second, fourth);
            const double c = orientation(third, fourth, first);
            const double d = orientation(third, fourth, second);
            return ((a > kEpsilon && b < -kEpsilon)
                    || (a < -kEpsilon && b > kEpsilon))
                && ((c > kEpsilon && d < -kEpsilon)
                    || (c < -kEpsilon && d > kEpsilon));
        }

        bool projectedTrianglesIntersect(
            const std::array<Eigen::Vector2d, 3>& first,
            const std::array<Eigen::Vector2d, 3>& second)
        {
            if(std::abs(orientation(first[0], first[1], first[2]))
                    <= kEpsilon
                || std::abs(orientation(second[0], second[1], second[2]))
                    <= kEpsilon) {
                return false;
            }
            for(int component = 0; component < 2; ++component) {
                const double firstMin = std::min({ first[0][component],
                    first[1][component], first[2][component] });
                const double firstMax = std::max({ first[0][component],
                    first[1][component], first[2][component] });
                const double secondMin = std::min({ second[0][component],
                    second[1][component], second[2][component] });
                const double secondMax = std::max({ second[0][component],
                    second[1][component], second[2][component] });
                if(firstMax < secondMin - kEpsilon
                    || secondMax < firstMin - kEpsilon) {
                    return false;
                }
            }
            for(int corner = 0; corner < 3; ++corner) {
                if(pointInsideTriangle(first[corner], second)
                    || pointInsideTriangle(second[corner], first)) {
                    return true;
                }
            }
            const Eigen::Vector2d firstCenter =
                (first[0] + first[1] + first[2]) / 3.0;
            const Eigen::Vector2d secondCenter =
                (second[0] + second[1] + second[2]) / 3.0;
            if(pointInsideTriangle(firstCenter, second)
                || pointInsideTriangle(secondCenter, first)) {
                return true;
            }
            for(int firstEdge = 0; firstEdge < 3; ++firstEdge) {
                for(int secondEdge = 0; secondEdge < 3; ++secondEdge) {
                    if(edgesIntersect(first[firstEdge], first[(firstEdge + 1) % 3],
                            second[secondEdge], second[(secondEdge + 1) % 3])) {
                        return true;
                    }
                }
            }
            return false;
        }

        double beamRadiusAt(
            const TzinavaParameters& parameters,
            double axialDistance)
        {
            if(parameters.beamKind == TzinavaBeamKind::Cylindrical
                || axialDistance <= parameters.cylindricalLengthMeters) {
                return parameters.beamRadiusMeters;
            }
            return parameters.beamRadiusMeters
                + (axialDistance - parameters.cylindricalLengthMeters)
                    * std::tan(parameters.coneHalfAngleRadians);
        }

        std::vector<Candidate> candidates(
            const TriangleMesh& mesh,
            const GunState& gun,
            const TzinavaInputModel& input,
            const TzinavaParameters& parameters,
            bool checkVisibility,
            std::size_t* visibilityQueryCount = nullptr)
        {
            const Eigen::Vector3d axis = gun.axis.normalized();
            Eigen::Vector3d basisX;
            Eigen::Vector3d basisY;
            projectionBasis(axis, basisX, basisY);
            std::vector<Candidate> output;
            output.reserve(mesh.faces.size());
            for(std::size_t faceIndex = 0; faceIndex < mesh.faces.size(); ++faceIndex) {
                const Eigen::Vector3d centroid = faceCentroid(mesh, faceIndex);
                const Eigen::Vector3d relative = centroid - gun.position;
                const double axialDistance = relative.dot(axis);
                if(axialDistance <= 0.0) {
                    continue;
                }
                const Eigen::Vector3d radial = relative - axialDistance * axis;
                const double beamRadius = beamRadiusAt(parameters, axialDistance);
                const auto& face = mesh.faces[faceIndex];
                if(radial.norm() > beamRadius) {
                    continue;
                }
                Candidate candidate;
                candidate.faceIndex = faceIndex;
                candidate.centroid = centroid;
                candidate.normal = faceNormal(mesh, faceIndex);
                candidate.facingCosine = -candidate.normal.dot(axis);
                candidate.standOffDistance = relative.norm();
                candidate.radialDistance = radial.norm();
                candidate.impactAngleDegrees = std::asin(std::clamp(
                    candidate.facingCosine, 0.0, 1.0)) * 180.0 / kPi;

                const Eigen::Vector3d rotationVelocity =
                    input.objectAngularSpeedRadiansPerSecond
                    * input.objectRotationAxis.normalized().cross(
                        centroid - input.objectRotationOrigin);
                const Eigen::Vector3d totalVelocity = gun.velocity - rotationVelocity;
                const Eigen::Vector3d planeNormal = totalVelocity.cross(axis);
                if(planeNormal.norm() > kEpsilon) {
                    const Eigen::Vector3d unitPlaneNormal = planeNormal.normalized();
                    const Eigen::Vector3d projectedNormal = candidate.normal
                        - candidate.normal.dot(unitPlaneNormal) * unitPlaneNormal;
                    candidate.projectedNormalDotAxis =
                        std::abs(projectedNormal.dot(axis));
                    if(candidate.projectedNormalDotAxis > kEpsilon) {
                        candidate.spotSpeedMetersPerSecond =
                            planeNormal.norm() / candidate.projectedNormalDotAxis;
                    }
                }
                for(int corner = 0; corner < 3; ++corner) {
                    const Eigen::Vector3d local =
                        mesh.vertices[face[corner]] - gun.position;
                    candidate.projected[corner] = Eigen::Vector2d(
                        local.dot(basisX), local.dot(basisY));
                }
                output.push_back(candidate);
            }
            if(!checkVisibility) {
                return output;
            }
            std::sort(output.begin(), output.end(),
                [](const Candidate& first, const Candidate& second) {
                    return first.standOffDistance < second.standOffDistance;
                });
            for(std::size_t targetIndex = 0;
                targetIndex < output.size(); ++targetIndex) {
                for(std::size_t otherIndex = 0;
                    otherIndex < targetIndex; ++otherIndex) {
                    if(visibilityQueryCount != nullptr) {
                        ++*visibilityQueryCount;
                    }
                    if(projectedTrianglesIntersect(
                            output[targetIndex].projected,
                            output[otherIndex].projected)) {
                        output[targetIndex].hidden = true;
                        break;
                    }
                }
            }
            return output;
        }

        double adaptiveTimeStep(
            const std::vector<Candidate>& visibleCandidates,
            const GunState& gun,
            const TzinavaParameters& parameters)
        {
            double timeStep = std::numeric_limits<double>::infinity();
            for(const Candidate& candidate : visibleCandidates) {
                if(candidate.hidden
                    || candidate.facingCosine <= kEpsilon
                    || candidate.spotSpeedMetersPerSecond <= kEpsilon
                    || candidate.projectedNormalDotAxis <= kEpsilon) {
                    continue;
                }
                timeStep = std::min(timeStep,
                    parameters.beamRadiusMeters
                        / (candidate.spotSpeedMetersPerSecond
                            * candidate.projectedNormalDotAxis));
            }
            if(gun.velocity.norm() > kEpsilon) {
                timeStep = std::min(timeStep,
                    parameters.beamRadiusMeters / gun.velocity.norm());
            }
            if(!std::isfinite(timeStep)) {
                return parameters.stationaryTimeStepSeconds;
            }
            return timeStep * parameters.timeStepOverlapFactor;
        }
    }

    void TzinavaLookupTable::validate() const
    {
        if(standOffDistancesMeters.size() < 2 || impactAnglesDegrees.size() < 2
            || thicknessMeters.size()
                != standOffDistancesMeters.size() * impactAnglesDegrees.size()
            || !std::is_sorted(standOffDistancesMeters.begin(),
                standOffDistancesMeters.end())
            || !std::is_sorted(impactAnglesDegrees.begin(),
                impactAnglesDegrees.end())) {
            throw std::invalid_argument(
                "Tzinava reproduction requires a complete, sorted SoD-angle "
                "thickness lookup table.");
        }
        for(const double value : thicknessMeters) {
            if(!std::isfinite(value) || value < 0.0) {
                throw std::invalid_argument(
                    "Tzinava thickness lookup table contains an invalid value.");
            }
        }
    }

    double TzinavaLookupTable::interpolate(
        double standOffDistanceMeters,
        double impactAngleDegrees) const
    {
        validate();
        if(standOffDistanceMeters < standOffDistancesMeters.front()
            || standOffDistanceMeters > standOffDistancesMeters.back()
            || impactAngleDegrees < impactAnglesDegrees.front()
            || impactAngleDegrees > impactAnglesDegrees.back()) {
            throw std::out_of_range(
                "Tzinava kinematics fall outside the supplied lookup table.");
        }
        const std::size_t distanceIndex = lowerInterval(
            standOffDistancesMeters, standOffDistanceMeters);
        const std::size_t angleIndex = lowerInterval(
            impactAnglesDegrees, impactAngleDegrees);
        const double distanceRatio =
            (standOffDistanceMeters - standOffDistancesMeters[distanceIndex])
            / (standOffDistancesMeters[distanceIndex + 1]
                - standOffDistancesMeters[distanceIndex]);
        const double angleRatio =
            (impactAngleDegrees - impactAnglesDegrees[angleIndex])
            / (impactAnglesDegrees[angleIndex + 1]
                - impactAnglesDegrees[angleIndex]);
        const auto value = [&](std::size_t distance, std::size_t angle) {
            return thicknessMeters[
                distance * impactAnglesDegrees.size() + angle];
        };
        const double low = (1.0 - angleRatio) * value(distanceIndex, angleIndex)
            + angleRatio * value(distanceIndex, angleIndex + 1);
        const double high =
            (1.0 - angleRatio) * value(distanceIndex + 1, angleIndex)
            + angleRatio * value(distanceIndex + 1, angleIndex + 1);
        return (1.0 - distanceRatio) * low + distanceRatio * high;
    }

    void TzinavaParameters::validate() const
    {
        thicknessTable.validate();
        if(!std::isfinite(beamRadiusMeters) || beamRadiusMeters <= 0.0
            || !std::isfinite(speedCoefficientB) || speedCoefficientB <= 0.0
            || !std::isfinite(speedCoefficientC) || speedCoefficientC <= 0.0
            || !std::isfinite(referenceSpotSpeedMillimetersPerSecond)
            || referenceSpotSpeedMillimetersPerSecond <= 0.0
            || !std::isfinite(timeStepOverlapFactor)
            || timeStepOverlapFactor <= 0.0
            || !std::isfinite(stationaryTimeStepSeconds)
            || stationaryTimeStepSeconds <= 0.0
            || !std::isfinite(lookupReferenceDwellSeconds)
            || lookupReferenceDwellSeconds <= 0.0) {
            throw std::invalid_argument(
                "Tzinava reproduction parameters are incomplete.");
        }
        if(gaussianRadialProfile
            && (!std::isfinite(gaussianSigmaMeters)
                || gaussianSigmaMeters <= 0.0)) {
            throw std::invalid_argument(
                "Tzinava Gaussian radial profile requires sigma.");
        }
        if(beamKind == TzinavaBeamKind::CylindricalConical
            && (cylindricalLengthMeters < 0.0
                || coneHalfAngleRadians <= 0.0)) {
            throw std::invalid_argument(
                "Tzinava cylindrical-conical beam parameters are invalid.");
        }
    }

    TzinavaResult TzinavaReproducer::run(
        const TzinavaInputModel& input,
        const TzinavaParameters& parameters,
        const ReproductionExecution& execution)
    {
        reportDiagnostic(execution, "Input validation",
            "mesh faces=" + std::to_string(input.initialMesh.faces.size())
                + ", gun poses=" + std::to_string(input.gunTrajectory.size()));
        parameters.validate();
        input.initialMesh.validate();
        if(input.gunTrajectory.size() < 2) {
            throw std::invalid_argument(
                "Tzinava reproduction requires at least two gun poses.");
        }
        if(std::adjacent_find(input.gunTrajectory.begin(),
                input.gunTrajectory.end(),
                [](const SprayPose& first, const SprayPose& second) {
                    return second.timeSeconds <= first.timeSeconds;
                }) != input.gunTrajectory.end()) {
            throw std::invalid_argument(
                "Tzinava gun trajectory times must be strictly increasing.");
        }
        if(input.objectRotationAxis.squaredNorm() <= kEpsilon) {
            throw std::invalid_argument("Tzinava object rotation axis is zero.");
        }

        const auto started = std::chrono::steady_clock::now();
        TzinavaResult result;
        reportDiagnostic(execution, "Mesh subdivision",
            "initial faces=" + std::to_string(input.initialMesh.faces.size()));
        result.subdividedMesh = subdivide(
            input.initialMesh, parameters.beamRadiusMeters);
        reportDiagnostic(execution, "Candidate search and deposition",
            "subdivided faces="
                + std::to_string(result.subdividedMesh.faces.size()));
        result.faceThicknessMeters.assign(
            result.subdividedMesh.faces.size(), 0.0);
        result.faceCentroids.resize(result.subdividedMesh.faces.size());
        std::vector<double> passagePeak(result.subdividedMesh.faces.size(), 0.0);
        std::vector<std::size_t> lastVisibleStep(
            result.subdividedMesh.faces.size(), 0);
        std::vector<std::size_t> lastVisiblePass(
            result.subdividedMesh.faces.size(), 0);
        result.statistics.evaluatedElementCount =
            result.subdividedMesh.faces.size();

        double time = input.gunTrajectory.front().timeSeconds;
        const double endTime = input.gunTrajectory.back().timeSeconds;
        std::size_t previousSegment = 0;
        std::size_t movingPass = 0;
        while(time < endTime - kEpsilon) {
            if(canceled(execution)) {
                result.canceled = true;
                break;
            }
            const std::size_t segment = trajectorySegmentAt(
                input.gunTrajectory, time);
            for(std::size_t boundary = previousSegment + 1;
                boundary <= segment; ++boundary) {
                const SprayPose& before = input.gunTrajectory[boundary - 1];
                const SprayPose& turn = input.gunTrajectory[boundary];
                const SprayPose& after = input.gunTrajectory[boundary + 1];
                const Eigen::Vector3d inbound = turn.position - before.position;
                const Eigen::Vector3d outbound = after.position - turn.position;
                if((!before.sprayEnabled && turn.sprayEnabled)
                    || (inbound.squaredNorm() > kEpsilon
                        && outbound.squaredNorm() > kEpsilon
                        && inbound.dot(outbound) < 0.0)) {
                    ++movingPass;
                }
            }
            previousSegment = segment;
            const GunState gun = gunStateAt(input.gunTrajectory, time);
            const TriangleMesh currentMesh = rotatedMesh(
                result.subdividedMesh,
                input.objectRotationOrigin,
                input.objectRotationAxis,
                input.objectAngularSpeedRadiansPerSecond
                    * (time - input.gunTrajectory.front().timeSeconds));
            const std::vector<Candidate> currentCandidates = candidates(
                currentMesh, gun, input, parameters, true,
                &result.statistics.visibilityQueryCount);
            double dt = adaptiveTimeStep(currentCandidates, gun, parameters);
            dt = std::min({ dt,
                input.gunTrajectory[segment + 1].timeSeconds - time,
                endTime - time });
            if(!std::isfinite(dt) || dt <= 0.0) {
                throw std::runtime_error(
                    "Tzinava adaptive time step is not positive.");
            }

            if(gun.sprayEnabled) {
                for(const Candidate& candidate : currentCandidates) {
                    ++result.statistics.candidatePairCount;
                    if(candidate.hidden) {
                        ++result.statistics.hiddenElementCount;
                        continue;
                    }
                    if(candidate.facingCosine <= kEpsilon) {
                        continue;
                    }
                    const bool stationary =
                        candidate.spotSpeedMetersPerSecond <= kEpsilon;
                    const double baseThickness = parameters.thicknessTable.interpolate(
                        candidate.standOffDistance,
                        candidate.impactAngleDegrees);
                    const double radialFactor = parameters.gaussianRadialProfile
                        ? std::exp(-0.5 * candidate.radialDistance
                            * candidate.radialDistance
                            / (parameters.gaussianSigmaMeters
                                * parameters.gaussianSigmaMeters))
                        : 1.0;
                    double increment = baseThickness * radialFactor;
                    if(stationary) {
                        increment *= dt / parameters.lookupReferenceDwellSeconds;
                    } else {
                        const double speedMillimetersPerSecond =
                            1000.0 * candidate.spotSpeedMetersPerSecond;
                        const double speedFunction = 1.0
                            / (parameters.speedCoefficientB
                                    * speedMillimetersPerSecond
                                + parameters.speedCoefficientC);
                        const double referenceSpeedFunction = 1.0
                            / (parameters.speedCoefficientB
                                    * parameters.referenceSpotSpeedMillimetersPerSecond
                                + parameters.speedCoefficientC);
                        increment *= speedFunction / referenceSpeedFunction;
                        const std::size_t face = candidate.faceIndex;
                        const std::size_t step =
                            result.statistics.trajectorySampleCount;
                        if(step == 0 || lastVisibleStep[face] != step - 1
                            || lastVisiblePass[face] != movingPass) {
                            passagePeak[face] = 0.0;
                        }
                        // Eq. (9) is one full moving pass. The paper does not specify
                        // which overlapping sample represents it; retain its peak.
                        const double additional =
                            std::max(0.0, increment - passagePeak[face]);
                        passagePeak[face] = std::max(passagePeak[face], increment);
                        lastVisibleStep[face] = step;
                        lastVisiblePass[face] = movingPass;
                        increment = additional;
                    }
                    result.faceThicknessMeters[candidate.faceIndex] += increment;
                }
            }
            ++result.statistics.trajectorySampleCount;
            time += dt;
            reportProgress(execution,
                (time - input.gunTrajectory.front().timeSeconds)
                    / (endTime - input.gunTrajectory.front().timeSeconds),
                "Tzinava 2020 reproduction");
        }

        reportDiagnostic(execution, "Result assembly",
            "evaluated time steps="
                + std::to_string(result.statistics.trajectorySampleCount));
        for(std::size_t faceIndex = 0;
            faceIndex < result.subdividedMesh.faces.size(); ++faceIndex) {
            result.faceCentroids[faceIndex] =
                faceCentroid(result.subdividedMesh, faceIndex);
        }
        result.statistics.elapsedMilliseconds = elapsedMilliseconds(started);
        return result;
    }
}
