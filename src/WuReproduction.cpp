#include <SprayThicknessPrediction/WuReproduction.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace spraythickness::published
{
    namespace
    {
        constexpr double kPi = 3.14159265358979323846;
        constexpr double kEpsilon = 1.0e-12;

        struct SurfaceHit
        {
            double distance{ std::numeric_limits<double>::infinity() };
            Eigen::Vector3d position = Eigen::Vector3d::Zero();
            Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
            std::size_t substrateFaceIndex{ static_cast<std::size_t>(-1) };
            bool depositedSurface{ false };
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

        std::vector<Eigen::Vector3d> rayDirections(
            const Eigen::Vector3d& axis,
            const WuParameters& parameters)
        {
            Eigen::Vector3d basisX;
            Eigen::Vector3d basisY;
            localBasis(axis, basisX, basisY);
            std::vector<Eigen::Vector3d> directions{ axis };
            for(double polar = parameters.rayAngularStepRadians;
                polar <= parameters.maximumDeflectionRadians + kEpsilon;
                polar += parameters.rayAngularStepRadians) {
                const std::size_t azimuthCount = std::max<std::size_t>(6,
                    static_cast<std::size_t>(std::ceil(
                        2.0 * kPi * std::sin(polar)
                        / parameters.rayAngularStepRadians)));
                for(std::size_t azimuthIndex = 0;
                    azimuthIndex < azimuthCount; ++azimuthIndex) {
                    const double azimuth = 2.0 * kPi
                        * static_cast<double>(azimuthIndex)
                        / static_cast<double>(azimuthCount);
                    directions.push_back((std::cos(polar) * axis
                        + std::sin(polar)
                            * (std::cos(azimuth) * basisX
                                + std::sin(azimuth) * basisY)).normalized());
                }
            }
            return directions;
        }

        bool intersectCylinder(
            const Eigen::Vector3d& origin,
            const Eigen::Vector3d& direction,
            const WuDepositedCylinder& cylinder,
            const Eigen::Vector3d& substrateNormal,
            SurfaceHit& hit)
        {
            const Eigen::Vector3d axis = cylinder.growthDirection.normalized();
            const double normalAlongAxis = substrateNormal.dot(axis);
            if(normalAlongAxis <= kEpsilon) {
                return false;
            }
            const Eigen::Vector3d delta = origin - cylinder.baseCenter;
            const Eigen::Vector3d directionPerpendicular =
                direction - direction.dot(substrateNormal) / normalAlongAxis * axis;
            const Eigen::Vector3d deltaPerpendicular =
                delta - delta.dot(substrateNormal) / normalAlongAxis * axis;
            const double a = directionPerpendicular.squaredNorm();
            const double b = 2.0 * directionPerpendicular.dot(deltaPerpendicular);
            const double c = deltaPerpendicular.squaredNorm()
                - cylinder.radiusMeters * cylinder.radiusMeters;
            double nearest = std::numeric_limits<double>::infinity();
            Eigen::Vector3d nearestNormal = Eigen::Vector3d::Zero();
            if(a > kEpsilon) {
                const double discriminant = b * b - 4.0 * a * c;
                if(discriminant >= 0.0) {
                    for(const double sign : { -1.0, 1.0 }) {
                        const double distance =
                            (-b + sign * std::sqrt(discriminant)) / (2.0 * a);
                        if(distance <= kEpsilon || distance >= nearest) {
                            continue;
                        }
                        const Eigen::Vector3d point = origin + distance * direction;
                        const double axial = (point - cylinder.baseCenter)
                            .dot(substrateNormal) / normalAlongAxis;
                        if(axial >= 0.0 && axial <= cylinder.heightMeters) {
                            nearest = distance;
                            const Eigen::Vector3d radial = point
                                - cylinder.baseCenter - axial * axis;
                            nearestNormal = (radial - substrateNormal
                                * radial.dot(axis) / normalAlongAxis).normalized();
                        }
                    }
                }
            }
            const double denominator = direction.dot(substrateNormal);
            if(std::abs(denominator) > kEpsilon) {
                for(const double axial : { 0.0, cylinder.heightMeters }) {
                    const Eigen::Vector3d center = cylinder.baseCenter + axial * axis;
                    const double distance = (center - origin)
                        .dot(substrateNormal) / denominator;
                    if(distance <= kEpsilon || distance >= nearest) {
                        continue;
                    }
                    const Eigen::Vector3d point = origin + distance * direction;
                    if((point - center).squaredNorm()
                        <= cylinder.radiusMeters * cylinder.radiusMeters) {
                        nearest = distance;
                        nearestNormal = axial > 0.0
                            ? substrateNormal : -substrateNormal;
                    }
                }
            }
            if(!std::isfinite(nearest)) {
                return false;
            }
            hit.distance = nearest;
            hit.position = origin + nearest * direction;
            hit.normal = nearestNormal;
            hit.substrateFaceIndex = cylinder.substrateFaceIndex;
            hit.depositedSurface = true;
            return true;
        }

        class DepositedCylinderIndex
        {
        public:
            explicit DepositedCylinderIndex(
                const TriangleMesh& substrate,
                const std::vector<WuDepositedCylinder>& cylinders)
                : m_substrate(substrate), m_cylinders(cylinders)
            {
            }

            void add()
            {
                const WuDepositedCylinder& cylinder = m_cylinders.back();
                const Eigen::Vector3d axis = cylinder.growthDirection.normalized();
                const Eigen::Vector3d normal = faceNormal(
                    m_substrate, cylinder.substrateFaceIndex);
                m_normals.push_back(normal);
                const Eigen::Vector3d top = cylinder.baseCenter
                    + cylinder.heightMeters * axis;
                Node bounds;
                for(int component = 0; component < 3; ++component) {
                    const double radius = cylinder.radiusMeters * std::sqrt(
                        std::max(0.0, 1.0 - normal[component] * normal[component]))
                        + kEpsilon;
                    bounds.lower[component] = std::min(
                        cylinder.baseCenter[component], top[component]) - radius;
                    bounds.upper[component] = std::max(
                        cylinder.baseCenter[component], top[component]) + radius;
                }
                m_bounds.push_back(bounds);
                if(m_cylinders.size() - m_indexedCount < 512) {
                    return;
                }
                m_indexedCount = m_cylinders.size();
                m_order.resize(m_indexedCount);
                std::iota(m_order.begin(), m_order.end(), std::size_t{ 0 });
                m_nodes.clear();
                m_nodes.reserve(m_indexedCount / 4 + 1);
                build(0, m_indexedCount);
            }

            void nearest(const Eigen::Vector3d& origin,
                const Eigen::Vector3d& direction, SurfaceHit& hit) const
            {
                std::size_t bestCylinder = 0;
                if(!m_nodes.empty()) {
                    search(0, origin, direction, hit, bestCylinder);
                }
                for(std::size_t index = m_indexedCount;
                    index < m_cylinders.size(); ++index) {
                    SurfaceHit candidate;
                    if(intersectCylinder(origin, direction,
                            m_cylinders[index], m_normals[index], candidate)) {
                        acceptHit(index, candidate, hit, bestCylinder);
                    }
                }
            }

        private:
            struct Node
            {
                Eigen::Vector3d lower = Eigen::Vector3d::Zero();
                Eigen::Vector3d upper = Eigen::Vector3d::Zero();
                std::size_t begin{ 0 };
                std::size_t end{ 0 };
                std::size_t left{ 0 };
                std::size_t right{ 0 };
            };

            std::size_t build(std::size_t begin, std::size_t end)
            {
                Node node;
                node.begin = begin;
                node.end = end;
                node.lower.setConstant(std::numeric_limits<double>::infinity());
                node.upper.setConstant(-std::numeric_limits<double>::infinity());
                for(std::size_t offset = begin; offset < end; ++offset) {
                    const Node& cylinder = m_bounds[m_order[offset]];
                    node.lower = node.lower.cwiseMin(cylinder.lower);
                    node.upper = node.upper.cwiseMax(cylinder.upper);
                }
                const std::size_t index = m_nodes.size();
                m_nodes.push_back(node);
                if(end - begin <= 8) {
                    return index;
                }
                Eigen::Index axis = 0;
                (node.upper - node.lower).maxCoeff(&axis);
                const std::size_t middle = begin + (end - begin) / 2;
                std::nth_element(m_order.begin() + begin,
                    m_order.begin() + middle, m_order.begin() + end,
                    [this, axis](std::size_t first, std::size_t second) {
                        return m_bounds[first].lower[axis]
                            + m_bounds[first].upper[axis]
                            < m_bounds[second].lower[axis]
                                + m_bounds[second].upper[axis];
                    });
                const std::size_t left = build(begin, middle);
                const std::size_t right = build(middle, end);
                m_nodes[index].left = left;
                m_nodes[index].right = right;
                return index;
            }

            static double rayBoxNear(const Node& node,
                const Eigen::Vector3d& origin,
                const Eigen::Vector3d& direction, double limit)
            {
                double near = 0.0;
                double far = limit;
                for(int component = 0; component < 3; ++component) {
                    if(std::abs(direction[component]) <= 1.0e-15) {
                        if(origin[component] < node.lower[component]
                            || origin[component] > node.upper[component]) {
                            return std::numeric_limits<double>::infinity();
                        }
                        continue;
                    }
                    const double first = (node.lower[component]
                        - origin[component]) / direction[component];
                    const double second = (node.upper[component]
                        - origin[component]) / direction[component];
                    near = std::max(near, std::min(first, second));
                    far = std::min(far, std::max(first, second));
                    if(near > far) {
                        return std::numeric_limits<double>::infinity();
                    }
                }
                return far > kEpsilon ? near
                    : std::numeric_limits<double>::infinity();
            }

            static void acceptHit(std::size_t cylinderIndex,
                const SurfaceHit& candidate, SurfaceHit& hit,
                std::size_t& bestCylinder)
            {
                if(candidate.distance < hit.distance
                    || (candidate.distance == hit.distance
                        && bestCylinder != 0
                        && cylinderIndex + 1 < bestCylinder)) {
                    hit = candidate;
                    bestCylinder = cylinderIndex + 1;
                }
            }

            void search(std::size_t index,
                const Eigen::Vector3d& origin,
                const Eigen::Vector3d& direction, SurfaceHit& hit,
                std::size_t& bestCylinder) const
            {
                const Node& node = m_nodes[index];
                if(rayBoxNear(node, origin, direction, hit.distance)
                    > hit.distance) {
                    return;
                }
                if(node.left == 0) {
                    for(std::size_t offset = node.begin; offset < node.end;
                        ++offset) {
                        SurfaceHit candidate;
                        const std::size_t cylinderIndex = m_order[offset];
                        if(intersectCylinder(origin, direction,
                                m_cylinders[cylinderIndex],
                                m_normals[cylinderIndex], candidate)) {
                            acceptHit(cylinderIndex, candidate, hit,
                                bestCylinder);
                        }
                    }
                    return;
                }
                const double left = rayBoxNear(m_nodes[node.left],
                    origin, direction, hit.distance);
                const double right = rayBoxNear(m_nodes[node.right],
                    origin, direction, hit.distance);
                if(left <= right) {
                    search(node.left, origin, direction, hit, bestCylinder);
                    search(node.right, origin, direction, hit, bestCylinder);
                } else {
                    search(node.right, origin, direction, hit, bestCylinder);
                    search(node.left, origin, direction, hit, bestCylinder);
                }
            }

            const TriangleMesh& m_substrate;
            const std::vector<WuDepositedCylinder>& m_cylinders;
            std::vector<Eigen::Vector3d> m_normals;
            std::vector<Node> m_bounds;
            std::vector<std::size_t> m_order;
            std::vector<Node> m_nodes;
            std::size_t m_indexedCount{ 0 };
        };

        bool nearestSurfaceHit(
            const TriangleMesh& substrate,
            const DepositedCylinderIndex& cylinders,
            const Eigen::Vector3d& origin,
            const Eigen::Vector3d& direction,
            SurfaceHit& hit)
        {
            RayHit meshHit;
            if(nearestRayHit(substrate, origin, direction, meshHit)) {
                hit.distance = meshHit.distance;
                hit.position = meshHit.position;
                hit.normal = meshHit.normal;
                hit.substrateFaceIndex = meshHit.faceIndex;
            }
            cylinders.nearest(origin, direction, hit);
            return std::isfinite(hit.distance);
        }

    }

    void WuParameters::validate() const
    {
        sprayAngleRelativeDepositionEfficiency.validate(5,
            "Wu spray-angle RDE");
        sprayDistanceRelativeDepositionEfficiency.validate(4,
            "Wu spray-distance RDE");
        traverseSpeedPeakCorrectionFactor.validate(2,
            "Wu traverse-speed PCF");
        if(!std::isfinite(peakCylinderHeightMeters)
            || peakCylinderHeightMeters <= 0.0
            || !std::isfinite(referencePoseDurationSeconds)
            || referencePoseDurationSeconds <= 0.0
            || !std::isfinite(gaussianSigmaMeters)
            || gaussianSigmaMeters <= 0.0
            || !std::isfinite(maximumDeflectionRadians)
            || maximumDeflectionRadians <= 0.0
            || !std::isfinite(rayAngularStepRadians)
            || rayAngularStepRadians <= 0.0
            || rayAngularStepRadians > maximumDeflectionRadians
            || !std::isfinite(cylinderRadiusMeters)
            || cylinderRadiusMeters <= 0.0) {
            throw std::invalid_argument(
                "Wu reproduction requires the measured Gaussian profile, ray "
                "discretization, cylinder radius, and all fitted coefficients.");
        }
    }

    WuResult WuReproducer::run(
        const WuInputModel& input,
        const WuParameters& parameters,
        const ReproductionExecution& execution)
    {
        const auto started = std::chrono::steady_clock::now();
        reportDiagnostic(execution, "Input validation",
            "substrate faces=" + std::to_string(input.substrate.faces.size())
                + ", nozzle poses="
                + std::to_string(input.nozzleTrajectory.size()));
        parameters.validate();
        input.substrate.validate();
        if(input.nozzleTrajectory.size() < 2) {
            throw std::invalid_argument(
                "Wu reproduction requires at least two nozzle poses.");
        }
        const std::vector<SprayPose>& poses = input.nozzleTrajectory;
        WuParameters sampling = parameters;
        if(input.generatedPlateStack) {
            sampling.rayAngularStepRadians = std::min(
                sampling.rayAngularStepRadians, 0.007);
        }
        const double topPlateZ = input.generatedPlateStack
            ? std::max_element(input.substrate.vertices.begin(),
                input.substrate.vertices.end(), [](const auto& first,
                    const auto& second) { return first.z() < second.z(); })->z()
            : 0.0;
        reportDiagnostic(execution, "Trajectory sampling",
            "original poses=" + std::to_string(poses.size())
                + ", inserted=0");
        WuResult result;
        result.substrate = input.substrate;
        DepositedCylinderIndex cylinderIndex(result.substrate,
            result.depositedCylinders);
        result.statistics.trajectorySampleCount = poses.size();
        result.statistics.evaluatedElementCount = input.substrate.faces.size();
        reportDiagnostic(execution, "Ray visibility and cylinder deposition",
            "original poses=" + std::to_string(poses.size()));
        std::size_t substrateHits = 0;
        std::size_t coatingHits = 0;
        std::size_t missedRays = 0;
        std::size_t behindNozzle = 0;
        std::size_t rejectedAngle = 0;
        std::size_t rejectedBackFace = 0;
        std::size_t rejectedDistance = 0;
        std::size_t rejectedHeight = 0;
        std::size_t disabledPoses = 0;
        std::size_t invalidAxisPoses = 0;
        bool reportedPlateSampling = false;

        double sprayDurationSeconds = 0.0;
        for(std::size_t poseIndex = 0; poseIndex + 1 < poses.size(); ++poseIndex) {
            if(canceled(execution)) {
                result.canceled = true;
                break;
            }
            const SprayPose& pose = poses[poseIndex];
            const double durationSeconds = poses[poseIndex + 1].timeSeconds
                - pose.timeSeconds;
            if(!std::isfinite(durationSeconds) || durationSeconds < 0.0) {
                throw std::invalid_argument(
                    "Wu reproduction requires finite, nondecreasing pose times.");
            }
            if(durationSeconds == 0.0) {
                continue;
            }
            if(!pose.sprayEnabled) {
                ++disabledPoses;
                continue;
            }
            if(pose.axis.squaredNorm() <= kEpsilon) {
                ++invalidAxisPoses;
                continue;
            }
            sprayDurationSeconds += durationSeconds;
            const double speedMillimetersPerSecond =
                pose.linearVelocity.norm() * 1000.0;
            const double speedPcf =
                parameters.traverseSpeedPeakCorrectionFactor.evaluate(
                    speedMillimetersPerSecond);
            if(!std::isfinite(speedPcf) || speedPcf <= 0.0) {
                throw std::out_of_range(
                    "Wu traverse-speed correction is non-positive at "
                    + std::to_string(speedMillimetersPerSecond)
                    + " mm/s (factor=" + std::to_string(speedPcf)
                    + "); use a speed supported by the calibration.");
            }
            const Eigen::Vector3d mainAxis = pose.axis.normalized();
            const std::vector<Eigen::Vector3d> directions =
                rayDirections(mainAxis, sampling);
            double cylinderRadius = parameters.cylinderRadiusMeters;
            if(input.generatedPlateStack && std::abs(mainAxis.z()) > 0.25) {
                const double standOff = std::abs(
                    (pose.position.z() - topPlateZ) / mainAxis.z());
                cylinderRadius = std::max(cylinderRadius,
                    0.75 * standOff * sampling.rayAngularStepRadians);
            }
            if(input.generatedPlateStack && !reportedPlateSampling) {
                reportDiagnostic(execution, "Plate sampling",
                    "rays/pose=" + std::to_string(directions.size())
                        + ", angular step="
                        + std::to_string(sampling.rayAngularStepRadians)
                        + " rad, physical cylinder radius="
                        + std::to_string(cylinderRadius * 1000.0) + " mm");
                reportedPlateSampling = true;
            }
            std::vector<WuDepositedCylinder> pending;
            pending.reserve(directions.size());
            for(const Eigen::Vector3d& direction : directions) {
                SurfaceHit hit;
                ++result.statistics.visibilityQueryCount;
                if(!nearestSurfaceHit(result.substrate,
                        cylinderIndex,
                        pose.position,
                        direction,
                        hit)) {
                    ++missedRays;
                    continue;
                }
                if(hit.depositedSurface) {
                    ++coatingHits;
                } else {
                    ++substrateHits;
                }
                const Eigen::Vector3d relative = hit.position - pose.position;
                const double axial = relative.dot(mainAxis);
                if(axial <= 0.0) {
                    ++behindNozzle;
                    continue;
                }
                const double radial = (relative - axial * mainAxis).norm();
                // The fitted angle efficiency describes the substrate, not the
                // microscopic side wall of a previously deposited cylinder.
                const Eigen::Vector3d substrateNormal = faceNormal(
                    result.substrate, hit.substrateFaceIndex);
                if(substrateNormal.dot(direction) >= -kEpsilon) {
                    ++rejectedBackFace;
                    continue;
                }
                const double sprayAngleDegrees = std::asin(std::clamp(
                    -substrateNormal.dot(direction),
                    0.0,
                    1.0)) * 180.0 / kPi;
                const double angleRde =
                    parameters.sprayAngleRelativeDepositionEfficiency.evaluate(
                        sprayAngleDegrees);
                const double distanceRde =
                    parameters.sprayDistanceRelativeDepositionEfficiency.evaluate(
                        relative.norm() * 1000.0);
                if(!std::isfinite(angleRde) || angleRde <= 0.0) {
                    ++rejectedAngle;
                    continue;
                }
                if(!std::isfinite(distanceRde) || distanceRde <= 0.0) {
                    ++rejectedDistance;
                    continue;
                }
                const double height = parameters.peakCylinderHeightMeters
                    * (durationSeconds / parameters.referencePoseDurationSeconds)
                    * std::exp(-0.5 * radial * radial
                        / (parameters.gaussianSigmaMeters
                            * parameters.gaussianSigmaMeters))
                    * angleRde * distanceRde * speedPcf;
                if(height <= 0.0 || !std::isfinite(height)) {
                    ++rejectedHeight;
                    continue;
                }
                pending.push_back({
                    hit.position,
                    -direction,
                    cylinderRadius,
                    height,
                    hit.substrateFaceIndex
                });
            }
            for(const WuDepositedCylinder& cylinder : pending) {
                result.depositedCylinders.push_back(cylinder);
                cylinderIndex.add();
                ++result.statistics.candidatePairCount;
            }
            reportProgress(execution,
                static_cast<double>(poseIndex + 1)
                    / static_cast<double>(poses.size() - 1),
                "Wu 2020 reproduction");
        }
        if(!result.canceled) {
            reportProgress(execution, 1.0, "Wu 2020 reproduction");
        }
        reportDiagnostic(execution, "Ray intersections",
            "rays=" + std::to_string(result.statistics.visibilityQueryCount)
                + ", substrate=" + std::to_string(substrateHits)
                + ", coating=" + std::to_string(coatingHits)
                + ", missed=" + std::to_string(missedRays));
        reportDiagnostic(execution, "Deposition decisions",
            "deposited=" + std::to_string(result.depositedCylinders.size())
                + ", behindNozzle=" + std::to_string(behindNozzle)
                + ", backFace=" + std::to_string(rejectedBackFace)
                + ", angle=" + std::to_string(rejectedAngle)
                + ", distance=" + std::to_string(rejectedDistance)
                + ", height=" + std::to_string(rejectedHeight)
                + ", sprayOffPoses=" + std::to_string(disabledPoses)
                + ", invalidAxisPoses=" + std::to_string(invalidAxisPoses));
        reportDiagnostic(execution, "Result assembly",
            "visibility rays="
                + std::to_string(result.statistics.visibilityQueryCount)
                + ", deposited cylinders="
                + std::to_string(result.depositedCylinders.size())
                + ", active spray duration="
                + std::to_string(sprayDurationSeconds) + " s");
        result.statistics.elapsedMilliseconds = elapsedMilliseconds(started);
        return result;
    }
}
