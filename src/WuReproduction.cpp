#include <SprayThicknessPrediction/WuReproduction.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
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
            SurfaceHit& hit)
        {
            const Eigen::Vector3d axis = cylinder.growthDirection.normalized();
            const Eigen::Vector3d delta = origin - cylinder.baseCenter;
            const Eigen::Vector3d directionPerpendicular =
                direction - direction.dot(axis) * axis;
            const Eigen::Vector3d deltaPerpendicular =
                delta - delta.dot(axis) * axis;
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
                        const double axial = (point - cylinder.baseCenter).dot(axis);
                        if(axial >= 0.0 && axial <= cylinder.heightMeters) {
                            nearest = distance;
                            nearestNormal = (point - cylinder.baseCenter
                                - axial * axis).normalized();
                        }
                    }
                }
            }
            const double denominator = direction.dot(axis);
            if(std::abs(denominator) > kEpsilon) {
                for(const double axial : { 0.0, cylinder.heightMeters }) {
                    const Eigen::Vector3d center = cylinder.baseCenter + axial * axis;
                    const double distance = (center - origin).dot(axis) / denominator;
                    if(distance <= kEpsilon || distance >= nearest) {
                        continue;
                    }
                    const Eigen::Vector3d point = origin + distance * direction;
                    if((point - center).squaredNorm()
                        <= cylinder.radiusMeters * cylinder.radiusMeters) {
                        nearest = distance;
                        nearestNormal = axial > 0.0 ? axis : -axis;
                    }
                }
            }
            if(!std::isfinite(nearest)) {
                return false;
            }
            hit.distance = nearest;
            hit.position = origin + nearest * direction;
            hit.normal = nearestNormal;
            return true;
        }

        bool nearestSurfaceHit(
            const TriangleMesh& substrate,
            const std::vector<WuDepositedCylinder>& cylinders,
            const Eigen::Vector3d& origin,
            const Eigen::Vector3d& direction,
            SurfaceHit& hit)
        {
            RayHit meshHit;
            if(nearestRayHit(substrate, origin, direction, meshHit)) {
                hit.distance = meshHit.distance;
                hit.position = meshHit.position;
                hit.normal = meshHit.normal;
            }
            for(const WuDepositedCylinder& cylinder : cylinders) {
                SurfaceHit cylinderHit;
                if(intersectCylinder(origin, direction, cylinder, cylinderHit)
                    && cylinderHit.distance < hit.distance) {
                    hit = cylinderHit;
                }
            }
            return std::isfinite(hit.distance);
        }

        std::vector<SprayPose> spatiallyResample(
            const std::vector<SprayPose>& poses,
            double maximumSpacing)
        {
            if(poses.size() < 2) {
                throw std::invalid_argument(
                    "Wu reproduction requires at least two nozzle poses.");
            }
            std::vector<SprayPose> result;
            for(std::size_t index = 0; index + 1 < poses.size(); ++index) {
                const SprayPose& first = poses[index];
                const SprayPose& second = poses[index + 1];
                const double distance = (second.position - first.position).norm();
                const std::size_t intervalCount = std::max<std::size_t>(1,
                    static_cast<std::size_t>(std::ceil(distance / maximumSpacing)));
                for(std::size_t interval = 0; interval < intervalCount; ++interval) {
                    const double ratio = static_cast<double>(interval)
                        / static_cast<double>(intervalCount);
                    SprayPose pose;
                    pose.timeSeconds = (1.0 - ratio) * first.timeSeconds
                        + ratio * second.timeSeconds;
                    pose.position = (1.0 - ratio) * first.position
                        + ratio * second.position;
                    pose.axis = ((1.0 - ratio) * first.axis
                        + ratio * second.axis).normalized();
                    const double duration = second.timeSeconds - first.timeSeconds;
                    if(duration > 0.0) {
                        pose.linearVelocity =
                            (second.position - first.position) / duration;
                    } else {
                        pose.linearVelocity.setZero();
                    }
                    pose.sprayEnabled = first.sprayEnabled;
                    result.push_back(pose);
                }
            }
            result.push_back(poses.back());
            return result;
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
        parameters.validate();
        input.substrate.validate();
        const std::vector<SprayPose> poses = spatiallyResample(
            input.nozzleTrajectory, 0.5 * parameters.gaussianSigmaMeters);
        const auto started = std::chrono::steady_clock::now();
        WuResult result;
        result.substrate = input.substrate;
        result.statistics.trajectorySampleCount = poses.size();
        result.statistics.evaluatedElementCount = input.substrate.faces.size();

        for(std::size_t poseIndex = 0; poseIndex + 1 < poses.size(); ++poseIndex) {
            if(canceled(execution)) {
                result.canceled = true;
                break;
            }
            const SprayPose& pose = poses[poseIndex];
            if(!pose.sprayEnabled || pose.axis.squaredNorm() <= kEpsilon) {
                continue;
            }
            const Eigen::Vector3d mainAxis = pose.axis.normalized();
            const std::vector<Eigen::Vector3d> directions =
                rayDirections(mainAxis, parameters);
            for(const Eigen::Vector3d& direction : directions) {
                SurfaceHit hit;
                ++result.statistics.visibilityQueryCount;
                if(!nearestSurfaceHit(result.substrate,
                        result.depositedCylinders,
                        pose.position,
                        direction,
                        hit)) {
                    continue;
                }
                const Eigen::Vector3d relative = hit.position - pose.position;
                const double axial = relative.dot(mainAxis);
                if(axial <= 0.0) {
                    continue;
                }
                const double radial = (relative - axial * mainAxis).norm();
                const double sprayAngleDegrees = std::asin(std::clamp(
                    std::abs(hit.normal.normalized().dot(direction)),
                    0.0,
                    1.0)) * 180.0 / kPi;
                const double angleRde =
                    parameters.sprayAngleRelativeDepositionEfficiency.evaluate(
                        sprayAngleDegrees);
                const double distanceRde =
                    parameters.sprayDistanceRelativeDepositionEfficiency.evaluate(
                        relative.norm() * 1000.0);
                const double speedPcf =
                    parameters.traverseSpeedPeakCorrectionFactor.evaluate(
                        pose.linearVelocity.norm() * 1000.0);
                const double height = parameters.peakCylinderHeightMeters
                    * std::exp(-0.5 * radial * radial
                        / (parameters.gaussianSigmaMeters
                            * parameters.gaussianSigmaMeters))
                    * angleRde * distanceRde * speedPcf;
                if(height <= 0.0 || !std::isfinite(height)) {
                    continue;
                }
                result.depositedCylinders.push_back({
                    hit.position,
                    -direction,
                    parameters.cylinderRadiusMeters,
                    height
                });
                ++result.statistics.candidatePairCount;
            }
            reportProgress(execution,
                static_cast<double>(poseIndex + 1)
                    / static_cast<double>(poses.size() - 1),
                "Wu 2020 reproduction");
        }
        result.statistics.elapsedMilliseconds = elapsedMilliseconds(started);
        return result;
    }
}
