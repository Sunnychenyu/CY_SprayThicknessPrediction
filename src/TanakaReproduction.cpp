#include <SprayThicknessPrediction/TanakaReproduction.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace spraythickness::published
{
    namespace
    {
        constexpr double kPi = 3.14159265358979323846;
        constexpr double kEpsilon = 1.0e-12;

        double elapsedMilliseconds(
            const std::chrono::steady_clock::time_point& started)
        {
            return std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - started).count();
        }

        double poseDuration(
            const std::vector<SprayPose>& poses,
            std::size_t index)
        {
            if(index + 1 >= poses.size()) {
                return 0.0;
            }
            return std::max(0.0,
                poses[index + 1].timeSeconds - poses[index].timeSeconds);
        }

        void sprayBasis(
            const Eigen::Vector3d& axis,
            Eigen::Vector3d& x,
            Eigen::Vector3d& y)
        {
            const Eigen::Vector3d unitAxis = axis.normalized();
            const Eigen::Vector3d reference = std::abs(unitAxis.z()) < 0.9
                ? Eigen::Vector3d::UnitZ()
                : Eigen::Vector3d::UnitY();
            x = reference.cross(unitAxis).normalized();
            y = unitAxis.cross(x).normalized();
        }
    }

    void TanakaParameters::validate() const
    {
        if(!std::isfinite(paintDischargeCubicMetersPerSecond)
            || paintDischargeCubicMetersPerSecond <= 0.0
            || !std::isfinite(referenceDistanceMeters)
            || referenceDistanceMeters <= 0.0
            || !std::isfinite(referenceSigmaXMeters)
            || referenceSigmaXMeters <= 0.0
            || !std::isfinite(referenceSigmaYMeters)
            || referenceSigmaYMeters <= 0.0
            || !std::isfinite(distanceExponent)
            || !std::isfinite(incidenceExponent)
            || incidenceExponent <= 0.0) {
            throw std::invalid_argument(
                "Tanaka reproduction requires Q, reference distance, both "
                "reference sigmas, alpha, and p.");
        }
    }

    TanakaResult TanakaReproducer::run(
        const TanakaInputModel& input,
        const TanakaParameters& parameters,
        const ReproductionExecution& execution)
    {
        parameters.validate();
        if(input.targetPoints.empty()) {
            throw std::invalid_argument("Tanaka target point set is empty.");
        }
        if(input.sprayPoses.size() < 2) {
            throw std::invalid_argument(
                "Tanaka reproduction requires at least two timed spray poses.");
        }

        const auto started = std::chrono::steady_clock::now();
        TanakaResult result;
        result.pointThicknessMeters.assign(input.targetPoints.size(), 0.0);
        result.statistics.trajectorySampleCount = input.sprayPoses.size();
        result.statistics.evaluatedElementCount = input.targetPoints.size();

        for(std::size_t poseIndex = 0;
            poseIndex + 1 < input.sprayPoses.size(); ++poseIndex) {
            if(canceled(execution)) {
                result.canceled = true;
                break;
            }
            const SprayPose& pose = input.sprayPoses[poseIndex];
            const double dt = poseDuration(input.sprayPoses, poseIndex);
            if(!pose.sprayEnabled || dt <= 0.0) {
                continue;
            }
            if(pose.axis.squaredNorm() <= kEpsilon) {
                throw std::invalid_argument("Tanaka spray axis is zero.");
            }
            const Eigen::Vector3d axis = pose.axis.normalized();
            Eigen::Vector3d basisX;
            Eigen::Vector3d basisY;
            sprayBasis(axis, basisX, basisY);

            for(std::size_t pointIndex = 0;
                pointIndex < input.targetPoints.size(); ++pointIndex) {
                const TanakaTargetPoint& point = input.targetPoints[pointIndex];
                const Eigen::Vector3d relative = point.position - pose.position;
                const double axialDistance = relative.dot(axis);
                if(axialDistance <= kEpsilon
                    || point.normal.squaredNorm() <= kEpsilon) {
                    continue;
                }
                const double surfaceCosine = std::clamp(
                    -point.normal.normalized().dot(relative.normalized()),
                    -1.0,
                    1.0);
                if(surfaceCosine <= 0.0) {
                    continue;
                }

                const double impactAngleDegrees =
                    std::asin(surfaceCosine) * 180.0 / kPi;
                const double incidenceFactor = std::pow(
                    impactAngleDegrees / 90.0,
                    parameters.incidenceExponent);
                const double distanceScale = std::pow(
                    axialDistance / parameters.referenceDistanceMeters,
                    parameters.distanceExponent);
                const double sigmaX =
                    distanceScale * parameters.referenceSigmaXMeters;
                const double sigmaY =
                    distanceScale * parameters.referenceSigmaYMeters;
                if(sigmaX <= kEpsilon || sigmaY <= kEpsilon) {
                    continue;
                }

                const double x = relative.dot(basisX);
                const double y = relative.dot(basisY);
                const double exponent = -0.5
                    * (x * x / (sigmaX * sigmaX)
                        + y * y / (sigmaY * sigmaY));
                const double thicknessRate =
                    parameters.paintDischargeCubicMetersPerSecond
                    / (2.0 * kPi * sigmaX * sigmaY)
                    * std::exp(exponent) * incidenceFactor;
                result.pointThicknessMeters[pointIndex] += thicknessRate * dt;
                ++result.statistics.candidatePairCount;
            }
            reportProgress(execution,
                static_cast<double>(poseIndex + 1)
                    / static_cast<double>(input.sprayPoses.size() - 1),
                "Tanaka 2024 reproduction");
        }
        result.statistics.elapsedMilliseconds = elapsedMilliseconds(started);
        return result;
    }
}
