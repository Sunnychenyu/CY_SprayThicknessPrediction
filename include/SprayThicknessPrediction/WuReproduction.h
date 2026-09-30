#pragma once

#include <SprayThicknessPrediction/PublishedReproductionCommon.h>

#include <limits>

namespace spraythickness::published
{
    struct WuDepositedCylinder
    {
        Eigen::Vector3d baseCenter = Eigen::Vector3d::Zero();
        Eigen::Vector3d growthDirection = Eigen::Vector3d::UnitZ();
        double radiusMeters{ 0.02e-3 };
        double heightMeters{ 0.0 };
        std::size_t substrateFaceIndex{ static_cast<std::size_t>(-1) };
    };

    struct WuInputModel
    {
        TriangleMesh substrate;
        std::vector<SprayPose> nozzleTrajectory;
        bool generatedPlateStack{ false };
    };

    struct WuParameters
    {
        double peakCylinderHeightMeters{
            std::numeric_limits<double>::quiet_NaN() };
        // Time represented by peakCylinderHeightMeters in the calibration.
        double referencePoseDurationSeconds{ 1.0 };
        double gaussianSigmaMeters{
            std::numeric_limits<double>::quiet_NaN() };
        double maximumDeflectionRadians{
            std::numeric_limits<double>::quiet_NaN() };
        double rayAngularStepRadians{
            std::numeric_limits<double>::quiet_NaN() };
        double cylinderRadiusMeters{ 0.02e-3 };
        Polynomial sprayAngleRelativeDepositionEfficiency;
        Polynomial sprayDistanceRelativeDepositionEfficiency;
        Polynomial traverseSpeedPeakCorrectionFactor;

        void validate() const;
    };

    struct WuResult
    {
        TriangleMesh substrate;
        std::vector<WuDepositedCylinder> depositedCylinders;
        ReproductionStatistics statistics;
        bool canceled{ false };
    };

    class WuReproducer
    {
    public:
        static WuResult run(
            const WuInputModel& input,
            const WuParameters& parameters,
            const ReproductionExecution& execution = {});
    };
}
