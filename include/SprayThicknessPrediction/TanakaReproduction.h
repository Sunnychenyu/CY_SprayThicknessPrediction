#pragma once

#include <SprayThicknessPrediction/PublishedReproductionCommon.h>

#include <limits>

namespace spraythickness::published
{
    struct TanakaTargetPoint
    {
        Eigen::Vector3d position = Eigen::Vector3d::Zero();
        Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
    };

    struct TanakaParameters
    {
        double paintDischargeCubicMetersPerSecond{
            std::numeric_limits<double>::quiet_NaN() };
        double referenceDistanceMeters{
            std::numeric_limits<double>::quiet_NaN() };
        double referenceSigmaXMeters{
            std::numeric_limits<double>::quiet_NaN() };
        double referenceSigmaYMeters{
            std::numeric_limits<double>::quiet_NaN() };
        double distanceExponent{ 1.0 };
        double incidenceExponent{ 0.701 };

        void validate() const;
    };

    struct TanakaInputModel
    {
        std::vector<TanakaTargetPoint> targetPoints;
        std::vector<SprayPose> sprayPoses;
    };

    struct TanakaResult
    {
        std::vector<double> pointThicknessMeters;
        ReproductionStatistics statistics;
        bool canceled{ false };
    };

    class TanakaReproducer
    {
    public:
        static TanakaResult run(
            const TanakaInputModel& input,
            const TanakaParameters& parameters,
            const ReproductionExecution& execution = {});
    };
}
