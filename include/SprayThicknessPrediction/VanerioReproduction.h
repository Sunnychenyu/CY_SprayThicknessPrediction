#pragma once

#include <SprayThicknessPrediction/PublishedReproductionCommon.h>

#include <limits>

namespace spraythickness::published
{
    struct VanerioInputModel
    {
        TriangleMesh initialStlSurface;
        std::vector<SprayPose> nozzleTrajectory;
    };

    struct VanerioParameters
    {
        double growthRateCoefficientMetersPerSecond{
            std::numeric_limits<double>::quiet_NaN() };
        double jetRadiusMeters{ std::numeric_limits<double>::quiet_NaN() };
        double jetShapeCoefficientK2{
            std::numeric_limits<double>::quiet_NaN() };
        double maximumMeshEdgeMeters{
            std::numeric_limits<double>::quiet_NaN() };
        double shadowGridStepMeters{
            std::numeric_limits<double>::quiet_NaN() };
        TabulatedCurve depositionEfficiencyByTangentAngle;
        TabulatedCurve depositionEfficiencyByDistance;
        TabulatedCurve profileStretchByDistance;

        void validate() const;
    };

    struct VanerioResult
    {
        TriangleMesh evolvedStlSurface;
        ReproductionStatistics statistics;
        std::vector<std::string> implementationNotes;
        bool canceled{ false };
    };

    class VanerioReproducer
    {
    public:
        static VanerioResult run(
            const VanerioInputModel& input,
            const VanerioParameters& parameters,
            const ReproductionExecution& execution = {});
    };
}
