#pragma once

#include <SprayCore/SprayModel.h>
#include <SprayTrajectoryCore/SprayTrajectory.h>
#include <WorkpieceCore/WorkpieceModel.h>

#include <cstddef>
#include <string>
#include <vector>

namespace spraythickness
{
    struct ThicknessSampleResult
    {
        size_t sampleIndex{ 0 };
        double thickness{ 0.0 };
        double targetThickness{ 0.0 };
        double error{ 0.0 };
    };

    class ThicknessField
    {
    public:
        std::vector<ThicknessSampleResult> results;

        bool empty() const;
        void resizeFromWorkpiece(const sprayworkpiece::WorkpieceModel& workpiece);
        void updateErrors();
    };

    struct ThicknessMetrics
    {
        double minThickness{ 0.0 };
        double maxThickness{ 0.0 };
        double averageThickness{ 0.0 };
        double meanError{ 0.0 };
        double maxAbsError{ 0.0 };
        double coverageRatio{ 0.0 };
        double underCoatedRatio{ 0.0 };
        double overCoatedRatio{ 0.0 };
    };

    struct ThicknessPredictionOptions
    {
        double timeStep{ 0.02 };
        double coverageTolerance{ 0.01 };
        double overCoatTolerance{ 0.01 };
        bool useSampleAreaWeight{ true };
    };

    struct ThicknessPredictionResult
    {
        ThicknessField field;
        ThicknessMetrics metrics;
        std::vector<std::string> warnings;
    };

    class ThicknessMetricsCalculator
    {
    public:
        static ThicknessMetrics calculate(
            const ThicknessField& field,
            const ThicknessPredictionOptions& options = {});
    };

    class SprayThicknessPredictor
    {
    public:
        static ThicknessPredictionResult predict(
            const sprayworkpiece::WorkpieceModel& workpiece,
            const spraycore::SprayTool& tool,
            const spraycore::SprayProcess& process,
            const spraytrajectory::SprayTrajectory& trajectory,
            const ThicknessPredictionOptions& options = {});
    };
}
