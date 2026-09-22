#include <SprayThicknessPrediction/ThicknessPrediction.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace spraythickness
{
    const char* thicknessModelId(ThicknessModelKind model)
    {
        switch(model) {
        case ThicknessModelKind::PaperGaussian:
            return "paper_gaussian";
        }
        return "unknown";
    }

    bool ThicknessField::empty() const
    {
        return results.empty();
    }

    void ThicknessField::resizeFromWorkpiece(const sprayworkpiece::WorkpieceModel& workpiece)
    {
        results.clear();
        results.reserve(workpiece.samples.size());

        for (size_t i = 0; i < workpiece.samples.size(); ++i)
        {
            ThicknessSampleResult result;
            result.sampleIndex = i;
            result.targetThickness = workpiece.samples[i].targetThickness;
            result.error = result.thickness - result.targetThickness;
            results.push_back(result);
        }
    }

    void ThicknessField::updateErrors()
    {
        for (auto& result : results)
            result.error = result.thickness - result.targetThickness;
    }

    ThicknessMetrics ThicknessMetricsCalculator::calculate(
        const ThicknessField& field,
        const ThicknessPredictionOptions& options)
    {
        ThicknessMetrics metrics;
        if (field.results.empty())
            return metrics;

        metrics.minThickness = std::numeric_limits<double>::max();
        metrics.maxThickness = std::numeric_limits<double>::lowest();

        double sumThickness = 0.0;
        double sumError = 0.0;
        size_t coveredCount = 0;
        size_t underCount = 0;
        size_t overCount = 0;

        for (const auto& result : field.results)
        {
            metrics.minThickness = std::min(metrics.minThickness, result.thickness);
            metrics.maxThickness = std::max(metrics.maxThickness, result.thickness);
            sumThickness += result.thickness;
            sumError += result.error;
            metrics.maxAbsError = std::max(metrics.maxAbsError, std::abs(result.error));

            const double lower = result.targetThickness - options.coverageTolerance;
            const double upper = result.targetThickness + options.overCoatTolerance;

            if (result.thickness >= lower && result.thickness <= upper)
                ++coveredCount;
            if (result.thickness < lower)
                ++underCount;
            if (result.thickness > upper)
                ++overCount;
        }

        const double count = static_cast<double>(field.results.size());
        metrics.averageThickness = sumThickness / count;
        metrics.meanError = sumError / count;
        metrics.coverageRatio = static_cast<double>(coveredCount) / count;
        metrics.underCoatedRatio = static_cast<double>(underCount) / count;
        metrics.overCoatedRatio = static_cast<double>(overCount) / count;
        return metrics;
    }

    ThicknessPredictionResult SprayThicknessPredictor::predict(
        const sprayworkpiece::WorkpieceModel& workpiece,
        const spraycore::SprayTool& tool,
        const spraycore::SprayProcess& process,
        const spraytrajectory::SprayTrajectory& trajectory,
        const ThicknessPredictionOptions& options)
    {
        ThicknessPredictionResult result;
        result.field.resizeFromWorkpiece(workpiece);

        if (workpiece.samples.empty())
        {
            result.warnings.push_back("Workpiece has no surface samples.");
            result.metrics = ThicknessMetricsCalculator::calculate(result.field, options);
            return result;
        }

        if (trajectory.empty())
        {
            result.warnings.push_back("Spray trajectory is empty.");
            result.metrics = ThicknessMetricsCalculator::calculate(result.field, options);
            return result;
        }

        if (options.trajectorySamplingMode == TrajectorySamplingMode::ResampleByTimeStep
            && options.timeStep <= 0.0)
        {
            result.warnings.push_back("Thickness prediction timeStep must be positive.");
            result.metrics = ThicknessMetricsCalculator::calculate(result.field, options);
            return result;
        }

        const auto samples = options.trajectorySamplingMode == TrajectorySamplingMode::OriginalPoints
            ? spraytrajectory::SprayTrajectorySampler::originalSamples(trajectory)
            : spraytrajectory::SprayTrajectorySampler::sample(trajectory, options.timeStep);

        if (samples.empty())
        {
            result.warnings.push_back("Spray trajectory sampling produced no samples.");
            result.metrics = ThicknessMetricsCalculator::calculate(result.field, options);
            return result;
        }

        bool sawDifferentProcessId = false;

        for (size_t sampleIndex = 0; sampleIndex < samples.size(); ++sampleIndex)
        {
            const auto& trajectorySample = samples[sampleIndex];
            if (!trajectorySample.sprayEnabled)
                continue;

            if (!trajectorySample.processId.empty() && !process.id.empty()
                && trajectorySample.processId != process.id)
            {
                sawDifferentProcessId = true;
                continue;
            }

            const double deltaTime = sampleIndex + 1 < samples.size()
                ? std::max(0.0, samples[sampleIndex + 1].time - trajectorySample.time)
                : 0.0;
            if(deltaTime <= 0.0)
                continue;

            const Eigen::Isometry3d toolPose = trajectorySample.tcpPose * tool.T_link_tool;

            for (size_t i = 0; i < workpiece.samples.size(); ++i)
            {
                const auto& surfaceSample = workpiece.samples[i];
                if (!surfaceSample.valid)
                    continue;

                const double depositRate = spraycore::SprayKernel::evaluateDepositRate(
                    tool,
                    process,
                    toolPose,
                    surfaceSample.position,
                    surfaceSample.normal);

                const double areaWeight = options.useSampleAreaWeight
                    ? std::max(0.0, surfaceSample.areaWeight)
                    : 1.0;

                result.field.results[i].thickness += depositRate * deltaTime * areaWeight;
            }
        }

        if (sawDifferentProcessId)
            result.warnings.push_back("Some spray trajectory samples used a different processId and were skipped.");

        result.field.updateErrors();
        result.metrics = ThicknessMetricsCalculator::calculate(result.field, options);
        return result;
    }
}
