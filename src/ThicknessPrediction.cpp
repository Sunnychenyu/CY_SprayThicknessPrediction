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
        ThicknessMetricsAccumulator accumulator(options);
        for (const auto& result : field.results)
            accumulator.add(result);
        return accumulator.metrics();
    }

    ThicknessMetricsAccumulator::ThicknessMetricsAccumulator(const ThicknessPredictionOptions& options)
        : m_coverageTolerance(options.coverageTolerance)
        , m_overCoatTolerance(options.overCoatTolerance)
    {
        m_metrics.minThickness = std::numeric_limits<double>::max();
        m_metrics.maxThickness = std::numeric_limits<double>::lowest();
    }

    ThicknessMetrics ThicknessMetricsAccumulator::metrics() const
    {
        if(m_count == 0) return {};
        ThicknessMetrics metrics = m_metrics;
        const double count = static_cast<double>(m_count);
        metrics.averageThickness = m_sumThickness / count;
        metrics.meanError = m_sumError / count;
        metrics.coverageRatio = static_cast<double>(m_coveredCount) / count;
        metrics.underCoatedRatio = static_cast<double>(m_underCount) / count;
        metrics.overCoatedRatio = static_cast<double>(m_overCount) / count;
        return metrics;
    }

    void ThicknessMetricsAccumulator::merge(const ThicknessMetricsAccumulator& other)
    {
        if(other.m_count == 0) return;
        m_metrics.minThickness = std::min(m_metrics.minThickness, other.m_metrics.minThickness);
        m_metrics.maxThickness = std::max(m_metrics.maxThickness, other.m_metrics.maxThickness);
        m_metrics.maxAbsError = std::max(m_metrics.maxAbsError, other.m_metrics.maxAbsError);
        m_sumThickness += other.m_sumThickness;
        m_sumError += other.m_sumError;
        m_count += other.m_count;
        m_coveredCount += other.m_coveredCount;
        m_underCount += other.m_underCount;
        m_overCount += other.m_overCount;
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
