#include <SprayThicknessPrediction/DemoThicknessPrediction.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace spraythickness
{
    namespace
    {
        constexpr double kEpsilon = 1.0e-12;
        constexpr double kTwoPi = 6.28318530717958647692;

        double clamp01(double value)
        {
            return std::max(0.0, std::min(1.0, value));
        }

        Eigen::Vector3d perpendicularAxis(const Eigen::Vector3d& direction)
        {
            Eigen::Vector3d basis = Eigen::Vector3d::UnitX();
            if(std::abs(direction.dot(basis)) > 0.85) {
                basis = Eigen::Vector3d::UnitY();
            }
            return direction.cross(basis).normalized();
        }
    }

    ThicknessPredictionResult DemoThicknessPredictor::predict(
        const sprayworkpiece::WorkpieceModel& workpiece,
        const DemoThicknessPredictionOptions& options)
    {
        ThicknessPredictionResult result;
        result.field.resizeFromWorkpiece(workpiece);
        if(workpiece.samples.empty()) {
            result.warnings.push_back("Workpiece has no surface samples.");
            result.metrics = ThicknessMetricsCalculator::calculate(result.field);
            return result;
        }

        Eigen::Vector3d direction = options.gradientDirection;
        if(!direction.allFinite() || direction.norm() <= kEpsilon) {
            direction = Eigen::Vector3d::UnitZ();
            result.warnings.push_back("Demo thickness gradient direction was invalid; +Z was used.");
        } else {
            direction.normalize();
        }
        const Eigen::Vector3d transverse = perpendicularAxis(direction);

        double minProjection = std::numeric_limits<double>::max();
        double maxProjection = std::numeric_limits<double>::lowest();
        double minTransverse = std::numeric_limits<double>::max();
        double maxTransverse = std::numeric_limits<double>::lowest();
        std::size_t validCount = 0;
        for(const sprayworkpiece::SurfaceSample& sample : workpiece.samples) {
            if(!sample.valid || !sample.position.allFinite()) {
                continue;
            }
            const double projection = sample.position.dot(direction);
            const double transverseProjection = sample.position.dot(transverse);
            minProjection = std::min(minProjection, projection);
            maxProjection = std::max(maxProjection, projection);
            minTransverse = std::min(minTransverse, transverseProjection);
            maxTransverse = std::max(maxTransverse, transverseProjection);
            ++validCount;
        }

        if(validCount == 0) {
            result.warnings.push_back("Workpiece has no valid finite surface samples.");
            result.metrics = ThicknessMetricsCalculator::calculate(result.field);
            return result;
        }

        const double lowThickness = std::min(options.minThicknessMeters, options.maxThicknessMeters);
        const double highThickness = std::max(options.minThicknessMeters, options.maxThicknessMeters);
        const double thicknessSpan = highThickness - lowThickness;
        const double projectionSpan = maxProjection - minProjection;
        const double transverseSpan = maxTransverse - minTransverse;
        const double waveRatio = std::max(0.0, std::min(0.45, options.transverseWaveRatio));

        if(projectionSpan <= kEpsilon) {
            result.warnings.push_back("Demo thickness projection range was degenerate; a uniform midpoint was used.");
        }

        for(std::size_t i = 0; i < workpiece.samples.size(); ++i) {
            const sprayworkpiece::SurfaceSample& sample = workpiece.samples[i];
            ThicknessSampleResult& sampleResult = result.field.results[i];
            if(!sample.valid || !sample.position.allFinite()) {
                sampleResult.thickness = lowThickness;
                continue;
            }

            double normalized = projectionSpan > kEpsilon
                ? (sample.position.dot(direction) - minProjection) / projectionSpan
                : 0.5;
            if(transverseSpan > kEpsilon && waveRatio > 0.0) {
                const double transverseNormalized =
                    (sample.position.dot(transverse) - minTransverse) / transverseSpan;
                normalized += waveRatio * std::sin(kTwoPi * transverseNormalized);
            }
            normalized = clamp01(normalized);
            sampleResult.thickness = lowThickness + thicknessSpan * normalized;
        }

        result.field.updateErrors();
        result.metrics = ThicknessMetricsCalculator::calculate(result.field);
        return result;
    }
}
