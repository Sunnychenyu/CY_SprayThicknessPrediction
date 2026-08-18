#include <SprayThicknessPrediction/ThicknessVisualization.h>

#include <SprayThicknessPrediction/DemoThicknessPrediction.h>

#include <algorithm>
#include <cmath>

namespace spraythickness
{
    namespace
    {
        double clamp01(double value)
        {
            return std::max(0.0, std::min(1.0, value));
        }

        double safeRange(double minValue, double maxValue)
        {
            const double range = maxValue - minValue;
            return std::abs(range) > 1.0e-9 ? range : 1.0;
        }

        Eigen::Vector3f lerp(const Eigen::Vector3f& a, const Eigen::Vector3f& b, float t)
        {
            return a + (b - a) * t;
        }
    }

    double ThicknessVisualization::evaluateDemoBurnerThickness(
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& normal,
        const Eigen::Vector3d& boundsMin,
        const Eigen::Vector3d& boundsMax)
    {
        sprayworkpiece::WorkpieceModel workpiece;
        for(const Eigen::Vector3d& samplePosition : { boundsMin, position, boundsMax }) {
            sprayworkpiece::SurfaceSample sample;
            sample.position = samplePosition;
            sample.normal = normal;
            workpiece.addSample(sample);
        }

        DemoThicknessPredictionOptions options;
        options.minThicknessMeters = 0.0;
        options.maxThicknessMeters = 1.0;
        const ThicknessPredictionResult result = DemoThicknessPredictor::predict(
            workpiece,
            options);
        return result.field.results.size() > 1
            ? clamp01(result.field.results[1].thickness)
            : 0.5;
    }

    Eigen::Vector3f ThicknessVisualization::mapThicknessToColor(
        double thickness,
        const ThicknessColorRange& range)
    {
        const double denom = safeRange(range.minThickness, range.maxThickness);
        const double t = clamp01((thickness - range.minThickness) / denom);

        if (t < 0.5)
        {
            return lerp(range.minColor, range.midColor, static_cast<float>(t * 2.0));
        }

        return lerp(range.midColor, range.maxColor, static_cast<float>((t - 0.5) * 2.0));
    }
}
