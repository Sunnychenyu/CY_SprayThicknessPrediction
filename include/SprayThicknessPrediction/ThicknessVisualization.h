#pragma once

#include <Eigen/Core>

namespace spraythickness
{
    struct ThicknessColorRange
    {
        double minThickness{ 0.0 };
        double maxThickness{ 1.0 };
        Eigen::Vector3f minColor{ 0.05f, 0.20f, 1.0f };
        Eigen::Vector3f midColor{ 0.10f, 0.85f, 0.25f };
        Eigen::Vector3f maxColor{ 1.0f, 0.08f, 0.02f };
    };

    class ThicknessVisualization
    {
    public:
        // Compatibility entry point for existing spray SDK consumers. New callers
        // should use DemoThicknessPredictor; remove after the next SDK deprecation window.
        static double evaluateDemoBurnerThickness(
            const Eigen::Vector3d& position,
            const Eigen::Vector3d& normal,
            const Eigen::Vector3d& boundsMin,
            const Eigen::Vector3d& boundsMax);

        // Kept for source compatibility without adding a Spray -> VisualizationSDK
        // dependency. New visualization code should use visualization::ScalarColorMap.
        static Eigen::Vector3f mapThicknessToColor(
            double thickness,
            const ThicknessColorRange& range = {});
    };
}
