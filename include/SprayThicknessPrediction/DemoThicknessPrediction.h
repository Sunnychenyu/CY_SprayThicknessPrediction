#pragma once

#include <SprayThicknessPrediction/ThicknessPrediction.h>

#include <Eigen/Core>

namespace spraythickness
{
    struct DemoThicknessPredictionOptions
    {
        Eigen::Vector3d gradientDirection = Eigen::Vector3d::UnitZ();
        double minThicknessMeters{ 20.0e-6 };
        double maxThicknessMeters{ 120.0e-6 };
        double transverseWaveRatio{ 0.10 };
    };

    class DemoThicknessPredictor
    {
    public:
        static ThicknessPredictionResult predict(
            const sprayworkpiece::WorkpieceModel& workpiece,
            const DemoThicknessPredictionOptions& options = {});
    };
}
