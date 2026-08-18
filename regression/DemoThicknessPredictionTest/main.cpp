#include <SprayThicknessPrediction/DemoThicknessPrediction.h>
#include <SprayThicknessPrediction/ThicknessVisualization.h>

#include <cmath>
#include <iostream>

namespace
{
    bool nearlyEqual(double lhs, double rhs, double tolerance = 1.0e-12)
    {
        return std::abs(lhs - rhs) <= tolerance;
    }
}

int main()
{
    sprayworkpiece::WorkpieceModel workpiece;
    for(double z : { 0.0, 0.5, 1.0 }) {
        sprayworkpiece::SurfaceSample sample;
        sample.position = Eigen::Vector3d(0.0, 0.0, z);
        workpiece.addSample(sample);
    }

    spraythickness::DemoThicknessPredictionOptions options;
    options.transverseWaveRatio = 0.0;
    const spraythickness::ThicknessPredictionResult first =
        spraythickness::DemoThicknessPredictor::predict(workpiece, options);
    const spraythickness::ThicknessPredictionResult second =
        spraythickness::DemoThicknessPredictor::predict(workpiece, options);
    if(first.field.results.size() != 3 || second.field.results.size() != 3) {
        std::cerr << "Unexpected result size.\n";
        return 1;
    }
    if(!nearlyEqual(first.field.results.front().thickness, options.minThicknessMeters) ||
        !nearlyEqual(first.field.results.back().thickness, options.maxThicknessMeters) ||
        !(first.field.results[0].thickness < first.field.results[1].thickness &&
            first.field.results[1].thickness < first.field.results[2].thickness)) {
        std::cerr << "Directional gradient contract failed.\n";
        return 2;
    }
    for(std::size_t i = 0; i < first.field.results.size(); ++i) {
        if(!std::isfinite(first.field.results[i].thickness) ||
            !nearlyEqual(first.field.results[i].thickness, second.field.results[i].thickness)) {
            std::cerr << "Prediction is not finite and deterministic.\n";
            return 3;
        }
    }

    sprayworkpiece::WorkpieceModel empty;
    const spraythickness::ThicknessPredictionResult emptyResult =
        spraythickness::DemoThicknessPredictor::predict(empty);
    if(!emptyResult.field.empty() || emptyResult.warnings.empty()) {
        std::cerr << "Empty input diagnostics contract failed.\n";
        return 4;
    }

    sprayworkpiece::WorkpieceModel flat;
    flat.addSample(sprayworkpiece::SurfaceSample());
    flat.addSample(sprayworkpiece::SurfaceSample());
    options.gradientDirection = Eigen::Vector3d::Zero();
    const spraythickness::ThicknessPredictionResult flatResult =
        spraythickness::DemoThicknessPredictor::predict(flat, options);
    if(flatResult.warnings.size() < 2 ||
        !std::isfinite(flatResult.field.results.front().thickness)) {
        std::cerr << "Degenerate input diagnostics contract failed.\n";
        return 5;
    }

    const double compatibilityValue =
        spraythickness::ThicknessVisualization::evaluateDemoBurnerThickness(
            Eigen::Vector3d(0.5, 0.5, 0.5),
            Eigen::Vector3d::UnitZ(),
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d::Ones());
    if(!std::isfinite(compatibilityValue) || compatibilityValue < 0.0 ||
        compatibilityValue > 1.0) {
        std::cerr << "Compatibility prediction wrapper failed.\n";
        return 6;
    }
    return 0;
}
