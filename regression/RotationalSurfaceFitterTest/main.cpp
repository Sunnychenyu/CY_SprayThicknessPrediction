#include <SprayThicknessPrediction/RotationalSurfaceFitter.h>

#include <cmath>
#include <iostream>

int main()
{
    constexpr std::size_t segmentCount = 32;
    constexpr double radius = 0.25;
    constexpr double height = 0.8;
    sprayworkpiece::WorkpieceModel workpiece;

    for(std::size_t segment = 0; segment < segmentCount; ++segment) {
        const double angle0 = 2.0 * 3.14159265358979323846
            * static_cast<double>(segment) / static_cast<double>(segmentCount);
        const double angle1 = 2.0 * 3.14159265358979323846
            * static_cast<double>(segment + 1) / static_cast<double>(segmentCount);
        const std::size_t base = workpiece.samples.size();
        for(const auto& point : {
            Eigen::Vector3d(radius * std::cos(angle0), radius * std::sin(angle0), 0.0),
            Eigen::Vector3d(radius * std::cos(angle1), radius * std::sin(angle1), 0.0),
            Eigen::Vector3d(radius * std::cos(angle1), radius * std::sin(angle1), height),
            Eigen::Vector3d(radius * std::cos(angle0), radius * std::sin(angle0), height) }) {
            sprayworkpiece::SurfaceSample sample;
            sample.position = point;
            workpiece.addSample(sample);
        }
        workpiece.triangleIndices.push_back(static_cast<std::uint32_t>(base));
        workpiece.triangleIndices.push_back(static_cast<std::uint32_t>(base + 1));
        workpiece.triangleIndices.push_back(static_cast<std::uint32_t>(base + 2));
        workpiece.triangleIndices.push_back(static_cast<std::uint32_t>(base));
        workpiece.triangleIndices.push_back(static_cast<std::uint32_t>(base + 2));
        workpiece.triangleIndices.push_back(static_cast<std::uint32_t>(base + 3));
    }

    const spraythickness::CylindricalSurfaceFitResult fit =
        spraythickness::fitCylindricalSurface(workpiece, 0);
    if(!fit.valid() || std::abs(fit.radius - radius) > 1.0e-6 ||
        std::abs(std::abs(fit.axisDirection.dot(Eigen::Vector3d::UnitZ())) - 1.0) > 1.0e-6 ||
        fit.rmsRadialError > 1.0e-6) {
        std::cerr << "Cylindrical surface fitting contract failed: "
            << fit.failureReason << "\n";
        return 1;
    }
    return 0;
}
