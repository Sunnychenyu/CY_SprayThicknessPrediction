#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace spraythickness::published
{
    struct TriangleMesh
    {
        std::vector<Eigen::Vector3d> vertices;
        std::vector<std::array<std::uint32_t, 3>> faces;

        bool empty() const;
        void validate() const;
    };

    struct SprayPose
    {
        double timeSeconds{ 0.0 };
        Eigen::Vector3d position = Eigen::Vector3d::Zero();
        Eigen::Vector3d axis = -Eigen::Vector3d::UnitZ();
        Eigen::Vector3d linearVelocity = Eigen::Vector3d::Zero();
        bool sprayEnabled{ true };
        int passIndex{ 0 };
    };

    struct ReproductionExecution
    {
        const std::atomic_bool* cancelRequested{ nullptr };
        std::function<void(double, const std::string&)> progress;
    };

    struct ReproductionStatistics
    {
        double elapsedMilliseconds{ 0.0 };
        std::size_t trajectorySampleCount{ 0 };
        std::size_t evaluatedElementCount{ 0 };
        std::size_t candidatePairCount{ 0 };
        std::size_t visibilityQueryCount{ 0 };
        std::size_t hiddenElementCount{ 0 };
    };

    struct RayHit
    {
        double distance{ 0.0 };
        std::size_t faceIndex{ 0 };
        Eigen::Vector3d position = Eigen::Vector3d::Zero();
        Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
    };

    struct TabulatedCurve
    {
        std::vector<double> arguments;
        std::vector<double> values;

        void validate(const std::string& name) const;
        double interpolate(double argument, const std::string& name) const;
    };

    struct Polynomial
    {
        std::vector<double> coefficients;

        void validate(std::size_t expectedCoefficientCount,
            const std::string& name) const;
        double evaluate(double argument) const;
    };

    bool canceled(const ReproductionExecution& execution);
    void reportProgress(
        const ReproductionExecution& execution,
        double value,
        const std::string& message);

    Eigen::Vector3d faceNormal(
        const TriangleMesh& mesh,
        std::size_t faceIndex);
    Eigen::Vector3d faceCentroid(
        const TriangleMesh& mesh,
        std::size_t faceIndex);
    bool nearestRayHit(
        const TriangleMesh& mesh,
        const Eigen::Vector3d& origin,
        const Eigen::Vector3d& direction,
        RayHit& hit,
        std::size_t ignoredFace = static_cast<std::size_t>(-1));
}
