#pragma once

#include <SprayThicknessPrediction/PublishedReproductionCommon.h>

#include <Eigen/Geometry>

#include <limits>

namespace spraythickness::published
{
    struct FukePolygon
    {
        std::vector<std::uint32_t> vertexIndices;
    };

    struct FukeInputModel
    {
        std::vector<Eigen::Vector3d> vertices;
        std::vector<FukePolygon> polygons;
        std::vector<double> timesSeconds;
        std::vector<Eigen::Isometry3d> workpiecePoses;
        Eigen::Vector3d vaporSourcePosition = Eigen::Vector3d::Zero();
        Eigen::Vector3d vaporSourceNormal = Eigen::Vector3d::UnitZ();
    };

    struct FukeParameters
    {
        double referenceThicknessRateMetersPerSecond{
            std::numeric_limits<double>::quiet_NaN() };
        double referenceDistanceMeters{
            std::numeric_limits<double>::quiet_NaN() };
        int plumeExponent{ 2 };

        void validate() const;
    };

    struct FukeResult
    {
        std::vector<Eigen::Vector3d> polygonCentroids;
        std::vector<double> polygonThicknessMeters;
        ReproductionStatistics statistics;
        bool canceled{ false };
    };

    class FukeReproducer
    {
    public:
        static FukeResult run(
            const FukeInputModel& input,
            const FukeParameters& parameters,
            const ReproductionExecution& execution = {});
    };
}
