#pragma once

#include <SprayThicknessPrediction/PublishedReproductionCommon.h>

#include <cstdint>
#include <limits>

namespace spraythickness::published
{
    struct DynamicSurfaceInputModel
    {
        TriangleMesh initialStlSurface;
        std::vector<SprayPose> nozzleTrajectory;
    };

    struct DynamicSurfaceParameters
    {
        double equivalentParticleDiameterMeters{
            std::numeric_limits<double>::quiet_NaN() };
        double depositedCylinderRadiusMeters{
            std::numeric_limits<double>::quiet_NaN() };
        double equivalentParticleRatePerSecond{
            std::numeric_limits<double>::quiet_NaN() };
        double angularMeanRadians{
            std::numeric_limits<double>::quiet_NaN() };
        double angularSigmaRadians{
            std::numeric_limits<double>::quiet_NaN() };
        double maximumEjectionAngleRadians{
            std::numeric_limits<double>::quiet_NaN() };
        std::size_t polarBinCount{ 0 };
        double depositionBatchDurationSeconds{
            std::numeric_limits<double>::quiet_NaN() };
        std::uint32_t randomSeed{ 0 };
        TabulatedCurve relativeBuildUpByInclinationRadians;

        double voxelLeafMeters{ std::numeric_limits<double>::quiet_NaN() };
        double uniformSamplingRadiusMeters{
            std::numeric_limits<double>::quiet_NaN() };
        double normalSearchRadiusMeters{
            std::numeric_limits<double>::quiet_NaN() };
        std::size_t normalMaximumNeighbors{ 0 };
        double boundaryGradientThresholdRadiansPerMeter{
            std::numeric_limits<double>::quiet_NaN() };
        double dbscanRadiusMeters{ std::numeric_limits<double>::quiet_NaN() };
        std::size_t dbscanMinimumPoints{ 0 };
        double alphaShapeRadiusMeters{
            std::numeric_limits<double>::quiet_NaN() };
        int poissonDepth{ 0 };
        std::size_t maximumHoleBoundaryEdges{ 0 };
        std::size_t smoothingIterations{ 0 };
        double smoothingRelaxation{
            std::numeric_limits<double>::quiet_NaN() };

        void validate() const;
    };

    struct DynamicDepositedCylinder
    {
        Eigen::Vector3d baseCenter = Eigen::Vector3d::Zero();
        Eigen::Vector3d growthDirection = Eigen::Vector3d::UnitZ();
        double radiusMeters{ 0.0 };
        double heightMeters{ 0.0 };
    };

    struct DynamicSurfaceResult
    {
        TriangleMesh evolvedStlSurface;
        std::vector<DynamicDepositedCylinder> depositedCylinders;
        ReproductionStatistics statistics;
        std::vector<std::string> implementationNotes;
        bool canceled{ false };
    };

    class DynamicSurfaceReproducer
    {
    public:
        static DynamicSurfaceResult run(
            const DynamicSurfaceInputModel& input,
            const DynamicSurfaceParameters& parameters,
            const ReproductionExecution& execution = {});
    };
}
