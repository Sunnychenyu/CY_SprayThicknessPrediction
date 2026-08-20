#pragma once

#include <SprayCore/SprayModel.h>
#include <SprayTrajectoryCore/SprayTrajectory.h>
#include <WorkpieceCore/WorkpieceModel.h>

#include <cstddef>
#include <cstdint>
#include <atomic>
#include <functional>
#include <string>
#include <vector>

namespace spraythickness
{
    enum class ThicknessModelKind
    {
        PaperGaussian
    };

    enum class TrajectorySamplingMode
    {
        OriginalPoints,
        ResampleByTimeStep
    };

    const char* thicknessModelId(ThicknessModelKind model);

    struct ThicknessSampleResult
    {
        size_t sampleIndex{ 0 };
        double thickness{ 0.0 };
        double targetThickness{ 0.0 };
        double error{ 0.0 };
    };

    class ThicknessField
    {
    public:
        std::vector<ThicknessSampleResult> results;

        bool empty() const;
        void resizeFromWorkpiece(const sprayworkpiece::WorkpieceModel& workpiece);
        void updateErrors();
    };

    struct ThicknessMetrics
    {
        double minThickness{ 0.0 };
        double maxThickness{ 0.0 };
        double averageThickness{ 0.0 };
        double meanError{ 0.0 };
        double maxAbsError{ 0.0 };
        double coverageRatio{ 0.0 };
        double underCoatedRatio{ 0.0 };
        double overCoatedRatio{ 0.0 };
    };

    struct ThicknessPredictionTiming
    {
        bool valid{ false };
        bool spatialFiltering{ false };
        bool spatialGridCacheHit{ false };
        bool axisymmetricMappingCacheHit{ false };
        bool manualSpatialGridCellSize{ false };
        double periodicMappingMilliseconds{ 0.0 };
        double axisymmetricMappingMilliseconds{ 0.0 };
        double spatialGridMilliseconds{ 0.0 };
        double bvhMilliseconds{ 0.0 };
        double uploadMilliseconds{ 0.0 };
        double dispatchMilliseconds{ 0.0 };
        double pureGpuMilliseconds{ 0.0 };
        double readbackMilliseconds{ 0.0 };
        double backendTotalMilliseconds{ 0.0 };
        double spatialGridCellSizeMeters{ 0.0 };
        std::size_t predictionVertexCount{ 0 };
        std::size_t spatialInputVertexCount{ 0 };
        std::size_t spatialSkippedVertexCount{ 0 };
        std::size_t sprayPointCount{ 0 };
        std::size_t spatialGridCellCount{ 0 };
        std::size_t spatialGridCandidateCellPairs{ 0 };
        std::size_t spatialGridCandidateVertexPairs{ 0 };
        std::size_t axisymmetricBindingCount{ 0 };
        std::size_t axisymmetricActiveBindingCount{ 0 };
        std::size_t axisymmetricMappedZeroCount{ 0 };
        std::size_t axisymmetricMappedNonzeroCount{ 0 };
        std::size_t axisymmetricMappedInvalidCount{ 0 };
        std::size_t axisymmetricActiveSegmentCount{ 0 };
        std::size_t axisymmetricActiveProfilePathCount{ 0 };
        int spatialGridDimensionX{ 0 };
        int spatialGridDimensionY{ 0 };
        int spatialGridDimensionZ{ 0 };
        std::size_t adaptiveBatchSize{ 0 };
    };

    struct ThicknessPredictionOptions
    {
        TrajectorySamplingMode trajectorySamplingMode{ TrajectorySamplingMode::OriginalPoints };
        double timeStep{ 0.02 };
        double coverageTolerance{ 0.01 };
        double overCoatTolerance{ 0.01 };
        bool useSampleAreaWeight{ true };
    };

    struct ThicknessPredictionResult
    {
        ThicknessField field;
        ThicknessMetrics metrics;
        ThicknessPredictionTiming timing;
        std::vector<std::string> warnings;
    };

    struct PaperGaussianParameters
    {
        double amplitudeMillimeters{ 0.5953252911 };
        double rotationRadians{ -0.04494047503 };
        double phiOffsetRadians{ 0.0002884115775 };
        double psiOffsetRadians{ -0.002083376641 };
        double sigmaPhiRadians{ 0.02793467718 };
        double sigmaPsiRadians{ 0.03265327731 };
        double referenceDistanceMeters{ 0.120 };
        double referenceAngleDegrees{ 90.0 };
        double referenceExposureSeconds{ 4.0 };
    };

    struct HeatHistoryParameters
    {
        double correctionAmplitude{ 4.164041379 };
        double historyScaleSeconds{ 62.61289749 };
        double referenceHistorySeconds{ 4.0 };
        double coolingTimeSeconds{ 0.144639594 };
        double activityThresholdRatio{ 0.05 };
    };

    struct PeriodicLocalPredictionOptions
    {
        bool enabled{ false };
        Eigen::Vector3d axisOrigin = Eigen::Vector3d::Zero();
        Eigen::Vector3d axisDirection = Eigen::Vector3d::UnitZ();
        Eigen::Vector3d referenceDirection = Eigen::Vector3d::UnitX();
        std::size_t sectorCount{ 1 };
        // Retained for source compatibility; exact sector selection ignores this value.
        double angularHaloRadians{ 0.0 };
        double maximumMappingDistanceMeters{ 1.0e-3 };
        double maximumNormalAngleDegrees{ 20.0 };
        double contributionCutoffRatio{ 1.0e-10 };
        bool reduceTrajectory{ false };
        bool fallbackToFullPrediction{ true };
    };

    struct AxisymmetricProfileBinding
    {
        std::uint32_t firstSampleIndex{ 0 };
        std::uint32_t secondSampleIndex{ 0 };
        // The selected segment is retained for branch diagnostics and cache
        // validation. It does not change the interpolation data path.
        std::uint32_t segmentIndex{ 0 };
        float interpolation{ 0.0f };
        bool active{ false };
    };

    struct AxisymmetricProfileSampleSegment
    {
        std::uint32_t firstSampleIndex{ 0 };
        std::uint32_t secondSampleIndex{ 0 };
        std::uint32_t profilePathIndex{ 0 };
        Eigen::Vector2d firstSectionPosition = Eigen::Vector2d::Zero();
        Eigen::Vector2d secondSectionPosition = Eigen::Vector2d::Zero();
    };

    struct AxisymmetricProfilePredictionOptions
    {
        bool enabled{ false };
        std::vector<sprayworkpiece::SurfaceSample> predictionSamples;
        std::vector<AxisymmetricProfileSampleSegment> sampleSegments;
        std::vector<AxisymmetricProfileBinding> fullVertexBindings;
        Eigen::Vector3d axisOrigin = Eigen::Vector3d::Zero();
        Eigen::Vector3d axisDirection = Eigen::Vector3d::UnitZ();
        Eigen::Vector3d radialDirection = Eigen::Vector3d::UnitX();
        Eigen::Vector2d selectionMinimum = Eigen::Vector2d::Zero();
        Eigen::Vector2d selectionMaximum = Eigen::Vector2d::Zero();
    };

    struct SpatialInfluenceFilteringOptions
    {
        bool enabled{ false };
        // Gaussian pattern values below this ratio are treated as zero by
        // the experimental spatial-filtered path.
        double contributionCutoffRatio{ 1.0e-10 };
        bool overrideGridCellSize{ false };
        double gridCellSizeMeters{ 0.0 };
        // Skip GPU dispatch for vertices whose spatial-grid cell has no
        // candidate spray points. Their thickness is exactly zero under the
        // same cutoff used by spatial filtering.
        bool filterCandidateVertices{ false };
        bool fallbackToFullPrediction{ true };
    };

    struct AdvancedThicknessPredictionOptions
    {
        ThicknessPredictionOptions base;
        PaperGaussianParameters deposition;
        HeatHistoryParameters history;
        bool enableBvhOcclusion{ true };
        bool enableHistoryCorrection{ true };
        double shadowBiasMeters{ 1.0e-5 };
        // Zero lets the GPU backend select and calibrate the batch size.
        std::size_t trajectoryBatchSize{ 0 };
        PeriodicLocalPredictionOptions periodicLocal;
        AxisymmetricProfilePredictionOptions axisymmetricProfile;
        SpatialInfluenceFilteringOptions spatialFiltering;
    };

    struct ThicknessPredictionTask
    {
        ThicknessModelKind model{ ThicknessModelKind::PaperGaussian };
        sprayworkpiece::WorkpieceModel workpiece;
        spraycore::SprayTool tool;
        spraycore::SprayProcess process;
        spraytrajectory::SprayTrajectory trajectory;
        AdvancedThicknessPredictionOptions options;
    };

    struct ThicknessPredictionExecution
    {
        const std::atomic_bool* cancelRequested{ nullptr };
        std::function<void(double, const std::string&)> progress;
    };

    class IThicknessPredictionBackend
    {
    public:
        virtual ~IThicknessPredictionBackend() = default;

        virtual ThicknessPredictionResult predict(
            const ThicknessPredictionTask& task,
            const ThicknessPredictionExecution& execution = {}) = 0;
    };

    class ThicknessMetricsCalculator
    {
    public:
        static ThicknessMetrics calculate(
            const ThicknessField& field,
            const ThicknessPredictionOptions& options = {});
    };

    class SprayThicknessPredictor
    {
    public:
        static ThicknessPredictionResult predict(
            const sprayworkpiece::WorkpieceModel& workpiece,
            const spraycore::SprayTool& tool,
            const spraycore::SprayProcess& process,
            const spraytrajectory::SprayTrajectory& trajectory,
            const ThicknessPredictionOptions& options = {});
    };
}
