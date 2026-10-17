#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace spraythickness
{
    enum class MeshSamplingStatus
    {
        Insufficient,
        Minimum,
        Recommended,
        Missed,
        BackFacing
    };

    struct MeshSamplingPointAssessment
    {
        // Index in the applied sampler output, including interpolated samples.
        std::size_t sampleIndex{ 0 };
        double timeSeconds{ 0.0 };
        double durationSeconds{ 0.0 };
        std::size_t visibleVertexCount{ 0 };
        MeshSamplingStatus status{ MeshSamplingStatus::Missed };
        bool centerSampled{ false };
        bool bothDirectionsSampled{ false };
    };

    struct MeshSamplingAssessment
    {
        bool applicable{ false };
        std::string unavailableReason;
        std::size_t minimumVertexCount{ 20 };
        std::size_t recommendedVertexCount{ 50 };
        double corePatternRatio{ 0.05 };
        std::size_t totalPointCount{ 0 };
        std::size_t airMovePointCount{ 0 };
        std::size_t ignoredPointCount{ 0 };
        std::size_t assessedPointCount{ 0 };
        std::size_t missedPointCount{ 0 };
        std::size_t backFacingPointCount{ 0 };
        std::size_t insufficientPointCount{ 0 };
        std::size_t recommendedPointCount{ 0 };
        std::size_t limitedDistributionPointCount{ 0 };
        std::size_t centerGapPointCount{ 0 };
        std::size_t directionalGapPointCount{ 0 };
        std::size_t minimumVisibleVertexCount{ 0 };
        double averageVisibleVertexCount{ 0.0 };
        double assessedDurationSeconds{ 0.0 };
        double insufficientDurationSeconds{ 0.0 };
        double recommendedDurationSeconds{ 0.0 };
        std::size_t worstSampleIndex{ 0 };
        double worstTimeSeconds{ 0.0 };
        double elapsedMilliseconds{ 0.0 };
        std::vector<MeshSamplingPointAssessment> points;
    };
}
