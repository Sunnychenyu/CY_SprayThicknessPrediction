#pragma once

#include <SprayThicknessPrediction/DynamicSurfaceReproduction.h>
#include <SprayThicknessPrediction/FukeReproduction.h>
#include <SprayThicknessPrediction/TanakaReproduction.h>
#include <SprayThicknessPrediction/TzinavaReproduction.h>
#include <SprayThicknessPrediction/VanerioReproduction.h>
#include <SprayThicknessPrediction/WuReproduction.h>
#include <SprayThicknessPrediction/ThicknessPrediction.h>

#include <variant>

namespace spraythickness
{
    enum class ReproductionAlgorithmKind
    {
        CurrentMethod,
        Tanaka2024,
        Tzinava2020,
        Wu2020,
        Fuke2005,
        Vanerio2021,
        DynamicSurface2026
    };

    const char* reproductionAlgorithmId(ReproductionAlgorithmKind algorithm);
    const char* reproductionAlgorithmName(ReproductionAlgorithmKind algorithm);

    using PublishedReproductionInput = std::variant<
        published::TanakaInputModel,
        published::TzinavaInputModel,
        published::WuInputModel,
        published::FukeInputModel,
        published::VanerioInputModel,
        published::DynamicSurfaceInputModel>;

    using PublishedReproductionParameters = std::variant<
        published::TanakaParameters,
        published::TzinavaParameters,
        published::WuParameters,
        published::FukeParameters,
        published::VanerioParameters,
        published::DynamicSurfaceParameters>;

    struct CurrentMethodReproductionResult
    {
        ThicknessPredictionResult prediction;
        published::ReproductionStatistics statistics;
        ThicknessModelKind thicknessModel{ ThicknessModelKind::PaperGaussian };
        bool bvhOcclusion{ false };
        bool historyCorrection{ false };
        bool canceled{ false };
    };

    using PublishedReproductionNativeResult = std::variant<
        CurrentMethodReproductionResult,
        published::TanakaResult,
        published::TzinavaResult,
        published::WuResult,
        published::FukeResult,
        published::VanerioResult,
        published::DynamicSurfaceResult>;

    struct AlgorithmReproductionTask
    {
        ReproductionAlgorithmKind algorithm{
            ReproductionAlgorithmKind::Tanaka2024 };
        PublishedReproductionInput input;
        PublishedReproductionParameters parameters;
    };

    struct AlgorithmReproductionResult
    {
        ReproductionAlgorithmKind algorithm{
            ReproductionAlgorithmKind::Tanaka2024 };
        PublishedReproductionNativeResult nativeResult;
    };

    const published::ReproductionStatistics& reproductionStatistics(
        const AlgorithmReproductionResult& result);
    bool reproductionCanceled(const AlgorithmReproductionResult& result);
    std::vector<std::string> reproductionImplementationNotes(
        const AlgorithmReproductionResult& result);

    class AlgorithmReproducer
    {
    public:
        static AlgorithmReproductionResult run(
            const AlgorithmReproductionTask& task,
            const published::ReproductionExecution& execution = {});
    };
}
