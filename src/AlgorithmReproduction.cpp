#include <SprayThicknessPrediction/AlgorithmReproduction.h>

#include <stdexcept>
#include <type_traits>

namespace spraythickness
{
    namespace
    {
        template<typename Input, typename Parameters, typename Reproducer>
        AlgorithmReproductionResult dispatch(
            const AlgorithmReproductionTask& task,
            const published::ReproductionExecution& execution)
        {
            const Input* input = std::get_if<Input>(&task.input);
            const Parameters* parameters =
                std::get_if<Parameters>(&task.parameters);
            if(input == nullptr || parameters == nullptr) {
                throw std::invalid_argument(
                    "The selected paper does not match the native input or parameter type.");
            }
            AlgorithmReproductionResult result;
            result.algorithm = task.algorithm;
            result.nativeResult = Reproducer::run(
                *input, *parameters, execution);
            return result;
        }
    }

    const char* reproductionAlgorithmId(ReproductionAlgorithmKind algorithm)
    {
        switch(algorithm) {
        case ReproductionAlgorithmKind::CurrentMethod:
            return "current_gpu";
        case ReproductionAlgorithmKind::Tanaka2024:
            return "tanaka_2024";
        case ReproductionAlgorithmKind::Tzinava2020:
            return "tzinava_2020";
        case ReproductionAlgorithmKind::Wu2020:
            return "wu_2020";
        case ReproductionAlgorithmKind::Fuke2005:
            return "fuke_2005";
        case ReproductionAlgorithmKind::Vanerio2021:
            return "vanerio_2021";
        case ReproductionAlgorithmKind::DynamicSurface2026:
            return "dynamic_surface_2026";
        }
        return "unknown";
    }

    const char* reproductionAlgorithmName(ReproductionAlgorithmKind algorithm)
    {
        switch(algorithm) {
        case ReproductionAlgorithmKind::CurrentMethod:
            return "Current BVH-GPU method";
        case ReproductionAlgorithmKind::Tanaka2024:
            return "Tanaka et al. (2024)";
        case ReproductionAlgorithmKind::Tzinava2020:
            return "Tzinava et al. (2020)";
        case ReproductionAlgorithmKind::Wu2020:
            return "Wu et al. (2020)";
        case ReproductionAlgorithmKind::Fuke2005:
            return "Fuke et al. (2005)";
        case ReproductionAlgorithmKind::Vanerio2021:
            return "Vanerio et al. (2021)";
        case ReproductionAlgorithmKind::DynamicSurface2026:
            return "Dynamic surface evolution (2026)";
        }
        return "Unknown reproduction";
    }

    const published::ReproductionStatistics& reproductionStatistics(
        const AlgorithmReproductionResult& result)
    {
        return std::visit([](const auto& native)
                -> const published::ReproductionStatistics& {
            return native.statistics;
        }, result.nativeResult);
    }

    bool reproductionCanceled(const AlgorithmReproductionResult& result)
    {
        return std::visit([](const auto& native) {
            return native.canceled;
        }, result.nativeResult);
    }

    std::vector<std::string> reproductionImplementationNotes(
        const AlgorithmReproductionResult& result)
    {
        return std::visit([](const auto& native) {
            using Result = std::decay_t<decltype(native)>;
            if constexpr(std::is_same_v<Result, published::VanerioResult>
                || std::is_same_v<Result, published::DynamicSurfaceResult>) {
                return native.implementationNotes;
            }
            return std::vector<std::string>{};
        }, result.nativeResult);
    }

    AlgorithmReproductionResult AlgorithmReproducer::run(
        const AlgorithmReproductionTask& task,
        const published::ReproductionExecution& execution)
    {
        switch(task.algorithm) {
        case ReproductionAlgorithmKind::Tanaka2024:
            return dispatch<published::TanakaInputModel,
                published::TanakaParameters,
                published::TanakaReproducer>(task, execution);
        case ReproductionAlgorithmKind::Tzinava2020:
            return dispatch<published::TzinavaInputModel,
                published::TzinavaParameters,
                published::TzinavaReproducer>(task, execution);
        case ReproductionAlgorithmKind::Wu2020:
            return dispatch<published::WuInputModel,
                published::WuParameters,
                published::WuReproducer>(task, execution);
        case ReproductionAlgorithmKind::Fuke2005:
            return dispatch<published::FukeInputModel,
                published::FukeParameters,
                published::FukeReproducer>(task, execution);
        case ReproductionAlgorithmKind::Vanerio2021:
            return dispatch<published::VanerioInputModel,
                published::VanerioParameters,
                published::VanerioReproducer>(task, execution);
        case ReproductionAlgorithmKind::DynamicSurface2026:
            return dispatch<published::DynamicSurfaceInputModel,
                published::DynamicSurfaceParameters,
                published::DynamicSurfaceReproducer>(task, execution);
        case ReproductionAlgorithmKind::CurrentMethod:
            throw std::invalid_argument(
                "The current GPU method is executed by ThicknessPrediction.");
        }
        throw std::invalid_argument("Unknown reproduction algorithm.");
    }
}
