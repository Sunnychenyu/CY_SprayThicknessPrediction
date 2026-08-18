#pragma once

#include <WorkpieceCore/WorkpieceModel.h>

#include <Eigen/Core>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace spraythickness
{
    struct CylindricalSurfaceFitOptions
    {
        double maximumAdjacentNormalAngleDegrees{ 20.0 };
        double maximumRmsRadialErrorMeters{ 0.002 };
        double maximumRadialErrorMeters{ 0.005 };
        std::size_t maximumTriangleCount{ 2000000 };
    };

    struct CylindricalSurfaceFitResult
    {
        std::size_t seedTriangleIndex{ 0 };
        std::vector<std::size_t> selectedTriangleIndices;
        std::vector<std::uint32_t> selectedVertexIndices;
        Eigen::Vector3d axisOrigin = Eigen::Vector3d::Zero();
        Eigen::Vector3d axisDirection = Eigen::Vector3d::UnitZ();
        double radius{ 0.0 };
        double rmsRadialError{ 0.0 };
        double maximumRadialError{ 0.0 };
        std::size_t triangleCount{ 0 };
        std::size_t vertexCount{ 0 };
        std::string failureReason;

        bool valid() const
        {
            return failureReason.empty() && triangleCount > 0 && radius > 0.0;
        }
    };

    CylindricalSurfaceFitResult fitCylindricalSurface(
        const sprayworkpiece::WorkpieceModel& workpiece,
        std::size_t seedTriangleIndex,
        const CylindricalSurfaceFitOptions& options = {});
}
