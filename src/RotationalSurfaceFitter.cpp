#include <SprayThicknessPrediction/RotationalSurfaceFitter.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <queue>
#include <unordered_set>
#include <vector>

namespace spraythickness
{
    namespace
    {
        using EdgeKey = std::array<long long, 6>;

        long long quantize(double value)
        {
            return static_cast<long long>(std::llround(value * 1.0e8));
        }

        std::array<long long, 3> quantizedPoint(const Eigen::Vector3d& point)
        {
            return { quantize(point.x()), quantize(point.y()), quantize(point.z()) };
        }

        EdgeKey makeEdgeKey(
            const Eigen::Vector3d& first,
            const Eigen::Vector3d& second)
        {
            const std::array<long long, 3> a = quantizedPoint(first);
            const std::array<long long, 3> b = quantizedPoint(second);
            if(a < b) {
                return { a[0], a[1], a[2], b[0], b[1], b[2] };
            }
            return { b[0], b[1], b[2], a[0], a[1], a[2] };
        }

        struct TriangleInfo
        {
            std::uint32_t indices[3] = { 0, 0, 0 };
            Eigen::Vector3d center = Eigen::Vector3d::Zero();
            Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
            bool valid = false;
        };

        bool finiteVector(const Eigen::Vector3d& value)
        {
            return value.allFinite();
        }

        void setFailure(
            CylindricalSurfaceFitResult& result,
            const char* message)
        {
            result.failureReason = message;
        }
    }

    CylindricalSurfaceFitResult fitCylindricalSurface(
        const sprayworkpiece::WorkpieceModel& workpiece,
        std::size_t seedTriangleIndex,
        const CylindricalSurfaceFitOptions& options)
    {
        CylindricalSurfaceFitResult result;
        const std::size_t triangleCount = workpiece.triangleIndices.size() / 3;
        if(workpiece.samples.empty() || triangleCount == 0) {
            setFailure(result, "The workpiece has no indexed triangles.");
            return result;
        }
        if(seedTriangleIndex >= triangleCount) {
            setFailure(result, "The picked triangle index is outside the workpiece.");
            return result;
        }

        std::vector<TriangleInfo> triangles(triangleCount);
        std::map<EdgeKey, std::vector<std::size_t>> edgeOwners;
        for(std::size_t triangle = 0; triangle < triangleCount; ++triangle) {
            TriangleInfo& info = triangles[triangle];
            info.indices[0] = workpiece.triangleIndices[triangle * 3];
            info.indices[1] = workpiece.triangleIndices[triangle * 3 + 1];
            info.indices[2] = workpiece.triangleIndices[triangle * 3 + 2];
            if(info.indices[0] >= workpiece.samples.size() ||
                info.indices[1] >= workpiece.samples.size() ||
                info.indices[2] >= workpiece.samples.size()) {
                continue;
            }
            const Eigen::Vector3d& a = workpiece.samples[info.indices[0]].position;
            const Eigen::Vector3d& b = workpiece.samples[info.indices[1]].position;
            const Eigen::Vector3d& c = workpiece.samples[info.indices[2]].position;
            const Eigen::Vector3d cross = (b - a).cross(c - a);
            if(!finiteVector(cross) || cross.squaredNorm() <= 1.0e-20) {
                continue;
            }
            info.center = (a + b + c) / 3.0;
            info.normal = cross.normalized();
            info.valid = true;
            edgeOwners[makeEdgeKey(a, b)].push_back(triangle);
            edgeOwners[makeEdgeKey(b, c)].push_back(triangle);
            edgeOwners[makeEdgeKey(c, a)].push_back(triangle);
        }
        if(!triangles[seedTriangleIndex].valid) {
            setFailure(result, "The picked triangle is degenerate or invalid.");
            return result;
        }

        std::vector<std::vector<std::size_t>> neighbors(triangleCount);
        for(const auto& entry : edgeOwners) {
            const std::vector<std::size_t>& owners = entry.second;
            for(std::size_t i = 0; i < owners.size(); ++i) {
                for(std::size_t j = i + 1; j < owners.size(); ++j) {
                    neighbors[owners[i]].push_back(owners[j]);
                    neighbors[owners[j]].push_back(owners[i]);
                }
            }
        }

        const double maximumAngleRadians =
            std::max(0.0, options.maximumAdjacentNormalAngleDegrees)
            * 3.14159265358979323846 / 180.0;
        const double minimumNormalDot = std::cos(maximumAngleRadians);
        std::vector<bool> included(triangleCount, false);
        std::queue<std::size_t> pending;
        pending.push(seedTriangleIndex);
        included[seedTriangleIndex] = true;
        std::vector<std::size_t> region;
        while(!pending.empty() && region.size() < options.maximumTriangleCount) {
            const std::size_t triangle = pending.front();
            pending.pop();
            region.push_back(triangle);
            for(const std::size_t neighbor : neighbors[triangle]) {
                if(included[neighbor] || !triangles[neighbor].valid) {
                    continue;
                }
                if(std::abs(triangles[triangle].normal.dot(triangles[neighbor].normal)) <
                    minimumNormalDot) {
                    continue;
                }
                included[neighbor] = true;
                pending.push(neighbor);
            }
        }
        if(region.size() < 3) {
            setFailure(result, "The selected surface did not expand to enough adjacent faces.");
            return result;
        }

        Eigen::Matrix3d normalCovariance = Eigen::Matrix3d::Zero();
        Eigen::Vector3d meanCenter = Eigen::Vector3d::Zero();
        for(const std::size_t triangle : region) {
            normalCovariance += triangles[triangle].normal * triangles[triangle].normal.transpose();
            meanCenter += triangles[triangle].center;
        }
        meanCenter /= static_cast<double>(region.size());
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> normalSolver(normalCovariance);
        if(normalSolver.info() != Eigen::Success) {
            setFailure(result, "Failed to estimate the cylindrical axis direction.");
            return result;
        }
        Eigen::Vector3d axis = normalSolver.eigenvectors().col(0);
        if(!finiteVector(axis) || axis.squaredNorm() <= 1.0e-16) {
            setFailure(result, "The selected surface has no stable axis direction.");
            return result;
        }
        axis.normalize();

        Eigen::Matrix3d system = axis * axis.transpose();
        Eigen::Vector3d rightHandSide = system * meanCenter;
        for(const std::size_t triangle : region) {
            const Eigen::Vector3d& normal = triangles[triangle].normal;
            const Eigen::Matrix3d projection =
                Eigen::Matrix3d::Identity() - normal * normal.transpose();
            system += projection;
            rightHandSide += projection * triangles[triangle].center;
        }
        if(std::abs(system.determinant()) <= 1.0e-14) {
            setFailure(result, "The selected surface has no stable axis position.");
            return result;
        }
        const Eigen::Vector3d axisOrigin = system.ldlt().solve(rightHandSide);
        if(!finiteVector(axisOrigin)) {
            setFailure(result, "Failed to estimate the cylindrical axis position.");
            return result;
        }

        std::unordered_set<std::uint32_t> vertexSet;
        for(const std::size_t triangle : region) {
            vertexSet.insert(triangles[triangle].indices[0]);
            vertexSet.insert(triangles[triangle].indices[1]);
            vertexSet.insert(triangles[triangle].indices[2]);
        }
        if(vertexSet.size() < 6) {
            setFailure(result, "The selected surface has too few unique vertices.");
            return result;
        }

        double radiusSum = 0.0;
        for(const std::uint32_t index : vertexSet) {
            const Eigen::Vector3d offset = workpiece.samples[index].position - axisOrigin;
            radiusSum += (offset - axis * offset.dot(axis)).norm();
        }
        const double radius = radiusSum / static_cast<double>(vertexSet.size());
        if(!std::isfinite(radius) || radius <= 1.0e-9) {
            setFailure(result, "The selected surface has an invalid fitted radius.");
            return result;
        }

        double squaredError = 0.0;
        double maximumError = 0.0;
        for(const std::uint32_t index : vertexSet) {
            const Eigen::Vector3d offset = workpiece.samples[index].position - axisOrigin;
            const double actualRadius = (offset - axis * offset.dot(axis)).norm();
            const double error = std::abs(actualRadius - radius);
            squaredError += error * error;
            maximumError = std::max(maximumError, error);
        }

        result.axisOrigin = axisOrigin;
        result.axisDirection = axis;
        result.seedTriangleIndex = seedTriangleIndex;
        result.selectedTriangleIndices = region;
        result.selectedVertexIndices.assign(vertexSet.begin(), vertexSet.end());
        result.radius = radius;
        result.rmsRadialError = std::sqrt(squaredError / static_cast<double>(vertexSet.size()));
        result.maximumRadialError = maximumError;
        result.triangleCount = region.size();
        result.vertexCount = vertexSet.size();
        if(result.rmsRadialError > options.maximumRmsRadialErrorMeters ||
            result.maximumRadialError > options.maximumRadialErrorMeters) {
            setFailure(result, "The expanded surface does not satisfy the cylindrical fit tolerance.");
        }
        return result;
    }
}
