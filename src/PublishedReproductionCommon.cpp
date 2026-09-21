#include <SprayThicknessPrediction/PublishedReproductionCommon.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace spraythickness::published
{
    namespace
    {
        constexpr double kEpsilon = 1.0e-12;

        bool rayTriangleIntersection(
            const Eigen::Vector3d& origin,
            const Eigen::Vector3d& direction,
            const Eigen::Vector3d& first,
            const Eigen::Vector3d& second,
            const Eigen::Vector3d& third,
            double& distance)
        {
            const Eigen::Vector3d edge1 = second - first;
            const Eigen::Vector3d edge2 = third - first;
            const Eigen::Vector3d p = direction.cross(edge2);
            const double determinant = edge1.dot(p);
            if(std::abs(determinant) <= kEpsilon) {
                return false;
            }
            const double inverse = 1.0 / determinant;
            const Eigen::Vector3d t = origin - first;
            const double u = t.dot(p) * inverse;
            if(u < 0.0 || u > 1.0) {
                return false;
            }
            const Eigen::Vector3d q = t.cross(edge1);
            const double v = direction.dot(q) * inverse;
            if(v < 0.0 || u + v > 1.0) {
                return false;
            }
            distance = edge2.dot(q) * inverse;
            return distance > kEpsilon;
        }
    }

    bool TriangleMesh::empty() const
    {
        return vertices.empty() || faces.empty();
    }

    void TriangleMesh::validate() const
    {
        if(empty()) {
            throw std::invalid_argument("Triangle mesh is empty.");
        }
        for(const auto& face : faces) {
            for(const std::uint32_t index : face) {
                if(index >= vertices.size()) {
                    throw std::invalid_argument("Triangle mesh contains an invalid vertex index.");
                }
            }
        }
    }

    bool canceled(const ReproductionExecution& execution)
    {
        return execution.cancelRequested != nullptr
            && execution.cancelRequested->load();
    }

    void reportProgress(
        const ReproductionExecution& execution,
        double value,
        const std::string& message)
    {
        if(execution.progress) {
            execution.progress(std::clamp(value, 0.0, 1.0), message);
        }
    }

    Eigen::Vector3d faceNormal(
        const TriangleMesh& mesh,
        std::size_t faceIndex)
    {
        const auto& face = mesh.faces.at(faceIndex);
        const Eigen::Vector3d normal =
            (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
        if(normal.squaredNorm() <= kEpsilon) {
            return Eigen::Vector3d::Zero();
        }
        return normal.normalized();
    }

    Eigen::Vector3d faceCentroid(
        const TriangleMesh& mesh,
        std::size_t faceIndex)
    {
        const auto& face = mesh.faces.at(faceIndex);
        return (mesh.vertices[face[0]] + mesh.vertices[face[1]]
            + mesh.vertices[face[2]]) / 3.0;
    }

    bool nearestRayHit(
        const TriangleMesh& mesh,
        const Eigen::Vector3d& origin,
        const Eigen::Vector3d& direction,
        RayHit& hit,
        std::size_t ignoredFace)
    {
        if(direction.squaredNorm() <= kEpsilon) {
            return false;
        }
        const Eigen::Vector3d unitDirection = direction.normalized();
        double nearestDistance = std::numeric_limits<double>::infinity();
        std::size_t nearestFace = 0;
        for(std::size_t faceIndex = 0; faceIndex < mesh.faces.size(); ++faceIndex) {
            if(faceIndex == ignoredFace) {
                continue;
            }
            const auto& face = mesh.faces[faceIndex];
            double distance = 0.0;
            if(rayTriangleIntersection(
                    origin,
                    unitDirection,
                    mesh.vertices[face[0]],
                    mesh.vertices[face[1]],
                    mesh.vertices[face[2]],
                    distance)
                && distance < nearestDistance) {
                nearestDistance = distance;
                nearestFace = faceIndex;
            }
        }
        if(!std::isfinite(nearestDistance)) {
            return false;
        }
        hit.distance = nearestDistance;
        hit.faceIndex = nearestFace;
        hit.position = origin + nearestDistance * unitDirection;
        hit.normal = faceNormal(mesh, nearestFace);
        return true;
    }

    void TabulatedCurve::validate(const std::string& name) const
    {
        if(arguments.size() < 2 || arguments.size() != values.size()
            || !std::is_sorted(arguments.begin(), arguments.end())) {
            throw std::invalid_argument(name + " must contain a sorted table.");
        }
        for(std::size_t index = 0; index < arguments.size(); ++index) {
            if(!std::isfinite(arguments[index]) || !std::isfinite(values[index])) {
                throw std::invalid_argument(name + " contains a non-finite value.");
            }
            if(index > 0 && arguments[index] <= arguments[index - 1]) {
                throw std::invalid_argument(name + " arguments must be unique.");
            }
        }
    }

    double TabulatedCurve::interpolate(
        double argument,
        const std::string& name) const
    {
        validate(name);
        if(argument < arguments.front() || argument > arguments.back()) {
            throw std::out_of_range(argument < arguments.front()
                ? name + " argument is below its calibrated range."
                : name + " argument is above its calibrated range.");
        }
        const auto upper = std::upper_bound(
            arguments.begin(), arguments.end(), argument);
        if(upper == arguments.end()) {
            return values.back();
        }
        const std::size_t high = static_cast<std::size_t>(
            upper - arguments.begin());
        const std::size_t low = high - 1;
        const double ratio = (argument - arguments[low])
            / (arguments[high] - arguments[low]);
        return (1.0 - ratio) * values[low] + ratio * values[high];
    }

    void Polynomial::validate(
        std::size_t expectedCoefficientCount,
        const std::string& name) const
    {
        if(coefficients.size() != expectedCoefficientCount) {
            throw std::invalid_argument(name + " has the wrong polynomial degree.");
        }
        for(const double coefficient : coefficients) {
            if(!std::isfinite(coefficient)) {
                throw std::invalid_argument(name + " contains a non-finite coefficient.");
            }
        }
    }

    double Polynomial::evaluate(double argument) const
    {
        double value = 0.0;
        for(auto coefficient = coefficients.rbegin();
            coefficient != coefficients.rend(); ++coefficient) {
            value = value * argument + *coefficient;
        }
        return value;
    }
}
