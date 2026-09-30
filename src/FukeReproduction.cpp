#include <SprayThicknessPrediction/FukeReproduction.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace spraythickness::published
{
    namespace
    {
        constexpr double kEpsilon = 1.0e-12;

        double elapsedMilliseconds(
            const std::chrono::steady_clock::time_point& started)
        {
            return std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - started).count();
        }

        TriangleMesh triangulate(
            const FukeInputModel& input,
            const Eigen::Isometry3d& pose)
        {
            TriangleMesh output;
            output.vertices.reserve(input.vertices.size());
            for(const Eigen::Vector3d& vertex : input.vertices) {
                output.vertices.push_back(pose * vertex);
            }
            for(std::size_t polygonIndex = 0;
                polygonIndex < input.polygons.size(); ++polygonIndex) {
                const auto& indices = input.polygons[polygonIndex].vertexIndices;
                if(indices.size() != 3 && indices.size() != 4) {
                    throw std::invalid_argument(
                        "Fuke polygons must be triangular or quadrilateral.");
                }
                for(const std::uint32_t index : indices) {
                    if(index >= input.vertices.size()) {
                        throw std::invalid_argument(
                            "Fuke polygon contains an invalid vertex index.");
                    }
                }
                output.faces.push_back({ indices[0], indices[1], indices[2] });
                if(indices.size() == 4) {
                    output.faces.push_back({ indices[0], indices[2], indices[3] });
                }
            }
            return output;
        }

        Eigen::Vector3d polygonCentroid(
            const FukeInputModel& input,
            const FukePolygon& polygon,
            const Eigen::Isometry3d& pose)
        {
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            for(const std::uint32_t index : polygon.vertexIndices) {
                centroid += pose * input.vertices.at(index);
            }
            return centroid / static_cast<double>(polygon.vertexIndices.size());
        }

        Eigen::Vector3d polygonNormal(
            const FukeInputModel& input,
            const FukePolygon& polygon,
            const Eigen::Isometry3d& pose)
        {
            const Eigen::Vector3d first = pose * input.vertices.at(
                polygon.vertexIndices.at(0));
            const Eigen::Vector3d second = pose * input.vertices.at(
                polygon.vertexIndices.at(1));
            const Eigen::Vector3d third = pose * input.vertices.at(
                polygon.vertexIndices.at(2));
            const Eigen::Vector3d normal = (second - first).cross(third - first);
            return normal.squaredNorm() > 1.0e-30
                ? normal.normalized()
                : Eigen::Vector3d::Zero();
        }

        bool visible(
            const TriangleMesh& mesh,
            const Eigen::Vector3d& source,
            const Eigen::Vector3d& centroid)
        {
            const Eigen::Vector3d ray = centroid - source;
            const double targetDistance = ray.norm();
            if(targetDistance <= kEpsilon) {
                return false;
            }
            RayHit hit;
            return !nearestRayHit(mesh, source, ray, hit)
                || hit.distance >= targetDistance - 1.0e-9;
        }
    }

    void FukeParameters::validate() const
    {
        if(!std::isfinite(referenceThicknessRateMetersPerSecond)
            || referenceThicknessRateMetersPerSecond <= 0.0
            || !std::isfinite(referenceDistanceMeters)
            || referenceDistanceMeters <= 0.0
            || plumeExponent < 2) {
            throw std::invalid_argument(
                "Fuke reproduction requires a reference deposition rate, "
                "reference source distance, and integer n >= 2.");
        }
    }

    FukeResult FukeReproducer::run(
        const FukeInputModel& input,
        const FukeParameters& parameters,
        const ReproductionExecution& execution)
    {
        reportDiagnostic(execution, "Input validation",
            "polygons=" + std::to_string(input.polygons.size())
                + ", poses=" + std::to_string(input.timesSeconds.size()));
        parameters.validate();
        if(input.vertices.empty() || input.polygons.empty()) {
            throw std::invalid_argument("Fuke polygon mesh is empty.");
        }
        if(input.timesSeconds.size() < 2
            || input.timesSeconds.size() != input.workpiecePoses.size()) {
            throw std::invalid_argument(
                "Fuke reproduction requires matching timed workpiece poses.");
        }
        if(!input.sprayEnabled.empty()
            && input.sprayEnabled.size() != input.timesSeconds.size()) {
            throw std::invalid_argument(
                "Fuke spray flags must match timed workpiece poses.");
        }
        if(input.vaporSourceNormal.squaredNorm() <= kEpsilon) {
            throw std::invalid_argument("Fuke vapor source normal is zero.");
        }

        const auto started = std::chrono::steady_clock::now();
        FukeResult result;
        result.polygonThicknessMeters.assign(input.polygons.size(), 0.0);
        result.polygonCentroids.resize(input.polygons.size());
        result.statistics.trajectorySampleCount = input.timesSeconds.size();
        result.statistics.evaluatedElementCount = input.polygons.size();
        const Eigen::Vector3d sourceAxis = input.vaporSourceNormal.normalized();
        reportDiagnostic(execution, "Polygon triangulation and visibility",
            "polygons=" + std::to_string(input.polygons.size()));

        for(std::size_t poseIndex = 0;
            poseIndex + 1 < input.workpiecePoses.size(); ++poseIndex) {
            if(canceled(execution)) {
                result.canceled = true;
                break;
            }
            const double dt = std::max(0.0,
                input.timesSeconds[poseIndex + 1] - input.timesSeconds[poseIndex]);
            if(dt <= 0.0) {
                continue;
            }
            if(!input.sprayEnabled.empty() && !input.sprayEnabled[poseIndex]) {
                continue;
            }
            const Eigen::Isometry3d& pose = input.workpiecePoses[poseIndex];
            const TriangleMesh mesh = triangulate(input, pose);
            for(std::size_t polygonIndex = 0;
                polygonIndex < input.polygons.size(); ++polygonIndex) {
                const Eigen::Vector3d centroid = polygonCentroid(
                    input, input.polygons[polygonIndex], pose);
                const Eigen::Vector3d normal = polygonNormal(
                    input, input.polygons[polygonIndex], pose);
                const Eigen::Vector3d sourceToPoint =
                    centroid - input.vaporSourcePosition;
                const double distance = sourceToPoint.norm();
                if(distance <= kEpsilon || normal.squaredNorm() <= kEpsilon) {
                    continue;
                }
                const Eigen::Vector3d direction = sourceToPoint / distance;
                const double plumeCosine = sourceAxis.dot(direction);
                const double incidenceCosine = -normal.dot(direction);
                if(plumeCosine <= 0.0 || incidenceCosine <= 0.0) {
                    continue;
                }
                ++result.statistics.visibilityQueryCount;
                if(!visible(mesh, input.vaporSourcePosition, centroid)) {
                    ++result.statistics.hiddenElementCount;
                    continue;
                }
                const double rate =
                    parameters.referenceThicknessRateMetersPerSecond
                    * std::pow(plumeCosine, parameters.plumeExponent)
                    * std::pow(parameters.referenceDistanceMeters / distance, 2.0)
                    * incidenceCosine;
                result.polygonThicknessMeters[polygonIndex] += rate * dt;
                result.polygonCentroids[polygonIndex] = centroid;
                ++result.statistics.candidatePairCount;
            }
            reportProgress(execution,
                static_cast<double>(poseIndex + 1)
                    / static_cast<double>(input.workpiecePoses.size() - 1),
                "Fuke 2005 reproduction");
        }
        reportDiagnostic(execution, "Result assembly",
            "visibility queries="
                + std::to_string(result.statistics.visibilityQueryCount)
                + ", accepted contributions="
                + std::to_string(result.statistics.candidatePairCount));
        result.statistics.elapsedMilliseconds = elapsedMilliseconds(started);
        return result;
    }
}
