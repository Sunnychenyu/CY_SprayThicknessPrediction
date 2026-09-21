#include <SprayThicknessPrediction/DynamicSurfaceReproduction.h>

#include <pcl/PolygonMesh.h>
#include <pcl/common/centroid.h>
#include <pcl/conversions.h>
#include <pcl/filters/uniform_sampling.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/surface/concave_hull.h>
#include <pcl/surface/poisson.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace spraythickness::published
{
    namespace
    {
        constexpr double kPi = 3.14159265358979323846;
        constexpr double kEpsilon = 1.0e-12;
        using Point = pcl::PointXYZ;
        using PointNormal = pcl::PointNormal;
        using Cloud = pcl::PointCloud<Point>;
        using NormalCloud = pcl::PointCloud<PointNormal>;

        struct DirectionSample
        {
            Eigen::Vector3d localDirection = Eigen::Vector3d::UnitZ();
            double weight{ 0.0 };
        };

        struct ReconstructedRegion
        {
            TriangleMesh mesh;
            std::vector<Eigen::Vector2d> boundary;
            Eigen::Vector3d origin = Eigen::Vector3d::Zero();
            Eigen::Vector3d axisX = Eigen::Vector3d::UnitX();
            Eigen::Vector3d axisY = Eigen::Vector3d::UnitY();
        };

        double elapsedMilliseconds(
            const std::chrono::steady_clock::time_point& started)
        {
            return std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - started).count();
        }

        Point toPoint(const Eigen::Vector3d& value)
        {
            return Point(static_cast<float>(value.x()),
                static_cast<float>(value.y()), static_cast<float>(value.z()));
        }

        Eigen::Vector3d toEigen(const Point& value)
        {
            return { value.x, value.y, value.z };
        }

        void localBasis(const Eigen::Vector3d& axis,
            Eigen::Vector3d& x, Eigen::Vector3d& y)
        {
            const Eigen::Vector3d reference = std::abs(axis.z()) < 0.9
                ? Eigen::Vector3d::UnitZ() : Eigen::Vector3d::UnitY();
            x = reference.cross(axis).normalized();
            y = axis.cross(x).normalized();
        }

        std::vector<DirectionSample> makeDirectionSamples(
            const DynamicSurfaceParameters& parameters)
        {
            std::vector<DirectionSample> samples;
            const double step = parameters.maximumEjectionAngleRadians
                / static_cast<double>(parameters.polarBinCount - 1);
            for(std::size_t polarIndex = 0;
                polarIndex < parameters.polarBinCount; ++polarIndex) {
                const double polar = step * static_cast<double>(polarIndex);
                const std::size_t azimuthCount = polarIndex == 0
                    ? 1 : std::max<std::size_t>(6,
                        static_cast<std::size_t>(std::ceil(
                            2.0 * kPi * std::sin(polar) / step)));
                const double normalized =
                    (polar - parameters.angularMeanRadians)
                    / parameters.angularSigmaRadians;
                const double ringWeight = std::exp(-0.5 * normalized * normalized);
                for(std::size_t azimuthIndex = 0;
                    azimuthIndex < azimuthCount; ++azimuthIndex) {
                    const double azimuth = 2.0 * kPi
                        * static_cast<double>(azimuthIndex)
                        / static_cast<double>(azimuthCount);
                    samples.push_back({ {
                        std::sin(polar) * std::cos(azimuth),
                        std::sin(polar) * std::sin(azimuth),
                        std::cos(polar) }, ringWeight });
                }
            }
            const double total = std::accumulate(samples.begin(), samples.end(),
                0.0, [](double sum, const DirectionSample& sample) {
                    return sum + sample.weight;
                });
            if(total <= 0.0 || !std::isfinite(total)) {
                throw std::invalid_argument(
                    "Dynamic-surface ejection probabilities are invalid.");
            }
            for(DirectionSample& sample : samples) {
                sample.weight /= total;
            }
            return samples;
        }

        double circleOverlapArea(double radius, double distance)
        {
            if(distance >= 2.0 * radius) {
                return 0.0;
            }
            if(distance <= kEpsilon) {
                return kPi * radius * radius;
            }
            const double ratio = std::clamp(distance / (2.0 * radius), 0.0, 1.0);
            return 2.0 * radius * radius * std::acos(ratio)
                - 0.5 * distance
                    * std::sqrt(std::max(0.0,
                        4.0 * radius * radius - distance * distance));
        }

        void redistributePairwiseOverlap(
            std::vector<DynamicDepositedCylinder>& cylinders)
        {
            if(cylinders.size() < 2) {
                return;
            }
            std::vector<double> addedVolume(cylinders.size(), 0.0);
            const Eigen::Vector3d axis =
                cylinders.front().growthDirection.normalized();
            for(std::size_t first = 0; first < cylinders.size(); ++first) {
                for(std::size_t second = first + 1;
                    second < cylinders.size(); ++second) {
                    const Eigen::Vector3d delta =
                        cylinders[second].baseCenter - cylinders[first].baseCenter;
                    const double radial = (delta - delta.dot(axis) * axis).norm();
                    const double area = circleOverlapArea(
                        cylinders[first].radiusMeters, radial);
                    if(area <= 0.0) {
                        continue;
                    }
                    const double firstStart =
                        cylinders[first].baseCenter.dot(axis);
                    const double secondStart =
                        cylinders[second].baseCenter.dot(axis);
                    const double overlapLength = std::max(0.0,
                        std::min(firstStart + cylinders[first].heightMeters,
                            secondStart + cylinders[second].heightMeters)
                        - std::max(firstStart, secondStart));
                    const double sharedVolume = area * overlapLength;
                    addedVolume[first] += 0.5 * sharedVolume;
                    addedVolume[second] += 0.5 * sharedVolume;
                }
            }
            for(std::size_t index = 0; index < cylinders.size(); ++index) {
                const double baseArea = kPi * cylinders[index].radiusMeters
                    * cylinders[index].radiusMeters;
                cylinders[index].heightMeters += addedVolume[index] / baseArea;
            }
        }

        Cloud::Ptr downsample(const std::vector<Eigen::Vector3d>& points,
            const DynamicSurfaceParameters& parameters)
        {
            Cloud::Ptr source(new Cloud);
            source->reserve(points.size());
            for(const Eigen::Vector3d& point : points) {
                source->push_back(toPoint(point));
            }
            Cloud::Ptr voxel(new Cloud);
            pcl::VoxelGrid<Point> voxelFilter;
            voxelFilter.setInputCloud(source);
            const float leaf = static_cast<float>(parameters.voxelLeafMeters);
            voxelFilter.setLeafSize(leaf, leaf, leaf);
            voxelFilter.filter(*voxel);

            Cloud::Ptr sampled(new Cloud);
            pcl::UniformSampling<Point> uniformFilter;
            uniformFilter.setInputCloud(voxel);
            uniformFilter.setRadiusSearch(
                parameters.uniformSamplingRadiusMeters);
            uniformFilter.filter(*sampled);
            return sampled;
        }

        std::vector<std::vector<int>> hybridNeighborhoods(
            const Cloud::Ptr& cloud,
            const DynamicSurfaceParameters& parameters)
        {
            pcl::KdTreeFLANN<Point> tree;
            tree.setInputCloud(cloud);
            std::vector<std::vector<int>> neighborhoods(cloud->size());
            for(std::size_t index = 0; index < cloud->size(); ++index) {
                std::vector<int> indices;
                std::vector<float> distances;
                tree.nearestKSearch((*cloud)[index],
                    static_cast<int>(parameters.normalMaximumNeighbors),
                    indices, distances);
                const double radiusSquared = parameters.normalSearchRadiusMeters
                    * parameters.normalSearchRadiusMeters;
                for(std::size_t neighbor = 0;
                    neighbor < indices.size(); ++neighbor) {
                    if(distances[neighbor] <= radiusSquared) {
                        neighborhoods[index].push_back(indices[neighbor]);
                    }
                }
            }
            return neighborhoods;
        }

        std::vector<Eigen::Vector3d> estimateNormals(
            const Cloud::Ptr& cloud,
            const std::vector<std::vector<int>>& neighborhoods,
            const Eigen::Vector3d& referencePoint)
        {
            std::vector<Eigen::Vector3d> normals(cloud->size());
            for(std::size_t index = 0; index < cloud->size(); ++index) {
                if(neighborhoods[index].size() < 3) {
                    throw std::runtime_error(
                        "Dynamic-surface normal neighborhood contains fewer than three points.");
                }
                Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
                for(const int neighbor : neighborhoods[index]) {
                    centroid += toEigen((*cloud)[neighbor]);
                }
                centroid /= static_cast<double>(neighborhoods[index].size());
                Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
                for(const int neighbor : neighborhoods[index]) {
                    const Eigen::Vector3d centered =
                        toEigen((*cloud)[neighbor]) - centroid;
                    covariance += centered * centered.transpose();
                }
                covariance /= static_cast<double>(neighborhoods[index].size());
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
                Eigen::Vector3d normal = solver.eigenvectors().col(0).normalized();
                if(normal.dot(referencePoint - toEigen((*cloud)[index])) < 0.0) {
                    normal = -normal;
                }
                normals[index] = normal;
            }
            return normals;
        }

        std::vector<bool> detectBoundary(
            const Cloud::Ptr& cloud,
            const std::vector<std::vector<int>>& neighborhoods,
            const std::vector<Eigen::Vector3d>& normals,
            double threshold)
        {
            std::vector<bool> boundary(cloud->size(), false);
            for(std::size_t index = 0; index < cloud->size(); ++index) {
                double angleSum = 0.0;
                double distanceSum = 0.0;
                std::size_t count = 0;
                for(const int neighbor : neighborhoods[index]) {
                    if(static_cast<std::size_t>(neighbor) == index) {
                        continue;
                    }
                    angleSum += std::acos(std::clamp(std::abs(
                        normals[index].dot(normals[neighbor])), 0.0, 1.0));
                    distanceSum += (toEigen((*cloud)[neighbor])
                        - toEigen((*cloud)[index])).norm();
                    ++count;
                }
                if(count == 0 || distanceSum <= kEpsilon) {
                    boundary[index] = true;
                    continue;
                }
                const double gradient = (angleSum / count)
                    / (distanceSum / count);
                boundary[index] = gradient > threshold;
            }
            return boundary;
        }

        std::vector<std::vector<int>> dbscan(
            const Cloud::Ptr& cloud,
            const std::vector<bool>& boundary,
            const DynamicSurfaceParameters& parameters)
        {
            pcl::KdTreeFLANN<Point> tree;
            tree.setInputCloud(cloud);
            std::vector<int> labels(cloud->size(), -1);
            int nextLabel = 0;
            for(std::size_t seed = 0; seed < cloud->size(); ++seed) {
                if(boundary[seed] || labels[seed] != -1) {
                    continue;
                }
                std::vector<int> neighbors;
                std::vector<float> distances;
                tree.radiusSearch((*cloud)[seed], parameters.dbscanRadiusMeters,
                    neighbors, distances);
                if(neighbors.size() < parameters.dbscanMinimumPoints) {
                    labels[seed] = -2;
                    continue;
                }
                const int label = nextLabel++;
                labels[seed] = label;
                std::deque<int> queue(neighbors.begin(), neighbors.end());
                while(!queue.empty()) {
                    const int point = queue.front();
                    queue.pop_front();
                    if(boundary[point]) {
                        continue;
                    }
                    if(labels[point] == -2) {
                        labels[point] = label;
                    }
                    if(labels[point] != -1) {
                        continue;
                    }
                    labels[point] = label;
                    neighbors.clear();
                    distances.clear();
                    tree.radiusSearch((*cloud)[point],
                        parameters.dbscanRadiusMeters, neighbors, distances);
                    if(neighbors.size() >= parameters.dbscanMinimumPoints) {
                        queue.insert(queue.end(), neighbors.begin(), neighbors.end());
                    }
                }
            }
            if(nextLabel == 0) {
                throw std::runtime_error(
                    "Dynamic-surface DBSCAN produced no reconstructable region.");
            }

            for(std::size_t index = 0; index < cloud->size(); ++index) {
                if(!boundary[index]) {
                    continue;
                }
                std::vector<int> neighbors;
                std::vector<float> distances;
                tree.nearestKSearch((*cloud)[index],
                    static_cast<int>(cloud->size()), neighbors, distances);
                for(const int neighbor : neighbors) {
                    if(labels[neighbor] >= 0) {
                        labels[index] = labels[neighbor];
                        break;
                    }
                }
            }

            std::vector<std::vector<int>> clusters(nextLabel);
            for(std::size_t index = 0; index < labels.size(); ++index) {
                if(labels[index] >= 0) {
                    clusters[labels[index]].push_back(static_cast<int>(index));
                }
            }
            clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                [](const auto& cluster) { return cluster.size() < 3; }),
                clusters.end());
            return clusters;
        }

        bool pointInPolygon(const Eigen::Vector2d& point,
            const std::vector<Eigen::Vector2d>& polygon)
        {
            bool inside = false;
            for(std::size_t first = 0, second = polygon.size() - 1;
                first < polygon.size(); second = first++) {
                const Eigen::Vector2d& a = polygon[first];
                const Eigen::Vector2d& b = polygon[second];
                if((a.y() > point.y()) != (b.y() > point.y())
                    && point.x() < (b.x() - a.x())
                            * (point.y() - a.y()) / (b.y() - a.y()) + a.x()) {
                    inside = !inside;
                }
            }
            return inside;
        }

        using MeshEdge = std::pair<std::uint32_t, std::uint32_t>;

        std::vector<std::vector<std::uint32_t>> boundaryLoops(
            const TriangleMesh& mesh)
        {
            std::map<MeshEdge, std::size_t> edgeCounts;
            for(const auto& face : mesh.faces) {
                for(int edge = 0; edge < 3; ++edge) {
                    ++edgeCounts[std::minmax(
                        face[edge], face[(edge + 1) % 3])];
                }
            }

            std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> graph;
            for(const auto& edge : edgeCounts) {
                if(edge.second != 1) {
                    continue;
                }
                graph[edge.first.first].push_back(edge.first.second);
                graph[edge.first.second].push_back(edge.first.first);
            }

            std::set<MeshEdge> visited;
            std::vector<std::vector<std::uint32_t>> loops;
            for(const auto& entry : graph) {
                const std::uint32_t start = entry.first;
                for(const std::uint32_t firstNeighbor : entry.second) {
                    const MeshEdge firstEdge = std::minmax(start, firstNeighbor);
                    if(visited.count(firstEdge) != 0) {
                        continue;
                    }
                    std::vector<std::uint32_t> loop{ start };
                    std::uint32_t previous = start;
                    std::uint32_t current = firstNeighbor;
                    bool closed = false;
                    while(true) {
                        visited.insert(std::minmax(previous, current));
                        if(current == start) {
                            closed = true;
                            break;
                        }
                        loop.push_back(current);
                        const auto found = graph.find(current);
                        if(found == graph.end() || found->second.size() != 2) {
                            break;
                        }
                        const std::uint32_t next = found->second[0] == previous
                            ? found->second[1] : found->second[0];
                        previous = current;
                        current = next;
                        if(loop.size() > graph.size()) {
                            break;
                        }
                    }
                    if(closed && loop.size() >= 3) {
                        loops.push_back(std::move(loop));
                    }
                }
            }
            return loops;
        }

        double projectedLoopArea(const TriangleMesh& mesh,
            const std::vector<std::uint32_t>& loop,
            const ReconstructedRegion& region)
        {
            double twiceArea = 0.0;
            for(std::size_t index = 0; index < loop.size(); ++index) {
                const Eigen::Vector3d first =
                    mesh.vertices[loop[index]] - region.origin;
                const Eigen::Vector3d second =
                    mesh.vertices[loop[(index + 1) % loop.size()]] - region.origin;
                const Eigen::Vector2d first2(
                    first.dot(region.axisX), first.dot(region.axisY));
                const Eigen::Vector2d second2(
                    second.dot(region.axisX), second.dot(region.axisY));
                twiceArea += first2.x() * second2.y()
                    - second2.x() * first2.y();
            }
            return 0.5 * twiceArea;
        }

        void fillSmallHoles(ReconstructedRegion& region,
            const Eigen::Vector3d& referencePoint,
            std::size_t maximumBoundaryEdges)
        {
            const auto loops = boundaryLoops(region.mesh);
            if(loops.size() < 2) {
                return;
            }

            std::size_t outerLoop = 0;
            double outerArea = 0.0;
            for(std::size_t index = 0; index < loops.size(); ++index) {
                const double area = std::abs(projectedLoopArea(
                    region.mesh, loops[index], region));
                if(area > outerArea) {
                    outerArea = area;
                    outerLoop = index;
                }
            }

            for(std::size_t loopIndex = 0; loopIndex < loops.size(); ++loopIndex) {
                const auto& loop = loops[loopIndex];
                if(loopIndex == outerLoop || loop.size() > maximumBoundaryEdges) {
                    continue;
                }
                Eigen::Vector3d center = Eigen::Vector3d::Zero();
                for(const std::uint32_t vertex : loop) {
                    center += region.mesh.vertices[vertex];
                }
                center /= static_cast<double>(loop.size());
                const std::uint32_t centerIndex =
                    static_cast<std::uint32_t>(region.mesh.vertices.size());
                region.mesh.vertices.push_back(center);
                for(std::size_t index = 0; index < loop.size(); ++index) {
                    std::array<std::uint32_t, 3> face{ loop[index],
                        loop[(index + 1) % loop.size()], centerIndex };
                    const Eigen::Vector3d normal =
                        (region.mesh.vertices[face[1]]
                            - region.mesh.vertices[face[0]]).cross(
                            region.mesh.vertices[face[2]]
                                - region.mesh.vertices[face[0]]);
                    if(normal.dot(referencePoint - center) < 0.0) {
                        std::swap(face[0], face[1]);
                    }
                    region.mesh.faces.push_back(face);
                }
            }
        }

        ReconstructedRegion reconstructCluster(
            const Cloud::Ptr& cloud,
            const std::vector<Eigen::Vector3d>& normals,
            const std::vector<bool>& sourceBoundary,
            const std::vector<int>& cluster,
            const Eigen::Vector3d& boundaryReferencePoint,
            const Eigen::Vector3d& normalReferencePoint,
            const DynamicSurfaceParameters& parameters)
        {
            Cloud::Ptr points(new Cloud);
            NormalCloud::Ptr pointsWithNormals(new NormalCloud);
            for(const int index : cluster) {
                points->push_back((*cloud)[index]);
                PointNormal value;
                value.x = (*cloud)[index].x;
                value.y = (*cloud)[index].y;
                value.z = (*cloud)[index].z;
                value.normal_x = static_cast<float>(normals[index].x());
                value.normal_y = static_cast<float>(normals[index].y());
                value.normal_z = static_cast<float>(normals[index].z());
                pointsWithNormals->push_back(value);
            }
            if(points->size() < 10) {
                throw std::runtime_error(
                    "Dynamic-surface reconstruction region contains fewer than ten points.");
            }

            Eigen::Vector4f centroid4;
            pcl::compute3DCentroid(*points, centroid4);
            Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
            const Eigen::Vector3d centroid = centroid4.head<3>().cast<double>();
            for(const Point& point : *points) {
                const Eigen::Vector3d centered = toEigen(point) - centroid;
                covariance += centered * centered.transpose();
            }
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
            ReconstructedRegion region;
            region.origin = centroid;
            region.axisX = solver.eigenvectors().col(2).normalized();
            region.axisY = solver.eigenvectors().col(1).normalized();

            pcl::ConcaveHull<Point> alphaShape;
            alphaShape.setInputCloud(points);
            alphaShape.setAlpha(parameters.alphaShapeRadiusMeters);
            alphaShape.setDimension(2);
            Cloud hullPoints;
            std::vector<pcl::Vertices> hullPolygons;
            alphaShape.reconstruct(hullPoints, hullPolygons);
            if(hullPolygons.empty()) {
                throw std::runtime_error(
                    "Dynamic-surface Alpha Shape produced no boundary.");
            }
            const Eigen::Vector3d referenceRelative =
                boundaryReferencePoint - region.origin;
            const Eigen::Vector2d reference2(
                referenceRelative.dot(region.axisX),
                referenceRelative.dot(region.axisY));
            double closestBoundarySquared = std::numeric_limits<double>::max();
            for(const pcl::Vertices& polygon : hullPolygons) {
                std::vector<Eigen::Vector2d> candidate;
                Eigen::Vector2d candidateCenter = Eigen::Vector2d::Zero();
                for(const std::uint32_t index : polygon.vertices) {
                    const Eigen::Vector3d relative =
                        toEigen(hullPoints[index]) - region.origin;
                    candidate.push_back({ relative.dot(region.axisX),
                        relative.dot(region.axisY) });
                    candidateCenter += candidate.back();
                }
                if(candidate.size() < 3) {
                    continue;
                }
                candidateCenter /= static_cast<double>(candidate.size());
                if(pointInPolygon(reference2, candidate)) {
                    region.boundary = std::move(candidate);
                    break;
                }
                const double distanceSquared =
                    (candidateCenter - reference2).squaredNorm();
                if(distanceSquared < closestBoundarySquared) {
                    closestBoundarySquared = distanceSquared;
                    region.boundary = std::move(candidate);
                }
            }
            if(region.boundary.size() < 3) {
                throw std::runtime_error(
                    "Dynamic-surface Alpha Shape produced no valid boundary loop.");
            }

            pcl::Poisson<PointNormal> poisson;
            poisson.setDepth(parameters.poissonDepth);
            poisson.setInputCloud(pointsWithNormals);
            pcl::PolygonMesh reconstructed;
            poisson.reconstruct(reconstructed);
            Cloud reconstructedPoints;
            pcl::fromPCLPointCloud2(reconstructed.cloud, reconstructedPoints);
            region.mesh.vertices.reserve(reconstructedPoints.size());
            for(const Point& point : reconstructedPoints) {
                region.mesh.vertices.emplace_back(point.x, point.y, point.z);
            }
            for(const pcl::Vertices& polygon : reconstructed.polygons) {
                if(polygon.vertices.size() < 3) {
                    continue;
                }
                for(std::size_t corner = 1;
                    corner + 1 < polygon.vertices.size(); ++corner) {
                    std::array<std::uint32_t, 3> face{
                        static_cast<std::uint32_t>(polygon.vertices[0]),
                        static_cast<std::uint32_t>(polygon.vertices[corner]),
                        static_cast<std::uint32_t>(
                            polygon.vertices[corner + 1]) };
                    const Eigen::Vector3d faceCenter =
                        (region.mesh.vertices[face[0]]
                            + region.mesh.vertices[face[1]]
                            + region.mesh.vertices[face[2]]) / 3.0;
                    const Eigen::Vector3d faceDirection =
                        (region.mesh.vertices[face[1]]
                            - region.mesh.vertices[face[0]]).cross(
                            region.mesh.vertices[face[2]]
                                - region.mesh.vertices[face[0]]);
                    if(faceDirection.dot(normalReferencePoint - faceCenter) < 0.0) {
                        std::swap(face[1], face[2]);
                    }
                    region.mesh.faces.push_back(face);
                }
            }
            if(region.mesh.faces.empty()) {
                throw std::runtime_error(
                    "Dynamic-surface Poisson reconstruction produced no faces.");
            }

            fillSmallHoles(region, normalReferencePoint,
                parameters.maximumHoleBoundaryEdges);

            std::vector<std::unordered_set<std::size_t>> adjacency(
                region.mesh.vertices.size());
            std::map<std::pair<std::uint32_t, std::uint32_t>, std::size_t> edges;
            for(const auto& face : region.mesh.faces) {
                for(int edge = 0; edge < 3; ++edge) {
                    const std::uint32_t first = face[edge];
                    const std::uint32_t second = face[(edge + 1) % 3];
                    adjacency[first].insert(second);
                    adjacency[second].insert(first);
                    ++edges[std::minmax(first, second)];
                }
            }
            std::vector<bool> fixed(region.mesh.vertices.size(), false);
            for(const auto& edge : edges) {
                if(edge.second == 1) {
                    fixed[edge.first.first] = true;
                    fixed[edge.first.second] = true;
                }
            }
            for(const int sourceIndex : cluster) {
                if(!sourceBoundary[sourceIndex]) {
                    continue;
                }
                const Eigen::Vector3d source = toEigen((*cloud)[sourceIndex]);
                std::size_t closest = 0;
                double closestSquared = std::numeric_limits<double>::max();
                for(std::size_t vertex = 0;
                    vertex < region.mesh.vertices.size(); ++vertex) {
                    const double squared =
                        (region.mesh.vertices[vertex] - source).squaredNorm();
                    if(squared < closestSquared) {
                        closestSquared = squared;
                        closest = vertex;
                    }
                }
                fixed[closest] = true;
            }
            for(std::size_t iteration = 0;
                iteration < parameters.smoothingIterations; ++iteration) {
                std::vector<Eigen::Vector3d> next = region.mesh.vertices;
                for(std::size_t vertex = 0;
                    vertex < region.mesh.vertices.size(); ++vertex) {
                    if(fixed[vertex] || adjacency[vertex].empty()) {
                        continue;
                    }
                    Eigen::Vector3d average = Eigen::Vector3d::Zero();
                    for(const std::size_t neighbor : adjacency[vertex]) {
                        average += region.mesh.vertices[neighbor];
                    }
                    average /= static_cast<double>(adjacency[vertex].size());
                    next[vertex] = (1.0 - parameters.smoothingRelaxation)
                            * region.mesh.vertices[vertex]
                        + parameters.smoothingRelaxation * average;
                }
                region.mesh.vertices = std::move(next);
            }

            std::vector<std::array<std::uint32_t, 3>> trimmedFaces;
            trimmedFaces.reserve(region.mesh.faces.size());
            for(const auto& face : region.mesh.faces) {
                const Eigen::Vector3d faceCenter =
                    (region.mesh.vertices[face[0]]
                        + region.mesh.vertices[face[1]]
                        + region.mesh.vertices[face[2]]) / 3.0;
                const Eigen::Vector3d relative = faceCenter - region.origin;
                if(pointInPolygon({ relative.dot(region.axisX),
                        relative.dot(region.axisY) }, region.boundary)) {
                    trimmedFaces.push_back(face);
                }
            }
            region.mesh.faces = std::move(trimmedFaces);
            if(region.mesh.faces.empty()) {
                throw std::runtime_error(
                    "Dynamic-surface Poisson mesh was removed by boundary trimming.");
            }
            return region;
        }

        TriangleMesh mergeRegions(const TriangleMesh& current,
            const std::vector<ReconstructedRegion>& regions)
        {
            TriangleMesh result;
            result.vertices = current.vertices;
            for(std::size_t faceIndex = 0;
                faceIndex < current.faces.size(); ++faceIndex) {
                const Eigen::Vector3d center = faceCentroid(current, faceIndex);
                bool replaced = false;
                for(const ReconstructedRegion& region : regions) {
                    const Eigen::Vector3d relative = center - region.origin;
                    if(pointInPolygon({ relative.dot(region.axisX),
                            relative.dot(region.axisY) }, region.boundary)) {
                        replaced = true;
                        break;
                    }
                }
                if(!replaced) {
                    result.faces.push_back(current.faces[faceIndex]);
                }
            }
            for(const ReconstructedRegion& region : regions) {
                const std::uint32_t offset =
                    static_cast<std::uint32_t>(result.vertices.size());
                result.vertices.insert(result.vertices.end(),
                    region.mesh.vertices.begin(), region.mesh.vertices.end());
                for(auto face : region.mesh.faces) {
                    face[0] += offset;
                    face[1] += offset;
                    face[2] += offset;
                    result.faces.push_back(face);
                }
            }
            result.validate();
            return result;
        }

        TriangleMesh reconstructSurface(const TriangleMesh& current,
            const std::vector<DynamicDepositedCylinder>& batch,
            const Eigen::Vector3d& nozzlePosition,
            const Eigen::Vector3d& nozzleAxis,
            const DynamicSurfaceParameters& parameters)
        {
            std::vector<Eigen::Vector3d> featurePoints;
            featurePoints.reserve(batch.size());
            for(const DynamicDepositedCylinder& cylinder : batch) {
                featurePoints.push_back(cylinder.baseCenter
                    + cylinder.heightMeters * cylinder.growthDirection);
            }
            Cloud::Ptr cloud = downsample(featurePoints, parameters);
            if(cloud->size() < 10) {
                throw std::runtime_error(
                    "Dynamic-surface downsampling retained fewer than ten feature points.");
            }
            const auto neighborhoods = hybridNeighborhoods(cloud, parameters);
            const auto normals = estimateNormals(
                cloud, neighborhoods, nozzlePosition);
            const auto boundary = detectBoundary(cloud, neighborhoods, normals,
                parameters.boundaryGradientThresholdRadiansPerMeter);
            const auto clusters = dbscan(cloud, boundary, parameters);
            std::vector<ReconstructedRegion> regions;
            regions.reserve(clusters.size());

            RayHit referenceHit;
            Eigen::Vector3d realReferencePoint = batch.front().baseCenter;
            if(nearestRayHit(current, nozzlePosition, nozzleAxis, referenceHit)) {
                realReferencePoint = referenceHit.position;
            } else {
                double closestAxisSquared = std::numeric_limits<double>::max();
                for(const DynamicDepositedCylinder& cylinder : batch) {
                    const Eigen::Vector3d delta =
                        cylinder.baseCenter - nozzlePosition;
                    const double axisDistanceSquared =
                        (delta - delta.dot(nozzleAxis) * nozzleAxis).squaredNorm();
                    if(axisDistanceSquared < closestAxisSquared) {
                        closestAxisSquared = axisDistanceSquared;
                        realReferencePoint = cylinder.baseCenter;
                    }
                }
            }
            std::size_t realReferenceCloudIndex = 0;
            double referenceDistanceSquared = std::numeric_limits<double>::max();
            for(std::size_t index = 0; index < cloud->size(); ++index) {
                const double distanceSquared =
                    (toEigen((*cloud)[index]) - realReferencePoint).squaredNorm();
                if(distanceSquared < referenceDistanceSquared) {
                    referenceDistanceSquared = distanceSquared;
                    realReferenceCloudIndex = index;
                }
            }
            for(const auto& cluster : clusters) {
                const bool containsRealReference = std::find(cluster.begin(),
                    cluster.end(), static_cast<int>(realReferenceCloudIndex))
                    != cluster.end();
                Eigen::Vector3d localReference = realReferencePoint;
                if(!containsRealReference) {
                    localReference = Eigen::Vector3d::Zero();
                    localReference.z() = -std::numeric_limits<double>::max();
                    for(const int index : cluster) {
                        const Eigen::Vector3d point = toEigen((*cloud)[index]);
                        localReference.x() += point.x();
                        localReference.y() += point.y();
                        localReference.z() = std::max(
                            localReference.z(), point.z());
                    }
                    localReference.x() /= static_cast<double>(cluster.size());
                    localReference.y() /= static_cast<double>(cluster.size());
                }
                regions.push_back(reconstructCluster(cloud, normals, boundary,
                    cluster, localReference, nozzlePosition, parameters));
            }
            return mergeRegions(current, regions);
        }
    }

    void DynamicSurfaceParameters::validate() const
    {
        relativeBuildUpByInclinationRadians.validate(
            "Dynamic-surface relative build-up rate Rbu");
        const std::vector<double> positiveValues{
            equivalentParticleDiameterMeters,
            depositedCylinderRadiusMeters,
            equivalentParticleRatePerSecond,
            angularSigmaRadians,
            maximumEjectionAngleRadians,
            depositionBatchDurationSeconds,
            voxelLeafMeters,
            uniformSamplingRadiusMeters,
            normalSearchRadiusMeters,
            boundaryGradientThresholdRadiansPerMeter,
            dbscanRadiusMeters,
            alphaShapeRadiusMeters
        };
        if(std::any_of(positiveValues.begin(), positiveValues.end(),
                [](double value) {
                    return !std::isfinite(value) || value <= 0.0;
                })
            || !std::isfinite(angularMeanRadians)
            || polarBinCount < 2
            || normalMaximumNeighbors < 3
            || dbscanMinimumPoints < 3
            || poissonDepth < 3
            || maximumHoleBoundaryEdges < 3
            || !std::isfinite(smoothingRelaxation)
            || smoothingRelaxation <= 0.0
            || smoothingRelaxation > 1.0) {
            throw std::invalid_argument(
                "Dynamic-surface reproduction requires all particle, batch, "
                "overlap, point-cloud, clustering, and reconstruction parameters.");
        }
    }

    DynamicSurfaceResult DynamicSurfaceReproducer::run(
        const DynamicSurfaceInputModel& input,
        const DynamicSurfaceParameters& parameters,
        const ReproductionExecution& execution)
    {
        parameters.validate();
        input.initialStlSurface.validate();
        if(input.nozzleTrajectory.size() < 2) {
            throw std::invalid_argument(
                "Dynamic-surface reproduction requires at least two timed poses.");
        }
        const auto started = std::chrono::steady_clock::now();
        DynamicSurfaceResult result;
        result.evolvedStlSurface = input.initialStlSurface;
        result.statistics.trajectorySampleCount = input.nozzleTrajectory.size();
        result.implementationNotes.push_back(
            "The paper does not publish the overlap-volume solver; parallel-cylinder "
            "pair intersections are redistributed equally to the participating cylinders.");
        result.implementationNotes.push_back(
            "The paper does not publish hole-repair and feature-smoothing numerical details; "
            "configured small boundary loops are closed by centroid triangle fans, and "
            "Laplacian smoothing fixes both the reconstructed boundary and mapped gradient features.");
        result.implementationNotes.push_back(
            "Subregions use the central-ray reference point or the paper's Eq. (19) virtual "
            "reference point; if no Alpha Shape loop contains it, the nearest loop centroid is selected.");

        const std::vector<DirectionSample> directions =
            makeDirectionSamples(parameters);
        std::vector<double> probabilities;
        probabilities.reserve(directions.size());
        for(const DirectionSample& direction : directions) {
            probabilities.push_back(direction.weight);
        }
        std::mt19937 generator(parameters.randomSeed);
        std::discrete_distribution<std::size_t> selectDirection(
            probabilities.begin(), probabilities.end());
        const double sphereVolume = kPi / 6.0
            * std::pow(parameters.equivalentParticleDiameterMeters, 3.0);
        const double baseCylinderHeight = sphereVolume
            / (kPi * parameters.depositedCylinderRadiusMeters
                * parameters.depositedCylinderRadiusMeters);

        for(std::size_t poseIndex = 0;
            poseIndex + 1 < input.nozzleTrajectory.size(); ++poseIndex) {
            if(canceled(execution)) {
                result.canceled = true;
                break;
            }
            const SprayPose& pose = input.nozzleTrajectory[poseIndex];
            const double interval =
                input.nozzleTrajectory[poseIndex + 1].timeSeconds
                - pose.timeSeconds;
            if(interval <= 0.0) {
                throw std::invalid_argument(
                    "Dynamic-surface trajectory times must be strictly increasing.");
            }
            if(!pose.sprayEnabled) {
                continue;
            }
            const Eigen::Vector3d axis = pose.axis.normalized();
            Eigen::Vector3d basisX;
            Eigen::Vector3d basisY;
            localBasis(axis, basisX, basisY);
            const std::size_t batchCount = std::max<std::size_t>(1,
                static_cast<std::size_t>(std::ceil(
                    interval / parameters.depositionBatchDurationSeconds)));
            const std::size_t totalParticles = static_cast<std::size_t>(
                std::llround(parameters.equivalentParticleRatePerSecond * interval));
            for(std::size_t batchIndex = 0;
                batchIndex < batchCount; ++batchIndex) {
                const std::size_t begin = totalParticles * batchIndex / batchCount;
                const std::size_t end = totalParticles * (batchIndex + 1) / batchCount;
                std::vector<DynamicDepositedCylinder> batch;
                batch.reserve(end - begin);
                for(std::size_t particle = begin; particle < end; ++particle) {
                    const DirectionSample& sample =
                        directions[selectDirection(generator)];
                    const Eigen::Vector3d direction =
                        (sample.localDirection.x() * basisX
                            + sample.localDirection.y() * basisY
                            + sample.localDirection.z() * axis).normalized();
                    RayHit hit;
                    ++result.statistics.visibilityQueryCount;
                    if(!nearestRayHit(result.evolvedStlSurface,
                            pose.position, direction, hit)) {
                        continue;
                    }
                    const double inclination = std::acos(std::clamp(
                        std::abs(hit.normal.dot(axis)), 0.0, 1.0));
                    const double buildUp =
                        parameters.relativeBuildUpByInclinationRadians.interpolate(
                            inclination,
                            "Dynamic-surface relative build-up rate Rbu");
                    if(buildUp <= 0.0) {
                        continue;
                    }
                    batch.push_back({ hit.position, -axis,
                        parameters.depositedCylinderRadiusMeters,
                        baseCylinderHeight * buildUp });
                }
                if(batch.empty()) {
                    continue;
                }
                redistributePairwiseOverlap(batch);
                result.evolvedStlSurface = reconstructSurface(
                    result.evolvedStlSurface, batch, pose.position, axis,
                    parameters);
                result.depositedCylinders.insert(
                    result.depositedCylinders.end(), batch.begin(), batch.end());
                result.statistics.candidatePairCount += batch.size();
            }
            reportProgress(execution,
                static_cast<double>(poseIndex + 1)
                    / static_cast<double>(input.nozzleTrajectory.size() - 1),
                "Dynamic surface evolution reproduction");
        }
        result.statistics.evaluatedElementCount =
            result.evolvedStlSurface.faces.size();
        result.statistics.elapsedMilliseconds = elapsedMilliseconds(started);
        return result;
    }
}
