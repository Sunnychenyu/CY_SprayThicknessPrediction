#include <SprayThicknessPrediction/DynamicSurfaceReproduction.h>
#include <SprayThicknessPrediction/FukeReproduction.h>
#include <SprayThicknessPrediction/TzinavaReproduction.h>
#include <SprayThicknessPrediction/VanerioReproduction.h>
#include <SprayThicknessPrediction/WuReproduction.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
    bool approximatelyEqual(double first, double second, double tolerance = 1.0e-10)
    {
        return std::abs(first - second) <= tolerance;
    }

    spraythickness::published::TriangleMesh makePlaneMesh()
    {
        spraythickness::published::TriangleMesh mesh;
        mesh.vertices = {
            Eigen::Vector3d(-0.02, -0.02, 0.0),
            Eigen::Vector3d(0.02, -0.02, 0.0),
            Eigen::Vector3d(0.02, 0.02, 0.0),
            Eigen::Vector3d(-0.02, 0.02, 0.0)
        };
        mesh.faces = { { 0, 1, 2 }, { 0, 2, 3 } };
        return mesh;
    }

    bool testFukePublishedEquation()
    {
        using namespace spraythickness::published;
        FukeInputModel input;
        input.vertices = {
            Eigen::Vector3d(-0.01, -0.005, 0.1),
            Eigen::Vector3d(0.01, -0.005, 0.1),
            Eigen::Vector3d(0.0, 0.01, 0.1)
        };
        input.polygons.push_back({ { 0, 2, 1 } });
        input.timesSeconds = { 0.0, 1.0 };
        input.workpiecePoses = {
            Eigen::Isometry3d::Identity(), Eigen::Isometry3d::Identity()
        };
        FukeParameters parameters;
        parameters.referenceThicknessRateMetersPerSecond = 2.0e-6;
        parameters.referenceDistanceMeters = 0.1;
        const FukeResult result = FukeReproducer::run(input, parameters);
        return result.polygonThicknessMeters.size() == 1
            && approximatelyEqual(result.polygonThicknessMeters.front(),
                2.0e-6, 1.0e-12);
    }

    bool testFukeMillimeterFace()
    {
        using namespace spraythickness::published;
        FukeInputModel input;
        input.vertices = {
            Eigen::Vector3d(0.0, 0.0, 0.1),
            Eigen::Vector3d(0.0, 0.001, 0.1),
            Eigen::Vector3d(0.001, 0.0, 0.1)
        };
        input.polygons.push_back({ { 0, 1, 2 } });
        input.timesSeconds = { 0.0, 1.0 };
        input.workpiecePoses = {
            Eigen::Isometry3d::Identity(), Eigen::Isometry3d::Identity()
        };
        FukeParameters parameters;
        parameters.referenceThicknessRateMetersPerSecond = 2.0e-6;
        parameters.referenceDistanceMeters = 0.1;
        const FukeResult result = FukeReproducer::run(input, parameters);
        return result.polygonThicknessMeters.size() == 1
            && result.polygonThicknessMeters.front() > 0.0;
    }

    bool testFukeSourceLineOfSight()
    {
        using namespace spraythickness::published;
        FukeInputModel input;
        input.vertices = {
            { -0.01, -0.01, 0.1 }, { 0.0, 0.01, 0.1 },
            { 0.01, -0.01, 0.1 }, { -0.01, -0.01, 0.2 },
            { 0.0, 0.01, 0.2 }, { 0.01, -0.01, 0.2 }
        };
        input.polygons = { { { 0, 1, 2 } }, { { 3, 4, 5 } } };
        input.timesSeconds = { 0.0, 1.0 };
        input.workpiecePoses = {
            Eigen::Isometry3d::Identity(), Eigen::Isometry3d::Identity()
        };
        FukeParameters parameters;
        parameters.referenceThicknessRateMetersPerSecond = 2.0e-6;
        parameters.referenceDistanceMeters = 0.1;
        const FukeResult result = FukeReproducer::run(input, parameters);
        const bool passed = result.polygonThicknessMeters.size() == 2
            && result.polygonThicknessMeters[0] > 0.0
            && approximatelyEqual(result.polygonThicknessMeters[1],
                0.0, 1.0e-12);
        if(!passed && result.polygonThicknessMeters.size() == 2) {
            std::cerr << "Fuke front/back thickness: "
                      << result.polygonThicknessMeters[0] << ", "
                      << result.polygonThicknessMeters[1]
                      << "; blocked=" << result.statistics.hiddenElementCount
                      << '\n';
        }
        return passed;
    }

    bool testFukeTimeAccumulation()
    {
        using namespace spraythickness::published;
        FukeInputModel input;
        input.vertices = {
            { -0.01, -0.01, 0.1 }, { 0.0, 0.02, 0.1 },
            { 0.01, -0.01, 0.1 }
        };
        input.polygons = { { { 0, 1, 2 } } };
        input.timesSeconds = { 0.0, 0.25, 1.0 };
        input.workpiecePoses.assign(3, Eigen::Isometry3d::Identity());
        FukeParameters parameters;
        parameters.referenceThicknessRateMetersPerSecond = 2.0e-6;
        parameters.referenceDistanceMeters = 0.1;
        const FukeResult result = FukeReproducer::run(input, parameters);
        const bool passed = result.polygonThicknessMeters.size() == 1
            && approximatelyEqual(result.polygonThicknessMeters[0], 2.0e-6,
                1.0e-12)
            && result.statistics.candidatePairCount == 2;
        if(!passed) {
            std::cerr << "Fuke accumulated thickness: "
                      << result.polygonThicknessMeters[0]
                      << "; contributions="
                      << result.statistics.candidatePairCount << '\n';
        }
        return passed;
    }

    bool testFukeDistanceAndAngles()
    {
        using namespace spraythickness::published;
        FukeInputModel input;
        input.vertices = {
            { -0.01, -0.01, 0.0 }, { 0.0, 0.02, 0.0 },
            { 0.01, -0.01, 0.0 }
        };
        input.polygons = { { { 0, 1, 2 } } };
        input.timesSeconds = { 0.0, 1.0 };
        Eigen::Isometry3d pose = Eigen::Isometry3d::Identity();
        pose.linear() = Eigen::AngleAxisd(
            3.14159265358979323846 / 3.0, Eigen::Vector3d::UnitY()).toRotationMatrix();
        pose.translation() = Eigen::Vector3d(0.0, 0.0, 0.2);
        input.workpiecePoses.assign(2, pose);
        FukeParameters parameters;
        parameters.referenceThicknessRateMetersPerSecond = 2.0e-6;
        parameters.referenceDistanceMeters = 0.1;
        const FukeResult tilted = FukeReproducer::run(input, parameters);
        pose = Eigen::Isometry3d::Identity();
        pose.translation() = Eigen::Vector3d(0.1, 0.0, 0.1);
        input.workpiecePoses.assign(2, pose);
        const FukeResult offAxis = FukeReproducer::run(input, parameters);
        return tilted.polygonThicknessMeters.size() == 1
            && offAxis.polygonThicknessMeters.size() == 1
            && approximatelyEqual(tilted.polygonThicknessMeters[0], 0.25e-6,
                1.0e-12)
            && approximatelyEqual(offAxis.polygonThicknessMeters[0],
                2.0e-6 / (4.0 * std::sqrt(2.0)), 1.0e-12);
    }

    bool testFukeDensePlateOcclusion()
    {
        using namespace spraythickness::published;
        FukeInputModel input;
        constexpr int cells = 12;
        constexpr double spacing = 0.002;
        for(double height : { 0.1, 0.2 }) {
            const std::uint32_t base =
                static_cast<std::uint32_t>(input.vertices.size());
            for(int y = 0; y <= cells; ++y) {
                for(int x = 0; x <= cells; ++x) {
                    input.vertices.emplace_back(
                        (x - cells / 2) * spacing,
                        (y - cells / 2) * spacing, height);
                }
            }
            for(int y = 0; y < cells; ++y) {
                for(int x = 0; x < cells; ++x) {
                    const std::uint32_t a = base + y * (cells + 1) + x;
                    input.polygons.push_back({ { a, a + cells + 1, a + 1 } });
                    input.polygons.push_back({ { a + 1, a + cells + 1,
                        a + cells + 2 } });
                }
            }
        }
        input.timesSeconds = { 0.0, 1.0 };
        input.workpiecePoses.assign(2, Eigen::Isometry3d::Identity());
        FukeParameters parameters;
        parameters.referenceThicknessRateMetersPerSecond = 2.0e-6;
        parameters.referenceDistanceMeters = 0.1;
        const FukeResult result = FukeReproducer::run(input, parameters);
        std::cout << "Fuke dense-plate elapsed: "
                  << result.statistics.elapsedMilliseconds << " ms\n";
        const std::size_t frontFaceCount = 2 * cells * cells;
        return result.polygonThicknessMeters.size() == 2 * frontFaceCount
            && std::all_of(result.polygonThicknessMeters.begin(),
                result.polygonThicknessMeters.begin() + frontFaceCount,
                [](double thickness) { return thickness > 0.0; })
            && std::all_of(result.polygonThicknessMeters.begin() + frontFaceCount,
                result.polygonThicknessMeters.end(),
                [](double thickness) { return thickness == 0.0; });
    }

    bool testFukeMovingDensePlateOcclusion()
    {
        using namespace spraythickness::published;
        FukeInputModel input;
        constexpr int cells = 14;
        constexpr double spacing = 0.001;
        for(double height : { 0.12, 0.14 }) {
            const std::uint32_t base =
                static_cast<std::uint32_t>(input.vertices.size());
            for(int y = 0; y <= cells; ++y) {
                for(int x = 0; x <= cells; ++x) {
                    input.vertices.emplace_back(
                        (x - cells / 2) * spacing,
                        (y - cells / 2) * spacing, height);
                }
            }
            for(int y = 0; y < cells; ++y) {
                for(int x = 0; x < cells; ++x) {
                    const std::uint32_t a = base + y * (cells + 1) + x;
                    input.polygons.push_back({ { a, a + cells + 1, a + 1 } });
                    input.polygons.push_back({ { a + 1, a + cells + 1,
                        a + cells + 2 } });
                }
            }
        }
        for(int step = -6; step <= 7; ++step) {
            Eigen::Isometry3d pose = Eigen::Isometry3d::Identity();
            pose.translation().x() = -step * 0.0005;
            input.timesSeconds.push_back((step + 6) * 0.01);
            input.workpiecePoses.push_back(pose);
        }
        FukeParameters parameters;
        parameters.referenceThicknessRateMetersPerSecond = 2.0e-6;
        parameters.referenceDistanceMeters = 0.12;
        const FukeResult result = FukeReproducer::run(input, parameters);
        const std::size_t topFaces = 2 * cells * cells;
        const std::size_t leakedFaces = std::count_if(
            result.polygonThicknessMeters.begin() + topFaces,
            result.polygonThicknessMeters.end(),
            [](double thickness) { return thickness > 0.0; });
        if(leakedFaces != 0) {
            std::cerr << "Fuke moving-plate shadow leaked onto "
                      << leakedFaces << " lower faces.\n";
        }
        const bool interiorOccluded =
            result.polygonThicknessMeters.size() == 2 * topFaces
            && leakedFaces == 0
            && std::all_of(result.polygonThicknessMeters.begin(),
                result.polygonThicknessMeters.begin() + topFaces,
                [](double thickness) { return thickness > 0.0; });
        if(!interiorOccluded) {
            return false;
        }
        Eigen::Isometry3d outsidePose = Eigen::Isometry3d::Identity();
        outsidePose.translation().x() = -0.017;
        input.timesSeconds = { 0.0, 1.0 };
        input.workpiecePoses.assign(2, outsidePose);
        const FukeResult exposedEdge = FukeReproducer::run(input, parameters);
        const std::size_t centerFace = topFaces
            + 2 * ((cells / 2) * cells + cells / 2);
        return exposedEdge.polygonThicknessMeters[centerFace] == 0.0
            && std::any_of(exposedEdge.polygonThicknessMeters.begin() + topFaces,
                exposedEdge.polygonThicknessMeters.end(),
                [](double thickness) { return thickness > 0.0; });
    }

    bool testFukeQuadrilateralVisibility()
    {
        using namespace spraythickness::published;
        FukeInputModel input;
        input.vertices = {
            { -0.01, -0.01, 0.1 }, { -0.01, 0.01, 0.1 },
            { 0.01, 0.01, 0.1 }, { 0.01, -0.01, 0.1 },
            { -0.01, -0.01, 0.2 }, { -0.01, 0.01, 0.2 },
            { 0.01, 0.01, 0.2 }, { 0.01, -0.01, 0.2 }
        };
        input.polygons = { { { 0, 1, 2, 3 } }, { { 4, 5, 6, 7 } } };
        input.timesSeconds = { 0.0, 1.0 };
        input.workpiecePoses.assign(2, Eigen::Isometry3d::Identity());
        FukeParameters parameters;
        parameters.referenceThicknessRateMetersPerSecond = 2.0e-6;
        parameters.referenceDistanceMeters = 0.1;
        const FukeResult result = FukeReproducer::run(input, parameters);
        return result.polygonThicknessMeters.size() == 2
            && approximatelyEqual(result.polygonThicknessMeters[0], 2.0e-6,
                1.0e-12)
            && result.polygonThicknessMeters[1] == 0.0
            && result.statistics.hiddenElementCount == 1;
    }

    bool testTzinavaStationaryBranch()
    {
        using namespace spraythickness::published;
        TzinavaInputModel input;
        input.initialMesh.vertices = {
            Eigen::Vector3d(-0.005, -0.005, 0.0),
            Eigen::Vector3d(0.005, -0.005, 0.0),
            Eigen::Vector3d(0.005, 0.005, 0.0),
            Eigen::Vector3d(-0.005, 0.005, 0.0)
        };
        input.initialMesh.faces = { { 0, 1, 2 }, { 0, 2, 3 } };
        input.gunTrajectory.push_back({ 0.0,
            Eigen::Vector3d(0.0, 0.0, 0.1),
            -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        input.gunTrajectory.push_back({ 1.0,
            Eigen::Vector3d(0.0, 0.0, 0.1),
            -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        TzinavaParameters parameters;
        parameters.beamRadiusMeters = 0.05;
        parameters.gaussianRadialProfile = false;
        parameters.stationaryTimeStepSeconds = 0.25;
        parameters.lookupReferenceDwellSeconds = 1.0;
        parameters.thicknessTable.standOffDistancesMeters = { 0.09, 0.11 };
        parameters.thicknessTable.impactAnglesDegrees = { 80.0, 90.0 };
        parameters.thicknessTable.thicknessMeters = {
            1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6
        };
        const TzinavaResult result = TzinavaReproducer::run(input, parameters);
        input.gunTrajectory.back().position.z() = 0.11;
        const TzinavaResult movingGun = TzinavaReproducer::run(input, parameters);
        return result.faceThicknessMeters.size() == 2
            && approximatelyEqual(result.faceThicknessMeters[0], 1.0e-6, 1.0e-12)
            && approximatelyEqual(result.faceThicknessMeters[1], 1.0e-6, 1.0e-12)
            && approximatelyEqual(movingGun.faceThicknessMeters[0], 1.0e-6,
                1.0e-12)
            && approximatelyEqual(movingGun.faceThicknessMeters[1], 1.0e-6,
                1.0e-12);
    }

    bool testTzinavaSpraySwitchBoundary()
    {
        using namespace spraythickness::published;
        TzinavaInputModel input;
        input.initialMesh.vertices = {
            { -0.001, -0.001, 0.0 },
            { 0.001, -0.001, 0.0 },
            { 0.0, 0.001, 0.0 }
        };
        input.initialMesh.faces = { { 0, 1, 2 } };
        input.gunTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), false },
            { 2.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), false }
        };
        TzinavaParameters parameters;
        parameters.beamRadiusMeters = 0.05;
        parameters.gaussianRadialProfile = false;
        parameters.stationaryTimeStepSeconds = 2.0;
        parameters.lookupReferenceDwellSeconds = 1.0;
        parameters.thicknessTable.standOffDistancesMeters = { 0.09, 0.11 };
        parameters.thicknessTable.impactAnglesDegrees = { 80.0, 90.0 };
        parameters.thicknessTable.thicknessMeters = {
            1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6
        };
        const TzinavaResult result = TzinavaReproducer::run(input, parameters);
        return result.faceThicknessMeters.size() == 1
            && approximatelyEqual(result.faceThicknessMeters[0], 1.0e-6,
                1.0e-12);
    }

    bool testTzinavaDensePlateOcclusion()
    {
        using namespace spraythickness::published;
        TzinavaInputModel input;
        constexpr int cells = 8;
        constexpr double spacing = 0.001;
        for(double height : { 0.0, -0.015 }) {
            const std::uint32_t base =
                static_cast<std::uint32_t>(input.initialMesh.vertices.size());
            for(int y = 0; y <= cells; ++y) {
                for(int x = 0; x <= cells; ++x) {
                    input.initialMesh.vertices.emplace_back(
                        (x - cells / 2) * spacing,
                        (y - cells / 2) * spacing, height);
                }
            }
            for(int y = 0; y < cells; ++y) {
                for(int x = 0; x < cells; ++x) {
                    const std::uint32_t a = base + y * (cells + 1) + x;
                    input.initialMesh.faces.push_back(
                        { a, a + 1, a + cells + 2 });
                    input.initialMesh.faces.push_back(
                        { a, a + cells + 2, a + cells + 1 });
                }
            }
        }
        input.gunTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        TzinavaParameters parameters;
        parameters.beamRadiusMeters = 0.05;
        parameters.gaussianRadialProfile = false;
        parameters.stationaryTimeStepSeconds = 0.25;
        parameters.lookupReferenceDwellSeconds = 1.0;
        parameters.thicknessTable.standOffDistancesMeters = { 0.09, 0.13 };
        parameters.thicknessTable.impactAnglesDegrees = { 0.0, 90.0 };
        parameters.thicknessTable.thicknessMeters = {
            0.0, 1.0e-6, 0.0, 1.0e-6
        };

        const TzinavaResult result = TzinavaReproducer::run(input, parameters);
        const std::size_t topFaceCount = cells * cells * 2;
        for(std::size_t face = 0; face < topFaceCount; ++face) {
            std::swap(input.initialMesh.faces[face][1],
                input.initialMesh.faces[face][2]);
        }
        const TzinavaResult reversedTop = TzinavaReproducer::run(input, parameters);
        return result.faceThicknessMeters.size() == 2 * topFaceCount
            && std::all_of(result.faceThicknessMeters.begin(),
                result.faceThicknessMeters.begin() + topFaceCount,
                [](double thickness) {
                    return approximatelyEqual(thickness, 1.0e-6, 1.0e-12);
                })
            && std::all_of(result.faceThicknessMeters.begin() + topFaceCount,
                result.faceThicknessMeters.end(), [](double thickness) {
                    return approximatelyEqual(thickness, 0.0, 1.0e-12);
                })
            && std::all_of(reversedTop.faceThicknessMeters.begin(),
                reversedTop.faceThicknessMeters.end(), [](double thickness) {
                    return approximatelyEqual(thickness, 0.0, 1.0e-12);
                });
    }

    bool testTzinavaOrthographicTriangleShadow()
    {
        using namespace spraythickness::published;
        TzinavaInputModel input;
        input.initialMesh.vertices = {
            { 0.004, -0.001, 0.09 },
            { 0.006, -0.001, 0.09 },
            { 0.005, 0.001, 0.09 },
            { 0.004, -0.001, 0.08 },
            { 0.006, -0.001, 0.08 },
            { 0.005, 0.001, 0.08 }
        };
        input.initialMesh.faces = { { 0, 1, 2 }, { 3, 4, 5 } };
        input.gunTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        TzinavaParameters parameters;
        parameters.beamRadiusMeters = 0.01;
        parameters.gaussianRadialProfile = false;
        parameters.stationaryTimeStepSeconds = 0.25;
        parameters.lookupReferenceDwellSeconds = 1.0;
        parameters.thicknessTable.standOffDistancesMeters = { 0.005, 0.025 };
        parameters.thicknessTable.impactAnglesDegrees = { 0.0, 90.0 };
        parameters.thicknessTable.thicknessMeters = {
            0.0, 1.0e-6, 0.0, 1.0e-6
        };

        const TzinavaResult result = TzinavaReproducer::run(input, parameters);
        return result.faceThicknessMeters.size() == 2
            && approximatelyEqual(result.faceThicknessMeters[0], 1.0e-6, 1.0e-12)
            && approximatelyEqual(result.faceThicknessMeters[1], 0.0, 1.0e-12);
    }

    bool testTzinavaMovingScanAccumulation()
    {
        using namespace spraythickness::published;
        TzinavaInputModel input;
        input.initialMesh.vertices = {
            Eigen::Vector3d(-0.0003, -0.0003, 0.0),
            Eigen::Vector3d(0.0003, -0.0003, 0.0),
            Eigen::Vector3d(0.0, 0.0006, 0.0)
        };
        input.initialMesh.faces = { { 0, 1, 2 } };
        input.gunTrajectory = {
            { 0.0, Eigen::Vector3d(-0.03, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 6.0, Eigen::Vector3d(0.03, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };

        TzinavaParameters parameters;
        parameters.beamRadiusMeters = 0.01;
        parameters.gaussianSigmaMeters = 0.003;
        parameters.gaussianRadialProfile = true;
        parameters.referenceSpotSpeedMillimetersPerSecond = 10.0;
        parameters.timeStepOverlapFactor = 0.1;
        parameters.stationaryTimeStepSeconds = 0.1;
        parameters.lookupReferenceDwellSeconds = 1.0;
        parameters.thicknessTable.standOffDistancesMeters = { 0.09, 0.12 };
        parameters.thicknessTable.impactAnglesDegrees = { 0.0, 90.0 };
        parameters.thicknessTable.thicknessMeters = {
            0.0, 10.0e-6, 0.0, 10.0e-6
        };

        const TzinavaResult onePass = TzinavaReproducer::run(input, parameters);
        parameters.timeStepOverlapFactor = 0.05;
        const TzinavaResult finerPass = TzinavaReproducer::run(input, parameters);
        input.gunTrajectory.push_back({ 12.0,
            Eigen::Vector3d(-0.03, 0.0, 0.1),
            -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        const TzinavaResult twoPasses = TzinavaReproducer::run(input, parameters);
        input.gunTrajectory.pop_back();
        input.gunTrajectory.back().timeSeconds = 3.0;
        const TzinavaResult fasterPass = TzinavaReproducer::run(input, parameters);
        input.gunTrajectory.back().timeSeconds = 6.0;
        for(Eigen::Vector3d& vertex : input.initialMesh.vertices) {
            vertex.y() += parameters.gaussianSigmaMeters;
        }
        const TzinavaResult offsetPass = TzinavaReproducer::run(input, parameters);
        const double expectedSpeedRatio =
            (parameters.speedCoefficientB * 10.0 + parameters.speedCoefficientC)
            / (parameters.speedCoefficientB * 20.0 + parameters.speedCoefficientC);

        const bool passed = onePass.faceThicknessMeters.size() == 1
            && finerPass.faceThicknessMeters.size() == 1
            && twoPasses.faceThicknessMeters.size() == 1
            && fasterPass.faceThicknessMeters.size() == 1
            && offsetPass.faceThicknessMeters.size() == 1
            && approximatelyEqual(onePass.faceThicknessMeters[0], 10.0e-6,
                0.2e-6)
            && approximatelyEqual(finerPass.faceThicknessMeters[0], 10.0e-6,
                0.2e-6)
            && approximatelyEqual(twoPasses.faceThicknessMeters[0], 20.0e-6,
                0.4e-6)
            && approximatelyEqual(fasterPass.faceThicknessMeters[0],
                10.0e-6 * expectedSpeedRatio, 0.2e-6)
            && approximatelyEqual(offsetPass.faceThicknessMeters[0],
                10.0e-6 * std::exp(-0.5), 0.2e-6);
        if(!passed) {
            std::cerr << "Tzinava scan thickness (one/fine/two/fast/offset, um): ";
            for(const TzinavaResult* result : {
                    &onePass, &finerPass, &twoPasses, &fasterPass, &offsetPass }) {
                for(double thickness : result->faceThicknessMeters) {
                    std::cerr << thickness * 1.0e6 << ' ';
                }
                std::cerr << "; ";
            }
            std::cerr << '\n';
        }
        return passed;
    }

    bool testTzinavaContinuousReversalAndBackface()
    {
        using namespace spraythickness::published;
        TzinavaInputModel input;
        input.initialMesh.vertices = {
            { -0.0003, -0.0003, 0.0 },
            { 0.0003, -0.0003, 0.0 },
            { 0.0, 0.0006, 0.0 }
        };
        input.initialMesh.faces = { { 0, 1, 2 } };
        input.gunTrajectory = {
            { 0.0, { -0.001, 0.0, 0.1 }, -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d::Zero(), true },
            { 1.0, { 0.001, 0.0, 0.1 }, -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d::Zero(), true }
        };
        TzinavaParameters parameters;
        parameters.beamRadiusMeters = 0.01;
        parameters.gaussianRadialProfile = false;
        parameters.referenceSpotSpeedMillimetersPerSecond = 2.0;
        parameters.timeStepOverlapFactor = 0.1;
        parameters.stationaryTimeStepSeconds = 0.1;
        parameters.lookupReferenceDwellSeconds = 1.0;
        parameters.thicknessTable.standOffDistancesMeters = { 0.09, 0.11 };
        parameters.thicknessTable.impactAnglesDegrees = { 0.0, 90.0 };
        parameters.thicknessTable.thicknessMeters = {
            0.0, 10.0e-6, 0.0, 10.0e-6
        };

        const TzinavaResult onePass = TzinavaReproducer::run(input, parameters);
        input.gunTrajectory.push_back({ 2.0, { -0.001, 0.0, 0.1 },
            -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        const TzinavaResult reversed = TzinavaReproducer::run(input, parameters);
        input.initialMesh.faces = { { 0, 2, 1 } };
        const TzinavaResult backface = TzinavaReproducer::run(input, parameters);
        return approximatelyEqual(onePass.faceThicknessMeters[0], 10.0e-6, 1.0e-12)
            && approximatelyEqual(reversed.faceThicknessMeters[0],
                20.0e-6, 1.0e-12)
            && approximatelyEqual(backface.faceThicknessMeters[0], 0.0, 1.0e-12);
    }

    bool testTzinavaMovingPlateOcclusion()
    {
        using namespace spraythickness::published;
        TzinavaInputModel input;
        constexpr int cells = 4;
        for(int y = 0; y <= cells; ++y) {
            for(int x = 0; x <= cells; ++x) {
                input.initialMesh.vertices.emplace_back(
                    0.01 * (static_cast<double>(x) / cells - 0.5),
                    0.01 * (static_cast<double>(y) / cells - 0.5), 0.0);
            }
        }
        for(int y = 0; y < cells; ++y) {
            for(int x = 0; x < cells; ++x) {
                const std::uint32_t a = y * (cells + 1) + x;
                input.initialMesh.faces.push_back(
                    { a, a + 1, a + cells + 2 });
                input.initialMesh.faces.push_back(
                    { a, a + cells + 2, a + cells + 1 });
            }
        }
        const std::uint32_t lowerBase =
            static_cast<std::uint32_t>(input.initialMesh.vertices.size());
        input.initialMesh.vertices.emplace_back(-0.0003, -0.0003, -0.01);
        input.initialMesh.vertices.emplace_back(0.0003, -0.0003, -0.01);
        input.initialMesh.vertices.emplace_back(0.0, 0.0006, -0.01);
        input.initialMesh.faces.push_back(
            { lowerBase, lowerBase + 1, lowerBase + 2 });
        input.gunTrajectory = {
            { 0.0, Eigen::Vector3d(-0.03, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 6.0, Eigen::Vector3d(0.03, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };

        TzinavaParameters parameters;
        parameters.beamKind = TzinavaBeamKind::CylindricalConical;
        parameters.beamRadiusMeters = 0.002;
        parameters.coneHalfAngleRadians = std::atan(0.08);
        parameters.gaussianSigmaMeters = 0.003;
        parameters.referenceSpotSpeedMillimetersPerSecond = 10.0;
        parameters.timeStepOverlapFactor = 0.1;
        parameters.stationaryTimeStepSeconds = 0.1;
        parameters.lookupReferenceDwellSeconds = 1.0;
        parameters.thicknessTable.standOffDistancesMeters = { 0.09, 0.12 };
        parameters.thicknessTable.impactAnglesDegrees = { 0.0, 90.0 };
        parameters.thicknessTable.thicknessMeters = {
            0.0, 10.0e-6, 0.0, 10.0e-6
        };

        const TzinavaResult result = TzinavaReproducer::run(input, parameters);
        const bool passed = result.faceThicknessMeters.size() > cells * cells * 2
            && *std::max_element(result.faceThicknessMeters.begin(),
                result.faceThicknessMeters.end() - 1) > 5.0e-6
            && approximatelyEqual(result.faceThicknessMeters.back(),
                0.0, 1.0e-12);
        if(!passed) {
            std::cerr << "Tzinava moving occlusion: faces="
                << result.faceThicknessMeters.size() << ", top nonzero="
                << std::count_if(result.faceThicknessMeters.begin(),
                    result.faceThicknessMeters.end() - 1,
                    [](double thickness) { return thickness > 0.0; })
                << ", bottom um="
                << result.faceThicknessMeters.back() * 1.0e6 << '\n';
        }
        return passed;
    }

    bool testWuDynamicCylinderDeposition()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 2.0e-6;
        parameters.gaussianSigmaMeters = 0.01;
        parameters.maximumDeflectionRadians = 0.01;
        parameters.rayAngularStepRadians = 0.01;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0, 0.0 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.0, 0.0 };

        const WuResult result = WuReproducer::run(input, parameters);
        if(result.depositedCylinders.empty()) {
            return false;
        }
        const WuDepositedCylinder& center = result.depositedCylinders.front();
        return approximatelyEqual(center.baseCenter.norm(), 0.0, 1.0e-12)
            && approximatelyEqual(center.heightMeters, 2.0e-6, 1.0e-12)
            && center.growthDirection.isApprox(Eigen::Vector3d::UnitZ());
    }

    bool testWuDepositionConservesSampleDuration()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 2.0e-6;
        parameters.gaussianSigmaMeters = 0.01;
        parameters.maximumDeflectionRadians = 0.01;
        parameters.rayAngularStepRadians = 0.01;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0, 0.0 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.0, 0.0 };
        const auto pose = [](double time) {
            return SprayPose{ time, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true };
        };
        input.nozzleTrajectory = { pose(0.0), pose(1.0) };
        const WuResult coarse = WuReproducer::run(input, parameters);
        input.nozzleTrajectory = { pose(0.0), pose(0.5), pose(1.0) };
        const WuResult fine = WuReproducer::run(input, parameters);
        const auto nominalVolume = [](const WuResult& result) {
            double volume = 0.0;
            for(const WuDepositedCylinder& cylinder : result.depositedCylinders) {
                volume += cylinder.radiusMeters * cylinder.radiusMeters
                    * cylinder.heightMeters;
            }
            return volume;
        };
        const double coarseVolume = nominalVolume(coarse);
        const double fineVolume = nominalVolume(fine);
        const bool stationaryConserved = !coarse.depositedCylinders.empty()
            && fine.depositedCylinders.size()
                == 2 * coarse.depositedCylinders.size()
            && approximatelyEqual(coarse.depositedCylinders.front().heightMeters,
                2.0e-6, 1.0e-12)
            && approximatelyEqual(fine.depositedCylinders.front().heightMeters,
                1.0e-6, 1.0e-12)
            && std::abs(fineVolume - coarseVolume) <= coarseVolume * 0.005;
        if(!stationaryConserved) {
            return false;
        }

        parameters.referencePoseDurationSeconds = 0.01;
        const auto scan = [](double step) {
            std::vector<SprayPose> poses;
            const int intervalCount = static_cast<int>(std::round(0.04 / step));
            for(int index = 0; index <= intervalCount; ++index) {
                const double time = index * step;
                poses.push_back({ time,
                    Eigen::Vector3d(0.06 * time, 0.0, 0.03),
                    -Eigen::Vector3d::UnitZ(),
                    Eigen::Vector3d(0.06, 0.0, 0.0), true });
            }
            return poses;
        };
        input.nozzleTrajectory = scan(0.01);
        const WuResult coarseScan = WuReproducer::run(input, parameters);
        input.nozzleTrajectory = scan(0.005);
        const WuResult fineScan = WuReproducer::run(input, parameters);
        const double coarseScanVolume = nominalVolume(coarseScan);
        const bool scanConserved = coarseScanVolume > 0.0
            && fineScan.depositedCylinders.size()
                == 2 * coarseScan.depositedCylinders.size()
            && std::abs(nominalVolume(fineScan) - coarseScanVolume)
                <= coarseScanVolume * 0.02;
        if(!scanConserved) {
            return false;
        }

        input.nozzleTrajectory = { pose(0.0), pose(0.04), pose(0.055) };
        const WuResult uneven = WuReproducer::run(input, parameters);
        const std::size_t raysPerPose =
            uneven.statistics.visibilityQueryCount / 2;
        return uneven.depositedCylinders.size() == 2 * raysPerPose
            && approximatelyEqual(uneven.depositedCylinders.front().heightMeters,
                8.0e-6, 1.0e-12)
            && approximatelyEqual(
                uneven.depositedCylinders[raysPerPose].heightMeters,
                3.0e-6, 1.0e-12);
    }

    bool testWuPlateRejectsBackFace()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        input.generatedPlateStack = true;
        for(double time : { 0.0, 1.0 }) {
            input.nozzleTrajectory.push_back({ time,
                Eigen::Vector3d(0.0, 0.0, -0.03),
                Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        }
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 5.0e-6;
        parameters.gaussianSigmaMeters = 0.004;
        parameters.maximumDeflectionRadians = 0.02;
        parameters.rayAngularStepRadians = 0.01;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0, 0.0 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.0, 0.0 };
        return WuReproducer::run(input, parameters).depositedCylinders.empty();
    }

    bool testWuSidewallUsesSubstrateAngle()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        for(const auto& position : { Eigen::Vector3d(0.0, 0.0, 0.1),
                Eigen::Vector3d(0.01045, 0.0, 0.1),
                Eigen::Vector3d(0.01045, 0.0, 0.1) }) {
            input.nozzleTrajectory.push_back({
                static_cast<double>(input.nozzleTrajectory.size()), position,
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        }
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 0.001;
        parameters.gaussianSigmaMeters = 0.1;
        parameters.maximumDeflectionRadians = 0.1;
        parameters.rayAngularStepRadians = 0.1;
        parameters.cylinderRadiusMeters = 0.0005;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients = {
            101.0, -400.0 / 90.0, 600.0 / (90.0 * 90.0),
            -400.0 / (90.0 * 90.0 * 90.0),
            100.0 / (90.0 * 90.0 * 90.0 * 90.0) };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients = { 1.0, 0.0 };

        const WuResult result = WuReproducer::run(input, parameters);
        double largestSidewallIncrement = 0.0;
        for(const WuDepositedCylinder& cylinder : result.depositedCylinders) {
            if(cylinder.baseCenter.z() > 1.0e-5
                && std::abs(cylinder.baseCenter.x()) < 0.001
                && std::abs(cylinder.baseCenter.y()) < 0.001) {
                largestSidewallIncrement = std::max(largestSidewallIncrement,
                    cylinder.heightMeters);
            }
        }
        return largestSidewallIncrement > 0.0
            && largestSidewallIncrement < 0.01;
    }

    bool testWuInvalidSpeedCorrection()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.12),
                -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d(0.5, 0.0, 0.0), true },
            { 0.001, Eigen::Vector3d(0.0005, 0.0, 0.12),
                -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d(0.5, 0.0, 0.0), true }
        };
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 5.342534219921746e-6;
        parameters.gaussianSigmaMeters = 0.003589657150328816;
        parameters.maximumDeflectionRadians = 0.09795983193;
        parameters.rayAngularStepRadians = 0.09795983193;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 16.57062, -0.9661822, 0.020547, -0.0001842, 5.98125e-7 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { -0.683333, 0.208855, -0.007281, 0.000072 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.239, -0.0054 };

        bool rejected = false;
        try {
            WuReproducer::run(input, parameters);
        } catch(const std::out_of_range& error) {
            rejected = std::string(error.what()).find("500")
                    != std::string::npos
                && std::string(error.what()).find("-1.461")
                    != std::string::npos;
        }
        input.nozzleTrajectory[0].position.z() = 0.03;
        input.nozzleTrajectory[1].position =
            Eigen::Vector3d(0.00005, 0.0, 0.03);
        input.nozzleTrajectory[0].linearVelocity.x() = 0.05;
        input.nozzleTrajectory[1].linearVelocity.x() = 0.05;
        const WuResult valid = WuReproducer::run(input, parameters);
        if(!rejected || valid.depositedCylinders.empty()) {
            std::cerr << "Wu speed: rejected=" << rejected
                      << ", valid cylinders="
                      << valid.depositedCylinders.size() << '\n';
        }
        return rejected && !valid.depositedCylinders.empty();
    }

    bool testWuCylinderStackAndPlateOcclusion()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        const std::size_t lowerOffset = input.substrate.vertices.size();
        for(const Eigen::Vector3d& vertex : makePlaneMesh().vertices) {
            input.substrate.vertices.push_back(
                vertex - 0.015 * Eigen::Vector3d::UnitZ());
        }
        input.substrate.faces.push_back({
            static_cast<std::uint32_t>(lowerOffset),
            static_cast<std::uint32_t>(lowerOffset + 1),
            static_cast<std::uint32_t>(lowerOffset + 2) });
        input.substrate.faces.push_back({
            static_cast<std::uint32_t>(lowerOffset),
            static_cast<std::uint32_t>(lowerOffset + 2),
            static_cast<std::uint32_t>(lowerOffset + 3) });
        for(double time : { 0.0, 1.0, 2.0 }) {
            input.nozzleTrajectory.push_back({ time,
                Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        }
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 2.0e-6;
        parameters.gaussianSigmaMeters = 0.01;
        parameters.maximumDeflectionRadians = 0.08;
        parameters.rayAngularStepRadians = 0.005;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0, 0.0 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.0, 0.0 };

        const WuResult result = WuReproducer::run(input, parameters);
        if(result.depositedCylinders.size() <= 512) {
            return false;
        }
        const WuDepositedCylinder& first = result.depositedCylinders.front();
        const auto second = std::find_if(result.depositedCylinders.begin() + 1,
            result.depositedCylinders.end(), [](const auto& cylinder) {
                return cylinder.baseCenter.x() == 0.0
                    && cylinder.baseCenter.y() == 0.0
                    && cylinder.baseCenter.z() > 0.0;
            });
        return second != result.depositedCylinders.end()
            && approximatelyEqual(first.baseCenter.z(), 0.0)
            && approximatelyEqual(second->baseCenter.z(),
                first.heightMeters, 1.0e-12)
            && std::all_of(result.depositedCylinders.begin(),
                result.depositedCylinders.end(), [](const auto& cylinder) {
                    return cylinder.substrateFaceIndex < 2
                        && cylinder.baseCenter.z() >= -1.0e-12;
                });
    }

    bool testWuRecommendedPlateScan()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        const std::size_t lowerOffset = input.substrate.vertices.size();
        for(const Eigen::Vector3d& vertex : makePlaneMesh().vertices) {
            input.substrate.vertices.push_back(
                vertex - 0.02 * Eigen::Vector3d::UnitZ());
        }
        input.substrate.faces.push_back({
            static_cast<std::uint32_t>(lowerOffset),
            static_cast<std::uint32_t>(lowerOffset + 1),
            static_cast<std::uint32_t>(lowerOffset + 2) });
        input.substrate.faces.push_back({
            static_cast<std::uint32_t>(lowerOffset),
            static_cast<std::uint32_t>(lowerOffset + 2),
            static_cast<std::uint32_t>(lowerOffset + 3) });
        for(int index = 0; index < 5; ++index) {
            input.nozzleTrajectory.push_back({ 0.01 * index,
                Eigen::Vector3d((index - 2) * 0.0006, 0.0, 0.03),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        }
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 5.342534219921746e-6;
        parameters.gaussianSigmaMeters = 0.003589657150328816;
        parameters.maximumDeflectionRadians = 0.09795983193;
        parameters.rayAngularStepRadians = 0.0081633193275;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 16.57062, -0.9661822, 0.020547, -0.0001842, 5.98125e-7 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { -0.683333, 0.208855, -0.007281, 0.000072 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.239, -0.0054 };

        std::vector<std::string> stages;
        ReproductionExecution execution;
        execution.diagnostic = [&stages](const std::string& stage,
            const std::string& details) {
            stages.push_back(stage + ": " + details);
        };
        const WuResult result = WuReproducer::run(input, parameters, execution);
        const auto hasStage = [&stages](const std::string& name) {
            return std::any_of(stages.begin(), stages.end(),
                [&name](const std::string& value) {
                    return value.find(name) == 0;
                });
        };
        const bool passed = result.statistics.trajectorySampleCount == 5
            && result.statistics.visibilityQueryCount > 0
            && !result.depositedCylinders.empty()
            && std::all_of(result.depositedCylinders.begin(),
                result.depositedCylinders.end(), [](const auto& cylinder) {
                    return cylinder.substrateFaceIndex < 2
                        && cylinder.baseCenter.z() >= -1.0e-12;
                })
            && hasStage("Ray intersections")
            && hasStage("Deposition decisions");
        if(!passed) {
            std::cerr << "Wu plate scan: poses="
                      << result.statistics.trajectorySampleCount
                      << ", rays=" << result.statistics.visibilityQueryCount
                      << ", cylinders=" << result.depositedCylinders.size()
                      << ", intersection log=" << hasStage("Ray intersections")
                      << ", deposition log=" << hasStage("Deposition decisions")
                      << '\n';
        }
        return passed;
    }

    bool testWuUsesOriginalPoses()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        for(int index : { 0, 4 }) {
            input.nozzleTrajectory.push_back({ 0.025 * index,
                Eigen::Vector3d(-0.003 + 0.0015 * index, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        }
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 2.0e-6;
        parameters.gaussianSigmaMeters = 0.003;
        parameters.maximumDeflectionRadians = 0.02;
        parameters.rayAngularStepRadians = 0.02;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0, 0.0 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.0, 0.0 };

        const WuResult sparse = WuReproducer::run(input, parameters);
        input.nozzleTrajectory.clear();
        for(int index = 0; index <= 4; ++index) {
            input.nozzleTrajectory.push_back({ 0.025 * index,
                Eigen::Vector3d(-0.003 + 0.0015 * index, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        }
        const WuResult dense = WuReproducer::run(input, parameters);
        if(sparse.statistics.trajectorySampleCount != 2
            || dense.statistics.trajectorySampleCount != 5
            || sparse.statistics.visibilityQueryCount * 4
                != dense.statistics.visibilityQueryCount
            || sparse.depositedCylinders.empty()
            || dense.depositedCylinders.empty()) {
            std::cerr << "Wu original poses: sparse poses="
                      << sparse.statistics.trajectorySampleCount
                      << ", dense poses="
                      << dense.statistics.trajectorySampleCount
                      << ", sparse rays="
                      << sparse.statistics.visibilityQueryCount
                      << ", dense rays="
                      << dense.statistics.visibilityQueryCount << '\n';
            return false;
        }
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 0.025, Eigen::Vector3d(0.0015001, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        const WuResult justAbove = WuReproducer::run(input, parameters);
        return justAbove.statistics.trajectorySampleCount == 2;
    }

    bool testWuGeneratedPlateCoverageAndPoseBatch()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        input.generatedPlateStack = true;
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.03),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.03),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), false }
        };
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 2.0e-6;
        parameters.gaussianSigmaMeters = 0.003;
        parameters.maximumDeflectionRadians = 0.098;
        parameters.rayAngularStepRadians = 0.0081633193275;
        parameters.cylinderRadiusMeters = 2.0e-5;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0, 0.0 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.0, 0.0 };

        const WuResult first = WuReproducer::run(input, parameters);
        if(first.depositedCylinders.size() < 600
            || std::any_of(first.depositedCylinders.begin(),
                first.depositedCylinders.end(), [](const auto& cylinder) {
                    return std::abs(cylinder.baseCenter.z()) > 1.0e-10;
                })) {
            return false;
        }
        for(double x = -0.0029; x <= 0.0029; x += 0.00005) {
            for(double y = -0.0029; y <= 0.0029; y += 0.00005) {
                if(x * x + y * y > 0.0029 * 0.0029) {
                    continue;
                }
                const bool covered = std::any_of(first.depositedCylinders.begin(),
                    first.depositedCylinders.end(), [x, y](const auto& cylinder) {
                        const double dx = cylinder.baseCenter.x() - x;
                        const double dy = cylinder.baseCenter.y() - y;
                        return dx * dx + dy * dy
                            <= cylinder.radiusMeters * cylinder.radiusMeters;
                    });
                if(!covered) {
                    std::cerr << "Wu generated plate uncovered at "
                              << x << ", " << y << '\n';
                    return false;
                }
            }
        }
        input.nozzleTrajectory[1].sprayEnabled = true;
        input.nozzleTrajectory.push_back({ 2.0,
            Eigen::Vector3d(0.0, 0.0, 0.03),
            -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), false });
        const WuResult stacked = WuReproducer::run(input, parameters);
        double centerHeight = 0.0;
        double outerHeight = 0.0;
        for(const WuDepositedCylinder& cylinder : stacked.depositedCylinders) {
            const double radius = cylinder.baseCenter.head<2>().norm();
            const double top = (cylinder.baseCenter
                + cylinder.heightMeters * cylinder.growthDirection).z();
            if(radius < 0.00005) {
                centerHeight = std::max(centerHeight, top);
            } else if(radius > 0.002 && radius < 0.0029) {
                outerHeight = std::max(outerHeight, top);
            }
        }
        return stacked.statistics.trajectorySampleCount == 3
            && stacked.statistics.visibilityQueryCount
                == 2 * first.statistics.visibilityQueryCount
            && centerHeight > outerHeight
            && std::any_of(stacked.depositedCylinders.begin(),
                stacked.depositedCylinders.end(), [](const auto& cylinder) {
                    return cylinder.baseCenter.z() > 0.0;
                });
    }

    bool testWuSixDiscreteSprayPoses()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        input.generatedPlateStack = true;
        for(int index = 0; index < 6; ++index) {
            input.nozzleTrajectory.push_back({
                static_cast<double>(index),
                Eigen::Vector3d(-0.015 + 0.006 * index, 0.0, 0.03),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        }
        input.nozzleTrajectory.push_back({ 6.0,
            Eigen::Vector3d(0.015, 0.0, 0.03),
            -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), false });
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 2.0e-6;
        parameters.gaussianSigmaMeters = 0.003;
        parameters.maximumDeflectionRadians = 0.03;
        parameters.rayAngularStepRadians = 0.007;
        parameters.cylinderRadiusMeters = 2.0e-5;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0, 0.0 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.0, 0.0 };
        const WuResult result = WuReproducer::run(input, parameters);
        if(result.statistics.trajectorySampleCount != 7
            || result.depositedCylinders.empty()) {
            return false;
        }
        for(int index = 0; index < 6; ++index) {
            const double centerX = -0.015 + 0.006 * index;
            if(!std::any_of(result.depositedCylinders.begin(),
                    result.depositedCylinders.end(), [centerX](const auto& cylinder) {
                        return std::abs(cylinder.baseCenter.x() - centerX)
                            < 1.0e-8;
                    })) {
                return false;
            }
        }
        for(int index = 0; index < 5; ++index) {
            const double gapX = -0.012 + 0.006 * index;
            if(std::any_of(result.depositedCylinders.begin(),
                    result.depositedCylinders.end(), [gapX](const auto& cylinder) {
                        return std::abs(cylinder.baseCenter.x() - gapX)
                            < cylinder.radiusMeters;
                    })) {
                return false;
            }
        }
        return true;
    }

    bool testWuZeroDepositionDiagnostic()
    {
        using namespace spraythickness::published;
        WuInputModel input;
        input.substrate = makePlaneMesh();
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 2.0e-6;
        parameters.gaussianSigmaMeters = 0.01;
        parameters.maximumDeflectionRadians = 0.02;
        parameters.rayAngularStepRadians = 0.02;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 0.0, 0.0, 0.0, 0.0, 0.0 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { 1.0, 0.0, 0.0, 0.0 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.0, 0.0 };

        std::string intersections;
        std::string decisions;
        ReproductionExecution execution;
        execution.diagnostic = [&intersections, &decisions](
            const std::string& stage, const std::string& details) {
            if(stage == "Ray intersections") {
                intersections = details;
            } else if(stage == "Deposition decisions") {
                decisions = details;
            }
        };
        const WuResult result = WuReproducer::run(input, parameters, execution);
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients[0] = -1.0;
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients[0] = -1.0;
        const WuResult invalidFactors = WuReproducer::run(input, parameters);
        return result.depositedCylinders.empty()
            && invalidFactors.depositedCylinders.empty()
            && intersections.find("substrate=8") != std::string::npos
            && intersections.find("missed=0") != std::string::npos
            && decisions.find("angle=8") != std::string::npos;
    }

    bool testVanerioDynamicSurfaceUpdate()
    {
        using namespace spraythickness::published;
        VanerioInputModel input;
        input.initialStlSurface.vertices = {
            Eigen::Vector3d(-0.005, -0.005, 0.0),
            Eigen::Vector3d(0.005, -0.005, 0.0),
            Eigen::Vector3d(0.0, 0.01, 0.0)
        };
        input.initialStlSurface.faces = { { 0, 1, 2 } };
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        VanerioParameters parameters;
        parameters.growthRateCoefficientMetersPerSecond = 3.0e-6;
        parameters.jetRadiusMeters = 0.05;
        parameters.jetShapeCoefficientK2 = 2.0;
        parameters.maximumMeshEdgeMeters = 0.1;
        parameters.shadowGridStepMeters = 0.001;
        parameters.depositionEfficiencyByTangentAngle = {
            { 0.0, 1.0 }, { 1.0, 1.0 }
        };
        parameters.depositionEfficiencyByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };
        parameters.profileStretchByDistance = {
            { 0.05, 0.15 }, { 0.5, 0.5 }
        };

        const VanerioResult result = VanerioReproducer::run(input, parameters);
        if(result.evolvedStlSurface.vertices.size() != 3) {
            return false;
        }
        const auto integral = [](double stretch) {
            const double k2 = 2.0;
            const double limit = std::min(1.0, 1.0 / stretch);
            return (-std::expm1(-k2 * stretch * stretch * limit * limit)
                    / (2.0 * k2 * stretch * stretch)
                    - 0.5 * std::exp(-k2) * limit * limit)
                / (-std::expm1(-k2));
        };
        const double expected = 3.0e-6 * integral(1.0) / integral(0.5);
        if(result.vertexThicknessMeters.size() != result.evolvedStlSurface.vertices.size()) {
            return false;
        }
        for(std::size_t index = 0; index < result.evolvedStlSurface.vertices.size();
            ++index) {
            if(!approximatelyEqual(result.evolvedStlSurface.vertices[index].z(),
                    expected, 1.0e-12)
                || !approximatelyEqual(result.vertexThicknessMeters[index],
                    expected, 1.0e-12)) {
                return false;
            }
        }
        return true;
    }

    bool testVanerioDensePlateDeposition()
    {
        using namespace spraythickness::published;
        VanerioInputModel input;
        constexpr int cells = 20;
        constexpr double spacing = 0.001;
        for(double height : { 0.0, -0.015 }) {
            const std::uint32_t base = static_cast<std::uint32_t>(
                input.initialStlSurface.vertices.size());
            for(int y = 0; y <= cells; ++y) {
                for(int x = 0; x <= cells; ++x) {
                    input.initialStlSurface.vertices.emplace_back(
                        (x - cells / 2) * spacing,
                        (y - cells / 2) * spacing, height);
                }
            }
            for(int y = 0; y < cells; ++y) {
                for(int x = 0; x < cells; ++x) {
                    const std::uint32_t a = base + y * (cells + 1) + x;
                    input.initialStlSurface.faces.push_back(
                        { a, a + 1, a + cells + 1 });
                    input.initialStlSurface.faces.push_back(
                        { a + 1, a + cells + 2, a + cells + 1 });
                }
            }
        }
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.12),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 0.001, Eigen::Vector3d(0.0, 0.0, 0.12),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        VanerioParameters parameters;
        parameters.growthRateCoefficientMetersPerSecond = 0.000148831322775;
        parameters.jetRadiusMeters = 0.010768971450986448;
        parameters.jetShapeCoefficientK2 = 4.36;
        parameters.maximumMeshEdgeMeters = 0.001;
        parameters.shadowGridStepMeters = 0.0005;
        parameters.depositionEfficiencyByTangentAngle = {
            { 0.0, 1.0 }, { 1.0, 1.0 }
        };
        parameters.depositionEfficiencyByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };
        parameters.profileStretchByDistance = {
            { 0.05, 0.12, 0.15 }, { 2.4, 1.0, 0.8 }
        };
        const VanerioResult result = VanerioReproducer::run(input, parameters);
        std::size_t raisedTop = 0;
        std::size_t raisedBottom = 0;
        for(const Eigen::Vector3d& vertex : result.evolvedStlSurface.vertices) {
            raisedTop += vertex.z() > 1.0e-12 ? 1 : 0;
            raisedBottom += vertex.z() > -0.015 + 1.0e-12
                && vertex.z() < -0.01 ? 1 : 0;
        }
        std::cout << "Vanerio dense plate: contributions="
                  << result.statistics.candidatePairCount
                  << ", raised top=" << raisedTop
                  << ", raised bottom=" << raisedBottom << '\n';
        return result.statistics.candidatePairCount > 0
            && raisedTop > 0 && raisedBottom == 0;
    }

    bool testVanerioNearestOutletVertex()
    {
        using namespace spraythickness::published;
        VanerioInputModel input;
        input.initialStlSurface.vertices = {
            { 0.0029, 0.0029, 0.02 }, { 0.0031, 0.0029, 0.02 },
            { 0.0029, 0.0031, 0.02 }, { 0.0005, 0.0005, 0.01995 },
            { 0.0007, 0.0005, 0.01995 }, { 0.0005, 0.0007, 0.01995 }
        };
        input.initialStlSurface.faces = { { 0, 1, 2 }, { 3, 4, 5 } };
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.12),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.12),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        VanerioParameters parameters;
        parameters.growthRateCoefficientMetersPerSecond = 3.0e-6;
        parameters.jetRadiusMeters = 0.05;
        parameters.jetShapeCoefficientK2 = 2.0;
        parameters.maximumMeshEdgeMeters = 0.1;
        parameters.shadowGridStepMeters = 0.01;
        parameters.depositionEfficiencyByTangentAngle = {
            { 0.0, 1.0 }, { 1.0, 1.0 }
        };
        parameters.depositionEfficiencyByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };
        parameters.profileStretchByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };

        const VanerioResult result = VanerioReproducer::run(input, parameters);
        const bool passed = result.evolvedStlSurface.vertices.size() == 6
            && approximatelyEqual(result.evolvedStlSurface.vertices[0].z(),
                0.02, 1.0e-12)
            && result.evolvedStlSurface.vertices[3].z() > 0.01995;
        if(!passed) {
            std::cerr << "Vanerio nearest outlet vertex: top="
                      << result.evolvedStlSurface.vertices[0].z()
                      << ", lower="
                      << result.evolvedStlSurface.vertices[3].z() << '\n';
        }
        return passed;
    }

    bool testVanerioCylinderAcrossGridBoundary()
    {
        using namespace spraythickness::published;
        VanerioInputModel input;
        input.initialStlSurface.vertices = {
            { -0.002, 0.0, 0.02 }, { -0.001, 0.0, 0.02 },
            { -0.0015, 0.001, 0.02 },
            { 0.001, 0.0, 0.0 }, { 0.002, 0.0, 0.0 },
            { 0.0015, 0.001, 0.0 }
        };
        input.initialStlSurface.faces = { { 0, 1, 2 }, { 3, 4, 5 } };
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.12),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.12),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        VanerioParameters parameters;
        parameters.growthRateCoefficientMetersPerSecond = 3.0e-6;
        parameters.jetRadiusMeters = 0.05;
        parameters.jetShapeCoefficientK2 = 2.0;
        parameters.maximumMeshEdgeMeters = 0.1;
        parameters.shadowGridStepMeters = 0.01;
        parameters.depositionEfficiencyByTangentAngle = {
            { 0.0, 1.0 }, { 1.0, 1.0 }
        };
        parameters.depositionEfficiencyByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };
        parameters.profileStretchByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };

        const VanerioResult result = VanerioReproducer::run(input, parameters);
        return result.evolvedStlSurface.vertices[0].z() > 0.02
            && approximatelyEqual(result.evolvedStlSurface.vertices[3].z(),
                0.0, 1.0e-12);
    }

    bool testVanerioDistanceProfileSupport()
    {
        using namespace spraythickness::published;
        VanerioInputModel input;
        input.initialStlSurface.vertices = {
            { 0.0145, -0.0005, 0.0 },
            { 0.0155, -0.0005, 0.0 },
            { 0.015, 0.0005, 0.0 }
        };
        input.initialStlSurface.faces = { { 0, 1, 2 } };
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        VanerioParameters parameters;
        parameters.growthRateCoefficientMetersPerSecond = 3.0e-6;
        parameters.jetRadiusMeters = 0.01;
        parameters.jetShapeCoefficientK2 = 2.0;
        parameters.maximumMeshEdgeMeters = 0.1;
        parameters.shadowGridStepMeters = 0.001;
        parameters.depositionEfficiencyByTangentAngle = {
            { 0.0, 1.0 }, { 1.0, 1.0 }
        };
        parameters.depositionEfficiencyByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };
        parameters.profileStretchByDistance = {
            { 0.05, 0.15 }, { 0.5, 0.5 }
        };

        const VanerioResult result = VanerioReproducer::run(input, parameters);
        return result.statistics.candidatePairCount == 0
            && approximatelyEqual(result.evolvedStlSurface.vertices[0].z(),
                0.0, 1.0e-12);
    }

    bool testVanerioInterPassRemeshing()
    {
        using namespace spraythickness::published;
        VanerioInputModel input;
        input.initialStlSurface.vertices = {
            { 0.0, 0.0, 0.0 }, { 0.005, 0.0, 0.0 },
            { 0.005, 0.005, 0.0 }, { 0.0, 0.005, 0.0 }
        };
        input.initialStlSurface.faces = { { 0, 1, 3 }, { 1, 2, 3 } };
        input.nozzleTrajectory = {
            { 0.0, { 0.001, 0.001, 0.1 }, -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d::Zero(), true },
            { 1.0, { 0.0015, 0.001, 0.1 }, -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d::Zero(), true },
            { 1.001, { 0.002, 0.001, 0.1 }, -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d::Zero(), false }
        };
        VanerioParameters parameters;
        parameters.growthRateCoefficientMetersPerSecond = 0.02;
        parameters.jetRadiusMeters = 0.003;
        parameters.jetShapeCoefficientK2 = 2.0;
        parameters.maximumMeshEdgeMeters = 0.008;
        parameters.shadowGridStepMeters = 0.001;
        parameters.depositionEfficiencyByTangentAngle = {
            { 0.0, 10.0 }, { 1.0, 1.0 }
        };
        parameters.depositionEfficiencyByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };
        parameters.profileStretchByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };

        const VanerioResult singlePass =
            VanerioReproducer::run(input, parameters);
        input.nozzleTrajectory[2].position.x() = 0.001;
        const VanerioResult reversingPass =
            VanerioReproducer::run(input, parameters);
        return singlePass.evolvedStlSurface.faces.size() == 2
            && reversingPass.evolvedStlSurface.faces.size() > 2
            && reversingPass.vertexThicknessMeters.size()
                == reversingPass.evolvedStlSurface.vertices.size();
    }

    bool testVanerioSharedEdgeRemeshing()
    {
        using namespace spraythickness::published;
        VanerioInputModel input;
        input.initialStlSurface.vertices = {
            { 0.0, 0.0, 0.0 }, { 0.01, 0.0, 0.0 },
            { 0.01, 0.01, 0.0 }, { 0.0, 0.01, 0.0 }
        };
        input.initialStlSurface.faces = { { 0, 1, 3 }, { 1, 2, 3 } };
        input.nozzleTrajectory = {
            { 0.0, { 0.0, 0.0, 0.1 }, -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d::Zero(), false },
            { 1.0, { 0.0, 0.0, 0.1 }, -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d::Zero(), false }
        };
        VanerioParameters parameters;
        parameters.growthRateCoefficientMetersPerSecond = 1.0e-6;
        parameters.jetRadiusMeters = 0.02;
        parameters.jetShapeCoefficientK2 = 2.0;
        parameters.maximumMeshEdgeMeters = 0.006;
        parameters.shadowGridStepMeters = 0.001;
        parameters.depositionEfficiencyByTangentAngle = {
            { 0.0, 1.0 }, { 1.0, 1.0 }
        };
        parameters.depositionEfficiencyByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };
        parameters.profileStretchByDistance = {
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };
        const VanerioResult result = VanerioReproducer::run(input, parameters);
        for(std::size_t first = 0;
            first < result.evolvedStlSurface.vertices.size(); ++first) {
            for(std::size_t second = first + 1;
                second < result.evolvedStlSurface.vertices.size(); ++second) {
                if((result.evolvedStlSurface.vertices[first]
                        - result.evolvedStlSurface.vertices[second]).norm()
                    < 1.0e-12) {
                    return false;
                }
            }
        }
        return result.evolvedStlSurface.faces.size() > 2;
    }

    double flatCylinderUnionVolume(
        const std::vector<spraythickness::published::DynamicDepositedCylinder>&
            cylinders)
    {
        constexpr int kCells = 800;
        constexpr double kPitch = 0.00005;
        constexpr double kOrigin = -0.02;
        std::vector<double> heights(kCells * kCells, 0.0);
        for(const auto& cylinder : cylinders) {
            const int firstX = std::max(0, static_cast<int>(std::floor(
                (cylinder.baseCenter.x() - cylinder.radiusMeters - kOrigin)
                    / kPitch)));
            const int lastX = std::min(kCells - 1,
                static_cast<int>(std::ceil(
                    (cylinder.baseCenter.x() + cylinder.radiusMeters - kOrigin)
                        / kPitch)));
            const int firstY = std::max(0, static_cast<int>(std::floor(
                (cylinder.baseCenter.y() - cylinder.radiusMeters - kOrigin)
                    / kPitch)));
            const int lastY = std::min(kCells - 1,
                static_cast<int>(std::ceil(
                    (cylinder.baseCenter.y() + cylinder.radiusMeters - kOrigin)
                        / kPitch)));
            for(int y = firstY; y <= lastY; ++y) {
                for(int x = firstX; x <= lastX; ++x) {
                    const double dx = kOrigin + (x + 0.5) * kPitch
                        - cylinder.baseCenter.x();
                    const double dy = kOrigin + (y + 0.5) * kPitch
                        - cylinder.baseCenter.y();
                    if(dx * dx + dy * dy
                        <= cylinder.radiusMeters * cylinder.radiusMeters) {
                        const std::size_t index = y * kCells + x;
                        heights[index] = std::max(heights[index],
                            cylinder.heightMeters);
                    }
                }
            }
        }
        double volume = 0.0;
        for(const double height : heights) {
            volume += height * kPitch * kPitch;
        }
        return volume;
    }

    bool testDynamicSurfaceBatchUpdate()
    {
        using namespace spraythickness::published;
        DynamicSurfaceInputModel input;
        input.initialStlSurface = makePlaneMesh();
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        DynamicSurfaceParameters parameters;
        parameters.equivalentParticleDiameterMeters = 100.0e-6;
        parameters.depositedCylinderRadiusMeters = 0.3e-3;
        parameters.equivalentParticleRatePerSecond = 600.0;
        parameters.angularMeanRadians = 0.0;
        parameters.angularSigmaRadians = 0.03;
        parameters.maximumEjectionAngleRadians = 0.06;
        parameters.polarBinCount = 12;
        parameters.depositionBatchDurationSeconds = 1.0;
        parameters.randomSeed = 7;
        parameters.relativeBuildUpByInclinationRadians = {
            { 0.0, 1.5707963267948966 }, { 1.0, 1.0 }
        };
        parameters.voxelLeafMeters = 0.15e-3;
        parameters.uniformSamplingRadiusMeters = 0.20e-3;
        parameters.normalSearchRadiusMeters = 2.0e-3;
        parameters.normalMaximumNeighbors = 20;
        parameters.boundaryGradientThresholdRadiansPerMeter = 1000.0;
        parameters.dbscanRadiusMeters = 1.5e-3;
        parameters.dbscanMinimumPoints = 3;
        parameters.alphaShapeRadiusMeters = 2.0e-3;
        parameters.poissonDepth = 5;
        parameters.maximumHoleBoundaryEdges = 20;
        parameters.smoothingIterations = 1;
        parameters.smoothingRelaxation = 0.1;

        std::vector<std::string> stages;
        ReproductionExecution execution;
        execution.diagnostic = [&stages](const std::string& stage,
                                   const std::string&) {
            stages.push_back(stage);
        };
        const DynamicSurfaceResult result =
            DynamicSurfaceReproducer::run(input, parameters, execution);
        const double particleVolume = 3.14159265358979323846 / 6.0
            * std::pow(parameters.equivalentParticleDiameterMeters, 3.0);
        const double targetVolume = result.depositedCylinders.size()
            * particleVolume;
        const double volumeError = std::abs(
            flatCylinderUnionVolume(result.depositedCylinders)
                - targetVolume) / targetVolume;
        const bool stagesReported = std::find(stages.begin(), stages.end(),
            "Point-cloud reconstruction") != stages.end()
            && std::find(stages.begin(), stages.end(), "Result assembly")
                != stages.end();
        if(!stagesReported) {
            return false;
        }

        parameters.equivalentParticleRatePerSecond = 20.0;
        parameters.voxelLeafMeters = 1.0;
        parameters.uniformSamplingRadiusMeters = 1.0;
        bool failureContainsCounts = false;
        try {
            DynamicSurfaceReproducer::run(input, parameters);
        } catch(const std::runtime_error& error) {
            const std::string message = error.what();
            failureContainsCounts = message.find("pose 1, batch 1")
                    != std::string::npos
                && message.find("local-density filtering retained")
                    != std::string::npos
                && message.find("deposited cylinders") != std::string::npos;
        }
        return failureContainsCounts
            && !result.depositedCylinders.empty()
            && volumeError < 0.05
            && !result.evolvedStlSurface.empty()
            && result.evolvedStlSurface.faces.size() > 2;
    }

    bool testDynamicSurfacePreservesOccludedSubstrate()
    {
        using namespace spraythickness::published;
        DynamicSurfaceInputModel input;
        input.initialStlSurface.vertices = {
            { -0.005, -0.005, 0.0 }, { 0.005, -0.005, 0.0 },
            { 0.005, 0.005, 0.0 }, { -0.005, 0.005, 0.0 },
            { -0.005, -0.005, -0.02 }, { 0.005, -0.005, -0.02 },
            { 0.005, 0.005, -0.02 }, { -0.005, 0.005, -0.02 }
        };
        input.initialStlSurface.faces = {
            { 0, 1, 2 }, { 0, 2, 3 }, { 4, 5, 6 }, { 4, 6, 7 }
        };
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 2.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        DynamicSurfaceParameters parameters;
        parameters.equivalentParticleDiameterMeters = 100.0e-6;
        parameters.depositedCylinderRadiusMeters = 0.3e-3;
        parameters.equivalentParticleRatePerSecond = 600.0;
        parameters.angularMeanRadians = 0.0;
        parameters.angularSigmaRadians = 0.03;
        parameters.maximumEjectionAngleRadians = 0.06;
        parameters.polarBinCount = 12;
        parameters.depositionBatchDurationSeconds = 1.0;
        parameters.randomSeed = 7;
        parameters.relativeBuildUpByInclinationRadians = {
            { 0.0, 1.5707963267948966 }, { 1.0, 1.0 }
        };
        parameters.voxelLeafMeters = 0.15e-3;
        parameters.uniformSamplingRadiusMeters = 0.20e-3;
        parameters.normalSearchRadiusMeters = 2.0e-3;
        parameters.normalMaximumNeighbors = 20;
        parameters.boundaryGradientThresholdRadiansPerMeter = 1000.0;
        parameters.dbscanRadiusMeters = 1.5e-3;
        parameters.dbscanMinimumPoints = 3;
        parameters.alphaShapeRadiusMeters = 2.0e-3;
        parameters.poissonDepth = 5;
        parameters.maximumHoleBoundaryEdges = 20;
        parameters.smoothingIterations = 1;
        parameters.smoothingRelaxation = 0.1;

        const DynamicSurfaceResult result =
            DynamicSurfaceReproducer::run(input, parameters);
        std::size_t lowerFaces = 0;
        for(const auto& face : result.evolvedStlSurface.faces) {
            if(std::all_of(face.begin(), face.end(), [&](std::uint32_t index) {
                    return std::abs(result.evolvedStlSurface.vertices[index].z()
                        + 0.02) < 1.0e-8;
                })) {
                ++lowerFaces;
            }
        }
        const bool lowerDeposited = std::any_of(
            result.depositedCylinders.begin(), result.depositedCylinders.end(),
            [](const DynamicDepositedCylinder& cylinder) {
                return cylinder.baseCenter.z() < -0.01;
            });
        const bool coatingHit = std::any_of(
            result.depositedCylinders.begin(), result.depositedCylinders.end(),
            [](const DynamicDepositedCylinder& cylinder) {
                return cylinder.baseCenter.z() > 1.0e-9;
            });
        if(lowerFaces != 2 || lowerDeposited) {
            std::cerr << "Dynamic-surface occluded substrate: lower faces="
                      << lowerFaces << ", lower deposited=" << lowerDeposited
                      << '\n';
        }
        return lowerFaces == 2 && !lowerDeposited && coatingHit
            && result.depositedCylinders.size() > 600;
    }

    bool testDynamicSurfaceDenseTrajectory()
    {
        using namespace spraythickness::published;
        DynamicSurfaceInputModel input;
        input.initialStlSurface = makePlaneMesh();
        for(std::size_t index = 0; index <= 1000; ++index) {
            input.nozzleTrajectory.push_back({
                static_cast<double>(index) * 0.001,
                Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        }
        DynamicSurfaceParameters parameters;
        parameters.equivalentParticleDiameterMeters = 100.0e-6;
        parameters.depositedCylinderRadiusMeters = 0.3e-3;
        parameters.equivalentParticleRatePerSecond = 600.0;
        parameters.angularMeanRadians = 0.0;
        parameters.angularSigmaRadians = 0.03;
        parameters.maximumEjectionAngleRadians = 0.06;
        parameters.polarBinCount = 12;
        parameters.depositionBatchDurationSeconds = 1.0;
        parameters.randomSeed = 7;
        parameters.relativeBuildUpByInclinationRadians = {
            { 0.0, 1.5707963267948966 }, { 1.0, 1.0 }
        };
        parameters.voxelLeafMeters = 0.15e-3;
        parameters.uniformSamplingRadiusMeters = 0.20e-3;
        parameters.normalSearchRadiusMeters = 2.0e-3;
        parameters.normalMaximumNeighbors = 20;
        parameters.boundaryGradientThresholdRadiansPerMeter = 1000.0;
        parameters.dbscanRadiusMeters = 1.5e-3;
        parameters.dbscanMinimumPoints = 3;
        parameters.alphaShapeRadiusMeters = 2.0e-3;
        parameters.poissonDepth = 5;
        parameters.maximumHoleBoundaryEdges = 20;
        parameters.smoothingIterations = 1;
        parameters.smoothingRelaxation = 0.1;

        const DynamicSurfaceInputModel denseInput = input;
        const DynamicSurfaceResult dense =
            DynamicSurfaceReproducer::run(input, parameters);
        input.nozzleTrajectory.erase(input.nozzleTrajectory.begin() + 1,
            input.nozzleTrajectory.end() - 1);
        const DynamicSurfaceResult sparse =
            DynamicSurfaceReproducer::run(input, parameters);
        input.nozzleTrajectory.push_back({ 1.002,
            Eigen::Vector3d(0.0, 0.0, 0.1),
            -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        const DynamicSurfaceResult shortTail =
            DynamicSurfaceReproducer::run(input, parameters);
        parameters.normalSearchRadiusMeters = 1.0e-3;
        std::vector<std::string> sparseFringeDetails;
        ReproductionExecution sparseFringeExecution;
        sparseFringeExecution.diagnostic = [&sparseFringeDetails](
                                               const std::string&,
                                               const std::string& details) {
            sparseFringeDetails.push_back(details);
        };
        const DynamicSurfaceResult sparseFringe =
            DynamicSurfaceReproducer::run(
                denseInput, parameters, sparseFringeExecution);
        const bool sparsePointsRemoved = std::any_of(
            sparseFringeDetails.begin(), sparseFringeDetails.end(),
            [](const std::string& details) {
                return details.find("sparse removed=") != std::string::npos
                    && details.find("sparse removed=0") == std::string::npos;
            });
        return !dense.evolvedStlSurface.empty()
            && dense.statistics.visibilityQueryCount == 600
            && dense.statistics.visibilityQueryCount
                == sparse.statistics.visibilityQueryCount
            && dense.depositedCylinders.size()
                == sparse.depositedCylinders.size()
            && shortTail.statistics.visibilityQueryCount == 601
            && shortTail.depositedCylinders.size()
                == sparse.depositedCylinders.size() + 1
            && !sparseFringe.evolvedStlSurface.empty()
            && sparsePointsRemoved;
    }

    bool testDynamicSurfaceSmallRegion()
    {
        using namespace spraythickness::published;
        DynamicSurfaceInputModel input;
        input.initialStlSurface = makePlaneMesh();
        input.nozzleTrajectory = {
            { 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.0, Eigen::Vector3d(0.014, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true },
            { 1.015, Eigen::Vector3d(0.014, 0.0, 0.1),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true }
        };
        DynamicSurfaceParameters parameters;
        parameters.equivalentParticleDiameterMeters = 100.0e-6;
        parameters.depositedCylinderRadiusMeters = 0.3e-3;
        parameters.equivalentParticleRatePerSecond = 600.0;
        parameters.angularMeanRadians = 0.0;
        parameters.angularSigmaRadians = 0.005;
        parameters.maximumEjectionAngleRadians = 0.015;
        parameters.polarBinCount = 12;
        parameters.depositionBatchDurationSeconds = 1.1;
        parameters.randomSeed = 7;
        parameters.relativeBuildUpByInclinationRadians = {
            { 0.0, 1.5707963267948966 }, { 1.0, 1.0 }
        };
        parameters.voxelLeafMeters = 0.15e-3;
        parameters.uniformSamplingRadiusMeters = 0.20e-3;
        parameters.normalSearchRadiusMeters = 2.0e-3;
        parameters.normalMaximumNeighbors = 20;
        parameters.boundaryGradientThresholdRadiansPerMeter = 1000.0;
        parameters.dbscanRadiusMeters = 1.5e-3;
        parameters.dbscanMinimumPoints = 3;
        parameters.alphaShapeRadiusMeters = 2.0e-3;
        parameters.poissonDepth = 5;
        parameters.maximumHoleBoundaryEdges = 20;
        parameters.smoothingIterations = 1;
        parameters.smoothingRelaxation = 0.1;

        bool skippedSmallRegion = false;
        ReproductionExecution execution;
        execution.diagnostic = [&skippedSmallRegion](
                                   const std::string&, const std::string& details) {
            skippedSmallRegion |= details.find("DBSCAN skipped ")
                != std::string::npos;
        };
        const DynamicSurfaceResult result =
            DynamicSurfaceReproducer::run(input, parameters, execution);
        return skippedSmallRegion && result.evolvedStlSurface.faces.size() > 2
            && result.depositedCylinders.size() > 10;
    }

    bool testDynamicSurfaceGeneratedPlateScan()
    {
        using namespace spraythickness::published;
        DynamicSurfaceInputModel input;
        input.initialStlSurface.vertices = {
            { -0.05, -0.05, 0.0 }, { 0.05, -0.05, 0.0 },
            { 0.05, 0.05, 0.0 }, { -0.05, 0.05, 0.0 },
            { -0.05, -0.05, -0.02 }, { 0.05, -0.05, -0.02 },
            { 0.05, 0.05, -0.02 }, { -0.05, 0.05, -0.02 }
        };
        input.initialStlSurface.faces = {
            { 0, 1, 2 }, { 0, 2, 3 }, { 4, 5, 6 }, { 4, 6, 7 }
        };
        for(std::size_t index = 0; index <= 560; ++index) {
            const double progress = index <= 280
                ? static_cast<double>(index) / 280.0
                : static_cast<double>(560 - index) / 280.0;
            input.nozzleTrajectory.push_back({
                static_cast<double>(index) * 0.001,
                Eigen::Vector3d(-0.07 + 0.14 * progress, 0.0, 0.12),
                -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(),
                index < 560 });
        }
        DynamicSurfaceParameters parameters;
        parameters.equivalentParticleDiameterMeters = 0.0002;
        parameters.depositedCylinderRadiusMeters = 0.0003;
        parameters.equivalentParticleRatePerSecond = 2890.3959983420814;
        parameters.angularMeanRadians = 0.0;
        parameters.angularSigmaRadians = 0.030385709702379613;
        parameters.maximumEjectionAngleRadians = 0.09795983192999999;
        parameters.polarBinCount = 25;
        parameters.depositionBatchDurationSeconds = 0.5;
        parameters.randomSeed = 2026;
        parameters.relativeBuildUpQuadraticPerDegreeSquared = 0.0001482;
        parameters.relativeBuildUpByInclinationRadians = {
            { 0.0, 0.6981317007977318 }, { 1.0, 0.76288 }
        };
        parameters.voxelLeafMeters = 0.00015;
        parameters.uniformSamplingRadiusMeters = 0.0002;
        parameters.normalSearchRadiusMeters = 0.002;
        parameters.normalMaximumNeighbors = 20;
        parameters.boundaryGradientThresholdRadiansPerMeter = 1000.0;
        parameters.dbscanRadiusMeters = 0.0015;
        parameters.dbscanMinimumPoints = 3;
        parameters.alphaShapeRadiusMeters = 0.002;
        parameters.poissonDepth = 5;
        parameters.maximumHoleBoundaryEdges = 20;
        parameters.smoothingIterations = 1;
        parameters.smoothingRelaxation = 0.1;

        try {
            const DynamicSurfaceResult result =
                DynamicSurfaceReproducer::run(input, parameters);
            std::size_t lowerHits = 0;
            std::size_t interiorLowerHits = 0;
            for(const auto& cylinder : result.depositedCylinders) {
                if(cylinder.baseCenter.z() < -0.01) {
                    ++lowerHits;
                    if(std::abs(cylinder.baseCenter.x()) < 0.049
                        && std::abs(cylinder.baseCenter.y()) < 0.049) {
                        ++interiorLowerHits;
                    }
                }
            }
            if(interiorLowerHits > 0) {
                std::cerr << "Dynamic-surface lower hits: " << lowerHits
                          << ", inside upper plate: " << interiorLowerHits
                          << '\n';
            }
            return result.depositedCylinders.size() > 100
                && result.evolvedStlSurface.faces.size() > 4
                && interiorLowerHits == 0;
        } catch(const std::exception& error) {
            std::cerr << "Generated plate scan: " << error.what() << '\n';
            return false;
        }
    }
    bool diagnoseWuPlateEdges()
    {
        using namespace spraythickness::published;
        constexpr double halfSide = 0.01;
        constexpr double overrun = 0.007;
        constexpr double speed = 0.06;
        constexpr double timeStep = 0.01;
        struct Knot
        {
            double time;
            double x;
        };
        const std::vector<Knot> knots = {
            { 0.0, -halfSide - overrun },
            { overrun / speed, -halfSide },
            { (overrun + 2.0 * halfSide) / speed, halfSide },
            { (2.0 * overrun + 2.0 * halfSide) / speed,
                halfSide + overrun },
            { (2.0 * overrun + 2.0 * halfSide) / speed,
                halfSide + overrun },
            { (3.0 * overrun + 2.0 * halfSide) / speed, halfSide },
            { (3.0 * overrun + 4.0 * halfSide) / speed, -halfSide },
            { (4.0 * overrun + 4.0 * halfSide) / speed,
                -halfSide - overrun }
        };
        const auto positionAt = [&knots](double time) {
            for(std::size_t index = 1; index < knots.size(); ++index) {
                if(time <= knots[index].time) {
                    const double duration = knots[index].time
                        - knots[index - 1].time;
                    const double alpha = duration > 0.0
                        ? (time - knots[index - 1].time) / duration : 0.0;
                    return knots[index - 1].x
                        + alpha * (knots[index].x - knots[index - 1].x);
                }
            }
            return knots.back().x;
        };
        const auto poses = [&knots, &positionAt, timeStep](int mode) {
            std::vector<Knot> samples;
            if(mode != 1) {
                const std::vector<Knot> controls = mode == 2
                    ? std::vector<Knot>{ knots.front(), knots[3], knots.back() }
                    : knots;
                samples.push_back(controls.front());
                for(std::size_t index = 1; index < controls.size(); ++index) {
                    const Knot first = controls[index - 1];
                    const Knot last = controls[index];
                    const double duration = last.time - first.time;
                    for(std::size_t step = 1;
                        step * timeStep < duration - 1.0e-12; ++step) {
                        const double alpha = step * timeStep / duration;
                        samples.push_back({ first.time + step * timeStep,
                            first.x + alpha * (last.x - first.x) });
                    }
                    samples.push_back(last);
                }
            } else {
                for(std::size_t step = 0;
                    step * timeStep < knots.back().time - 1.0e-12;
                    ++step) {
                    const double time = step * timeStep;
                    samples.push_back({ time, positionAt(time) });
                }
                samples.push_back(knots.back());
            }
            std::vector<SprayPose> result;
            result.reserve(samples.size());
            for(std::size_t index = 0; index < samples.size(); ++index) {
                Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
                if(index + 1 < samples.size()) {
                    const double duration = samples[index + 1].time
                        - samples[index].time;
                    if(duration > 0.0) {
                        velocity.x() = (samples[index + 1].x
                            - samples[index].x) / duration;
                    }
                }
                result.push_back({ samples[index].time,
                    Eigen::Vector3d(samples[index].x, 0.0, 0.03),
                    -Eigen::Vector3d::UnitZ(), velocity,
                    index + 1 < samples.size() });
            }
            return result;
        };
        WuInputModel input;
        input.substrate.vertices = {
            { -halfSide, -halfSide, 0.0 },
            { halfSide, -halfSide, 0.0 },
            { halfSide, halfSide, 0.0 },
            { -halfSide, halfSide, 0.0 }
        };
        input.substrate.faces = { { 0, 1, 2 }, { 0, 2, 3 } };
        input.generatedPlateStack = true;
        WuParameters parameters;
        parameters.peakCylinderHeightMeters = 5.342534219921746e-6;
        parameters.gaussianSigmaMeters = 0.003589657150328816;
        parameters.maximumDeflectionRadians = 0.09795983193;
        parameters.rayAngularStepRadians = 0.0081633193275;
        parameters.cylinderRadiusMeters = 2.0e-5;
        parameters.sprayAngleRelativeDepositionEfficiency.coefficients =
            { 16.57062, -0.9661822, 0.020547, -0.0001842, 5.98125e-7 };
        parameters.sprayDistanceRelativeDepositionEfficiency.coefficients =
            { -0.683333, 0.208855, -0.007281, 0.000072 };
        parameters.traverseSpeedPeakCorrectionFactor.coefficients =
            { 1.239, -0.0054 };
        for(const int mode : { 0, 1, 2 }) {
            input.nozzleTrajectory = poses(mode);
            const WuResult result = WuReproducer::run(input, parameters);
            std::cerr << (mode == 0 ? "segmented"
                : mode == 1 ? "uniform" : "leg-only")
                      << " poses=" << input.nozzleTrajectory.size()
                      << " cylinders=" << result.depositedCylinders.size()
                      << '\n';
            for(const double x : { -halfSide, -0.008, 0.0, 0.008,
                    halfSide }) {
                std::size_t count = 0;
                std::size_t coatingHits = 0;
                double maxHeight = 0.0;
                double pointHeight = 0.0;
                for(const WuDepositedCylinder& cylinder :
                    result.depositedCylinders) {
                    const Eigen::Vector3d top = cylinder.baseCenter
                        + cylinder.heightMeters * cylinder.growthDirection;
                    if(std::abs(cylinder.baseCenter.x() - x) <= 0.0005
                        && std::abs(cylinder.baseCenter.y()) <= 0.0005) {
                        ++count;
                        coatingHits += cylinder.baseCenter.z() > 1.0e-10;
                        maxHeight = std::max(maxHeight, top.z());
                    }
                    const double dx = top.x() - x;
                    if(dx * dx + top.y() * top.y()
                        <= cylinder.radiusMeters * cylinder.radiusMeters) {
                        pointHeight = std::max(pointHeight, top.z());
                    }
                }
                std::cerr << "  x=" << x * 1000.0
                          << " mm, band cylinders=" << count
                          << ", on-coating=" << coatingHits
                          << ", band peak=" << maxHeight * 1.0e6
                          << " um, centerline=" << pointHeight * 1.0e6
                          << " um\n";
            }
        }
        return true;
    }
}

int main(int argc, char** argv)
{
    if(argc == 2 && std::string(argv[1]) == "--wu-duration-only") {
        return testWuDepositionConservesSampleDuration() ? 0 : 41;
    }
    if(argc == 2 && std::string(argv[1]) == "--wu-edge-diagnostic") {
        return diagnoseWuPlateEdges() ? 0 : 39;
    }
    if(argc == 2 && std::string(argv[1]) == "--wu-only") {
        bool passed = true;
        const auto check = [&passed](bool result, const char* name) {
            if(!result) {
                std::cerr << name << " failed.\n";
                passed = false;
            }
        };
        check(testWuDynamicCylinderDeposition(), "Wu dynamic deposition");
        check(testWuDepositionConservesSampleDuration(), "Wu time-weighted deposition");
        check(testWuPlateRejectsBackFace(), "Wu plate back-face rejection");
        check(testWuSidewallUsesSubstrateAngle(), "Wu sidewall angle");
        check(testWuInvalidSpeedCorrection(), "Wu speed correction");
        check(testWuCylinderStackAndPlateOcclusion(), "Wu plate occlusion");
        check(testWuRecommendedPlateScan(), "Wu recommended scan");
        check(testWuUsesOriginalPoses(), "Wu original poses");
        check(testWuGeneratedPlateCoverageAndPoseBatch(), "Wu plate coverage");
        check(testWuSixDiscreteSprayPoses(), "Wu six discrete poses");
        check(testWuZeroDepositionDiagnostic(), "Wu zero-deposition log");
        if(!passed) {
            std::cerr << "Wu reproduction tests failed.\n";
            return 38;
        }
        return 0;
    }
    if(argc == 2 && std::string(argv[1]) == "--fuke-only") {
        if(!testFukePublishedEquation()) {
            std::cerr << "Fuke published-equation test failed.\n";
            return 11;
        }
        if(!testFukeMillimeterFace()) {
            std::cerr << "Fuke millimeter-face test failed.\n";
            return 19;
        }
        if(!testFukeSourceLineOfSight()) {
            std::cerr << "Fuke line-of-sight test failed.\n";
            return 28;
        }
        if(!testFukeTimeAccumulation()) {
            std::cerr << "Fuke time-accumulation test failed.\n";
            return 34;
        }
        if(!testFukeDistanceAndAngles()) {
            std::cerr << "Fuke distance-and-angle test failed.\n";
            return 35;
        }
        if(!testFukeDensePlateOcclusion()) {
            std::cerr << "Fuke dense-plate test failed.\n";
            return 36;
        }
        if(!testFukeMovingDensePlateOcclusion()) {
            std::cerr << "Fuke moving dense-plate test failed.\n";
            return 39;
        }
        if(!testFukeQuadrilateralVisibility()) {
            std::cerr << "Fuke quadrilateral-visibility test failed.\n";
            return 37;
        }
        return 0;
    }
    if(argc == 2 && std::string(argv[1]) == "--vanerio-only") {
        return testVanerioDynamicSurfaceUpdate()
            && testVanerioDensePlateDeposition()
            && testVanerioNearestOutletVertex()
            && testVanerioCylinderAcrossGridBoundary()
            && testVanerioDistanceProfileSupport()
            && testVanerioInterPassRemeshing()
            && testVanerioSharedEdgeRemeshing() ? 0 : 29;
    }
    if(argc == 2 && std::string(argv[1]) == "--tzinava-only") {
        if(!testTzinavaSpraySwitchBoundary()) {
            std::cerr << "Tzinava spray-switch boundary test failed.\n";
            return 30;
        }
        if(!testTzinavaDensePlateOcclusion()) {
            std::cerr << "Tzinava dense-plate occlusion test failed.\n";
            return 18;
        }
        if(!testTzinavaOrthographicTriangleShadow()) {
            std::cerr << "Tzinava orthographic shadow test failed.\n";
            return 26;
        }
        if(!testTzinavaMovingScanAccumulation()) {
            std::cerr << "Tzinava moving-scan accumulation test failed.\n";
            return 24;
        }
        if(!testTzinavaContinuousReversalAndBackface()) {
            std::cerr << "Tzinava reversal and backface test failed.\n";
            return 31;
        }
        if(!testTzinavaMovingPlateOcclusion()) {
            std::cerr << "Tzinava moving-plate occlusion test failed.\n";
            return 25;
        }
        return 0;
    }
    if(argc == 2 && std::string(argv[1]) == "--dynamic-only") {
        if(!testDynamicSurfaceBatchUpdate()) {
            std::cerr << "Dynamic-surface batch test failed.\n";
            return 15;
        }
        if(!testDynamicSurfacePreservesOccludedSubstrate()) {
            std::cerr << "Dynamic-surface substrate-preservation test failed.\n";
            return 37;
        }
        if(!testDynamicSurfaceDenseTrajectory()) {
            std::cerr << "Dynamic-surface dense-trajectory test failed.\n";
            return 17;
        }
        if(!testDynamicSurfaceSmallRegion()) {
            std::cerr << "Dynamic-surface small-region test failed.\n";
            return 22;
        }
        if(!testDynamicSurfaceGeneratedPlateScan()) {
            std::cerr << "Dynamic-surface generated-plate scan test failed.\n";
            return 23;
        }
        return 0;
    }
    if(!testFukePublishedEquation()) {
        std::cerr << "Fuke published-equation test failed.\n";
        return 11;
    }
    if(!testFukeMillimeterFace()) {
        std::cerr << "Fuke millimeter-face test failed.\n";
        return 19;
    }
    if(!testFukeSourceLineOfSight()) {
        std::cerr << "Fuke line-of-sight test failed.\n";
        return 28;
    }
    if(!testFukeTimeAccumulation()
        || !testFukeDistanceAndAngles()
        || !testFukeDensePlateOcclusion()
        || !testFukeQuadrilateralVisibility()
        || !testFukeMovingDensePlateOcclusion()) {
        std::cerr << "Fuke accumulation, geometry, or dense-plate test failed.\n";
        return 33;
    }
    if(!testTzinavaStationaryBranch()) {
        std::cerr << "Tzinava stationary-branch test failed.\n";
        return 12;
    }
    if(!testTzinavaSpraySwitchBoundary()) {
        std::cerr << "Tzinava spray-switch boundary test failed.\n";
        return 30;
    }
    if(!testTzinavaDensePlateOcclusion()) {
        std::cerr << "Tzinava dense-plate occlusion test failed.\n";
        return 18;
    }
    if(!testTzinavaOrthographicTriangleShadow()) {
        std::cerr << "Tzinava orthographic shadow test failed.\n";
        return 26;
    }
    if(!testTzinavaMovingScanAccumulation()) {
        std::cerr << "Tzinava moving-scan accumulation test failed.\n";
        return 24;
    }
    if(!testTzinavaContinuousReversalAndBackface()) {
        std::cerr << "Tzinava reversal and backface test failed.\n";
        return 31;
    }
    if(!testTzinavaMovingPlateOcclusion()) {
        std::cerr << "Tzinava moving-plate occlusion test failed.\n";
        return 25;
    }
    if(!testWuDynamicCylinderDeposition()) {
        std::cerr << "Wu dynamic-cylinder test failed.\n";
        return 13;
    }
    if(!testWuDepositionConservesSampleDuration()) {
        std::cerr << "Wu time-weighted deposition test failed.\n";
        return 41;
    }
    if(!testWuPlateRejectsBackFace()) {
        std::cerr << "Wu plate back-face test failed.\n";
        return 40;
    }
    if(!testWuSidewallUsesSubstrateAngle()) {
        std::cerr << "Wu sidewall-angle test failed.\n";
        return 39;
    }
    if(!testWuInvalidSpeedCorrection()) {
        std::cerr << "Wu invalid-speed correction test failed.\n";
        return 20;
    }
    if(!testWuCylinderStackAndPlateOcclusion()) {
        std::cerr << "Wu cylinder-stack and plate-occlusion test failed.\n";
        return 16;
    }
    if(!testWuRecommendedPlateScan()) {
        std::cerr << "Wu recommended plate-scan test failed.\n";
        return 34;
    }
    if(!testWuUsesOriginalPoses()) {
        std::cerr << "Wu original-pose test failed.\n";
        return 1;
    }
    if(!testWuGeneratedPlateCoverageAndPoseBatch()) {
        std::cerr << "Wu generated-plate coverage test failed.\n";
        return 35;
    }
    if(!testWuSixDiscreteSprayPoses()) {
        std::cerr << "Wu six-discrete-poses test failed.\n";
        return 1;
    }
    if(!testWuZeroDepositionDiagnostic()) {
        std::cerr << "Wu zero-deposition diagnostic test failed.\n";
        return 36;
    }
    if(!testVanerioDynamicSurfaceUpdate()) {
        std::cerr << "Vanerio dynamic-surface test failed.\n";
        return 14;
    }
    if(!testVanerioDensePlateDeposition()) {
        std::cerr << "Vanerio dense-plate deposition test failed.\n";
        return 21;
    }
    if(!testVanerioNearestOutletVertex()) {
        std::cerr << "Vanerio nearest-outlet-vertex test failed.\n";
        return 27;
    }
    if(!testVanerioCylinderAcrossGridBoundary()) {
        std::cerr << "Vanerio cylinder boundary test failed.\n";
        return 31;
    }
    if(!testVanerioDistanceProfileSupport()) {
        std::cerr << "Vanerio distance-profile support test failed.\n";
        return 32;
    }
    if(!testVanerioInterPassRemeshing()) {
        std::cerr << "Vanerio inter-pass remeshing test failed.\n";
        return 28;
    }
    if(!testVanerioSharedEdgeRemeshing()) {
        std::cerr << "Vanerio shared-edge remeshing test failed.\n";
        return 29;
    }
    if(!testDynamicSurfaceBatchUpdate()) {
        std::cerr << "Dynamic-surface batch test failed.\n";
        return 15;
    }
    if(!testDynamicSurfacePreservesOccludedSubstrate()) {
        std::cerr << "Dynamic-surface substrate-preservation test failed.\n";
        return 37;
    }
    if(!testDynamicSurfaceDenseTrajectory()) {
        std::cerr << "Dynamic-surface dense-trajectory test failed.\n";
        return 17;
    }
    if(!testDynamicSurfaceSmallRegion()) {
        std::cerr << "Dynamic-surface small-region test failed.\n";
        return 22;
    }
    if(!testDynamicSurfaceGeneratedPlateScan()) {
        std::cerr << "Dynamic-surface generated-plate scan test failed.\n";
        return 23;
    }
    return 0;
}
