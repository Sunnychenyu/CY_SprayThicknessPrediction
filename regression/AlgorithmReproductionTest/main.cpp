#include <SprayThicknessPrediction/DynamicSurfaceReproduction.h>
#include <SprayThicknessPrediction/FukeReproduction.h>
#include <SprayThicknessPrediction/TanakaReproduction.h>
#include <SprayThicknessPrediction/TzinavaReproduction.h>
#include <SprayThicknessPrediction/VanerioReproduction.h>
#include <SprayThicknessPrediction/WuReproduction.h>

#include <cmath>
#include <iostream>

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

    bool testTanakaPublishedEquation()
    {
        using namespace spraythickness::published;
        TanakaInputModel input;
        input.targetPoints.push_back({ Eigen::Vector3d::Zero(),
            Eigen::Vector3d::UnitZ() });
        input.sprayPoses.push_back({ 0.0, Eigen::Vector3d(0.0, 0.0, 0.1),
            -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        input.sprayPoses.push_back({ 1.0, Eigen::Vector3d(0.0, 0.0, 0.1),
            -Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero(), true });
        TanakaParameters parameters;
        parameters.paintDischargeCubicMetersPerSecond = 1.0e-9;
        parameters.referenceDistanceMeters = 0.1;
        parameters.referenceSigmaXMeters = 0.01;
        parameters.referenceSigmaYMeters = 0.02;
        const TanakaResult result = TanakaReproducer::run(input, parameters);
        const double expected = 1.0e-9
            / (2.0 * 3.14159265358979323846 * 0.01 * 0.02);
        return result.pointThicknessMeters.size() == 1
            && approximatelyEqual(result.pointThicknessMeters.front(),
                expected, 1.0e-12);
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
        return result.faceThicknessMeters.size() == 2
            && approximatelyEqual(result.faceThicknessMeters[0], 1.0e-6, 1.0e-12)
            && approximatelyEqual(result.faceThicknessMeters[1], 1.0e-6, 1.0e-12);
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
            { 0.05, 0.15 }, { 1.0, 1.0 }
        };

        const VanerioResult result = VanerioReproducer::run(input, parameters);
        if(result.evolvedStlSurface.vertices.size() != 3) {
            return false;
        }
        for(const Eigen::Vector3d& vertex : result.evolvedStlSurface.vertices) {
            if(!approximatelyEqual(vertex.z(), 3.0e-6, 1.0e-12)) {
                return false;
            }
        }
        return true;
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

        const DynamicSurfaceResult result =
            DynamicSurfaceReproducer::run(input, parameters);
        return !result.depositedCylinders.empty()
            && !result.evolvedStlSurface.empty()
            && result.evolvedStlSurface.faces.size() > 2;
    }
}

int main()
{
    if(!testTanakaPublishedEquation()) {
        std::cerr << "Tanaka published-equation test failed.\n";
        return 10;
    }
    if(!testFukePublishedEquation()) {
        std::cerr << "Fuke published-equation test failed.\n";
        return 11;
    }
    if(!testTzinavaStationaryBranch()) {
        std::cerr << "Tzinava stationary-branch test failed.\n";
        return 12;
    }
    if(!testWuDynamicCylinderDeposition()) {
        std::cerr << "Wu dynamic-cylinder test failed.\n";
        return 13;
    }
    if(!testVanerioDynamicSurfaceUpdate()) {
        std::cerr << "Vanerio dynamic-surface test failed.\n";
        return 14;
    }
    if(!testDynamicSurfaceBatchUpdate()) {
        std::cerr << "Dynamic-surface batch test failed.\n";
        return 15;
    }
    return 0;
}
