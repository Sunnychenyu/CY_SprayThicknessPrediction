#pragma once

#include <SprayThicknessPrediction/PublishedReproductionCommon.h>

#include <limits>

namespace spraythickness::published
{
    enum class TzinavaBeamKind
    {
        Cylindrical,
        CylindricalConical
    };

    struct TzinavaLookupTable
    {
        std::vector<double> standOffDistancesMeters;
        std::vector<double> impactAnglesDegrees;
        std::vector<double> thicknessMeters;

        void validate() const;
        double interpolate(
            double standOffDistanceMeters,
            double impactAngleDegrees) const;
    };

    struct TzinavaInputModel
    {
        TriangleMesh initialMesh;
        std::vector<SprayPose> gunTrajectory;
        Eigen::Vector3d objectRotationOrigin = Eigen::Vector3d::Zero();
        Eigen::Vector3d objectRotationAxis = Eigen::Vector3d::UnitZ();
        double objectAngularSpeedRadiansPerSecond{ 0.0 };
    };

    struct TzinavaParameters
    {
        TzinavaBeamKind beamKind{ TzinavaBeamKind::Cylindrical };
        double beamRadiusMeters{ std::numeric_limits<double>::quiet_NaN() };
        double cylindricalLengthMeters{ 0.0 };
        double coneHalfAngleRadians{ 0.0 };
        double gaussianSigmaMeters{ std::numeric_limits<double>::quiet_NaN() };
        bool gaussianRadialProfile{ true };
        double speedCoefficientB{ 4.8e-6 };
        double speedCoefficientC{ 3.5e-4 };
        double referenceSpotSpeedMillimetersPerSecond{ 502.0 };
        double timeStepOverlapFactor{ 0.5 };
        double stationaryTimeStepSeconds{
            std::numeric_limits<double>::quiet_NaN() };
        double lookupReferenceDwellSeconds{
            std::numeric_limits<double>::quiet_NaN() };
        TzinavaLookupTable thicknessTable;

        void validate() const;
    };

    struct TzinavaResult
    {
        TriangleMesh subdividedMesh;
        std::vector<Eigen::Vector3d> faceCentroids;
        std::vector<double> faceThicknessMeters;
        ReproductionStatistics statistics;
        bool canceled{ false };
    };

    class TzinavaReproducer
    {
    public:
        static TzinavaResult run(
            const TzinavaInputModel& input,
            const TzinavaParameters& parameters,
            const ReproductionExecution& execution = {});
    };
}
