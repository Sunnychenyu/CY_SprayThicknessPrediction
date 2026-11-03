#pragma once
#include <string>
#include <vector>
namespace spraythickness {
// Values are evaluated in mm / degrees; coefficients follow descending polynomial order.
enum class EfficiencyCurveKind { Unit=0,LegacyDistance,LegacyAngle,Poly1,Poly2,Poly3,Poly4,
    Exp1,Exp2,Power1,Power2,Gauss1,Gauss2,Rat11,Rat21,Rat22,LinearNodes };
struct CalibratedEfficiencyCurve {
    EfficiencyCurveKind kind=EfficiencyCurveKind::Unit;
    std::vector<double> parameters,x,y;
    // Optional linear continuation from (0,0) to the curve at this input value.
    double lowerRampEnd=0.0;
    bool operator==(const CalibratedEfficiencyCurve& other) const {
        return kind==other.kind && parameters==other.parameters && x==other.x && y==other.y
            && lowerRampEnd==other.lowerRampEnd;
    }
};
double evaluateCalibratedEfficiency(const CalibratedEfficiencyCurve& curve,double value);
struct PaperGaussianParameters;
struct HeatHistoryParameters;
// elapsedSeconds is between interval ends and includes cooling gaps; rates are mm/s.
double applyPaperGaussianHistory(double rateMillimetersPerSecond,double pattern,double dtSeconds,
    double elapsedSeconds,double& tau,const HeatHistoryParameters& parameters,double densityCorrection=1.0);
bool samePaperGaussianParameters(const PaperGaussianParameters& a,const PaperGaussianParameters& b);
double paperGaussianDensityCorrection(const PaperGaussianParameters& p,double distanceMm,double angleDegrees);
double paperGaussianPattern(const PaperGaussianParameters& p,double phiRadians,double psiRadians);
double paperGaussianSupportAngle(const PaperGaussianParameters& p,double gaussianRadius);
double paperGaussianDistanceEfficiency(const PaperGaussianParameters& p,double distanceMillimeters);
double paperGaussianAngleEfficiency(const PaperGaussianParameters& p,double angleDegrees);
// CPU reference of the exact base deposition law, before optional heat history / occlusion.
double paperGaussianRateMillimeters(const PaperGaussianParameters& p,double distanceMillimeters,
    double angleDegrees,double phiRadians,double psiRadians);
}
