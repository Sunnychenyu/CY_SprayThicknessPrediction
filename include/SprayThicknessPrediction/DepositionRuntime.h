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
    bool operator==(const CalibratedEfficiencyCurve& other) const {
        return kind==other.kind && parameters==other.parameters && x==other.x && y==other.y;
    }
};
double evaluateCalibratedEfficiency(const CalibratedEfficiencyCurve& curve,double value);
struct PaperGaussianParameters;
bool samePaperGaussianParameters(const PaperGaussianParameters& a,const PaperGaussianParameters& b);
double paperGaussianDistanceEfficiency(const PaperGaussianParameters& p,double distanceMillimeters);
double paperGaussianAngleEfficiency(const PaperGaussianParameters& p,double angleDegrees);
// CPU reference of the exact base deposition law, before optional heat history / occlusion.
double paperGaussianRateMillimeters(const PaperGaussianParameters& p,double distanceMillimeters,
    double angleDegrees,double phiRadians,double psiRadians);
}
