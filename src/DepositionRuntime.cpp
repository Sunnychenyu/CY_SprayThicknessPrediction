#include <SprayThicknessPrediction/ThicknessPrediction.h>
#include <algorithm>
#include <cmath>
#include <limits>
namespace spraythickness {
bool samePaperGaussianParameters(const PaperGaussianParameters& a,const PaperGaussianParameters& b) {
    const auto& x=a.calibratedHistory;const auto& y=b.calibratedHistory;
    return a.gaussianComponents==b.gaussianComponents && x.continuousInput==y.continuousInput
        && x.referencePeakRateMillimetersPerSecond==y.referencePeakRateMillimetersPerSecond
        && a.useCalibratedHistory==b.useCalibratedHistory
        && x.correctionAmplitude==y.correctionAmplitude && x.historyScaleSeconds==y.historyScaleSeconds
        && x.referenceHistorySeconds==y.referenceHistorySeconds && x.coolingTimeSeconds==y.coolingTimeSeconds
        && x.activityThresholdRatio==y.activityThresholdRatio
        && a.calibrationModelId==b.calibrationModelId && a.useCalibratedEfficiencies==b.useCalibratedEfficiencies
        && a.densityMultiplier==b.densityMultiplier && a.useCalibratedDensity==b.useCalibratedDensity
        && a.distanceDensity==b.distanceDensity && a.angleDensity==b.angleDensity
        && a.densityDistanceMinMm==b.densityDistanceMinMm && a.densityDistanceMaxMm==b.densityDistanceMaxMm
        && a.densityAngleMinDegrees==b.densityAngleMinDegrees && a.densityAngleMaxDegrees==b.densityAngleMaxDegrees
        && a.efficiencyIncludesGeometry==b.efficiencyIncludesGeometry
        && a.distanceEfficiency==b.distanceEfficiency && a.angleEfficiency==b.angleEfficiency
        && a.amplitudeMillimeters==b.amplitudeMillimeters && a.rotationRadians==b.rotationRadians
        && a.phiOffsetRadians==b.phiOffsetRadians && a.psiOffsetRadians==b.psiOffsetRadians
        && a.sigmaPhiRadians==b.sigmaPhiRadians && a.sigmaPsiRadians==b.sigmaPsiRadians
        && a.referenceDistanceMeters==b.referenceDistanceMeters && a.referenceAngleDegrees==b.referenceAngleDegrees
        && a.referenceExposureSeconds==b.referenceExposureSeconds;
}
double evaluateCalibratedEfficiency(const CalibratedEfficiencyCurve& c,double x){
    if(c.lowerRampEnd>0 && x<c.lowerRampEnd)
        return std::max(0.0,evaluateCalibratedEfficiency(c,c.lowerRampEnd))*std::max(0.0,x)/c.lowerRampEnd;
    const auto p=[&](int i){return i<int(c.parameters.size())?c.parameters[i]:0.0;};
    const auto gauss=[&](int i){const double s=std::max(std::abs(p(i+2)),1e-9);return p(i)*std::exp(-std::pow((x-p(i+1))/s,2));};
    const auto denominator=[](double v){return std::abs(v)<1e-9?1e-9:v;};
    switch(c.kind){
    case EfficiencyCurveKind::Unit:return 1;
    case EfficiencyCurveKind::LegacyDistance:return 1/(.001443*x*x-.3073*x+28.31);
    case EfficiencyCurveKind::LegacyAngle:x=std::clamp(x,0.0,180.0);return .0002564*x*x+.01058*x+.6358;
    case EfficiencyCurveKind::Poly1:case EfficiencyCurveKind::Poly2:case EfficiencyCurveKind::Poly3:case EfficiencyCurveKind::Poly4:{double v=0;for(double a:c.parameters)v=v*x+a;return v;}
    case EfficiencyCurveKind::Exp1:return p(0)*std::exp(p(1)*x);
    case EfficiencyCurveKind::Exp2:return p(0)*std::exp(p(1)*x)+p(2)*std::exp(p(3)*x);
    case EfficiencyCurveKind::Power1:return p(0)*std::pow(std::max(x,1e-9),p(1));
    case EfficiencyCurveKind::Power2:return p(0)*std::pow(std::max(x,1e-9),p(1))+p(2);
    case EfficiencyCurveKind::Gauss1:return gauss(0);
    case EfficiencyCurveKind::Gauss2:return gauss(0)+gauss(3);
    case EfficiencyCurveKind::Rat11:return (p(0)*x+p(1))/denominator(x+p(2));
    case EfficiencyCurveKind::Rat21:return (p(0)*x*x+p(1)*x+p(2))/denominator(x+p(3));
    case EfficiencyCurveKind::Rat22:return (p(0)*x*x+p(1)*x+p(2))/denominator(x*x+p(3)*x+p(4));
    case EfficiencyCurveKind::LinearNodes:
        if(c.x.empty() || c.x.size()!=c.y.size())return std::numeric_limits<double>::quiet_NaN();
        if(x<=c.x.front())return c.y.front();if(x>=c.x.back())return c.y.back();
        {const auto hi=std::upper_bound(c.x.begin(),c.x.end(),x)-c.x.begin();const auto lo=hi-1;
        const double t=(x-c.x[lo])/(c.x[hi]-c.x[lo]);return c.y[lo]*(1-t)+c.y[hi]*t;}
    }
    return 1;
}
double paperGaussianDistanceEfficiency(const PaperGaussianParameters& p,double d){
    if(p.useCalibratedEfficiencies)return std::max(evaluateCalibratedEfficiency(p.distanceEfficiency,d),1e-9);
    const double v=std::log(std::max(d,1e-6)/106.8615)/.0839;return .9789+.1456*std::exp(-v*v);
}
double paperGaussianAngleEfficiency(const PaperGaussianParameters& p,double a){
    if(p.useCalibratedEfficiencies)return std::max(evaluateCalibratedEfficiency(p.angleEfficiency,a),0.0);
    return std::max(1e-6,(40.81414+.73808*a-.000838598*a*a)*.01);
}
double paperGaussianPattern(const PaperGaussianParameters& p,double phi,double psi){
    if(!p.gaussianComponents.empty()){
        double result=0;for(const auto& c:p.gaussianComponents){const double u=(phi-c.phi0)/c.sigmaPhi,v=(psi-c.psi0)/c.sigmaPsi;
            result+=c.coefficient*std::exp(-.5*(u*u+v*v));}return result;
    }
    const double x=phi-p.phiOffsetRadians,y=psi-p.psiOffsetRadians,c=std::cos(p.rotationRadians),s=std::sin(p.rotationRadians);
    const double u=(c*x+s*y)/p.sigmaPhiRadians,v=(-s*x+c*y)/p.sigmaPsiRadians;return std::exp(-.5*(u*u+v*v));
}
double paperGaussianSupportAngle(const PaperGaussianParameters& p,double radius){
    constexpr double limit=1.5533430342749532;
    if(p.gaussianComponents.empty())return std::clamp(
        std::max(std::abs(p.phiOffsetRadians),std::abs(p.psiOffsetRadians))
        +radius*std::max(std::abs(p.sigmaPhiRadians),std::abs(p.sigmaPsiRadians)),1e-6,limit);
    double coefficientSum=0,angle=0;
    for(const auto& c:p.gaussianComponents)coefficientSum+=c.coefficient;
    // Spatial filtering compares absolute mixture density, not a unit-peak Gaussian.
    const double supportRadius=std::sqrt(radius*radius+2*std::log(std::max(1.0,coefficientSum)));
    for(const auto& c:p.gaussianComponents){
        const double u=std::abs(c.phi0)+supportRadius*c.sigmaPhi,v=std::abs(c.psi0)+supportRadius*c.sigmaPsi;
        if(u>=limit || v>=limit)return limit;
        angle=std::max(angle,std::atan(std::hypot(std::tan(u),std::tan(v))));
    }
    return std::clamp(angle,1e-6,limit);
}

double paperGaussianDensityCorrection(const PaperGaussianParameters& p,double d,double a){
    if(!p.useCalibratedDensity)return p.densityMultiplier;
    const double rd=evaluateCalibratedEfficiency(p.distanceDensity,p.referenceDistanceMeters*1000);
    const double ra=evaluateCalibratedEfficiency(p.angleDensity,p.referenceAngleDegrees);
    const double localD=evaluateCalibratedEfficiency(p.distanceDensity,std::clamp(d,p.densityDistanceMinMm,p.densityDistanceMaxMm));
    const double localA=evaluateCalibratedEfficiency(p.angleDensity,std::clamp(a,p.densityAngleMinDegrees,p.densityAngleMaxDegrees));
    return p.densityMultiplier*rd*ra/std::max(localD*localA,1e-12);
}
double paperGaussianRateMillimeters(const PaperGaussianParameters& p,double d,double angle,double phi,double psi){
    if(d<=0 || angle<=0)return 0;
    const double ref=p.referenceDistanceMeters*1000,pi=3.14159265358979323846;
    const double projection=std::cos((90-angle)*pi/180)/std::max(std::cos((90-p.referenceAngleDegrees)*pi/180),1e-9);
    const double efficiency=paperGaussianDistanceEfficiency(p,d)/paperGaussianDistanceEfficiency(p,ref)
        *paperGaussianAngleEfficiency(p,angle)/paperGaussianAngleEfficiency(p,p.referenceAngleDegrees);
    return p.amplitudeMillimeters/std::max(p.referenceExposureSeconds,1e-9)*(p.efficiencyIncludesGeometry?1.0:std::pow(ref/d,2)*projection)*efficiency
        *paperGaussianPattern(p,phi,psi)*paperGaussianDensityCorrection(p,d,angle);
}
double applyPaperGaussianHistory(double rate,double pattern,double dt,double elapsed,double& tau,const HeatHistoryParameters& h,double densityCorrection)
{
    tau*=std::exp(-std::max(0.0,elapsed)/h.coolingTimeSeconds);
    if(h.continuousInput){
        const double d=std::clamp(rate/(densityCorrection*h.referencePeakRateMillimetersPerSecond),0.0,1.0);
        tau+=h.coolingTimeSeconds*(-std::expm1(-dt/h.coolingTimeSeconds))*d;
    }else if(pattern>=h.activityThresholdRatio)tau+=dt;
    const double f=1+h.correctionAmplitude*(std::exp(-tau/h.historyScaleSeconds)-std::exp(-h.referenceHistorySeconds/h.historyScaleSeconds));
    return rate*dt*(h.continuousInput?f:std::max(1.0,f));
}

}
