#include <SprayThicknessPrediction/ThicknessPrediction.h>
#include <algorithm>
#include <cmath>
#include <limits>
namespace spraythickness {
bool samePaperGaussianParameters(const PaperGaussianParameters& a,const PaperGaussianParameters& b) {
    const auto& x=a.calibratedHistory;const auto& y=b.calibratedHistory;
    return x.continuousInput==y.continuousInput
        && x.referencePeakRateMillimetersPerSecond==y.referencePeakRateMillimetersPerSecond
        && a.useCalibratedHistory==b.useCalibratedHistory
        && x.correctionAmplitude==y.correctionAmplitude && x.historyScaleSeconds==y.historyScaleSeconds
        && x.referenceHistorySeconds==y.referenceHistorySeconds && x.coolingTimeSeconds==y.coolingTimeSeconds
        && x.activityThresholdRatio==y.activityThresholdRatio
        && a.calibrationModelId==b.calibrationModelId && a.useCalibratedEfficiencies==b.useCalibratedEfficiencies
        && a.efficiencyIncludesGeometry==b.efficiencyIncludesGeometry
        && a.distanceEfficiency==b.distanceEfficiency && a.angleEfficiency==b.angleEfficiency
        && a.amplitudeMillimeters==b.amplitudeMillimeters && a.rotationRadians==b.rotationRadians
        && a.phiOffsetRadians==b.phiOffsetRadians && a.psiOffsetRadians==b.psiOffsetRadians
        && a.sigmaPhiRadians==b.sigmaPhiRadians && a.sigmaPsiRadians==b.sigmaPsiRadians
        && a.referenceDistanceMeters==b.referenceDistanceMeters && a.referenceAngleDegrees==b.referenceAngleDegrees
        && a.referenceExposureSeconds==b.referenceExposureSeconds;
}
double evaluateCalibratedEfficiency(const CalibratedEfficiencyCurve& c,double x){
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
    if(p.useCalibratedEfficiencies)return std::max(evaluateCalibratedEfficiency(p.angleEfficiency,a),1e-9);
    return std::max(1e-6,(40.81414+.73808*a-.000838598*a*a)*.01);
}
double paperGaussianRateMillimeters(const PaperGaussianParameters& p,double d,double angle,double phi,double psi){
    if(d<=0 || angle<=0)return 0;
    const double a=std::cos(p.rotationRadians)*(phi-p.phiOffsetRadians)+std::sin(p.rotationRadians)*(psi-p.psiOffsetRadians);
    const double b=-std::sin(p.rotationRadians)*(phi-p.phiOffsetRadians)+std::cos(p.rotationRadians)*(psi-p.psiOffsetRadians);
    const double ref=p.referenceDistanceMeters*1000,pi=3.14159265358979323846;
    const double projection=std::cos((90-angle)*pi/180)/std::max(std::cos((90-p.referenceAngleDegrees)*pi/180),1e-9);
    const double efficiency=paperGaussianDistanceEfficiency(p,d)/paperGaussianDistanceEfficiency(p,ref)
        *paperGaussianAngleEfficiency(p,angle)/paperGaussianAngleEfficiency(p,p.referenceAngleDegrees);
    return p.amplitudeMillimeters/std::max(p.referenceExposureSeconds,1e-9)*(p.efficiencyIncludesGeometry?1.0:std::pow(ref/d,2)*projection)*efficiency
        *std::exp(-.5*(std::pow(a/p.sigmaPhiRadians,2)+std::pow(b/p.sigmaPsiRadians,2)));
}
double applyPaperGaussianHistory(double rate,double pattern,double dt,double elapsed,double& tau,const HeatHistoryParameters& h)
{
    tau*=std::exp(-std::max(0.0,elapsed)/h.coolingTimeSeconds);
    if(h.continuousInput){
        const double d=std::clamp(rate/h.referencePeakRateMillimetersPerSecond,0.0,1.0);
        tau+=h.coolingTimeSeconds*(-std::expm1(-dt/h.coolingTimeSeconds))*(d>=h.activityThresholdRatio?d:0.0);
    }else if(pattern>=h.activityThresholdRatio)tau+=dt;
    const double f=1+h.correctionAmplitude*(std::exp(-tau/h.historyScaleSeconds)-std::exp(-h.referenceHistorySeconds/h.historyScaleSeconds));
    return rate*dt*(h.continuousInput?f:std::max(1.0,f));
}

}
