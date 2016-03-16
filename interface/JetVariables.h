//Track builder infos
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TLorentzVector.h"
//#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h" 
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
KalmanVertexFitter vtxFitter(true);
//Namespaces
using namespace reco;
using namespace edm;
using namespace std;
//Check that the track is a good track
bool is_goodtrk(Track trk,const reco::Vertex& vtx){
 bool isgoodtrk = false;
 if(trk.pt()>1 &&
   trk.hitPattern().numberOfValidHits()>=8 && 
   trk.hitPattern().numberOfValidPixelHits()>=2 &&
   trk.normalizedChi2()<5 &&
   std::abs(trk.dxy(vtx.position()))<0.2 &&
   std::abs(trk.dz(vtx.position()))<17    
   ) isgoodtrk = true;
 return isgoodtrk;
}
//Get transient tracks from track
TransientTrack get_ttrk(Track trk, const TransientTrackBuilder& ttrkbuilder){
 TransientTrack ttrk;
 ttrk = ttrkbuilder.build(&trk);
 return ttrk;
}
vector<TransientTrack> get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
 vector<TransientTrack> ttrks;
 for(uint tr=0; tr<trks.size(); tr++){
  TransientTrack ttrk = ttrkbuilder.build(&trks[tr]);
  ttrks.push_back(ttrk);
 }
 return ttrks;
}
//Get transient pv from tracks
TransientVertex get_tv(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
 TransientVertex tv;
 vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
 if(ttrks.size()>=2) tv = vtxFitter.vertex(ttrks);
 return tv;
}
TransientVertex get_ttv(vector<TransientTrack> ttrks, const TransientTrackBuilder& ttrkbuilder){
 TransientVertex tv;
 if(ttrks.size()>=2) tv = vtxFitter.vertex(ttrks);
 return tv;
}
//Get chi2 information from trk vertex
void get_chi2info(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& chi2tot, double& chi2ndf, double& chi2pval){
 if(trks.size()>=2){
  TransientVertex trks_tv = get_tv(trks, ttrkbuilder);
  if(trks_tv.isValid()){
   chi2tot  = trks_tv.totalChiSquared();
   chi2ndf  = trks_tv.degreesOfFreedom();
   chi2pval = TMath::Prob(trks_tv.totalChiSquared(),trks_tv.degreesOfFreedom());
  }
 }
}
pair<double,double> dca2trks(Track tkA, Track tkB, const TransientTrackBuilder& ttrkbuilder){
 double dca3d2trks_sig = 0;
 double dca2d2trks_sig = 0;
 TransientTrack ttkA = get_ttrk(tkA, ttrkbuilder);
 TransientTrack ttkB = get_ttrk(tkB, ttrkbuilder);
 if(ttkA.impactPointTSCP().isValid() && ttkB.impactPointTSCP().isValid()){
  FreeTrajectoryState state1 = ttkA.impactPointTSCP().theState();
  FreeTrajectoryState state2 = ttkB.impactPointTSCP().theState();
  TwoTrackMinimumDistance minDist;
  minDist.calculate(state1, state2);
  if(minDist.status()){
   std::pair<GlobalPoint,GlobalPoint> pcas = minDist.points();
   GlobalPoint pca1 = pcas.first;
   GlobalPoint pca2 = pcas.second;
   ROOT::Math::SVector<double, 3> distanceVector(pca1.x()-pca2.x(), pca1.y()-pca2.y(), pca1.z()-pca2.z());
   const float twoTkDis3D = ROOT::Math::Mag(distanceVector);
   float mass     = 0.139526;
   float massigma = mass*1e-6;
   float chi2 = 0.0f, ndf = 0.0f;
   KinematicParticleFactoryFromTransientTrack pFactory;
   RefCountedKinematicParticle tkAParticle = pFactory.particle(ttkA, mass, chi2, ndf, massigma);
   RefCountedKinematicParticle tkBParticle = pFactory.particle(ttkB, mass, chi2, ndf, massigma);
   float sig[6];
   sig[0] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,0);
   sig[1] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,1);
   sig[2] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,1);
   sig[3] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,2);
   sig[4] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,2);
   sig[5] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(2,2);
   ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca1Cov(sig, sig+6);
   sig[0] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,0);
   sig[1] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,1);
   sig[2] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,1);
   sig[3] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,2);
   sig[4] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,2);
   sig[5] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(2,2);
   ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca2Cov(sig, sig+6);
   ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > totCov = pca1Cov + pca2Cov;
   double twoTkDist3DEr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTkDis3D;
   dca3d2trks_sig = twoTkDis3D/twoTkDist3DEr;
   distanceVector(2) = 0.0;
   double twoTauDis2D    = ROOT::Math::Mag(distanceVector);
   double twoTauDist2DEr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTauDis2D;
   dca2d2trks_sig = twoTauDis2D/twoTauDist2DEr;
    }//if(minDist.status())
 }//ttkA.impactPointTSCP().isValid() && ttkB.impactPointTSCP().isValid()
 return make_pair(dca3d2trks_sig,dca2d2trks_sig);
}

//Get Num of two-trk vertices
void get_2trksinfo(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& num2v, double& numno2v,double& dca3d2t, double& dca3dno2t, double& dca2d2t, double& dca2dno2t){
 for(uint t=0; t<trks.size(); t++){
  for(uint t2=t+1; t2<trks.size(); t2++){
   vector<Track> twotrks;
   twotrks.push_back(trks[t]);  
   twotrks.push_back(trks[t2]);    
   TransientVertex tv = get_tv(twotrks, ttrkbuilder);
   if(tv.isValid() && TMath::Prob(tv.totalChiSquared(),tv.degreesOfFreedom())>0.05){
    num2v++;
    pair<double,double> dca2trks3d2d = dca2trks(trks[t], trks[t2], ttrkbuilder);
    dca3d2t += dca2trks3d2d.first;
    dca2d2t += dca2trks3d2d.second;
   }else{
    numno2v++;
    pair<double,double> dca2trks3d2d = dca2trks(trks[t], trks[t2], ttrkbuilder);
    dca3dno2t += dca2trks3d2d.first;
    dca2dno2t += dca2trks3d2d.second;
   }
  }
 }
}
//Get jet tracks
//void get_jettrks(const pat::Jet& jet, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, int& jet_ndaus, int& jet_chtrks, int& jet_chtrkspv, int& jet_chtrksnpv, vector<Track>& jetchtrks, vector<Track>& jetchtrksnpv){
 void get_jettrks(const pat::Jet& jet, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, int& jet_ndaus, int& jet_chtrks, int& jet_chtrkspv, int& jet_chtrksnpv, vector<Track>& jetchtrks, vector<Track>& jetchtrkspv,vector<Track>& jetchtrksnpv, int& jet_nleptons,double& jet_chtrks_pt,double& jet_chtrks_M,double& jet_chtrkspv_pt,double& jet_chtrkspv_M,double& jet_chtrksnpv_pt,double& jet_chtrksnpv_M){
 //Access jet daughters
 vector<CandidatePtr> jdaus(jet.daughterPtrVector());
 sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
 jet_ndaus = jdaus.size(); 
 math::PtEtaPhiMLorentzVector chtrks_LV(0,0,0,0);
 math::PtEtaPhiMLorentzVector chtrkspv_LV(0,0,0,0);
 math::PtEtaPhiMLorentzVector chtrksnpv_LV(0,0,0,0);
 for(int jd=0; jd<jet_ndaus; jd++){
  const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
  if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
  Track trk = Track(jcand.pseudoTrack());
  bool isgoodtrk = is_goodtrk(trk,vtx);
  //Minimal conditions for a track 
  //if(jcand.charge()!=0 && jcand.numberOfHits()>0)
  if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1)
  {
   if(fabs(jcand.pdgId())==11 || fabs(jcand.pdgId())==13) jet_nleptons++;
   math::PtEtaPhiMLorentzVector chtrks_LV2(jcand.pt(),jcand.eta(),jcand.phi(),jcand.mass());
   chtrks_LV += chtrks_LV2;
   jet_chtrks++; 
   jetchtrks.push_back(trk);
   //Other conditions on jet daughters
   //Using fromPV method
   if(jcand.fromPV()==pat::PackedCandidate::PVUsedInFit){
    jet_chtrkspv++;
    jetchtrkspv.push_back(trk);
    math::PtEtaPhiMLorentzVector chtrkspv_LV2(jcand.pt(),jcand.eta(),jcand.phi(),jcand.mass());
    chtrkspv_LV += chtrkspv_LV2;
   }else{
    math::PtEtaPhiMLorentzVector chtrksnpv_LV2(jcand.pt(),jcand.eta(),jcand.phi(),jcand.mass());
    chtrksnpv_LV += chtrksnpv_LV2;
    jet_chtrksnpv++;
    jetchtrksnpv.push_back(trk);
   }
  }//Ch trks 
 }//Loop on jet daus 
 jet_chtrks_pt    = chtrks_LV.pt();
 jet_chtrks_M     = chtrks_LV.M();
 jet_chtrkspv_pt  = chtrkspv_LV.pt();
 jet_chtrkspv_M   = chtrkspv_LV.M();
 jet_chtrksnpv_pt = chtrksnpv_LV.pt();
 jet_chtrksnpv_M  = chtrksnpv_LV.M();
}
bool is_loosePOG_jetmuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h){
  bool ismu = false;
  for(const pat::Muon &mu : *muon_h){
    if(deltaR(jcand.p4(),mu.p4())<0.1 && fabs(jcand.pt()-mu.pt())/mu.pt()<0.05){
     if(mu.isLooseMuon()) ismu = true;
     if(ismu) break;
    }
  }
  return ismu;
}
bool is_softLep_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<std::vector<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(std::vector<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele;
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05 ){
      const HitPattern &hitPattern = lele.gsfTrack().get()->hitPattern();
      uint32_t hit = hitPattern.getHitPattern(HitPattern::TRACK_HITS, 0);
      bool hitCondition = !(HitPattern::validHitFilter(hit) && ((HitPattern::pixelBarrelHitFilter(hit) && HitPattern::getLayer(hit) < 3) || HitPattern::pixelEndcapHitFilter(hit)));
      if(!hitCondition && lele.passConversionVeto()) isele = true;
      if(isele) break;
    }
  }
  return isele;
}
bool is_vetoPOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<std::vector<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(std::vector<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele;
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05){
      double ooEmooP = 999;
      if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
      else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
      else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy() );
      if(lele.full5x5_sigmaIetaIeta()<0.012 && fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.0126 && fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.107
         && lele.hcalOverEcal()<0.186 && ooEmooP<0.239
&& lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=2 && lele.passConversionVeto()
        ) isele = true;
      if(isele) break;
    }
  }
  return isele;
}
bool is_loosePOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<std::vector<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(std::vector<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele;
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05) {
      double ooEmooP = 999;
      if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
      else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
      else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy() );
      if(lele.full5x5_sigmaIetaIeta()<0.0105 && fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.00976 && fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.0929
         && lele.hcalOverEcal()<0.0765 && ooEmooP<0.184
 && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=2 && lele.passConversionVeto()
        ) isele = true;
      if(isele) break;
  }
  }
  return isele;
}
void get_jetlep(const pat::Jet& jet, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder,int& jet_ndaus,edm::Handle<std::vector<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,double& bjetness_num_pdgid_eles, double& bjetness_num_pdgid_mus, double& bjetness_num_soft_eles, double& bjetness_num_vetonoipnoiso_eles, double& bjetness_num_loosenoipnoiso_eles, double& bjetness_num_loose_mus){
 vector<CandidatePtr> jdaus(jet.daughterPtrVector());
 sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
 jet_ndaus = jdaus.size();
 for(int jd=0; jd<jet_ndaus; jd++){
  const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
  if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
  Track trk = Track(jcand.pseudoTrack());
  bool isgoodtrk = is_goodtrk(trk,vtx);
  if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
   if(fabs(jcand.pdgId())==13){
    bjetness_num_pdgid_mus++;
    if(is_loosePOG_jetmuon(jcand,muon_h)) bjetness_num_loose_mus++;
   }
   if(fabs(jcand.pdgId())==11){
    bjetness_num_pdgid_eles++;
    if(is_softLep_jetelectron(jcand,electron_pat,vtx)) bjetness_num_soft_eles++;
    if(is_vetoPOGNoIPNoIso_jetelectron(jcand,electron_pat,vtx)) bjetness_num_vetonoipnoiso_eles++;
    if(is_loosePOGNoIPNoIso_jetelectron(jcand,electron_pat,vtx)) bjetness_num_loosenoipnoiso_eles++;
   }
  }   
 }
}
//Get the 3D IP val, sig, sigIP 
void IPToolsValues3D(const TransientTrack ttrk, reco::Vertex vtx, GlobalVector gv3D, double& trk_IP3D_val,double& trk_IP3D_sig, double& trk_sIP3D_val,double& trk_sIP3D_sig, double& trk_IP3D_err, double& trk_sIP3D_err){
 //3D
 trk_IP3D_val = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.value();
 trk_IP3D_err = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.error();
 trk_IP3D_sig = IPTools::absoluteImpactParameter3D(ttrk,vtx).second.significance();
 //s3D
 trk_sIP3D_val = IPTools::signedImpactParameter3D(ttrk,gv3D,vtx).second.value();
 trk_sIP3D_err = IPTools::signedImpactParameter3D(ttrk,gv3D,vtx).second.error();
 trk_sIP3D_sig = IPTools::signedImpactParameter3D(ttrk,gv3D,vtx).second.significance();
}
 //Get the 2DIP value, sig, signIP
void IPToolsValues2D(const TransientTrack ttrk, reco::Vertex vtx, GlobalVector gv2D,double& trk_IP2D_val,double& trk_IP2D_sig,double& trk_sIP2D_val,double& trk_sIP2D_sig,double& trk_IP2D_err,double& trk_sIP2D_err ){
 //2D
 trk_IP2D_val = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.value();
 trk_IP2D_err = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.error();
 trk_IP2D_sig = IPTools::absoluteTransverseImpactParameter(ttrk,vtx).second.significance();
 //s2D
 trk_sIP2D_val = IPTools::signedTransverseImpactParameter(ttrk,gv2D,vtx).second.value();
 trk_sIP2D_err = IPTools::signedTransverseImpactParameter(ttrk,gv2D,vtx).second.error();
 trk_sIP2D_sig = IPTools::signedTransverseImpactParameter(ttrk,gv2D,vtx).second.significance();
}
 // Get the 1DIP value, sig and signedIP
void IPToolsValues1D(const TransientTrack ttrk, reco::Vertex vtx, GlobalVector gv, double& trk_IP1D_val,double& trk_IP1D_sig, double& trk_sIP1D_val,double& trk_sIP1D_sig,double& trk_IP1D_err, double& trk_sIP1D_err){ 
 GlobalPoint vert(vtx.position().x(), vtx.position().y(), vtx.position().z());
 TrajectoryStateClosestToPoint traj = ttrk.trajectoryStateClosestToPoint(vert);
 trk_IP1D_val = fabs(traj.perigeeParameters().longitudinalImpactParameter());
 trk_IP1D_err = traj.perigeeError().longitudinalImpactParameterError();
 trk_IP1D_sig = trk_IP1D_val/trk_IP1D_err;
 //Get the Sign
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint impactPoint = closestIn3DSpaceState.globalPosition();
 GlobalVector IPVec(0,0,impactPoint.z()-vtx.position().z());
 double prod   = IPVec.dot(gv);
 double sign   = (prod>=0) ? 1. : -1.;
 trk_sIP1D_val = sign*trk_IP1D_val;
 trk_sIP1D_err = sign*trk_IP1D_err;
 trk_sIP1D_sig = sign*trk_IP1D_sig;
}
////
//Get the IP3D and IP2D using the defination of IPTools 
////
void get_IP_TE(const TransientTrack ttrk, reco::Vertex vtx,GlobalVector gv ,double& trk_IP1D_val1, double& trk_IP1D_sig1 ,double& trk_sIP1D_val1,double& trk_sIP1D_sig1, int dimension){
 GlobalPoint vert(vtx.position().x(), vtx.position().y(), vtx.position().z());
 //AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 TransverseImpactPointExtrapolator extrapolator(ttrk.field()); 
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector IPVec(0,0,refPoint.z()-vtx.position().z());
 AlgebraicSymMatrix33 error = refPointErr.matrix()+ vertexPositionErr.matrix();
 GlobalVector diff          = refPoint - vertexPosition;
 AlgebraicVector3 vDiff;
//  AlgebraicSymMatrix33 error = refPointErr.matrix()+ vertexPositionErr.matrix();
 // GlobalVector diff          = refPoint - vertexPosition;
 double  x;
 double  y;
 double  z;
 if(dimension==1){
  x        =0;
  y        =0;
  z        =diff.z();
  vDiff[0] = 0;
  vDiff[1] = 0;
  vDiff[2] = diff.z();
 }else if(dimension==2){
   x        = diff.x();
   y        = diff.y();
   z        = 0;
   vDiff[0] = diff.x();
   vDiff[1] = diff.y();
   vDiff[2] = 0;
  }else{
    x        = diff.x();
    y        = diff.y();
    z        = diff.z();
    vDiff[0] = diff.x();
    vDiff[1] = diff.y();
    vDiff[2] = diff.z();
  }
 trk_IP1D_val1 = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
 double err2   = ROOT::Math::Similarity(error,vDiff);
 double err    = 0;
 if (trk_IP1D_val1 != 0)  err  =  sqrt(err2)/trk_IP1D_val1;
 trk_IP1D_sig1  = trk_IP1D_val1/err;
 double prod    = IPVec.dot(gv);
 double sign    = (prod>=0) ? 1. : -1.;
 trk_sIP1D_val1 = sign*trk_IP1D_val1;
 trk_sIP1D_sig1 = sign*trk_IP1D_sig1;
}

void get_IP_AE(const TransientTrack ttrk, reco::Vertex vtx,GlobalVector gv,double& trk_IP1D_val,double& trk_sIP1D_val,double& trk_IP1D_sig,double& trk_sIP1D_sig,double& trk_IP1D_err, double& trk_sIP1D_err, int dimension){
 GlobalPoint vert(vtx.position().x(), vtx.position().y(), vtx.position().z());
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector IPVec(refPoint.x()-vtx.position().x(),refPoint.y()-vtx.position().y(),refPoint.z()-vtx.position().z());
 AlgebraicSymMatrix33 error = refPointErr.matrix()+ vertexPositionErr.matrix();
 GlobalVector diff          = refPoint - vertexPosition;
 AlgebraicVector3 vDiff;
 double  x;
 double  y;
 double  z;
 if(dimension==1){
  x        = 0;
  y        = 0;
  z        = diff.z();
  vDiff[0] = 0;
  vDiff[1] = 0;
  vDiff[2] = diff.z();
 }else if(dimension==2){
   x        = diff.x();
   y        = diff.y();
   z        = 0;
   vDiff[0] = diff.x();
   vDiff[1] = diff.y();
   vDiff[2] = 0;
  }else{
    x        = diff.x();
    y        = diff.y();
    z        = diff.z();
    vDiff[0] = diff.x();
    vDiff[1] = diff.y();
    vDiff[2] = diff.z();
  }
 trk_IP1D_val = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
 double err2  = ROOT::Math::Similarity(error,vDiff);
 if (trk_IP1D_val != 0) trk_IP1D_err  =  sqrt(err2)/trk_IP1D_val;
 trk_IP1D_sig  = trk_IP1D_val/trk_IP1D_err;
 double prod   = IPVec.dot(gv);
 double sign   = (prod>=0) ? 1. : -1.;
 trk_sIP1D_val = sign*trk_IP1D_val;
 trk_sIP1D_err = sign*trk_IP1D_err;
 trk_sIP1D_sig = sign*trk_IP1D_sig;
}
void zIPValues(const TransientTrack ttrk, reco::Vertex vtx, GlobalVector gv,double& zIPvalue){
 SignedTransverseImpactParameter stip;
 zIPvalue = stip.zImpactParameter(ttrk,gv,vtx).second.value();
}
 //new function to calculate IP1D and investigations for sign
void get_IP_AE1D_x(const TransientTrack ttrk, reco::Vertex vtx,GlobalVector gv,double& trk_sIP1D_val_x,double& diff_RP_PV_x,double& prod_x){
 GlobalPoint vert(vtx.position().x(), vtx.position().y(), vtx.position().z());
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 //double diff_x,diff_y,diff_z;
 GlobalVector diff             = refPoint - vertexPosition;
 GlobalVector IPVec_x(refPoint.x()-vtx.position().x(),0,0);
 diff_RP_PV_x          = diff.x();
 double trk_IP1D_val_x = sqrt(pow(diff.x(),2));
 prod_x                = IPVec_x.dot(gv);
 double sign           = (prod_x>=0) ? 1. : -1.;
 trk_sIP1D_val_x       = sign*trk_IP1D_val_x;
} 
void get_IP_AE1D_y(const TransientTrack ttrk, reco::Vertex vtx,GlobalVector gv,double& trk_sIP1D_val_y,double& diff_RP_PV_y,double& prod_y){
 GlobalPoint vert(vtx.position().x(), vtx.position().y(), vtx.position().z());
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector diff             = refPoint - vertexPosition;
 //Calculate difference between refrence point and vtx position
 GlobalVector IPVec_y(0,refPoint.y()-vtx.position().y(),0);
 diff_RP_PV_y          = diff.y();
 double trk_IP1D_val_y = sqrt(pow(diff.y(),2));
 prod_y                = IPVec_y.dot(gv);
 double sign           = (prod_y>=0) ? 1. : -1.;
 trk_sIP1D_val_y       = sign*trk_IP1D_val_y;
}
void get_IP_AE1D_z(const TransientTrack ttrk, reco::Vertex vtx,GlobalVector gv,double& trk_sIP1D_val_z,double& diff_RP_PV_z,double& prod_z){
 GlobalPoint vert(vtx.position().x(), vtx.position().y(), vtx.position().z());
 AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
 TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(),vert);
 GlobalPoint refPoint          = tsos.globalPosition();
 GlobalError refPointErr       = tsos.cartesianError().position();
 GlobalPoint vertexPosition    = RecoVertex::convertPos(vtx.position());
 GlobalError vertexPositionErr = RecoVertex::convertError(vtx.error());
 GlobalVector diff             = refPoint - vertexPosition;
 GlobalVector IPVec_z(0,0,refPoint.z()-vtx.position().z());
 diff_RP_PV_z          = diff.z();
 double trk_IP1D_val_z = sqrt(pow(diff.z(),2));
 prod_z                = IPVec_z.dot(gv);
 double sign           = (prod_z>=0) ? 1. : -1.;
 trk_sIP1D_val_z       = sign*trk_IP1D_val_z;
}

//Get the jet tracks
void get_jettrks_PV(const pat::Jet& jet, const TransientTrackBuilder& ttrkbuilder, vector<Track>& jetchtrks_PV){
 vector<CandidatePtr> jdaus(jet.daughterPtrVector());
 sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
 int jet_ndaus = jdaus.size();
 for(int jd=0; jd<jet_ndaus; jd++){
  const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
  //if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
  Track trk = Track(jcand.pseudoTrack());
  if(jcand.charge()!=0 && jcand.numberOfHits()>0 && jcand.fromPV()==pat::PackedCandidate::PVUsedInFit) jetchtrks_PV.push_back(trk);
 }//loop on jdaus
}
//Get the jet tracks not from PV
void get_jettrks_nPV(const pat::Jet& jet, const TransientTrackBuilder& ttrkbuilder, vector<Track>& jetchtrks_nPV){
 vector<CandidatePtr> jdaus(jet.daughterPtrVector());
 sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
 int jet_ndaus = jdaus.size();
 for(int jd=0; jd<jet_ndaus; jd++){
  const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
  //if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
  Track trk = Track(jcand.pseudoTrack());
  if(jcand.charge()!=0 && jcand.numberOfHits()>0) jetchtrks_nPV.push_back(trk);
 }
}
//Measure Flight Distance 
void get_fd(reco::Vertex pv,reco::Vertex sv, GlobalVector gv,double& fd3dval, double& fd3dsig, double& fd2dval, double& fd2dsig){
 fd3dval = SecondaryVertex::computeDist3d(pv,sv,gv,true).value();
 fd3dsig = SecondaryVertex::computeDist3d(pv,sv,gv,true).significance();
 fd2dval = SecondaryVertex::computeDist2d(pv,sv,gv,true).value();
 fd2dsig = SecondaryVertex::computeDist2d(pv,sv,gv,true).significance();
}
//Measure Decay Length
void DecayLength3D(const TransientTrack ttrk,reco::Vertex vtx,GlobalVector gv,double& DL_val,double& DL_error,double& DL_sig){
 TrajectoryStateOnSurface closestToJetState = IPTools::closestApproachToJet(ttrk.impactPointState(), vtx, gv, ttrk.field());
 if(!closestToJetState.isValid()) return;
 GlobalVector jetDirection  = gv.unit();
 GlobalVector jetDirection2(jetDirection.x(),jetDirection.y(),jetDirection.z());
 GlobalPoint vertexPosition(vtx.position().x(),vtx.position().y(),vtx.position().z());
 GlobalPoint Point1(closestToJetState.globalPosition().x(),closestToJetState.globalPosition().y(),closestToJetState.globalPosition().z());
 DL_val = jetDirection2.dot(Point1-vertexPosition);
 AlgebraicVector3 j;
 j[0] = jetDirection.x();
 j[1] = jetDirection.y();
 j[2] = jetDirection.z();
 AlgebraicVector6 jj;
 jj[0] = jetDirection.x();
 jj[1] = jetDirection.y();
 jj[2] = jetDirection.z();
 jj[3] = 0.;
 jj[4] = 0.;
 jj[5] = 0.;
 reco::Vertex vtx2=vtx;
 double trackError2 = ROOT::Math::Similarity(jj,closestToJetState.cartesianError().matrix());
 double vertexError2 = ROOT::Math::Similarity(j,vtx2.covariance());
 DL_error = sqrt(trackError2+vertexError2);
 DL_sig   = DL_val/DL_error;
}
//Measure Decay Length 2D
void DecayLength2D(const TransientTrack ttrk,reco::Vertex vtx,GlobalVector gv,double& DL_val,double& DL_error,double& DL_sig){
 TrajectoryStateOnSurface closestToJetState = IPTools::closestApproachToJet(ttrk.impactPointState(), vtx, gv, ttrk.field());
 if(!closestToJetState.isValid()) return;
 GlobalVector jetDirection  = gv.unit();
 GlobalVector jetDirection2(jetDirection.x(),jetDirection.y(),0);
 GlobalPoint vertexPosition(vtx.position().x(),vtx.position().y(),0);
 GlobalPoint Point1(closestToJetState.globalPosition().x(),closestToJetState.globalPosition().y(),0);
 DL_val = jetDirection2.dot(Point1-vertexPosition);
 AlgebraicVector3 j;
 j[0] = jetDirection.x();
 j[1] = jetDirection.y();
 j[2] = 0;
 AlgebraicVector6 jj;
 jj[0] = jetDirection.x();
 jj[1] = jetDirection.y();
 jj[2] = 0;
 jj[3] = 0.;
 jj[4] = 0.;
 jj[5] = 0.;
 reco::Vertex vtx2=vtx;
 double trackError2 = ROOT::Math::Similarity(jj,closestToJetState.cartesianError().matrix());
 double vertexError2 = ROOT::Math::Similarity(j,vtx2.covariance());
 DL_error = sqrt(trackError2+vertexError2);
 DL_sig = DL_val/DL_error;
}
//Decay Length 1D
void DecayLength1D(const TransientTrack ttrk,reco::Vertex vtx,GlobalVector gv,double& DL_val,double& DL_error,double& DL_sig){
 TrajectoryStateOnSurface closestToJetState = IPTools::closestApproachToJet(ttrk.impactPointState(), vtx, gv, ttrk.field());
 if(!closestToJetState.isValid()) return;
 GlobalVector jetDirection  = gv.unit();
 GlobalVector jetDirection2(0,0,jetDirection.z());
 GlobalPoint vertexPosition(0,0,vtx.position().z());
 GlobalPoint Point1(0,0,closestToJetState.globalPosition().z());
 DL_val = jetDirection2.dot(Point1-vertexPosition);
 AlgebraicVector3 j;
 j[0] = 0;
 j[1] = 0;
 j[2] = jetDirection.z();
 AlgebraicVector6 jj;
 jj[0] = 0;
 jj[1] = 0;
 jj[2] = jetDirection.z();
 jj[3] = 0.;
 jj[4] = 0.;
 jj[5] = 0.;
 reco::Vertex vtx2=vtx;
 double trackError2 = ROOT::Math::Similarity(jj,closestToJetState.cartesianError().matrix());
 double vertexError2 = ROOT::Math::Similarity(j,vtx2.covariance());
 DL_error = sqrt(trackError2+vertexError2);
 DL_sig = DL_val/DL_error;
}


