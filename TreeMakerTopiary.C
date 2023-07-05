#define TreeMakerTopiary_cxx
#include "TreeMakerTopiary.h"
#include "RestFrames/RestFrames.hh"
#include "JetCorrectionUncertainty.h"
#include "JetCorrectorParameters.h"
#include "JetResolution.h"
#include "XYMETCorrection_withUL17andUL18andUL16.h"
#include <TLeaf.h>
#include <TVector.h>
#include <TH1.h>
#include <TProfile.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include "GEScaleSyst.cc"

using std::vector;
using std::string;

RestFrames::RFKey ensure_autoload(1);
using namespace RestFrames;

double deltaR(TLorentzVector T1, TLorentzVector T2){
  double p1 = T1.Phi(); 
  double p2 = T2.Phi(); 
  double e1 = T1.Eta(); 
  double e2 = T2.Eta(); 
  auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);
  return std::sqrt((e1-e2)*(e1-e2) + dp*dp);
}

TLorentzVector closestgenjetParticle(TLorentzVector jet, vector<TLorentzVector>* genjets){
  double deltarmin = std::numeric_limits<double>::infinity();
  TLorentzVector next;
  for(unsigned int i=0; i<genjets->size(); ++i) {
    TLorentzVector pi = genjets->at(i);
    double dr = deltaR(genjets->at(i), jet);
    //std::cout<<"dr is :  "<< dr << " Pt : Eta : Phi:  " << pi.Pt() << " : " << pi.Eta() << " : " << pi.Phi() <<std::endl;
    if(dr < deltarmin && pi != jet) {
      deltarmin = dr;
      next = pi;
    }
  }
  //std::cout<<"next pt : eta : phi :  "<< next.Pt() << " : " << next.Eta() << " : " << next.Phi() << std::endl;
  return next;
}

TVector2 getShiftedPtvec(double pt,double ptunc, double phi,int dir){
  double vecx = (pt+dir*ptunc)*std::cos(phi);//dir needs to -1 for down, +1 for up
  double vecy = (pt+dir*ptunc)*std::sin(phi);
  TVector2 vec(vecx,vecy);
  return vec;
}

float getPDFUncertainty(vector<float>* LHEPdfWeight){
  // http://nnpdf.mi.infn.it/wp-content/uploads/2019/03/NNPDFfits_Positivity_PhysicalCrossSections_v2.pdf
  std::sort(LHEPdfWeight->begin(), LHEPdfWeight->begin()+100);
  float err = (LHEPdfWeight->at(83) - LHEPdfWeight->at(16))/2;
  return err;
}

std::pair<float,float> getQCDScaleUpDwn(vector<float>* scales){
  // [mur=1, muf=1], [mur=1, muf=2], [mur=1, muf=0.5], [mur=2, muf=1], [mur=2, muf=2], [mur=2, muf=0.5], [mur=0.5, muf=1], [mur=0.5, muf=2], [mur=0.5, muf=0.5]
  //Remove the 0.5 and 2 pairs
  //Take the highest and lowest as the up and down scale
  vector<float>* redscale = scales;
  
  //for (std::size_t i = 0; i < redscale->size();++i){
  //  std::cout<<i<<" "<<redscale->at(i)<<std::endl;
  //}
  redscale->erase(std::next(redscale->begin(),5));
  redscale->erase(std::next(redscale->begin(),6));
  redscale->erase(redscale->begin());//deletes [1,1]

  //std::cout<<"The new scales"<<std::endl;
  float max = -99999;
  float min = 999999;
  for (std::size_t i = 0; i < redscale->size();++i){
    //std::cout<<i<<" "<<redscale->at(i)<<std::endl;
    if (redscale->at(i) > max){
      max = redscale->at(i);
    }
    if (redscale->at(i) < min){
      min = redscale->at(i);
    }
  }
  
  std::pair<float,float> updwn(min,max);
  //std::cout<<"The dwn/up weights"<<std::endl;
  //std::cout<<"       "<<updwn.first<<std::endl;
  //std::cout<<"       "<<updwn.second<<std::endl;
  return updwn;

}

			      

TLorentzVector doHEMshiftJet(TLorentzVector jet,bool fid){
  if ((jet.Eta() < -1.3) && (jet.Eta() > -2.5) && (jet.Phi() > -1.57) && (jet.Phi() < -0.87) && fid) {
    jet = jet*0.80;//scale down the energy by 20%
    }
  if ((jet.Eta() < -2.5) && (jet.Eta() > -3.0) && (jet.Phi() > -1.57) && (jet.Phi() < -0.87) && fid) {
    jet = jet*0.70;//scale down the energy by 20%
    }
  return jet;
}
  
std::pair<double,double> getMETCorrHEM(TLorentzVector ojet, TLorentzVector njet){
  //double  getMETCorrHEM(TLorentzVector ojet, TLorentzVector njet){
  double oldphi = ojet.Phi();
  double newphi = njet.Phi();
  double oldpt = ojet.Pt();
  double newpt = njet.Pt();
  TLorentzVector diff = ojet-njet;
  double diffx = diff.Pt()*std::cos(newphi);
  double diffy = diff.Pt()*std::sin(newphi);
  std::pair<double,double> metshift(diffx,diffy);

  /*
  std::cout<<"    Original Jet Pt "<<ojet.Pt()<<std::endl;
  std::cout<<"   Original Jet Phi "<<ojet.Phi()<<std::endl;
  std::cout<<"        new Jet Pt  "<<njet.Pt()<<std::endl;
  std::cout<<"        new Jet Phi "<<njet.Phi()<<std::endl;
  std::cout<<"       diff Jet Pt  "<<diff.Pt()<<std::endl;
  std::cout<<"       diff Jet Phi "<<diff.Phi()<<std::endl;

  std::cout<<"       diff ptx     "<<metshift.first<<std::endl;
  std::cout<<"       diff pty     "<<metshift.second<<std::endl;
  //*/  
  return metshift;
}
  


std::vector<double> GetMuonTriggerSF(std::vector<double> sfvec1,std::vector<double> sfvec2,TLorentzVector obj1,TLorentzVector obj2,float highbinedge,TH2 *heff) {
  std::vector<double> sfvec;
  double leptonsf  = 1;
  double lsfup     = 0;
  double lsfdwn    = 0;
  int leptonbin1   = -1;
  int leptonbin2   = -1;
  double ptcheck1 = obj1.Pt();
  double ptcheck2 = obj2.Pt();
  double obj1eff = -1;
  double obj2eff = -1;
  if (ptcheck1 >= highbinedge) {
    ptcheck1 =  highbinedge - 20.0;//safely within last bin, but a hack
  }
  if (ptcheck2 >= highbinedge) {
    ptcheck2 =  highbinedge - 20.0;//safely within last bin, but a hack
  }
  leptonbin1 = heff->FindBin(ptcheck1,std::fabs(obj1.Eta()));
  leptonbin2 = heff->FindBin(ptcheck2,std::fabs(obj2.Eta()));
  obj1eff = heff->GetBinContent(leptonbin1);
  obj2eff = heff->GetBinContent(leptonbin2);
  leptonsf = (1-(1-sfvec1[0]*obj1eff)*(1-sfvec2[0]*obj2eff))/(1-(1-obj1eff)*(1-obj1eff));
  //The unc will only include the sf unc, as the eff should be small
  double dsf1 = obj1eff*(1-sfvec2[0]*obj2eff)/(obj2eff+obj1eff-obj2eff*obj1eff);
  double dsf2 = obj2eff*(1-sfvec1[0]*obj1eff)/(obj2eff+obj1eff-obj2eff*obj1eff);
  lsfup = sqrt(dsf1*dsf1*sfvec1[1]*sfvec1[1]+dsf2*dsf2*sfvec2[1]*sfvec2[1]);
  lsfdwn = lsfup;
  //write the output
  sfvec.push_back(leptonsf);
  sfvec.push_back(lsfup);
  sfvec.push_back(lsfdwn);

  /*
  std::cout<<"    The scale factor of the leading muon is "<<sfvec1[0]<<std::endl;
  std::cout<<"    The scale factor of the sleading muon is "<<sfvec2[0]<<std::endl;
  std::cout<<"    The scale factor of the event  is "<<sfvec[0]<<std::endl;
  std::cout<<"    The scale unc    of the event  is "<<sfvec[1]<<std::endl;
  */
  return sfvec;
}

std::vector<double> combineTheLeptonSF(std::vector<double> sfvec1,std::vector<double> sfvec2) {
  std::vector<double> outv;
  double sf = 1;
  double upunc = 0;
  double dnunc = 0;
  sf = sfvec1[0]*sfvec2[0];
  outv.push_back(sf);
  upunc = sf*sqrt(pow(sfvec1[1],2)/sfvec1[0]+pow(sfvec2[1],2)/sfvec2[0]);
  dnunc = sf*sqrt(pow(sfvec1[2],2)/sfvec1[0]+pow(sfvec2[2],2)/sfvec2[0]);
  outv.push_back(upunc);
  outv.push_back(dnunc);
  //std::cout<<"Value of unc sfvec1:  "<<sfvec1[1]<<std::endl;
  //std::cout<<"Value of unc sfvec2:  "<<sfvec2[1]<<std::endl;
  //std::cout<<"combounc val:  "<<upunc<<std::endl;
  //std::cout<<"Value of second vec sf: "<<sfvec2[0]<<std::endl;
  //std::cout<<"Value of combined sf:   "<<sf<<std::endl;
  return outv;
}

float getMagnitude(std::vector<int> vec){
  float sum2 = 0;
  for (std::size_t i = 0;i < vec.size();++i){
    sum2 += vec[i]*vec[i];
  }
  float mag = sqrt(sum2);
  return mag;
}

std::vector<double> GetElectronPtEtaSF(int year,TH2 *hist,TLorentzVector obj, bool isTrig) {
  std::vector<double> sfvec;
  double leptonsf  = 1;
  double lsfup  = 0;
  double lsfdwn = 0;
  int leptonbin    = -1;
  double ptcheck = obj.Pt();

  if (isTrig) {
      if (year == 18  && ptcheck >= 2000.0){
	ptcheck = 1990.0;
      }
      if (year != 18 && ptcheck >= 1000.0){
	ptcheck == 990.0;
      }
    }

  else if (ptcheck >= 500.0) {
    ptcheck = 400.0;//safely within last bin, but a hack
  }

  //Gather the electron's scale factors
  leptonbin = hist->FindBin(obj.Eta(),ptcheck);
  leptonsf  = hist->GetBinContent(leptonbin);
  lsfup     = hist->GetBinErrorUp(leptonbin);
  lsfdwn    = hist->GetBinErrorLow(leptonbin);

  if (isTrig){
    lsfup = 0.005;
    lsfdwn = 0.005;
  }
  sfvec.push_back(leptonsf);
  sfvec.push_back(lsfup);
  sfvec.push_back(lsfdwn);
  return sfvec;
}

std::vector<double> GetMuonPtEtaSF(int year,TH2 *hist,TLorentzVector obj,float highbinedge, bool isID) {
  //std::cout<<"This is the hist "<<hist<<std::endl;
  std::vector<double> sfvec;
  double leptonsf  = 1;
  double lsfup  = 0;
  double lsfdwn = 0;
  int leptonbin    = -1;
  double ptcheck = obj.Pt();
  if (ptcheck >= highbinedge) {
    ptcheck =  highbinedge - 20.0;//safely within last bin, but a hack
  }
  if (year == 16 && isID) {
    leptonbin = hist->FindBin(ptcheck,obj.Eta());
  }
  else{
    leptonbin = hist->FindBin(ptcheck,std::abs(obj.Eta()));
  }
  leptonsf = hist->GetBinContent(leptonbin);
  lsfup    = hist->GetBinErrorUp(leptonbin);
  lsfdwn   = hist->GetBinErrorLow(leptonbin);
  sfvec.push_back(leptonsf);
  sfvec.push_back(lsfup);
  sfvec.push_back(lsfdwn);
  return sfvec;
}

void TreeMakerTopiary::Loop(std::string outputFileName, float totalOriginalEvents, int sampleType,int year,int anchan, TVector metsys)
//void TreeMakerTopiary::Loop(std::string outputFileName, float totalOriginalEvents, int sampleType,int year, int anchan, std::vector<int>  *metsys)
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   

   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("RunNum",1);
   fChain->SetBranchStatus("LumiBlockNum",1);
   fChain->SetBranchStatus("EvtNum",1);
   //fChain->SetBranchStatus("GenMET",1);
   //fChain->SetBranchStatus("GenMETPhi",1);
   fChain->SetBranchStatus("TriggerPass",1);
   fChain->SetBranchStatus("JetsAK8Clean*",1);
   fChain->SetBranchStatus("Muons*",1);
   //fChain->SetBranchStatus("Electrons*",1);
   fChain->SetBranchStatus("MET",1);
   fChain->SetBranchStatus("METPhi",1);
   fChain->SetBranchStatus("METUp",1);
   fChain->SetBranchStatus("METPhiUp",1);
   fChain->SetBranchStatus("METDown",1);
   fChain->SetBranchStatus("METPhiDown",1);
   fChain->SetBranchStatus("SelectedMuons",1);
   //fChain->SetBranchStatus("SelectedMuonsTunepMomenErr",1);
   fChain->SetBranchStatus("ZCandidates*",1);
   fChain->SetBranchStatus("SelectedElectrons*",1);
   fChain->SetBranchStatus("eeBadScFilter",1);
   fChain->SetBranchStatus("ZCandidatesMuMu",1);
   fChain->SetBranchStatus("ZCandidatesEE",1);
   fChain->SetBranchStatus("ZCandidatesEU",1);
   fChain->SetBranchStatus("NVtx",1);
   //fChain->SetBranchStatus("JetsAK8*",1);
   //fChain->SetBranchStatus("GenJetsAK8*",1);
   //fChain->SetBranchStatus("Jets*",1);

   if (sampleType > 0){
     fChain->SetBranchStatus("GenParticles*",1);
     fChain->SetBranchStatus("ScaleWeights",1);
     fChain->SetBranchStatus("PDFweights",1);
     //fChain->SetBranchStatus("puSysDown",1);
     //fChain->SetBranchStatus("puSysUp",1);
     //fChain->SetBranchStatus("puWeight",1);
   }

   if (year == 16 or year == 17 && sampleType > 0){
     fChain->SetBranchStatus("NonPrefiringProb*",1);
   }


   TFile qcdnnloFile("../DYCorrection/lindert_qcd_nnlo_sf.root","READ");
   TH1D *hqcdnnlosf  = (TH1D*)qcdnnloFile.Get("eej");
   hqcdnnlosf->SetDirectory(0);
   qcdnnloFile.Close();
   TFile ewknloFile("../DYCorrection/merged_kfactors_zjets.root","READ");
   TH1F *hewknlosf = (TH1F*)ewknloFile.Get("kfactor_monojet_ewk");
   hewknlosf->SetDirectory(0);
   ewknloFile.Close();


   //Initialize Stuff
   TLorentzVector hCandidate;
   TLorentzVector ZCandidate;
   double ZCandidate_pt;
   double ZCandidate_phi;
   double ZCandidate_eta;
   double ZCandidate_m;
   double ZupCandidate_pt;
   double ZCorrCandidate_pt;
   double ZupCandidate_m;
   double ZCorrCandidate_m;
   double ZdnCandidate_pt;
   double ZdnCandidate_m;
   double hCandidate_pt;
   double hCandidate_phi;
   double hCandidate_eta;
   double hCandidate_m;
   double hCandidate_sd;
   double hCandidate_dmdhbbvqcd;
   double hCandidate_dmdzbbvqcd;
   double hCandidate_dmdzhbbvqcd;
   double hCandidate_middb;
   double mEstZp;
   double mEstND;
   double mEstNS;
   double evntwhem;
   double  evntwkf;
   double  evntwbtag;
   double  evntwmusf;
   double evntwmutrig;
   double  muiduncup;
   double  muiduncdwn;
   double  mutriguncup;
   double  mutriguncdwn;
   double  evntweltrig;
   double  eltriguncup;
   double  eltriguncdwn;
   double  evntwelidsf;
   double  eliduncup;
   double  eliduncdwn;
   double  evntwelrecosf;
   double  elrecouncup;
   double  elrecouncdwn;
   double  btaguncup;
   double  btaguncdwn;
   TLorentzVector LMuCandidate;
   double LMuCandidate_pt;
   double LMuCandidate_phi;
   double LMuCandidate_eta;
   double LMuCandidate_m;
   double LMuCandidate_ptunc;
   TLorentzVector sLMuCandidate;
   double sLMuCandidate_pt;
   double sLMuCandidate_phi;
   double sLMuCandidate_eta;
   double sLMuCandidate_m;
   double sLMuCandidate_ptunc;
   TLorentzVector LEleCandidate;
   double LEleCandidate_pt;
   double LEleCandidate_phi;
   double LEleCandidate_eta;
   double LEleCandidate_m;
   TLorentzVector sLEleCandidate;
   double sLEleCandidate_pt;
   double sLEleCandidate_phi;
   double sLEleCandidate_eta;
   double sLEleCandidate_m;
   double ghCandidate_pt;
   double ghCandidate_phi;
   double ghCandidate_eta;
   double ghCandidate_m;
   double gzCandidate_pt;
   double gzCandidate_phi;
   double gzCandidate_eta;
   double gzCandidate_m;
   double channelflag;
   double metusable;
   double metxycorr;
   double metmuup;
   double metmudn;
   double metxyphicorr;
   double metphiusable;
   TLorentzVector eventleade;
   TLorentzVector eventleadmu;
   float pdfweightup;
   float pdfweightdn;
   float qcdweightup;
   float qcdweightdn;
   float puweightnom;
   float puweightup;
   float puweightdn;
   float puweightvtxup;
   float puweightvtxdn;

   //Define the skimmed skim  output file and tree
   TFile* trimFile = new TFile(outputFileName.c_str(),"recreate");
   TTree* trimTree = fChain->CloneTree(0);
   TH1F*  hnskimed = new TH1F("hnskimed","number of events at skim level",1,0,1);
   TH1F*  hnskimedup = new TH1F("hnskimedup","weighed qcd scale events up sum",1,0,1);
   TH1F*  hnskimeddwn = new TH1F("hnskimeddwn","weighed qcd scale events dwn sum",1,0,1);
   TH1F*  hnorigevnts = new TH1F("hnorigevnts","original number of events, preskim",1,0,1);
   TH1F*  htrigpass = new TH1F("htrigpass","trigger pass",1,0,1);
   TH1F*  hZpass = new TH1F("hZpass","Z pass",1,0,1);
   TH1F*  hHpass = new TH1F("hHpass","h pass",1,0,1);
   TH1F*  hpass  = new TH1F("hpass","passing all req",1,0,1);
   TH1F*  hnpu = new TH1F("hnpu","number pileup weighted",1,0,1);
   TH1F*  hnpuup = new TH1F("hnpuup","number pileup weighted up",1,0,1);
   TH1F*  hnpudwn = new TH1F("hnpudwn","number pileup weighted dwn",1,0,1);
   TH1F*  hnpunumup = new TH1F("hnpunumup","number pileup weighted nvtx up",1,0,1);
   TH1F*  hnpunumdwn = new TH1F("hnpunumdwn","number pileup weighted nvtx dwn",1,0,1);


   //Add these horrible plots for efficiency
   TH1F* hzpasstrig_pt = new TH1F("hzpasstrig_pt","Z pt > 100 and passing triggers",74,60,800);
   TH1F* hzbuild_pt = new TH1F("hzbuild_pt","Z pt > 100",74,60,800);
   TH1F* hzpasstrig_eta = new TH1F("hzpasstrig_eta","Z pt > 100 and passing triggers",24,-2.4,2.4);
   TH1F* hzbuild_eta = new TH1F("hzbuild_eta","Z pt > 100",24,-2.4,2.4);
   
   TH1F* hlmupasstrig_pt = new TH1F("hlmupasstrig_pt","Z pt > 100 and passing triggers",74,60,800);
   TH1F* hlmubuild_pt = new TH1F("hlmubuild_pt","Z pt > 100",74,60,800);
   TH1F* hlmupasstrig_eta = new TH1F("hlmupasstrig_eta","Z pt > 100and passing triggers",24,-2.4,2.4);
   TH1F* hlmubuild_eta = new TH1F("hlmubuild_eta","Z pt > 100",24,-2.4,2.4);

   TH1F* hjetpt = new TH1F("hjetpt","All Jet Pt",48,0,1200);
   TH1F* hjeteta = new TH1F("hjeteta","All Jet eta",30,-2.8,2.8);
   TH1F* hjetphi = new TH1F("hjetphi","All Jet phi",30,-3.14159,3.14159);
   TH1F* hjetptJER = new TH1F("hjetptJER","All Jet Pt After JER",48,0,1200);
   TH1F* hjetetaJER = new TH1F("hjetetaJER","All Jet eta after JER",30,-2.8,2.8);
   TH1F* hjetphiJER = new TH1F("hjetphiJER","All Jet phi after JER",30,-3.14159,3.14159);
   
   TBranch *hCand     = trimTree->Branch("hCandidate","TLorentzVector",&hCandidate);
   TBranch *hCand_pt  = trimTree->Branch("hCandidate_pt",&hCandidate_pt,"hCandidate_pt/D");
   TBranch *hCand_phi = trimTree->Branch("hCandidate_phi",&hCandidate_phi,"hCandidate_phi/D");
   TBranch *hCand_eta = trimTree->Branch("hCandidate_eta",&hCandidate_eta,"hCandidate_eta/D");
   TBranch *hCand_m   = trimTree->Branch("hCandidate_m",&hCandidate_m,"hCandidate_m/D");
   TBranch *hCand_sd  = trimTree->Branch("hCandidate_sd",&hCandidate_sd,"hCandidate_sd/D");
   TBranch *hCand_dmdhbbvqcd  = trimTree->Branch("hCandidate_DeepMassDecorrelTagHbbvsQCD",&hCandidate_dmdhbbvqcd,"hCandidate_dmdhbbvqcd/D");
   TBranch *hCand_dmdzbbvqcd  = trimTree->Branch("hCandidate_DeepMassDecorrelTagZbbvsQCD",&hCandidate_dmdzbbvqcd,"hCandidate_dmdzbbvqcd/D");
   TBranch *hCand_dmdzhbbvqcd = trimTree->Branch("hCandidate_DeepMassDecorrelTagZHbbvsQCD",&hCandidate_dmdzhbbvqcd,"hCandidate_dmdzhbbvqcd/D");
   TBranch *hCand_middb       = trimTree->Branch("hCandidate_pfMassIndependentDeepDoubleBvLJetTagsProbHbb",&hCandidate_middb,"hCandidate_middb/D");
   TBranch *ZCand     = trimTree->Branch("ZCandidate","TLorentzVector",&ZCandidate);
   TBranch *ZCand_pt  = trimTree->Branch("ZCandidate_pt",&ZCandidate_pt,"ZCandidate_pt/D");
   TBranch *ZCorrCand_pt  = trimTree->Branch("ZCorrCandidate_pt",&ZCorrCandidate_pt,"ZCorrCandidate_pt/D");
   TBranch *ZupCand_pt  = trimTree->Branch("ZupCandidate_pt",&ZupCandidate_pt,"ZupCandidate_pt/D");
   TBranch *ZdnCand_pt  = trimTree->Branch("ZdnCandidate_pt",&ZdnCandidate_pt,"ZdnCandidate_pt/D");
   TBranch *ZCand_phi = trimTree->Branch("ZCandidate_phi",&ZCandidate_phi,"ZCandidate_phi/D");
   TBranch *ZCand_eta = trimTree->Branch("ZCandidate_eta",&ZCandidate_eta,"ZCandidate_eta/D");
   TBranch *ZCand_m   = trimTree->Branch("ZCandidate_m",&ZCandidate_m,"ZCandidate_m/D");
   TBranch *ZdnCand_m   = trimTree->Branch("ZdnCandidate_m",&ZdnCandidate_m,"ZdnCandidate_m/D");
   TBranch *ZCorrCand_m   = trimTree->Branch("ZCorrCandidate_m",&ZCorrCandidate_m,"ZCorrCandidate_m/D");
   TBranch *ZupCand_m   = trimTree->Branch("ZupCandidate_m",&ZupCandidate_m,"ZupCandidate_m/D");
   TBranch *ZpMest    = trimTree->Branch("ZPrime_mass_est",&mEstZp,"mEstZp/D");
   TBranch *NDMest    = trimTree->Branch("ND_mass_est",&mEstND,"mEstND/D");
   TBranch *NSMest    = trimTree->Branch("NS_mass_est",&mEstNS,"mEstNS/D");
   TBranch *evntweighthem = trimTree->Branch("event_weight_hem",&evntwhem,"evntwhem/D");
   TBranch *evntweightkf = trimTree->Branch("event_weight_kf",&evntwkf,"evntwkf/D");
   TBranch *evntweightbtag = trimTree->Branch("event_weight_btag",&evntwbtag,"evntwbtag/D");
   TBranch *evntbtaguncup = trimTree->Branch("event_weight_btaguncup",&btaguncup,"btaguncup/D");
   TBranch *evntbtaguncdwn = trimTree->Branch("event_weight_btaguncdwn",&btaguncdwn,"btaguncdwn/D");
   TBranch *evntweightmuid = trimTree->Branch("event_weight_muid",&evntwmusf,"evntwmusf/D");
   TBranch *evntweightmutrig = trimTree->Branch("event_weight_mutrig",&evntwmutrig,"evntwmutrig/D");
   TBranch *evntweighteltrig = trimTree->Branch("event_weight_eltrig",&evntweltrig,"evntweltrig/D");
   TBranch *evntmuiduncup = trimTree->Branch("event_weight_muiduncup",&muiduncup,"muiduncup/D");
   TBranch *evntmuiduncdwn = trimTree->Branch("event_weight_muiduncdwn",&muiduncdwn,"muiduncdwn/D");
   TBranch *evntmutriguncup = trimTree->Branch("event_weight_mutriguncup",&mutriguncup,"mutriguncup/D");
   TBranch *evntmutriguncdwn = trimTree->Branch("event_weight_mutriguncdwn",&mutriguncdwn,"mutriguncdwn/D");
   TBranch *evnteltriguncup = trimTree->Branch("event_weight_eltriguncup",&eltriguncup,"eltriguncup/D");
   TBranch *evnteltriguncdwn = trimTree->Branch("event_weight_eltriguncdwn",&eltriguncdwn,"eltriguncdwn/D");
   TBranch *evntweightelid = trimTree->Branch("event_weight_elid",&evntwelidsf,"evntwelidsf/D");
   TBranch *evnteliduncup = trimTree->Branch("event_weight_eliduncup",&eliduncup,"eliduncup/D");
   TBranch *evnteliduncdwn = trimTree->Branch("event_weight_eliduncdwn",&eliduncdwn,"eliduncdwn/D");
   TBranch *evntweightelreco = trimTree->Branch("event_weight_elreco",&evntwelrecosf,"evntwelrecosf/D");
   TBranch *evntelrecouncup = trimTree->Branch("event_weight_elrecouncup",&elrecouncup,"elrecouncup/D");
   TBranch *evntelrecouncdwn = trimTree->Branch("event_weight_elrecouncdwn",&elrecouncdwn,"elrecouncdwn/D");
   TBranch *qcdweightsdwn = trimTree->Branch("qcdweight_dwn",&qcdweightdn,"qcdweightdn/F");
   TBranch *qcdweightsup = trimTree->Branch("qcdweight_up",&qcdweightup,"qcdweightup/F");
   TBranch *pdfweightsdwn = trimTree->Branch("pdfweight_dwn",&pdfweightdn,"pdfweightdn/F");
   TBranch *pdfweightsup = trimTree->Branch("pdfweight_up",&pdfweightup,"pdfweightup/F");
   TBranch *puweights    = trimTree->Branch("puweight",&puweightnom,"puweightnom/F");
   TBranch *puweightsup    = trimTree->Branch("puweight_up",&puweightup,"puweightup/F");
   TBranch *puweightsdn    = trimTree->Branch("puweight_dwn",&puweightdn,"puweightdn/F");
   TBranch *puweightsvtxup    = trimTree->Branch("puweightvtx_up",&puweightvtxup,"puweightvtxup/F");
   TBranch *puweightsvtxdn    = trimTree->Branch("puweightvtx_dwn",&puweightvtxdn,"puweightvtxdn/F");
   TBranch *channelf   = trimTree->Branch("channel_flag",&channelflag,"channelflag/D");
   TBranch *LMuCand     = trimTree->Branch("LMuCandidate","TLorentzVector",&LMuCandidate);
   TBranch *LMuCand_pt  = trimTree->Branch("LMuCandidate_pt",&LMuCandidate_pt,"LMuCandidate_pt/D");
   TBranch *LMuCand_ptunc  = trimTree->Branch("LMuCandidate_ptunc",&LMuCandidate_ptunc,"LMuCandidate_ptunc/D");
   TBranch *LMuCand_phi = trimTree->Branch("LMuCandidate_phi",&LMuCandidate_phi,"LMuCandidate_phi/D");
   TBranch *LMuCand_eta = trimTree->Branch("LMuCandidate_eta",&LMuCandidate_eta,"LMuCandidate_eta/D");
   TBranch *sLMuCand     = trimTree->Branch("sLMuCandidate","TLorentzVector",&sLMuCandidate);
   TBranch *sLMuCand_pt  = trimTree->Branch("sLMuCandidate_pt",&sLMuCandidate_pt,"sLMuCandidate_pt/D");
   TBranch *sLMuCand_ptunc  = trimTree->Branch("sLMuCandidate_ptunc",&sLMuCandidate_ptunc,"sLMuCandidate_ptunc/D");
   TBranch *sLMuCand_phi = trimTree->Branch("sLMuCandidate_phi",&sLMuCandidate_phi,"sLMuCandidate_phi/D");
   TBranch *sLMuCand_eta = trimTree->Branch("sLMuCandidate_eta",&sLMuCandidate_eta,"sLMuCandidate_eta/D");
   TBranch *LEleCand     = trimTree->Branch("LEleCandidate","TLorentzVector",&LEleCandidate);
   TBranch *LEleCand_pt  = trimTree->Branch("LEleCandidate_pt",&LEleCandidate_pt,"LEleCandidate_pt/D");
   TBranch *LEleCand_phi = trimTree->Branch("LEleCandidate_phi",&LEleCandidate_phi,"LEleCandidate_phi/D");
   TBranch *LEleCand_eta = trimTree->Branch("LEleCandidate_eta",&LEleCandidate_eta,"LEleCandidate_eta/D");
   TBranch *sLEleCand     = trimTree->Branch("sLEleCandidate","TLorentzVector",&sLEleCandidate);
   TBranch *sLEleCand_pt  = trimTree->Branch("sLEleCandidate_pt",&sLEleCandidate_pt,"sLEleCandidate_pt/D");
   TBranch *sLEleCand_phi = trimTree->Branch("sLEleCandidate_phi",&sLEleCandidate_phi,"sLEleCandidate_phi/D");
   TBranch *sLEleCand_eta = trimTree->Branch("sLEleCandidate_eta",&sLEleCandidate_eta,"sLEleCandidate_eta/D");
   TBranch *ghCand_pt  = trimTree->Branch("ghCandidate_pt",&ghCandidate_pt,"ghCandidate_pt/D");
   TBranch *ghCand_phi = trimTree->Branch("ghCandidate_phi",&ghCandidate_phi,"ghCandidate_phi/D");
   TBranch *ghCand_eta = trimTree->Branch("ghCandidate_eta",&ghCandidate_eta,"ghCandidate_eta/D");
   TBranch *ghCand_m   = trimTree->Branch("ghCandidate_m",&ghCandidate_m,"ghCandidate_m/D");
   TBranch *gzCand_pt  = trimTree->Branch("gzCandidate_pt",&gzCandidate_pt,"gzCandidate_pt/D");
   TBranch *gzCand_phi = trimTree->Branch("gzCandidate_phi",&gzCandidate_phi,"gzCandidate_phi/D");
   TBranch *gzCand_eta = trimTree->Branch("gzCandidate_eta",&gzCandidate_eta,"gzCandidate_eta/D");
   TBranch *gzCand_m   = trimTree->Branch("gzCandidate_m",&gzCandidate_m,"gzCandidate_m/D");
   TBranch *metu       = trimTree->Branch("metsuable",&metusable,"metusable/D");
   TBranch *metu_phi   = trimTree->Branch("metphiusable",&metphiusable,"metphiusable/D");
   TBranch *metuxy       = trimTree->Branch("metxycorr",&metxycorr,"metxycorr/D");
   TBranch *metmuupb       = trimTree->Branch("metmuup",&metmuup,"metmuup/D");
   TBranch *metmudnb       = trimTree->Branch("metmudn",&metmudn,"metmudn/D");
   TBranch *metuxy_phi   = trimTree->Branch("metxyphicorr",&metxyphicorr,"metxyphicorr/D");

   //Uncertainty + Scale Factor Stuff
   //The index of the nonzero place is the ssytematic explored
   //The sign of it is up/down
   int systs = metsys.GetNoElements();
   int systidx = -999;
   if (metsys.Norm1() == 1.0) {
     for (int i = 0;i < systs; ++i) {
       if (metsys[i] != 0) {
	 systidx = i;
	 break;
       }
     }
   std::cout<<"The index for systematics is: "<<systidx<<std::endl;
   std::cout<<"The value we are applying is: "<<metsys[systidx]<<std::endl;
   }
   
   //
   int jecsys = metsys[1];//metsys is a 6D vector. JER, JES ...
   int jeRsys = metsys[0];
   std::cout<<"The value JER we are applying is: "<< jeRsys <<std::endl;


   
   TString uncfile = "badstring";
   TString btagfile = "badstring";
   TString muonsffile = "badstring";
   TString muonsfhname = "badstring";
   TString muontrigsffile = "badstring";
   TString muontrigsfhname = "badstring";
   TString muontrigefhname = "badstring";
   TString electronsffile = "badstring";
   TString electronsfhname = "badstring";
   TString electronIDsffile = "badstring";
   TString electronIDsfhname = "badstring";
   TString eltrigsffile = "badstring";
   TString eltrigsfhname = "badstring";
   TString yearstr = "badstring";
   TString unjerfile = "badstring";
   TString unjerfilesf = "badstring";
   TString pufileName = "badstring";

   TH1F *hbtagsf = 0;
   TH1F *hbtagsfuncup = 0;
   TH1F *hbtagsfuncdwn = 0;
   TH2D *hmuonsf = 0;
   TH2F *helectronsf = 0;
   TH2F *helectronIDsf = 0;
   TH2D *helectronTRIGsf = 0;
   TH2F *hmuontrigsf = 0;
   TH2F *hmuontrigef = 0;
   TH1F *hpuweight = 0;
   TH1F *hpuweightup = 0;
   TH1F *hpuweightdwn = 0;
   TProfile *hnvtxtrans = 0;
   
   if (sampleType == 3) {//ttbar
     btagfile = "btagsf/DeepAK8MassDecorrelZHbbvQCD_ttscalefactors_Zptcut_Hptcut_metcut_btagwp.root";
   }
   else {
     btagfile = "btagsf/DeepAK8MassDecorrelZHbbvQCD_scalefactors_Zptcut_Hptcut_metcut_btagwp.root";
   }

   TFile btagsf(btagfile,"READ");
   
   std::cout<<"Looking at year  "<<year<<std::endl;
   if (year == 180){
     std::cout<<"In Run2018RunA"<<std::endl;
     yearstr = "2018";
     uncfile = "JEC/Autumn18_RunA_V19_DATA/Autumn18_RunA_V19_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Autumn18_V7b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Autumn18_V7b_DATA_SF_AK8PFPuppi.txt";
   }
   if (year == 181){
     yearstr = "2018";
     uncfile = "JEC/Autumn18_RunB_V19_DATA/Autumn18_RunB_V19_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Autumn18_V7b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Autumn18_V7b_DATA_SF_AK8PFPuppi.txt";
   }
   if (year == 182){
     yearstr = "2018";
     uncfile = "JEC/Autumn18_RunC_V19_DATA/Autumn18_RunC_V19_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Autumn18_V7b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Autumn18_V7b_DATA_SF_AK8PFPuppi.txt";
   }
   if (year == 183){
     yearstr = "2018";
     uncfile = "JEC/Autumn18_RunD_V19_DATA/Autumn18_RunD_V19_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Autumn18_V7b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Autumn18_V7b_DATA_SF_AK8PFPuppi.txt";
   }
   if (year == 18) {
     yearstr = "2018";
     //JEC
     std::cout<<"In Autumn18"<<std::endl;
     uncfile = "JEC/Autumn18_V19_MC_UncertaintySources_AK4PFPuppi.txt";
     //JER
     unjerfile = "JER/Autumn18_V7b_MC_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Autumn18_V7b_MC_SF_AK8PFPuppi.txt";
     //Muon ID SF
     muonsffile = "leptonsf/Run2018ABCD_muon_SF_ID.root";
     muonsfhname = "NUM_TightID_DEN_TrackerMuons_pt_abseta";
     TFile muonsff(muonsffile,"READ");
     hmuonsf = new TH2D(*((TH2D*)muonsff.Get(muonsfhname)));
     hmuonsf->SetDirectory(0);
     muonsff.Close();
     //Muon Trigger SF
     //std::cout<<"About to try and do the muon trigger sf"<<std::endl;
     muontrigsffile = "leptonsf/Run2018_muon_SF_TRIGGER.root";
     muontrigsfhname = "NUM_Mu50_TkMu100_DEN_TightID_abseta_pt";
     muontrigefhname = "NUM_Mu50_TkMu100_DEN_TightID_abseta_pt_efficiencyMC";
     TFile muontrigsff(muontrigsffile,"READ");
     hmuontrigsf = new TH2F(*((TH2F*)muontrigsff.Get(muontrigsfhname)));
     hmuontrigef = new TH2F(*((TH2F*)muontrigsff.Get(muontrigefhname)));
     hmuontrigsf->SetDirectory(0);
     hmuontrigef->SetDirectory(0);
     muontrigsff.Close();
     //Electron reco SF
     electronsffile = "leptonsf/Run2018_electron_SF.root";
     electronsfhname = "EGamma_SF2D";
     TFile elecsff(electronsffile,"READ");
     helectronsf = new TH2F(*((TH2F*)elecsff.Get(electronsfhname)));
     helectronsf->SetDirectory(0);
     elecsff.Close();
     //Electron ID SF
     electronIDsffile = "leptonsf/Run2018_electron_SFID.root";
     electronIDsfhname = "EGamma_SF2D";
     TFile elecIDsff(electronIDsffile,"READ");
     helectronIDsf = new TH2F(*((TH2F*)elecIDsff.Get(electronIDsfhname)));
     helectronIDsf->SetDirectory(0);
     elecIDsff.Close();
     //Electron Tigger SF
     eltrigsffile = "leptonsf/Ele115orEleIso32orPho200_SF_2018.root";
     eltrigsfhname = "SF_TH2F";
     TFile elecTRIGsff(eltrigsffile,"READ");
     helectronTRIGsf = new TH2D(*((TH2D*)elecTRIGsff.Get(eltrigsfhname)));
     helectronTRIGsf->SetDirectory(0);
     elecTRIGsff.Close();
     //Btag SF
     hbtagsf = new TH1F(*((TH1F*)btagsf.Get("2018sf")));
     hbtagsfuncup = new TH1F(*((TH1F*)btagsf.Get("2018uncUp")));
     hbtagsfuncdwn = new TH1F(*((TH1F*)btagsf.Get("2018uncDown")));
     hbtagsfuncup->SetDirectory(0);
     hbtagsfuncdwn->SetDirectory(0);
     hbtagsf->SetDirectory(0);
     btagsf.Close();
     //PileUp ReWeighting
     pufileName = "pilueuphists/PileupHistograms_2018_69mb_pm5.root";
     TFile pufile(pufileName,"READ");
     hpuweight = new TH1F(*((TH1F*)pufile.Get("pu_weights_central")));
     hpuweightup = new TH1F(*((TH1F*)pufile.Get("pu_weights_up")));
     hpuweightdwn = new TH1F(*((TH1F*)pufile.Get("pu_weights_down")));
     hpuweight->SetDirectory(0);
     hpuweightup->SetDirectory(0);
     hpuweightdwn->SetDirectory(0);
     pufile.Close();
     TFile numintfile("pilueuphists/Run2_161718_pileUpTransferFunctions_Averages_Zptcut_Hptcut_metcut_btagwp.root");
     hnvtxtrans = new TProfile(*((TProfile*)numintfile.Get("tprofile_allmc_2018")));
     hnvtxtrans->SetDirectory(0);
     numintfile.Close();
   }
   if (year == 17) {
     yearstr = "2017";
     //JEC
     uncfile = "JEC/Fall17_17Nov2017_V32_MC.tar-1/Fall17_17Nov2017_V32_MC_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Fall17_V3b_MC_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Fall17_V3b_MC_SF_AK8PFPuppi.txt";
     //Muon ID SF
     muonsffile = "leptonsf/Run2017BCDEF_muon_SF_ID.root";
     muonsfhname = "NUM_TightID_DEN_genTracks_pt_abseta";
     TFile muonsff(muonsffile,"READ");
     hmuonsf = new TH2D(*((TH2D*)muonsff.Get(muonsfhname)));
     hmuonsf->SetDirectory(0);
     muonsff.Close();
     //Muon Trigger SF
     std::cout<<"About to try and do the muon trigger sf"<<std::endl;
     muontrigsffile = "leptonsf/Run2017_muon_SF_TRIGGER.root";
     muontrigsfhname = "NUM_Mu50_TkMu100_DEN_TightID_abseta_pt";
     muontrigefhname = "NUM_Mu50_TkMu100_DEN_TightID_abseta_pt_efficiencyMC";
     TFile muontrigsff(muontrigsffile,"READ");
     hmuontrigsf = new TH2F(*((TH2F*)muontrigsff.Get(muontrigsfhname)));
     hmuontrigef = new TH2F(*((TH2F*)muontrigsff.Get(muontrigefhname)));
     hmuontrigsf->SetDirectory(0);
     hmuontrigef->SetDirectory(0);
     muontrigsff.Close();
     //Electron reco SF
     electronsffile = "leptonsf/Run2017BCDEF_electron_SF.root";
     electronsfhname = "EGamma_SF2D";
     TFile elecsff(electronsffile,"READ");
     helectronsf = new TH2F(*((TH2F*)elecsff.Get(electronsfhname)));
     helectronsf->SetDirectory(0);
     elecsff.Close();
     //Electron ID SF
     electronIDsffile = "leptonsf/Run2017_electron_SFID.root";
     electronIDsfhname = "EGamma_SF2D";
     TFile elecIDsff(electronIDsffile,"READ");
     helectronIDsf = new TH2F(*((TH2F*)elecIDsff.Get(electronIDsfhname)));
     helectronIDsf->SetDirectory(0);
     elecIDsff.Close();
     //Electron Tigger SF
     eltrigsffile = "leptonsf/Ele115orEleIso35orPho200_SF_2017.root";
     eltrigsfhname = "SF_TH2F";
     TFile elecTRIGsff(eltrigsffile,"READ");
     helectronTRIGsf = new TH2D(*((TH2D*)elecTRIGsff.Get(eltrigsfhname)));
     helectronTRIGsf->SetDirectory(0);
     elecTRIGsff.Close();
     //Btag SF
     hbtagsf = new TH1F(*((TH1F*)btagsf.Get("2017sf")));
     hbtagsfuncup = new TH1F(*((TH1F*)btagsf.Get("2017uncUp")));
     hbtagsfuncdwn = new TH1F(*((TH1F*)btagsf.Get("2017uncDown")));
     hbtagsfuncup->SetDirectory(0);
     hbtagsfuncdwn->SetDirectory(0);
     hbtagsf->SetDirectory(0);
     hbtagsf->SetDirectory(0);
     btagsf.Close();
     //PileUp ReWeighting
     pufileName = "pilueuphists/PileupHistograms_2017_69mb_pm5.root";
     TFile pufile(pufileName,"READ");
     hpuweight = new TH1F(*((TH1F*)pufile.Get("pu_weights_central")));
     hpuweightup = new TH1F(*((TH1F*)pufile.Get("pu_weights_up")));
     hpuweightdwn = new TH1F(*((TH1F*)pufile.Get("pu_weights_down")));
     hpuweight->SetDirectory(0);
     hpuweightup->SetDirectory(0);
     hpuweightdwn->SetDirectory(0);
     pufile.Close();
     TFile numintfile("pilueuphists/Run2_161718_pileUpTransferFunctions_Averages_Zptcut_Hptcut_metcut_btagwp.root");
     hnvtxtrans = new TProfile(*((TProfile*)numintfile.Get("tprofile_allmc_2017")));
     hnvtxtrans->SetDirectory(0);
     numintfile.Close();

   }
   if (year == 170) {
     yearstr = "2017";
     uncfile = "JEC/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Fall17_V3b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Fall17_V3b_DATA_SF_AK8PFPuppi.txt";
       }
   if (year == 171) {
     yearstr = "2017";
     uncfile = "JEC/Fall17_17Nov2017C_V32_DATA/Fall17_17Nov2017C_V32_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Fall17_V3b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Fall17_V3b_DATA_SF_AK8PFPuppi.txt";
   }
   if (year == 172) {
     yearstr = "2017";
     uncfile = "JEC/Fall17_17Nov2017DE_V32_DATA/Fall17_17Nov2017DE_V32_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Fall17_V3b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Fall17_V3b_DATA_SF_AK8PFPuppi.txt";
   }
   if (year == 173) {
     yearstr = "2017";
     uncfile = "JEC/Fall17_17Nov2017F_V32_DATA/Fall17_17Nov2017F_V32_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Fall17_V3b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Fall17_V3b_DATA_SF_AK8PFPuppi.txt";
   }
   if (year == 16) {
     yearstr = "2016";
     //JEC
     uncfile = "JEC/Summer16_07Aug2017_V11_MC_UncertaintySources_AK8PFPuppi.txt";
     //JER
     unjerfile = "JER/Summer16_25nsV1b_MC_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Summer16_25nsV1b_MC_SF_AK8PFPuppi.txt";
     //Muon ID SF
     muonsffile = "leptonsf/Run2016_muon_SF_ID.root";
     muonsfhname = "hflip";
     TFile muonsff(muonsffile,"READ");
     hmuonsf = new TH2D(*((TH2D*)muonsff.Get(muonsfhname)));
     hmuonsf->SetDirectory(0);
     muonsff.Close();
     //Muon Trigger SF
     std::cout<<"About to try and do the muon trigger sf"<<std::endl;
     muontrigsffile = "leptonsf/Run2016_muon_SF_TRIGGER.root";
     muontrigsfhname = "NUM_Mu50_TkMu50_DEN_TightID_abseta_pt";
     muontrigefhname = "NUM_Mu50_TkMu50_DEN_TightID_abseta_pt_efficiencyMC";
     TFile muontrigsff(muontrigsffile,"READ");
     hmuontrigsf = new TH2F(*((TH2F*)muontrigsff.Get(muontrigsfhname)));
     hmuontrigef = new TH2F(*((TH2F*)muontrigsff.Get(muontrigefhname)));
     hmuontrigsf->SetDirectory(0);
     hmuontrigef->SetDirectory(0);
     muontrigsff.Close();
     //Electron reco SF
     electronsffile = "leptonsf/Run2016BCDEFGH_electron_SF.root";
     electronsfhname = "EGamma_SF2D";
     TFile elecsff(electronsffile,"READ");
     helectronsf = new TH2F(*((TH2F*)elecsff.Get(electronsfhname)));
     helectronsf->SetDirectory(0);
     elecsff.Close();
     //Electron ID SF
     electronIDsffile = "leptonsf/Run2016_electron_SFID.root";
     electronIDsfhname = "EGamma_SF2D";
     TFile elecIDsff(electronIDsffile,"READ");
     helectronIDsf = new TH2F(*((TH2F*)elecIDsff.Get(electronIDsfhname)));
     helectronIDsf->SetDirectory(0);
     elecIDsff.Close();
     //Electron Tigger SF
     eltrigsffile = "leptonsf/Ele115orEleIso27orPho175_SF_2016.root";
     eltrigsfhname = "SF_TH2F";
     TFile elecTRIGsff(eltrigsffile,"READ");
     helectronTRIGsf = new TH2D(*((TH2D*)elecTRIGsff.Get(eltrigsfhname)));
     helectronTRIGsf->SetDirectory(0);
     elecTRIGsff.Close();
     //Btag SF
     hbtagsf = new TH1F(*((TH1F*)btagsf.Get("2016sf")));
     hbtagsfuncup = new TH1F(*((TH1F*)btagsf.Get("2016uncUp")));
     hbtagsfuncdwn = new TH1F(*((TH1F*)btagsf.Get("2016uncDown")));
     hbtagsfuncup->SetDirectory(0);
     hbtagsfuncdwn->SetDirectory(0);
     hbtagsf->SetDirectory(0);
     hbtagsf->SetDirectory(0);
     btagsf.Close();
     //PileUp ReWeighting
     pufileName = "pilueuphists/PileupHistograms_2016_69mb_pm5.root";
     TFile pufile(pufileName,"READ");
     hpuweight = new TH1F(*((TH1F*)pufile.Get("pu_weights_central")));
     hpuweightup = new TH1F(*((TH1F*)pufile.Get("pu_weights_up")));
     hpuweightdwn = new TH1F(*((TH1F*)pufile.Get("pu_weights_down")));
     hpuweight->SetDirectory(0);
     hpuweightup->SetDirectory(0);
     hpuweightdwn->SetDirectory(0);
     pufile.Close();
     TFile numintfile("pilueuphists/Run2_161718_pileUpTransferFunctions_Averages_Zptcut_Hptcut_metcut_btagwp.root");
     hnvtxtrans = new TProfile(*((TProfile*)numintfile.Get("tprofile_allmc_2016")));
     hnvtxtrans->SetDirectory(0);
     numintfile.Close();
   }
   if (year == 160) {
     yearstr = "2016";
     uncfile = "JEC/Summer16_07Aug2017BCD_V11_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Summer16_25nsV1b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Summer16_25nsV1b_DATA_SF_AK8PFPuppi.txt";
   }
   if (year == 161) {
     yearstr = "2016";
     uncfile = "JEC/Summer16_07Aug2017EF_V11_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Summer16_25nsV1b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Summer16_25nsV1b_DATA_SF_AK8PFPuppi.txt";
   }
   if (year == 162) {
     yearstr = "2016";
     uncfile = "JEC/Summer16_07Aug2017GH_V11_DATA_UncertaintySources_AK8PFPuppi.txt";
     unjerfile = "JER/Summer16_25nsV1b_DATA_PtResolution_AK8PFPuppi.txt";
     unjerfilesf = "JER/Summer16_25nsV1b_DATA_SF_AK8PFPuppi.txt";
   }
   std::cout<<"Using JEC uncertainty file "<<uncfile<<std::endl;
   std::cout<<"Using JER uncertainty file "<<unjerfile<<std::endl;
   std::cout<<"Using JER SF  file "<<unjerfilesf<<std::endl;
   std::cout<<"Using btagsf file  "<<btagfile<<std::endl;
   std::cout<<"Using muon ID sf file  "<<muonsffile<<std::endl;
   std::cout<<"Using muon trigger sf file  "<<muontrigsffile<<std::endl;
   std::cout<<"Using electron sf file  "<<electronsffile<<std::endl;

   JetCorrectionUncertainty* jec_unc = new JetCorrectionUncertainty(*(new JetCorrectorParameters(uncfile.Data(),"Total")));
   // JER from the txt files : https://github.com/cms-jet/JRDatabase/tree/master/textFiles
   JME::JetResolution resolution_ = JME::JetResolution(unjerfile.Data());
   JME::JetResolutionScaleFactor res_sf_ = JME::JetResolutionScaleFactor(unjerfilesf.Data());

   //bring data era encoding back to normal numbers
   //Data era added as a digit on the end of the year
   //dividing by ten makes the year the integer when rounded
   if (year > 100){
     year = year/10;
     year = round(year);
   }

   //Info and holders
   hnskimed->SetBinContent(1,nentries);
   hnorigevnts->SetBinContent(1,totalOriginalEvents);
   float zmwinlow = 70.;
   float zmwinhi  = 110.;
   float hptcut   = 250.;

   ///*
   //Recursive Jigsaw Part
   LabRecoFrame         LABcontra("LABcontra","LABcontra");
   DecayRecoFrame       Zp("Zp","Z'");
   DecayRecoFrame       ND("ND","N_{D}");
   DecayRecoFrame       NDbar("NDbar","N_{Dbar}");
   VisibleRecoFrame     Z("Z","Z");
   InvisibleRecoFrame   NS("NS","N_{S}");
   VisibleRecoFrame     h("h","h");
   InvisibleRecoFrame   NSbar("NSbar","Z_{Sbar}");

   LABcontra.SetChildFrame(Zp);
   Zp.AddChildFrame(ND);
   Zp.AddChildFrame(NDbar);
   ND.AddChildFrame(Z);
   ND.AddChildFrame(NS);
   NDbar.AddChildFrame(h);
   NDbar.AddChildFrame(NSbar);

   LABcontra.InitializeTree();
   
   // Invisible Group
   InvisibleGroup INVcontra("INVcontra","NS NS Jigsaws");
   INVcontra.AddFrame(NS);
   INVcontra.AddFrame(NSbar);
   
   // Set NS NS~ mass equal to Z h mass
   SetMassInvJigsaw NSNSM("NSNSM", "M_{NSNS} = m_{Zh}");
   INVcontra.AddJigsaw(NSNSM);
   
   SetRapidityInvJigsaw NSNSR("NSNSR", "#eta_{NSNS} = #eta_{ZH}");
   INVcontra.AddJigsaw(NSNSR);
   NSNSR.AddVisibleFrames(LABcontra.GetListVisibleFrames());
   
   //MinMassesSqInvJigsaw MinMND("MinMND","min M_{D}, M_{ND}= M_{NDbar}",2);
   ContraBoostInvJigsaw MinMND("MinMND","min M_{ND}, M_{ND}= M_{NDbar}");
   INVcontra.AddJigsaw(MinMND);
   MinMND.AddVisibleFrame(Z, 0);
   MinMND.AddVisibleFrame(h, 1);
   MinMND.AddInvisibleFrame(NS, 0);
   MinMND.AddInvisibleFrame(NSbar, 1);

   LABcontra.InitializeAnalysis();
   //*/

   //Trigger Stuff
   string trgtit;
   string delim  = ",";
   string ourtrg;
   std::vector<string> thetrigs;
   int trgidx    = -1;
   int trgval;
   std::vector<int> trgvals;
   std::vector<int> trgidxs;
   TFile * fthen = 0;
   TFile * fnow = 0;
   //std::vector<string> trig18mu = {"HLT_Mu55_v","HLT_TkMu100_v"};
   std::vector<string> trig18mu = {"HLT_Mu50_v","HLT_TkMu100_v"};
   std::vector<string> trig17mu = {"HLT_Mu50_v","HLT_TkMu100_v"};
   std::vector<string> trig16mu = {"HLT_Mu50_v","HLT_TkMu50_v"};
   std::vector<string> trig18e = {"HLT_Ele32_WPTight_Gsf_v","HLT_Photon200_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v"};
   std::vector<string> trig17e = {"HLT_Ele35_WPTight_Gsf_v","HLT_Photon200_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v"};
   std::vector<string> trig16e = {"HLT_Ele27_WPTight_Gsf_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v","HLT_Photon175_v"};
   std::vector<string> trig18emu = {"HLT_Mu50_v","HLT_TkMu100_v","HLT_Ele32_WPTight_Gsf_v","HLT_Photon200_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v"};
   std::vector<string> trig17emu = {"HLT_Mu50_v","HLT_TkMu100_v","HLT_Ele35_WPTight_Gsf_v","HLT_Photon200_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v"};
   std::vector<string> trig16emu = {"HLT_Mu50_v","HLT_TkMu50_v","HLT_Ele27_WPTight_Gsf_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v","HLT_Photon175_v"};


   std::vector<std::vector<string>> trig18 = {{"no"},trig18emu,trig18e,{"no"},trig18mu};
   std::vector<std::vector<string>> trig17 = {{"no"},trig17emu,trig17e,{"no"},trig17mu};
   std::vector<std::vector<string>> trig16 = {{"no"},trig16emu,trig16e,{"no"},trig16mu};
   std::vector<std::vector<std::vector<string>>> trigs = {trig16,trig17,trig18};
   
   //counters
   float qcdwskimup = 0;
   float qcdwskimdwn = 0;
   float puweighttot = 0;
   float puweightuptot = 0;
   float puweightdwntot = 0;
   float puweightvtxuptot = 0;
   float puweightvtxdntot = 0;
   int emcounter = 0;
   int counttrigpass = 0;
   int countzpass    = 0;
   int counthpass    = 0;
   int countmetpass  = 0;
   int countpass     = 0;
   int zmumu = 0;
   int zmumuzee = 0;
   int zmumuzemu = 0;
   int zmumuzeezemu = 0;
   int zee = 0;
   int zeezemu = 0;
   int zemu = 0;
   int znorecfat = 0;
   int znounrecfat = 0;
   int znofat = 0;
   int countzcand = 0;
   int countfat = 0;
   int emumulead = 0;
   int emuelead  = 0;

   //std::cout<<"Number of files that make up the TChain: "<<fChain->GetListOfFiles()->GetSize()<<std::endl;this is weird, does not match files

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //Define some boos to signal events to write
      bool passZ    = false;
      bool passh    = false;
      bool passMET  = false;
      bool passTrig = false;
      bool passFil  = false;
      bool mumuchan = false;
      bool killelHEM = false;
      double channel = -1.0;

      //A counter, for my sanity
      if (jentry%25000 == 0) {
      	std::cout<<"    analyzing event "<<jentry<<std::endl;
      }

      
      //debug
      //std::cout<<"    analyzing event "<<jentry<<std::endl;
      //for (std::size_t i = 0; i < PDFweights->size();++i) {
      //	std::cout<<PDFweights->at(i)<<std::endl;
      //}
      //if (jentry%20 == 0) {
      //std::cout<<"    analyzing event "<<jentry<<std::endl;
      //}
      //if (jentry == 20) {
      //break;
      //}
     

      //Trigger decisions
      size_t pos = 0;
      string token;
      std::vector<int> checktrgidxs;
      thetrigs = trigs[year-16][anchan];
	for (std::size_t i = 0; i < thetrigs.size();++i){
	  ourtrg = thetrigs[i];
	  trgtit = fChain->GetBranch("TriggerPass")->GetTitle();
	  fthen = fChain->GetFile();
	  while ((pos = trgtit.find(delim)) != std::string::npos+1 && token != ourtrg) {
	    token = trgtit.substr(0,pos);
	    trgidx += 1;
	    trgtit.erase(0,pos+delim.length());
	  }
	  checktrgidxs.push_back(trgidx);
	  trgidx = -1;
	}
      if (checktrgidxs.size() != 0) {//This keeps the trg indexs in scope
	trgidxs = checktrgidxs;
      }
      trgvals.clear();
      for (std::size_t i = 0; i < trgidxs.size();++i){
	trgval = TriggerPass->at(trgidxs[i]);
	trgvals.push_back(trgval);
      }


      //Do the trigger for non-emu samples
      if ((trgvals[0] == 1 || trgvals[1] == 1) && anchan != 1) {//if not emu channel !!! Only works for mumu channerl at the moment
	passTrig = true;
	counttrigpass += 1;
      }

      //eeBadScFilter
      if (sampleType < 0) {
	if (fChain->GetLeaf("eeBadScFilter")->GetValue() == 1.0) {
	  passFil = true;
	}
      }

      //DY+Jets k-factors+GenParticleInfo
      float evntwkf_hold = 1.;
      TLorentzVector theGenZ;
      TLorentzVector theGenH;
      std::pair<float,float> qcdupdwn(1.0,1.0);
      float pdfwup = 1.0;
      float pdfwdn = 1.0;
      float puwnom = 1;
      float puwup =1;
      float puwdwn =1;
      float puwnumup =1;
      float puwnumdwn =1;
      if (sampleType > 0) {//Not Data
	int gpid;
	unsigned long ngen = GenParticles->size();
	for (unsigned long i = 0; i < ngen; ++i) {
	  int gpid = GenParticles_PdgId->at(i);
	  if (gpid == 23) {
	    theGenZ = GenParticles->at(i);
	  }
	  if (gpid == 25) {
	    theGenH = GenParticles->at(i);
	  }
	}

	//Calculate the k-factor
	float qcdnlosf   = 1;
	double qcdnnlosf = 1;
	float ewknlosf   = 1;
	float ntruenom = NVtx;
	float ntrueup =  NVtx;
	float ntruedwn = NVtx; 
	if (sampleType == 2) {//DY+Jets
	  qcdnlosf = 1.423*exp(-0.002257*theGenZ.Pt())+0.451;
	  if (qcdnlosf <= 0.0) {
	    qcdnlosf = 1.;
	  }
	  int zptbinqcd  = hqcdnnlosf->FindBin(theGenZ.Pt());
	  qcdnnlosf = hqcdnnlosf->GetBinContent(zptbinqcd);
	  if (qcdnnlosf <= 0.0) {
	    qcdnnlosf = 1.;
	  }
	  int zptbinewk = hewknlosf->FindBin(theGenZ.Pt());
	  ewknlosf = hewknlosf->GetBinContent(zptbinewk);
	  if (ewknlosf <= 0.0) {
	    ewknlosf = 1.;
	  }
	  evntwkf_hold = ewknlosf*qcdnnlosf*qcdnlosf;
	}
	//Do PDF uncertainty
	pdfwup = 1.0 + getPDFUncertainty(PDFweights);
	pdfwdn = 1.0 - getPDFUncertainty(PDFweights);
	
	//std::cout<<"PDF Up "<<pdfwup<<std::endl;
	//std::cout<<"PDF Dn "<<pdfwdn<<std::endl;

	//Do QCD Uncertainty
	qcdupdwn = getQCDScaleUpDwn(ScaleWeights);
	qcdwskimup+= qcdupdwn.second;
	qcdwskimdwn+= qcdupdwn.first;

	//Get Pileup Weight
	ntruenom  = hnvtxtrans->GetBinContent(hnvtxtrans->FindBin(NVtx));
	ntrueup   = ntruenom+hnvtxtrans->GetBinError(hnvtxtrans->FindBin(NVtx));
	ntruedwn  = ntruenom-hnvtxtrans->GetBinError(hnvtxtrans->FindBin(NVtx));
	puwnom    = hpuweight->GetBinContent(hpuweight->FindBin(ntruenom));
	puwup     = hpuweightup->GetBinContent(hpuweightup->FindBin(ntruenom));
	puwdwn    = hpuweightdwn->GetBinContent(hpuweightdwn->FindBin(ntruenom));
	puwnumup  = hpuweight->GetBinContent(hpuweight->FindBin(ntrueup));
	puwnumdwn = hpuweight->GetBinContent(hpuweight->FindBin(ntruedwn));
	puweighttot+=puwnom;
	puweightuptot+=puwup;
	puweightdwntot+=puwdwn;
	puweightvtxuptot+=puwnumup;
	puweightvtxdntot+=puwnumdwn;
	//std::cout<<"    NVTx       "<<NVtx<<std::endl;
	//std::cout<<"    True Num   "<<ntruenom<<std::endl;
	//std::cout<<"    True NumUp "<<ntrueup<<std::endl;
	//std::cout<<"    True NumDn "<<ntruedwn<<std::endl;
	//std::cout<<"       puw nom "<<puwnom<<std::endl;
	//std::cout<<"        puw up "<<puwup<<std::endl;
	//std::cout<<"        puw dn "<<puwdwn<<std::endl;			    
	//std::cout<<"    puw num up "<<puwnumup<<std::endl;
	//std::cout<<"    puw num dn "<<puwnumdwn<<std::endl;			    
	
      }

      
      //Only deal with exact events
      unsigned int nmutest = Muons->size();
      unsigned int neltest = Electrons->size();
      
      //Z exploration
      unsigned int nselmu = SelectedMuons->size();
      unsigned int nselel = SelectedElectrons->size();
      unsigned int nZmumu = ZCandidatesMuMu->size();
      unsigned int nZee = ZCandidatesEE->size();
      //unsigned int nZeu = 0;  
      unsigned int nZeu = ZCandidatesEU->size();//Does not work on old DY ntuples

      emcounter += nZeu;
      
      TLorentzVector leadmu;
      TLorentzVector subleadmu;
      double leadmupt_unc = 0;
      double subleadmupt_unc = 0;
      TLorentzVector leade;
      TLorentzVector subleade;
      TLorentzVector testmu;
      TLorentzVector teste;
      double lptmax = 0;
      unsigned int nZs = ZCandidates->size();
      TLorentzVector theZ;
      double baseZdiff = 99999;
      int muld = 0;
      int elld = 0;
      double evntwhem_hold = 1.0;

      //Build shifted Z
      TLorentzVector upleadmu;
      TLorentzVector dnleadmu;
      TLorentzVector upsubleadmu;
      TLorentzVector dnsubleadmu;
      TLorentzVector upZ;
      TLorentzVector dnZ;
      TLorentzVector corrZ;

      //Channel Flags
     ///*
      //std::cout<<"At the part for before the Z"<<std::endl;
      if (nZmumu > 0 && nZee == 0 && nZeu == 0 && anchan == 4){
	if (passTrig) countzcand += 1;
	//in binary 100, 4 in decimal
	channel = 4.;//4 in decimal
	mumuchan = true;
	//zmumu += 1;
	std::vector<TLorentzVector>::iterator muit;
	for (muit = SelectedMuons->begin(); muit != SelectedMuons->end();++muit) { //might need to move after Z
	  if (muit->Pt() > lptmax) {
	    lptmax = muit->Pt();
	    leadmu.SetPtEtaPhiM(muit->Pt(),muit->Eta(),muit->Phi(),muit->M());
	  }
	  else {
	    subleadmu.SetPtEtaPhiM(muit->Pt(),muit->Eta(),muit->Phi(),muit->M());;
	 }
	}

	//This omitts a check, but the leading muon should always be first
	//std::cout<<"About to check the uncertainty"<<std::endl;
	//leadmupt_unc = *SelectedMuonsTunepMomenErr->begin();
	//subleadmupt_unc = SelectedMuonsTunepMomenErr->back();
	//std::cout<<"The leading muon's error "<<leadmupt_unc<<std::endl;
	//std::cout<<"The subleading muon's error "<<subleadmupt_unc<<std::endl;

	//upleadmu.SetPtEtaPhiM(leadmu.Pt()+leadmupt_unc,leadmu.Eta(),leadmu.Phi(),leadmu.M());
	//dnleadmu.SetPtEtaPhiM(leadmu.Pt()-leadmupt_unc,leadmu.Eta(),leadmu.Phi(),leadmu.M());
	//upsubleadmu.SetPtEtaPhiM(subleadmu.Pt()+subleadmupt_unc,subleadmu.Eta(),subleadmu.Phi(),subleadmu.M());
	//dnsubleadmu.SetPtEtaPhiM(subleadmu.Pt()-subleadmupt_unc,subleadmu.Eta(),subleadmu.Phi(),subleadmu.M());
	//upZ = upleadmu+upsubleadmu;
	//dnZ = dnleadmu+dnsubleadmu;

	//Checks the muon momentum scale crap`
	//Prepare to do muon momentum tuning
	TLorentzVector leadmucorr;
	TLorentzVector leadmucorrup;
	TLorentzVector leadmucorrdn;
	TLorentzVector subleadmucorr;
	TLorentzVector subleadmucorrup;
	TLorentzVector subleadmucorrdown;


	if (sampleType > 0) {
	   std::string muera = "badstring";
	  if (year == 18){
	    muera = "2018";
	  }
	  if (year == 17){
	    muera = "2017";
	  }
	  if (year == 16){
	    muera = "2016";
	  }
	GEScaleSyst GE(muera);
	GE.SetVerbose(0);

	int qlmu = 0;
	int qslmu = 0;
	//std::cout<<Muons->size()<<std::endl;
	for (muit = Muons->begin(); muit != Muons->end();++muit){
	  if (deltaR(leadmu,*muit) < 0.2 && deltaR(subleadmu,*muit) > 0.2) {
	    qlmu = Muons_charge->at(muit-Muons->begin());
	  }
	  else if (deltaR(subleadmu,*muit) < 0.2 && deltaR(leadmu,*muit) > 0.2) {
	    qslmu = Muons_charge->at(muit-Muons->begin());
	  }
	  else if (deltaR(subleadmu,*muit) < 0.2 && deltaR(leadmu,*muit) < 0.2) {
	    if (deltaR(subleadmu,*muit) < deltaR(leadmu,*muit)){
	      qslmu = Muons_charge->at(muit-Muons->begin());
	    }
	    else if (deltaR(subleadmu,*muit) > deltaR(leadmu,*muit)){
	      qlmu = Muons_charge->at(muit-Muons->begin());
	    }
	  }
	  //else {
	  //  std::cout<<"!!!!!!!!!!!!!!!!!No muons matching the Z muons found"<<std::endl;
	  //  std::cout<<" size of muon array: "<<Muons->size()<<std::endl;
      
	    
	  //}
	}
	//std::cout<<"Leading muon charge: "<<qlmu <<std::endl;
	//std::cout<<"subLeading muon charge: "<<qslmu <<std::endl;
	leadmucorr = GE.GEScaleCorrLvec(leadmu.Pt(),leadmu.Eta(),leadmu.Phi(),qlmu,0,0);
	leadmucorrup = GE.GEScaleCorrLvec(leadmu.Pt(),leadmu.Eta(),leadmu.Phi(),qlmu,0,1);
	leadmucorrdn = GE.GEScaleCorrLvec(leadmu.Pt(),leadmu.Eta(),leadmu.Phi(),qlmu,0,2);
	subleadmucorr = GE.GEScaleCorrLvec(subleadmu.Pt(),subleadmu.Eta(),subleadmu.Phi(),qslmu,0,0);
	subleadmucorrup = GE.GEScaleCorrLvec(subleadmu.Pt(),subleadmu.Eta(),subleadmu.Phi(),qslmu,0,1);
	subleadmucorrdown = GE.GEScaleCorrLvec(subleadmu.Pt(),subleadmu.Eta(),subleadmu.Phi(),qslmu,0,2);
	//std::cout<<"Original leading muon pT: "<<leadmu.Pt()<<std::endl;
	//std::cout<<"corr leading muon pT:     "<<leadmucorr.Pt()<<std::endl;

	corrZ = leadmucorr+subleadmucorr;
	upZ   = leadmucorrup+subleadmucorrup;
	dnZ   = leadmucorrdn+subleadmucorrdown;
	}

	
	std::vector<TLorentzVector>::iterator zit;
	for (zit = ZCandidatesMuMu->begin(); zit != ZCandidatesMuMu->end(); ++zit) {
	  double massZdiff = std::abs(91.1876 - zit->M());
	  if ((massZdiff < baseZdiff) && (zit->M() >= zmwinlow) && (zit->M() <= zmwinhi)) {
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Eta(),zit->Phi(),zit->M());
	    passZ = true;
	    //if (passTrig){
	      zmumu += 1;
	      //}

	  }
	}
      }
      if (nZmumu > 0 && nZee > 0 && nZeu == 0  && anchan == 6) {
	//110 in binary, 6 in decimal
	channel = 6.;//
	//if (passTrig){
	zmumuzee += 1;
	//}
      }
      if (nZmumu > 0 && nZee == 0 && nZeu > 0 && anchan == 5) {
	//101 in binary, 5 in decimal
	channel = 5.;//
	int matchmu = 0;
	//if (passTrig) {
	  std::vector<TLorentzVector>::iterator muit;
	  std::vector<TLorentzVector>::iterator muit2;
	  for (muit = SelectedMuons->begin(); muit != SelectedMuons->end();++muit) { //might need to move after Z
	    for (muit2 = muit+1;muit2 != SelectedMuons->end();++muit2){
	      if (*muit2 == *muit) {
		matchmu += 1;}
	    }
	    //}
	  if (matchmu == 0){
	    zmumuzemu +=1;
	  }
	}
      }
      if (nZmumu > 0 && nZee > 0 && nZeu > 0  && anchan == 7) {
	//111 in binary, 7 in decimal
	channel = 7.;
	//if (passTrig) {
	zmumuzeezemu +=1;
	//}
      }
      if (nZmumu == 0 && nZee > 0 && nZeu == 0  && anchan == 2) {
	//010
	channel = 2.;
	//if (passTrig) {
	zee +=1;
	//}
       ///*
	std::vector<TLorentzVector>::iterator eit;
	for (eit = SelectedElectrons->begin(); eit != SelectedElectrons->end();++eit) { 
	  if (eit->Pt() > lptmax) {
	    lptmax = eit->Pt();
	    leade.SetPtEtaPhiM(eit->Pt(),eit->Eta(),eit->Phi(),eit->M());
	  }
	  else {
	    subleade.SetPtEtaPhiM(eit->Pt(),eit->Eta(),eit->Phi(),eit->M());;
	 }
	}
	std::vector<TLorentzVector>::iterator zit;
	for (zit = ZCandidatesEE->begin(); zit != ZCandidatesEE->end(); ++zit) {
	  double massZdiff = std::abs(91.1876 - zit->M());
	  if ((massZdiff < baseZdiff) && (zit->M() > zmwinlow) && (zit->M() < zmwinhi)) {
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Eta(),zit->Phi(),zit->M());
	    passZ = true;
	  }
	}
	}
      if (nZmumu == 0 && nZee > 0 && nZeu > 0 && anchan == 3) {
	//011
	channel = 3.;
	//if (passTrig) {
	zeezemu +=1;
	//}
      }
      //std::cout<<"At the part for before the Z"<<std::endl;
      if (nZmumu == 0 && nZee == 0 && nZeu > 0  && anchan == 1) {
	//001
	//std::cout<<"At least we are in the emu channel part!"<<std::endl;
	channel = 1.;
	std::vector<TLorentzVector>::iterator zit;
	for (zit = ZCandidatesEU->begin(); zit != ZCandidatesEU->end(); ++zit) {
	  double massZdiff = std::abs(91.1876 - zit->M());
	  if ((massZdiff < baseZdiff) && (zit->M() > zmwinlow) && (zit->M() < zmwinhi)) {
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Eta(),zit->Phi(),zit->M());
	    passZ = true;
	    zemu += 1;
	    testmu = SelectedMuons->at(0);
	    teste = SelectedElectrons->at(0);
	    //HEM Test
	    if (teste.Eta() > -3.0 && teste.Eta() < -1.3 && teste.Phi() > -1.57 && teste.Phi() < -0.87 && ((year >= 180) or year == 18)){
	      killelHEM = true;
	      evntwhem_hold = 0.36;
	    }
	    if (testmu.Pt() > teste.Pt()) {
	      leadmu = testmu;
	      subleade = teste;
	      //std::cout<<"    The muon is leading, the electron is subleading!"<<std::endl;
	      muld = 1;
	    }
	    else {
	      leade = teste;
	      subleadmu = testmu;
	      //std::cout<<"    The electron is leading, the muon is subleading!"<<std::endl;
	      elld = 1;
	    }
	  }
	}
      }
     //*/
      //Do the emu trigger
      //emutrig
      //passTrig = true;
      bool passmu = false;
      bool passel = false;
      if (Electrons->size() > 0){
	eventleade = Electrons->at(0);
      }
      if (Muons->size() > 0){
	eventleadmu = Muons->at(0);
      }
      if (sampleType > 0 && anchan == 1 ) {//if emu channel
	if (Electrons->size() > 0 && Muons->size() > 0) {
	  if (eventleade.Pt() > 27 ) {
	    if (trgvals[2] == 1 || trgvals[3] == 1 || trgvals[4] == 1) {
	      passTrig = true;
	      counttrigpass += true;
	      passel = true;
	    }
	  }
	  if (eventleadmu.Pt() > 50 ) {
	    if (trgvals[0] == 1 || trgvals[1] == 1) {
	      passTrig = true;
	      counttrigpass += true;
	      passmu = true;
	    }
	  }
	}
      }
      //if (passTrig){
      //std::cout<<"Passed the trigger"<<std::endl;
      //}
      //std::cout<<"Magnitude of trig vals "<<getMagnitude(trgvals)<<std::endl;
      if (sampleType < 0 && anchan == 1 ) {//if emu channel data
	if (Electrons->size() > 0 && Muons->size() > 0) {
	  //eventleade = Electrons->at(0);
	  //eventleadmu = Muons->at(0);
	  if (getMagnitude(trgvals) > 0.0){
	    passTrig = true;
	    counttrigpass += 1;
	  }
	}

      }
      //Z Candidate Build
      //For old ntuples
      /*
	leadmupt_unc = 0;
	subleadmupt_unc = 0;
      if (nselmu > 0 && nselel == 0 && anchan == 4) {
      	mumuchan = true;
	channel = 4;
	std::vector<TLorentzVector>::iterator muit;
	for (muit = SelectedMuons->begin(); muit != SelectedMuons->end();++muit) { //might need to move after Z
	  if (muit->Pt() > lptmax) {
	    lptmax = muit->Pt();
	    leadmu.SetPtEtaPhiM(muit->Pt(),muit->Eta(),muit->Phi(),muit->M());
	  }
	  else {
	    subleadmu.SetPtEtaPhiM(muit->Pt(),muit->Eta(),muit->Phi(),muit->M());;
	 }
	}
      }
      if (nZs > 0) {
	std::vector<TLorentzVector>::iterator zit;
	for (zit = ZCandidates->begin(); zit != ZCandidates->end(); ++zit) {
	  double massZdiff = std::abs(91.1876 - zit->M());
	  if ((massZdiff < baseZdiff) && (zit->M() > zmwinlow) && (zit->M() < zmwinhi)) {
	    baseZdiff = massZdiff;
	    theZ.SetPtEtaPhiM(zit->Pt(),zit->Eta(),zit->Phi(),zit->M());
	    passZ = true;
	  }
	}
      }
      //*/


      //Do Lepton SF
      std::vector<double> muidsfvec = {1,0,0};//total muon weight
      std::vector<double> elidsfvec = {1,0,0};//total electron ID weight
      std::vector<double> elrecosfvec = {1,0,0};//total electron Reco weight
      std::vector<double> muidsfleadv = {1,0,0};//muon channel leadmuon
      std::vector<double> elidsfleadv = {1,0,0};//electron channel leadelectron
      std::vector<double> elrecoslleadv = {1,0,0};//electron channel leadelectron
      std::vector<double> muidsfsublv = {1,0,0};//muon channel subleadmuon
      std::vector<double> elidsfsublv = {1,0,0};//electron channel subleadelectron
      std::vector<double> elrecosublv = {1,0,0};//electron channel subleadelectron
      std::vector<double> mutrigsfvec= {1,0,0};//muon trigger sf
      std::vector<double> eltrigsfvec= {1,0,0};//electron trigger
      
      if (sampleType > 0 && anchan == 4) {//not data, mumuchannel
	//std::cout<<"Doing Muon ID sf"<<std::endl;
	muidsfleadv = GetMuonPtEtaSF(year,hmuonsf,leadmu,120.0,true);
	muidsfsublv = GetMuonPtEtaSF(year,hmuonsf,subleadmu,120.0,true);
	muidsfvec = combineTheLeptonSF(muidsfleadv,muidsfsublv);
	//std::cout<<"Doing Muon Trigger sf"<<std::endl;
	//mutrigsfvec = GetMuonPtEtaSF(year,hmuontrigsf,eventleadmu,200.0,false);//old way
	std::vector<double> lmutrigsfv = GetMuonPtEtaSF(year,hmuontrigsf,leadmu,500.0,false);
	std::vector<double> slmutrigsfv = GetMuonPtEtaSF(year,hmuontrigsf,subleadmu,500.0,false);
	mutrigsfvec = GetMuonTriggerSF(lmutrigsfv,slmutrigsfv,leadmu,subleadmu,500.0,hmuontrigef);
      }
      else if (sampleType > 0 && anchan == 2) {//not data, ee channel
	elidsfleadv = GetElectronPtEtaSF(year,helectronIDsf,leade,false);
	elidsfsublv = GetElectronPtEtaSF(year,helectronIDsf,subleade,false);
	elrecoslleadv = GetElectronPtEtaSF(year,helectronsf,leade,false);
	elrecosublv = GetElectronPtEtaSF(year,helectronIDsf,subleade,false);
	elidsfvec = combineTheLeptonSF(elidsfleadv,elidsfsublv);
	elrecosfvec = combineTheLeptonSF(elrecoslleadv,elrecosublv);
      }
      else if (sampleType > 0 && anchan == 1) {//not data, emu channel
	//Trigger scale factors should be applied based on the trigger object
	//Done higher in the code, for maximal confusion.
	//ID scale factors can be based off of the objects, done here
	if (muld == 1) {
	  muidsfvec = GetMuonPtEtaSF(year,hmuonsf,leadmu,120.0,true);
	  elrecosfvec = GetElectronPtEtaSF(year,helectronsf,subleade,false);
	  elidsfvec   = GetElectronPtEtaSF(year,helectronIDsf,subleade,false);
	}
	if (elld == 1) {
	  elrecosfvec = GetElectronPtEtaSF(year,helectronsf,leade,false);
	  elidsfvec   = GetElectronPtEtaSF(year,helectronIDsf,leade,false);
	  muidsfvec   = GetMuonPtEtaSF(year,hmuonsf,subleadmu,120.0,true);
	}
	if (passmu && not passel) {
	  mutrigsfvec = GetMuonPtEtaSF(year,hmuontrigsf,eventleadmu,500.0,false);
	}
	else if (passel && not passmu){
	  eltrigsfvec = GetElectronPtEtaSF(year,helectronTRIGsf,eventleade,true);
	}
	else {//using this background for the muon channel
	  mutrigsfvec = GetMuonPtEtaSF(year,hmuontrigsf,eventleadmu,500.0,false);
	}
      }
 
      //Higgs Candidate Build
      //JetBranch = JetsAK8Clean;
      // unsigned long nfat = JetsAK8->size();
      //unsigned long nfat = JetBranch->size();
      unsigned long nfat = JetsAK8Clean->size();
      TLorentzVector theh;
      theh.SetPtEtaPhiM(0.0,0.0,0.0,0.0);;
      TLorentzVector fat;
      double basehdiff = 99999;
      double basept = 0;
      double hsd = 0;
      double fsd = 0;
      double hsf = 1.0;
      double hbtagunup = 0;
      double hbtagundwn = 0;
      int    sfidx = 0;
      double hdmdhbbvqcd = 0;
      double hdmdzbbvqcd = 0;
      double hdmdzhbbvqcd = 0;
      double hmiddb = 0;
      bool   fid = 0;
      std::pair<double,double> metshift;
      bool hemshift = 0;
      //double metshift = 0;
      double jerfactor = 0;

      /*
      if (nZs > 0 && nfat == 0){
	//std::cout<<"Found an event with a Z Candidate but no reclustered fat jets"<<std::endl;
	znorecfat +=1;
      }

      if (nZs > 0 && nunfat == 0){
	//std::cout<<"Found an event with a Z Candidate but no reclustered fat jets"<<std::endl;
	znounrecfat +=1;
      }

      if (nZs > 0 && nfat == 0 && nunfat == 0){
	//std::cout<<"Found an event with a Z Candidate but no reclustered fat jets"<<std::endl;
	znofat +=1;
      }
      */

      ///*
      //reclustered jets
      if (nfat > 0) {
	if (passTrig) countfat += 1;
	for (unsigned long i =0; i < nfat; ++i) {
	  fat = JetsAK8Clean->at(i);
	  fsd = JetsAK8Clean_softDropMass->at(i);
	  fid = JetsAK8Clean_ID->at(i);
	  //std::cout<<"right before jer factor"<<std::endl;
	  jerfactor = JetsAK8Clean_jerFactor->at(i);
	  hjetpt->Fill(fat.Pt(),evntwkf_hold);
	  hjeteta->Fill(fat.Eta(),evntwkf_hold);
	  hjetphi->Fill(fat.Phi(),evntwkf_hold);
	  jec_unc->setJetEta(fat.Eta());
	  jec_unc->setJetPt(fat.Pt());
	  double unc = 0;
	  //std::cout<<"right before jec factor"<<std::endl;
	  unc = std::abs(jec_unc->getUncertainty(true));
	  double jecsysfac = 1 + jecsys*unc;
	  fat = fat*jecsysfac;
	  /*
	  TLorentzVector hemfat = doHEMshiftJet(fat,fid);
	  if (hemfat.M() != fat.M()){
	    //std::cout<<"Found a shifted jet, need a MET shift"<<std::endl;
	    //std::cout<<"                   Nominal METx Shift "<<metshift.first<<std::endl;
	    //std::cout<<"                   Nominal METy Shift "<<metshift.second<<std::endl;
	    hemshift = 1;
	    std::pair<double,double> metshiftj = getMETCorrHEM(fat,hemfat);
	    double metshiftxtot = metshift.first+metshiftj.first;
	    double metshiftytot = metshift.second+metshiftj.second;
	    metshift = std::make_pair(metshiftxtot,metshiftytot);
	    //std::cout<<"                       new METx Shift "<<metshift.first<<std::endl;
	    //std::cout<<"                       new METy Shift "<<metshift.second<<std::endl;
	    fat = hemfat;
	  }
	  //*/

	  if (jeRsys != 0){//Should do this only if we want to look at JEr systematics
	    TLorentzVector closest_genjet = closestgenjetParticle(fat, GenJetsAK8);
	    fat = fat * (1./jerfactor); // undo smearing 
	    // redo smearing with new smearing factor from txt file
	    float recopt = fat.Pt();
	    float recoeta = fat.Eta();
	    float abseta = fabs(recoeta);
	    float rho = fat.Rho();
	    //
	    float resolution = resolution_.getResolution({{JME::Binning::JetPt, recopt}, {JME::Binning::JetEta, recoeta}, {JME::Binning::Rho, rho}});
	    //
	    //
	    float genpt = -1;
	    if(!(closest_genjet.Pt() == 0) && deltaR(closest_genjet, fat) < 0.5*0.8){
	      genpt = closest_genjet.Pt();
	    }
	    if( fabs(genpt-recopt) > 3*resolution*recopt){
	      genpt=-1;
	    }
	    if(genpt < 15.0f) {
	      genpt=-1.;
	    }
	  //
	  //
	  float c = -1;
	  if (jeRsys == 1) { 
	    c = res_sf_.getScaleFactor({{JME::Binning::JetPt, recopt}, {JME::Binning::JetEta, recoeta}}, Variation::UP);
	  } else if (jeRsys == -1) {
	    c = res_sf_.getScaleFactor({{JME::Binning::JetPt, recopt}, {JME::Binning::JetEta, recoeta}}, Variation::DOWN);
	  } else{
	    c = res_sf_.getScaleFactor({{JME::Binning::JetPt, recopt}, {JME::Binning::JetEta, recoeta}});
	  }
	  //
	  // Calculate the new pt
	  float new_pt = -1.;

	  //An alternative approach, which does not require the presence of a matching particle-level jet, is the stochastic smearing: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution

	  //TRandom rand((int)(1000*abseta));
	  //float random_gauss = rand.Gaus(0, resolution);
	  if(genpt > 0){
	    new_pt = std::max(0.0f, genpt + c * (recopt - genpt));
	  } else{
	    float random_gauss = gRandom->Gaus(0, resolution);

	    new_pt = recopt * (1 + random_gauss*sqrt(std::max(c*c-1, 0.0f)));
	  }

	  fat *= new_pt / recopt;
	  //
	  //
	  hjetptJER->Fill(fat.Pt(),evntwkf_hold);
	  hjetetaJER->Fill(fat.Eta(),evntwkf_hold);
	  hjetphiJER->Fill(fat.Phi(),evntwkf_hold);
	  }

	  double masshdiff = std::abs(125.18 - fsd);
	  if ((masshdiff < basehdiff) && (fat.Pt() > hptcut) && fid && std::abs(fat.Eta()) < 2.4 && (fsd > 10)) {
	    basehdiff = masshdiff;
	    theh = fat;
	    hsd = fsd;
	    hdmdhbbvqcd  = JetsAK8Clean_DeepMassDecorrelTagHbbvsQCD->at(i);
	    hdmdzbbvqcd  = JetsAK8Clean_DeepMassDecorrelTagZbbvsQCD->at(i);
	    hdmdzhbbvqcd = JetsAK8Clean_DeepMassDecorrelTagZHbbvsQCD->at(i);
	    hmiddb = JetsAK8Clean_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(i);
	    passh = true;
	  }
	}
      }

      if (sampleType > 0) {
	sfidx = hbtagsf->FindBin(theh.Pt());
	if (sfidx > 0) {
	  hsf = hbtagsf->GetBinContent(sfidx);
	  hbtagunup = hbtagsfuncup->GetBinContent(sfidx);
	  hbtagundwn = hbtagsfuncdwn->GetBinContent(sfidx);
	}
      }
      //*/
      
      //btag sf debug
      //if (passh) {
      //std::cout<<"The jet pT "<<theh.Pt()<<std::endl;
      //std::cout<<"The sf "<<hsf<<std::endl;	
      //std::cout<<"The sfidx "<<sfidx<<std::endl;
      //std::cout<<"The sferrup  "<<hbtagunup<<std::endl;
      //std::cout<<"The sferrdwn "<<hbtagundwn<<std::endl;
      //}
	    
      //unreclustered jets
      /*
      if (nfat > 0) {
      for (unsigned long i =0; i < nfat; ++i) {
	fat = JetsAK8->at(i);
        fsd = JetsAK8_softDropMass->at(i);
        fid = JetsAK8_ID->at(i);
        double masshdiff = std::abs(125.18 - fsd);
        if ((masshdiff < basehdiff) && (fat.Pt() > hptcut) && fid && std::abs(fat.Eta()) < 2.4 && (fsd > 10)) {
	  basehdiff = masshdiff;
          theh = fat;
          hsd = fsd;
          hdmdhbbvqcd  = JetsAK8_DeepMassDecorrelTagHbbvsQCD->at(i);
          hdmdzbbvqcd  = JetsAK8_DeepMassDecorrelTagZbbvsQCD->at(i);
          hdmdzhbbvqcd = JetsAK8_DeepMassDecorrelTagZHbbvsQCD->at(i);
          hmiddb = JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(i);
          passh = true;
        }
      }
      }
      //*/

      //MET
      //double ptmiss     = fChain->GetLeaf("METclean")->GetValue(0);
      //double ptmiss_phi = fChain->GetLeaf("METPhiclean")->GetValue();
      double ptmiss     = fChain->GetLeaf("MET")->GetValue(0);
      double ptmiss_phi = fChain->GetLeaf("METPhi")->GetValue(0);

      //std::cout<<"Original MET "<<ptmiss<<std::endl;

      if (systidx > -1 && sampleType > 0){//if  a met systematic call has been made (otherwise the idx is initiallize to -999)
	if (metsys[systidx] > 0) {//Checks up or down
	  ptmiss = METUp->at(systidx);//The index of the corresponding met in the ntuple is the same
	  ptmiss_phi = METPhiUp->at(systidx);
	}
	else if (metsys[systidx] < 0) {
	  ptmiss = METDown->at(systidx);
	  ptmiss_phi = METPhiDown->at(systidx);
	}
	else {
	  continue;
	}
      }

      //std::cout<<"MET after systematics "<<ptmiss<<std::endl;
      
      //HEM15/16 test
      if (hemshift) {
	double ptmiss_pxbeta  = ptmiss*std::cos(ptmiss_phi);
	double ptmiss_pybeta  = ptmiss*std::sin(ptmiss_phi);
	double smetx = ptmiss_pxbeta+metshift.first;
	double smety = ptmiss_pybeta+metshift.second;
	double retmet = std::sqrt(smetx*smetx+smety*smety);
	TVector2 metv(ptmiss_pxbeta,ptmiss_pybeta);
	TVector2 metshiftv(metshift.first,metshift.second);
	TVector2 shiftedmetv = metv+metshiftv;
	double metshiftphi = shiftedmetv.Phi();
	if (metshiftphi > M_PI){
	  metshiftphi = metshiftphi - 2*M_PI;
	}

	ptmiss = shiftedmetv.Mod();
	ptmiss_phi = metshiftphi;

	
	/*
	std::cout<<"                   Original MET:     "<<ptmiss<<std::endl;
	std::cout<<"                   Ori_phi  MET:     "<<ptmiss_phi<<std::endl;
	std::cout<<"      Check of MET nom TVector2, phi "<<metv.Phi()<<std::endl;
	std::cout<<" Check phi conversion, 2pi + ori phi "<<2*M_PI+ptmiss_phi<<std::endl;
	std::cout<<" Check phi conversion, back to   val "<<metv.Phi()-2*M_PI<<std::endl;
	
	std::cout<<"resummed shifted met:      "<<retmet<<std::endl;
	std::cout<<"Check of magnitude of shifted met "<<shiftedmetv.Mod()<<std::endl;
	
	std::cout<<"Phi of the shifted MET   , phi "<<shiftedmetv.Phi()<<std::endl;
	std::cout<<"Phi of the shifted MET   , phi "<<metshiftphi<<std::endl;
	//*/
      }

      //muon momentum shift
      /*
      double ptmiss_x  = ptmiss*std::cos(ptmiss_phi);
      double ptmiss_y  = ptmiss*std::sin(ptmiss_phi);
      TVector2 metvec0(ptmiss_x,ptmiss_y);
      double leadmu_x  = leadmu.Pt()*std::cos(leadmu.Phi());
      double leadmu_y  = leadmu.Pt()*std::sin(leadmu.Phi());
      TVector2 leadmuptvec(leadmu_x,leadmu_y);
      double subleadmu_x  = subleadmu.Pt()*std::cos(subleadmu.Phi());
      double subleadmu_y  = subleadmu.Pt()*std::sin(subleadmu.Phi());
      TVector2 subleadmuptvec(subleadmu_x,subleadmu_y);
      TVector2 leadmuupvec = getShiftedPtvec(leadmu.Pt(),leadmupt_unc,leadmu.Phi(),1);
      TVector2 leadmudnvec = getShiftedPtvec(leadmu.Pt(),leadmupt_unc,leadmu.Phi(),-1);
      TVector2 subleadmuupvec = getShiftedPtvec(subleadmu.Pt(),subleadmupt_unc,subleadmu.Phi(),1);
      TVector2 subleadmudnvec = getShiftedPtvec(subleadmu.Pt(),subleadmupt_unc,subleadmu.Phi(),-1);

      TVector2 cleanedmet = metvec0+leadmuptvec+subleadmuptvec;
      TVector2 metmuonsup = cleanedmet-leadmuupvec-subleadmuupvec;
      TVector2 metmuonsdn = cleanedmet-leadmudnvec-subleadmudnvec;

      double muupmet = metmuonsup.Mod();
      double mudnmet = metmuonsdn.Mod();
      */
      
      //std::cout<<"Initial MET                             : "<<ptmiss<<std::endl;
      //std::cout<<"    MET with the muons no longer counted: "<<cleanedmet.Mod()<<std::endl;
      //std::cout<<"                       MET with up muons: "<<metmuonsup.Mod()<<std::endl;
      //std::cout<<"                       MET with dn muons: "<<metmuonsdn.Mod()<<std::endl;



      //Do MET Correction
      bool isMC = true;
      if (sampleType < 0) {
	isMC = false;
      }

      
      std::pair<double,double> metxycorrpair = METXYCorr_Met_MetPhi(ptmiss,ptmiss_phi,fChain->GetLeaf("RunNum")->GetValue(0),yearstr,isMC,fChain->GetLeaf("NVtx")->GetValue(0));

      //std::cout<<"MET after xy shift "<<ptmiss<<std::endl;
      
      /*
      std::cout<<"Original MET:     "<<ptmiss<<std::endl;
      std::cout<<"correct  MET:     "<<metxycorrpair.first<<std::endl;
      std::cout<<"Ori_phi  MET:     "<<ptmiss_phi<<std::endl;
      std::cout<<"correct  MET_phi: "<<metxycorrpair.second<<std::endl;
      //*/
      double ptmiss_px  = ptmiss*std::cos(ptmiss_phi);
      double ptmiss_py  = ptmiss*std::sin(ptmiss_phi);
      TVector3 met3     = TVector3(ptmiss_px,ptmiss_py,0.0);
      TVector3 met3xy   = TVector3(metxycorrpair.first*std::cos(metxycorrpair.second),metxycorrpair.first*std::sin(metxycorrpair.second),0.0);

      //met jes syst debug
      //std::cout<<"This is the jec syst multiplier: "<<jecsys<<std::endl;
      //std::cout<<"This is the value of the nominal MET:       "<<fChain->GetLeaf("MET")->GetValue(0)<<std::endl;
      //std::cout<<"This is the value of the METUp:             "<<METUp->at(1)<<std::endl;
      //std::cout<<"This is the value of the METDown:           "<<METDown->at(1)<<std::endl;
      //std::cout<<"This is the value of the MET you are using: "<<ptmiss<<std::endl;
      
      //recursive jigsaw
      //  /*
      LABcontra.ClearEvent();
      //INVcontra.SetLabFrameThreeVector(met3);
      INVcontra.SetLabFrameThreeVector(met3xy);
      Z.SetLabFrameFourVector(theZ);
      h.SetLabFrameFourVector(theh);
      LABcontra.AnalyzeEvent();
      mEstZp = Zp.GetMass();
      mEstND = ND.GetMass();
      mEstNS = NS.GetMass();
      //*/

      //Just for the efficiency plots
      if (passZ && (channel == anchan) && theZ.Pt() > 100.0) {
	hzbuild_pt->Fill(theZ.Pt());
	hzbuild_eta->Fill(theZ.Eta());
	hlmubuild_pt->Fill(leade.Pt());
	hlmubuild_eta->Fill(leade.Eta());
      }

      if (passZ && passTrig && (channel == anchan) && theZ.Pt() > 100.0) {
	countzpass +=1 ;
	hzpasstrig_pt->Fill(theZ.Pt());
	hzpasstrig_eta->Fill(theZ.Eta());
	hlmupasstrig_pt->Fill(leade.Pt());
	hlmupasstrig_eta->Fill(leade.Eta());
      }

      //if (passh && passZ && passTrig && (channel == anchan)) {//not mucmuchan, but if channel == anchan
        if (passZ && passTrig && (channel == anchan)) {//not mucmuchan, but if channel == anchan 
	//if (passh && passZ ) {//Removed Z Channel Requirement
	//if (channel == anchan && passZ) {//id'd lepton and gen higgs plots
	hCandidate = theh;
	hCandidate_pt  = theh.Pt();
	hCandidate_phi = theh.Phi();
	hCandidate_eta = theh.Eta();
	hCandidate_m   = theh.M();
	hCandidate_sd  = hsd;
	hCandidate_dmdhbbvqcd = hdmdhbbvqcd;
	hCandidate_dmdzbbvqcd = hdmdzbbvqcd;
	hCandidate_dmdzhbbvqcd = hdmdzhbbvqcd;
	hCandidate_middb = hmiddb;
	ZCandidate = theZ;
	ZCandidate_pt  = theZ.Pt();
	ZCorrCandidate_pt  = corrZ.Pt();
	ZupCandidate_pt  = upZ.Pt();
	ZdnCandidate_pt  = dnZ.Pt();
	ZCandidate_phi = theZ.Phi();
	ZCandidate_eta = theZ.Eta();
	ZCorrCandidate_m   = corrZ.M();
	ZCandidate_m   = theZ.M();
	ZupCandidate_m   = upZ.M();
	ZdnCandidate_m   = dnZ.M();
	//event weight ones
	evntwhem = evntwhem_hold;
	evntwkf = evntwkf_hold;
	evntwbtag = hsf;
	evntwmusf = muidsfvec[0];
	evntwmutrig = mutrigsfvec[0];
	evntweltrig = eltrigsfvec[0];
	evntwelidsf = elidsfvec[0];
	evntwelrecosf = elrecosfvec[0];
	//The scale factor uncertainties
	btaguncup  = hbtagunup;
	btaguncdwn = hbtagundwn;
	muiduncup  = muidsfvec[1];
	muiduncdwn = muidsfvec[2];
	mutriguncup  = mutrigsfvec[1];
	mutriguncdwn = mutrigsfvec[2];
	eltriguncup  = eltrigsfvec[1];
	eltriguncdwn = eltrigsfvec[2];
	eliduncup  = elidsfvec[1];
	eliduncdwn = elidsfvec[2];
	elrecouncup = elrecosfvec[1];
	elrecouncdwn = elrecosfvec[2];
	//Theory Uncs
	pdfweightup = pdfwup;
	pdfweightdn = pdfwdn;
	qcdweightup = qcdupdwn.second;
	qcdweightdn = qcdupdwn.first;
	//pileup
	puweightnom = puwnom;
	puweightup = puwup;
	puweightdn = puwdwn;
	puweightvtxup = puwnumup;
	puweightvtxdn = puwnumdwn;

	LMuCandidate = leadmu;
	LMuCandidate_pt  = leadmu.Pt();
	LMuCandidate_ptunc  = leadmupt_unc;
	LMuCandidate_phi = leadmu.Phi();
	LMuCandidate_eta = leadmu.Eta();
	LMuCandidate_m   = leadmu.M();
	sLMuCandidate = subleadmu;
	sLMuCandidate_pt  = subleadmu.Pt();
	sLMuCandidate_ptunc  = subleadmupt_unc;
	sLMuCandidate_phi = subleadmu.Phi();
	sLMuCandidate_eta = subleadmu.Eta();
	sLMuCandidate_m   = subleadmu.M();
	LEleCandidate = leade;
	LEleCandidate_pt  = leade.Pt();
	LEleCandidate_phi = leade.Phi();
	LEleCandidate_eta = leade.Eta();
	LEleCandidate_m   = leade.M();
	sLEleCandidate = subleade;
	sLEleCandidate_pt  = subleade.Pt();
	sLEleCandidate_phi = subleade.Phi();
	sLEleCandidate_eta = subleade.Eta();
	sLEleCandidate_m   = subleade.M();
	ghCandidate_pt  = theGenH.Pt();
	ghCandidate_phi = theGenH.Phi();
	ghCandidate_eta = theGenH.Eta();
	ghCandidate_m   = theGenH.M();
	gzCandidate_pt  = theGenZ.Pt();
	gzCandidate_phi = theGenZ.Phi();
	gzCandidate_eta = theGenZ.Eta();
	gzCandidate_m   = theGenZ.M();
	metusable = ptmiss;
	metphiusable = ptmiss_phi;
	metxyphicorr = metxycorrpair.second;
	metxycorr= metxycorrpair.first;
	//metmuup = muupmet;
	//metmudn = mudnmet;
	channelflag = channel;
	counthpass += 1;
	if (muld > 0) emumulead+=1;
	if (elld > 0) emuelead+=1; 
      }

	//std::cout<<"Pass Z "<<passZ<<std::endl;
	//std::cout<<"Pass trig "<<passTrig<<std::endl;
      
      //Fill the Tree
      if (Cut(ientry) < 0) continue;
      //if (passZ && passh && passTrig && sampleType > 0 && (channel == anchan)) {//usual

      if (passZ && passh && passTrig && sampleType > 0 && (channel == anchan)) {
	//if (passZ && passh && sampleType > 0) {//for Zee channel checks
	//if (passZ && (sampleType > 0) && (channel == anchan)){
	//std::cout<<"This is where I think I am, in this passing place"<<std::endl;
	trimTree->Fill();
	countpass += 1;
      }

	/////ucomment!!!
	//else if (passZ && passh && passTrig && sampleType < 0 && passFil && (channel == anchan)) {
      else if (passZ && passh && passTrig && sampleType < 0 && passFil && (channel == anchan) && not killelHEM) {
	//if (passZ && passh && sampleType == 0 && passFil) {//for Zee channel checks
	trimTree->Fill();
	countpass += 1;
      }

   }

   
   std::cout<<"Passing Trigger req:        "<<counttrigpass<<std::endl;
   std::cout<<"Passing Z  req, w/ pT cut:  "<<countzpass<<std::endl;
   std::cout<<"Passing h  req:             "<<counthpass<<std::endl;
   std::cout<<"Passing    req:             "<<countpass<<std::endl;
   std::cout<<"Passing number of fat jets: "<<countfat<<std::endl;
   //std::cout<<"Had a Z candidate:          "<<countzcand<<std::endl;
   std::cout<<"had a leading muon:         "<<emumulead<<std::endl;
   std::cout<<"had a leading electron:     "<<emuelead<<std::endl;
   //std::cout<<"Straight saving of Zemu "<<emcounter<<std::endl;
   

   ///*
   std::cout<<"Events with zmumu        "<<zmumu<<std::endl;
   std::cout<<"Events with zmumuzee     "<<zmumuzee<<std::endl;
   std::cout<<"Events with zmumuzemu    "<<zmumuzemu<<std::endl;
   std::cout<<"Events with zmumuzeezemu "<<zmumuzeezemu<<std::endl;
   std::cout<<"Events with zee          "<<zee<<std::endl;
   std::cout<<"Events with zeezemu      "<<zeezemu<<std::endl;

   std::cout<<"Number of events in skims:                  "<<nentries<<std::endl;
   std::cout<<"Number of events in skims qcdweighted up:   "<<qcdwskimup<<std::endl;
   std::cout<<"Number of events in skims qcdweighted dn:   "<<qcdwskimdwn<<std::endl;
//std::cout<<"Events with zemu         "<<zemu<<std::endl;

   /*
   std::cout<<"Events with no reclustered fat jets, but a Z "<<znorecfat<<std::endl;
   std::cout<<"Events with no orignal  fat jets, but a Z    "<<znounrecfat<<std::endl;
   std::cout<<"Events with no fat jets, but a Z             "<<znofat<<std::endl;
   */
   htrigpass->SetBinContent(1,counttrigpass);
   hZpass->SetBinContent(1,countzpass);
   hHpass->SetBinContent(1,counthpass);
   hpass->SetBinContent(1,countpass);
   hnskimedup->SetBinContent(1,qcdwskimup);
   hnskimeddwn->SetBinContent(1,qcdwskimdwn);
   hnpu->SetBinContent(1,puweighttot);
   hnpuup->SetBinContent(1,puweightuptot);
   hnpudwn->SetBinContent(1,puweightdwntot);
   hnpunumup->SetBinContent(1,puweightvtxuptot);
   hnpunumdwn->SetBinContent(1,puweightvtxdntot);
   
   trimFile->Write();
   trimFile->Close();   
   //std::cout<<"trimmed to "<< passEvents <<" events"<<std::endl;
   std::cout<<"Completed your topiary garden, hopefully your tastes have not changed."<<std::endl;

   
}

