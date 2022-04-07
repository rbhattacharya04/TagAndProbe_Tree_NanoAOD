import ROOT
import os
import numpy as np
#ROOT.ROOT.EnableImplicitMT()

#verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)


#tag_Muon_code = '''
#using namespace ROOT::VecOps;
#int getTagMuon(UInt_t& nMuon, ROOT::VecOps::RVec<Float_t> &pt, ROOT::VecOps::RVec<Float_t> &eta, ROOT::VecOps::RVec<Bool_t> &isGlobal){
#  ROOT::VecOps::RVec<Int_t> tagMuons = -1; 
#  for(int i=0; i<nMuon; i++){
#    if(isGlobal[i] && pt[i]> 20 && eta[i] <2){
#	tagMuons.push_back(i);
#    }
#  }
#  return tag_Muon_index;
#}
#'''
#ROOT.gInterpreter.Declare(tag_Muon_code)


#probe_Muon_code = '''
#using namespace ROOT::VecOps;
#int getProbeMuons(UInt_t& nMuon, int tag_Muon, ROOT::VecOps::RVec<Float_t> &pt, ROOT::VecOps::RVec<Bool_t> &isTracker, ROOT::VecOps::RVec<Int_t> &charge){
#   ROOT::VecOps::RVec<Int_t> probeMuons;
#   for(int i=0; i<nMuon; i++){
#     if (i == tag_Muon) continue;
#     //if ((charge[tag_Muon] + charge[i]) != 0) continue;
#     if (pt[i] > 10 && isTracker[i]){
#	 probeMuons.push_back(i);
#     }
#   }
#   return probeMuons;
#}
#'''
#ROOT.gInterpreter.Declare(probe_Muon_code)

create_TP_pair_code = '''
using namespace ROOT::VecOps;


struct TPPair{
  ROOT::Math::PtEtaPhiMVector LorentzVector;
  
  int tag_charge;
  float tag_pt;
  float tag_ptErr;        
  float tag_eta;
  float tag_phi;
  UChar_t tag_cleanmask;
  float tag_dxy;
  float tag_dxyErr;
  float tag_dz;
  float tag_dzErr;
  UChar_t tag_genPartFlav;
  int tag_genPartIdx;
  UChar_t tag_highPtId;
  bool tag_inTimeMuon;
  float tag_ip3d;
  bool tag_isGlobal;
  bool tag_isPFcand;
  bool tag_isTracker;
  int tag_jetIdx;
  float tag_jetPtRelv2;
  float tag_jetRelIso;
  float tag_mass;
  bool tag_mediumId;
  bool tag_mediumPromptId;
  UChar_t tag_miniIsoId;
  float tag_miniPFRelIso_all;
  float tag_miniPFRelIso_chg;
  UChar_t tag_multiIsoId;
  UChar_t tag_mvaId;
  float tag_mvaTTH;
  int tag_nStations;
  int tag_nTrackerLayers;
  int tag_pdgId;
  UChar_t tag_pfIsoId;
  float tag_pfRelIso03_all;
  float tag_pfRelIso03_chg;
  float tag_pfRelIso04_all;
  float tag_segmentComp;
  float tag_sip3d;
  bool tag_softId;
  bool tag_softMvaId;
  int tag_tightCharge;
  bool tag_tightId;
  UChar_t tag_tkIsoId;
  bool tag_triggerIdLoose;


  int probe_charge;
  float probe_pt;
  float probe_ptErr;        
  float probe_eta;
  float probe_phi;
  UChar_t probe_cleanmask;
  float probe_dxy;
  float probe_dxyErr;
  float probe_dz;
  float probe_dzErr;
  UChar_t probe_genPartFlav;
  int probe_genPartIdx;
  UChar_t probe_highPtId;
  bool probe_inTimeMuon;
  float probe_ip3d;
  bool probe_isGlobal;
  bool probe_isPFcand;
  bool probe_isTracker;
  int probe_jetIdx;
  float probe_jetPtRelv2;
  float probe_jetRelIso;
  float probe_mass;
  bool probe_mediumId;
  bool probe_mediumPromptId;
  UChar_t probe_miniIsoId;
  float probe_miniPFRelIso_all;
  float probe_miniPFRelIso_chg;
  UChar_t probe_multiIsoId;
  UChar_t probe_mvaId;
  float probe_mvaTTH;
  int probe_nStations;
  int probe_nTrackerLayers;
  int probe_pdgId;
  UChar_t probe_pfIsoId;
  float probe_pfRelIso03_all;
  float probe_pfRelIso03_chg;
  float probe_pfRelIso04_all;
  float probe_segmentComp;
  float probe_sip3d;
  bool probe_softId;
  bool probe_softMvaId;
  int probe_tightCharge;
  bool probe_tightId;
  UChar_t probe_tkIsoId;
  bool probe_triggerIdLoose;
};

ROOT::VecOps::RVec<TPPair> createTPPair(ROOT::VecOps::RVec<Int_t> &tag_Muon_charge,
					ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector> &tag_4_vector,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_ptErr,
                                        ROOT::VecOps::RVec<UChar_t> &tag_Muon_cleanmask,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_dxy,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_dxyErr,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_dz,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_dzErr,
                                        ROOT::VecOps::RVec<UChar_t> &tag_Muon_genPartFlav,
                                        ROOT::VecOps::RVec<Int_t> &tag_Muon_genPartIdx,
                                        ROOT::VecOps::RVec<UChar_t> &tag_Muon_highPtId,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_inTimeMuon,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_ip3d,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_isGlobal,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_isPFcand,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_isTracker,
                                        ROOT::VecOps::RVec<Int_t> &tag_Muon_jetIdx,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_jetPtRelv2,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_jetRelIso,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_mediumId,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_mediumPromptId,
                                        ROOT::VecOps::RVec<UChar_t> &tag_Muon_miniIsoId,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_miniPFRelIso_all,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_miniPFRelIso_chg,
                                        ROOT::VecOps::RVec<UChar_t> &tag_Muon_multiIsoId,
                                        ROOT::VecOps::RVec<UChar_t> &tag_Muon_mvaId,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_mvaTTH,
                                        ROOT::VecOps::RVec<Int_t> &tag_Muon_nStations,
                                        ROOT::VecOps::RVec<Int_t> &tag_Muon_nTrackerLayers,
                                        ROOT::VecOps::RVec<Int_t> &tag_Muon_pdgId,
                                        ROOT::VecOps::RVec<UChar_t> &tag_Muon_pfIsoId,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_pfRelIso03_all,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_pfRelIso03_chg,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_pfRelIso04_all,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_segmentComp,
                                        ROOT::VecOps::RVec<Float_t> &tag_Muon_sip3d,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_softId,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_softMvaId,
                                        ROOT::VecOps::RVec<Int_t> &tag_Muon_tightCharge,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_tightId,
                                        ROOT::VecOps::RVec<UChar_t> &tag_Muon_tkIsoId,
                                        ROOT::VecOps::RVec<Bool_t> &tag_Muon_triggerIdLoose,
					ROOT::VecOps::RVec<Int_t> &probe_Muon_charge,
                                        ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector> &probe_4_vector,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_ptErr,
                                        ROOT::VecOps::RVec<UChar_t> &probe_Muon_cleanmask,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_dxy,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_dxyErr,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_dz,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_dzErr,
                                        ROOT::VecOps::RVec<UChar_t> &probe_Muon_genPartFlav,
                                        ROOT::VecOps::RVec<Int_t> &probe_Muon_genPartIdx,
                                        ROOT::VecOps::RVec<UChar_t> &probe_Muon_highPtId,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_inTimeMuon,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_ip3d,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_isGlobal,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_isPFcand,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_isTracker,
                                        ROOT::VecOps::RVec<Int_t> &probe_Muon_jetIdx,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_jetPtRelv2,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_jetRelIso,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_mediumId,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_mediumPromptId,
                                        ROOT::VecOps::RVec<UChar_t> &probe_Muon_miniIsoId,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_miniPFRelIso_all,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_miniPFRelIso_chg,
                                        ROOT::VecOps::RVec<UChar_t> &probe_Muon_multiIsoId,
                                        ROOT::VecOps::RVec<UChar_t> &probe_Muon_mvaId,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_mvaTTH,
                                        ROOT::VecOps::RVec<Int_t> &probe_Muon_nStations,
                                        ROOT::VecOps::RVec<Int_t> &probe_Muon_nTrackerLayers,
                                        ROOT::VecOps::RVec<Int_t> &probe_Muon_pdgId,
                                        ROOT::VecOps::RVec<UChar_t> &probe_Muon_pfIsoId,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_pfRelIso03_all,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_pfRelIso03_chg,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_pfRelIso04_all,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_segmentComp,
                                        ROOT::VecOps::RVec<Float_t> &probe_Muon_sip3d,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_softId,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_softMvaId,
                                        ROOT::VecOps::RVec<Int_t> &probe_Muon_tightCharge,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_tightId,
                                        ROOT::VecOps::RVec<UChar_t> &probe_Muon_tkIsoId,
                                        ROOT::VecOps::RVec<Bool_t> &probe_Muon_triggerIdLoose){
  int n_tag = tag_Muon_charge.size();
  int n_probe = probe_Muon_charge.size();
  ROOT::VecOps::RVec<TPPair> TPPairs;
  for(size_t i=0;i<n_tag;i++){
    for(size_t j=0;j<n_probe;j++){
      if(tag_Muon_charge[i] == probe_Muon_charge[j]) continue;
      ROOT::Math::PtEtaPhiMVector TPPair_LorentzVector = tag_4_vector[i] + probe_4_vector[j];
      

      TPPair tp;
      tp.LorentzVector = TPPair_LorentzVector;
      
      tp.tag_charge = tag_Muon_charge[i];
      tp.tag_pt = tag_4_vector[i].Pt();
      tp.tag_ptErr = tag_Muon_ptErr[i];        
      tp.tag_eta = tag_4_vector[i].Eta();
      tp.tag_phi = tag_4_vector[i].Phi();
      tp.tag_cleanmask = tag_Muon_cleanmask[i];
      tp.tag_dxy = tag_Muon_dxy[i];
      tp.tag_dxyErr = tag_Muon_dxyErr[i];
      tp.tag_dz = tag_Muon_dz[i];
      tp.tag_dzErr = tag_Muon_dzErr[i];
      tp.tag_genPartFlav = tag_Muon_genPartFlav[i];
      tp.tag_genPartIdx = tag_Muon_genPartIdx[i];
      tp.tag_highPtId = tag_Muon_highPtId[i];
      tp.tag_inTimeMuon = tag_Muon_inTimeMuon[i];
      tp.tag_ip3d = tag_Muon_ip3d[i];
      tp.tag_isGlobal = tag_Muon_isGlobal[i];
      tp.tag_isPFcand = tag_Muon_isPFcand[i];
      tp.tag_isTracker = tag_Muon_isTracker[i];
      tp.tag_jetIdx = tag_Muon_jetIdx[i];
      tp.tag_jetPtRelv2 = tag_Muon_jetPtRelv2[i];
      tp.tag_jetRelIso = tag_Muon_jetRelIso[i];
      tp.tag_mass = tag_4_vector[i].M();
      tp.tag_mediumId = tag_Muon_mediumId[i];
      tp.tag_mediumPromptId = tag_Muon_mediumPromptId[i];
      tp.tag_miniIsoId = tag_Muon_miniIsoId[i];
      tp.tag_miniPFRelIso_all = tag_Muon_miniPFRelIso_all[i];
      tp.tag_miniPFRelIso_chg = tag_Muon_miniPFRelIso_chg[i] ;
      tp.tag_multiIsoId = tag_Muon_multiIsoId[i];
      tp.tag_mvaId = tag_Muon_mvaId[i];
      tp.tag_mvaTTH = tag_Muon_mvaTTH[i];
      tp.tag_nStations = tag_Muon_nStations[i];
      tp.tag_nTrackerLayers = tag_Muon_nTrackerLayers[i];
      tp.tag_pdgId = tag_Muon_pdgId[i];
      tp.tag_pfIsoId = tag_Muon_pfIsoId[i];
      tp.tag_pfRelIso03_all = tag_Muon_pfRelIso03_all[i];
      tp.tag_pfRelIso03_chg = tag_Muon_pfRelIso03_chg[i];
      tp.tag_pfRelIso04_all = tag_Muon_pfRelIso04_all[i];
      tp.tag_segmentComp = tag_Muon_segmentComp[i];
      tp.tag_sip3d = tag_Muon_sip3d[i];
      tp.tag_softId = tag_Muon_softId[i];
      tp.tag_softMvaId = tag_Muon_softMvaId[i];
      tp.tag_tightCharge = tag_Muon_tightCharge[i];
      tp.tag_tightId = tag_Muon_tightId[i];
      tp.tag_tkIsoId = tag_Muon_tkIsoId[i];
      tp.tag_triggerIdLoose = tag_Muon_triggerIdLoose[i];


      tp.probe_charge = probe_Muon_charge[j];
      tp.probe_pt = probe_4_vector[j].Pt();
      tp.probe_ptErr = probe_Muon_ptErr[j];        
      tp.probe_eta = probe_4_vector[j].Eta();
      tp.probe_phi = probe_4_vector[j].Phi();
      tp.probe_cleanmask = probe_Muon_cleanmask[j];
      tp.probe_dxy = probe_Muon_dxy[j];
      tp.probe_dxyErr = probe_Muon_dxyErr[j];
      tp.probe_dz = probe_Muon_dz[j];
      tp.probe_dzErr = probe_Muon_dzErr[j];
      tp.probe_genPartFlav = probe_Muon_genPartFlav[j];
      tp.probe_genPartIdx = probe_Muon_genPartIdx[j];
      tp.probe_highPtId = probe_Muon_highPtId[j];
      tp.probe_inTimeMuon = probe_Muon_inTimeMuon[j];
      tp.probe_ip3d = probe_Muon_ip3d[j];
      tp.probe_isGlobal = probe_Muon_isGlobal[j];
      tp.probe_isPFcand = probe_Muon_isPFcand[j];
      tp.probe_isTracker = probe_Muon_isTracker[j];
      tp.probe_jetIdx = probe_Muon_jetIdx[j];
      tp.probe_jetPtRelv2 = probe_Muon_jetPtRelv2[j];
      tp.probe_jetRelIso = probe_Muon_jetRelIso[j];
      tp.probe_mass = probe_4_vector[j].M();
      tp.probe_mediumId = probe_Muon_mediumId[j];
      tp.probe_mediumPromptId = probe_Muon_mediumPromptId[j];
      tp.probe_miniIsoId = probe_Muon_miniIsoId[j];
      tp.probe_miniPFRelIso_all = probe_Muon_miniPFRelIso_all[j];
      tp.probe_miniPFRelIso_chg = probe_Muon_miniPFRelIso_chg[j] ;
      tp.probe_multiIsoId = probe_Muon_multiIsoId[j];
      tp.probe_mvaId = probe_Muon_mvaId[j];
      tp.probe_mvaTTH = probe_Muon_mvaTTH[j];
      tp.probe_nStations = probe_Muon_nStations[j];
      tp.probe_nTrackerLayers = probe_Muon_nTrackerLayers[j];
      tp.probe_pdgId = probe_Muon_pdgId[j];
      tp.probe_pfIsoId = probe_Muon_pfIsoId[j];
      tp.probe_pfRelIso03_all = probe_Muon_pfRelIso03_all[j];
      tp.probe_pfRelIso03_chg = probe_Muon_pfRelIso03_chg[j];
      tp.probe_pfRelIso04_all = probe_Muon_pfRelIso04_all[j];
      tp.probe_segmentComp = probe_Muon_segmentComp[j];
      tp.probe_sip3d = probe_Muon_sip3d[j];
      tp.probe_softId = probe_Muon_softId[j];
      tp.probe_softMvaId = probe_Muon_softMvaId[j];
      tp.probe_tightCharge = probe_Muon_tightCharge[j];
      tp.probe_tightId = probe_Muon_tightId[j];
      tp.probe_tkIsoId = probe_Muon_tkIsoId[j];
      tp.probe_triggerIdLoose = probe_Muon_triggerIdLoose[j];
      

      // std::pair<size_t,size_t> TPPair = std::make_pair(i,j);
      //TPPairs.push_back(std::make_pair(TPPair,TPPair_LorentzVector));  
      TPPairs.push_back(tp);  
    }
  }
  //return ROOT::VecOps::Map(tag_Muon_charge,tag_4_vector,probe_Muon_charge,probe_4_vector,TPPairs);
  return TPPairs;
} 

'''



ROOT.gInterpreter.Declare(create_TP_pair_code)

test_code = '''

using namespace ROOT::VecOps;
double test(ROOT::VecOps::RVec<Int_t> &tag_Muon_charge,ROOT::VecOps::RVec<Int_t> &probe_Muon_charge){
  std::cout<<"tag_Muon_charge[0] = "<<tag_Muon_charge[0]<<std::endl;  
  return 5;
}
'''

ROOT.gInterpreter.Declare(test_code)

get_mass_code = '''
using namespace ROOT::VecOps;
using TPSave = std::pair<std::pair<size_t,size_t>,ROOT::Math::PtEtaPhiMVector>;
ROOT::VecOps::RVec<Float_t> getMass(ROOT::VecOps::RVec<TPSave> TPPairs){
  ROOT::VecOps::RVec<Float_t> mass;
  for (auto m:TPPairs) mass.push_back(m.second.M());
  return mass; 
}

'''
ROOT.gInterpreter.Declare(get_mass_code)


get_Tags_code = '''
using namespace ROOT::VecOps;
using TPSave = std::pair<std::pair<size_t,size_t>,ROOT::Math::PtEtaPhiMVector>;
ROOT::VecOps::RVec<Int_t> getTags(ROOT::VecOps::RVec<TPSave> TPPairs){
  ROOT::VecOps::RVec<Int_t> Tags;
  for (auto m:TPPairs) Tags.push_back(m.first.first);
  return Tags; 
}

'''
ROOT.gInterpreter.Declare(get_Tags_code)


get_Probes_code = '''
using namespace ROOT::VecOps;
using TPSave = std::pair<std::pair<size_t,size_t>,ROOT::Math::PtEtaPhiMVector>;
ROOT::VecOps::RVec<Int_t> getProbes(ROOT::VecOps::RVec<TPSave> TPPairs){
  ROOT::VecOps::RVec<Int_t> Probes;
  for (auto m:TPPairs) Probes.push_back(m.first.second);
  return Probes; 
}

'''
ROOT.gInterpreter.Declare(get_Probes_code)

get_pt_code = '''
using namespace ROOT::VecOps;
float getPt(ROOT::VecOps::RVec<Int_t> &tag_Muon_charge, ROOT::VecOps::RVec<Int_t> &probe_Muon_charge, int nMuon){
   std::cout<<"getPt"<<std::endl;
   float pt = 0;
   int n_tag = tag_Muon_charge.size();
   int n_probe = probe_Muon_charge.size();
   std::cout<<"n_tag = "<<n_tag<<" n_probe = "<<n_probe<<" nMuon = "<<nMuon<<std::endl;
   //std::cout<<"tag_Muon_charge[0] = "<<tag_Muon_charge[0]<<std::endl; 
   for(int i=0;i<n_tag;i++){ 
     std::cout<<"i = "<<i<<" tag_Muon_charge[i] = " <<tag_Muon_charge[i]<<std::endl;
     pt = tag_Muon_charge[i];
   } 
   return pt; 
}

'''
ROOT.gInterpreter.Declare(get_pt_code)

flatten_code = '''
using namespace ROOT::VecOps;
template<typename T>
ROOT::VecOps::RVec<T> flatten(std::vector<ROOT::VecOps::RVec<T>> const &vec)
{
  int size = 0;
  for (auto &v: vec) {
    size += v.size();
  }

  ROOT::VecOps::RVec<T> flattened;
  flattened.reserve(size);
  for (auto const &v: vec) {
    for (auto &e: v) flattened.push_back(e);        
  }
  return flattened;
}

'''
ROOT.gInterpreter.Declare(flatten_code)

makeTPDF_code = '''
using namespace ROOT::VecOps;
template<typename T>
ROOT::RDataFrame makeTPDF(ROOT::RDataFrame df){
  auto TPPairs = df.Take<ROOT::VecOps::RVec<TPPair>>("TPPairs");
  auto TPPairs_vector = *TPPairs;
  auto TPPairs_flatten = flatten(TPPairs_vector);
  auto getTP = [=] (int i) {return TPPairs_flatten[i];};
  ROOT::RDataFrame TPdf(TPPairs_flatten.size());
  TPdf.Define("TPPairs",getTP);
  return TPdf;   
}
'''
ROOT.gInterpreter.Declare(makeTPDF_code)

dirname="/afs/cern.ch/work/r/rbhattac/public/TnPTool_Tutorial/CMSSW_10_6_2/src/MuonAnalysis/TagAndProbe/test/zmumu/Sample_NanoAodFiles"

from os import listdir
from os.path import isfile,join

files = [join(dirname, f) for f in listdir(dirname) if f.endswith(".root")]


filenames = ROOT.std.vector('string')()

for name in files: filenames.push_back(name)


d = ROOT.RDataFrame("Events", filenames)

tag_selection = "Muon_isGlobal && Muon_pt >20 && abs(Muon_eta)<2"
probe_selection = "Muon_isTracker && Muon_pt >10"
masscut = "40 < mass < 200"


colNames = d.GetColumnNames();

for colName in colNames:
        colName = str(colName)
	if (colName.find('Muon', 0, 4) != -1):
             d = d.Define("tag_"+colName,colName+"[("+tag_selection+")]")
             d = d.Define("probe_"+colName,colName+"[("+probe_selection+")]")

d = d.Define("n_tag","tag_Muon_pt.size()")
d = d.Define("n_probe","probe_Muon_pt.size()")
d = d.Define("tag_4_vector","ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(tag_Muon_pt,tag_Muon_eta,tag_Muon_phi,tag_Muon_mass)").Define("probe_4_vector","ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(probe_Muon_pt,probe_Muon_eta,probe_Muon_phi,probe_Muon_mass)")
d = d.Define("TPPairs","createTPPair(tag_Muon_charge,tag_4_vector,tag_Muon_ptErr,tag_Muon_cleanmask,tag_Muon_dxy,tag_Muon_dxyErr,tag_Muon_dz,tag_Muon_dzErr,tag_Muon_genPartFlav,tag_Muon_genPartIdx,tag_Muon_highPtId,tag_Muon_inTimeMuon,tag_Muon_ip3d,tag_Muon_isGlobal,tag_Muon_isPFcand,tag_Muon_isTracker,tag_Muon_jetIdx,tag_Muon_jetPtRelv2,tag_Muon_jetRelIso,tag_Muon_mediumId,tag_Muon_mediumPromptId,tag_Muon_miniIsoId,tag_Muon_miniPFRelIso_all,tag_Muon_miniPFRelIso_chg,tag_Muon_multiIsoId,tag_Muon_mvaId,tag_Muon_mvaTTH,tag_Muon_nStations,tag_Muon_nTrackerLayers,tag_Muon_pdgId,tag_Muon_pfIsoId,tag_Muon_pfRelIso03_all,tag_Muon_pfRelIso03_chg,tag_Muon_pfRelIso04_all,tag_Muon_segmentComp,tag_Muon_sip3d,tag_Muon_softId,tag_Muon_softMvaId,tag_Muon_tightCharge,tag_Muon_tightId,tag_Muon_tkIsoId,tag_Muon_triggerIdLoose,probe_Muon_charge,probe_4_vector,probe_Muon_ptErr,probe_Muon_cleanmask,probe_Muon_dxy,probe_Muon_dxyErr,probe_Muon_dz,probe_Muon_dzErr,probe_Muon_genPartFlav,probe_Muon_genPartIdx,probe_Muon_highPtId,probe_Muon_inTimeMuon,probe_Muon_ip3d,probe_Muon_isGlobal,probe_Muon_isPFcand,probe_Muon_isTracker,probe_Muon_jetIdx,probe_Muon_jetPtRelv2,probe_Muon_jetRelIso,probe_Muon_mediumId,probe_Muon_mediumPromptId,probe_Muon_miniIsoId,probe_Muon_miniPFRelIso_all,probe_Muon_miniPFRelIso_chg,probe_Muon_multiIsoId,probe_Muon_mvaId,probe_Muon_mvaTTH,probe_Muon_nStations,probe_Muon_nTrackerLayers,probe_Muon_pdgId,probe_Muon_pfIsoId,probe_Muon_pfRelIso03_all,probe_Muon_pfRelIso03_chg,probe_Muon_pfRelIso04_all,probe_Muon_segmentComp,probe_Muon_sip3d,probe_Muon_softId,probe_Muon_softMvaId,probe_Muon_tightCharge,probe_Muon_tightId,probe_Muon_tkIsoId,probe_Muon_triggerIdLoose)")

#print(type(d))

TPPairs = d.AsNumpy(columns=["TPPairs"])
TPPairs_v = TPPairs["TPPairs"]

TPPairs_np = np.array([])

for e in TPPairs_v:
	for f in e:
        	TPPairs_np = np.append(TPPairs_np,f)

print(TPPairs_np)

mass = np.array([e.LorentzVector.M() for e in TPPairs_np])
pt = np.array([e.probe_pt for e in TPPairs_np])
TPDF = ROOT.RDF.MakeNumpyDataFrame({'mass': mass, 'pt': pt})

#TPDF = ROOT.RDF.MakeNumpyDataFrame(data)
#TPDF = ROOT.makeTPDF(d)
#TPDF = TPDF.Define("mass","TPPairs.LorentzVector.M()")

#print(d.GetColumnType("TPPairs"))

#test = d.Take<ROOT.VecOps.RVec<ROOT.TPPair>>("TPPairs");
#print(type(test))

#d = d.Define("mass","getMass(TPPairs)")
#d = d.Filter("All(mass > 50.) && All(mass < 140.)")
#d = d.Define("tags","getTags(TPPairs)").Define("tag_pt","tag_Muon_pt[tags]")
#d = d.Define("n_probe_2","getProbes(TPPairs).size()")
#d = d.Define("pt","getPt(tag_Muon_charge,probe_Muon_charge,nMuon)")
#d = d.Define("test","test(tag_Muon_charge,probe_Muon_charge)")

#d = d.Define("mass","ROOT::VecOps::InvariantMasses(tag_Muon_pt,tag_Muon_eta,tag_Muon_phi,tag_Muon_mass,probe_Muon_pt,probe_Muon_eta,probe_Muon_phi,probe_Muon_mass)")





#tagmuons=d.Define('tag_Muon',"getTagMuon(nMuon,Muon_pt,Muon_eta,Muon_isGlobal)").Filter("tag_Muon > -1")
#probemuons=tagmuons.Define('probe_Muon',"getProbeMuons(nMuon,tag_Muon,Muon_pt,Muon_isTracker,Muon_charge)")

#tagmuons=tagmuons.Define('tag_Muon_pt','Muon_pt[tag_Muon]')


#tagmuons = d.Filter("All(Muon_isGlobal) && All(Muon_pt>20.) && All(abs(Muon_eta)<2.)", "tagmuon selection")
#d = d.Filter('All(Muon_pt>20.)')
#probemuons = d.Filter("All(Muon_isTracker) && All(Muon_pt>10)")

print("A")
h1 = TPDF.Histo1D(("mass","maass",200,0,200),"mass");
#h2 = d.Histo1D(("mass2","maass2",200,0,200),"mass2")
#h4 = d.Histo1D(("n_tag","n_tag",10,0,10),"n_tag")
#h5 = d.Histo1D(("n_probe","n_probe",10,0,10),"n_probe")
#h6 = d.Histo1D(("n_tag_2","n_tag_2",10,0,10),"n_tag_2")
#h7 = d.Histo1D(("n_probe_2","n_probe_2",10,0,10),"n_probe_2")
h2=d.Histo1D(("pt","pt",100,0,100),"tag_Muon_pt")
#h2 = d1.Histo1D(("n_Muons_2","n_Muons",5,0,5),"nMuon");
#h2 = tagmuons.Histo1D(("tag_Muon","tag_Muon",20,-10,10),"tag_Muon");
#h3 = d.Histo1D(("probe_Muon_charge","probe_Muon_charge",10,-5,5),"probe_Muon_charge")
outHistFile = ROOT.TFile("test_histogram_2.root","RECREATE")
print("B")
h1.Write()
h2.Write()
#h3.Write()
#h4.Write()
#h5.Write()
#h6.Write()
#h7.Write()
#branchList = ROOT.vector('string')()
#for branchName in ["mass"]:
#    branchList.push_back(branchName)
#d = d.Snapshot<ROOT.Float_t>("test","test_tree.root",branchList) 
#test = d.Take<vector<float>("nMuon");
#print(Rtest.GetPtr()))
print("C")
outHistFile.Close()
TPDF.Snapshot("myTree","test_TP_tree.root")
