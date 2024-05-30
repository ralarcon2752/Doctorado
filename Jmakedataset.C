#include <iostream>
#include <vector>
#include "TLorentzVector.h"

#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>

//#include <TError.h>
//#include <ROOT/RDataFrame.hxx>

using namespace std;

const float muon_mass = 0.105700;  //0.1057128;
const float pion_mass = 0.1395263;

// dimuon cuts: 0=jpsi, 1=psi2s, 2=upsns
const float   m_min_[3] = {    2.9,    3.35,    8.5 };
const float   m_max_[3] = {    3.3,    4.05,   11.5 };
const int     m_bin_[3] = {     40,      70,    300 };
const float   t_cut_[3] = {    25.,     14.,    10. };
const float   y_cut_[3] = {    1.2,     1.2,    1.2 };
const int     i_hlt_[3] = {   4096,      32,      2 };
//const TString i_tag_[3] = { 'jpsi', 'psi2s', 'upsns'};
//my triggers = { 'HLT_Dimuon25_Jpsi_v*', 'HLT_Dimuon14_PsiPrime_v*', 'HLT_Dimuon10_Upsilon_y1p4_v*'};
//https://github.com/cms-sw/cmssw/tree/CMSSW_12_4_X/HLTrigger/Configuration/python

// single muon acc cuts
const float muon_eta_cut  = 2.4;  
const float muon_pt_cut   = 3.0;

const UInt_t Max_nMuon = 25;
const UInt_t Max_nMM = 250;

// charge correlation is assumed, e.g. charge_j+cherge_l = 0
bool is_seagull_(Float_t phi_j, Float_t phi_l, Float_t charge_j) {
  bool seagull_   = true;
  Double_t phiPN = phi_j - phi_l;
  if (charge_j < 0)  phiPN = phi_l - phi_j;
  if (phiPN > TMath::Pi()) phiPN -= 2.*TMath::Pi();
  else if (phiPN <= -TMath::Pi()) phiPN += 2.*TMath::Pi();
  if (phiPN>0.) seagull_ = false;
  return seagull_;
}

void Jmakedataset(Long64_t nbreak=0, UInt_t iset=0) {
  
  bool usePV =false;
  //gErrorIgnoreLevel = 2001;
  //gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");// Ignore Warnings
  //ROOT::EnableImplicitMT(8);// Tell ROOT you want to go parallel
  
  TChain *tree = new TChain("mm_tree");
  //tree->Add("DataMiniAOD/*.root/rootuple/mm_tree");
  //test
  //tree->Add("DataMiniAOD/ParkingDoubleMuonLowMass0_1.root/rootuple/mm_tree");
  tree->Add("../Data/ParkingDoubleMuonLowMass0_1.root/rootuple/mm_tree");

  Long64_t nentries = tree->GetEntries();
  std::cout << "Total entries in Events "  << nentries << std::endl;

  TLorentzVector* dimuon_p4 = 0;
  TLorentzVector* muonP_p4 = 0;
  TLorentzVector* muonM_p4 = 0;
  Float_t         vProb;
  Float_t         ppdlPV;
  Float_t         ppdlErrPV;
  UInt_t          trigger;
  Int_t           charge;
  Float_t         ppdl;
  Float_t         ppdlErr;
  
  tree->SetBranchAddress("dimuon_p4",&dimuon_p4);
  tree->SetBranchAddress("muonP_p4",&muonP_p4);
  tree->SetBranchAddress("muonM_p4",&muonM_p4);
  tree->SetBranchAddress("vProb",&vProb);
  tree->SetBranchAddress("trigger",&trigger);
  tree->SetBranchAddress("charge",&charge);
  
  // setting proper ppdl branches
  if (usePV) {
    tree->SetBranchAddress("ppdlPV",&ppdl);
    tree->SetBranchAddress("ppdlErrPV",&ppdlErr);
  } else {
    tree->SetBranchAddress("ppdlBS",&ppdl);
    tree->SetBranchAddress("ppdlErrBS",&ppdlErr);
  }
  
  Long64_t nevents = 0;
  Long64_t npassed = 0;
  
  Float_t mass,pt,y,ct,cterr;
  Float_t mu1pt,mu2pt,mu1y,mu2y;
  
  TString ofile = "jpsi-rootuple.root";
  if (iset==1) ofile = "psi2s-rootuple.root";
  if (iset==2) ofile = "upsns-rootuple.root";
  TFile *outfile = new TFile(ofile,"recreate");
  TTree *dimuontree = new TTree("dimuonTree", "HLT Dimuon Tree");
  
  dimuontree->Branch("mass", &mass, "mass/F");
  dimuontree->Branch("pt",   &pt,   "pt/F");
  dimuontree->Branch("y",    &y,    "y/F");
  dimuontree->Branch("ct",   &ct,   "ct/F");
  dimuontree->Branch("cterr",&cterr,"cterr/F");
  dimuontree->Branch("mu1pt",   &mu1pt,   "mu1pt/F");
  dimuontree->Branch("mu1y",    &mu1y,    "mu1y/F");
  dimuontree->Branch("mu2pt",   &mu2pt,   "mu2pt/F");
  dimuontree->Branch("mu2y",    &mu2y,    "mu2y/F");
  
  TH1F *control_m = new TH1F("control_m", "control plot;M(#mu^{+}#mu^{-}) [GeV];Events per 10 MeV",m_bin_[iset],m_min_[iset],m_max_[iset]);
  TH2F *h2mu1 = new TH2F("h2mu1", "h2mu1",48,-2.4,2.4,86,3.0,50.0);
  TH2F *h2mu2 = new TH2F("h2mu2", "h2mu2",48,-2.4,2.4,86,3.0,50.0);
  
  //Long64_t is_good = 0;
  for (Long64_t k=0;k<nentries; k++) {
    tree->GetEntry(k);
    if ( (k%(nentries/10)) == 0 ) std::cout << k << " out of " << nentries << " ..." << std::endl;
    
    nevents++;
    if(nevents==nbreak)  break;

    if(TMath::IsNaN(ppdl)==1)    continue;
    if(TMath::IsNaN(ppdlErr)==1) continue;
    
    if (!((trigger&i_hlt_[iset])>0))          continue;
    if (dimuon_p4->Pt() < t_cut_[iset])       continue;
    if (dimuon_p4->M()<m_min_[iset])          continue;
    if (dimuon_p4->M()>m_max_[iset])          continue;
    if (charge!=0)                            continue;
    if (vProb<0.01)                           continue;
    if (TMath::Abs(dimuon_p4->Rapidity())>2.5) continue;
    if (muonP_p4->Pt()<muon_pt_cut || TMath::Abs(muonP_p4->Eta())>muon_eta_cut) continue;
    if (muonM_p4->Pt()<muon_pt_cut || TMath::Abs(muonM_p4->Eta())>muon_eta_cut) continue;
    
    // ALGUNA CONDICION PARA "SEAGUL"????????????????????????
    Float_t muchargej = 0;
    Float_t muphij = 0.0; Float_t muphil=0.0;
    if (muonP_p4->Pt() > muonM_p4->Pt() ) {
      mu1pt  = muonP_p4->Pt();
      mu1y   = muonP_p4->Eta();
      mu2pt = muonM_p4->Pt();
      mu2y   = muonM_p4->Eta();
      muchargej = 1;
      muphij = muonP_p4->Phi();
      muphil = muonM_p4->Phi();
    } else {
      mu1pt  = muonM_p4->Pt();
      mu1y   = muonM_p4->Eta();
      mu2pt  = muonP_p4->Pt();
      mu2y   = muonP_p4->Eta();
      muchargej = -1;
      muphij = muonM_p4->Phi();
      muphil = muonP_p4->Phi();
    }
    if (!is_seagull_(muphij,muphil,muchargej))  continue;
    
    mass  = dimuon_p4->M();
    pt    = dimuon_p4->Pt();
    y     = dimuon_p4->Rapidity();
    ct    = ppdl;
    cterr = ppdlErr;
    
    h2mu1->Fill(muonP_p4->Eta(),muonP_p4->Pt());
    h2mu2->Fill(muonM_p4->Eta(),muonM_p4->Pt());
    control_m->Fill(mass);
    dimuontree->Fill();
    npassed++;
  }
  
  h2mu1->Write();
  h2mu2->Write();
  control_m->Write();
  dimuontree->Write();
  outfile->Close();
  delete outfile;
  
  printf("total of analyzed events: %lli, total of stored: %lli \n", nevents, npassed);
}
