#include <iostream>
#include <vector>
#include "TLorentzVector.h"

#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include "RooFit.h"

//#include <TError.h>
#include <ROOT/RDataFrame.hxx>

using namespace std;
using namespace ROOT;
//using namespace RooFit;

//Particle masses
const float muon_mass = 0.105700;  //0.1057128
const float pion_mass = 0.1395263;
const float ks_mass   = 0.497611;
const float Jpsi_mass = 3.096900;
//muon cuts
const float muon_eta_cut  = 2.2;  
const float muon_pt_cut   = 4.0;

// dimuon cuts: 0=jpsi, 1=psi2s, 2=upsns
const float   m_min_[3] = {    2.9,    3.35,    8.5 };
const float   m_max_[3] = {    3.3,    4.05,   11.5 };
const int     m_bin_[3] = {     40,      70,    300 };
const float   t_cut_[3] = {    25.,     14.,    10. };
const float   y_cut_[3] = {    1.2,     1.2,    1.2 };
const int     i_hlt_[3] = {   4096,      32,      2 };



const UInt_t Max_nMuon = 25;
const UInt_t Max_nMM = 250;


void Jmakedataset_RDF(UInt_t iset=0)
{  
  bool usePV =false;
  //gErrorIgnoreLevel = 2001;
  //gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");// Ignore Warnings
  ROOT::EnableImplicitMT(4);// Tell ROOT you want to go parallel
  
  TChain tree("treeS");
  tree.Add("ROOT_slim_Bdks0_2023MiniAOD.root/treeS");

  ROOT::RDataFrame dDF(tree);
  auto nentries = dDF.Count();
  cout << " Total entries in Events " << *nentries <<endl;

////////////////////////////J/Psi filters
  // Estos filtros me parece que ya están establecidos; sin embargo, me parecio buena idea agregarlos por práctica.

  auto soft_req = "softmu1==1 && softmu2 ==1";
  auto Jpsi_window = Form("massJ<%1.6f+0.010000 && massJ>%1.6f-0.010000",Jpsi_mass,Jpsi_mass);
  auto muon_filter = Form("mu1pt>%1.1f && mu2pt>%1.1f && TMath::Abs(mu1eta)<%1.1f && TMath::Abs(mu2eta)<%1.1f",muon_pt_cut,muon_pt_cut,muon_eta_cut,muon_eta_cut);
 
  auto dDFJpsi = dDF.Filter(soft_req,"soft_req").Filter(Jpsi_window,"Jpsi_window").Filter(muon_filter,"muon_filter");


////////////////////////////////////// B0d filters

  auto trigger = "tri_DMu4_LM_Displaced==1";
  auto ks_window = Form("masslamb<%1.6f+0.010000 && masslamb>%1.6f-0.010000",ks_mass,ks_mass);
  auto Jptfilter = "Jpsipt>8.0";
  auto Bptfilter = "Bpt>10.0";
  auto probfilter = "Bprob > 0.1";

  auto dDFBd = dDFJpsi.Filter(trigger,"trigger").Filter(ks_window,"ks_window").Filter(Jptfilter,"Jptfilter").Filter(Bptfilter,"Bptfilter").Filter(probfilter,"probfilter");
  
  

//  auto selDimuon = Form("mu1pt>%1.2f && mu2pt>%1.2f && mu1pt<%1.2f && mu2pt<2.5",t_cut_[iset],m_min_[iset],m_max_[iset]);
 
//  auto seldata = dDF.Filter(selDimuon, "selDimuon");
 
  const vector<string> columns  = {"mass", "pt", "y","ct","cterr","mu1pt","mu1y","mu2pt","mu2y"};
  TString ofile = "jpsi-rootuple_RDF.root";
  if (iset==1) ofile = "psi2s-rootuple_RDF.root";
  if (iset==2) ofile = "upsns-rootuple_RDF.root";
  dDFBd.Snapshot("dimuonTree", ofile, columns);
  
}
