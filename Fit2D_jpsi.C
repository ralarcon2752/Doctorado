#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooDecay.h"
#include "RooGaussModel.h"

#include "TFile.h"
#include "TTree.h"


#include "RooNumIntConfig.h"

#include "MyUtilities.h"

void save_result(ofstream& salida, RooRealVar Ns, RooRealVar Nb)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError();
  cout << " el archivo se escribio bien" << endl; 
  return;
}

void Fit2D_jpsi(float ptmin=25.0, float ptmax=30.0, float rapmin=0.0, float rapmax=0.3){
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);
  
  //TChain *_nt = new TChain("dimuonTree","");
  //_nt->Add("jpsi-rootuple.root");
  TFile *fin = new TFile("jpsi-rootuple.root");
  //TFile *fin = new TFile("jpsi-rootuplemy.root");
  //TFile *fin = new TFile("MyPyRootTry/jpsi-rootuple.root");
  //TFile *fin = new TFile("MyPyRootTry/jpsi-rootuple_RDF.root");
  //TFile *fin = new TFile("MyPyRootTry/jpsi-rootuple_RDF_pyroot.root");
  TTree *_nt = (TTree*)fin->Get("dimuonTree");
  /*
  std::cout << "These are the columns:" << std::endl;
  for (auto branch : *_nt->GetListOfBranches()) {
    std::cout << "Branch: " << branch->GetName() << std::endl;
  }
  */
  
  Double_t Mmin = 2.915; 
  Double_t Mmax = 3.295;
  Double_t ctmin = -0.05; 
  Double_t ctmax = 0.25;

  char cut[198];
  int nn_c;
  //nn_c = sprintf(cut,"pt >= %5.1f && pt < %5.1f && TMath::Abs(y) >= %4.2f && TMath::Abs(y) < %4.2f && mu1pt>4.0 && mu2pt>4.0 && TMath::IsNaN(ct)!=1 && TMath::IsNaN(cterr)!=1",ptmin,ptmax,rapmin,rapmax);
  nn_c = sprintf(cut,"pt >= %5.1f && pt < %5.1f && TMath::Abs(y) >= %4.2f && TMath::Abs(y) < %4.2f && mu1pt>4.0 && mu2pt>4.0 && TMath::Abs(mu1y) < 1.2 && TMath::Abs(mu2y) < 1.2",ptmin,ptmax,rapmin,rapmax);
  std::cout << "INFO: cut -> " << cut << "  with size = " << nn_c << std::endl;

  //RooRealVar mass("mass","#mu^{+}#mu^{-} invariant mass [GeV]",Mmin,Mmax,"GeV");
  RooRealVar mass("mass","#mu^{+}#mu^{#font[122]{\55}} invariant mass [GeV]",Mmin,Mmax);
  RooRealVar ct("ct","Decay length [cm]",-0.05,0.35);
  RooRealVar cterr("cterr","cterr",0.0005,0.005,"cm");
  RooRealVar pt("pt","pt",ptmin,ptmax,"GeV");
  RooRealVar y("y","y",-rapmax,rapmax);
  RooRealVar mu1pt("mu1pt","mu1pt",4.0,200.0);
  RooRealVar mu2pt("mu2pt","mu2pt",4.0,200.0);
  RooRealVar mu1y("mu1y","mu1y",-1.2,1.2);
  RooRealVar mu2y("mu2y","mu2y",-1.2,1.2);

  RooDataSet *data = new RooDataSet("data","data",_nt,RooArgSet(mass,ct,cterr,pt,y,mu1pt,mu2pt,mu1y,mu2y),cut);
  Double_t dataentries = data->numEntries();
  //Double_t ptmean = data->mean(pt);
  std::cout << "INFO: Total entries in Tree " << _nt->GetEntries() << std::endl;
  std::cout << "INFO: Total entries in dataset " << dataentries << std::endl;
  std::cout << "dataset range mean " << data->mean(pt) << std::endl;
  //return;

  //**********************************
  //   Define Mass Model
  //**********************************
  
  Double_t n_val = 2.15; //1.95196;
  Double_t alpha_val = 1.8; // 1.54004;
  Double_t a1_val = -2.;
  Double_t m_val = 3.091;
  Double_t MASS_PEAK = 3.091;
  
  Double_t n_err = n_val*0.05; //1.32741e-01;
  Double_t alpha_err = alpha_val*0.05*5; //4.76148e-02;
  Double_t a1_err = a1_val*0.05*(-1.0);
  Double_t m_err = 1.0e-3;

  //if (rapmin>=1.2) m_val = 3.092;
  //if (rapmin>=1.5) m_val = 3.090;
  if (ptmin>=32.0) m_val = 3.090;
  if (ptmin>=42.0) m_val = 3.089;
  if (ptmin>=60.0) m_val = 3.087;

  RooRealVar m_mean("mean","m_mean",m_val,MASS_PEAK-.03,MASS_PEAK+.03);
  //RooRealVar m_mean("mean","m_mean",m_val);//,MASS_PEAK-.03,MASS_PEAK+.03);  
  RooRealVar m_sigma1("sig1"," Mass width 1",0.020,0.010,0.070);
  
  RooRealVar alpha("alpha","alpha",alpha_val,0.,5.);
  RooRealVar n("n","n",n_val);//,0.,100.);
  //RooRealVar n("n","n",n_val,0.,100.);
  RooRealVar a1("a1","A1",a1_val,-10.,0.);

  //RooGaussian a1_cons("a1_cons","a1cons",a1,RooConst(a1_val),RooConst(a1_err));
  RooRealVar a1fix("a1fix","a1fix",a1_val);
  RooRealVar a1fix_err("a1fix_err","a1fix_err",a1_err);
  RooGaussian a1_cons("a1_cons","a1cons",a1,a1fix,a1fix_err);
  RooRealVar nfix("nfix","nfix",n_val);
  RooRealVar nfix_err("nfix_err","nfix_err",n_err);
  RooGaussian n_cons("n_cons","ncons",n,nfix,nfix_err);
  RooRealVar alphafix("alphafix","alphafix",alpha_val);
  RooRealVar alphafix_err("alphafix_err","alphafix_err",alpha_err);
  RooGaussian alpha_cons("alpha_cons","alphacons",alpha,alphafix,alphafix_err);
  //RooGaussian m_cons("m_cons","mcons",m_mean,RooConst(m_val),RooConst(m_err));
  
  //Double_t fsig1_val = 0.55 - 0.55*(ptmin-10.)/70.-rapmin*0.05/0.3; 
  //if (fsig1_val < 0.12) fsig1_val = 0.0;
  //RooRealVar m_fraction_sigma1("fsig1","m_fraction_sigma1",fsig1_val,0.0,1.0); 
  RooRealVar m_fraction_sigma1("fsig1","m_fraction_sigma1",0.5,0.0,1.0); 

  Double_t dsig_val = 0.005;
  Double_t dsig_error = 0.0005;
  
  RooRealVar m_sigma2("dsig","m_d_sigma1",dsig_val,0.,.100);
  RooRealVar dsigfix("dsigfix","dsigfix",dsig_val);
  RooRealVar dsigfix_err("dsigfix_err","dsigfix_err",dsig_error);
  RooGaussian dsig_cons("dsig_cons","dsigcons",m_sigma2,dsigfix,dsigfix_err);
  
  RooFormulaVar m_d1("sig2", "sigma2", "@0+@1", RooArgList(m_sigma1,m_sigma2));
  RooCBShape m_cb1("m_cb1","m_cb1",mass, m_mean, m_sigma1, alpha, n);
  RooGaussian m_cb2("m_cb2","m_cb2",mass, m_mean,m_d1);
  //RooGaussian m_cb2("m_cb2","m_cb2",mass, m_mean,m_sigma2);
  //RooAddPdf pdf_m_signal("pdf_m_signal","signalshape",RooArgList(m_cb1,m_cb2),m_fraction_sigma1);
  RooAbsPdf *pdf_m_signal = 0;
  pdf_m_signal = new RooAddPdf("pdf_m_signal","pdf_m_signal",RooArgList(m_cb1,m_cb2),RooArgList(m_fraction_sigma1));

  //if (fsig1_val == 0.0 || ptmin>=40.) {
  if (ptmin>=40.) {
      m_fraction_sigma1.setVal(1.);
      m_fraction_sigma1.setConstant(kTRUE);
      m_sigma2.setVal(0.);
      m_sigma2.setConstant(kTRUE);
    }

  //-----------------------------------------------------------------
  // combinatorial background PDF (prompt or non-prompt J/psi)
  //RooExponential pdf_m_combinatorial("pdf_m_combinatorial","background mass PDF",mass,a1);   
  //RooChebychev pdf_m_combinatorial("pdf_m_combinatorial","Background for Mass",mass,RooArgList(a1));
  RooAbsPdf *pdf_m_combinatorial = 0;
  pdf_m_combinatorial = new RooExponential("pdf_m_combinatorial","background mass PDF",mass,a1);

  //-----------------------------------------------------------------
  // full model
  
  RooRealVar n_signal("n_sig","n_signal",dataentries*.9,dataentries*.5,dataentries*1.3);
  RooRealVar n_combinatorial("n_bck","n_combinatorial",dataentries*.1,dataentries*.01,dataentries*.5);

  RooRealVar n_prompt("n_p","n_signal_prompt",dataentries*.9*.5,dataentries*.1,dataentries);
  RooRealVar n_nonprompt("n_np","n_signal_nonprompt",dataentries*.9*.5,dataentries*.1,dataentries);
  
  //RooAddPdf m_model("m_model","m_model", RooArgList(pdf_m_signal, pdf_m_combinatorial), RooArgList(n_signal, n_combinatorial));
  RooAbsPdf * m_model = 0;
  //m_model = new RooAddPdf("m_model","m_model", RooArgList(*pdf_m_signal, *pdf_m_combinatorial), RooArgList(n_signal, n_combinatorial));
  m_model = new RooAddPdf("m_model","m_model", 
                    RooArgList(*pdf_m_signal,*pdf_m_signal, *pdf_m_combinatorial),
                    RooArgList(n_nonprompt,n_prompt, n_combinatorial));

  //-----------------------------------------------------------------
  //------------ tmp mass fit procedure -------------------
  
  RooAbsPdf * m_modeltmp = 0;
  m_modeltmp = new RooAddPdf("m_modeltmp","m_modeltmp",
			  RooArgList(*pdf_m_signal, *pdf_m_combinatorial),
			  RooArgList(n_signal, n_combinatorial));
  
  RooArgSet* extfitconstraintsmass = new RooArgSet(a1_cons,alpha_cons,dsig_cons);
  //RooFitResult* m_Fit = m_modeltmp->fitTo(*data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
  RooFitResult* m_Fit = m_modeltmp->fitTo(*data,Extended(),ExternalConstraints(*extfitconstraintsmass),Minos(kFALSE),Save(kTRUE), NumCPU(4));
  data->Print("v"); 
  m_Fit->Print("v");
  if ( (Int_t)m_Fit->covQual()!=3 ) {
    std::cout << "*** Fit fail *** " << std::endl;
  }
  /*
  RooPlot *xframe2 = mass.frame();
  data->plotOn(xframe2);
  m_modeltmp->plotOn(xframe2);
  xframe2->Draw();
  */
  m_mean.setConstant(kTRUE);
  a1.setConstant(kTRUE);
  //m_sigma1.setConstant(kTRUE);
  //m_fraction_sigma1.setConstant(kTRUE);
  //return;

  
  Double_t EtaL = rapmin*100.0;
  Double_t EtaH = rapmax*100.0;

 //-----------------------------------------------------------------
 //------------ ctau model -------------------

  // Signal prompt model
  RooRealVar res_sig_mean("res_m","res_m",0.0);
  Double_t res_sig_val = 0.95;
  Double_t res_sig_error = .02;
  RooRealVar res_sig_sigma("res_s","res_s",res_sig_val,0.3,2.0);
  RooGaussian res_sig_cons("res_sig_cons","res_sig_cons",res_sig_sigma,RooConst(res_sig_val),RooConst(res_sig_error));

  RooGaussModel res_signal1("res_signal1","res_signal1",ct,res_sig_mean,res_sig_sigma,cterr);
  Double_t res_d_val = 0.7;
  Double_t res_d_error = res_d_val*.05;
  RooRealVar res_sig_diff("res_d","res_s",res_d_val,.0,5.5);
  RooGaussian res_d_cons("res_d_cons","res_d_cons",res_sig_diff,RooConst(res_d_val),RooConst(res_d_error));
  RooFormulaVar res_sig_sigma2("res_sig_sigma2", "sigma2", "@0+@1", RooArgList(res_sig_sigma,res_sig_diff));
  RooGaussModel res_signal2("res_signal2","res_signal2",ct,res_sig_mean,res_sig_sigma2,cterr);
  //RooRealVar  res_frac_sigma1("res_f1","res_f1",0.88);//,0.,1.);
  RooRealVar  res_frac_sigma1("res_f1","res_f1",0.88,0.,1.);
  //RooAddModel res_signal("res_signal",  "Resolution ",RooArgList(res_signal1,res_signal2),RooArgList(res_frac_sigma1));
  RooAddModel *res_signal = 0;
  res_signal = new RooAddModel("res_signal",  "Resolution ",RooArgList(res_signal1,res_signal2),RooArgList(res_frac_sigma1));


  if(ptmin>= 40) {
    res_frac_sigma1.setVal(1.);
    res_frac_sigma1.setConstant(kTRUE);
    res_sig_diff.setVal(0.);
    res_sig_diff.setConstant(kTRUE);
  }

  // Signal nonprompt model
  RooRealVar ctau("#lambda","ct",0.035,0.010,0.070);
  //RooDecay pdf_t_nonprompt("pdf_t_nonprompt","pdf_t_nonprompt",ct,ctau,*res_signal,RooDecay::SingleSided);
  RooAbsPdf  *pdf_t_nonprompt = 0;
  pdf_t_nonprompt = new RooDecay("pdf_t_nonprompt","pdf_t_nonprompt",ct,ctau,*res_signal,RooDecay::SingleSided);

  RooRealVar nonprompt_fraction("npf","npf",0.7,0.0,1.0);
  //RooAddPdf pdf_t_signal("pdf_t_signal","pdf_t_signal",RooArgList(*pdf_t_nonprompt,*res_signal),RooArgList(n_nonprompt,n_prompt));
  RooAbsPdf *pdf_t_signal = 0;
  pdf_t_signal = new RooAddPdf("pdf_t_signal","pdf_t_signal",RooArgList(*pdf_t_nonprompt,*res_signal),RooArgList(n_nonprompt,n_prompt));

  // Bkg model
  RooRealVar res_bck_mean("res_bck_mean","res_bck_mean",0.0);
  
  Double_t bck_s1_val = 1.7;
  //if (ptmin>15) bck_s1_val = 1.6;
  //if (ptmin>20) bck_s1_val = 1.5;
  if (ptmin>30) bck_s1_val = 1.4;
  if (ptmin>40) bck_s1_val = 1.3;
  if (ptmin>60) bck_s1_val = 1.2;
  if (ptmin>80) bck_s1_val = 1.;
  
  Double_t bck_s1_error = 0.05;
  RooRealVar res_bck_sigma1("bck_s1","bck_s1",bck_s1_val,0.1,5.);
  RooGaussian bck_s1_cons("bck_s1_cons","bck_ct_cons",res_bck_sigma1,RooConst(bck_s1_val),RooConst(bck_s1_error));
  
  Double_t bck_d_val = 4.5;
  Double_t bck_d_error = 0.2;
  RooRealVar res_bck_sigma2("bck_d","bck_d",bck_d_val,0.,40.0);
  RooGaussian bck_d_cons("bck_d_cons","bck_d_cons",res_bck_sigma2,RooConst(bck_d_val),RooConst(bck_d_error));
  
  RooFormulaVar res_bck_dsig("bck_s2", "bck_s2", "@0+@1", RooArgList(res_bck_sigma1,res_bck_sigma2));
  
  RooGaussModel res_back1("res_back1","res_back1",ct,res_sig_mean,res_bck_sigma1,cterr);
  RooGaussModel res_back2("res_back2","res_back2",ct,res_sig_mean,res_bck_dsig,cterr);

  Double_t bck_fs1_val = 0.8;
  Double_t bck_fs1_error = 0.03;
  RooRealVar res_b_frac_sigma1("bck_fs1","bck_fs1",bck_fs1_val);//,0.0,1.);
  RooGaussian bck_fs1_cons("bck_fs1_cons","bck_fs1_cons",res_b_frac_sigma1,RooConst(bck_fs1_val),RooConst(bck_fs1_error));
  
  RooAddModel res_back("res_back",  "Resolution b",RooArgList(res_back1,res_back2),RooArgList(res_b_frac_sigma1));

  Double_t bck_ct_val = 0.039;
  Double_t bck_ct_error = 0.001;
  RooRealVar ctau_b_nonprompt("bck_ct","bck_ctnp",bck_ct_val);//, 0.010, .070);
  RooGaussian bck_ct_cons("bck_ct_cons","bck_ct_cons",ctau_b_nonprompt,RooConst(bck_ct_val),RooConst(bck_ct_error));
  //use the same resolution function as the signal for now
  //RooDecay pdf_t_b_nonprompt("pdf_t_b_nonprompt","pdf_t_b_nonprompt",ct,ctau_b_nonprompt,*res_signal,RooDecay::SingleSided);
  RooAbsPdf  *pdf_t_b_nonprompt = 0;
  pdf_t_b_nonprompt = new RooDecay("pdf_t_b_nonprompt","pdf_t_b_nonprompt",ct,ctau_b_nonprompt,*res_signal,RooDecay::SingleSided);

  Double_t bck_npf_val = 0.78;
  if (ptmin < 30) bck_npf_val = 0.74;
  if (ptmin < 25) bck_npf_val = 0.70;
  if (ptmin < 20) bck_npf_val = 0.66;
  if (ptmin < 15) bck_npf_val = 0.64;
  Double_t bck_npf_val_error = 0.03;
  RooRealVar nonprompt_b_fraction("bck_npf","bck_fnp",bck_npf_val,0.0,1.0);
  RooGaussian bck_npf_cons("bck_npf_cons","bck_ct_cons",nonprompt_b_fraction,RooConst(bck_npf_val),RooConst(bck_npf_val_error));
  
  //RooAddPdf pdf_t_combinatorial("pdf_t_combinatorial","pdf_t_combinatorial0", RooArgList(*pdf_t_b_nonprompt,res_back),RooArgList(nonprompt_b_fraction));
  RooAbsPdf  *pdf_t_combinatorial = 0;
  pdf_t_combinatorial= new RooAddPdf("pdf_t_combinatorial","pdf_t_combinatorial0", RooArgList(*pdf_t_b_nonprompt,res_back),RooArgList(nonprompt_b_fraction));

  //Total ct Model

  //RooAddPdf t_model("t_model","Total time Model",
  //              RooArgList(pdf_t_signal, pdf_t_combinatorial),
  //               RooArgList(n_signal, n_combinatorial));
  RooAddPdf t_model("t_model","Total time Model",RooArgList(*pdf_t_nonprompt,*res_signal, *pdf_t_combinatorial),
		    RooArgList(n_nonprompt,n_prompt, n_combinatorial));


  
  // Total 2D Model
  
  RooProdPdf pdf_mxt_combinatorial("pdf_mxt_combinatorial","pdf_mt_combinatorial",RooArgSet(*pdf_m_combinatorial,*pdf_t_combinatorial));
  //RooProdPdf pdf_mxt_signal("pdf_mxt_signal","pdf_mt_signal",RooArgSet(pdf_m_signal,pdf_t_signal));
  RooProdPdf pdf_mxt_signal_np("pdf_mxt_signal_np","pdf_mt_signal_np",RooArgSet(*pdf_m_signal,*pdf_t_nonprompt));
  RooProdPdf pdf_mxt_signal_p("pdf_mxt_signal_p","pdf_mt_signal_p",RooArgSet(*pdf_m_signal,*res_signal));
  
  RooAddPdf mxt_model_u("mxt_model_u","Total mass,time Model",
			RooArgList(pdf_mxt_signal_np,pdf_mxt_signal_p, pdf_mxt_combinatorial),
			RooArgList(n_nonprompt,n_prompt, n_combinatorial));
  
  n_signal.setConstant(kFALSE);
  n_combinatorial.setConstant(kFALSE);
  
  /*   
  RooProdPdf mxt_model("mxt_model","model with constraint",
		       RooArgSet(mxt_model_u,
				 a1_cons,
				 //... n free and alpha cte. n_cons,alpha_cons,
				 alpha_cons,
				 dsig_cons,
				 //bck_ct_cons,
				 res_sig_cons,
				 res_d_cons,
				 bck_npf_cons,
				 bck_s1_cons,
				 //bck_fs1_cons,
				 bck_d_cons
				 ));
  */

  RooArgSet* extfitconstraints = new RooArgSet(alpha_cons,dsig_cons,res_sig_cons,res_d_cons,bck_npf_cons,bck_s1_cons,bck_d_cons);
  
  //RooFitResult* t_Fit = mxt_model.fitTo(*data,Extended(),Minos(kFALSE),NumCPU(4),Offset(kTRUE),Save(),ConditionalObservables(RooArgSet(cterr)));
  RooFitResult* t_Fit = mxt_model_u.fitTo(*data,Extended(),ExternalConstraints(*extfitconstraints),Minos(kFALSE),NumCPU(4),Offset(kTRUE),Save(),ConditionalObservables(RooArgSet(cterr)));
  //RooFitResult* t_Fit = mxt_model_u.fitTo(*data,Minos(kFALSE),NumCPU(8),Offset(kTRUE),Save(),ConditionalObservables(RooArgSet(cterr)));
  data->Print("v"); 
  t_Fit->Print("v");
  if ( (Int_t)t_Fit->covQual()!=3 ) {
    std::cout << "*** Fit fail *** " << std::endl;
  }

  //-----------------------------------------------------------------
  //------------ Mass plot -------------------
  m_mean.setConstant(kTRUE);
  n_signal.setConstant(kFALSE);
  n_combinatorial.setConstant(kFALSE);
  
  TCanvas* canv_nominal_pull= CreateCanvasNomPull("canv_nominal_pull", t_Fit, data, mass, Mmax, Mmin, m_model, pdf_m_signal, pdf_m_combinatorial, n_nonprompt, n_prompt, n_combinatorial, m_sigma1, m_sigma2, m_fraction_sigma1, m_mean, ptmin, ptmax, rapmin, rapmax);  
  canv_nominal_pull->Print(Form("plots_ptandeta2D/mass2D_jpsiFit_ptetabins_%1.0f_%1.0f_%1.0f_%1.0f_Pull.png",ptmin, ptmax, EtaL, EtaH));
  canv_nominal_pull->Print(Form("plots_ptandeta2D/mass2D_jpsiFit_ptetabins_%1.0f_%1.0f_%1.0f_%1.0f_Pull.pdf",ptmin, ptmax, EtaL, EtaH));

  //-----------------------------------------------------------------
  //------------ Lifetime plot -------------------
  m_mean.setConstant(kTRUE);
  m_sigma1.setConstant(kTRUE);
  m_fraction_sigma1.setConstant(kTRUE);
  m_sigma2.setConstant(kTRUE);

  TCanvas* canvct_nominal_pull= CreateCanvasctPull("canvct_nominal_pull", t_Fit, data, ct, cterr, ctmax, ctmin, &t_model, pdf_t_nonprompt, res_signal, pdf_t_combinatorial, n_signal, n_nonprompt, n_prompt, n_combinatorial, ptmin, ptmax, rapmin, rapmax);  
  canvct_nominal_pull->Print(Form("plots_ptandeta2D/ct2D_jpsiFit_ptetabins_%1.0f_%1.0f_%1.0f_%1.0f_Pull.png",ptmin, ptmax, EtaL, EtaH));
  canvct_nominal_pull->Print(Form("plots_ptandeta2D/ct2D_jpsiFit_ptetabins_%1.0f_%1.0f_%1.0f_%1.0f_Pull.pdf",ptmin, ptmax, EtaL, EtaH));

  //https://root-forum.cern.ch/t/how-to-close-or-delete-a-canvas-inside-a-macro/23619
  if (canvct_nominal_pull) {
    canvct_nominal_pull->Close();
    gSystem->ProcessEvents();
    delete canvct_nominal_pull;
    canvct_nominal_pull = 0;
    std::cout<<"it has been close\n";
  }
  if(canvct_nominal_pull)std::cout<<"Didn't close\n";

 //parameters yieldas output file
 ofstream salida_nominal(Form("txtfiles_ptandeta2D/output_jpsiFit_nominal_%1.0f_%1.0f_%1.0f_%1.0f.txt",ptmin, ptmax, EtaL, EtaH));
 salida_nominal.is_open();
 save_result(salida_nominal, n_prompt, n_combinatorial);
  
 std:cout << std::endl << "Done with plotting" << std::endl << std::endl;
  
}
  

