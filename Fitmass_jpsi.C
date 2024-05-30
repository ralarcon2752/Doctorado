#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

#include "RooFit.h"
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

#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooNumIntConfig.h"

using namespace RooFit;
using namespace std;

void save_result(ofstream& salida, RooRealVar Ns, RooRealVar Nb)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError();
  cout << " el archivo se escribio bien" << endl; 
  return;
}

TCanvas* CreateCanvasNomPull(TString cname, RooFitResult* result, RooDataSet *data, RooRealVar M, Double_t supM, Double_t infM,  RooAbsPdf* MassModel, RooAbsPdf* sumgau, RooAbsPdf* bkg1, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth, Double_t etal, Double_t etah) 
{
  //Double_t nbin = ((supM-infM)/0.010);
  //Double_t nbin = ((supM-infM)/0.001)+1;
  Double_t nbin = ((supM-infM)/0.004)+1;

  int H = 600;
  int W = 800;
  
  TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
  //TCanvas *c1 = new TCanvas(cname,cname,W,H);
  //c1->Divide(1,2);
  //c1->cd(1) ;
  c1->cd() ;  
  c1->SetLeftMargin(0.005);
  c1->SetRightMargin(0.01);
  c1->SetTopMargin(0.09);
  c1->SetBottomMargin(0.1);
  //gPad->SetLogy();


  TPad *pad1 = new TPad("pad1", "padi",0.01,0.411,0.9903769, 0.99 );
  pad1->SetLeftMargin(0.09);   
  pad1->SetRightMargin(0.019);
  pad1->SetTopMargin(0.09);
  pad1->SetBottomMargin(0.0);
 

  TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.41);
  pad2->SetLeftMargin(0.09);
  pad2->SetRightMargin(0.019);  
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.25);
  //pad2->SetTickx(0);
  pad2->SetFillColor(0);
  pad2->SetGridx(0);
  pad2->SetGridy(0);

  pad1->Draw();
  pad2->Draw();
  pad1->cd(); 

  RooPlot* Mframe = M.frame(infM,supM,nbin);
  data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0));
  //MassModel->plotOn(Mframe);
  MassModel->plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
  
  RooHist* hpullm2 = Mframe->pullHist() ;
  
  Double_t nfloatpars = result->floatParsFinal().getSize();
  Double_t ndof = nbin - nfloatpars;
  Double_t chi2tmp = Mframe->chiSquare()*ndof;
  Double_t probChi2 = TMath::Prob(chi2tmp, ndof);
  
  cout<<" Ndof : "<<ndof<< endl;
  cout<<" Chi2 : "<<chi2tmp<< endl;
  cout<<" Chi2/ndof : "<<Mframe->chiSquare()<<" Chi2Prob : "<< probChi2 <<  endl;


  //MassModel->plotOn(Mframe,Components(*sumgau),LineColor(kRed),LineWidth(2),Name("Signal"));
  MassModel->plotOn(Mframe,Components(*bkg1),LineColor(kGreen),LineStyle(kDashed),LineWidth(2),Name("bkg")); 
  data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
  MassModel->plotOn(Mframe);
 
  //Mframe->SetYTitle("Events / 1.0 MeV");
  Mframe->SetYTitle("Events / 4.0 MeV");  
  Mframe->SetLabelSize(0.07,"XY");
  Mframe->SetTitleSize(0.08,"XY");
  Mframe->GetYaxis()->CenterTitle();   
  Mframe->GetXaxis()->CenterTitle();
  Mframe->GetYaxis()->SetNdivisions(505,1);
  Mframe->GetXaxis()->SetNdivisions(505,1);
  Mframe->GetXaxis()->SetTickLength(0.06);    
  Mframe->GetXaxis()->SetDecimals(1); 
  Mframe->SetTitleOffset(0.8,"X");
  Mframe->SetTitleOffset(0.6,"Y");
  Mframe->SetMinimum(1.0); 
  Mframe->Draw();
  
  //TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
  TLegend *leg = new TLegend(0.18,0.50,0.38,0.88); 
  leg->SetTextSize(0.06);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
  leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
  //leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
  leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l");
  leg->Draw();

  Double_t Mpsi = mean.getVal()*1000.0;
  Double_t MpsiE = mean.getError()*1000.0;
 
  TLegend *legpar = new TLegend(0.6,0.58,0.8,0.88);
  legpar->SetTextSize(0.06);
  legpar->SetTextFont(42);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("M(J/#psi) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  //legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{J/#psi} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
  legpar->Draw();
  
  TLegend *legMass = new TLegend(0.6,0.3,0.8,0.5);
  legMass->SetTextFont(42); 
  legMass->SetTextSize(0.06);  
  legMass->SetFillColor(0); 
  legMass->SetBorderSize(0);
  legMass->SetFillStyle(0); 
  //legMass->SetHeader(Form(" %1.1f #leq p_{T}^{#mu^{+}#mu^{-}} < %1.1f GeV ",ptl,pth));
  //legMass->SetHeader(Form(" %1.2f #leq |y^{#mu^{+}#mu^{-}}| < %1.2f",etal,etah));
  legMass->AddEntry("",Form("%1.1f < p_{T}^{#mu^{+}#mu^{-}} < %1.1f GeV ",ptl,pth),"");
  legMass->AddEntry("",Form("%1.1f < |y^{#mu^{+}#mu^{-}}| < %1.1f",etal,etah),"");
  legMass->Draw();
 
  pad2->cd();
  
  // Create a new frame to draw the pull distribution 
  RooPlot* framem2 = M.frame(infM,supM,nbin) ;
  //framem2->addPlotable(hpullm2,"P") ;
  framem2->addPlotable(hpullm2,"XP") ;// setting Y error to cero
  
  framem2->SetYTitle(" (Data-Fit)/#sigma");
  framem2->SetLabelSize(0.1,"XY");
  framem2->SetTitleSize(0.13,"X");
  framem2->SetTitleSize(0.11,"Y");  
  framem2->GetYaxis()->CenterTitle();   
  framem2->GetXaxis()->CenterTitle();
  framem2->GetYaxis()->SetNdivisions(505,1);
  framem2->GetXaxis()->SetNdivisions(505,1);
  framem2->GetXaxis()->SetTickLength(0.07);  
  framem2->SetTitleOffset(0.9,"X");
  framem2->SetTitleOffset(0.4,"Y");
  framem2->SetMaximum(4.9);
  framem2->SetMinimum(-4.9);
  framem2->Draw();
 
  TLine *line1 = new TLine(infM,0.0,supM,0.0);
  line1->SetLineColor(1);
  line1->SetLineWidth(1);
  line1->Draw();
  TLine *line2 = new TLine(infM,4.0,supM,4.0);
  line2->SetLineColor(2);
  line2->SetLineWidth(2);
  line2->SetLineStyle(2);
  line2->Draw();
  TLine *line3 = new TLine(infM,-4.0,supM,-4.0);
  line3->SetLineColor(2);
  line3->SetLineWidth(2);
  line3->SetLineStyle(2);
  line3->Draw();
  
  c1->cd();
  
  TLatex *   tex1 = new TLatex(0.97,0.95,"L = zzz fb^{-1} (#sqrt{s} = 13.6 TeV)");
  tex1->SetNDC();
  tex1->SetTextAlign(31);
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.04); 
  tex1->SetLineWidth(2);
  
  TLatex *tex2 = new TLatex(0.11,0.95,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.04); 
  tex2->SetLineWidth(2);
  
  TLatex *tex3 = new TLatex(0.19,0.95,"Preliminary");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.04); 
  tex3->SetLineWidth(2);
  
  tex1->Draw();  
  tex2->Draw();
  tex3->Draw();
  
  c1->Modified();
  /*
  gPad->Update();
  gPad->RedrawAxis();
  TLine l;
  l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
  */
  return c1; 
  
}



void Fitmass_jpsi(float ptmin=25.0, float ptmax=26.0, float rapmin=0.0, float rapmax=1.2){
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);
  
  //TChain *_nt = new TChain("dimuonTree","");
  //_nt->Add("jpsi-rootuple.root");
  TFile *fin = new TFile("jpsi-rootuple.root");
  TTree *_nt = (TTree*)fin->Get("dimuonTree");
  /*
  std::cout << "These are the columns:" << std::endl;
  for (auto branch : *_nt->GetListOfBranches()) {
    std::cout << "Branch: " << branch->GetName() << std::endl;
  }
  return;
  */
  
  Double_t Mmin = 2.915; 
  Double_t Mmax = 3.295;

  char cut[128];
  int nn_c;
  nn_c = sprintf(cut,"pt >= %5.1f && pt < %5.1f && TMath::Abs(y) >= %4.2f && TMath::Abs(y) < %4.2f && mu1pt>4.0 && mu2pt>4.0",ptmin,ptmax,rapmin,rapmax);
  //nn_c = sprintf(cut,"pt >= %5.1f && pt < %5.1f && TMath::Abs(y) >= %4.2f && TMath::Abs(y) < %4.2f && mu1pt>4.0 && mu2pt>4.0 && TMath::IsNaN(ct)!=1 && TMath::IsNaN(cterr)!=1",ptmin,ptmax,rapmin,rapmax);
  std::cout << "INFO: cut -> " << cut << "  with size = " << nn_c << std::endl;

  //RooRealVar mass("mass","#mu^{+}#mu^{-} invariant mass [GeV]",Mmin,Mmax,"GeV");
  RooRealVar mass("mass","#mu^{+}#mu^{-} invariant mass [GeV]",Mmin,Mmax);
  RooRealVar ct("ct","ct",-0.05,0.35,"cm");
  RooRealVar cterr("cterr","cterr",0.0005,0.005,"cm");
  RooRealVar pt("pt","pt",ptmin,ptmax,"GeV");
  RooRealVar y("y","y",-rapmax,rapmax);
  RooRealVar mu1pt("mu1pt","mu1pt",4.0,200.0);
  RooRealVar mu2pt("mu2pt","mu2pt",4.0,200.0);

  RooDataSet *data = new RooDataSet("data","data",_nt,RooArgSet(mass,ct,cterr,pt,y,mu1pt,mu2pt),cut);
  Double_t dataentries = data->numEntries();
  //Double_t ptmean = data->mean(pt);
  std::cout << "INFO: Total entries in Tree " << _nt->GetEntries() << std::endl;
  std::cout << "INFO: Total entries in dataset " << dataentries << std::endl;
  std::cout << "dataet range mean " << data->mean(pt) << std::endl;

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

  RooRealVar m_mean("mean","m_mean",m_val,MASS_PEAK-.03,MASS_PEAK+.03);
  //RooRealVar m_mean("mean","m_mean",m_val);//,MASS_PEAK-.03,MASS_PEAK+.03);  
  RooRealVar m_sigma1("sig1"," Mass width 1",0.020,0.010,0.070);
  
  RooRealVar alpha("alpha","alpha",alpha_val,0.,5.);
  RooRealVar n("n","n",n_val);//,0.,100.);
  //RooRealVar n("n","n",n_val,0.,100.);
  RooRealVar a1("a1","A1",a1_val,-10.,0.);

  RooGaussian a1_cons("a1_cons","a1cons",a1,RooConst(a1_val),RooConst(a1_err));
  RooGaussian alpha_cons("alpha_cons","alphacons",alpha,RooConst(alpha_val),RooConst(alpha_err));
  //RooGaussian n_cons("n_cons","ncons",n,RooConst(n_val),RooConst(n_err));
  //RooGaussian m_cons("m_cons","mcons",m_mean,RooConst(m_val),RooConst(m_err));

  //Double_t fsig1_val = 0.55 - 0.55*(ptmin-10.)/70.-rapmin*0.05/0.3; 
  //if (fsig1_val < 0.12) fsig1_val = 0.0;
  //RooRealVar m_fraction_sigma1("fsig1","m_fraction_sigma1",fsig1_val,0.0,1.0); 
  RooRealVar m_fraction_sigma1("fsig1","m_fraction_sigma1",0.5,0.0,1.0); 

  Double_t dsig_val = 0.005;
  Double_t dsig_error = 0.0005;
  
  RooRealVar m_sigma2("dsig","m_d_sigma1",dsig_val,0.,.100);
  RooGaussian dsig_cons("dsig_cons","dsigcons",m_sigma2,RooConst(dsig_val),RooConst(dsig_error));
  
  RooFormulaVar m_d1("sig2", "sigma2", "@0+@1", RooArgList(m_sigma1,m_sigma2));
  RooCBShape m_cb1("m_cb1","m_cb1",mass, m_mean,m_sigma1, alpha,n);
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
  //RooAddPdf m_model("m_model","m_model", RooArgList(pdf_m_signal, pdf_m_combinatorial), RooArgList(n_signal, n_combinatorial));
  RooAbsPdf * m_model = 0;
  m_model = new RooAddPdf("m_model","m_model",
			  RooArgList(*pdf_m_signal, *pdf_m_combinatorial),
			  RooArgList(n_signal, n_combinatorial));

  //-----------------------------------------------------------------
  //------------ Mass fit procedure -------------------
  RooArgSet* extfitconstraintsmass = new RooArgSet(a1_cons,alpha_cons,dsig_cons);
  
  //RooFitResult* m_Fit = m_model->fitTo(*data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
  RooFitResult* m_Fit = m_model->fitTo(*data,Extended(),ExternalConstraints(*extfitconstraintsmass),Minos(kFALSE),Save(kTRUE), NumCPU(4));
  data->Print("v"); 
  m_Fit->Print("v");
  if ( (Int_t)m_Fit->covQual()!=3 ) {
    std::cout << "*** Fit fail *** " << std::endl;
  }

  Double_t EtaL = rapmin*100.0;
  Double_t EtaH = rapmax*100.0; 
  TCanvas* canv_nominal_pull= CreateCanvasNomPull("canv_nominal_pull", m_Fit, data, mass, Mmax, Mmin, m_model, pdf_m_signal, pdf_m_combinatorial, n_signal, n_combinatorial, m_sigma1, m_sigma2, m_fraction_sigma1, m_mean, ptmin, ptmax, rapmin, rapmax);  
 canv_nominal_pull->Print(Form("plots_ptandeta/mass_jpsiFit_ptetabins_%1.0f_%1.0f_%1.0f_%1.0f_Pull.png",ptmin, ptmax, EtaL, EtaH));
 canv_nominal_pull->Print(Form("plots_ptandeta/mass_jpsiFit_ptetabins_%1.0f_%1.0f_%1.0f_%1.0f_Pull.pdf",ptmin, ptmax, EtaL, EtaH));

 //parameters yieldas output file
 ofstream salida_nominal(Form("txtfiles_ptandeta/output_jpsiFit_nominal_%1.0f_%1.0f_%1.0f_%1.0f.txt",ptmin, ptmax, EtaL, EtaH));
 salida_nominal.is_open();
 save_result(salida_nominal, n_signal, n_combinatorial);

}
  
