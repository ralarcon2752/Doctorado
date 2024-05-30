#ifndef MYUTILITIES_H
#define MYUTILITIES_H

#include <RooPlot.h>
#include "RooFit.h"
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

using namespace RooFit;
using namespace std;

TCanvas* CreateCanvasNomPull(TString cname, RooFitResult* result, RooDataSet *data, RooRealVar M, Double_t supM, Double_t infM,  RooAbsPdf* MassModel, RooAbsPdf* sumgau, RooAbsPdf* bkg1, RooRealVar Nsnp, RooRealVar Nsp, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth, Double_t etal, Double_t etah) 
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
  //Double_t probChi2 = TMath::Prob(chi2tmp, ndof)*100.0;
  Double_t probChi2 = TMath::Prob(chi2tmp, ndof);
  cout<<" Ndof : "<<ndof<< endl;
  cout<<" Chi2 : "<<chi2tmp<< endl;
  cout<<" Chi2/ndof : "<<Mframe->chiSquare()<<" Chi2Prob : "<< probChi2 <<  endl;


  //MassModel->plotOn(Mframe,Components(*sumgau),LineColor(kRed),LineWidth(2),Name("Signal"));
  MassModel->plotOn(Mframe,Components(*bkg1),LineColor(kGreen),LineWidth(2),LineStyle(kDashed),Name("bkg")); 
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
  //legpar->AddEntry("",Form("M(J/#psi) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  //legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{J/#psi} np = %1.0f #pm %1.0f",Nsnp.getVal(),Nsnp.getError()),"");
  legpar->AddEntry("",Form("N_{J/#psi} p = %1.0f #pm %1.0f",Nsp.getVal(),Nsp.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
  legpar->Draw();
  
  TLegend *legMass = new TLegend(0.6,0.23,0.8,0.45);
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
  /*
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
  */
  
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

TCanvas* CreateCanvasctPull(TString cname, RooFitResult* result, RooDataSet *data, RooRealVar myct, RooRealVar mycterr, Double_t supct, Double_t infct,  RooAbsPdf* ctmodel, RooAbsPdf* npmod, RooAbsPdf* pmod, RooAbsPdf* bkgmod, RooRealVar Nst, RooRealVar Nsnp, RooRealVar Nsp, RooRealVar Nb, Double_t ptl, Double_t pth, Double_t etal, Double_t etah){

  int H = 600;
  int W = 800;
  TCanvas* canvt = new TCanvas("canvt","canvt",50,50,W,H);
  canvt->cd() ;
  canvt->SetLeftMargin(0.005);
  canvt->SetRightMargin(0.01);
  canvt->SetTopMargin(0.09);
  canvt->SetBottomMargin(0.1);
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
  gPad->SetLogy();

  Int_t ctbins = 600;//150 cada 20 mum
  RooPlot* frame_t = myct.frame(Range(infct,supct),Bins(ctbins));
  data->plotOn(frame_t,DataError(RooAbsData::SumW2),MarkerSize(0.7));
  ctmodel->plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),LineWidth(3),Name("fittotalct"));
  //ctmodel->plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),Normalization(data->numEntries(),RooAbsReal::NumEvent),LineWidth(3),Name("fittotalct"));
  Double_t nfp = result->floatParsFinal().getSize();
  Double_t ndof = ctbins - nfp;
  Float_t chi2 = frame_t->chiSquare();
  Float_t chi2_n = frame_t->chiSquare(nfp);
  cout << endl << "ndof = " << ndof << endl;
  cout << endl << "Chi2()/ndof = " << chi2 << endl;
  cout << "Chi2() = " << chi2*(ctbins-1-nfp)<< ", ndof = "  << ctbins-1-nfp << ", Chi2Prob = " << TMath::Prob(chi2*(ctbins-1-nfp),ctbins-1-nfp)<< endl;
  cout << "Chi2(" << nfp << ") = " << chi2_n << endl;
  RooHist* pull_t = frame_t->pullHist();

  Double_t n_tot = Nsp.getVal()+Nsnp.getVal();
  Double_t npf = Nsnp.getVal()/n_tot;
  Double_t n_tot_error = TMath::Sqrt(Nsp.getError()*Nsp.getError() + Nsnp.getError()*Nsnp.getError());
  Double_t npf_error = npf * TMath::Sqrt( pow(Nsnp.getError()/Nsnp.getVal(),2) + pow(n_tot_error/n_tot,2) );
  Double_t n_back =Nb.getVal();
  
  //Double_t np= Nst.getVal()*(1. - nonprompt_fraction.getVal());
  //Double_t nnp= Nst.getVal()*nonprompt_fraction.getVal();
  //ctmodel.plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),Components(bkgmod),LineColor(kGreen),LineWidth(3),LineStyle(kDashed));
  //npmod.plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),LineColor(kBlue),LineWidth(3),Normalization(nnp, RooAbsReal::NumEvent));
  //pmod.plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),LineColor(kBlack),LineWidth(3),Normalization(np, RooAbsReal::NumEvent));
  
  Double_t np= Nsp.getVal();                                                    
  Double_t nnp= Nsnp.getVal();                                                       
  //ctmodel.plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),Components(bkgmod),LineColor(9),LineWidth(3),LineStyle(kDashed));      
  bkgmod->plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),LineColor(kGreen),LineWidth(3),LineStyle(kDashed),Normalization(n_back, RooAbsReal::NumEvent),Name("bkgct"));
  npmod->plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),LineColor(kRed),LineWidth(3),Normalization(nnp, RooAbsReal::NumEvent),Name("signalnp"));
  pmod->plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),LineColor(kBlack),LineWidth(3),Normalization(np, RooAbsReal::NumEvent),Name("signalp"));
  data->plotOn(frame_t,DataError(RooAbsData::SumW2),MarkerSize(0.7),XErrorSize(0),Name("Datact"));
  ctmodel->plotOn(frame_t,ProjWData(mycterr,*data,kTRUE),LineWidth(3));			     
  
  frame_t->SetTitle("");
  frame_t->SetMinimum(0.2);
  Double_t ppdl_max_ = Nst.getVal() * 2.;
  frame_t->SetMaximum(ppdl_max_);
  
  frame_t->GetXaxis()->SetTitle("Decay length [cm]");
  frame_t->GetYaxis()->SetTitle("Events / 5 #mum");
  frame_t->GetXaxis()->SetLabelFont(42);
  frame_t->GetXaxis()->SetTitleFont(42);
  frame_t->GetYaxis()->SetLabelFont(42);
  frame_t->GetYaxis()->SetTitleFont(42);
  frame_t->GetXaxis()->SetLabelOffset(0.004);
  frame_t->GetXaxis()->SetTitleOffset(0.9);
  frame_t->GetXaxis()->SetTitleSize(0.055);
  frame_t->GetXaxis()->SetLabelSize(0.050); //0.045);
  frame_t->GetYaxis()->SetLabelOffset(0.004);
  frame_t->GetYaxis()->SetTitleOffset(0.6);
  frame_t->GetYaxis()->SetTitleSize(0.08);
  frame_t->GetYaxis()->SetLabelSize(0.07); //0.045);
  frame_t->Draw();

  TLegend *leg = new TLegend(0.7,0.45,0.9,0.88); 
  leg->SetTextSize(0.06);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(frame_t->findObject("Datact")," Data","ep"); 
  leg->AddEntry(frame_t->findObject("fittotalct")," Fit result","l");
  leg->AddEntry(frame_t->findObject("signalp")," Prompt signal","l");
  leg->AddEntry(frame_t->findObject("signalnp")," Nonprompt signal","l");			     
  leg->AddEntry(frame_t->findObject("bkgct"),"Comb. backg.","l");
  leg->Draw();			     

  TLegend *legct = new TLegend(0.35,0.66,0.55,0.88);
  legct->SetTextFont(42);
  legct->SetTextSize(0.06);
  legct->SetFillColor(0);
  legct->SetBorderSize(0);
  legct->SetFillStyle(0);
  legct->AddEntry("",Form("%1.1f < p_{T}^{#mu^{+}#mu^{-}} < %1.1f GeV ",ptl,pth),"");
  legct->AddEntry("",Form("%1.1f < |y^{#mu^{+}#mu^{-}}| < %1.1f",etal,etah),"");
  legct->Draw();


  pad2->cd();

  // Create a new frame to draw the pull distribution
  //pull_t->SetMarkerSize(.7);
  RooPlot* framem2_t = myct.frame(Range(infct,supct),Bins(ctbins)) ;
  //framem2_t->addPlotable(pull_t,"P") ;
  framem2_t->addPlotable(pull_t,"XP") ;// setting Y error to cero
  framem2_t->SetYTitle(" (Data-Fit)/#sigma");
  framem2_t->SetLabelSize(0.1,"XY");
  framem2_t->SetTitleSize(0.13,"X");
  framem2_t->SetTitleSize(0.11,"Y");
  framem2_t->GetYaxis()->CenterTitle();
  framem2_t->GetXaxis()->CenterTitle();
  framem2_t->GetYaxis()->SetNdivisions(505,1);
  //framem2_t->GetXaxis()->SetNdivisions(305,1);
  framem2_t->GetXaxis()->SetTickLength(0.07);
  framem2_t->SetTitleOffset(0.9,"X");
  framem2_t->SetTitleOffset(0.4,"Y");
  framem2_t->SetMaximum(4.9);
  framem2_t->SetMinimum(-4.9);
  framem2_t->Draw();

  canvt->cd();
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

  canvt->Update();
  //data->Print("v"); 
  //result->Print("v");
  return canvt; 

}



#endif
