#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TKey.h"
#include <iostream>
#include <TStyle.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

//#include "Macro.h"
void setTDRStyle();


void Opt_isocut(){
 
  
  setTDRStyle();
  //gStyle->SetPalette(1);
  
  std::vector<TString> masses;
  masses.push_back("40");
  masses.push_back("50");
  masses.push_back("60");
  masses.push_back("100");
  masses.push_back("200");
  masses.push_back("500");
  masses.push_back("1100");
  masses.push_back("1500");

  TCanvas* c1 = new TCanvas("Plot", "Plot", 800, 600);
      
  TLegend *legend = new TLegend(.7, 0.7, .9, 0.9);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);

  TString path= "/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/HNDiElectronOptimisation/periodBtoH/";
  TFile * fnp = new TFile(path+"HNDiElectronOptimisation_DoubleEG_SKnonprompt_dilep_cat_v8-0-6.root");
  TFile * fcf = new TFile(path+"HNDiElectronOptimisation_SKchargeflip_dilep_cat_v8-0-6.root");
  TFile * fmc = new TFile(path+"HNDiElectronOptimisation_mc_dilep_cat_v8-0-6.root");


  TH1F* h_npcutflow= (TH1F*)fnp->Get(("ISOcutflowdxysig__d0"));
  TH1F* h_cfcutflow= (TH1F*)fcf->Get(("ISOcutflowdxysig_"));
  TH1F* h_mccutflow= (TH1F*)fmc->Get(("ISOcutflowdxysig_"));

  for(unsigned int im=0; im < masses.size(); im++){
    TFile * fsig = new TFile((path+"HNDiElectronOptimisation_HNEmEm_"+ masses.at(im)+"_cat_v8-0-6.root"));
    TFile * fsig2 = new TFile((path+"HNDiElectronOptimisation_HNEpEp_"+ masses.at(im)+"_cat_v8-0-6.root"));

    
    cout << "\n ----------------" << endl;
    cout <<  "New mass : " << masses.at(im) << endl;
    TH1F* h_cutflow= (TH1F*)fsig->Get(("ISOcutflowdxysig_"));
    TH1F* h_cutflow2= (TH1F*)fsig2->Get(("ISOcutflowdxysig_"));
    
    h_cutflow->Add(h_cutflow2);
    
    vector<TString> isocuts;
    isocuts.push_back("050");
    isocuts.push_back("0525");
    isocuts.push_back("055");
    isocuts.push_back("060");
    isocuts.push_back("065");
    isocuts.push_back("075");
    isocuts.push_back("100");
    isocuts.push_back("125");



    vector<TString> binlabels;
    for(unsigned int iiso = 0 ; iiso < isocuts.size(); iiso++){
      for(unsigned int iiso2 = 0 ; iiso2 < isocuts.size(); iiso2++){
	binlabels.push_back(isocuts[iiso]+"_"+isocuts[iiso2]);
      }
    }


    TH1F* h_mass_isoopt = new TH1F(("isocut_opt_HN"+ masses.at(im)).Data(),("isocut_opt_HN"+ masses.at(im)).Data(), binlabels.size()+1,0.,double(binlabels.size()+1));
    int iisocut(0);
    for(unsigned int iiso = 0 ; iiso < isocuts.size(); iiso++){
      for(unsigned int iiso2 = 0 ; iiso2 < isocuts.size(); iiso2++, iisocut++){
	h_mass_isoopt->GetXaxis()->SetBinLabel(iisocut+1, binlabels[iisocut]);
	cout << iisocut+1 << " " << binlabels[iisocut] << endl;
      }
    }

    TH1F* h_cutflow_ref= (TH1F*)fsig->Get(("ISOREF"));

    
    int ibin(0);
    float maxcont(0.);
    for(int i= 0; i < binlabels.size()+1;i++){
      float sig_eff = h_cutflow->GetBinContent(i+1)/ h_cutflow_ref->GetBinContent(1);
      float tot_bkg = h_npcutflow->GetBinContent(i+1) + h_cfcutflow->GetBinContent(i+1)+ h_mccutflow->GetBinContent(i+1);
      float bkgtmp = tot_bkg + (0.3*h_npcutflow->GetBinContent(i+1))*(0.3*h_npcutflow->GetBinContent(i+1));
      float denom= 1. + sqrt(bkgtmp);
      float punzi = sig_eff/denom;

      if(punzi  > maxcont){
        maxcont=punzi;
	ibin=i;
	cout << sig_eff << " " << h_npcutflow->GetBinContent(ibin+1)  <<endl;

      }
      h_mass_isoopt->SetBinContent(i+1,punzi);
    } 
    cout << "Max punzi = " << maxcont << " bin = " << binlabels[ibin] <<  " " << endl;
    h_mass_isoopt->GetYaxis()->SetTitle("Punzi"); 
    h_mass_isoopt->GetXaxis()->SetTitle("d_{xy} cut");
    
    h_mass_isoopt->SetLineWidth(3.);
    if(im==0)  h_mass_isoopt->SetLineColor(kRed);
    if(im==1)  h_mass_isoopt->SetLineColor(kBlue);
    if(im==2)  h_mass_isoopt->SetLineColor(kGreen+4);
    if(im==3)  h_mass_isoopt->SetLineColor(kCyan);
    if(im==4)  h_mass_isoopt->SetLineColor(kSpring-2);
    if(im==5)  h_mass_isoopt->SetLineColor(kOrange);

    if(im==4)  h_mass_isoopt->SetLineStyle(2.);
    if(im==5)  h_mass_isoopt->SetLineStyle(4.);    

    //h_mass_isoopt->GetYaxis()->SetRangeUser(0.001, 0.005);
    h_mass_isoopt->GetYaxis()->SetTitleSize(0.05);
    h_mass_isoopt->GetYaxis()->SetTitleOffset(1.5);
    h_mass_isoopt->GetYaxis()->SetLabelSize(0.04);
    
    h_mass_isoopt->GetXaxis()->SetTitleSize(0.05);
    h_mass_isoopt->GetXaxis()->SetLabelSize(0.04);

    if(im==0)h_mass_isoopt->Draw("hist");
    else h_mass_isoopt->Draw("histsame");
    
    legend->AddEntry(h_mass_isoopt,("m_{N}= "+ masses.at(im)).Data(),"l");
  }
  
  legend->Draw("same");
  c1->SaveAs(("iso_dxysigopt.pdf"));
}



void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);


  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);


  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);
  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  //  tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);

  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);

  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.4);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset



  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.4);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}
