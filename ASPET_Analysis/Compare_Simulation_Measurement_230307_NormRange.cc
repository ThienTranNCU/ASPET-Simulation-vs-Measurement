//Analysis is ased on the Singles

#include "Riostream.h"
#include <string.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>

#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TRint.h"
#include "TObject.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TSpectrum.h"


#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TBrowser.h"

#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <complex>
#include <TApplication.h>
#include <TMath.h>


void Compare_Simulation_Measurement_230307_NormRange(){
	
Double_t dif=0;
TH1D *h1m=new TH1D("h1m"," ",64,0,32);
TH1D *h1b=new TH1D("h1b"," ",64,0,32);
   
 h1m->GetXaxis()->SetTitle("X index");
Double_t M1D[64]={0},B1D[64]={0};

TFile f1("230228_Data_Simulation_ASPET_CGMH220312/Phantom5cm/PMMA/DETECTOR_DATA_100MeV_2mins_230220_Coincidences_DAQ/Coincidences_1DMaps.root","read");
TCanvas *c1 = (TCanvas*)f1.Get("c4;1");
TH1D *h1 = (TH1D*)c1->GetPrimitive("hc");
TH1D *h2 = (TH1D*)c1->GetPrimitive("hk");

TFile f2("Data_CGMH_220312/100104.298_A3-2MINS_Processed/All_outputTree.root","read");
TH2D *h3 = (TH2D*)f2.Get("h2Map;1");

TFile f3("Data_CGMH_220312/101740.889_C1BG_Processed/All_outputTree.root","read");
TH2D *hb = (TH2D*)f3.Get("h2Map;1");

hb->Scale(0.172166);

for(Int_t i=0;i<64;i++){
		for(Int_t j=0;j<32;j++){
			M1D[i] += h3->GetBinContent(i+1,j+1);
			B1D[i] += hb->GetBinContent(i+1,j+1);
		}
}
Double_t sim;
dif=0;
for(Int_t i=0;i<64;i++){
	sim=h2->GetBinContent(i+1);
	M1D[i]=M1D[i]-B1D[i];
	h1m->SetBinContent(i+1,M1D[i]);
	h1b->SetBinContent(i+1,B1D[i]);
	if(M1D[i]!=0){
	dif=dif+ (abs(sim*1.0-M1D[i]*1.0))/(M1D[i]*1.0);
	}
	//cout<<M1D[i]<<" "<<sim<<" "<<dif<<endl;
}
dif=dif/63.0;

Double_t sc= 1+dif;
h2->Scale(sc);

         TF1 *fGs1 = new TF1("fGs1","[0]+([1]-[0])/(1+exp([2]*(x-[3])))",18.0,25.);
//    fG->SetParameters(-1.5,10.0);
		fGs1->SetParameters(5000.0,1.0,0.3,11.0);

 //   h2->Fit(fGs1,"RQ");
    Double_t p1=fGs1->GetParameter(3);
    Double_t e1=fGs1->GetParError(3);
    Double_t p0=fGs1->GetParameter(2);
	Double_t e0=fGs1->GetParError(2);
 cout<<p0<<" "<<e0<<" "<<p1<<" "<<e1<<endl;

TF1 *fGs2 = new TF1("fGs2","[0]+([1]-[0])/(1+exp([2]*(x-[3])))",0,6);
//    fG->SetParameters(-1.5,10.0);
		fGs2->SetParameters(5000.0,1.0,0.3,11.0);

 //   h2->Fit(fGs2,"RQ+");
    p1=fGs2->GetParameter(3);
    e1=fGs2->GetParError(3);
    p0=fGs2->GetParameter(2);
	e0=fGs2->GetParError(2);
 cout<<p0<<" "<<e0<<" "<<p1<<" "<<e1<<endl;


         TF1 *fGm1 = new TF1("fGm1","[0]+([1]-[0])/(1+exp([2]*(x-[3])))",18,25.);
//    fG->SetParameters(-1.5,10.0);
		fGm1->SetParameters(5000.0,1.0,0.3,11.0);

//    h1m->Fit(fGm1,"RQ");
    p1=fGm1->GetParameter(3);
    e1=fGm1->GetParError(3);
    p0=fGm1->GetParameter(2);
	e0=fGm1->GetParError(2);
 cout<<p0<<" "<<e0<<" "<<p1<<" "<<e1<<endl;
 
          TF1 *fGm2 = new TF1("fGm2","[0]+([1]-[0])/(1+exp([2]*(x-[3])))",0,6);
//    fG->SetParameters(-1.5,10.0);
		fGm2->SetParameters(5000.0,1.0,0.3,11.0);

//    h1m->Fit(fGm2,"RQ+");
    p1=fGm2->GetParameter(3);
    e1=fGm2->GetParError(3);
    p0=fGm2->GetParameter(2);
	e0=fGm2->GetParError(2);
 cout<<p0<<" "<<e0<<" "<<p1<<" "<<e1<<endl;
 

TCanvas* c11 = new TCanvas("c11","1D projection",100,10,1000,500);//700-rong, 900-dai window hien thi
c11->SetGrid(1);
h2->Draw();
h1m->SetLineColor(2);
h1m->SetLineWidth(2);
//h1->Draw("same");
h1m->Draw("l same");
//h1b->Draw("same");
h1->SetLineColor(kGreen+2);
h2->SetLineColor(4);
h1b->SetLineColor(1);
h2->SetLineWidth(2);
h1->SetLineWidth(2);
h1b->SetLineWidth(2);

TLegend *leg = new TLegend(0.75,0.75,0.98,0.95);
leg->AddEntry(h1m,"Measurement","L");
leg->AddEntry(h2,"Simulation (mod 1 - scaled)","L");
//leg->AddEntry(h2,"Simulation dead channels killed","L");
//leg->AddEntry(hb,"Background","L");

leg->Draw();
//c11->SaveAs("Compare_230308_NormRange/1D_projection.png");

//TCanvas* c12 = new TCanvas("c12","2D reconstruction",100,10,1000,500);//700-rong, 900-dai window hien thi
//gStyle->SetPalette(55);
//h3->Draw("colz");
//c12->SaveAs("Compare_230308_NormRange/2DReconstructed_measurement.png");

}
