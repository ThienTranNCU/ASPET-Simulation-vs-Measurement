//This program is writtern to run on Root to analysis simulation output from Gate ASPET 220312 CGMH
//This program will give: 2D back to back reconstruction, 1D projection, time vs detection, and total eneryg spectrum.
//Put all the analysis results in a root file and create a folder in the same directory as the input file
// to compile: g++ Read_SIM_ASPET_Coin_DAQ_v3.cc -o Read_SIM_ASPET_Coin_DAQ_v3 `root-config --cflags --glibs` -lSpectrum
// there are 6 arguments: 1. file path (required), 2. file name (required), 3. start time [sec](optional), 4. stop time [sec](optional), ...
//5. Energy threshold [keV](optional), 6. Efficiency correction text file (optional)
// Tran Cong Thien 230224
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
#include "TRandom2.h"
#include "TH1D.h"
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
#include "TLegend.h"
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
#include <complex>
#include <TApplication.h>
#include <TMath.h>


using namespace std;
////ChannelID converting

Int_t bin_X(Int_t cry,Int_t mod, Int_t res){
	Int_t M = mod/2;
	Int_t ch = cry/8;
	Int_t bx = M*8+ch;
	if(res==1){
			bx=31-bx;
	}
	return bx;
}

Int_t bin_Y(Int_t cry,Int_t mod, Int_t res){
	Int_t M = mod%2;
	Int_t ch = cry%8;
	Int_t by = M*8+ch;

	return by;
}
///////////////////////////////////////////////////////////// Main program
int main(int argc, char *argv[])
{
	
char fname[4][1000];
sprintf(fname[0],argv[2]);   															 //file name
sprintf(fname[1],argv[1]);																// file directory
sprintf(fname[2],"%s/%s.root",fname[1],fname[0]);							//open file address
cout<<"Open file:  "<<fname[2]<<endl;

/////////////////////////////////////////DAQ detector response

Int_t board0,px0,py0,termi=1;
Double_t DAQ,DAQco[2][32][16]={0},Kill[2][32][16]={1},kv;
 
for(Int_t b=0;b<2;b++){
	for(Int_t i=0;i<32;i++){
		for(Int_t j=0;j<16;j++){
			DAQco[b][i][j]=1;
			Kill[b][i][j]=1;
		}
	}
} 
 if(argc>=6){
sprintf(fname[3],argv[6]);																				//	DAQ efficiency file
ifstream inDAQ;
inDAQ.open(Form(fname[3]));	 
	 
 while(1){
	  if (!inDAQ.good()) break;  
	  inDAQ>>board0>>px0>>py0>>DAQ>>kv; 
	  DAQco[board0][px0][py0]=DAQ;	 
	  Kill[board0][px0][py0]=kv;	   
} 
inDAQ.close();
}
	
///////////////////////////// Time and Energy threshold

	
const Int_t  xb=32, yb=16;
	Float_t start =20;            //Start time Defaul
	Float_t stop =140;		// Stop time Defaul
	Float_t Eth= 400;    		// energy threshold keV Defaul

	if(argc>=4){start = strtod(argv[3],NULL);stop = strtod(argv[4],NULL);} //DAQ TIME
    if(argc>=5){Eth=strtod(argv[5],NULL);} //Energy threshold
    
	
 Int_t x1,x2,y1,y2,t,to=0,count=0,count1,count2;
 Double_t pX,pY,pmX,pmY;

 Int_t const tbin=((stop-start)+40);
 TH1D *ht1=new TH1D("ht1","Detection rate",tbin,start-20,stop+20);
 ht1->GetXaxis()->SetTitle("Time (sec)");
 

 TH1D *he1=new TH1D("he1","Total energy spectrum",1000,0,1000);
he1->GetXaxis()->SetTitle("Energy (keV)");

 TH1D *ht2=new TH1D("ht2","Detection rate",tbin,start-20,stop+20);
 ht2->GetXaxis()->SetTitle("Time (sec)");
 

 TH1D *he2=new TH1D("he2","Total energy spectrum",1000,0,1000);
he2->GetXaxis()->SetTitle("Energy (keV)");

TH1D *hc=new TH1D("hc"," ",xb*2,0,32);
hc->GetXaxis()->SetTitle("X channel");
TH1D *hk=new TH1D("hk"," ",xb*2,0,32);
hk->GetXaxis()->SetTitle("X channel");
TH1D *hg=new TH1D("hg"," ",xb*2,0,32);

Int_t thc[64]={0};
Int_t thk[64]={0};
Int_t thg[64]={0};

 TH2D *h1=new TH2D("h1","Board 0",xb,0,32,yb,0,16);
 h1->GetYaxis()->SetTitle("Y channel");
h1->GetXaxis()->SetTitle("X channel");
 TH2D *h2=new TH2D("h2","Board 1",xb,0,32,yb,0,16);
h2->GetYaxis()->SetTitle("Y channel");
h2->GetXaxis()->SetTitle("X channel");
Double_t dis=0;	

 TH2D *hc1=new TH2D("hc1","DAQ adapted reconstruction",xb*2,0,32,yb*2,0,16);
hc1->GetYaxis()->SetTitle("Y channel");
hc1->GetXaxis()->SetTitle("X channel");

 TH2D *hg1=new TH2D("hg1","Global position reconstructed",xb*2,0,102.4,yb*2,-25.6,25.6);
hg1->GetYaxis()->SetTitle("Y channel");
hg1->GetXaxis()->SetTitle("X channel");

 TH2D *hk1=new TH2D("hk1","Dead modules killed reconstructtion",xb*2,0,32,yb*2,0,16);
hk1->GetYaxis()->SetTitle("Y channel");
hk1->GetXaxis()->SetTitle("X channel");

	///////////////////////out put
//ofstream out;
//out.open(Form("DETECTOR_DATA_PMMA_80MeV_first300events.txt"));

////##############################################################################################################simulation

	Int_t        nbytes = 0;
	Float_t          energy1;
    Double_t         time1;
    Int_t            runID;
    Int_t            eventID1;
    Float_t          globalPosX1,globalPosY1,globalPosZ1;
    Int_t            gantryID1,crystalID1,moduleID1,rsectorID1;

Float_t          energy2;
Double_t         time2;

    Int_t            eventID2;
    Float_t          globalPosX2,globalPosY2,globalPosZ2;
    Int_t            gantryID2,crystalID2,moduleID2,rsectorID2;

      
    //TFile *f = new TFile("Gamma_220105.root");
    TFile *f = new TFile(fname[2]);
    TTree *Coin = (TTree*)gDirectory->Get("Coincidences;1");
    Coin->SetBranchAddress("eventID1",&eventID1);
    Coin->SetBranchAddress("crystalID1",&crystalID1);
    Coin->SetBranchAddress("gantryID1",&gantryID1);
    Coin->SetBranchAddress("energy1",&energy1);
    Coin->SetBranchAddress("moduleID1",&moduleID1);
    Coin->SetBranchAddress("time1",&time1);
    Coin->SetBranchAddress("rsectorID1",&rsectorID1);
    Coin->SetBranchAddress("runID",&runID);
    Coin->SetBranchAddress("globalPosX1",&globalPosX1);
    Coin->SetBranchAddress("globalPosY1",&globalPosY1);
    Coin->SetBranchAddress("globalPosZ1",&globalPosZ1);
    
    Coin->SetBranchAddress("eventID2",&eventID2);
    Coin->SetBranchAddress("crystalID2",&crystalID2);
    Coin->SetBranchAddress("gantryID2",&gantryID2);
    Coin->SetBranchAddress("energy2",&energy2);
    Coin->SetBranchAddress("moduleID2",&moduleID2);
    Coin->SetBranchAddress("time2",&time2);
    Coin->SetBranchAddress("rsectorID2",&rsectorID2);
    Coin->SetBranchAddress("globalPosX2",&globalPosX2);
    Coin->SetBranchAddress("globalPosY2",&globalPosY2);
    Coin->SetBranchAddress("globalPosZ2",&globalPosZ2);
    
    Int_t px,py,mx,my,cx,cy,bx1,by1,bx2,by2;
    Float_t pgx,pgy,pcx,pcy,r,cr;
    
    TRandom2 *rand = new TRandom2(0);              ////random coin
    
    Int_t nentries = Coin->GetEntries();
 //  Int_t Zmean,ymean;
     for (Int_t i=0; i<nentries;i++) {
        nbytes += Coin->GetEntry(i);
        
        if((time1>=start)&&(time2>=start)&&(time1<=stop)&&(time2<=stop)&&(energy1*1000>=Eth)&&(energy2*1000>=Eth)){    ////Time and energy threshold
       
       	bx1= bin_X(crystalID1,moduleID1,rsectorID1);
		by1= bin_Y(crystalID1,moduleID1,rsectorID1);
   
		bx2= bin_X(crystalID2,moduleID2,rsectorID2);
		by2= bin_Y(crystalID2,moduleID2,rsectorID2);
		
		pcx=(bx1+bx2)/2.0;
		pcy=(by1+by2)/2.0;
		
		termi=Kill[rsectorID1][bx1][by1]*Kill[rsectorID2][bx2][by2];
		
		if(termi==1){
			 hk1->Fill(pcx,pcy);
			 he2->Fill(energy1*1000);
		     he2->Fill(energy2*1000);
		     ht2->Fill(time1);
	        ht2->Fill(time2);
			 
		}
		
		
		r = rand->Rndm();
		cr=DAQco[rsectorID1][bx1][by1]*DAQco[rsectorID2][bx2][by2]*Kill[rsectorID1][bx1][by1]*Kill[rsectorID2][bx2][by2];	
				
		if(r<=cr){     // DAQ conditon
			
		he1->Fill(energy1*1000);
		he1->Fill(energy2*1000);
		ht1->Fill(time1);
		ht1->Fill(time2);
			                                                                                    
		if(rsectorID1==0){
			h1->Fill(bx1,by1);
		}else{
			h2->Fill(bx1,by1);
		}
		
		if(rsectorID2==0){
			h1->Fill(bx2,by2);
		}else{
			h2->Fill(bx2,by2);
		}
		
		
        hc1->Fill(pcx,pcy);
	}
	
	}
}
//out.close();

for(Int_t i=0;i<64;i++){
	for(Int_t j=0;j<32;j++){
		thc[i]+=hc1->GetBinContent(i+1,j+1);
		thk[i]+=hk1->GetBinContent(i+1,j+1);
	}
}

for(Int_t i=0;i<64;i++){
	hc->SetBinContent(i+1,thc[i]);
	hk->SetBinContent(i+1,thk[i]);
}

char tname[1000];
char dname[200];

if(argc>=6){
	sprintf(dname,"Coincidences_DAQ");
}else{
	sprintf(dname,"Coincidences");
}

sprintf(tname,"mkdir %s/%s_%s",fname[1], fname[0],dname);
system(tname);

TCanvas* c1 = new TCanvas("c1","histogram",100,10,1000,500);//700-rong, 900-dai window hien thi
gStyle->SetPalette(55);

hc1->Draw("colz");
sprintf(tname,"%s/%s_%s/Coincidences_DAQ_2DMaps.png",fname[1], fname[0],dname);
c1->SaveAs(tname);
sprintf(tname,"%s/%s_%s/Coincidences_DAQ_2DMaps.root", fname[1],fname[0],dname);
c1->SaveAs(tname);

TCanvas* k1 = new TCanvas("k1","histogram",100,10,1000,500);//700-rong, 900-dai window hien thi
gStyle->SetPalette(55);

hk1->Draw("colz");
sprintf(tname,"%s/%s_%s/Coincidences_Kill_2DMaps.png",fname[1], fname[0],dname);
k1->SaveAs(tname);
sprintf(tname,"%s/%s_%s/Coincidences_Kill_2DMaps.root", fname[1],fname[0],dname);
k1->SaveAs(tname);

TCanvas* c1d = new TCanvas("c1d","histogram",100,10,1000,1000);//700-rong, 900-dai window hien thi
gStyle->SetPalette(55);
c1d->Divide(1,2);
c1d->cd(1);
h1->Draw("colz");
c1d->cd(2);
h2->Draw("colz");

sprintf(tname,"%s/%s_%s/LOR_2DMaps.png",fname[1], fname[0],dname);
c1d->SaveAs(tname);
sprintf(tname,"%s/%s_%s/LOR_2DMaps.root",fname[1], fname[0],dname);
c1d->SaveAs(tname);

TCanvas* c2 = new TCanvas("c2","histogram",100,10,1000,500);
c2->SetGrid();
ht2->Draw();
ht1->SetLineColor(2);
ht1->SetLineWidth(2);
ht1->Draw("same");
ht2->SetLineColor(4);
ht2->SetLineWidth(2);

TLegend *legt = new TLegend(0.75,0.75,0.95,0.95);
legt->AddEntry(ht1,"DAQ adapted","L");
legt->AddEntry(ht2,"Dead modules killed","L");
legt->Draw();
sprintf(tname,"%s/%s_%s/Det_Time.png", fname[1],fname[0],dname);
c2->SaveAs(tname);
sprintf(tname,"%s/%s_%s/Det_Time.root",fname[1], fname[0],dname);
c2->SaveAs(tname);

TCanvas* c3 = new TCanvas("c3","histogram",100,10,1000,500);
c3->SetGrid();
he2->Draw();
he1->SetLineColor(2);
he1->SetLineWidth(2);
he1->Draw("same");
he2->SetLineColor(4);
he2->SetLineWidth(2);

TLegend *lege = new TLegend(0.75,0.75,0.95,0.95);
lege->AddEntry(he1,"DAQ adapted","L");
lege->AddEntry(he2,"Dead modules killed","L");
lege->Draw();
sprintf(tname,"%s/%s_%s/Energy_Spectrum.png",fname[1], fname[0],dname);
c3->SaveAs(tname);
sprintf(tname,"%s/%s_%s/Energy_Spectrum.root", fname[1],fname[0],dname);
c3->SaveAs(tname);


TCanvas* c4 = new TCanvas("c4","histogram",100,10,1000,500);
c4->SetGrid();
hk->Draw();
hc->SetLineColor(2);
hc->SetLineWidth(2);
hc->Draw("same");
hk->SetLineColor(4);
hk->SetLineWidth(2);

TLegend *leg = new TLegend(0.75,0.75,0.95,0.95);
leg->AddEntry(hc,"DAQ adapted","L");
leg->AddEntry(hk,"Dead modules killed","L");
leg->Draw();

sprintf(tname,"%s/%s_%s/Coincidences_1DMaps.png",fname[1], fname[0],dname);
c4->SaveAs(tname);
sprintf(tname,"%s/%s_%s/Coincidences_1DMaps.root", fname[1],fname[0],dname);
c4->SaveAs(tname);

return 0;

}
