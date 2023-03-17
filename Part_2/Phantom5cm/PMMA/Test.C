
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TH1F.h"
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
#include "TPaveStats.h"


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


void Test(){
    
    Float_t          energy1,energy2;
    Double_t         time1,time2;
    Int_t            runID;
    Int_t            eventID1,eventID2;
    Float_t          globalPosX1,globalPosY1,globalPosZ1,globalPosX2,globalPosY2,globalPosZ2;
    
    
    Int_t            nbytes = 0;
    Int_t            moduleID1,rsectorID1;
    Int_t            crystalID1;
    Int_t            gantryID1;
    
    
    //TFile *f = new TFile("Gamma_220105.root");
    TFile *f = new TFile("DETECTOR_DATA_PMMA_100MeV_C11.root");
    TTree *Coin = (TTree*)gDirectory->Get("Coincidences");
    Coin->SetBranchAddress("eventID1",&eventID1);
    Coin->SetBranchAddress("crystalID1",&crystalID1);
    Coin->SetBranchAddress("gantryID1",&gantryID1);
    Coin->SetBranchAddress("energy1",&energy1);
    Coin->SetBranchAddress("moduleID1",&moduleID1);
    Coin->SetBranchAddress("time1",&time1);
    Coin->SetBranchAddress("rsectorID1",&rsectorID1);
    
    Coin->SetBranchAddress("runID",&runID);
    Coin->SetBranchAddress("eventID1",&eventID1);
    Coin->SetBranchAddress("energy1",&energy1);
    Coin->SetBranchAddress("time1",&time1);
    Coin->SetBranchAddress("globalPosX1",&globalPosX1);
    Coin->SetBranchAddress("globalPosY1",&globalPosY1);
    Coin->SetBranchAddress("globalPosZ1",&globalPosZ1);
    Coin->SetBranchAddress("eventID2",&eventID2);
    Coin->SetBranchAddress("energy2",&energy2);
    Coin->SetBranchAddress("time2",&time2);
    Coin->SetBranchAddress("globalPosX2",&globalPosX2);
    Coin->SetBranchAddress("globalPosY2",&globalPosY2);
    Coin->SetBranchAddress("globalPosZ2",&globalPosZ2);
    
    
    
    int divisor=2;int quotient; int remainder;
    Float_t Eth=0.200;
    Float_t Ethresh[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5}; //11 values of threshold specified
    //const int Nthresh=sizeof(Ethresh)/sizeof(Ethresh[0]);
    int Nthresh=11;
    cout<<"Number of E thresholds="<<Nthresh<<endl;
    Float_t nCoin[11]; for(int i=0; i<Nthresh; i++){nCoin[i]=0;}
    //Float Ethreshmean=0;
    
    
    
    TH2F *hist1 = new TH2F("hist1", "0 cm 100 MeV",64,0, 102.4, 32, -25.6, 25.6);
    TH1F *hcoin1 = new TH1F("hcoin1","hcoin1",2000,-1000,1000);
    TH1F *E1 = new   TH1F("E1","Gamma Spectrum",60,0,0.6);
    TH1F *E2 = new   TH1F("E2","Gamma Spectrum",60,0,0.6);
    char hname[100];
    TH1F *Eall[11];
    for(int i=0;i<11; i++)
    {sprintf(hname,"Eth%d",i);
        Eall[i] = new TH1F(hname,hname,120,0,0.6);
    }
    
    Int_t nentries = Coin->GetEntries();
    int Zmean,ymean;
    float tdelay=0;
    Double_t ns = 1*pow(10,-12);
    cout<<"Total coincidence events="<<nentries<<endl;
    cout<<"List of the first ten coincidences:";
    cout<<"time1"<<"\t"<<"time2"<<"\t"<<"tdelay"<<"\t"<<"energy1"<<"\t"<<"energy2"<<endl;
    for (Int_t i=0; i<nentries;i++) {
        nbytes += Coin->GetEntry(i);
        tdelay = (time1-time2)/ns;
        //coincidences are printed
        if(i<10){cout<<time1<<"\t"<<time2<<"\t"<<tdelay<<"\t"<<energy1<<"\t"<<energy2<<endl;}
        float Ymean=(globalPosY1+globalPosY2)/(2);
        float Zmean=(globalPosZ1+globalPosZ2)/(2);
        
        //Espectra for module1 and module2 for Eth threshold
        //Fill PROJECTED SOURCE ONLY FOR SET THRESHOLD
        if(energy1>Eth&&energy2>Eth)
        {
            hist1->Fill(Zmean,Ymean);
            hcoin1->Fill(tdelay);
            E1->Fill(energy1);
            E2->Fill(energy2);
        }
        //E spectra for various thresholds
        for(int i=0; i<Nthresh; i++)
        {
            if(energy1>=Ethresh[i]&&energy2>=Ethresh[i])
            { nCoin[i]=nCoin[i]+1;
                Eall[i]->Fill(energy1);  //Fills only one module of the detector
            }
        }
        //hcoin1->Fill(-tdelay);
        //}
    }
    
    TCanvas *C0=new TCanvas("C0", "Plane0",1000,500);
    C0->cd();
    
    gStyle->SetPalette(kRainBow);
    
    hist1->GetXaxis()->SetTitle("Z Axis (mm)");
    hist1->GetYaxis()->SetTitle("Y Axis (mm)");
    hist1->Draw();
    hist1->Draw("colz");
    C0->Print("Source Projection.png");
    C0->SaveAs("Source Projection.root");
    
    
    TCanvas *CE1=new TCanvas("CE1", "Plane0",500,500);
    //    E1->SetLineColor(kBlue);E1->Draw();
    //    E1->SetLineColor(kRed);E2->Draw("SAME");
    for(int i=0; i<Nthresh; i++)
    {
        Eall[i]->SetLineColor(2*i);
        Eall[i]->Draw("SAME");
    }
    CE1->cd();
    
    TCanvas *Ccoin=new TCanvas("Ccoin","Ccoin",500,500);
    hcoin1->Draw();
    
    TCanvas *Cthresh = new TCanvas("Cthresh","Coin vs Thresh",500,500);
    TGraph* grCthresh = new TGraph(Nthresh,Ethresh,nCoin);
    grCthresh->Draw("AC*");
    cout<<"Ethresh(MeV)"<<"\t No of coincidences"<<endl;
    for(int i=0; i<Nthresh; i++)
    {
        cout<<Ethresh[i]<<"\t"<<nCoin[i]<<endl;
    }
    
    
    
}



