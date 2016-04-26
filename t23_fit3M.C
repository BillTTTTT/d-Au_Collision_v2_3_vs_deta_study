//-----------------------------------------------
//Code to graph v1 and v2 vs deta
//for AMPT Model for d+Au at 200GeV
//
//Author: P. Yin
//-----------------------------------------------


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <iterator>
#include <cmath>
#include <TLegend.h>
#include <TLatex.h>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TF1.h"
using namespace std;

const int nog = 4;

vector<TH1F*> dhis;
vector<TF1*> fit;

void t23_fit3M()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(111);

	TFile *fin = new TFile("d+Au_on_200GeV.root");
    //TFile *fref= new TFile("p+p_200GeV_reference.root");
    TFile *fref= new TFile("ampt_off_200GeV.root");

    TH1F *refs;
    TH1F *refb;
    TH1F *refcount;
    TH1F *refC;

	//Define 20 histograms
    for(int i=0; i<nog; i++)
    {
    	dhis.push_back((TH1F*)fin->Get(Form("d_%i",i)));
    }

    refs = (TH1F*)fref->Get("s");
    refb = (TH1F*)fref->Get("b");
    refcount = (TH1F*)fref->Get("hcount");
    refC = (TH1F*)refs->Clone("refC");

    int count = refcount->GetBinContent(2);

    refb->Scale(1.0/count);
    refC->Divide(refb);
    refC->Scale(refb->Integral()/refs->Integral());

    //refC->Draw();
    
    TF1 *fitref = new TF1("fitref","[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi());

    refC->Fit(fitref,"R");

    //Get the parameters from the function
    fitref = (TF1*)refC->GetFunction("fitref");
    float c = fitref->GetParameter(0);
    float x = fitref->GetParameter(1);
    float y = fitref->GetParameter(2);
    float z = fitref->GetParameter(3);
    
	//Define the fit graphs
    for(int i=0; i<nog; i++)
    {
        fit.push_back(new TF1(Form("f_%i",i),"[0]*(0.999996*(1 + 2*-0.0308293*cos(x) + 2*0.00434637*cos(2*x) + 2*-0.00176899*cos(3*x)))+[1]", -0.5*TMath::Pi(), 1.5*TMath::Pi()));
        //fit.push_back(new TF1(Form("f_%i",i),"[4]*([0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x)))+[5]", -0.5*TMath::Pi(), 1.5*TMath::Pi()));
        //fit[i]->SetParameter(0,c);
        //fit[i]->SetParameter(1,x);
        //fit[i]->SetParameter(2,y);
        //fit[i]->SetParameter(3,z);
    }
    

    TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
    c1->Divide(2,2);

    //TCanvas *c2 = new TCanvas("c1","c1",1600,1200);
    //c1->Divide(2,2);

    //Draw graphs
    for(int i=0; i<nog; i++)
    {   
        c1->cd(i+1);
        dhis[i]->Draw();

        //Fit the graphs
        dhis[i]->Fit(fit[i],"R");

        fit[i] = (TF1*)dhis[i]->GetFunction(Form("f_%i",i));

        dhis[i]->SetXTitle("#Delta#phi");

        //dhis[i]->Add(fit[i],-1);

        //c2->cd(i+1);
        //dhis[i]->Draw();
    }
	
    //c1->Print("222d+Au_off_20GeV.pdf");
    
}






















