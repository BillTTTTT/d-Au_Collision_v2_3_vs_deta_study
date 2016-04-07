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

	TFile *fin = new TFile("d+Au_off_20GeV.root");

	//Define 20 histograms
    for(int i=0; i<nog; i++)
    {
    	dhis.push_back((TH1F*)fin->Get(Form("d_%i",i)));
    }
    
	//Define the fit graphs
    for(int i=0; i<nog; i++)
    {
        fit.push_back(new TF1(Form("f_%i",i),"[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi()));
    }
	
    TNtuple *tuple;
    tuple = (TNtuple*)fin->Get("v12v22v32");

    float v12;
    float v22;
    float v32;
    float sv12[4];
    float sv22[4];
    float sv32[4];

    tuple->SetBranchAddress("v12", &v12);
    tuple->SetBranchAddress("v22", &v22);
    tuple->SetBranchAddress("v32", &v32);

    for(int i=0; i<4; i++)
    {
        for(int j=i; j<tuple->GetEntries(); j=j+4)
        {
            tuple->GetEntry(j);
            sv12[i] = sv12[i] + v12;
            sv22[i] = sv22[i] + v22;
            sv32[i] = sv32[i] + v32;
        }
    }

    TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
    c1->Divide(2,2);

    TLatex latex;

    //Show the number of v12, v22, v32
    for(int i=0; i<nog; i++)
    {
        sv12[i] = sv12[i] / 1000.0;
        char strv12[10];
        sprintf(strv12, "%f", sv12[i]);
        sv22[i] = sv22[i] / 1000.0;
        char strv22[10];
        sprintf(strv22, "%f", sv22[i]);
        sv32[i] = sv32[i] / 1000.0;
        char strv32[10];
        sprintf(strv32, "%f", sv32[i]);

        int markx = i%2;
        int marky = i/2;
        float x1 = 0.08 + 0.5 * markx;
        float x2 = 0.15 + 0.5 * markx;
        float y1 = 0.90 - 0.5 * marky;
        float y2 = 0.86 - 0.5 * marky;
        float y3 = 0.82 - 0.5 * marky;
        
        latex.SetTextSize(0.02);
        latex.SetTextAlign(11);
        latex.DrawLatex(x1, y1,"#LTcos(1#Delta#phi)#GT = ");
        latex.DrawLatex(x2, y1,strv12);
        latex.DrawLatex(x1, y2,"#LTcos(2#Delta#phi)#GT = ");
        latex.DrawLatex(x2, y2,strv22);
        latex.DrawLatex(x1, y3,"#LTcos(3#Delta#phi)#GT = ");
        latex.DrawLatex(x2, y3,strv32);
    }

    //Draw graphs
    for(int i=0; i<nog; i++)
    {   
        c1->cd(i+1);
        dhis[i]->Draw();

        //Fit the graphs
        dhis[i]->Fit(fit[i],"R");

        //Get the parameters from the function
        fit[i] = (TF1*)dhis[i]->GetFunction(Form("f_%i",i));
        float c = fit[i]->GetParameter(0);
        float x = fit[i]->GetParameter(1);
        float y = fit[i]->GetParameter(2);
        float z = fit[i]->GetParameter(3);

        //Draw only v12 fit function
        TF1 *f1 = new TF1("f1","[0]*(1 + 2*[1]*cos(x))",-0.5*TMath::Pi(), 1.5*TMath::Pi()); 
        f1->SetParameter(0,c);
        f1->SetParameter(1,x);
        f1->Draw("same");
        f1->SetLineColor(kGreen);

        //Draw only v22 fit function
        TF1 *f2 = new TF1("f2","[0]*(1 + 2*[1]*cos(2*x))",-0.5*TMath::Pi(), 1.5*TMath::Pi()); 
        f2->SetParameter(0,c);
        f2->SetParameter(1,y);
        f2->Draw("same");
        f2->SetLineColor(kBlack);

        //Draw only v32 fit function
        TF1 *f3 = new TF1("f3","[0]*(1 + 2*[1]*cos(3*x))",-0.5*TMath::Pi(), 1.5*TMath::Pi()); 
        f3->SetParameter(0,c);
        f3->SetParameter(1,z);
        f3->Draw("same");
        f3->SetLineColor(kOrange);

        TF1 *f4 = new TF1("f4","[0]*1",-0.5*TMath::Pi(), 1.5*TMath::Pi()); 
        f4->SetParameter(0,c);
        f4->SetLineStyle(7);
        f4->SetLineWidth(1);
        f4->Draw("same");
        f4->SetLineColor(kBlue);

        
        dhis[i]->SetXTitle("#Delta#phi");

        TLegend *leg = new TLegend(0.65,0.15,0.85,0.35);
        leg->AddEntry(fit[i],"fit line","L");
        leg->AddEntry(f1,"cos(1#Delta#phi)","L");
        leg->AddEntry(f2,"cos(2#Delta#phi)","L");
        leg->AddEntry(f3,"cos(3#Delta#phi)","L");
        leg->AddEntry(f4,"Horizantal line of c1","L");
        leg->Draw("same");
    }
	
    c1->Print("d+Au_off_20GeV.pdf");

}






















