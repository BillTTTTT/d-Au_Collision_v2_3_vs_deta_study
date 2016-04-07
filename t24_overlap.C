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

void t24_overlap()
{
	gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0.0001);
    gStyle->SetLegendTextSize(0.03);

	TFile *fin1 = new TFile("d+Au_on_20GeV.root");
	TFile *fin2 = new TFile("d+Au_off_20GeV.root");

	TH1F *h11;
	TH1F *h12;
	TH1F *h21;
	TH1F *h22;

	h11 = (TH1F*)fin1->Get("t120");
	h21 = (TH1F*)fin2->Get("t120");

	h12 = (TH1F*)fin1->Get("t220");
	h22 = (TH1F*)fin2->Get("t220");

	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	h11->SetMarkerStyle(20);
	h11->Draw();
	h21->SetMarkerStyle(24);
	h21->Draw("same");

	TLegend *leg1 = new TLegend(0.65,0.15,0.85,0.35);
	leg1->AddEntry(h11,"with cascade","lep");
    leg1->AddEntry(h21,"without cascade","lep");
    leg1->Draw("same");

	c1->Print("v1_20GeV.pdf");

	TCanvas *c2 = new TCanvas("c2","c2",800,600);
	h12->SetMarkerStyle(20);
	h12->GetYaxis()->SetRangeUser(0,0.005);
	h12->Draw();
	h22->SetMarkerStyle(24);
	h22->Draw("same");

	TLegend *leg2 = new TLegend(0.15,0.65,0.35,0.85);
	leg2->AddEntry(h12,"with cascade","lep");
    leg2->AddEntry(h22,"without cascade","lep");
    leg2->Draw("same");

    c2->Print("v2_20GeV.pdf");

}

















