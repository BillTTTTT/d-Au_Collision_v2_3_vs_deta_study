//-----------------------------------------------
//Code to graph v3 and v2 vs deta
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

//------------------------
//Variables
//------------------------

struct particle
{
    int   id;
    float px;
    float py;
    float pz;
    float x;
    float y;
    float z;
    float eta;
    float phi;
    float pT;
};

//Number of nucleons in system
// --> p+Au = 198
// --> d+Au = 199
// --> d+Pb = 210
// --> p+Pb = 209
const int NUCL = 199;

const int nog = 4;

vector<particle> particles;
int n = 0;

vector<TH1F*> dhis2;
vector<TF1*>  fit;

TH1F *hcount;

int   ct[10];
float v12[10];
float v22[10];
float v32[10];

TProfile* t220;
TProfile* t120;

void calvs(int index, float dphi)
{
    v12[index] = v12[index] + TMath::Cos(1 * dphi);
    v22[index] = v22[index] + TMath::Cos(2 * dphi);
    v32[index] = v32[index] + TMath::Cos(3 * dphi);
    ct[index]++;
}

void processEvent()
{
    if(particles.size() == 0) return;

    //Calculate the pseudo-rapidity gap and store particles into right bin
    for(unsigned int i=0; i<particles.size(); i++)
    {
        for(unsigned int j=0; j<particles.size(); j++)
        {
            if(i<=j) continue;

            float deta;
            float dphi;
            int   index;

            //Calculate the pseudo-rapidity gap
            deta = abs(particles[i].eta - particles[j].eta);

            //Get rad of deta>2
            if(deta>2) continue;

            dphi = particles[i].phi - particles[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            index = deta / (2.0 / nog);

            calvs(index, dphi);

            dhis2[index]->Fill(dphi);

            t220->Fill(deta,TMath::Cos(2 * dphi));
            t120->Fill(deta,TMath::Cos(1 * dphi));

            int i1;

            i1 = deta / (2.0 / 10);

            hcount->Fill(i1);
        }
    }

}

 void parseampt()
{
    //Read in data file
    ifstream dataFile;
    //dataFile.open("ana/ampt.dat");
    dataFile.open("/Users/Bill/Desktop/Lab/d-Au_Collision_v2_3_vs_deta_study/10K_Data/ampt_0.dat");

    //Skip the job if not dataFile
    if (!dataFile)  return;

    //In this while loop, program will read the data file line by line
    while(dataFile)
    {
        int    evtnumber;
        int    testnum;
        int    nlist;
        double impactpar;
        int    npartproj;
        int    nparttarg;
        int    npartprojelas;
        int    npartprojinelas;
        int    nparttargelas;
        int    nparttarginelas;
        double junk;

        //Get the header of each event
        dataFile >> evtnumber >> testnum >> nlist >> impactpar >> npartproj >> nparttarg >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas >> junk;

        if (!dataFile) break;

        //Analysis each particle in the event
        for (int i=0; i<nlist; i++)
        {
            int partid;
            float pv[3];
            float mass;
            double space[4];

            dataFile >> partid >> pv[0] >> pv[1] >> pv[2] >> mass >> space[0] >> space[1] >> space[2] >> space[3];

            //Skip non-charged particles that we are not interested in
            //+-211 are pions,  +-321 are kaons, +-2212 are protons
            if(abs(partid) != 211 && abs(partid) != 321 && abs(partid) != 2212) continue;

            if(pv[2] > 99.9) continue;

            //Calculate the energy
            //float energy = TMath::Sqrt(pv[0]*pv[0] + pv[1]*pv[1] + pv[2]*pv[2]);

            float pt = TMath::Sqrt(pv[0]*pv[0] + pv[1]*pv[1]);

            if(pt==0) continue;

            //Make Lorentz vector
            TLorentzVector ev(pv[0], pv[1], pv[2], 0);

            //Get pT, phi, pseudorapidity, particle id, px, py, and pz. Store them into p.
            particle p;

            p.eta = ev.Eta();
            p.pT  = pt;
            p.phi = ev.Phi();
            p.px  = pv[0];
            p.py  = pv[1];
            p.pz  = pv[2];
            p.x   = space[0];
            p.y   = space[1];
            p.z   = space[2];

            //if(p.pT == 0) continue;
            if(p.pT <= 0.2) continue;
            if(abs(p.eta) >= 1) continue;

            particles.push_back(p);
        }

        processEvent();
        particles.clear();

        if (!dataFile) break;
    }
}


void t22_dphi_ampe_condor()
{
    //Define 4 histograms
    for(int i=0; i<nog; i++)
    {
        dhis2.push_back(new TH1F(Form("d_%i",i),  Form("v2: dphi in range of [%f,%f)",i*0.1,(i+1)*0.1), 50, -0.5*TMath::Pi(), 1.5*TMath::Pi()));
    }

    //Define count histgram
    hcount = new TH1F("hcount", "count", 10, 0, 10);

    t220 = new TProfile("t220","t220",10,-0.01,1.99,1,+1);
    t120 = new TProfile("t120","t120",10,-0.01,1.99,1,+1);

    parseampt();

    TNtuple *tuple = new TNtuple("v12v22v32", "v12v22v32", "v12:v22:v32", 32000);

    //Show the number of v12, v22, v32
    for(int i=0; i<nog; i++)
    {
        v12[i] = v12[i] / ct[i];
        v22[i] = v22[i] / ct[i];
        v32[i] = v32[i] / ct[i];
        tuple->Fill(v12[i],v22[i],v32[i]);
    }

    //Make a file to store outputs
    TFile *fout = new TFile("out.root","RECREATE");

    //Save graphs
    for(int i=0; i<nog; i++)
    {   
        dhis2[i]->Write();
    }
    t220->Write();
    t120->Write();
    tuple->Write();
    hcount->Write();
    fout->Close();
}





















