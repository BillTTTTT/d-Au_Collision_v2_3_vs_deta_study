//-----------------------------------------------
//Code to graph dphi distribution
//for AMPT Model for d+Au at 200GeV
//with new reference fitting method
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
    float rsquare;
    float psi2;
    float psi3;
};

//Number of nucleons in system
// --> p+Au = 198
// --> d+Au = 199
// --> d+Pb = 210
// --> p+Pb = 209
const int NUCL = 199;

const int nog = 4;

vector<particle> particles;
vector<particle> p1;
vector<particle> p2;
vector<particle> p3;
vector<particle> p4;

vector<TH1F*> dhiss;
vector<TH1F*> dhisb;
TH1F *hcount;

void move()
{
    //move particles from p2 to p1
    p1.clear();
    for(unsigned int i=0; i<p2.size(); i++)
    {
        p1.push_back(p2[i]);
    }

    //move particles from p3 to p2
    p2.clear();
    for(unsigned int i=0; i<p3.size(); i++)
    {
        p2.push_back(p3[i]);
    }

    //move particles from p4 to p3
    p3.clear();
    for(unsigned int i=0; i<p4.size(); i++)
    {
        p3.push_back(p4[i]);
    }
}

void processRealPair()
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

            dhiss[index]->Fill(dphi);
        }
    }

}

void processMixPair()
{
    if(p1.size() == 0) return;

    //Calculate the pseudo-rapidity gap and store particles into right bin
    for(unsigned int i=0; i<p1.size(); i++)
    {
        for(unsigned int j=0; j<p2.size(); j++)
        {
            float deta;
            float dphi;
            int   index;

            //Calculate the pseudo-rapidity gap
            deta = abs(p1[i].eta - p2[j].eta);

            //Get rad of deta>2
            if(deta>2) continue;

            dphi = p1[i].phi - p2[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            index = deta / (2.0 / nog);

            dhisb[index]->Fill(dphi);
            hcount->Fill(index);
        }

        for(unsigned int j=0; j<p3.size(); j++)
        {
            float deta;
            float dphi;
            int   index;

            //Calculate the pseudo-rapidity gap
            deta = abs(p1[i].eta - p3[j].eta);

            //Get rad of deta>2
            if(deta>2) continue;

            dphi = p1[i].phi - p3[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            index = deta / (2.0 / nog);

            dhisb[index]->Fill(dphi);
            hcount->Fill(index);
        }

        for(unsigned int j=0; j<p4.size(); j++)
        {
            float deta;
            float dphi;
            int   index;

            //Calculate the pseudo-rapidity gap
            deta = abs(p1[i].eta - p4[j].eta);

            //Get rad of deta>2
            if(deta>2) continue;

            dphi = p1[i].phi - p4[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            index = deta / (2.0 / nog);

            dhisb[index]->Fill(dphi);
            hcount->Fill(index);
        }
    }
}

 void parseampt()
{
    //Read in data file
    ifstream dataFile;
    dataFile.open("ana/ampt.dat");
    //dataFile.open("/Users/Bill/Desktop/Lab/d-Au_Collision_v2_3_vs_deta_study/10K_Data/ampt_0.dat");

    //Skip the job if not dataFile
    if(!dataFile)  return;

    //count event
    int count = 0;

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

        count = count + 1;

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

        if(count == 1)
        {
            for(unsigned int i=0; i<particles.size(); i++)
            {
                p1.push_back(particles[i]);
            }
        }
        else if(count == 2)
        {
            for(unsigned int i=0; i<particles.size(); i++)
            {
                p2.push_back(particles[i]);
            }
        }
        else if(count == 3)
        {
            for(unsigned int i=0; i<particles.size(); i++)
            {
                p3.push_back(particles[i]);
            }
        }
        else if(count == 4)
        {
            for(unsigned int i=0; i<particles.size(); i++)
            {
                p4.push_back(particles[i]);
            }
        }
        else
        {
            processMixPair();

            move();

            p4.clear();
            for(unsigned int i=0; i<particles.size(); i++)
            {
                p4.push_back(particles[i]);
            }
        }

        processRealPair();
        particles.clear();

        if (!dataFile) break;
    }
}


void t24_new_reference_fitting()
{
    //Define 4 real pair histograms
    for(int i=0; i<nog; i++)
    {
        dhiss.push_back(new TH1F(Form("s_%i",i),  Form("v2: dphi in range of [%f,%f)",i*0.1,(i+1)*0.1), 50, -0.5*TMath::Pi(), 1.5*TMath::Pi()));
    }

    //Define 4 real pair histograms
    for(int i=0; i<nog; i++)
    {
        dhisb.push_back(new TH1F(Form("b_%i",i),  Form("v2: dphi in range of [%f,%f)",i*0.1,(i+1)*0.1), 50, -0.5*TMath::Pi(), 1.5*TMath::Pi()));
    }

    //Define count histgram
    hcount = new TH1F("hcount", "count", 4, 0, 4);

    parseampt();

    //Make a file to store outputs
    TFile *fout = new TFile("out.root","RECREATE");

    //Save graphs
    for(int i=0; i<nog; i++)
    {   
        dhiss[i]->Write();
        dhisb[i]->Write();
    }
    hcount->Write();
    fout->Close();
}





















