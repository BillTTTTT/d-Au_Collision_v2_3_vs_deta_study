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

//Number of bins in histogram
const int nob = 20;

vector<particle> particles;
int n = 0;

vector<TH1F*> dhis2;

int   ct[20];
float v12[20];
float v22[20];
float v32[30];

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

    //Calculate centroid
    float cmx=0;
    float cmy=0;
    float xsum=0;
    float ysum=0;
    float avercos2=0;
    float avercos3=0;
    float aversin2=0;
    float aversin3=0;
    float aver2=0;

    //Calculate centroid
    for (unsigned int i=0; i<particles.size(); i++)
    {
        xsum = xsum + particles[i].x;
        ysum = ysum + particles[i].y;
    }

    cmx = xsum/particles.size();
    cmy = ysum/particles.size();

    //Calculate each useful value of collision particles
    //Store them in vectors
    for (unsigned int i=0; i<particles.size(); i++)
    {
        //Shift to center of mass frame
        particles[i].x = particles[i].x - cmx;
        particles[i].y = particles[i].y - cmy;

        particles[i].rsquare = particles[i].x*particles[i].x + particles[i].y*particles[i].y;
    }

    //Calculate the average values for computing epsilon_2
    for (unsigned int i=0; i<particles.size(); i++)
    {
        avercos2 = avercos2 + particles[i].rsquare * TMath::Cos(2*particles[i].phi);
        aversin2 = aversin2 + particles[i].rsquare * TMath::Sin(2*particles[i].phi);
    }

    avercos2 = avercos2 / particles.size();
    aversin2 = aversin2 / particles.size();

    //Calculate the average values for computing epsilon_3
    for (unsigned int i=0; i<particles.size(); i++)
    {
        avercos3 = avercos3 + particles[i].rsquare * TMath::Cos(3*particles[i].phi);
        aversin3 = aversin3 + particles[i].rsquare * TMath::Sin(3*particles[i].phi);
    }

    avercos3 = avercos3 / particles.size();
    aversin3 = aversin3 / particles.size();

    //Calculate psi2 and psi3
    float psi2;
    float psi3;

    psi2 = (TMath::ATan2(aversin2,avercos2) + TMath::Pi())/2.0;
    psi3 = (TMath::ATan2(aversin3,avercos3) + TMath::Pi())/3.0;

    for (unsigned int i=0; i<particles.size(); i++)
    {
        particles[i].psi2 = psi2;
        particles[i].psi3 = psi3;
    }

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

            index = deta / (2.0 / nob);

            calvs(index, dphi);

            dhis2[index]->Fill(dphi);
        }
    }

}

 void parseampt()
{
    for(int i=0; i<5; i++)
    {
        //Read in data file
        ifstream dataFile;
        dataFile.open(Form("/Users/air/Desktop/Lab/d-Au_Collision_v2_3_vs_deta_study/10K_Data/ampt_%i.dat",i));

        //Print some comment as the file is/isn't successfully opened
        if (!dataFile)
        {
            cout << Form("--> File %i does not exist\n",i+1) << endl << endl;
            return;
        }
        else
        {
            cout << Form("--> Successfully opened file number %i\n",i+1) << endl << endl;
        }

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
                float energy = TMath::Sqrt(pv[0]*pv[0] + pv[1]*pv[1] + pv[2]*pv[2]);

                float pt = TMath::Sqrt(pv[0]*pv[0] + pv[1]*pv[1]);

                if(pt==0) continue;

                //Make Lorentz vector
                TLorentzVector ev(pv[0], pv[1], pv[2], energy);

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
}

void t21_distribution_20bins()
{
    for(int i=0; i<20; i++)
    {
        dhis2.push_back(new TH1F(Form("d_%i",i),  Form("v2: dphi in range of [%f,%f)",i*0.1,(i+1)*0.1), 50, -0.5*TMath::Pi(), 1.5*TMath::Pi()));
    }

    parseampt();
    
    TCanvas *c1 = new TCanvas("c1","c1",2000,1600);
    c1->Divide(5,4);

    TLatex latex;

    for(int i=0; i<20; i++)
    {
        v12[i] = v12[i] / ct[i];
        char strv12[10];
        sprintf(strv12, "%f", v12[i]);
        v22[i] = v22[i] / ct[i];
        char strv22[10];
        sprintf(strv22, "%f", v22[i]);
        v32[i] = v32[i] / ct[i];
        char strv32[10];
        sprintf(strv32, "%f", v32[i]);

        int markx = i%5;
        int marky = i/5;
        float x1 = 0.040 + 0.2 * markx;
        float x2 = 0.075 + 0.2 * markx;
        float y1 = 0.95 - 0.25 * marky;
        float y2 = 0.94 - 0.25 * marky;
        float y3 = 0.93 - 0.25 * marky;
        
        latex.SetTextSize(0.01);
        latex.SetTextAlign(11);
        latex.DrawLatex(x1, y1,"#LTcos(1#Delta#phi)#GT = ");
        latex.DrawLatex(x2, y1,strv12);
        latex.DrawLatex(x1, y2,"#LTcos(2#Delta#phi)#GT = ");
        latex.DrawLatex(x2, y2,strv22);
        latex.DrawLatex(x1, y3,"#LTcos(3#Delta#phi)#GT = ");
        latex.DrawLatex(x2, y3,strv32);
    }

    for(int i=0; i<20; i++)
    {   
        c1->cd(i+1);
        dhis2[i]->Draw();
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(111);
        dhis2[i]->SetXTitle("#Delta#phi");
    }

    c1->Print("AMPT_20dis.pdf");
}





















