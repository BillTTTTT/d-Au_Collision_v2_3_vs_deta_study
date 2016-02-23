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
    int   mark;
};

//Number of nucleons in system
// --> p+Au = 198
// --> d+Au = 199
// --> d+Pb = 210
// --> p+Pb = 209
const int NUCL = 199;

vector<particle> particles;
int n = 0;

vector<float> v2bins;
vector<float> v3bins;

TProfile* t2;
TProfile* t3;
TProfile* t220;
TProfile* t320;

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
        for(unsigned int j=1; j<particles.size(); j++)
        {
            if(i==j) continue;

            float deta;
            float dphi;

            //Calculate the pseudo-rapidity gap
            deta = abs(particles[i].eta - particles[j].eta);

            //Get rad of deta>2
            if(deta>2) continue;

            dphi = particles[i].phi - particles[j].phi;

            t2->Fill(deta,TMath::Cos(2 * dphi));
            t3->Fill(deta,TMath::Cos(3 * dphi));
            t220->Fill(deta,TMath::Cos(2 * dphi));
            t320->Fill(deta,TMath::Cos(3 * dphi));
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
      		printf(Form("--> File %i does not exist\n",i+1));
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
                p.mark= 0;

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

void t19_ampt_v3_deta()
{
    t2 = new TProfile("t2","t2",100,-0.01,1.99,-1,+1);
    t3 = new TProfile("t3","t3",100,-0.01,1.99,-1,+1);
    t220 = new TProfile("t220","t220",20,-0.01,1.99,-1,+1);
    t320 = new TProfile("t320","t320",20,-0.01,1.99,-1,+1);

    parseampt();

    TCanvas *c1 = new TCanvas("c1","v2",700,600);
    t2->SetMarkerStyle(24);
    t2->Draw("p,e");
    t220->SetMarkerStyle(20);
    t220->Draw("p,e,same");
    TCanvas *c2 = new TCanvas("c2","v3",700,600);
    t2->SetMarkerStyle(24);
    t2->Draw("p,e");
    t320->SetMarkerStyle(20);
    t320->Draw("p,e,same");

    //c1->Print("v2.pdf");
    //c2->Print("v3.pdf");
}






























