//-----------------------------------------------
//Code to graph v3 and v2 vs deta
//for UrQMD Model for d+Au at 200GeV
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
#include "TF1.h"
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
            if(i>=j) continue;

            float deta;
            float dphi;

            //Calculate the pseudo-rapidity gap
            deta = abs(particles[i].eta - particles[j].eta);

            //Get rad of deta>2
            if(deta>2) continue;

            dphi = particles[i].phi - particles[j].phi;

            if(dphi > TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            t2->Fill(deta,TMath::Cos(2 * dphi));
            t3->Fill(deta,TMath::Cos(3 * dphi));
            t220->Fill(deta,TMath::Cos(2 * dphi));
            t320->Fill(deta,TMath::Cos(3 * dphi));
        }
    }

}

void parseurqmd()
{
	//Make a tokens as string format to store each line read in from data file
	string linestr;
	vector<string> tokens;

	//Since there are a set of data files storing the data, we need a for loop to read and calculate them one by one
	for(int nfile=0; nfile<100; nfile++)
	{
		//---------------------------------------------------------------------
		//Read in file
		ifstream dataFile;

		//There are 4 folders named as "UrQMD_20GeV", "UrQMD_39GeV", "UrQMD_62GeV", "UrQMD_200GeV"
		//In each folder, data file is named as "test.f20_0", "test.f20_1", "test.f20_2", and so on
		dataFile.open(Form("/gpfs/mnt/gpfs02/phenix/plhf3/pengqi/UrQMD_d+Au/UrQMD_200GeV/test.f20_%i",nfile));

		//Print some comment as the file is/isn't successfully opened
		if (!dataFile)
		{
			cout << Form("--> File test.f20_%i does not exist\n",nfile+1) << endl << endl;
			continue;
		}
		else
		{
			cout << Form("--> Successfully opened file number %i\n",nfile+1) << endl << endl;
		}

		//In this while loop, program will read the data file line by line
		while(dataFile)
		{
			//Clear the tokens before get any new line
			tokens.clear();
			//Get a line from data file
			std::getline(dataFile,linestr);
			//Breaking up tokens into individual words
			istringstream iss(linestr);
			//Put breaked up words back to the tokens
			copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

	  		//Ignore file header starting with #
			if(tokens[0] == "#") continue;

	  		//Find new event header
			if(atoi(tokens[0].c_str()) == 0 && atoi(tokens[1].c_str()) == NUCL)
			{
	    		//Skip the following lines
				for(int i=0; i<NUCL; i++)
				{
					std::getline(dataFile,linestr);
				}

				continue;
			}

	 		//Find collision header in the current event for 2 -> many scatterings
			if(tokens[0] == "2")
			{
				//Get the number of particals produced
				int numcollproducts = atoi(tokens[1].c_str());

				//Skip over all of the lines corresponding to collision products
				for(int i=0; i<numcollproducts+2; i++)
				{
					std::getline(dataFile,linestr);	
				}
				continue;
			}

			//Identify final-state particles and store them
			//The number of final state particles for the event (which has at least one collision) should greater than initial.
			//In this case, it should be greater than 199
			if(atoi(tokens[0].c_str()) > 199 && tokens[1] == "0" && tokens.size() == 2)
			{
				//Do calculation for each particle
				for(int i=0; i<atoi(tokens[0].c_str()); i++)
				{
					vector<string> finalstate_tokens;
					std::getline(dataFile,linestr);
					istringstream iss4(linestr);
					copy(istream_iterator<string>(iss4), istream_iterator<string>(), back_inserter(finalstate_tokens));				

					particle pf;

					//Get Lorentz vector
					TLorentzVector ev(atof(finalstate_tokens[3].c_str()), atof(finalstate_tokens[4].c_str()), atof(finalstate_tokens[5].c_str()), atof(finalstate_tokens[6].c_str()));
					
					//Skip the particle that we are not interested in
					//+-211 is pions,  +-321 is kaons, +-2212 is protons
					if(abs(atof(finalstate_tokens[1].c_str())) != 211 && abs(atof(finalstate_tokens[1].c_str())) != 321 && abs(atof(finalstate_tokens[1].c_str())) != 2212) continue;

					//Get pT, phi, pseudorapidity, particle id, px, py, and pz. Store them into pf.
					pf.pT  = ev.Pt();
					pf.phi = ev.Phi();
					pf.eta = ev.Eta();
					pf.id  = atof(finalstate_tokens[1].c_str());
					pf.px  = atof(finalstate_tokens[3].c_str());
					pf.py  = atof(finalstate_tokens[4].c_str());
					pf.pz  = atof(finalstate_tokens[5].c_str());

					if(pf.pT <= 0.2) continue;
                	if(abs(pf.eta) >= 1) continue;

                	particles.push_back(p);
				}

				processEvent();
				particles.clear();
			}

			//Skip the final state of no-collision event
			if(atoi(tokens[0].c_str()) == NUCL)
			{
				for(int i=0; i<NUCL; i++)
				{
					std::getline(dataFile,linestr);
				}

				continue;
			}

			//Skip "0    0" line at the end of each event
			if(tokens[0] == "0" && tokens[1] == "0") continue;

			if (!dataFile) break;
		}

		dataFile.close();
	}
	
}

void t20_UrQMD_v2_deta()
{
    t2 = new TProfile("t2","t2",100,-0.01,1.99,-1,+1);
    t3 = new TProfile("t3","t3",100,-0.01,1.99,-1,+1);
    t220 = new TProfile("t220","t220",20,-0.01,1.99,-1,+1);
    t320 = new TProfile("t320","t320",20,-0.01,1.99,-1,+1);
    TF1 *fit1 = new TF1("fit1","pol0",0,2); 
    TF1 *fit2 = new TF1("fit2","pol0",0,2); 

    parseurqmd();


    TCanvas *c1 = new TCanvas("c1","v2",700,600);
    t220->Fit(fit1,"R");
    t2->SetMarkerStyle(24);
    t2->Draw("p,e");
    t220->SetMarkerStyle(20);
    t220->Draw("p,e,same");

    TCanvas *c2 = new TCanvas("c2","v3",700,600);
    t320->Fit(fit2,"R");
    t3->SetMarkerStyle(24);
    t3->Draw("p,e");
    t320->SetMarkerStyle(20);
    t320->Draw("p,e,same");

    c1->Print("v2_UrQMD.pdf");
    c2->Print("v3_UrQMD.pdf");
}





























