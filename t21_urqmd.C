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
const int nob = 4;

const int nog = 4;

vector<particle> particles;
int n = 0;

vector<TH1F*> dhis2;
vector<TF1*>  fit;

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

void parseurqmd()
{
    //Make a tokens as string format to store each line read in from data file
    string linestr;
    vector<string> tokens;

    //Since there are a set of data files storing the data, we need a for loop to read and calculate them one by one
    for(int nfile=0; nfile<200; nfile++)
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

                    particles.push_back(pf);
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

void t21_urqmd()
{
    //Define 4 histograms
    for(int i=0; i<nog; i++)
    {
        dhis2.push_back(new TH1F(Form("d_%i",i),  Form("v2: dphi in range of [%f,%f)",i*0.1,(i+1)*0.1), 50, -0.5*TMath::Pi(), 1.5*TMath::Pi()));
    }

    //Define the fit graphs
    for(int i=0; i<nog; i++)
    {
        fit.push_back(new TF1(Form("f_%i",i),"[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi()));
    }

    parseurqmd();
    
    TCanvas *c1 = new TCanvas("c1","c1",2000,1600);
    c1->Divide(2,2);

    TLatex latex;

    //Show the number of v12, v22, v32
    for(int i=0; i<nog; i++)
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

        int markx = i%2;
        int marky = i/2;
        float x1 = 0.08 + 0.5 * markx;
        float x2 = 0.18 + 0.5 * markx;
        float y1 = 0.90 - 0.5 * marky;
        float y2 = 0.87 - 0.5 * marky;
        float y3 = 0.84 - 0.5 * marky;
        
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
        dhis2[i]->Draw();

        //Fit the graphs
        dhis2[i]->Fit(fit[i],"R");

        //Get the parameters from the function
        fit[i] = (TF1*)dhis2[i]->GetFunction(Form("f_%i",i));
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

        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(111);
        dhis2[i]->SetXTitle("#Delta#phi");

        TLegend *leg = new TLegend(0.65,0.15,0.85,0.35);
        leg->AddEntry(fit[i],"fit line","L");
        leg->AddEntry(f1,"cos(1#Delta#phi)","L");
        leg->AddEntry(f2,"cos(2#Delta#phi)","L");
        leg->AddEntry(f3,"cos(3#Delta#phi)","L");
        leg->Draw("same");
    }

    c1->Print("AMPT_4dis_20GeV.pdf");
}





















