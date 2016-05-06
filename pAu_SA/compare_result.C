#include <cmath>
#include <string.h>
#include <iostream.h>

#include "TPavesText.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TPDF.h"
#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGraphErrors.h"

void compare_result()
{
	TString infile1 = "save_fitting_AN_fitting_method1.root";
	TString infile2 = "save_fitting_AN_fitting_method2.root";

	TMultiGraph *Multi_G = new TMultiGraph();
	Multi_G->SetNameTitle("Multi_G","Forward Arm A_{N}");

	TFile *file1 = TFile::Open(infile1);
	TFile *file2 = TFile::Open(infile2);
	
	TGraphErrors *G1_1 = (TGraphErrors *)file1->Get("arm0_beam0_method1");
	TGraphErrors *G1_2 = (TGraphErrors *)file1->Get("arm0_beam1_method1");
	
	TGraphErrors *G2_1 = (TGraphErrors *)file2->Get("arm0_beam0_method2");
        TGraphErrors *G2_2 = (TGraphErrors *)file2->Get("arm0_beam1_method2"); 	
//	TMultiGraph *G1 = (TMultiGraph*)file1->Get("graph0"); 
//	TMultiGraph *G2 = (TMultiGraph*)file2->Get("graph0");
	TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
//	Multi_G->GetXaxis()->SetTitle("p_{T}(MeV)");
//	Multi_G->GetYaxis()->SetTitle("A_{N}");

	TLine *l1 = new TLine(0,0,10.0,0);
	
	TCanvas *c0 = new TCanvas("c0","c0",800,600);
	c0->cd();
	c0->SetGridy();
	Multi_G->Add(G1_1);
	Multi_G->Add(G1_2);
	Multi_G->Add(G2_1);
	Multi_G->Add(G2_2);
	Multi_G->Draw("AP");
	Multi_G->GetXaxis()->SetLimits(0,10.0);
	Multi_G->SetMaximum(0.1);
	Multi_G->SetMinimum(-0.1);
	Multi_G->GetXaxis()->SetTitle("p_{T}(MeV)");
	Multi_G->GetYaxis()->SetTitle("A_{N}");
	leg->AddEntry(G1_1,"North arm Blue beam fit with 3 parameters","lep");
	leg->AddEntry(G1_2,"North arm Yellow beam fit with 3 paratmeters","lep");
	leg->AddEntry(G2_1,"North arm Blue beam fit with 2 parameters","lep");
        leg->AddEntry(G2_2,"North arm Yellow beam fit with 2 paratmeters","lep");
	leg->Draw("SAME");
	l1->Draw("SAME");

	ifstream fin1("AN_method1");
	ifstream fin2("AN_method2");	
	double nb1[2]; 
	double enb1[2];
	double ny1[2]; 
	double eny1[2];
	double sb1[2]; 
	double esb1[2];
        double sy1[2]; 
	double esy1[2];

        double nb2[2]; 
	double enb2[2];
        double ny2[2];
	double eny2[2];
        double sb2[2]; 
	double esb2[2];
        double sy2[2]; 
	double esy2[2];

        double bnb1[2]; 
	double ebnb1[2];
        double bny1[2]; 
	double ebny1[2];
        double bsb1[2]; 
	double ebsb1[2];
        double bsy1[2]; 
	double ebsy1[2];
        
        double bnb2[2];
	double ebnb2[2];
        double bny2[2];
	double  ebny2[2];
        double bsb2[2];
	double  ebsb2[2];
        double bsy2[2];
	double  ebsy2[2];

	int count1=0;
	int count2=0;

	while (fin1>>nb1[count1]>>enb1[count1]>>ny1[count1]>>eny1[count1]>>sb1[count1]>>esb1[count1]>>sy1[count1]>>esy1[count1]>>bnb1[count1]>>ebnb1[count1]>>bny1[count1]>>ebny1[count1]>>bsb1[count1]>>ebsb1[count1]>>bsy1[count1]>>ebsy1[count1]) count1++;

	while (fin2>>nb2[count2]>>enb2[count2]>>ny2[count2]>>eny2[count2]>>sb2[count2]>>esb2[count2]>>sy2[count2]>>esy2[count2]>>bnb2[count2]>>ebnb2[count2]>>bny2[count2]>>ebny2[count2]>>bsb2[count2]>>ebsb2[count2]>>bsy2[count2]>>ebsy2[count2]) count2++;

	for (int i = 0; i< 2; i++)
	{
		cout<<  (nb1[i]-nb2[i])/enb1[i]<<"  "<<
			(ny1[i]-ny2[i])/eny1[i]<<"  "<<
			(sb1[i]-nb2[i])/esb1[i]<<"  "<<
                        (sy1[i]-nb2[i])/esy1[i]<<"  "<<
			(bnb1[i]-bnb2[i])/ebnb1[i]<<"  "<<
                        (bny1[i]-bny2[i])/ebny1[i]<<"  "<<
                        (bsb1[i]-bnb2[i])/ebsb1[i]<<"  "<<
                        (bsy1[i]-bnb2[i])/ebsy1[i]<<"  "<<endl;
	}
}
