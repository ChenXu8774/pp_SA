#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>
#include <string.h>

#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGraphErrors.h"
void result_lessbin( int check = 6)
{

	TString infile="432639_hist.root";
//	TString infile="diMuon_his.lst.root";
	gStyle->SetOptFit();
	
	TFile *file = TFile::Open(infile);
	if (check == 0)
	{
		TH2F *_sig_s_uu = (TH2F*)file->Get("_sah_pT_os_sig_n_PU");
		TH2F *_sig_s_du = (TH2F*)file->Get("_sah_pT_os_sig_n_MU");
//		TH2F *_sig_s_dd = (TH2F*)file->Get("_sah_pT_os_lsb_n_MM");
	
		ofstream fout("sig_n_countsbin.txt");

		for (int i = 0; i<16; i++)
		{
			fout<<//"event "<<setw(3)<<i<<" : " 
			_sig_s_uu->GetBinContent(1,i+1)<<"      "<<
			_sig_s_du->GetBinContent(1,i+1)
//			<<" uu "<<_sig_s_uu->GetBinContent(1,i+1)
//			<<" ud "<<_sig_s_ud->GetBinContent(1,i+1)
//			<<" du "<<_sig_s_du->GetBinContent(1,i+1)
//			<<" dd "<<_sig_s_dd->GetBinContent(1,i+1)
			<<endl;
		}
	}

	if (check == 1)
	{
		TH2 *_sig_n_al_blue = (TH2F*)file->Get("_sah_pT_os_sig_n_A_L_Blue_SINGLE");
		TH2 *_sig_s_al_blue = (TH2F*)file->Get("_sah_pT_os_sig_s_A_L_Blue_SINGLE");

		TH2 *_lsb_n_al_blue = (TH2F*)file->Get("_sah_pT_os_lsb_n_A_L_Blue_SINGLE");
                TH2 *_lsb_s_al_blue = (TH2F*)file->Get("_sah_pT_os_lsb_s_A_L_Blue_SINGLE");

		ofstream fout("al.txt");
		for (int i = 0; i<16; i++)
			fout<<"bin "<<setw(3)<<i<<" : "
			<<setw(15)<<_sig_n_al_blue->GetBinContent(1,i+1)
//			<<setw(15)<<_sig_s_al_blue->GetBinContent(1,i+1)

			<<setw(15)<<_lsb_n_al_blue->GetBinContent(1,i+1)
//                        <<setw(15)<<_lsb_s_al_blue->GetBinContent(1,i+1)
			

			<<endl;
	}

	if (check == 6)
	{
		TH2 *_sig_n_al_blue = (TH2F*)file->Get("_sah_xF_os_sig_s_A_L_Blue_SINGLE");

                TH2 *_lsb_n_al_blue = (TH2F*)file->Get("_sah_xF_os_lsb_s_A_L_Blue_SINGLE");
		ofstream fout("xf_al.txt");
		for (int i = 0; i<16; i++)
			fout<<"bin "<<setw(3)<<i<<" : "
			<<setw(15)<<_sig_n_al_blue->GetBinContent(1,i+1)
			<<setw(15)<<_lsb_n_al_blue->GetBinContent(1,i+1)
			<<endl;	
	}
	if (check == 2)
	{
		int fill;
		double polb, poly;	
		double p1,p2, total;
		ifstream fin1("fill_pol.txt");
		ofstream fout("fill_lumi.txt");
		while (fin1>>fill>>polb>>poly)
		{
			char *n1 = Form("RelLumi_FILL_0%d_SUM",fill);
//			char *n2 = Form("RelLumi_FILL_0%d_SUM",fill);
//			char *n3 = Form("RelLumi_FILL_0%d_Pol_Ylue",fill);
			TH1D *_Rel_SUM = (TH1D*)file->Get(n1);
//			TH1D *_Rel_Pol_B=(TH1D*)file->Get(n2);
//			TH1D *_Rel_Pol_Y=(TH1D*)file->Get(n3);
// 			p1 = _Rel_Pol_B->GetBinContent(1)/_Rel_SUM->GetBinContent(1);
			total +=_Rel_SUM->GetBinContent(1);
			fout <<" fill " << fill<< "       "<<"Lumi:   "<< _Rel_SUM->GetBinContent(1) <<endl;
		}
		cout<<"total lumi: " <<total<<endl;
	}
	
	if (check == 3)
        {
		int run;
		double lumi;
		ifstream fin("spinQA.txt");
		ofstream fout("lumi_run.txt");
	
		while (fin>>run)
		{
			char *n = Form("RelLumi_RUN_0000%d_SUM",run);
			TH1D *h = (TH1D*)file->Get(n);
			if (h==NULL) continue;
			fout<<"run "<<run<<"  lumi:   "<<h->GetBinContent(1) <<endl;
		}
	}

	if (check == 4)
	{
		ofstream fout("compare_A_L.txt");
		double P = 0.5928415;
		double R = 1.0060345;
		double N_P = 0;
		double N_M = 0;
		TH1D *h1 = (TH1D*)file->Get("_sah_pT_os_sig_n_PU");
		TH1D *h2 = (TH1D*)file->Get("_sah_pT_os_sig_n_MU");
		TH1D *h3 = (TH1D*)file->Get("_sah_pT_os_sig_n_A_L_Blue_SINGLE");
		for(int i =0; i < 16; i++)
		fout << (1./P*(h1->GetBinContent(1,i+1)-R*h2->GetBinContent(1,i+1))/(h1->GetBinContent(1,i+1)+R*h2->GetBinContent(1,i+1)))<<"	"<<h3->GetBinContent(1,i+1)<<endl;
		
	}
	
	if (check == 5)
	{
		ofstream fout("xf_n_hits_sig.txt");
		TH1D *h1 = (TH1D*)file->Get("_sah_xF_os_sig_s_PU");
		TH1D *h2 = (TH1D*)file->Get("_sah_xF_os_sig_s_MU");
		for (int i = 0; i <16; i ++)
		fout<<"event:"<<setw(3)<<i<<" : "
		<<" UP "<<h1->GetBinContent(1,i+1)
		<<" DN "<<h2->GetBinContent(1,i+1)
		<<endl;
	}
}
