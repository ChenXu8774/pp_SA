#include <cmath>
#include <string.h>
#include <iostream.h>
#include <fstream>

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

void result(int method = 2)
{
	
	TFile *f_out = new TFile(Form("save_fitting_AN_fitting_method%d.root",method),"RECREATE");
	TString infile ="422070_hist.root";
	double r[2][2];//[arm][n_pt_bin]

	double average_pt[2][2][2];//[arm][beam][pt_bin]
	average_pt[0][0][0] = 1.081;
	average_pt[0][1][0] = 1.081;
	average_pt[0][0][1] = 2.98;
	average_pt[0][1][1] = 2.98;

	average_pt[1][0][0] = 1.068;
	average_pt[1][1][0] = 1.068;
        average_pt[1][0][1] = 2.964;
	average_pt[1][1][1] = 2.964;

	average_pt[0][0][0] = average_pt[0][0][0]+0.5*(method-1);
        average_pt[0][1][0] = average_pt[0][1][0]+0.5*(method-1)+0.25;
	average_pt[0][0][1] = average_pt[0][0][1]+0.5*(method-1);
        average_pt[0][1][1] = average_pt[0][1][1]+0.5*(method-1)+0.25;


        average_pt[1][0][0] = average_pt[1][0][0]+0.5*(method-1);
        average_pt[1][1][0] = average_pt[1][1][0]+0.5*(method-1)+0.25;
	average_pt[1][0][1] = average_pt[1][0][1]+0.5*(method-1);
        average_pt[1][1][1] = average_pt[1][1][1]+0.5*(method-1)+0.25;


	r[0][0] = 0.141;
	r[0][1] = 0.106;

	r[1][0] = 0.165;
	r[1][1] = 0.117;

	double AN_inc[2][2][2]; //[arm][beam][pt_bin]	
	double AN_lsb[2][2][2]; //[arm][beam][pt_bin]  
	double AN_phy[2][2][2]; //[arm][beam][pt_bin] 
	
	double eAN_inc[2][2][2]; //[arm][beam][pt_bin]   
        double eAN_lsb[2][2][2]; //[arm][beam][pt_bin]  
        double eAN_phy[2][2][2]; //[arm][beam][pt_bin]
	gStyle->SetOptFit();
	
	TFile *file = TFile::Open(infile);
	//signal hitogram
	TH2F *sig_n_A_L_Blue = (TH2F*)file->Get("_sah_pT_os_sig_n_A_L_Blue");
	TH2F *sig_n_A_L_Yellow = (TH2F*)file->Get("_sah_pT_os_sig_n_A_L_Yellow");
	TH2F *sig_s_A_L_Blue = (TH2F*)file->Get("_sah_pT_os_sig_s_A_L_Blue");
        TH2F *sig_s_A_L_Yellow = (TH2F*)file->Get("_sah_pT_os_sig_s_A_L_Yellow");

	//background histogram
	TH2F *bgr_n_A_L_Blue = (TH2F*)file->Get("_sah_pT_os_lsb_n_A_L_Blue");
        TH2F *bgr_n_A_L_Yellow = (TH2F*)file->Get("_sah_pT_os_lsb_n_A_L_Yellow");
        TH2F *bgr_s_A_L_Blue = (TH2F*)file->Get("_sah_pT_os_lsb_s_A_L_Blue");
        TH2F *bgr_s_A_L_Yellow = (TH2F*)file->Get("_sah_pT_os_lsb_s_A_L_Yellow");

	const Int_t n_pt_bin = 2;
	const Int_t n_phi_bin = 16;
// sig_n_A_L_Blue->GetNbinsX();

	//signal variables
	float al_b_n[n_pt_bin][n_phi_bin];
	float al_error_b_n[n_pt_bin][n_phi_bin];
	float phi_b_n[n_phi_bin];

	float al_y_n[n_pt_bin][n_phi_bin];
        float al_error_y_n[n_pt_bin][n_phi_bin];
        float phi_y_n[n_phi_bin];

	float al_b_s[n_pt_bin][n_phi_bin];
        float al_error_b_s[n_pt_bin][n_phi_bin];
        float phi_b_s[n_phi_bin];

        float al_y_s[n_pt_bin][n_phi_bin];
        float al_error_y_s[n_pt_bin][n_phi_bin];
        float phi_y_s[n_phi_bin];

	//background variables
	float bal_b_n[n_pt_bin][n_phi_bin];
        float bal_error_b_n[n_pt_bin][n_phi_bin];
        float bphi_b_n[n_phi_bin];

        float bal_y_n[n_pt_bin][n_phi_bin];
        float bal_error_y_n[n_pt_bin][n_phi_bin];
        float bphi_y_n[n_phi_bin];

        float bal_b_s[n_pt_bin][n_phi_bin];
        float bal_error_b_s[n_pt_bin][n_phi_bin];
        float bphi_b_s[n_phi_bin];

        float bal_y_s[n_pt_bin][n_phi_bin];
        float bal_error_y_s[n_pt_bin][n_phi_bin];
        float bphi_y_s[n_phi_bin];

	for (int j = 0; j<n_pt_bin; j++)
		for (int i = 0; i<n_phi_bin; i++)
		{

		//give value to signal variable
			al_b_n[j][i] = sig_n_A_L_Blue->GetBinContent(j+1,i+1);
			al_error_b_n[j][i] = sig_n_A_L_Blue->GetBinError(j+1,i+1);
			phi_b_n[i] = sig_n_A_L_Blue->GetYaxis()->GetBinCenter(i+1);

			al_y_n[j][i] = sig_n_A_L_Yellow->GetBinContent(j+1,i+1);
                	al_error_y_n[j][i] = sig_n_A_L_Yellow->GetBinError(j+1,i+1);
                	phi_y_n[i] = sig_n_A_L_Yellow->GetYaxis()->GetBinCenter(i+1);		

			al_b_s[j][i] = sig_s_A_L_Blue->GetBinContent(j+1,i+1);
                	al_error_b_s[j][i] = sig_s_A_L_Blue->GetBinError(j+1,i+1);
                	phi_b_s[i] = sig_s_A_L_Blue->GetYaxis()->GetBinCenter(i+1);
	
        	        al_y_s[j][i] = sig_s_A_L_Yellow->GetBinContent(j+1,i+1);
               	 	al_error_y_s[j][i] = sig_s_A_L_Yellow->GetBinError(j+1,i+1);
                	phi_y_s[i] = sig_s_A_L_Yellow->GetYaxis()->GetBinCenter(i+1); 
	
			//give value to background variable
			bal_b_n[j][i] = bgr_n_A_L_Blue->GetBinContent(j+1,i+1);
                	bal_error_b_n[j][i] = bgr_n_A_L_Blue->GetBinError(j+1,i+1);
                	bphi_b_n[i] = bgr_n_A_L_Blue->GetYaxis()->GetBinCenter(i+1);
	
        	        bal_y_n[j][i] = bgr_n_A_L_Yellow->GetBinContent(j+1,i+1);
               	 	bal_error_y_n[j][i] = bgr_n_A_L_Yellow->GetBinError(j+1,i+1);
                	bphi_y_n[i] = bgr_n_A_L_Yellow->GetYaxis()->GetBinCenter(i+1);

                	bal_b_s[j][i] = bgr_s_A_L_Blue->GetBinContent(j+1,i+1);
                	bal_error_b_s[j][i] = bgr_s_A_L_Blue->GetBinError(j+1,i+1);
                	bphi_b_s[i] = bgr_s_A_L_Blue->GetYaxis()->GetBinCenter(i+1);

                	bal_y_s[j][i] = bgr_s_A_L_Yellow->GetBinContent(j+1,i+1);
                	bal_error_y_s[j][i] = bgr_s_A_L_Yellow->GetBinError(j+1,i+1);
                	bphi_y_s[i] = bgr_s_A_L_Yellow->GetYaxis()->GetBinCenter(i+1);
		}
	
	TGraphErrors *phi_sig_n_A_L_Blue[n_pt_bin];
	TGraphErrors *phi_sig_n_A_L_Yellow[n_pt_bin];
	TGraphErrors *phi_sig_s_A_L_Blue[n_pt_bin];
	TGraphErrors *phi_sig_s_A_L_Yellow[n_pt_bin];
	TGraphErrors *phi_bgr_n_A_L_Blue[n_pt_bin];
	TGraphErrors *phi_bgr_n_A_L_Yellow[n_pt_bin];
	TGraphErrors *phi_bgr_s_A_L_Blue[n_pt_bin];
	TGraphErrors *phi_bgr_s_A_L_Yellow[n_pt_bin];


	for (int i = 0; i < n_pt_bin; i++)
	{
		//Graph for signal
		phi_sig_n_A_L_Blue[i] = new TGraphErrors(n_phi_bin,phi_b_n,al_b_n[i],0,al_error_b_n[i]); 
		phi_sig_n_A_L_Yellow[i] = new TGraphErrors(n_phi_bin,phi_y_n,al_y_n[i],0,al_error_y_n[i]);
		phi_sig_s_A_L_Blue[i] = new TGraphErrors(n_phi_bin,phi_b_s,al_b_s[i],0,al_error_b_s[i]);
        	phi_sig_s_A_L_Yellow[i] = new TGraphErrors(n_phi_bin,phi_y_s,al_y_s[i],0,al_error_y_s[i]);
		
		//Graph for background 
		phi_bgr_n_A_L_Blue[i] = new TGraphErrors(n_phi_bin,bphi_b_n,bal_b_n[i],0,bal_error_b_n[i]);
        	phi_bgr_n_A_L_Yellow[i] = new TGraphErrors(n_phi_bin,bphi_y_n,bal_y_n[i],0,bal_error_y_n[i]);
        	phi_bgr_s_A_L_Blue[i] = new TGraphErrors(n_phi_bin,bphi_b_s,bal_b_s[i],0,bal_error_b_s[i]);
        	phi_bgr_s_A_L_Yellow[i] = new TGraphErrors(n_phi_bin,bphi_y_s,bal_y_s[i],0,bal_error_y_s[i]);

		phi_sig_n_A_L_Blue[i]->SetTitle(Form("North_incl_Blue_pTBin_%d",i));
		phi_sig_n_A_L_Blue[i]->GetXaxis()->SetTitle("phi");
		phi_sig_n_A_L_Blue[i]->GetYaxis()->SetTitle("AN");
		
		phi_sig_n_A_L_Yellow[i]->SetTitle(Form("North_incl_Yellow_pTBin_%d",i));
        	phi_sig_n_A_L_Yellow[i]->GetXaxis()->SetTitle("phi");
        	phi_sig_n_A_L_Yellow[i]->GetYaxis()->SetTitle("AN");	
	
		phi_sig_s_A_L_Blue[i]->SetTitle(Form("South_incl_Blue_pTBin_%d",i));
        	phi_sig_s_A_L_Blue[i]->GetXaxis()->SetTitle("phi");
       	 	phi_sig_s_A_L_Blue[i]->GetYaxis()->SetTitle("AN");
        
        	phi_sig_s_A_L_Yellow[i]->SetTitle(Form("South_incl_Yellow_pTBin_%d",i));
        	phi_sig_s_A_L_Yellow[i]->GetXaxis()->SetTitle("phi");
        	phi_sig_s_A_L_Yellow[i]->GetYaxis()->SetTitle("AN");

        	phi_bgr_n_A_L_Blue[i]->SetTitle(Form("North_bgr_Blue_pTBin_%d",i));
        	phi_bgr_n_A_L_Blue[i]->GetXaxis()->SetTitle("phi");
        	phi_bgr_n_A_L_Blue[i]->GetYaxis()->SetTitle("AN");

        	phi_bgr_n_A_L_Yellow[i]->SetTitle(Form("North_bgr_Yellow_pTBin_%d",i));
        	phi_bgr_n_A_L_Yellow[i]->GetXaxis()->SetTitle("phi");
        	phi_bgr_n_A_L_Yellow[i]->GetYaxis()->SetTitle("AN");

        	phi_bgr_s_A_L_Blue[i]->SetTitle(Form("South_bgr_Blue_pTBin_%d",i));
        	phi_bgr_s_A_L_Blue[i]->GetXaxis()->SetTitle("phi");
        	phi_bgr_s_A_L_Blue[i]->GetYaxis()->SetTitle("AN");

        	phi_bgr_s_A_L_Yellow[i]->SetTitle(Form("South_bgr_Yellow_pTBin_%d",i));
        	phi_bgr_s_A_L_Yellow[i]->GetXaxis()->SetTitle("phi");
        	phi_bgr_s_A_L_Yellow[i]->GetYaxis()->SetTitle("AN");
	}	//end for phi_Al TGraphErrors construction

	TF1 *f_nb[n_pt_bin];
	TF1 *f_ny[n_pt_bin];
	TF1 *f_sb[n_pt_bin];
	TF1 *f_sy[n_pt_bin];
	TF1 *f_bnb[n_pt_bin];
	TF1 *f_bny[n_pt_bin];
	TF1 *f_bsb[n_pt_bin];
	TF1 *f_bsy[n_pt_bin];

	for(int i =0; i<n_pt_bin;i++)
	{
		if (method == 1)
		{
			f_nb[i] = new TF1("f_nb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-3.1415926,3.1415926);
			f_ny[i] = new TF1("f_ny","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-3.1415926,3.1415926);
			f_sb[i] = new TF1("f_sb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-3.1415926,3.1415926);
			f_sy[i] = new TF1("f_sy","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-3.1415926,3.1415926);

			f_bnb[i] = new TF1("f_bnb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-3.1415926,3.1415926);
        		f_bny[i] = new TF1("f_bny","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-3.1415926,3.1415926);
        		f_bsb[i] = new TF1("f_bsb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-3.1415926,3.1415926);
        		f_bsy[i] = new TF1("f_bsy","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-3.1415926,3.1415926);
		}

		if (method == 2)
		{
			f_nb[i] = new TF1("f_nb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        		f_ny[i] = new TF1("f_ny","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        		f_sb[i] = new TF1("f_sb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        		f_sy[i] = new TF1("f_sy","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
	

			f_bnb[i] = new TF1("f_bnb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        		f_bny[i] = new TF1("f_bny","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        		f_bsb[i] = new TF1("f_bsb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        		f_bsy[i] = new TF1("f_bsy","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
		}
	}   //end for tf construction

	ofstream fout1(Form("AN_method%d",method));
	
	//fit 
	for (int i =0; i< n_pt_bin; i++)
	{
		f_nb[i]->SetParName(0,"AN");
		f_ny[i]->SetParName(0,"AN");
		f_sb[i]->SetParName(0,"AN");
		f_sy[i]->SetParName(0,"AN");
		f_bnb[i]->SetParName(0,"AN");
		f_bny[i]->SetParName(0,"AN");
		f_bsb[i]->SetParName(0,"AN");
		f_bsy[i]->SetParName(0,"AN");

		phi_sig_n_A_L_Blue[i]->Fit(f_nb[i]);
		phi_sig_n_A_L_Yellow[i]->Fit(f_ny[i]);
        	phi_sig_s_A_L_Blue[i]->Fit(f_sb[i]);
        	phi_sig_s_A_L_Yellow[i]->Fit(f_sy[i]);

		phi_bgr_n_A_L_Blue[i]->Fit(f_bnb[i]);
        	phi_bgr_n_A_L_Yellow[i]->Fit(f_bny[i]);
        	phi_bgr_s_A_L_Blue[i]->Fit(f_bsb[i]);
        	phi_bgr_s_A_L_Yellow[i]->Fit(f_bsy[i]);
		fout1<<f_nb[i]->GetParameter(0)<<"  "<<f_nb[i]->GetParError(0)<<"  "<<
	        f_ny[i]->GetParameter(0)<<"  "<<f_ny[i]->GetParError(0)<<"  "<<
	        f_sb[i]->GetParameter(0)<<"  "<<f_sb[i]->GetParError(0)<<"  "<<
	        f_sy[i]->GetParameter(0)<<"  "<<f_sy[i]->GetParError(0)<<"  "<<
	        f_bnb[i]->GetParameter(0)<<"  "<<f_bnb[i]->GetParError(0)<<"  "<<
                f_bny[i]->GetParameter(0)<<"  "<<f_bny[i]->GetParError(0)<<"  "<<
                f_bsb[i]->GetParameter(0)<<"  "<<f_bsb[i]->GetParError(0)<<"  "<<
                f_bsy[i]->GetParameter(0)<<"  "<<f_bsy[i]->GetParError(0)<<endl;
	}   //end for fit


	for (int i = 0; i< n_pt_bin; i++)
	{
		AN_inc[0][0][i] = f_nb[i]->GetParameter(0);
		eAN_inc[0][0][i]= f_nb[i]->GetParError(0);

		AN_lsb[0][0][i] = f_bnb[i]->GetParameter(0);
		eAN_lsb[0][0][i]= f_bnb[i]->GetParError(0);		 

		AN_inc[0][1][i] = f_ny[i]->GetParameter(0);
                eAN_inc[0][1][i]= f_ny[i]->GetParError(0);

                AN_lsb[0][1][i] = f_bny[i]->GetParameter(0);
                eAN_lsb[0][1][i]= f_bny[i]->GetParError(0);
		
		AN_inc[1][0][i] = f_sb[i]->GetParameter(0);
                eAN_inc[1][0][i]= f_sb[i]->GetParError(0);

                AN_lsb[1][0][i] = f_bsb[i]->GetParameter(0);
                eAN_lsb[1][0][i]= f_bsb[i]->GetParError(0);

		AN_inc[1][1][i] = f_sy[i]->GetParameter(0);
                eAN_inc[1][1][i]= f_sy[i]->GetParError(0);

                AN_lsb[1][1][i] = f_bsy[i]->GetParameter(0);
                eAN_lsb[1][1][i]= f_bsy[i]->GetParError(0);
	}

	ofstream fout(Form("an_check_method%d.txt",method));
	for (int i = 0; i < 2; i++)
		for(int j = 0; j <2; j++)
			for(int k = 0; k <n_pt_bin; k++)
			{
				AN_phy[i][j][k] = (AN_inc[i][j][k] - r[i][k]*AN_lsb[i][j][k])/(1. - r[i][k]);
				eAN_phy[i][j][k]= sqrt(eAN_inc[i][j][k]*eAN_inc[i][j][k] + 
				r[i][k]*r[i][k]*eAN_lsb[i][j][k]*eAN_lsb[i][j][k]) / (1. - r[i][k]);
				fout <<"arm "<<i <<"  beam  "<<j<<"  pt_bin "<<k<<"  inc	"<<AN_inc[i][j][k]<<"    lsb:"<<AN_lsb[i][j][k] <<"	AN_phy: " << AN_phy[i][j][k]<<"	eAN_phy: "<<eAN_phy[i][j][k]<<endl;		
				
			}
	TMultiGraph *mg[2];
	for (int i =0; i<2; i++) 
		mg[i] = new TMultiGraph(Form("graph%d",i),Form("graph%d",i));		

//	TGraphAsymmErrors *AN_vs_pT[2][2];//[arm][beam]
	TGraphErrors *AN_vs_pT[2][2];
	TLine *l1 = new TLine(0,0,10.0,0);

	char title[2][2];
	title[0][0] = "A_N_vs_pt_North_Blue";
	title[0][1] = "A_N_vs_pt_North_Yellow";
	title[1][0] = "A_N_vs_pt_South_Blue";
        title[1][1] = "A_N_vs_pt_South_Yellow";	
	
	f_out->cd();
	double pt_edge[3]={0.,2.,10.};
	double ex_high[2][2][2];
	double ex_low[2][2][2];
	for (int i = 0; i < 2; i++) //arm
                for(int j = 0; j <2; j++)  //beam
		for(int k = 0; k <2; k++)
		{
			ex_high[i][j][k] = pt_edge[k+1] - average_pt[i][j][k];
			ex_low[i][j][k]  = average_pt[i][j][k] - pt_edge[k];
		}
	for (int i = 0; i < 2; i++) //arm
                for(int j = 0; j <2; j++)  //beam
		{
//			AN_vs_pT[i][j] = new TGraphAsymmErrors(n_pt_bin,average_pt[i],AN_phy[i][j],ex_low[i][j],ex_high[i][j],eAN_phy[i][j],eAN_phy[i][j]);
			AN_vs_pT[i][j] = new TGraphErrors(n_pt_bin,average_pt[i][j],AN_phy[i][j],0,eAN_phy[i][j]);
			AN_vs_pT[i][j]->SetLineWidth(2);
		//	AN_vs_pT[i][j]->SetTitle(title[i][j]);
			AN_vs_pT[i][j]->SetNameTitle(Form("arm%d_beam%d_method%d",i,j,method),Form("arm%d_beam%d_method%d",i,j,method));
			AN_vs_pT[i][j]->GetXaxis()->SetTitle("p_{T}(MeV)");
			AN_vs_pT[i][j]->GetXaxis()->SetLimits(0,10.0);

			AN_vs_pT[i][j]->GetYaxis()->SetTitle("A_{N}");
			AN_vs_pT[i][j]->SetMaximum(0.15);
			AN_vs_pT[i][j]->SetMinimum(-0.15);

			AN_vs_pT[i][j]->SetMarkerSize(2);
			AN_vs_pT[i][j]->SetMarkerStyle(2+3*j);
			AN_vs_pT[i][j]->SetLineColor(2+method);
			AN_vs_pT[i][j]->Write();
		}	
	TCanvas *c0[2];
	
	c0[0] = new TCanvas("forward_plot","forward_plot",1000,800);
	c0[0]->cd();
	mg[0]->Add(AN_vs_pT[0][0]);
	mg[0]->Add(AN_vs_pT[1][1]);
	mg[0]->Draw("AP");
	l1->Draw("SAME");
	mg[0]->GetXaxis()->SetTitle("p_{T}");
        mg[0]->GetYaxis()->SetTitle("A_{N}");

	gPad->Modified();
	mg[0]->GetXaxis()->SetLimits(0,10.0);
        mg[0]->SetMinimum(-0.15);
        mg[0]->SetMaximum(0.15);
	leg1 = new TLegend(0.1,0.7,0.48,0.9);
	leg1->AddEntry(AN_vs_pT[0][0],"North arm_blue beam","lep");
	leg1->AddEntry(AN_vs_pT[1][1],"South arm_yellow beam","lep");
	leg1->Draw("SAME");

	c0[1] = new TCanvas("backward_plot","backward_plot",1000,800);
        c0[1]->cd();
	mg[1]->Add(AN_vs_pT[0][1]);
        mg[1]->Add(AN_vs_pT[1][0]);
	mg[1]->Draw("AP");
	l1->Draw("SAME");
        mg[1]->GetXaxis()->SetTitle("p_{T}");
        mg[1]->GetYaxis()->SetTitle("A_{N}");

        gPad->Modified();
        mg[1]->GetXaxis()->SetLimits(0,10.0);
        mg[1]->SetMinimum(-0.15);
        mg[1]->SetMaximum(0.15);
	leg2 = new TLegend(0.1,0.7,0.48,0.9);
        leg2->AddEntry(AN_vs_pT[0][1],"North arm_yellow beam","lep");
        leg2->AddEntry(AN_vs_pT[1][0],"South arm_blue beam","lep");
        leg2->Draw("SAME");

	
	//Draw it
//	TCanvas *c0[2][2];  //[forward/backward]
	TCanvas *c1[n_pt_bin];
	TCanvas *c2[n_pt_bin];
	mg[0]->Write();
	mg[1]->Write();
//	c0[0]->Write();
//	c0[1]->Write();

/*	for (int i = 0; i < 2; i++)
		{
			c0[i][j] = new TCanvas(Form("A_N_vs_pt_arm%d_beam%d",i,j), Form("A_N_vs_pt_arm%d_beam%d",i,j),800,600);
			c0[i][j]->cd();
			AN_vs_pT[i][j]->Draw("AP");
			l1->Draw("SAME");
			c0[i][j]->Write();
		}
*/	

	for(int i = 0; i<n_pt_bin; i++)
	{
	c1[i] = new TCanvas(Form("Canvas1_%d",i),Form("Canvas1_%d",i),1000,800);
	c1[i]->Divide(2,2);
	c1[i]->cd(1);
	phi_sig_n_A_L_Blue[i]->Draw("AP");
	f_nb[i]->Draw("SAME");
	f_nb[i]->Write();
	c1[i]->cd(3);
        phi_sig_n_A_L_Yellow[i]->Draw("AP");
	f_ny[i]->Draw("SAME");
	c1[i]->cd(2);
        phi_sig_s_A_L_Blue[i]->Draw("AP");
	f_sb[i]->Draw("SAME");
	c1[i]->cd(4);
        phi_sig_s_A_L_Yellow[i]->Draw("AP");
	f_sy[i]->Draw("SAME");

	c1[i]->Write();

	c2[i] = new TCanvas(Form("Canvas2_%d",i),Form("Canvas2_%d",i),1000,800);
        c2[i]->Divide(2,2);
        c2[i]->cd(1);
        phi_bgr_n_A_L_Blue[i]->Draw("AP");
        f_bnb[i]->Draw("SAME");
        c2[i]->cd(3);
        phi_bgr_n_A_L_Yellow[i]->Draw("AP");
        f_bny[i]->Draw("SAME");
        c2[i]->cd(2);
        phi_bgr_s_A_L_Blue[i]->Draw("AP");
        f_bsb[i]->Draw("SAME");
        c2[i]->cd(4);
        phi_bgr_s_A_L_Yellow[i]->Draw("AP");
        f_bsy[i]->Draw("SAME");

	c2[i]->Write();
	}	//end for drawing
}
