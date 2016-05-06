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
	
	TFile *f_out = new TFile("save_fitting_AN.root","RECREATE");
	TString infile ="432639_hist.root";

	double r[2][2];//[arm][pt]
	double xF_r[2][1];//[arm][xF]
	double average_pt[2][2];//[arm][pt_bin]
	double average_xf[2][1];//[arm][xf_bin]

	average_pt[0][0] = 1.081;
	average_pt[0][1] = 2.98;
	average_pt[1][0] = 1.068;
	average_pt[1][1] = 2.964;

	average_xf[0][0] = 0.5;
	average_xf[1][0] = 0.5;

	r[0][0] = 0.096;
        r[0][1] = 0.066;

        r[1][0] = 0.219;
        r[1][1] = 0.090;

	xF_r[0][0] = 0.105;
	xF_r[1][0] = 0.192;

	double AN_inc[2][2]; //[arm][pt_bin]	
	double AN_lsb[2][2]; //[arm][pt_bin]  
	double AN_phy[2][2]; //[arm][pt_bin] 

	double eAN_inc[2][2]; //[arm][pt_bin]   
        double eAN_lsb[2][2]; //[arm][pt_bin]  
        double eAN_phy[2][2]; //[arm][pt_bin]

	double xF_AN_inc[2][1]; //[arm][xF]    
        double xF_AN_lsb[2][1]; //[arm][xF]  
        double xF_AN_phy[2][1]; //[arm][xF] 

        double xF_eAN_inc[2][1]; //[arm][xF]   
        double xF_eAN_lsb[2][1]; //[arm][xF]  
        double xF_eAN_phy[2][1]; //[arm][xF]	

	gStyle->SetOptFit();
	
	TFile *file = TFile::Open(infile);
	//pT histogram
	//signal histogram
	TH2F *sig_n_A_L_Blue = (TH2F*)file->Get("_sah_pT_os_sig_n_A_L_Blue_SINGLE");
	TH2F *sig_s_A_L_Blue = (TH2F*)file->Get("_sah_pT_os_sig_s_A_L_Blue_SINGLE");

	//background histogram
	TH2F *bgr_n_A_L_Blue = (TH2F*)file->Get("_sah_pT_os_lsb_n_A_L_Blue_SINGLE");
        TH2F *bgr_s_A_L_Blue = (TH2F*)file->Get("_sah_pT_os_lsb_s_A_L_Blue_SINGLE");

	//xF histogram
	TH2F *xF_sig_n_A_L_Blue = (TH2F*)file->Get("_sah_xF_os_sig_n_A_L_Blue_SINGLE");
	TH2F *xF_sig_s_A_L_Blue = (TH2F*)file->Get("_sah_xF_os_sig_s_A_L_Blue_SINGLE");
	TH2F *xF_bgr_n_A_L_Blue = (TH2F*)file->Get("_sah_xF_os_lsb_n_A_L_Blue_SINGLE");
        TH2F *xF_bgr_s_A_L_Blue = (TH2F*)file->Get("_sah_xF_os_lsb_s_A_L_Blue_SINGLE");

	const Int_t n_pt_bin = 2;
	const Int_t n_phi_bin = 16;
// sig_n_A_L_Blue->GetNbinsX();

	//signal variables
	float al_b_n[n_pt_bin][n_phi_bin];
	float al_error_b_n[n_pt_bin][n_phi_bin];
	float phi_b_n[n_phi_bin];

	float al_b_s[n_pt_bin][n_phi_bin];
        float al_error_b_s[n_pt_bin][n_phi_bin];
        float phi_b_s[n_phi_bin];

	float xF_al_b_n[n_pt_bin][n_phi_bin];
        float xF_al_error_b_n[n_pt_bin][n_phi_bin];
        float xF_phi_b_n[n_phi_bin];

        float xF_al_b_s[n_pt_bin][n_phi_bin];
        float xF_al_error_b_s[n_pt_bin][n_phi_bin];
	float xF_phi_b_s[n_phi_bin];
	//background variables
	float bal_b_n[n_pt_bin][n_phi_bin];
        float bal_error_b_n[n_pt_bin][n_phi_bin];
        float bphi_b_n[n_phi_bin];

        float bal_b_s[n_pt_bin][n_phi_bin];
        float bal_error_b_s[n_pt_bin][n_phi_bin];
        float bphi_b_s[n_phi_bin];

	float xF_bal_b_n[n_pt_bin][n_phi_bin];
        float xF_bal_error_b_n[n_pt_bin][n_phi_bin];
        float xF_bphi_b_n[n_phi_bin];

        float xF_bal_b_s[n_pt_bin][n_phi_bin];
        float xF_bal_error_b_s[n_pt_bin][n_phi_bin];
        float xF_bphi_b_s[n_phi_bin];

	for (int j = 0; j<n_pt_bin; j++)
	{
		for (int i = 0; i<n_phi_bin; i++)
		{

		//give value to signal variable
			al_b_n[j][i] = sig_n_A_L_Blue->GetBinContent(j+1,i+1);
			al_error_b_n[j][i] = sig_n_A_L_Blue->GetBinError(j+1,i+1);
			phi_b_n[i] = sig_n_A_L_Blue->GetYaxis()->GetBinCenter(i+1);

			al_b_s[j][i] = sig_s_A_L_Blue->GetBinContent(j+1,i+1);
        	        al_error_b_s[j][i] = sig_s_A_L_Blue->GetBinError(j+1,i+1);
        	        phi_b_s[i] = sig_s_A_L_Blue->GetYaxis()->GetBinCenter(i+1);
	
		//give value to background variable
			bal_b_n[j][i] = bgr_n_A_L_Blue->GetBinContent(j+1,i+1);
        	        bal_error_b_n[j][i] = bgr_n_A_L_Blue->GetBinError(j+1,i+1);
        	        bphi_b_n[i] = bgr_n_A_L_Blue->GetYaxis()->GetBinCenter(i+1);

        	        bal_b_s[j][i] = bgr_s_A_L_Blue->GetBinContent(j+1,i+1);
        	        bal_error_b_s[j][i] = bgr_s_A_L_Blue->GetBinError(j+1,i+1);
        	        bphi_b_s[i] = bgr_s_A_L_Blue->GetYaxis()->GetBinCenter(i+1);
	
			if (j == 0)
			{
				xF_al_b_n[j][i] = xF_sig_n_A_L_Blue->GetBinContent(j+1,i+1);
                        	xF_al_error_b_n[j][i] = xF_sig_n_A_L_Blue->GetBinError(j+1,i+1);
                        	xF_phi_b_n[i] = xF_sig_n_A_L_Blue->GetYaxis()->GetBinCenter(i+1);

                        	xF_al_b_s[j][i] = xF_sig_s_A_L_Blue->GetBinContent(j+1,i+1);
                        	xF_al_error_b_s[j][i] = xF_sig_s_A_L_Blue->GetBinError(j+1,i+1);
                        	xF_phi_b_s[i] = xF_sig_s_A_L_Blue->GetYaxis()->GetBinCenter(i+1);

				xF_bal_b_n[j][i] = xF_bgr_n_A_L_Blue->GetBinContent(j+1,i+1);
                        	xF_bal_error_b_n[j][i] = xF_bgr_n_A_L_Blue->GetBinError(j+1,i+1);
                        	xF_bphi_b_n[i] = xF_bgr_n_A_L_Blue->GetYaxis()->GetBinCenter(i+1);

                        	xF_bal_b_s[j][i] = xF_bgr_s_A_L_Blue->GetBinContent(j+1,i+1);
                        	xF_bal_error_b_s[j][i] = xF_bgr_s_A_L_Blue->GetBinError(j+1,i+1);
                        	xF_bphi_b_s[i] = xF_bgr_s_A_L_Blue->GetYaxis()->GetBinCenter(i+1);
			}
		}
	}
	TGraphErrors *phi_sig_n_A_L_Blue[n_pt_bin];
	TGraphErrors *phi_sig_s_A_L_Blue[n_pt_bin];
	TGraphErrors *phi_bgr_n_A_L_Blue[n_pt_bin];
	TGraphErrors *phi_bgr_s_A_L_Blue[n_pt_bin];

	TGraphErrors *xF_phi_sig_n_A_L_Blue[1];
        TGraphErrors *xF_phi_sig_s_A_L_Blue[1];
        TGraphErrors *xF_phi_bgr_n_A_L_Blue[1];
        TGraphErrors *xF_phi_bgr_s_A_L_Blue[1];

	for (int i = 0; i < n_pt_bin; i++)
	{
	//Graph for signal
	phi_sig_n_A_L_Blue[i] = new TGraphErrors(n_phi_bin,phi_b_n,al_b_n[i],0,al_error_b_n[i]); 
	phi_sig_s_A_L_Blue[i] = new TGraphErrors(n_phi_bin,phi_b_s,al_b_s[i],0,al_error_b_s[i]);
	
	//Graph for background 
	phi_bgr_n_A_L_Blue[i] = new TGraphErrors(n_phi_bin,bphi_b_n,bal_b_n[i],0,bal_error_b_n[i]);
        phi_bgr_s_A_L_Blue[i] = new TGraphErrors(n_phi_bin,bphi_b_s,bal_b_s[i],0,bal_error_b_s[i]);

	phi_sig_n_A_L_Blue[i]->SetTitle(Form("North_incl_Blue_pTBin_%d",i));
	phi_sig_n_A_L_Blue[i]->GetXaxis()->SetTitle("phi");
	phi_sig_n_A_L_Blue[i]->GetYaxis()->SetTitle("AN");
	phi_sig_n_A_L_Blue[i]->SetMinimum(-0.4);
	phi_sig_n_A_L_Blue[i]->SetMaximum(0.4);

	phi_sig_s_A_L_Blue[i]->SetTitle(Form("South_incl_Blue_pTBin_%d",i));
        phi_sig_s_A_L_Blue[i]->GetXaxis()->SetTitle("phi");
        phi_sig_s_A_L_Blue[i]->GetYaxis()->SetTitle("AN");
	phi_sig_s_A_L_Blue[i]->SetMinimum(-0.4);
        phi_sig_s_A_L_Blue[i]->SetMaximum(0.4);

        phi_bgr_n_A_L_Blue[i]->SetTitle(Form("North_bgr_Blue_pTBin_%d",i));
        phi_bgr_n_A_L_Blue[i]->GetXaxis()->SetTitle("phi");
        phi_bgr_n_A_L_Blue[i]->GetYaxis()->SetTitle("AN");
	phi_bgr_n_A_L_Blue[i]->SetMinimum(-0.4);
        phi_bgr_n_A_L_Blue[i]->SetMaximum(0.4);

        phi_bgr_s_A_L_Blue[i]->SetTitle(Form("South_bgr_Blue_pTBin_%d",i));
        phi_bgr_s_A_L_Blue[i]->GetXaxis()->SetTitle("phi");
        phi_bgr_s_A_L_Blue[i]->GetYaxis()->SetTitle("AN");
	phi_bgr_s_A_L_Blue[i]->SetMinimum(-0.4);
        phi_bgr_s_A_L_Blue[i]->SetMaximum(0.4);
	}	//end for phi_Al TGraphErrors construction

	xF_phi_sig_n_A_L_Blue[0] = new TGraphErrors(n_phi_bin,xF_phi_b_n,xF_al_b_n[0],0,xF_al_error_b_n[0]);
        xF_phi_sig_s_A_L_Blue[0] = new TGraphErrors(n_phi_bin,xF_phi_b_s,xF_al_b_s[0],0,xF_al_error_b_s[0]);

	xF_phi_bgr_n_A_L_Blue[0] = new TGraphErrors(n_phi_bin,xF_bphi_b_n,xF_bal_b_n[0],0,xF_bal_error_b_n[0]);
        xF_phi_bgr_s_A_L_Blue[0] = new TGraphErrors(n_phi_bin,xF_bphi_b_s,xF_bal_b_s[0],0,xF_bal_error_b_s[0]);
	
	xF_phi_sig_n_A_L_Blue[0]->SetTitle("North_incl_Blue");
        xF_phi_sig_n_A_L_Blue[0]->GetXaxis()->SetTitle("phi");
        xF_phi_sig_n_A_L_Blue[0]->GetYaxis()->SetTitle("AN");
        xF_phi_sig_n_A_L_Blue[0]->SetMinimum(-0.4);
        xF_phi_sig_n_A_L_Blue[0]->SetMaximum(0.4);

        xF_phi_sig_s_A_L_Blue[0]->SetTitle("South_incl_Blue");
        xF_phi_sig_s_A_L_Blue[0]->GetXaxis()->SetTitle("phi");
        xF_phi_sig_s_A_L_Blue[0]->GetYaxis()->SetTitle("AN");
        xF_phi_sig_s_A_L_Blue[0]->SetMinimum(-0.4);
        xF_phi_sig_s_A_L_Blue[0]->SetMaximum(0.4);

        xF_phi_bgr_n_A_L_Blue[0]->SetTitle("North_bgr_Blue");
        xF_phi_bgr_n_A_L_Blue[0]->GetXaxis()->SetTitle("phi");
        xF_phi_bgr_n_A_L_Blue[0]->GetYaxis()->SetTitle("AN");
        xF_phi_bgr_n_A_L_Blue[0]->SetMinimum(-0.4);
        xF_phi_bgr_n_A_L_Blue[0]->SetMaximum(0.4);

        xF_phi_bgr_s_A_L_Blue[0]->SetTitle("South_bgr_Blue");
        xF_phi_bgr_s_A_L_Blue[0]->GetXaxis()->SetTitle("phi");
        xF_phi_bgr_s_A_L_Blue[0]->GetYaxis()->SetTitle("AN");
        xF_phi_bgr_s_A_L_Blue[0]->SetMinimum(-0.4);
        xF_phi_bgr_s_A_L_Blue[0]->SetMaximum(0.4);
	

	TF1 *f_nb[n_pt_bin];
	TF1 *f_sb[n_pt_bin];
	TF1 *f_bnb[n_pt_bin];
	TF1 *f_bsb[n_pt_bin];

	TF1 *xF_f_nb[1];
        TF1 *xF_f_sb[1];
        TF1 *xF_f_bnb[1];
        TF1 *xF_f_bsb[1];

	for(int i =0; i<n_pt_bin;i++)
	{

		if (method == 2)
		{
			f_nb[i] = new TF1("f_nb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        		f_sb[i] = new TF1("f_sb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
	

			f_bnb[i] = new TF1("f_bnb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        		f_bsb[i] = new TF1("f_bsb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
		}
	}   //end for tf construction

	xF_f_nb[0] = new TF1("xF_f_nb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        xF_f_sb[0] = new TF1("xF_f_sb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);

        xF_f_bnb[0] = new TF1("xF_f_bnb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);
        xF_f_bsb[0] = new TF1("xF_f_bsb","[1]*TMath::Sin(x)+[0]*TMath::Cos(x)",-3.1415926,3.1415926);

	for (int i =0; i< n_pt_bin; i++)
        {
		f_nb[i]->SetParName(0,"AN");
		f_sb[i]->SetParName(0,"AN");
		f_bnb[i]->SetParName(0,"AN");	
		f_bsb[i]->SetParName(0,"AN");
	}	
	//fit 
	
	xF_f_nb[0]->SetParName(0,"AN");
        xF_f_sb[0]->SetParName(0,"AN");
        xF_f_bnb[0]->SetParName(0,"AN");
        xF_f_bsb[0]->SetParName(0,"AN");
	
	for (int i =0; i< n_pt_bin; i++)
	{
	
		phi_sig_n_A_L_Blue[i]->Fit(f_nb[i]);
        	phi_sig_s_A_L_Blue[i]->Fit(f_sb[i]);

		phi_bgr_n_A_L_Blue[i]->Fit(f_bnb[i]);
        	phi_bgr_s_A_L_Blue[i]->Fit(f_bsb[i]);

		if (i == 0)
		{
			xF_phi_sig_n_A_L_Blue[i]->Fit(xF_f_nb[i]);
                	xF_phi_sig_s_A_L_Blue[i]->Fit(xF_f_sb[i]);

                	xF_phi_bgr_n_A_L_Blue[i]->Fit(xF_f_bnb[i]);
                	xF_phi_bgr_s_A_L_Blue[i]->Fit(xF_f_bsb[i]);
		}

	}   //end for fit


	for (int i = 0; i< n_pt_bin; i++)
        {
		AN_inc[0][i] = f_nb[i]->GetParameter(0);
		eAN_inc[0][i] = f_nb[i]->GetParError(0);
	
		AN_lsb[0][i] = f_bnb[i]->GetParameter(0);
		eAN_lsb[0][i] = f_bnb[i]->GetParError(0);

		AN_inc[1][i] = f_sb[i]->GetParameter(0);
		eAN_inc[1][i]= f_sb[i]->GetParError(0);
	
		AN_lsb[1][i] = f_bsb[i]->GetParameter(0);
                eAN_lsb[1][i]= f_bsb[i]->GetParError(0);
	}

	xF_AN_inc[0][0] = xF_f_nb[0]->GetParameter(0);
        xF_eAN_inc[0][0] = xF_f_nb[0]->GetParError(0);

        xF_AN_lsb[0][0] = xF_f_bnb[0]->GetParameter(0);
        xF_eAN_lsb[0][0] = xF_f_bnb[0]->GetParError(0);

        xF_AN_inc[1][0] = xF_f_sb[0]->GetParameter(0);
        xF_eAN_inc[1][0]= xF_f_sb[0]->GetParError(0);

        xF_AN_lsb[1][0] = xF_f_bsb[0]->GetParameter(0);
        xF_eAN_lsb[1][0]= xF_f_bsb[0]->GetParError(0);

	
	for (int i = 0; i<2; i++)
		for (int j = 0; j<n_pt_bin; j++)
		{
			AN_phy[i][j] = (AN_inc[i][j] - r[i][j]*AN_lsb[i][j])/(1. - r[i][j]);

			eAN_phy[i][j]= sqrt(eAN_inc[i][j]*eAN_inc[i][j] +
                                       r[i][j]*r[i][j]*eAN_lsb[i][j]*eAN_lsb[i][j]) / (1. - r[i][j]);
			if (j == 0)
			{
				xF_AN_phy[i][j] = (xF_AN_inc[i][j] - xF_r[i][j]*xF_AN_lsb[i][j])/(1. - xF_r[i][j]);

                        	xF_eAN_phy[i][j]= sqrt(xF_eAN_inc[i][j]*xF_eAN_inc[i][j] +
                                	       xF_r[i][j]*xF_r[i][j]*xF_eAN_lsb[i][j]*xF_eAN_lsb[i][j]) / (1. - xF_r[i][j]);	
				cout<<"xF_AN_phy: "<<xF_AN_phy[i][j]<< "   Error: "<< xF_eAN_phy[i][j]<<endl;
			}
		}

	TGraphErrors *AN_vs_pT[2];
	TGraphErrors *AN_vs_xF[2];
	TLegend *leg[2];
	TLegend *xF_leg[2];
	TLine *l1 = new TLine(0,0,10.0,0);
	TLine *l2 = new TLine(0,0,1.0,0);	

	TCanvas *c;
	c = new TCanvas("c","c",1000,400);
	c->Divide(2,1);

	TCanvas *c_xF;
	c_xF = new TCanvas("c_xF","c_xF",1000,400);
	c_xF->Divide(2,1);

	f_out->cd();	
	for (int i = 0; i < 2; i++) //arm
	{
		AN_vs_pT[i] = new TGraphErrors(n_pt_bin,average_pt[i],AN_phy[i],0,eAN_phy[i]);
                AN_vs_pT[i]->SetLineWidth(2);
		AN_vs_pT[i]->SetNameTitle(Form("AN vs p_{T} Arm%d",i),Form("AN vs p_{T} Arm%d",i));
		AN_vs_pT[i]->GetXaxis()->SetTitle("p_{T}(MeV)");
                AN_vs_pT[i]->GetXaxis()->SetLimits(0,10.0);
                AN_vs_pT[i]->GetYaxis()->SetTitle("A_{N}");
                AN_vs_pT[i]->SetMaximum(0.15);
                AN_vs_pT[i]->SetMinimum(-0.15);	
		AN_vs_pT[i]->SetMarkerSize(2);
		AN_vs_pT[i]->SetMarkerStyle(2+3*i);

		AN_vs_xF[i] = new TGraphErrors(1,average_xf[i],xF_AN_phy[i],0,xF_eAN_phy[i]);
                AN_vs_xF[i]->SetLineWidth(2);
                AN_vs_xF[i]->SetNameTitle(Form("Arm%d",i),Form("Arm%d",i));
                AN_vs_xF[i]->GetXaxis()->SetTitle("p_{T}(MeV)");
                AN_vs_xF[i]->GetXaxis()->SetLimits(0,1.0);
                AN_vs_xF[i]->GetYaxis()->SetTitle("A_{N}");
                AN_vs_xF[i]->SetMaximum(0.15);
                AN_vs_xF[i]->SetMinimum(-0.15);
                AN_vs_xF[i]->SetMarkerSize(2);
                AN_vs_xF[i]->SetMarkerStyle(2+3*i);
	}
	for (int i = 0; i < 2; i++) //arm
        {
		c->cd(i+1);
		AN_vs_pT[i]->Draw("AP");
		leg[i] = new TLegend(0.1,0.7,0.48,0.9);
		if(i == 0)
		leg[i]->AddEntry(AN_vs_pT[i],"North Arm A_{N} vs p_{T}","lep");
		else if (i ==1)
		leg[i]->AddEntry(AN_vs_pT[i],"South Arm A_{N} vs p_{T}","lep");
		l1->Draw("SAME");
		leg[i]->Draw("SAME");

		c_xF->cd(i+1);
		AN_vs_xF[i]->Draw("AP");
                xF_leg[i] = new TLegend(0.1,0.7,0.48,0.9);
                if(i == 0)
                xF_leg[i]->AddEntry(AN_vs_xF[i],"North Arm A_{N} vs. x_{F} ","lep");
                else if (i ==1)
                xF_leg[i]->AddEntry(AN_vs_xF[i],"South Arm A_{N} vs. x_{F} ","lep");
                l2->Draw("SAME");
                xF_leg[i]->Draw("SAME");
		
		c->Write();
	}

	
	TCanvas *c1[n_pt_bin];
	for(int i = 0; i<n_pt_bin; i++)
	{
		c1[i] = new TCanvas(Form("Canvas1_%d",i),Form("Canvas1_%d",i),1000,800);
		c1[i]->Divide(2,2);
		c1[i]->cd(1);
		phi_sig_n_A_L_Blue[i]->Draw("AP");
		f_nb[i]->Draw("SAME");
		c1[i]->cd(2);
        	phi_sig_s_A_L_Blue[i]->Draw("AP");
		f_sb[i]->Draw("SAME");
		c1[i]->cd(3);
		phi_bgr_n_A_L_Blue[i]->Draw("AP");
		f_bnb[i]->Draw("SAME");
		c1[i]->cd(4);
		phi_bgr_s_A_L_Blue[i]->Draw("AP");
		f_bsb[i]->Draw("SAME");
		c1[i]->Write();

	}	//end for pT drawing

	TCanvas *c2;
	c2 = new TCanvas("xF vs AN","xF vs AN",1000,800);
	c2->Divide(2,2);
	c2->cd(1);
	xF_phi_sig_n_A_L_Blue[0]->Draw("AP");
        xF_f_nb[0]->Draw("SAME");

	c2->cd(2);
	xF_phi_sig_s_A_L_Blue[0]->Draw("AP");
        xF_f_sb[0]->Draw("SAME");

	c2->cd(3);
	xF_phi_bgr_n_A_L_Blue[0]->Draw("AP");
        xF_f_bnb[0]->Draw("SAME");

	c2->cd(4);
	xF_phi_bgr_s_A_L_Blue[0]->Draw("AP");
        xF_f_bsb[0]->Draw("SAME");
	c2->Write();
}
