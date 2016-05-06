#include <iostream>
#include <iomanip>
#include <fstream>


#include "TFile.h"
#include "TTree.h"
static const int NCROSS=120;

#define OUTPUT "GL1Pdata.root"

void createSpinTree
(
 //int run1 = 365513, int run2 = 365513//TEST RUN 12
 //int run1 = 364822, int run2 = 368798//ALL RUN 12
 int run1 = 422070, int run2 = 432008//ALL RUN 15 pp200
//   int run1 = 422, int run2 = 425011
 )
{
	int runNum = 0;
	int filNum = 0;
	int XingShift = 0;
	float polB = 0.;
	float polY = 0.;
	float polBstat = 0.;
	float polBsyst = 0.;
	float polYstat = 0.;
	float polYsyst = 0.;
	short spinB[120];
	short spinY[120];
	short spinB_PHENIX[120];
	short spinY_PHENIX[120];
	float bbc_in[120];
	//float bbc_in_corr[120];
	float bbc_out[120];
	float zdc_in[120];
	float zdc_out[120];
	int crossing_PHENIX[120];
	int crossing_RHIC[120];

	for (int i=0; i<120; i++)
	{
		spinB[i] = 0;
		spinY[i] = 0;
		spinB_PHENIX[i] = 0;
		spinY_PHENIX[i] = 0;
		bbc_in[i] = 0;
		//bbc_in_corr[i] = 0;
		bbc_out[i] = 0;
		zdc_in[i] = 0;
		zdc_out[i] = 0;
	}

	TFile *newF = new TFile(OUTPUT,"recreate");
	TTree *newT = new TTree("T","polarization and BBCLL1 data");

	newT->Branch("RunNumber",  &runNum,   "RunNumber/I");
	newT->Branch("FillNumber", &filNum,   "FillNumber/I");
	newT->Branch("XingShift",  &XingShift,"XingShift/I");
	newT->Branch("PolB",       &polB,     "PolB/F");
	newT->Branch("ePolBstat",  &polBstat, "ePolBstat/F");
	newT->Branch("ePolBsyst",  &polBsyst, "ePolBsyst/F");
	newT->Branch("PolY",       &polY,     "PolY/F");
	newT->Branch("ePolYstat",  &polYstat, "ePolYstat/F");
	newT->Branch("ePolYsyst",  &polYsyst, "ePolYsyst/F");
	newT->Branch("SpinB_PHENIX",      spinB_PHENIX,     "SpinB_PHENIX[120]/S");
	newT->Branch("SpinY_PHENIX",      spinY_PHENIX,     "SpinY_PHENIX[120]/S");
	newT->Branch("SpinB",      spinB,     "SpinB[120]/S");
	newT->Branch("SpinY",      spinY,     "SpinY[120]/S");
	newT->Branch("BBCin",      bbc_in,    "BBCin[120]/F");
	//newT->Branch("BBCin_corr",      bbc_in_corr,    "BBCin_corr[120]/F");
	newT->Branch("BBCout",      bbc_out,    "BBCout[120]/F");
	newT->Branch("ZDCin",      zdc_in,    "ZDCin[120]/F");
	newT->Branch("ZDCout",      zdc_out,    "ZDCout[120]/F");
	newT->Branch("crossing_PHENIX",      crossing_PHENIX,    "crossing_PHENIX[120]/I");
	newT->Branch("crossing_RHIC",      crossing_RHIC,    "crossing_RHIC[120]/I");

	gSystem->Load("libuspin.so");
	SpinDBOutput spin_out("phnxrc");
	SpinDBContent spin_cont;

  //gSystem->Load("/direct/phenix+spin2/chenxu/CVS_PHENIX/install/lib/libScalerCorrection.so");
 // ScalerCorrection* scaler = new ScalerCorrection();

	cout<<"begin loop"<<endl;

	int badrunNum = -999;
	for(int run=run1; run<run2+1; run++)
	{
		spin_out.StoreDBContent(run,run);
		if(spin_out.CheckRunRowStore(run)!=1){continue;}
		spin_out.GetDBContentStore(spin_cont,run);

		runNum = spin_cont.GetRunNumber();
		filNum = spin_cont.GetFillNumber();
		XingShift = spin_cont.GetCrossingShift();

		if(runNum==badrunNum)continue;

		if(XingShift<0)
		{
			cout
				<<"XingShift Error: "
				<<" runNum: "<<runNum
				<<" filNum: "<<filNum
				<<" XingShift: "<<XingShift
				<<endl;
			continue;
		}

		//Long64_t Scaler_BBC_Residual[120];
		//scaler->GetScalerCount(run, "BBC", "Residual", Scaler_BBC_Residual);

		for(int i=0; i<NCROSS; i++)
		{
			double bpol,bpolerr,ypol,ypolerr;

			int Xing_PHENIX = i;
			int Xing_RHIC = (i+XingShift)%NCROSS;

			crossing_PHENIX[i] = Xing_PHENIX;
			crossing_RHIC[i] = Xing_RHIC;

			//if(Xing_RHIC<0||Xing_RHIC>=120||Xing_PHENIX<0||Xing_PHENIX>=120)
			if(XingShift!=5)
			{
				if(Xing_PHENIX == 0)
				{
					cout
						<<"XingShift Abnormal: "
						<<" runNum: "<<runNum
						<<" filNum: "<<filNum
						<<" XingShift: "<<XingShift
						<<" Xing_PHENIX: "<<Xing_PHENIX
						<<" Xing_RHIC: "<<Xing_RHIC
						<<endl;
				}
			}

			spinB_PHENIX[Xing_PHENIX] = spin_cont.GetSpinPatternBluePHENIX(Xing_PHENIX);
			spinY_PHENIX[Xing_PHENIX] = spin_cont.GetSpinPatternYellowPHENIX(Xing_PHENIX);
			//bbc_in[Xing_PHENIX]  = spin_cont.GetScalerBbcVertexCutPHENIX(Xing_PHENIX);
			//bbc_out[Xing_PHENIX] = spin_cont.GetScalerBbcNoCutPHENIX(Xing_PHENIX);
			//zdc_in[Xing_PHENIX]  = spin_cont.GetScalerZdcNarrowPHENIX(Xing_PHENIX);
			//zdc_out[Xing_PHENIX] = spin_cont.GetScalerZdcWidePHENIX(Xing_PHENIX);
			//spin_cont.GetPolarizationBluePHENIX(Xing_PHENIX,bpol,bpolerr);
			//spin_cont.GetPolarizationYellowPHENIX(Xing_PHENIX,ypol,ypolerr);

			spinB[Xing_RHIC]   = spin_cont.GetSpinPatternBlue(Xing_RHIC);
			spinY[Xing_RHIC]   = spin_cont.GetSpinPatternYellow(Xing_RHIC);
			bbc_in[Xing_RHIC]  = spin_cont.GetScalerBbcVertexCut(Xing_RHIC);
			//if(Scaler_BBC_Residual[Xing_RHIC]>0)
			//	bbc_in_corr[Xing_RHIC] = Scaler_BBC_Residual[Xing_RHIC];
			//else
			//	bbc_in_corr[Xing_RHIC] = bbc_in[Xing_RHIC];
			bbc_out[Xing_RHIC] = spin_cont.GetScalerBbcNoCut(Xing_RHIC);
			zdc_in[Xing_RHIC]  = spin_cont.GetScalerZdcNarrow(Xing_RHIC);
			zdc_out[Xing_RHIC] = spin_cont.GetScalerZdcWide(Xing_RHIC);
			spin_cont.GetPolarizationBlue(Xing_RHIC,bpol,bpolerr);
			spin_cont.GetPolarizationYellow(Xing_RHIC,ypol,ypolerr);

			if(spinB[Xing_RHIC]!=spinB_PHENIX[Xing_PHENIX])
			{
				cout
					<<"Method Error: "
					<<" runNum: "<<runNum
					<<" filNum: "<<filNum
					<<" spinB[Xing_PHENIX]: "<<Xing_PHENIX<<"\t"<<spinB[Xing_PHENIX]
					<<" spinB[Xing_RHIC]: "<<Xing_RHIC<<"\t"<<spinB[Xing_RHIC]
					<<endl;
				continue;
			}


			polB = bpol;
			polY = ypol;
			polBstat = bpolerr;
			polBsyst = 0.;
			polYstat = ypolerr;
			polYsyst = 0.;

			if(runNum!=badrunNum)
			{
				if(spinB[Xing_RHIC]<-100)  {cout<<"Pattern Error: spinB[Xing_RHIC]<-100  "<<setw(10)<<spinB[Xing_RHIC]  <<" "<<"runNum: "<<runNum<<" filNum: "<<filNum<<" Xing_RHIC: "<<Xing_RHIC<<endl; badrunNum = runNum; continue;} 
				if(spinY[Xing_RHIC]<-100)  {cout<<"Pattern Error: spinY[Xing_RHIC]<-100  "<<setw(10)<<spinY[Xing_RHIC]  <<" "<<"runNum: "<<runNum<<" filNum: "<<filNum<<" Xing_RHIC: "<<Xing_RHIC<<endl; badrunNum = runNum; continue;}

				if(bbc_in[Xing_RHIC]<0)    {cout<<"Scaler Error: bbc_in[Xing_RHIC]<0     "<<setw(10)<<bbc_in[Xing_RHIC] <<" "<<"runNum: "<<runNum<<" filNum: "<<filNum<<" Xing_RHIC: "<<Xing_RHIC<<endl; badrunNum = runNum; continue;}
				if(bbc_out[Xing_RHIC]<0)   {cout<<"Scaler Error: bbc_out[Xing_RHIC]<0    "<<setw(10)<<bbc_out[Xing_RHIC]<<" "<<"runNum: "<<runNum<<" filNum: "<<filNum<<" Xing_RHIC: "<<Xing_RHIC<<endl; badrunNum = runNum; continue;}
				if(zdc_in[Xing_RHIC]<0)    {cout<<"Scaler Error: zdc_in[Xing_RHIC]<0     "<<setw(10)<<zdc_in[Xing_RHIC] <<" "<<"runNum: "<<runNum<<" filNum: "<<filNum<<" Xing_RHIC: "<<Xing_RHIC<<endl; badrunNum = runNum; continue;}
				if(zdc_out[Xing_RHIC]<0)   {cout<<"Scaler Error: zdc_out[Xing_RHIC]<0    "<<setw(10)<<zdc_out[Xing_RHIC]<<" "<<"runNum: "<<runNum<<" filNum: "<<filNum<<" Xing_RHIC: "<<Xing_RHIC<<endl; badrunNum = runNum; continue;}
				if(XingShift<0)    {cout<<"Scaler Error: XingShift[Xing_RHIC]<0  "<<setw(10)<<XingShift <<" "<<"runNum: "<<runNum<<" filNum: "<<filNum<<" Xing_RHIC: "<<Xing_RHIC<<endl; badrunNum = runNum; continue;}
			}
		}

		double tc_x_blue,tc_x_blueerr,tc_y_blue,tc_y_blueerr;
		double tc_x_yellow,tc_x_yellowerr,tc_y_yellow,tc_y_yellowerr;
		spin_cont.GetTransCompBlueX(tc_x_blue,tc_x_blueerr);
		spin_cont.GetTransCompBlueY(tc_y_blue,tc_y_blueerr);
		spin_cont.GetTransCompYellowX(tc_x_yellow,tc_x_yellowerr);
		spin_cont.GetTransCompYellowY(tc_y_yellow,tc_y_yellowerr);

		newT->Fill();
	}

	cout<<"Finished!!"<<endl;

	newF->cd();
	newT->Write();
	newF->Close();
}
