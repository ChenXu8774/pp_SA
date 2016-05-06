#include "JpsiAN.hh"

#include <iostream>
#include <iomanip>
#include <fstream>


#include "TF1.h"
#include "TPDF.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TAxis.h"
#include "TString.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TPavesText.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAttText.h"
#include "TStyle.h"


#define  MinMs      1.5
#define  MaxMs      7.5
#define  MsRange    6.0    // MaxMs - MinMs
#define  SIG        2.0

#define  n_phi_bin  16
#define  n_pt_bin   2
// background mass range
#define  Bg1_MinMs  1.5
#define  Bg1_MaxMs  2.4


// values from PDG  [GeV]
#define  MuMass     0.105658367
#define  JpsiMass   0.3096916
#define  PsiMass    0.368609


#define  Deg2Rad    0.017453293
#define  PI         3.14159265358979323846
#define  TwoPI      6.28318530717958647692
#define  PIby2      1.57079632679489661923
#define  SqrtTwoPI  2.50662827463100050241


#define  NBins      120
#define  FitLim1    2.0
#define  FitLim2    4.5

#define  XINGs           120
#define  BlueSpinShift  -0.24 // different in run 15
#define  YeloSpinShift   0.   // different in run 15

#define  NCUT  20    // what is NCUT? each arm must have more than NCUT j/psi per fill

//#define  RUNLIST   "goodrun.list"
//#define  RUNLIST   "bbc_abortgap.list"
#define RUNLIST   "goodrun.list"

using namespace std;

JpsiAN::JpsiAN()
{
	TFile *f_out = new TFile("save_histo.root", "RECREATE");
  	ReadGoodRunList();
  	ReadFillTable();

 	SetHistos();

  	_file = new TFile("Run15_pp_goodrun.root");
//  	_t = (TTree*)_file->Get("dimuon");

  	_t = (TTree*)_file->Get("T");
  	SetTree();

  	set<int> fill;

  	for (int evt=0; evt<_t->GetEntries(); evt++)
    	{
      		_t->GetEntry(evt);

      		if (Cut()) continue;
      		if (charge != 0) continue;

      		fill.insert(_fillTable.find(Run_Number)->second);
    	}


  	_nFills = fill.size();
  	cout << "#fills  " << _nFills << endl;
  	_Fill = new double[_nFills];


  	int count=0;
  	set<int>::iterator itr;
  	for (itr=fill.begin(); itr!=fill.end(); itr++) _Fill[count++] = *itr;


  	SetArrays();    // this requires _nFills
  	GetSpinInfo();

  	RelLumi();
 	GetYeilds();
	GetAN();
	f_out->cd();
	c0->Write();
	c_yeild[0]->Write();
	c_yeild[1]->Write();
	cPt->Write();
}
//-----------------------------------------------------------------------------

JpsiAN::~JpsiAN()
{
  	delete [] _Fill;
  	delete [] _eFill;

  	for (int i=0; i<2; i++)
    	{


      		delete [] _R1[i];
      		delete [] _R2[i];
      		delete [] _R3[i];

      		delete [] _eR1[i];
      		delete [] _eR2[i];
      		delete [] _eR3[i];
    	}

  	map<int, short*>::iterator itr;

  	for (itr = _spinB.begin(); itr != _spinB.end(); itr++)
    		delete [] itr->second;

  	for (itr = _spinY.begin(); itr != _spinY.end(); itr++)
    		delete [] itr->second;


  	map<int, double*>::iterator itr2;
  	for (itr2 = _BBCinB.begin(); itr2 != _BBCinB.end(); itr2++)
    		delete [] itr2->second;
  	for (itr2 = _BBCinY.begin(); itr2 != _BBCinY.end(); itr2++)
    		delete [] itr2->second;

  	_file->Close();
  	delete _file;
}
//-----------------------------------------------------------------------------


bool JpsiAN::Cut()
{
  	// returns true if cut is applied

  	if (_runList.count(Run_Number) == 0) return true;

  	if (same_event != true) return true;
  	if (mass<=1.0 || mass>=4.0) return true;
  	if (Tr0_pz*Tr1_pz < 0.) return true;
  	if (charge !=0) return true;
  	if (fabs(Evt_bbcZ) == 0.) return true;
  	if (fabs(Evt_bbcZ) >= 30.) return true;

  	if (Tr0_pz >0 && (fabs(Tr0_DG0) >= 25. || fabs(Tr1_DG0) >= 25.)) return true;
  	if (Tr0_pz <0 && (fabs(Tr0_DG0) >= 30. || fabs(Tr1_DG0) >= 30.)) return true;

 	if (Tr0_DDG0 >= 10. || Tr1_DDG0 >= 10.) return true;
  	if (Tr0_ntrhits<=9 || Tr1_ntrhits<=9) return true;
  	if (Tr0_nidhits<=5 || Tr1_nidhits<=5) return true;
  	if (Tr0_lastgap<=2 || Tr1_lastgap<=2) return true;
  	if (Tr0_dca_r>5. || Tr1_dca_r>5.) return true;

  	if ( pT >= 6.) return true;

  	if (pz >= 100.) return true;
  	if (fabs(rapidity) <= 1.2 || fabs(rapidity) >= 2.2) return true;

  	if (Tr0_pz >0 && (!((lvl1_trigscaled & 0x00400000)||(lvl1_trigscaled &0x00100000)))) return true;
  	if (Tr0_pz <0 && (!((lvl1_trigscaled & 0x00800000)||(lvl1_trigscaled &0x00200000)))) return true;

  	if (Evt_vtxchi2 >= 5.) return true;

  	return false;
}
//-----------------------------------------------------------------------------


void JpsiAN::SetHistos()
{
  	cout<<"Sethistos................"<<endl;
	
	_hMass[North] = new TH1D("hMassNorth","hMassNorth", NBins, 2, 4.5);
  	_hMass[South] = new TH1D("hMassSouth","hMassSouth", NBins, 2, 4.5);

  	_hXF[North]   = new TH1D("hXFnorth","hXFnorth",100, 0., 0.5);
  	_hXF[South]   = new TH1D("hXFsouth","hXFsouth",100, 0., 0.5);

  	_hPt[North][0]   = new TH1D("hPtNorth_bin0","hPtNorth_bin0",60, 0., 6.);
  	_hPt[South][0]   = new TH1D("hPtSouth_bin0","hPtSouth_bin0",60, 0., 6.);
  
	_hPt[North][1]   = new TH1D("hPtNorth_bin1","hPtNorth_bin1",60, 0., 6.);
        _hPt[South][1]   = new TH1D("hPtSouth_bin1","hPtSouth_bin1",60, 0., 6.);
}
//-----------------------------------------------------------------------------

void JpsiAN::GetYeilds()
{
  	cout << "JpsiAN::GetYeilds" << endl;

  	int arm  = -9;

  	TVector3 *vDimu = new TVector3(0.,0.,0.);
  	double phi;
	
  	int fill, Xing;

  	short spinB, spinY;

	double pT_MinJpsiMs[2][n_pt_bin];
	double pT_MaxJpsiMs[2][n_pt_bin];

	pT_MinJpsiMs[North][0] =2.88121;
        pT_MaxJpsiMs[North][0] =3.43513;
        pT_MinJpsiMs[North][1] =2.87472;
        pT_MaxJpsiMs[North][1] =3.45992;

        pT_MinJpsiMs[South][0] =2.84384;
        pT_MaxJpsiMs[South][0] =3.39732;
        pT_MinJpsiMs[South][1] =2.84184;
        pT_MaxJpsiMs[South][1] =3.40799;

        const Double_t pt_edge[n_pt_bin + 1] = {0,2,10};
        for (int i =0; i <2; i++)
        {
                _hIcl_u_u_B[i] = new TH2D(Form("incl.arm%d_u_u",i),Form("incl.arm%d_u_u",i),n_pt_bin,pt_edge,n_phi_bin,-PI,PI);
                _hIcl_u_d_B[i] = new TH2D(Form("incl.arm%d_u_d",i),Form("incl.arm%d_u_d",i),n_pt_bin,pt_edge,n_phi_bin,-PI,PI);
                _hIcl_d_u_B[i] = new TH2D(Form("incl.arm%d_d_u",i),Form("incl.arm%d_d_u",i),n_pt_bin,pt_edge,n_phi_bin,-PI,PI);
                _hIcl_d_d_B[i] = new TH2D(Form("incl.arm%d_d_d",i),Form("incl.arm%d_d_d",i),n_pt_bin,pt_edge,n_phi_bin,-PI,PI);

                _hBgr_u_u_B[i] = new TH2D(Form("bgrd.arm%d_u_u",i),Form("bgrd.arm%d_u_u",i),n_pt_bin,pt_edge,n_phi_bin,-PI,PI);
                _hBgr_u_d_B[i] = new TH2D(Form("bgrd.arm%d_u_d",i),Form("bgrd.arm%d_u_d",i),n_pt_bin,pt_edge,n_phi_bin,-PI,PI);
                _hBgr_d_u_B[i] = new TH2D(Form("bgrd.arm%d_d_u",i),Form("bgrd.arm%d_d_u",i),n_pt_bin,pt_edge,n_phi_bin,-PI,PI);
                _hBgr_d_d_B[i] = new TH2D(Form("bgrd.arm%d_d_d",i),Form("bgrd.arm%d_d_d",i),n_pt_bin,pt_edge,n_phi_bin,-PI,PI);
        }

	ofstream event("event_output.txt");


	double JPSI_SIG_MASS_MAX = -1.;
	double JPSI_SIG_MASS_MIN = -1.;

  	for (int evt=0; evt<_t->GetEntries(); evt++)
    	{
      		_t->GetEntry(evt);

      		if (Cut()) continue;

      		arm = rapidity > 0 ? North : South;      

      		vDimu->SetXYZ(px, py, pz);

      		phi = double(vDimu->Phi());

      		fill = _fillTable.find(Run_Number)->second;
		
      		Xing = int(SpinX_ID);

  		// NOTE: blue/yellow beam goes towards north/south.
      		spinB = _spinB.find(fill)->second[Xing];
      		spinY = _spinY.find(fill)->second[Xing];

		int pt_bin = -1;
		for(int i=0;i<2;i++)
			if (pT >pt_edge[i] && pT<pt_edge[i+1])
			{
				JPSI_SIG_MASS_MAX = pT_MaxJpsiMs[arm][i];
				JPSI_SIG_MASS_MIN = pT_MinJpsiMs[arm][i];
				pt_bin = i;		
			}
		if (mass > JPSI_SIG_MASS_MIN && mass<JPSI_SIG_MASS_MAX && charge==0)
//      	if (mass > _MinJpsiMs[arm] && mass < _MaxJpsiMs[arm] && charge == 0)
		{
	  	//fout << Run_Number << endl;
			if(spinB == 1 && spinY == 1) _hIcl_u_u_B[arm]->Fill(pT,phi);
			if(spinB == 1 && spinY ==-1) _hIcl_u_d_B[arm]->Fill(pT,phi);
			if(spinB ==-1 && spinY == 1) _hIcl_d_u_B[arm]->Fill(pT,phi);
			if(spinB ==-1 && spinY ==-1) _hIcl_d_d_B[arm]->Fill(pT,phi);
			
			_hit_plots[arm]->Fill(SpinX_ID,_RunTable.find(Run_Number)->second);
	  		_hXF[arm]->Fill(fabs(xF));
	  		_hPt[arm][pt_bin]->Fill(pT);

		}

      		if (mass > Bg1_MinMs && mass < Bg1_MaxMs && charge == 0)
		{
			if(spinB == 1 && spinY == 1) _hBgr_u_u_B[arm]->Fill(pT,phi);
                        if(spinB == 1 && spinY ==-1) _hBgr_u_d_B[arm]->Fill(pT,phi);
                        if(spinB ==-1 && spinY == 1) _hBgr_d_u_B[arm]->Fill(pT,phi);
                        if(spinB ==-1 && spinY ==-1) _hBgr_d_d_B[arm]->Fill(pT,phi);
			_hit_plots[arm]->Fill(SpinX_ID,_RunTable.find(Run_Number)->second);
		}

                event<<
                        setw(15)<<" run num: "<<setw(10)<<Run_Number
                        <<setw(10)<<" beamclk0: "<<setw(15)<<beamclk0
                        <<setw(10)<<" mass: "<<setw(10)<<mass
                        <<setw(6)<<" pT: "<<setw(8)<<pT
                        <<setw(8)<<" phi: "<<setw(8)<<phi
			<<setw(8)<<" Xing: "<<setw(8)<<SpinX_ID
                        <<endl;

	

    	}
	
	ofstream n_hit0("hit_plots_s_sig.txt");
	ofstream n_hit1("hit_plots_n_sig.txt");	
	for (int k=0; k<n_phi_bin; k++)
	{
		n_hit0<<"event:"<<setw(3)<<k<<" : "
		<<" uu "<<_hIcl_u_u_B[1]->GetBinContent(1,k+1)
		<<" ud "<<_hIcl_u_d_B[1]->GetBinContent(1,k+1)
		<<" du "<<_hIcl_d_u_B[1]->GetBinContent(1,k+1)
		<<" dd "<<_hIcl_d_d_B[1]->GetBinContent(1,k+1)
		<<endl;
		n_hit1<<"event:"<<setw(3)<<k<<" : "
		<<" uu "<<_hIcl_u_u_B[0]->GetBinContent(1,k+1)
                <<" ud "<<_hIcl_u_d_B[0]->GetBinContent(1,k+1)
                <<" du "<<_hIcl_d_u_B[0]->GetBinContent(1,k+1)
                <<" dd "<<_hIcl_d_d_B[0]->GetBinContent(1,k+1)
		<<endl;
	}			

        c0 = new TCanvas("Xing_hits_check","Xing_hits_check",1000,800);
        c0->Divide(2,1);
        for(int i = 0; i <2; i++)
        {
               c0->cd(i+1);
               _hit_plots[i]->Draw("colz");
        }

	
	for(int i =0; i<2; i++)
	{
		c_yeild[i] = new TCanvas(Form("Arm%d_hits",i), Form("Arm%d_hits",i),1200,1000);
		c_yeild[i]->Divide(4,2);
		c_yeild[i]->cd(1);
		_hIcl_u_u_B[i]->Draw("colz");
		c_yeild[i]->cd(2);
		_hIcl_u_d_B[i]->Draw("colz");
		c_yeild[i]->cd(3);
		_hIcl_d_u_B[i]->Draw("colz");
                c_yeild[i]->cd(4);
		_hIcl_d_d_B[i]->Draw("colz");
                c_yeild[i]->cd(5);
		_hBgr_u_u_B[i]->Draw("colz");
                c_yeild[i]->cd(6);
		_hBgr_u_d_B[i]->Draw("colz");
                c_yeild[i]->cd(7);
		_hBgr_d_u_B[i]->Draw("colz");
                c_yeild[i]->cd(8);
		_hBgr_d_d_B[i]->Draw("colz");
	}
  	//fout.close();
  	gStyle->SetOptFit(0);
  	TCanvas *cXF = new TCanvas("cXF","cXF",1000,500);
 	_vCanvas.push_back(cXF);
  	cXF->Divide(2,1);
  	cXF->cd(1);
  	_hXF[North]->Draw();
  	_xF[North] = _hXF[North]->GetMean();
  	cXF->cd(2);
  	_hXF[South]->Draw();
  	_xF[South] = _hXF[South]->GetMean();

  	cPt = new TCanvas("cPt","cPt",1000,1000);
 	_vCanvas.push_back(cPt);
  	cPt->Divide(2,2);
  	cPt->cd(1);
  	_hPt[North][0]->Draw();
	cPt->cd(2);
        _hPt[North][1]->Draw();
  	cPt->cd(3);
  	_hPt[South][0]->Draw();
	cPt->cd(4);
        _hPt[South][1]->Draw();

  	cout << endl;
}
//-----------------------------------------------------------------------------

void JpsiAN::GetAN()
{
	cout<<"JpsiAN::GetAN()................"<<endl;	
	double n_inc_uu_B, n_inc_ud_B, n_inc_du_B,n_inc_dd_B; 
	double n_bgr_uu_B, n_bgr_ud_B, n_bgr_du_B,n_bgr_dd_B; 

	double n_inc_uu_Y, n_inc_ud_Y, n_inc_du_Y,n_inc_dd_Y;
	double n_bgr_uu_Y, n_bgr_ud_Y, n_bgr_du_Y,n_bgr_dd_Y;

	double al_inc_B[2][n_pt_bin][n_phi_bin],al_inc_Y[2][n_pt_bin][n_phi_bin];	//[arm][pt_bin][phi_bin]  
	double al_bgr_B[2][n_pt_bin][n_phi_bin],al_bgr_Y[2][n_pt_bin][n_phi_bin]; 	//[arm][pt_bin][phi_bin] 

	double eal_inc_B[2][n_pt_bin][n_phi_bin],eal_inc_Y[2][n_pt_bin][n_phi_bin];
	double eal_bgr_B[2][n_pt_bin][n_phi_bin],eal_bgr_Y[2][n_pt_bin][n_phi_bin];

	ofstream test1("test_hits.txt");

	double phi_edge[n_phi_bin];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j<n_pt_bin; j++)  //2 for arm
		{
			for (int k = 0; k < n_phi_bin; k++)
			{
				if(i==0&&j==0) phi_edge[k] = _hIcl_u_u_B[i]->GetYaxis()->GetBinLowEdge(k+1);

				n_inc_uu_B= _hIcl_u_u_B[i]->GetBinContent(j+1,k+1);	
				n_inc_ud_B= _hIcl_u_d_B[i]->GetBinContent(j+1,k+1);
				n_inc_du_B= _hIcl_d_u_B[i]->GetBinContent(j+1,k+1);
                                n_inc_dd_B= _hIcl_d_d_B[i]->GetBinContent(j+1,k+1);

				if (i == 0 && j==0)
					test1<<"event "<<setw(3)<<i<<" : "
					<<" uu "<<n_inc_uu_B
					<<" ud "<<n_inc_ud_B
					<<" du "<<n_inc_du_B
					<<" dd "<<n_inc_dd_B
					<<endl;

				n_bgr_uu_B= _hBgr_u_u_B[i]->GetBinContent(j+1,k+1);
                                n_bgr_ud_B= _hBgr_u_d_B[i]->GetBinContent(j+1,k+1);
                                n_bgr_du_B= _hBgr_d_u_B[i]->GetBinContent(j+1,k+1);
                                n_bgr_dd_B= _hBgr_d_d_B[i]->GetBinContent(j+1,k+1); 
 
				n_inc_uu_Y= _hIcl_u_u_B[i]->GetBinContent(j+1,k+1);
                                n_inc_ud_Y= _hIcl_d_u_B[i]->GetBinContent(j+1,k+1);
                                n_inc_du_Y= _hIcl_u_d_B[i]->GetBinContent(j+1,k+1);
                                n_inc_dd_Y= _hIcl_d_d_B[i]->GetBinContent(j+1,k+1);

                                n_bgr_uu_Y= _hBgr_u_u_B[i]->GetBinContent(j+1,k+1);
                                n_bgr_ud_Y= _hBgr_d_u_B[i]->GetBinContent(j+1,k+1);
                                n_bgr_du_Y= _hBgr_u_d_B[i]->GetBinContent(j+1,k+1);
                                n_bgr_dd_Y= _hBgr_d_d_B[i]->GetBinContent(j+1,k+1);
				
				al_inc_B[i][j][k] = 1. / _PB * (n_inc_uu_B +  R1B * n_inc_ud_B - R2B * n_inc_du_B - R3B * n_inc_dd_B) / (n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B);
				al_bgr_B[i][j][k] = 1. / _PB * (n_bgr_uu_B +  R1B * n_bgr_ud_B - R2B * n_bgr_du_B - R3B * n_bgr_dd_B) / (n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B);
				    cout<<" uu "<<n_inc_uu_B
                                        <<" ud "<<n_inc_ud_B
                                        <<" du "<<n_inc_du_B
                                        <<" dd "<<n_inc_dd_B
                                  	<<_PB<<"	"<<R1B<<"	"<<R2B<<"	"<<R3B
					<<endl;
				al_inc_Y[i][j][k] = 1./_PY*(n_inc_uu_Y+R1Y*n_inc_ud_Y-R2Y*n_inc_du_Y-R3Y*n_inc_dd_Y)/
                                                           (n_inc_uu_Y+R1Y*n_inc_ud_Y+R2Y*n_inc_du_Y+R3Y*n_inc_dd_Y);
                                al_bgr_Y[i][j][k] = 1./_PY*(n_bgr_uu_Y+R1Y*n_bgr_ud_Y-R2Y*n_bgr_du_Y-R3Y*n_bgr_dd_Y)/
                                                           (n_bgr_uu_Y+R1Y*n_bgr_ud_Y+R2Y*n_bgr_du_Y+R3Y*n_bgr_dd_Y);

				eal_inc_B[i][j][k]=sqrt(
				pow(1./_PB/_PB*(n_inc_uu_B +  R1B * n_inc_ud_B - R2B * n_inc_du_B - R3B * n_inc_dd_B) / (n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B)*_ePB,2) +
				pow(2./_PB*(R2B * n_inc_du_B + R3B * n_inc_dd_B)/(n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B)/(n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B),2)*n_inc_uu_B+
				pow(2./_PB*R1B*(R2B * n_inc_du_B + R3B * n_inc_dd_B)/(n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B)/(n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B),2)*n_inc_du_B+
				pow(2./_PB*R2B*(n_inc_uu_B +  R1B * n_inc_ud_B)/(n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B)/(n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B),2)*n_inc_ud_B+
				pow(2./_PB*R3B*(n_inc_uu_B +  R1B * n_inc_ud_B)/(n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B)/(n_inc_uu_B +  R1B * n_inc_ud_B + R2B * n_inc_du_B + R3B * n_inc_dd_B),2)*n_inc_dd_B);


				eal_bgr_B[i][j][k]=sqrt(
                                pow(1./_PB/_PB*(n_bgr_uu_B +  R1B * n_bgr_ud_B - R2B * n_bgr_du_B - R3B * n_bgr_dd_B) / (n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B)*_ePB,2) +
                                pow(2./_PB*(R2B * n_bgr_du_B + R3B * n_bgr_dd_B)/(n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B)/(n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B),2)*n_bgr_uu_B+
                                pow(2./_PB*R1B*(R2B * n_bgr_du_B + R3B * n_bgr_dd_B)/(n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B)/(n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B),2)*n_bgr_du_B+
                                pow(2./_PB*R2B*(n_bgr_uu_B +  R1B * n_bgr_ud_B)/(n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B)/(n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B),2)*n_bgr_ud_B+
		              	pow(2./_PB*R3B*(n_bgr_uu_B +  R1B * n_bgr_ud_B)/(n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B)/(n_bgr_uu_B +  R1B * n_bgr_ud_B + R2B * n_bgr_du_B + R3B * n_bgr_dd_B),2)*n_bgr_dd_B);

                                eal_inc_Y[i][j][k]=sqrt(
                                pow(1./_PY/_PY*(n_inc_uu_Y +  R1Y * n_inc_ud_Y - R2Y * n_inc_du_Y - R3Y * n_inc_dd_Y) / (n_inc_uu_Y +  R1Y * n_inc_ud_Y + R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y)*_ePY,2) +
                                pow(2./_PY*(R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y)/(n_inc_uu_Y +  R1Y * n_inc_ud_Y + R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y)/(n_inc_uu_Y +  R1Y * n_inc_ud_Y + R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y),2)*n_inc_uu_Y+
                                pow(2./_PY*R1Y*(R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y)/(n_inc_uu_Y +  R1Y * n_inc_ud_Y + R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y)/(n_inc_uu_Y +  R1Y * n_inc_ud_Y + R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y),2)*n_inc_du_Y+
                                pow(2./_PY*R2Y*(n_inc_uu_Y +  R1Y * n_inc_ud_Y)/(n_inc_uu_Y +  R1Y * n_inc_ud_Y + R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y)/(n_inc_uu_Y +  R1Y * n_inc_ud_Y + R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y),2)*n_inc_ud_Y+                                
				pow(2./_PY*R3Y*(n_inc_uu_Y +  R1Y * n_inc_ud_Y)/(n_inc_uu_Y +  R1Y * n_inc_ud_Y + R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y)/(n_inc_uu_Y +  R1Y * n_inc_ud_Y + R2Y * n_inc_du_Y + R3Y * n_inc_dd_Y),2)*n_inc_dd_Y);

                                eal_bgr_Y[i][j][k]=sqrt(
                                pow(1./_PY/_PY*(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y - R2Y * n_bgr_du_Y - R3Y * n_bgr_dd_Y) / (n_bgr_uu_Y +  R1Y * n_bgr_ud_Y + R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y)*_ePY,2) +
                                pow(2./_PY*(R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y)/(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y + R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y)/(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y + R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y),2)*n_bgr_uu_Y+
                                pow(2./_PY*R1Y*(R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y)/(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y + R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y)/(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y + R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y),2)*n_bgr_du_Y+
                                pow(2./_PY*R2Y*(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y)/(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y + R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y)/(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y + R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y),2)*n_bgr_ud_Y+
                                pow(2./_PY*R3Y*(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y)/(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y + R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y)/(n_bgr_uu_Y +  R1Y * n_bgr_ud_Y + R2Y * n_bgr_du_Y + R3Y * n_bgr_dd_Y),2)*n_bgr_dd_Y);
			}
		}
	}

	ofstream al_check("al_check.txt");
	for (int i = 0; i <n_phi_bin; i++)
	{
		al_check<<"bin "<<setw(3)<<i<<" : "
		<<setw(15)<< al_inc_B[0][0][i]
                <<setw(15)<< al_inc_B[1][0][i]
                <<setw(15)<< al_inc_Y[0][0][i]
                <<setw(15)<< al_inc_Y[1][0][i]

                <<setw(15)<< al_bgr_B[0][0][i]
                <<setw(15)<< al_bgr_B[1][0][i]
                <<setw(15)<< al_bgr_Y[0][0][i]
                <<setw(15)<< al_bgr_Y[1][0][i]
		<<endl;
	}
	al_check<<_PY<<"	"<<_PB<<"	"<< R1B<<"	"<<R2B<<"	"<<R3B <<endl; 
	al_check.close();

	TGraphErrors *Graph_al_inc_B[2][n_pt_bin];
	TGraphErrors *Graph_al_bgr_B[2][n_pt_bin];

        TGraphErrors *Graph_al_inc_Y[2][n_pt_bin];
        TGraphErrors *Graph_al_bgr_Y[2][n_pt_bin];

	for (int i = 0; i < 2; i++)
	{	
		for(int j = 0; j < n_pt_bin; j++)
		{
			Graph_al_inc_B[i][j] = new TGraphErrors(n_phi_bin,phi_edge,al_inc_B[i][j],0,eal_inc_B[i][j]);
			Graph_al_inc_B[i][j]->SetTitle(Form("AN_inclusive_Blue_arm%d_ptbin%d",i,j));
			Graph_al_inc_B[i][j]->GetXaxis()->SetTitle("#phi");
			Graph_al_inc_B[i][j]->GetYaxis()->SetTitle("A_{L}");

			Graph_al_bgr_B[i][j] = new TGraphErrors(n_phi_bin,phi_edge,al_bgr_B[i][j],0,eal_bgr_B[i][j]);
                        Graph_al_bgr_B[i][j]->SetTitle(Form("AN_background_Blue_arm%d_ptbin%d",i,j));
			Graph_al_bgr_B[i][j]->GetXaxis()->SetTitle("#phi");
                        Graph_al_bgr_B[i][j]->GetYaxis()->SetTitle("A_{L}");			

			Graph_al_inc_Y[i][j] = new TGraphErrors(n_phi_bin,phi_edge,al_inc_Y[i][j],0,eal_inc_Y[i][j]);
                        Graph_al_inc_Y[i][j]->SetTitle(Form("AN_inclusive_Yole_arm%d_ptbin%d",i,j));
			Graph_al_inc_Y[i][j]->GetXaxis()->SetTitle("#phi");
                        Graph_al_inc_Y[i][j]->GetYaxis()->SetTitle("A_{L}");

                        Graph_al_bgr_Y[i][j] = new TGraphErrors(n_phi_bin,phi_edge,al_bgr_Y[i][j],0,eal_bgr_Y[i][j]);
                        Graph_al_bgr_Y[i][j]->SetTitle(Form("AN_background_Yole_arm%d_ptbin%d",i,j));
			Graph_al_bgr_Y[i][j]->GetXaxis()->SetTitle("#phi");
                        Graph_al_bgr_Y[i][j]->GetYaxis()->SetTitle("A_{L}");	

		}		
	}
	TF1 *f_inc_B[2][n_pt_bin];
        TF1 *f_bgr_B[2][n_pt_bin];
	TF1 *f_inc_Y[2][n_pt_bin];
        TF1 *f_bgr_Y[2][n_pt_bin];
	
	for (int i = 0; i<2; i++)
	{
		for (int j = 0; j < n_pt_bin; j++)
		{
			f_inc_B[i][j] = new TF1(Form("inc_arm%d_bin%d_B",i,j),"[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-PI,PI);
			f_bgr_B[i][j] = new TF1(Form("bgr_arm%d_bin%d_B",i,j),"[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-PI,PI);
		
			f_inc_Y[i][j] = new TF1(Form("inc_arm%d_bin%d_Y",i,j),"[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-PI,PI);
                        f_bgr_Y[i][j] = new TF1(Form("bgr_arm%d_bin%d_Y",i,j),"[1]*TMath::Sin(x)+[0]*TMath::Cos(x)+[2]",-PI,PI);
		}
	}		
	gStyle->SetOptFit();	

	for (int i = 0; i<2; i++)
        {
	        for (int j = 0; j < n_pt_bin; j++)
                {
			Graph_al_inc_B[i][j]->Fit(f_inc_B[i][j]);
			Graph_al_bgr_B[i][j]->Fit(f_bgr_B[i][j]);
			Graph_al_inc_Y[i][j]->Fit(f_inc_Y[i][j]);
			Graph_al_bgr_Y[i][j]->Fit(f_bgr_Y[i][j]);
		}
	}
	TCanvas *c[2][n_pt_bin];


        for (int i = 0; i<2; i++)
        { 
	        for (int j = 0; j < n_pt_bin; j++)
                {
			c[i][j] = new TCanvas(Form("A_{N}_arm%d_bin%d",i,j),Form("A_{N}_arm%d_bin%d",i,j),1000,800);
                        c[i][j]->Divide(2,2);

                        c[i][j]->cd(1);
			Graph_al_inc_B[i][j]->Draw("AP");
			f_inc_B[i][j]->Draw("SAME");

			c[i][j]->cd(2);
			Graph_al_bgr_B[i][j]->Draw("AP");
			f_bgr_B[i][j]->Draw("SAME");

			c[i][j]->cd(3);
			Graph_al_inc_Y[i][j]->Draw("AP");	
			f_inc_Y[i][j]->Draw("SAME");

			c[i][j]->cd(4);
			Graph_al_bgr_Y[i][j]->Draw("AP");
			f_bgr_Y[i][j]->Draw("SAME");
		}
	}
}


void JpsiAN::ReadGoodRunList()
{
  	cout << "JpsiAN::ReadGoodRunList" << endl;

  	ifstream fin(RUNLIST);
  	int run;

  	while (fin >> run) _runList.insert(run);

  	_nRuns = _runList.size();

  	cout << "Number of good runs: " << _nRuns << "\n" << endl;

	for (int i =0; i<2; i++)
	_hit_plots[i] = new TH2D(Form("counts_arm%d",i),Form("counts_arm%d",i),120,0,119,_nRuns,0,_nRuns-1);

	int count;
	ifstream fin2("run_list.txt");
	while (fin2>>count>>run) _RunTable[run] = count;
}
//-----------------------------------------------------------------------------


void JpsiAN::ReadFillTable()
{
  	cout << "JpsiAN::ReadFillTable" << endl;

   	ifstream fin("runfilltable.list");
//  	ifstream fin("halfplane_runfilltable.list");
//  	ifstream fin("pp_goodrunfill.list");
  	int run, fill;
  	while (fin >> run >> fill) _fillTable[run] = fill;

  	fin.close();
}
//-----------------------------------------------------------------------------


void JpsiAN::GetSpinInfo()
{
  	cout << "JpsiAN::GetSpinInfo............." << endl;
	
  	TFile *fin = new TFile("GL1Pdata.root");
  	TTree *t = (TTree*)fin->Get("T");

  	cout << "Reading ... " << fin->GetName() << endl;

  	int    runNum, filNum;
	int    XingShift;

  	float  PolB, ePolBstat, ePolBsyst;
  	float  PolY, ePolYstat, ePolYsyst;

  	short spinB[XINGs];
  	short spinY[XINGs];
  	float bbc_in[XINGs];


  	t->SetBranchAddress("RunNumber",  &runNum);
  	t->SetBranchAddress("FillNumber", &filNum);
	t->SeaBranchAddress("XingShift",  &XingShift);		

  	t->SetBranchAddress("PolB",       &PolB);
  	t->SetBranchAddress("ePolBstat",  &ePolBstat);
  	t->SetBranchAddress("ePolBsyst",  &ePolBsyst);
  	t->SetBranchAddress("PolY",       &PolY);
  	t->SetBranchAddress("ePolYstat",  &ePolYstat);
  	t->SetBranchAddress("ePolYsyst",  &ePolYsyst);

  	t->SetBranchAddress("SpinB",      spinB);
  	t->SetBranchAddress("SpinY",      spinY);
  	t->SetBranchAddress("BBCin",      bbc_in);

	ofstream run_xing_lumi("run_xing_lumi.txt");
	
	int xing;
  	int fill=0, oldFill=0;
  	bool exists = false;
  	for (int i=0; i<t->GetEntries(); i++)
    	{
      		t->GetEntry(i);
      		if (_runList.count(runNum)==0) continue;
      		if (PolB == 0 || PolY == 0) cout << runNum << endl;

      		fill = filNum;
      		exists = false;
      		for (int k=0; k<_nFills; k++)
		if (fill == _Fill[k]) {exists = true; break;}

      		if (!exists)
		{
	  		cout << fill << "  " << runNum << endl;
	  		continue;
		}

      		if (fill != oldFill)
		{
	  		oldFill = fill;

	  		_PolB[fill]  = PolB;
	  		_PolY[fill]  = PolY;
	  		_ePolB[fill] =
	    		float(sqrt(ePolBstat*ePolBstat + ePolBsyst*ePolBsyst));
	  		_ePolY[fill] =
	    		float(sqrt(ePolYstat*ePolYstat + ePolYsyst*ePolYsyst));
	  		_spinB[fill] = new short[XINGs];
	  		_spinY[fill] = new short[XINGs];

	  		_BBCinB[fill] = new double[XINGs];
	  		_BBCinY[fill] = new double[XINGs];

	  		for (int j=0; j<XINGs; j++)
	    		{
				if(j >= XingShift) xing = j-XingShift;
				else xing = (j+120-XingShift); 
	      			_spinB[fill][xing] = spinB[j];
	      			_spinY[fill][xing] = spinY[j];

	      			_BBCinB[fill][j] = 0.;
	      			_BBCinY[fill][j] = 0.;
	    		}
		}

		for (int j=0; j<XINGs; j++)
		{
			if(j >= XingShift) xing = j-XingShift;
                        else xing = (j+120-XingShift);
	    		_BBCinB[fill][xing] += double(bbc_in[j]);
	    		_BBCinY[fill][xing] += double(bbc_in[j]);
			run_xing_lumi<<runNum<<"	"<<setw(5)<<j<<setw(15)<<bbc_in[j]<<endl;	

		}
	}    

  	fin->Close();

  	cout << _nFills << "  " << _BBCinB.size() << "  " << _PolB.size() << endl;

  	double *pol[2];
  	double *pole[2];
  	pol[Blue]  = new double[_nFills];
  	pol[Yelo]  = new double[_nFills];
  	pole[Blue] = new double[_nFills];
 	pole[Yelo] = new double[_nFills];

  	for (int i=0; i<_nFills; i++)
    	{
      		fill = int(_Fill[i]);
      		pol [Blue][i] =  _PolB.find(fill)->second;
      		pol [Yelo][i] =  _PolY.find(fill)->second;
      		pole[Blue][i] = _ePolB.find(fill)->second;
      		pole[Yelo][i] = _ePolY.find(fill)->second;

      		cout << fill << "   " << pol[Blue][i] << "   " << pol[Yelo][i] << "	"<<pole[Blue][i]<<"	"<<pole[Yelo][i]<<endl;
    	}
  	cout << endl;


  	delete [] pol [Blue];
  	delete [] pole[Blue];
  	delete [] pol [Yelo];
  	delete [] pole[Yelo];
}
//-----------------------------------------------------------------------------


void JpsiAN::RelLumi()
{
  	cout << "JpsiAN::RelLumi" << endl;
  	int fill= -99;
  	short spinB, spinY;
  	double L_up_up, L_up_dn, L_dn_up, L_dn_dn, L_total;
  	double weight[_nFills];
	ofstream lumi_fill("lumi_fill.txt");
	
  	double bbcin;
	_L_uuB = 0;
	_L_udB = 0;
	_L_duB = 0;
	_L_ddB = 0;

	_L_uuY = 0;
        _L_udY = 0;
        _L_duY = 0;
        _L_ddY = 0;
  	cout << _nFills << endl;
  	for (int i=0; i<_nFills; i++)
    	{
      		fill = int(_Fill[i]);
      		weight[i] = 0.;
      		L_up_up = 0.;
      		L_up_dn = 0.;
      		L_dn_up = 0.;
      		L_dn_dn = 0.;

      	for (int xing=0; xing<XINGs; xing++)
		{
	  		spinB = _spinB.find(fill)->second[xing];
	  		spinY = _spinY.find(fill)->second[xing];
	  		bbcin = _BBCinB.find(fill)->second[xing];

	  		if (spinB == UP && spinY == UP) L_up_up += bbcin;
	  		if (spinB == UP && spinY == DN) L_up_dn += bbcin;
	  		if (spinB == DN && spinY == UP) L_dn_up += bbcin;
	  		if (spinB == DN && spinY == DN) L_dn_dn += bbcin;

	  		weight[i] =weight[i] + bbcin; 
			
		}
		_L_uuB = _L_uuB + L_up_up;
		_L_udB = _L_udB + L_up_dn;
		_L_duB = _L_duB + L_dn_up;
		_L_ddB = _L_ddB + L_dn_dn;

		_L_uuY = _L_uuY + L_up_up; 
		_L_udY = _L_udY + L_dn_up;
    	        _L_duY = _L_duY + L_up_dn;
        	_L_ddY = _L_ddY + L_dn_dn;   
 
        }

	cout <<"test for luminosity "<<_L_uuB <<"	"<<_L_udB<<"	"<<_L_duB<<"	"<<_L_ddB<<endl;

	R1B = _L_uuB/_L_udB;
	R2B = _L_uuB/_L_duB;
	R3B = _L_uuB/_L_ddB;

	R1Y = _L_uuY/_L_udY;
        R2Y = _L_uuY/_L_duY;
        R3Y = _L_uuY/_L_ddY;
	L_total = 0;
	_PB=0;
        _ePB=0;
        _PY=0;
        _ePY=0;
	ofstream polar_fill("polar_fill.txt");
	
	double polB, polY, epolB, epolY;
	for (int i=0; i<_nFills; i++)
	{
		fill = int(_Fill[i]);
		lumi_fill<<setw(6)<<"fill: "<<fill<<setw(8)<< "weight: "<<weight[i]<<"  polB  "<<_PolB.find(fill)->second<<"polY  "<<_PolY.find(fill)->second<<endl;
		polB = _PolB.find(fill)->second;
		polY = _PolY.find(fill)->second;
		epolB= _ePolB.find(fill)->second;
		epolY= _ePolY.find(fill)->second;

		polar_fill << setw(8)<<weight[i]<<setw(8)<<"polB "<<setw(8)<<polB<<setw(8)<<"polY "<<setw(8)<<polY<<setw(8)<<"epolB "<<setw(8)<<epolB<<setw(8)<<"epolY "<<setw(10)<<epolY<<setw(10)<<_PY<<endl;

		_PB = _PB+ weight[i] * polB;
		_ePB=_ePB+ weight[i] * epolB;
		_PY = _PY+ weight[i] * polY;
                _ePY=_ePY+ weight[i] * epolY;	
		L_total += weight[i]; 
	}
	

	_PB = _PB/L_total;
	_ePB = _ePB/L_total;
	_PY = _PY/L_total;
        _ePY = _ePY/L_total;
	cout << "pol for Blue: "<< _PB << "	pol error for blue: " <<_ePB<<endl;
	cout << "pol for yellow: "<< _PY << "     pol error for yellow: " <<_ePY<<endl;

}
//-----------------------------------------------------------------------------


void JpsiAN::SetArrays()
{
  	cout << "JpsiAN::SetArrays" << endl;

  	_eFill = new double[_nFills];
  	for (int i=0; i<_nFills; i++) _eFill[i] = 0.;

  		for (int i=0; i<2; i++)
    		{
      			_R1[i] = new double[_nFills];
      			_R2[i] = new double[_nFills];
      			_R3[i] = new double[_nFills];

      			_eR1[i] = new double[_nFills];
      			_eR2[i] = new double[_nFills];
      			_eR3[i] = new double[_nFills];

      				for (int j=0; j<_nFills; j++)
				{
	  				_R1[i][j] = 0.;
	  				_R2[i][j] = 0.;
	  				_R3[i][j] = 0.;

	  				_eR1[i][j] = 0.;
	  				_eR2[i][j] = 0.;
	  				_eR3[i][j] = 0.;

				}

    		}		

}
//-----------------------------------------------------------------------------


void JpsiAN::SetTree()
{
  	_t->SetBranchAddress("Run_Number",        &Run_Number);
  	_t->SetBranchAddress("Evt_Number",        &Evt_Number);
  	_t->SetBranchAddress("Evt_Nmu",           &Evt_Nmu);
  	_t->SetBranchAddress("Evt_bbcZ",          &Evt_bbcZ);
  	_t->SetBranchAddress("Evt_vtxchi2",       &Evt_vtxchi2);
  	_t->SetBranchAddress("Tr0_DDG0",          &Tr0_DDG0);
  	_t->SetBranchAddress("Tr0_DG0",           &Tr0_DG0);
  	_t->SetBranchAddress("Tr0_DS3",           &Tr0_DS3);
  	_t->SetBranchAddress("Tr0_idchi2",        &Tr0_idchi2);
  	_t->SetBranchAddress("Tr0_lastgap",       &Tr0_lastgap);
  	_t->SetBranchAddress("Tr0_ntrhits",       &Tr0_ntrhits);
  	_t->SetBranchAddress("Tr1_nidhits",       &Tr1_nidhits);
  	_t->SetBranchAddress("Tr0_dca_r",         &Tr0_dca_r);  
  	_t->SetBranchAddress("Tr0_px",            &Tr0_px);
  	_t->SetBranchAddress("Tr0_py",            &Tr0_py);
  	_t->SetBranchAddress("Tr0_pz",            &Tr0_pz);
  	_t->SetBranchAddress("Tr0_trhits",        &Tr0_trhits);
  	_t->SetBranchAddress("Tr1_DDG0",          &Tr1_DDG0);
  	_t->SetBranchAddress("Tr1_DG0",           &Tr1_DG0);
  	_t->SetBranchAddress("Tr1_DS3",           &Tr1_DS3);
  	_t->SetBranchAddress("Tr1_ntrhits",       &Tr1_ntrhits);
  	_t->SetBranchAddress("Tr1_lastgap",       &Tr1_lastgap);
  	_t->SetBranchAddress("Tr1_dca_r",         &Tr1_dca_r);
  	_t->SetBranchAddress("Tr1_idchi2",        &Tr1_idchi2);
  	_t->SetBranchAddress("Tr0_nidhits",       &Tr0_nidhits);
  	_t->SetBranchAddress("Tr1_px",            &Tr1_px);
  	_t->SetBranchAddress("Tr1_py",            &Tr1_py);
  	_t->SetBranchAddress("Tr1_pz",            &Tr1_pz);
  	_t->SetBranchAddress("charge",            &charge);
  	_t->SetBranchAddress("mass",              &mass);
  	_t->SetBranchAddress("p",                 &p);
  	_t->SetBranchAddress("pT",                &pT);
  	_t->SetBranchAddress("rapidity",          &rapidity);
  	_t->SetBranchAddress("x1",                &x1);
  	_t->SetBranchAddress("x2",                &x2);
  	_t->SetBranchAddress("xF",                &xF);
  	_t->SetBranchAddress("SpinX_ID",          &SpinX_ID);
  	_t->SetBranchAddress("px",		    &px);
  	_t->SetBranchAddress("py",                &py);
  	_t->SetBranchAddress("pz",                &pz);
  	_t->SetBranchAddress("lvl1_trigscaled",   &lvl1_trigscaled);
  	_t->SetBranchAddress("same_event",        &same_event);
	_t->SetBranchAddress("beamclk0",	  &beamclk0);
}
