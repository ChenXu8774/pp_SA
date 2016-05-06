#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TSystem.h>
#include <TMath.h>
#include <TLorentzVector.h>
	
//#include <RpcSingleMuonContainer.h>
//#include <RpcDiMuonContainer.h>
//#include <DiMuonContainer.h>
//#include <RunHeader.h>
//#include <TrigLvl1.h>
//#include <SyncObject.h>
	
static const int NCROSS=120;
	
using namespace std;
	
int prepare_ntuples
(
 char *f_in_list_name = "picoDST.list", char *f_out_root_name = "Run15_pp_goodrun.root"
)
{
//	gSystem->Load("/phenix/spin/phnxsp01/yuhw/CVS_PHENIX/install/lib/libpicodst_object.so");
	gSystem->Load("libpicodst_object.so");
	//input chain
	TChain chain_in_dimuon("T");
	TChain chain_in_header("T1");
	std::ifstream h_inf(f_in_list_name);
	if(!h_inf.is_open())
	{   
		cout << "The File " << f_in_list_name << " can not be opened!" << endl;
	}   
	while(!h_inf.eof())
	{   
		char filename[1000];
		h_inf.getline(filename,1000);
		if(h_inf.eof()) continue;
		chain_in_dimuon.Add(filename);
		chain_in_header.Add(filename);
	}
	DiMuonContainer *dimuoncontainer = NULL;
	chain_in_dimuon.SetBranchAddress("DST/DiMuonContainer", &dimuoncontainer);
	TrigLvl1 *triglvl1 = NULL;
	chain_in_dimuon.SetBranchAddress("DST/TrigLvl1",&triglvl1);
	SyncObject * syncobject = NULL;
	chain_in_dimuon.SetBranchAddress("DST/Sync",&syncobject);
	RunHeader *runheader = NULL;
	chain_in_header.SetBranchAddress("RUN/RunHeader",&runheader);
	RpcSingleMuonContainer *rpcsinglemuoncontainer = NULL;
	RpcDiMuonContainer * rpcdimuoncontainer = NULL;
	chain_in_dimuon.SetBranchAddress("DST/RpcDiMuonContainer",&rpcdimuoncontainer);

	//output root file
	TFile * f_out = new TFile(f_out_root_name,"RECREATE");
	TTree * t_out = new TTree("T","DiMuon NTUPLE");

	float Run_Number = -999;
	bool  same_event = false;
	float Evt_Number = -999;
	float Evt_Nmu = -999;
	float Evt_bbcCentrality = -999;
	float Evt_bbcZ = -999;
	float Evt_vtxchi2 = -999;
	float Evt_vtxoor = -999;
	float Tr0_DDG0 = -999;
	float Tr0_DG0 = -999;
	float Tr0_DS3 = -999;
        float Tr0_idchi2 = -999;
        float Tr0_trchi2 = -999;
        float Tr0_idhits = -999;
        float Tr0_nidhits = -999;
        float Tr0_px = -999;
        float Tr0_py = -999;
        float Tr0_pz = -999;
	float Tr0_trhits = -999;
	float Tr0_ntrhits = -999;
	float Tr1_DDG0 = -999;
        float Tr1_DG0 = -999;
        float Tr1_DS3 = -999;
        float Tr1_idchi2 = -999;
        float Tr1_trchi2 = -999;
        float Tr1_idhits = -999;
        float Tr1_nidhits = -999;
        float Tr1_px = -999;
        float Tr1_py = -999;
        float Tr1_pz = -999;
        float Tr1_trhits = -999;
        float Tr1_ntrhits = -999;
        float charge = -999;
        float costhCS = -999;
	float phiHX = -999;
	float phiGJ = -999;
	float phiCS = -999;
        float mass = -999;
        float mass_fvtxmutr = -999;
        float X0 = -999;
        float Y0 = -999;
        float Z0 = -999;
	float R0 = -999;
        float p = -999;
        float pT = -999;
        float px = -999;
        float py = -999;
        float pz = -999;
        float rapidity = -999;
        float x1 = -999;
        float x2 = -999;
        float xF = -999;
        UInt_t lvl1_trigraw = 0;
        UInt_t lvl1_triglive = 0;
        UInt_t lvl1_trigscaled = 0;
        unsigned int SpinX_ID = -999;
	unsigned int beamclk0 = 0;
        unsigned int beamclk1 = 0;
        float Evt_fvtxX = -999;
        float Evt_fvtxY = -999;
        float Evt_fvtxZ = -999;

        float Tr0_Rpc1St1Time = -999;
        float Tr1_Rpc1St1Time = -999;
        float Tr0_Rpc3St3Time = -999;
        float Tr1_Rpc3St3Time = -999;
	
        short Tr0_lastgap = -999;
        short Tr1_lastgap = -999;


	//fvtx info
	float dca_r = -999;
        float Tr0_dca_r = -999;
        float Tr1_dca_r = -999;
        float Tr0_dca_phi = -999;
        float Tr1_dca_phi = -999;
        float Tr0_dca_z = -999;
        float Tr1_dca_z = -999;
        float Tr0_chi2_fvtxmutr = -999;
        float Tr1_chi2_fvtxmutr = -999;
        float Tr0_dr_fvtx = -999;
        float Tr1_dr_fvtx = -999;
        float Tr0_dtheta_fvtx = -999;
        float Tr1_dtheta_fvtx = -999;
        float Tr0_dphi_fvtx = -999;
        float Tr1_dphi_fvtx = -999;
	
	
        t_out->Branch("lvl1_trigraw",&lvl1_trigraw,"lvl1_trigraw/i");
        t_out->Branch("lvl1_triglive",&lvl1_triglive,"lvl1_triglive/i");
        t_out->Branch("lvl1_trigscaled",&lvl1_trigscaled,"lvl1_trigscaled/i");
	
        t_out->Branch("Evt_fvtxX",&Evt_fvtxX,"Evt_fvtxX/F");
        t_out->Branch("Evt_fvtxY",&Evt_fvtxY,"Evt_fvtxY/F");
        t_out->Branch("Evt_fvtxZ",&Evt_fvtxZ,"Evt_fvtxZ/F");
  //      t_out->Branch("Tr0_Rpc1St1Time",&Tr0_Rpc1St1Time,"Tr0_Rpc1St1Time/F");
    //    t_out->Branch("Tr1_Rpc1St1Time",&Tr1_Rpc1St1Time,"Tr1_Rpc1St1Time/F");
      //  t_out->Branch("Tr0_Rpc3St3Time",&Tr0_Rpc3St3Time,"Tr0_Rpc3St3Time/F");
        //t_out->Branch("Tr1_Rpc3St3Time",&Tr1_Rpc3St3Time,"Tr1_Rpc3St3Time/F");

        t_out->Branch("Tr0_lastgap",&Tr0_lastgap,"Tr0_lastgap/S");
        t_out->Branch("Tr1_lastgap",&Tr1_lastgap,"Tr1_lastgap/S");
        t_out->Branch("dca_r",&dca_r,"dca_r/F");
        t_out->Branch("Tr0_dca_r",&Tr0_dca_r,"Tr0_dca_r/F");
        t_out->Branch("Tr1_dca_r",&Tr1_dca_r,"Tr1_dca_r/F");
        t_out->Branch("Tr0_dca_phi",&Tr0_dca_phi,"Tr0_dca_phi/F");
        t_out->Branch("Tr1_dca_phi",&Tr1_dca_phi,"Tr1_dca_phi/F");
        t_out->Branch("Tr0_dca_z",&Tr0_dca_z,"Tr0_dca_z/F");
        t_out->Branch("Tr1_dca_z",&Tr1_dca_z,"Tr1_dca_z/F");
        t_out->Branch("Tr0_chi2_fvtxmutr",&Tr0_chi2_fvtxmutr,"Tr0_chi2_fvtxmutr/F");
        t_out->Branch("Tr1_chi2_fvtxmutr",&Tr1_chi2_fvtxmutr,"Tr1_chi2_fvtxmutr/F");
        t_out->Branch("Tr0_dr_fvtx",&Tr0_dr_fvtx,"Tr0_dr_fvtx/F");
        t_out->Branch("Tr1_dr_fvtx",&Tr1_dr_fvtx,"Tr1_dr_fvtx/F");
        t_out->Branch("Tr0_dtheta_fvtx",&Tr0_dtheta_fvtx,"Tr0_dtheta_fvtx/F");
        t_out->Branch("Tr1_dtheta_fvtx",&Tr1_dtheta_fvtx,"Tr1_dtheta_fvtx/F");
        t_out->Branch("Tr0_dphi_fvtx",&Tr0_dphi_fvtx,"Tr0_dphi_fvtx/F");
        t_out->Branch("Tr1_dphi_fvtx",&Tr1_dphi_fvtx,"Tr1_dphi_fvtx/F");
	
        t_out->Branch("Run_Number",&Run_Number,"Run_Number/F");
        t_out->Branch("same_event",&same_event,"same_event/O");
        t_out->Branch("Evt_Number",&Evt_Number,"Evt_Number/F");
        t_out->Branch("Evt_Nmu",&Evt_Nmu,"Evt_Nmu/F");
        t_out->Branch("Evt_bbcCentrality",&Evt_bbcCentrality,"Evt_bbcCentrality/F");
        t_out->Branch("Evt_bbcZ",&Evt_bbcZ,"Evt_bbcZ/F");
        t_out->Branch("Evt_vtxchi2",&Evt_vtxchi2,"Evt_vtxchi2/F");
        t_out->Branch("Evt_vtxoor",&Evt_vtxoor,"Evt_vtxoor/F");
        t_out->Branch("Tr0_DDG0",&Tr0_DDG0,"Tr0_DDG0/F");
        t_out->Branch("Tr0_DG0",&Tr0_DG0,"Tr0_DG0/F");
        t_out->Branch("Tr0_DS3",&Tr0_DS3,"Tr0_DS3/F");
        t_out->Branch("Tr0_idchi2",&Tr0_idchi2,"Tr0_idchi2/F");
        t_out->Branch("Tr0_trchi2",&Tr0_trchi2,"Tr0_trchi2/F");
        t_out->Branch("Tr0_idhits",&Tr0_idhits,"Tr0_idhits/F");
        t_out->Branch("Tr0_nidhits",&Tr0_nidhits,"Tr0_nidhits/F");
        t_out->Branch("Tr0_px",&Tr0_px,"Tr0_px/F");
        t_out->Branch("Tr0_py",&Tr0_py,"Tr0_py/F");
        t_out->Branch("Tr0_pz",&Tr0_pz,"Tr0_pz/F");
        t_out->Branch("Tr0_trhits",&Tr0_trhits,"Tr0_trhits/F");
        t_out->Branch("Tr0_ntrhits",&Tr0_ntrhits,"Tr0_ntrhits/F");
        t_out->Branch("Tr1_DDG0",&Tr1_DDG0,"Tr1_DDG0/F");
        t_out->Branch("Tr1_DG0",&Tr1_DG0,"Tr1_DG0/F");
        t_out->Branch("Tr1_DS3",&Tr1_DS3,"Tr1_DS3/F");
        t_out->Branch("Tr1_idchi2",&Tr1_idchi2,"Tr1_idchi2/F");
        t_out->Branch("Tr1_trchi2",&Tr1_trchi2,"Tr1_trchi2/F");
        t_out->Branch("Tr1_idhits",&Tr1_idhits,"Tr1_idhits/F");
        t_out->Branch("Tr1_nidhits",&Tr1_nidhits,"Tr1_nidhits/F");
        t_out->Branch("Tr1_px",&Tr1_px,"Tr1_px/F");
        t_out->Branch("Tr1_py",&Tr1_py,"Tr1_py/F");
        t_out->Branch("Tr1_pz",&Tr1_pz,"Tr1_pz/F");
        t_out->Branch("Tr1_trhits",&Tr1_trhits,"Tr1_trhits/F");
        t_out->Branch("Tr1_ntrhits",&Tr1_ntrhits,"Tr1_ntrhits/F");
        t_out->Branch("charge",&charge,"charge/F");
        t_out->Branch("costhCS",&costhCS,"costhCS/F");
	t_out->Branch("phiHX",&phiHX,"phiHX/F");
	t_out->Branch("phiCS",&phiCS,"phiCS/F");
	t_out->Branch("phiGJ",&phiGJ,"phiGJ/F");

        t_out->Branch("mass",&mass,"mass/F");
        t_out->Branch("mass_fvtxmutr",&mass_fvtxmutr,"mass_fvtxmutr/F");
        t_out->Branch("X0",&X0,"X0/F");
        t_out->Branch("Y0",&Y0,"Y0/F");
        t_out->Branch("Z0",&Z0,"Z0/F");
	t_out->Branch("R0",&R0,"R0/F");
        t_out->Branch("p",&p,"p/F");
        t_out->Branch("pT",&pT,"pT/F");
        t_out->Branch("px",&px,"px/F");
        t_out->Branch("py",&py,"py/F");
        t_out->Branch("pz",&pz,"pz/F");
        t_out->Branch("rapidity",&rapidity,"rapidity/F");
        t_out->Branch("x1",&x1,"x1/F");
        t_out->Branch("x2",&x2,"x2/F");
        t_out->Branch("xF",&xF,"xF/F");
        t_out->Branch("SpinX_ID",&SpinX_ID,"SpinX_ID/i");
        t_out->Branch("beamclk0",&beamclk0,"beamclk0/i");
        t_out->Branch("beamclk1",&beamclk1,"beamclk1/i");

        Long64_t Nentries_dimuons = chain_in_dimuon.GetEntries();
        cout<<"Nentries_dimuons = "<<Nentries_dimuons<<endl;
	
        Long64_t Nentries_header = chain_in_header.GetEntries();
        cout<<"Nentries_header = "<<Nentries_header<<endl;
        Long64_t ientry_header = 0;
	
        Int_t oldTreeNumber = -999;
        for(Int_t ientry_dimuon = 0; ientry_dimuon < Nentries_dimuons; ientry_dimuon++)
        {   
                if(ientry_dimuon%100000==0)cout<<"processing "<<100.0*ientry_dimuon/Nentries_dimuons<<"%"<<endl;
	
                if((chain_in_dimuon.LoadTree(ientry_dimuon))<0) break;

                chain_in_dimuon.GetEntry(ientry_dimuon);
	
                ientry_header = chain_in_dimuon.GetTreeNumber();
                chain_in_header.GetEntry(ientry_header);
                if(runheader) Run_Number = runheader->get_RunNumber();
                else Run_Number = -999;

                if(dimuoncontainer==NULL) continue;
                UInt_t nDiMuons = dimuoncontainer->get_nDiMuons();
                if(nDiMuons<1) continue;
                Evt_Nmu = nDiMuons;
                //Evt_Number = syncobject->EventNumber();
                Evt_Number = 0;
	
                Evt_bbcCentrality = dimuoncontainer->get_Evt_Cent();
                Evt_bbcZ = dimuoncontainer->get_Evt_bbcZ();

                Evt_fvtxX = dimuoncontainer->get_Evt_fvtxX();
                Evt_fvtxY = dimuoncontainer->get_Evt_fvtxY();
                Evt_fvtxZ = dimuoncontainer->get_Evt_fvtxZ();

	for (Int_t i1=0;i1<nDiMuons;i1++)
        { 
/*		if(rpcdimuoncontainer)
                        {
                               if(rpcdimuoncontainer->get_nRpcDiMuons()==nDiMuons)
                                {
                                        RpcDiMuon *rpcdimuon = rpcdimuoncontainer->get_RpcDiMuon(i1);

                                        Tr0_Rpc1St1Time = rpcdimuon->get_Tr0_Rpc1St1Time();
                                        Tr0_Rpc3St3Time = rpcdimuon->get_Tr0_Rpc3St3Time();
                                        Tr1_Rpc1St1Time = rpcdimuon->get_Tr1_Rpc1St1Time();
                                        Tr1_Rpc3St3Time = rpcdimuon->get_Tr1_Rpc3St3Time();
                                }
                        }
                        else
                        {
                                Tr0_Rpc1St1Time = -999; 
                                Tr0_Rpc3St3Time = -999; 
                                Tr1_Rpc1St1Time = -999; 
                                Tr1_Rpc3St3Time = -999; 
                        }*/

                        DiMuon *dimuon = (DiMuon*)dimuoncontainer->get_DiMuon(i1);

                        same_event = dimuon->get_same_event();
	
                        Tr0_lastgap = dimuon->get_Tr0_lastgap();
                        Tr1_lastgap = dimuon->get_Tr1_lastgap();
	
                        dca_r = dimuon->get_dca_r();
                        Tr0_dca_r = dimuon->get_Tr0_dca_r();
                        Tr1_dca_r = dimuon->get_Tr1_dca_r();
                        Tr0_dca_phi = dimuon->get_Tr0_dca_phi();
                        Tr1_dca_phi = dimuon->get_Tr1_dca_phi();
                        Tr0_dca_z = dimuon->get_Tr0_dca_z();
                        Tr1_dca_z = dimuon->get_Tr1_dca_z();
                        Tr0_chi2_fvtxmutr = dimuon->get_Tr0_chi2_fvtxmutr();
                        Tr1_chi2_fvtxmutr = dimuon->get_Tr1_chi2_fvtxmutr();
                        Tr0_dr_fvtx = dimuon->get_Tr0_dr_fvtx();
                        Tr1_dr_fvtx = dimuon->get_Tr1_dr_fvtx();
                        Tr0_dtheta_fvtx = dimuon->get_Tr0_dtheta_fvtx();
                        Tr1_dtheta_fvtx = dimuon->get_Tr1_dtheta_fvtx();
                        Tr0_dphi_fvtx = dimuon->get_Tr0_dphi_fvtx();
                        Tr1_dphi_fvtx = dimuon->get_Tr1_dphi_fvtx();

                        Evt_vtxchi2 = dimuon->get_Evt_vtxchi2();
                        Evt_vtxoor = dimuon->get_Evt_vtxoor();

                        Tr0_DDG0 = dimuon->get_Tr0_DDG0();
                        Tr0_DG0 = dimuon->get_Tr0_DG0();
                        Tr0_DS3 = dimuon->get_Tr0_DS3();
                        Tr0_idchi2 = dimuon->get_Tr0_idchi2();
                        Tr0_trchi2 = dimuon->get_Tr0_trchi2();
                        Tr0_idhits = dimuon->get_Tr0_idhits();
                        Tr0_nidhits = dimuon->get_Tr0_nidhits();
                        Tr0_px = dimuon->get_Tr0_px();
                        Tr0_py = dimuon->get_Tr0_py();
                        Tr0_pz = dimuon->get_Tr0_pz();
                        Tr0_trhits = dimuon->get_Tr0_trhits();
                        Tr0_ntrhits = dimuon->get_Tr0_ntrhits();
                        Tr1_DDG0 = dimuon->get_Tr1_DDG0();
                        Tr1_DG0 = dimuon->get_Tr1_DG0();
                        Tr1_DS3 = dimuon->get_Tr1_DS3();
                        Tr1_idchi2 = dimuon->get_Tr1_idchi2();
                        Tr1_trchi2 = dimuon->get_Tr1_trchi2();
                        Tr1_idhits = dimuon->get_Tr1_idhits();
                        Tr1_nidhits = dimuon->get_Tr1_nidhits();
                        Tr1_px = dimuon->get_Tr1_px();
                        Tr1_py = dimuon->get_Tr1_py();
                        Tr1_pz = dimuon->get_Tr1_pz();
                        Tr1_trhits = dimuon->get_Tr1_trhits();
                        Tr1_ntrhits = dimuon->get_Tr1_ntrhits();
                        charge = dimuon->get_charge();
                        costhCS = dimuon->get_costhCS();
			phiHX = dimuon->get_phiHX();
			phiGJ = dimuon->get_phiGJ();
			phiCS = dimuon->get_phiCS();
                        mass = dimuon->get_mass();
                        mass_fvtxmutr = dimuon->get_mass_fvtxmutr();
                        X0 = dimuon->get_X0();
                        Y0 = dimuon->get_Y0();
                        Z0 = dimuon->get_Z0();
			R0 = TMath::Sqrt(X0*X0+Y0*Y0);

			px = dimuon->get_px();
			py = dimuon->get_py();
			pz = dimuon->get_pz();
			pT = dimuon->get_pT();

			p = TMath::Sqrt(pT*pT+pz*pz);
                        rapidity = dimuon->get_rapidity();

                        x1 = mass/200*TMath::Exp(rapidity);
                        x2 = mass/200*TMath::Exp(-rapidity);
                        xF = x1-x2;

                        if(triglvl1)
                        {
                                lvl1_trigraw    = triglvl1->get_lvl1_trigraw();
                                //cout<<"lvl1_trigraw= "<<lvl1_trigraw<<endl;
                                lvl1_triglive   = triglvl1->get_lvl1_triglive();
                                lvl1_trigscaled = triglvl1->get_lvl1_trigscaled();
                                SpinX_ID = triglvl1->get_lvl1_clock_cross();
                                beamclk0 = triglvl1->get_lvl1_beam_clk(0);
                                beamclk1 = triglvl1->get_lvl1_beam_clk(1);
                        }
                        else
                        {
                                lvl1_trigraw    = 0;
                                lvl1_triglive   = 0;
                                lvl1_trigscaled = 0;
                                SpinX_ID = -999;
                                beamclk0 = -999;
                                beamclk1 = -999;
                        }

                        t_out->Fill();
                }
        }
	
        f_out->cd();
        t_out->Write();
        f_out->Close();
        return 1;
}			



