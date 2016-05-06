// $Id: Test_SpinAna.C,v 1.1 2013/05/12 07:07:36 jinhuang Exp $

/*!
 * \file Test_SpinAna.C
 * \brief Example how to draw the result of spin analyzer
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2013/05/12 07:07:36 $
 */

#include "SaveCanvas.C"
#include "SetOKStyle.C"
#include <string.h>
#include <cassert>
#include <fstream>
#include <iostream>
using namespace std;

void
Test_SpinAna()
{
  SetOKStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  TVirtualFitter::SetDefaultFitter("Minuit2");
  gSystem->Load("/gpfs/mnt/gpfs02/phenix/spin/spin2/chenxu/spinAnalyzer/install/lib/libspin_analyzer.so");
//  gSystem->Load("/phenix/spin/phnxsp01/yuhw/CVS_PHENIX/install/lib/libspin_analyzer.so");

//  MakeGoodRunList();
//  CrossingCheck();
//  AsymmetryCheck();
//  CrossingCheck("/direct/phenix+spin/phnxsp01/mxliu/run13/ana_pp510//dimuon_run13pp510-short.list_hist.root","RelLumi");
//  AsymmetryCheck("/direct/phenix+spin/phnxsp01/mxliu/run13/ana_pp510//dimuon_run13pp510-short.list_hist.root","DiMuonPT","A_LL",2);
//  AsymmetryCheck(
//      "/phenix/u/jinhuang/links/spin_data/data/Run13_SpinAna1/SngMuons_hist.lst.root",
//      "SngMuonPT", "A_L_Yellow", 2);

//  ////////////////////////////////////////////////////////////
//  // Run13 Dimuons
//  ////////////////////////////////////////////////////////////
//
//  const int n_hists = 22;
//  const char * hist_names[] =
//    { "RelLumi", "DiMuonB2BQm", "DiMuonB2BQp", "DiMuonB2B", "DiMuonMassN",
//        "DiMuonMassQmN", "DiMuonMassQmS", "DiMuonMassQm", "DiMuonMassQpN",
//        "DiMuonMassQpS", "DiMuonMassQp", "DiMuonMassS", "DiMuonMass",
//        "DiMuonPTN", "DiMuonPTS", "DiMuonPT", "DiMuonPTmN", "DiMuonPTmS",
//        "DiMuonPTm", "DiMuonPTpN", "DiMuonPTpS", "DiMuonPTp" };
//  //RelLumi have to be the first one for asymmetry calculation
////  TString list = "/phenix/u/jinhuang/links/spin_data/data/Run13_SpinAna1/DiMuons_hist_10.lst";
//  TString list = "/phenix/u/jinhuang/links/spin_data/data/Run13_SpinAna1/DiMuons_hist.lst";
//  MergeHist(list,"spinQA.txt", n_hists, hist_names);
//  CrossingCheck(list+".root","RelLumi");
//  ////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Run13 Sngmuons
////////////////////////////////////////////////////////////

  const int n_hists = 5;
  const char * hist_names[] =
    { "RelLumi", // Histogram container for RelLumi_RUN_0000386773
        "_sah_pT_os_lsb_n", //  Histogram container for SngHadronEta_RUN_0000386773
        "_sah_pT_os_lsb_s", // Histogram container for SngHadronEtam_RUN_0000386773
        "_sah_pT_os_sig_n", // Histogram container for SngHadronEtap_RUN_0000386773
        "_sah_pT_os_sig_s", //  Histogram container for SngHadronPTN_RUN_0000386773
     /*   "SngHadronPTS", //  Histogram container for SngHadronPTS_RUN_0000386773
        "SngHadronPT", // Histogram container for SngHadronPT_RUN_0000386773
        "SngHadronPTmN", // Histogram container for SngHadronPTmN_RUN_0000386773
        "SngHadronPTmS", // Histogram container for SngHadronPTmS_RUN_0000386773
        "SngHadronPTm", //  Histogram container for SngHadronPTm_RUN_0000386773
        "SngHadronPTpN", // Histogram container for SngHadronPTpN_RUN_0000386773
        "SngHadronPTpS", // Histogram container for SngHadronPTpS_RUN_0000386773
        "SngHadronPTp", //  Histogram container for SngHadronPTp_RUN_0000386773
        "SngMuonEtaW10", // Histogram container for SngMuonEtaW10_RUN_0000386773
        "SngMuonEtaW10m", //  Histogram container for SngMuonEtaW10m_RUN_0000386773
        "SngMuonEtaW10p", //  Histogram container for SngMuonEtaW10p_RUN_0000386773
        "SngMuonEtaW15", // Histogram container for SngMuonEtaW15_RUN_0000386773
        "SngMuonEtaW15m", //  Histogram container for SngMuonEtaW15m_RUN_0000386773
        "SngMuonEtaW15p", //  Histogram container for SngMuonEtaW15p_RUN_0000386773
        "SngMuonEtaW20", // Histogram container for SngMuonEtaW20_RUN_0000386773
        "SngMuonEtaW20m", //  Histogram container for SngMuonEtaW20m_RUN_0000386773
        "SngMuonEtaW20p", //  Histogram container for SngMuonEtaW20p_RUN_0000386773
        "SngMuonEta", //  Histogram container for SngMuonEta_RUN_0000386773
        "SngMuonEtam", // Histogram container for SngMuonEtam_RUN_0000386773
        "SngMuonEtap", // Histogram container for SngMuonEtap_RUN_0000386773
        "SngMuonPTN", //  Histogram container for SngMuonPTN_RUN_0000386773
        "SngMuonPTS", //  Histogram container for SngMuonPTS_RUN_0000386773
        "SngMuonPT", // Histogram container for SngMuonPT_RUN_0000386773
        "SngMuonPTmS", // Histogram container for SngMuonPTmS_RUN_0000386773
        "SngMuonPTm", //  Histogram container for SngMuonPTm_RUN_0000386773
        "SngMuonPTpN", // Histogram container for SngMuonPTpN_RUN_0000386773
        "SngMuonPTp" // Histogram container for SngMuonPTp_RUN_0000386773*/
      };
  //RelLumi have to be the first one for asymmetry calculation
//  TString list =
//  "/phenix/u/jinhuang/links/spin_data/data/Run13_SpinAna1/SngMuons_hist_10.lst";
  TString list = "/gpfs/mnt/gpfs02/phenix/spin/spin2/chenxu/spinAnalyzer/diMuon_his.lst";
  MergeHist(list, "spinQA.txt", n_hists, hist_names);
  CrossingCheck(list + ".root", "RelLumi");
  ////////////////////////////////////////////////////////////

//  Test();
//  ReadTest();

//  GetBeamPos();
//  GetBeamPos_Sum();

//  GetBeamPos_FitOffSet("data/taxi_1264.filelist_hist_oldVertex.root");
//  GetBeamPos_FitOffSet_FVTX("data/taxi_1264.filelist_hist_oldVertex.root");
//    GetBeamPos_FitOffSet("data/taxi_1264.filelist_hist_newVertex.root");

//    DCA_Pz("data/taxi_1264.filelist_hist_oldVertex.root");
  //  DCA_Pz("data/taxi_1264.filelist_hist_newVertex.root");

//    DCA_SinglePhiSlice("data/taxi_1264.filelist_hist_oldVertex.root");
//    DCA_SinglePhiSlice("data/taxi_1264.filelist_hist_newVertex.root");

//  DCA_SinglePhiSlice_Fill("data/taxi_1264.filelist_hist_oldVertex.root");
//  DCA_SinglePhiSlice_Fill("data/taxi_1264.filelist_hist_newVertex.root");
}

void
Test(
    TString infile =
        "/phenix/u/jinhuang/links/spin_data/data/Run13_SpinAna1/386828_0.root_DiMuons_hist.root")
{

  TFile *_file0 = TFile::Open(infile);

  saHist * h = RelLumi_saHist;

  h->AutoLoad(3);

  h->Print();

  delete _file0;

  h->Print();

}

void
MergeHist(
    TString infile =
        "/gpfs/mnt/gpfs02/phenix/spin/spin2/chenxu/spinAnalyzer/dimuon_his.lst",
    TString spinQA_file = "spinQA.txt", const int n_hists = 0,
    const char * hist_names[] = NULL)
{
  //load good runs
  fstream f(spinQA_file, ios_base::in);
  assert(f.is_open());

  saHist::n_list good_runs;

  string line;
  while (!f.eof())
    {
      getline(f, line);
      stringstream ss(line);

      int run = 0;

      ss >> run;

      if (run > 0)
        good_runs.push_back(run);
    }

  cout << good_runs.size() << " Good Runs loaded from " << spinQA_file << endl;

  //load source files
  TFile *_file1 = TFile::Open(infile + ".root", "recreate");

  fstream f_lst(infile, ios_base::in);

  string line;

  saHist * hs[1000] =
    { NULL };
  for (int i = 0; i < n_hists; i++)
    hs[i] = NULL;

  while (!f_lst.eof())
    {
      getline(f_lst, line);

      cout << "Process " << line << endl;
      if (line.length() == 0)
        continue;

      TFile *_file0 = TFile::Open(line.c_str());
      assert(_file0);

      for (int i = 0; i < n_hists; i++)
        {

          _file0->cd();

          const TString hist_name = hist_names[i];

//          cout << "Loading " << hist_name << endl;

          saHist * h_src = (saHist *) _file0->GetObjectChecked(
              hist_name + "_saHist", "saHist");
          assert(h_src);
          h_src->Verbosity(0, true);
          h_src->AutoLoad(0);
          h_src->Verbosity(0, true);

          h_src->RemoveAsymmetryHisto();

          if (!(hs[i]))
            {
              _file1->cd();
              cout
                  << "====================  h_src->MakeTemplate();  ====================="
                  << endl;
              (hs[i]) = h_src->MakeTemplate();
              assert((hs[i]));

//              cout << "Template Made " << hist_name << " : "
//                  << (hs[i])->GetName() << " : " << hs[i]->GetName() << endl;
            }

          (hs[i])->AdoptSubHist(h_src, true, true, good_runs);

          delete h_src;

//          cout << "Processed " << hist_name << " : " << (hs[i])->GetName()
//              << " : " << hs[0] << ":" << hs[1] << endl;
        }

      delete _file0;
      _file1->cd();
    }

  for (int i = 0; i < n_hists; i++)
    {

      const TString hist_name = hist_names[i];

      cout << "Final touch on " << hist_name << " : " << (hs[i])->GetName()
          << endl;

      if (i > 0)
        {

          cout
              << "====================  (hs[i])->CalcAsymmetry(hs[0], saHist::BbcNoCutPHENIX); =====================";
          (hs[i])->CalcAsymmetry(hs[0], saHist::BbcVertexCut);

          (hs[i])->Verbosity(1);

//          (hs[i])->CalcAsymmetry_Chi2Fit();
        }

      cout << "====================  h_new -> Write()  ====================="
          << endl;

      (hs[i])->Write();
    }

//  for (int i = 0; i < n_hists; i++)
//    {
//      delete (hs[i]);
//    }
  delete _file1;
}

void
CrossingCheck(TString infile =
    "data/sngmuon_run12pp510_taxi.list_hist_wrongnuchmask.root",
    TString histname = "SngMuonPT")
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *_file0 = TFile::Open(infile);

  saHist * RelLumi_saHist = (saHist *) _file0->GetObjectChecked(
      histname + "_saHist", "saHist");
  assert(RelLumi_saHist);
  RelLumi_saHist->AutoLoad();

  vector<int> runlist = RelLumi_saHist->get_run_list();

  cout << "Have " << runlist.size() << " runs" << endl;

  TH2F * h_run = new TH2F("h_run", "Crossing Check;RHIC Crossing;Run ID", 120,
      .5, 119.5, runlist.size(), -.5, runlist.size() - .5);
  TH1F * h_run_list = new TH1F("h_run_list",
      "Run look up table;Run ID;Run Number", runlist.size(), -.5,
      runlist.size() - .5);

  fstream f("run_list.txt", ios_base::out);

  for (int runid = 0; runid < runlist.size(); runid++)
    {
      saHist * RelLumi_run = RelLumi_saHist->getsaHist_run(runlist[runid]);

      assert(RelLumi_run);

      TH2F * h_crossing_check = RelLumi_run->getHisto("CROSSING_CHECKS");
      assert(h_crossing_check);

      assert(h_crossing_check->GetNbinsX() == 120);

      h_run_list->SetBinContent(runid + 1, runlist[runid]);
      f << runid << "\t" << runlist[runid] << endl;

      for (int c = 1; c <= h_crossing_check->GetNbinsX(); c++)
        {

          h_run->SetBinContent(c, runid + 1,
              h_crossing_check->GetBinContent(c, 1));

        }

    }

  f.close();

  TCanvas *c1 = new TCanvas("CrossingCheck", "CrossingCheck", 1000, 900);

  c1->Divide(2, 1);
  int idx = 1;
  TPad * p;

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  h_run->Draw("colz");
  p->SetLogz();

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  h_run_list->Draw();

  SaveCanvas(c1,
      TString(_file0->GetName()) + TString("_") + TString(c1->GetName()) + "_"
          + histname, kFALSE);

}

void
AsymmetryCheck(TString infile =
    "data/sngmuon_run12pp510_taxi.list_hist_wrongnuchmask.root",
    TString histname = "SngMuonPT", TString hist_type = "A_LL", int bin = 4)
{

//  gStyle->SetOptStat(111);
  gStyle->SetOptFit(1111);

  TFile *_file0 = TFile::Open(infile);

  saHist * RelLumi_saHist = (saHist *) _file0->GetObjectChecked(
      histname + "_saHist", "saHist");
  assert(RelLumi_saHist);
  RelLumi_saHist->AutoLoad(2);

  vector<int> runlist = RelLumi_saHist->get_run_list();

  cout << "Have " << runlist.size() << " runs" << endl;

  TGraphErrors * ge = new TGraphErrors(runlist.size());

  TH1F * h_run_list = new TH1F("h_run_list",
      "Run look up table;Run ID;Run Number", runlist.size(), -.5,
      runlist.size() - .5);

  fstream f("run_list.txt", ios_base::out);

  TH1 * h_type = NULL;
  for (int runid = 0; runid < runlist.size(); runid++)
    {
      saHist * RelLumi_run = RelLumi_saHist->getsaHist_run(runlist[runid]);

      assert(RelLumi_run);

      h_type = RelLumi_run->getHisto(hist_type.Data());
      assert(h_type);

      h_run_list->SetBinContent(runid + 1, runlist[runid]);
      f << runid << "\t" << runlist[runid] << endl;

      (ge->GetX())[runid] = runlist[runid];
      (ge->GetY())[runid] = h_type->GetBinContent(bin);
      (ge->GetEY())[runid] = h_type->GetBinError(bin);

    }

  f.close();

  TCanvas *c1 = new TCanvas("AsymmetryCheck", "AsymmetryCheck", 1000, 900);

//  c1->Divide(2, 1);
  int idx = 1;
  TPad * p;

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  ge->SetMarkerColor(kRed);
  ge->SetLineColor(kRed - 7);

  ge->Draw("A*E");

  ge->GetHistogram()->SetTitle(
      Form("%s, %s: bin %d;Run Number;%s", h_type->GetTitle(), hist_type.Data(),
          bin, hist_type.Data()));
  double range = TMath::Max(abs(ge->GetHistogram()->GetYaxis()->GetXmin()),
      abs(ge->GetHistogram()->GetYaxis()->GetXmax()));
  range = TMath::Min(1., TMath::Abs(range));
  ge->GetHistogram()->GetYaxis()->SetRangeUser(-1, 1);

  ge->Fit("pol0", "M");
  ge->SetName(histname + "_" + hist_type + Form("_bin%d_AsymCheck", bin));

  TLine * zero = new TLine(ge->GetHistogram()->GetXaxis()->GetXmin(), 0,
      ge->GetHistogram()->GetXaxis()->GetXmax(), 0);
  zero->Draw();

  SaveCanvas(c1,
      TString(_file0->GetName()) + TString("_") + TString(c1->GetName()) + "_"
          + histname + "_" + hist_type + Form("_bin%d", bin), kFALSE);

}

void
MakeGoodRunList(
    TString infile =
        "/direct/phenix+spin/phnxsp01/mxliu/run13/ana_pp510/Run13pp500_SpinData_26Jul2013.root")
{

  TFile *_file0 = TFile::Open(infile);
  assert(_file0);

  TTree * T = (TTree *) _file0->GetObjectChecked("T", "TTree");
  assert(T);

  const Long64_t n = T->Draw("RunNumber:ScalerA[1]:ScalerA[0]:ScalerA[2]", "1",
      "goff");
  assert(n >0);

  fstream f("spinQA.txt", ios_base::out);

  for (Long64_t i = 0; i < n; i++)
    {
      const double * runs = T->GetV1();
      const double * gl1gaps = T->GetV2();
      const double * crossing0s = T->GetV3();
      const double * crossing2s = T->GetV4();

      int run = (int) (runs[i]);
      Long64_t gl1gap = (Long64_t) (gl1gaps[i]);
      Long64_t crossing0 = (Long64_t) (crossing0s[i]);
      Long64_t crossing2 = (Long64_t) (crossing2s[i]);

      cout << "Run = " << run << ", gap = " << gl1gap << " - " << crossing0
          << " - " << crossing2 << endl;

      if (gl1gap == 0 && crossing0 > 0 && crossing2 > 0)
        {
          f << run;

          for (int j = 0; j < 120; j++)
            f << "\t" << 1;

          f << endl;
        }
    }

  f.close();
}

void
ReadTest()
{

  //  TFile *_file0 = TFile::Open("data/366636_0.lst_hist.root");
  //  TFile *_file0 = TFile::Open("data/taxi_1264.filelist_hist.root");
  //  TFile *_file0 = TFile::Open("data/taxi_1264_100.filelist_hist.root");
  TFile *_file0 = TFile::Open("data/taxi_1264_2.filelist_hist.root");

  saHist * InvMass_saHist = (saHist *) _file0->GetObjectChecked(
      "InvMass_saHist", "saHist");
  assert(InvMass_saHist);
  saHist * RelLumi_saHist = (saHist *) _file0->GetObjectChecked(
      "RelLumi_saHist", "saHist");
  assert(RelLumi_saHist);

  InvMass_saHist->AutoLoad();

  RelLumi_saHist->AutoLoad();

  ///////////////////////

  TCanvas *c1 = new TCanvas("Test_SpinAna_ReadTest", "Test_SpinAna_ReadTest",
      1000, 900);
  c1->Divide(2, 2);
  int idx = 1;
  TPad * p;

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  p->SetLogy();

  InvMass_saHist->getHisto("YIELD")->Draw();
  InvMass_saHist->getHisto("YIELD")->SetYTitle("Yield");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  InvMass_saHist->getHisto("A_LL")->Draw();
  InvMass_saHist->getHisto("A_LL")->SetYTitle("Double spin asymmetry");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  InvMass_saHist->getHisto("A_L_Blue")->Draw();
  InvMass_saHist->getHisto("A_L_Blue")->SetYTitle(
      "Single spin asymmetry (Blue)");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  InvMass_saHist->getHisto("A_L_Yellow")->Draw();
  InvMass_saHist->getHisto("A_L_Yellow")->SetYTitle(
      "Single spin asymmetry (Yelllow)");

  c1->Paint();

  SaveCanvas(c1,
      TString(_file0->GetName()) + TString("_") + TString(c1->GetName()),
      kFALSE);

}

void
DCA_Pz(TString infile = "data/taxi_1264.filelist_hist.root")
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *_file0 = TFile::Open(infile);

//  return;

  DCA_Pz_Arm_SUM->GetZaxis()->SetRange(1, 1);
  TH2D * DCA_Pz_Arm_SUM_xy_arm0 = (TH2D *) DCA_Pz_Arm_SUM->Project3D("xy_arm0");
  DCA_Pz_Arm_SUM->GetZaxis()->SetRange(2, 2);
  TH2D * DCA_Pz_Arm_SUM_xy_arm1 = (TH2D *) DCA_Pz_Arm_SUM->Project3D("xy_arm1");

  const TH2D * DCA_Pz_Arm_SUM_xy_arms[2] =
    { DCA_Pz_Arm_SUM_xy_arm0, DCA_Pz_Arm_SUM_xy_arm1 };
  const int bins[][2] =
    {
      { 1, 2 },
      { 3, 4 },
      { 5, 8 },
      { 9, 40 } };

  TCanvas *c1 = new TCanvas("DCA_Pz", "DCA_Pz", 1800, 1000);
  c1->Divide(4, 2);
  int idx = 1;

  for (int arm = 0; arm < 2; arm++)
    {
      const TH2D * DCA_Pz_Arm_SUM_xy_arm = DCA_Pz_Arm_SUM_xy_arms[arm];
      assert(DCA_Pz_Arm_SUM_xy_arm);

      for (int bin = 0; bin < 4; bin++)
        {
          p = (TPad *) c1->cd(idx++);
          c1->Update();

          p->SetLogy();

          const int bin1 = bins[bin][0];
          const int bin2 = bins[bin][1];

          const TString name = Form("_pz%d", bin);

          const double low_pz =
              DCA_Pz_Arm_SUM_xy_arm->GetXaxis()->GetBinLowEdge(bin1);
          const double high_pz =
              DCA_Pz_Arm_SUM_xy_arm->GetXaxis()->GetBinUpEdge(bin2);

          TH1D * h1 = DCA_Pz_Arm_SUM_xy_arm->ProjectionY(
              DCA_Pz_Arm_SUM_xy_arm->GetName() + name, bin1, bin2);
          assert(h1);

          TF1 * f1 = new TF1(Form("gaus1_%d_%d", arm, bin), "gaus", -1, 1);
          TF1 * f2 = new TF1(Form("dgaus1_%d_%d", arm, bin),
              "gaus + [3]*exp(-0.5*((x-[1])/[4])**2) + [5]", -1, 1);

          f1->SetLineColor(kGreen);
          f2->SetLineColor(kRed);
          f2->SetNpx(1000);

          h1->Fit(f1, "MR");

          f2->SetParameters(f1->GetParameter(0) / 2, f1->GetParameter(1),
              f1->GetParameter(2), f1->GetParameter(0) / 2,
              f1->GetParameter(2) / 4, 0);

          f2->SetParName(1, "Mean");
          f2->SetParName(5, "Const");

          f2->SetParName(2, "#sigma (W)");
          f2->SetParName(4, "#sigma (N)");

          f2->SetParName(0, "Amp (W)");
          f2->SetParName(3, "Amp (N)");

          f2->SetParName(2, "#sigma (W)");
          f2->SetParName(4, "#sigma (N)");

          h1->Fit(f2, "MR");

          h1->SetTitle(
              Form(
                  "Arm %d DCA for %.1f<P_{z}<%.1f GeV/c;vertex (Kalman fit - VTX) #upoint #hat{p_{T}} + cor. (cm)",
                  arm, low_pz, high_pz));

          h1->Draw();
        }

    }

  SaveCanvas(c1, TString(infile) + TString("_") + TString(c1->GetName()),
      kFALSE);
}

void
DCA_SinglePhiSlice(TString infile = "data/taxi_1264.filelist_hist.root")
{

  const int bin1 = 20;
  const int bin2 = 23;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *_file0 = TFile::Open(infile);

  DCA_Phi_Arm_SUM->GetZaxis()->SetRange(1, 1);
  TH2D * DCA_Phi_Arm_SUM_xy_arm0 = (TH2D *) DCA_Phi_Arm_SUM->Project3D(
      "xy_arm0");
  DCA_Phi_Arm_SUM->GetZaxis()->SetRange(2, 2);
  TH2D * DCA_Phi_Arm_SUM_xy_arm1 = (TH2D *) DCA_Phi_Arm_SUM->Project3D(
      "xy_arm1");

  DCA_Phi_Arm_SUM_xy_arm0->FitSlicesY();
  DCA_Phi_Arm_SUM_xy_arm1->FitSlicesY();

  const TString name = Form("_phi%d_%d", bin1, bin2);

  const double low_pz = DCA_Phi_Arm_SUM_xy_arm0->GetXaxis()->GetBinLowEdge(
      bin1);
  const double high_pz = DCA_Phi_Arm_SUM_xy_arm0->GetXaxis()->GetBinUpEdge(
      bin2);

  TH1D * h1 = DCA_Phi_Arm_SUM_xy_arm0->ProjectionY(
      "DCA_Phi_Arm_SUM_xy_arm0" + name, bin1, bin2);

  const int arm = 0;
  const int bin = 0;

  TF1 * f1 = new TF1(Form("gaus1_%d_%d", arm, bin), "gaus", -1, 1);
  TF1 * f2 = new TF1(Form("dgaus1_%d_%d", arm, bin),
      "gaus + [3]*exp(-0.5*((x-[1])/[4])**2) + [5]", -1, 1);

  f1->SetLineColor(kGreen);
  f2->SetLineColor(kRed);
  f2->SetNpx(1000);

  h1->Fit(f1, "MR0");

  f2->SetParameters(f1->GetParameter(0) / 2, f1->GetParameter(1),
      f1->GetParameter(2), f1->GetParameter(0) / 2, f1->GetParameter(2) / 4, 0);

  f2->SetParName(1, "Mean");
  f2->SetParName(5, "Const");

  f2->SetParName(2, "#sigma (W)");
  f2->SetParName(4, "#sigma (N)");

  f2->SetParName(0, "Amp (W)");
  f2->SetParName(3, "Amp (N)");

  f2->SetParName(2, "#sigma (W)");
  f2->SetParName(4, "#sigma (N)");

  h1->Fit(f2, "MR0");

  h1->SetTitle(
      Form(
          "Arm %d DCA for P_{z}>10 GeV/c;vertex (Kalman fit - VTX) #upoint #hat{p_{T}} + cor. (cm)",
          arm));

  TCanvas *c1 = new TCanvas("DCA_SinglePhiSlice", "DCA_SinglePhiSlice", 1800,
      900 * .8);
  c1->Divide(2, 1);
  int idx = 1;

  c1->cd(idx++);
  c1->Update();

  DCA_Phi_Arm_SUM_xy_arm0->GetXaxis()->SetRangeUser(low_pz, high_pz);
  DCA_Phi_Arm_SUM_xy_arm0->Draw("COLZ");
  DCA_Phi_Arm_SUM_xy_arm0_1->Draw("same");

  TF1 * f = new TF1("f1", "[0] * cos(x) + [1] * sin(x)", -10, 10);
  f->SetLineColor(kMagenta);
  DCA_Phi_Arm_SUM_xy_arm0_1->Fit(f, "M0");
  f->Draw("same");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  p->SetLogy();

  h1->Draw();
//  f1->Draw("same");
  f2->Draw("same");

  SaveCanvas(c1, TString(infile) + TString("_") + TString(c1->GetName()),
      kFALSE);
}

void
DCA_SinglePhiSlice_Fill(TString infile = "data/taxi_1264.filelist_hist.root")
{

  const int bin1 = 20;
  const int bin2 = 23;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *_file0 = TFile::Open(infile);
  assert(_file0);

  TH3F * DCA_Phi_Arm_SUM = (TH3F *) _file0->GetObjectChecked("DCA_Phi_Arm_SUM",
      "TH3F");
  saHist * sahist = (saHist *) _file0->GetObjectChecked("DCA_Phi_Arm_saHist",
      "saHist");
  assert(sahist);
//  sahist->Verbosity(0);
  sahist->AutoLoad();

  TH2F * DCA_Phi_Arm_FillSum = new TH2F("DCA_Phi_Arm_FillSum",
      "DCA VS Fill;Fill ID;DCA", //
      50, .5, 50.5, DCA_Phi_Arm_SUM->GetXaxis()->GetNbins(),
      DCA_Phi_Arm_SUM->GetXaxis()->GetXmin(),
      DCA_Phi_Arm_SUM->GetXaxis()->GetXmax());

  int nfill = 0;
  int fill = 0;
  TH3F * h3 = 0;
  for (int i = 0; i < sahist->nsaHists(); i++)
    {
      saHist * sub_sahist = sahist->getsaHist(i);

      string name = sub_sahist->get_name();

      cout << "Process " << name << endl;

      int n = sscanf(name.c_str(), "DCA_Phi_Arm_FILL_%d", &fill);

      if (n == 1)
        {
          h3 = (TH3F *) sub_sahist->getHisto("SUM");

          assert(h3);

          nfill++;
          cout << "Fill = " << fill << " -> " << nfill << endl;

          h3->GetZaxis()->SetRange(1, 1);
          TH2D * DCA_Phi_Arm_SUM_xy_arm0 = (TH2D *) h3->Project3D("xy_arm0");

          const TString name2 = Form("_phi%d_%d", bin1, bin2);

          const double low_pz =
              DCA_Phi_Arm_SUM_xy_arm0->GetXaxis()->GetBinLowEdge(bin1);
          const double high_pz =
              DCA_Phi_Arm_SUM_xy_arm0->GetXaxis()->GetBinUpEdge(bin2);

          TH1D * h1 = DCA_Phi_Arm_SUM_xy_arm0->ProjectionY(
              DCA_Phi_Arm_SUM_xy_arm0->GetName() + name2, bin1, bin2);

          for (int bin = 1; bin < h1->GetNbinsX(); bin++)
            {
              DCA_Phi_Arm_FillSum->SetBinContent(nfill, bin,
                  h1->GetBinContent(bin));
            }

        }
    } //      for (int i = 0; i < sahist->nsaHists(); i++)

  TH1D * h1 = DCA_Phi_Arm_FillSum->ProjectionY("DCA_Phi_Arm_FillSum_Fill26", 21,
      28);
  assert(DCA_Phi_Arm_FillSum_Fill26);
  h1->Sumw2();

  const int arm = 0;
  const int bin = 0;

  TF1 * f1 = new TF1(Form("gaus1_%d_%d", arm, bin), "gaus", -1, 1);
  TF1 * f2 = new TF1(Form("dgaus1_%d_%d", arm, bin),
      "gaus + [3]*exp(-0.5*((x-[1])/[4])**2) + [5]", -1, 1);

  f1->SetLineColor(kGreen);
  f2->SetLineColor(kRed);
  f2->SetNpx(1000);

  h1->Fit(f1, "MR0");

  f2->SetParameters(f1->GetParameter(0) / 2, f1->GetParameter(1),
      f1->GetParameter(2), f1->GetParameter(0) / 2, f1->GetParameter(2) / 4, 0);

  f2->SetParName(1, "Mean");
  f2->SetParName(5, "Const");

  f2->SetParName(2, "#sigma (W)");
  f2->SetParName(4, "#sigma (N)");

  f2->SetParName(0, "Amp (W)");
  f2->SetParName(3, "Amp (N)");

  f2->SetParName(2, "#sigma (W)");
  f2->SetParName(4, "#sigma (N)");

  h1->Fit(f2, "MR0");

  h1->SetTitle(
      Form(
          "Arm %d DCA for P_{z}>10 GeV/c;vertex (Kalman fit - VTX) #upoint #hat{p_{T}} + cor. (cm)",
          arm));

  TCanvas *c1 = new TCanvas("DCA_SinglePhiSlice_Fill",
      "DCA_SinglePhiSlice_Fill", 1800, 900);
  c1->Divide(2, 1);
  int idx = 1;

  c1->cd(idx++);
  c1->Update();

  DCA_Phi_Arm_FillSum->Draw("colz");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  p->SetLogy();

  h1->Draw();
//  f1->Draw("same");
  f2->Draw("same");

  SaveCanvas(c1, TString(infile) + TString("_") + TString(c1->GetName()),
      kFALSE);
}

void
GetBeamPos(TString data = "data/taxi_1627.list")
{

//  int verbosity = 1;
  int verbosity = 0;

  TVirtualFitter::SetDefaultFitter("Minuit2");

  fstream f(data, ios_base::in);
  string line;
  assert(f.is_open());

  fstream flog(data + ".GetBeamPos", ios_base::out);
  assert(flog.is_open());

  while (!f.eof())
    {
      getline(f, line);

      cout << "Process " << line << endl;

      int run = 0;
      int fill = 0;

      TH2F * h_XY_run = 0;
      TH2F * h_XY_fill = 0;

      TFile *_file0 = TFile::Open(line.c_str());
      if (!_file0)
        continue;

      saHist * sahist = (saHist *) _file0->GetObjectChecked("VTX_XY_saHist",
          "saHist");
      assert(sahist);
      sahist->Verbosity(verbosity);
      sahist->AutoLoad();

      for (int i = 0; i < sahist->nsaHists(); i++)
        {
          saHist * sub_sahist = sahist->getsaHist(i);

          string name = sub_sahist->get_name();

          if (run == 0)
            {
              int n = sscanf(name.c_str(), "VTX_XY_RUN_%d", &run);

              if (n == 1)
                {
                  cout << "Run = " << run << endl;
                  h_XY_run = (TH2F *) sub_sahist->getHisto("SUM");

                  assert(h_XY_run);
                }
            }

          if (fill == 0)
            {
              int n = sscanf(name.c_str(), "VTX_XY_FILL_%d", &fill);

              if (n == 1)
                {
                  cout << "Fill = " << fill << endl;
                  h_XY_fill = (TH2F *) sub_sahist->getHisto("SUM");

                  assert(h_XY_fill);
                }
            }
        } //      for (int i = 0; i < sahist->nsaHists(); i++)

      if (!h_XY_run || !h_XY_fill)
        continue;

      TH1D * projX = h_XY_run->ProjectionX();
      TH1D * projY = h_XY_run->ProjectionY();

      pair<double, double> x = DrawAndFit(projX);
      pair<double, double> y = DrawAndFit(projY);

      flog << run << "\t" << x.first << "\t" << y.first << "\t" << x.second
          << "\t" << y.second << "\t" << fill << endl;

      if (verbosity)
        {

          TCanvas *c1 = new TCanvas("Test_SpinAna_GetBeamPos",
              "Test_SpinAna_GetBeamPos", 1800 * .5, 900 * .5);
          c1->Divide(2, 1);
          int idx = 1;

          c1->cd(idx++);
          c1->Update();
          projX->Draw();

          c1->cd(idx++);
          c1->Update();
          projY->Draw();

          SaveCanvas(c1,
              TString(_file0->GetName()) + TString("_")
                  + TString(c1->GetName()), kFALSE);

          delete c1;
        }

      delete sahist;
      delete _file0;
    }
//      TProfile h_X =

  flog.close();

  return;

//
}

pair<double, double>
DrawAndFit(TH1D * h)
{
  assert(h);

  const double default_signma = .02;
  const double peak = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());

  TF1 * f1 = new TF1("gaus", "gaus", peak - 2 * default_signma,
      peak + 2 * default_signma);
  TF1 * f2 = new TF1("dgaus", "gaus + [3]*exp(-0.5*((x-[1])/[4])**2) + [5]",
      peak - 4 * default_signma, peak + 4 * default_signma);
  TF1 * f3 =
      new TF1("tgaus",
          "gaus + [3]*exp(-0.5*((x-[1])/[4])**2) + [5]+ [6]*exp(-0.5*((x-[1])/[7])**2)",
          peak - 4 * default_signma, peak + 4 * default_signma);

  f1->SetLineColor(kGreen);
  f2->SetLineColor(kRed);
  f3->SetLineColor(kBlue);

  f1->SetParameters(12, peak, default_signma);
  h->Fit(f1, "MRQ");

  f2->SetParameters(f1->GetParameter(0) / 2, f1->GetParameter(1),
      f1->GetParameter(2), f1->GetParameter(0) / 2, f1->GetParameter(2) * 2, 0);
  h->Fit(f2, "MRQ");

//  f3->SetParameters(f2->GetParameter(0), f2->GetParameter(1),
//      f2->GetParameter(2), f2->GetParameter(3), f2->GetParameter(4),
//      f2->GetParameter(5), f2->GetParameter(0), f2->GetParameter(2) / 2);
//
//  h->Fit(f3, "MR");

  return pair<double, double>(f2->GetParameter(1), f2->GetParError(1));
}

void
GetBeamPos_Sum(TString infile = "data/taxi_1627.list")
{

  int n = 0;
  double id[10000] =
    { 0 };
  double run[10000] =
    { 0 };
  double fill_id[10000] =
    { 0 };

  double x[10000] =
    { 0 };
  double ex[10000] =
    { 0 };
  double y[10000] =
    { 0 };
  double ey[10000] =
    { 0 };

  int min_run = 1e7;
  int max_run = 0;

  int last_fill = 0;

  int n_fill = 0;
  int fill_first_run_id[10000] =
    { 0 };

  // read in
  fstream flog(infile + ".GetBeamPos", ios_base::in);
  string line;
  while (getline(flog, line))
    {

      stringstream sline(line);

//      cout << line <<endl;

      int fill;

      sline >> run[n] >> x[n] >> y[n] >> ex[n] >> ey[n] >> fill;

      min_run = min_run > run[n] ? run[n] : min_run;
      max_run = max_run > run[n] ? max_run : run[n];

//        cout << n << "  -  " << run[n] << "\t" << x[n] << "\t" << y[n] << "\t"
//            << ex[n] << "\t" << ey[n]<<" - "
//            <<(fabs(x[n]) > .1) << endl;

      id[n] = n;

      if (ex[n] < .5 && ey[n] < .5 && ex[n] > 0.00001 && ey[n] > .00001)
        n++;

    }
  flog.close();

  TGraphErrors * gex = new TGraphErrors(n, run, x, 0, ex);
  gex->SetName("gex");
  TGraphErrors * gey = new TGraphErrors(n, run, y, 0, ey);
  gey->SetName("gey");

  TGraph * grun = new TGraph(n, id, run);
  grun->SetName("grun");

  // interpolate

  const int min_run_src = min_run;

  n = 0;
  fstream flog(infile + ".ExportRunList", ios_base::in);
  while (getline(flog, line))
    {

      stringstream sline(line);
      int tmp_run = 0;

      sline >> tmp_run;

      if (tmp_run < 365e3)
        continue;

      run[n] = tmp_run;

      if (n == 0)
        {
          cout << run[n] << " : " << line << endl;
        }

      min_run = min_run > run[n] ? run[n] : min_run;
      max_run = max_run > run[n] ? max_run : run[n];

      if (run[n] < min_run_src)
        {

          TF1 f2("fit2", "pol0", min_run_src, min_run_src + 100);

          gex->Fit(&f2, "MRQ");

          x[n] = f2.Eval(run[n]);
          ex[n] = 0;

          gey->Fit(&f2, "MRQ");
          y[n] = f2.Eval(run[n]);
          ey[n] = 0;

        }
      else
        {

          TF1 f("fit", "pol2", run[n] - 4, run[n] + 4);

          TFitResultPtr ptr = gex->Fit(&f, "MRQS");
          TFitResult* fit_res = ptr.Get();

          if (fit_res)
            {

              x[n] = f.Eval(run[n]);
              ex[n] = 0;

              gey->Fit(&f, "MRQ");
              y[n] = f.Eval(run[n]);
              ey[n] = 0;

            }
          else
            {

              cout << "Expand range for run " << run[n] << endl;

              TF1 f2("fit2", "pol0", run[n] - 1000, run[n] + 1000);

              gex->Fit(&f2, "MRQ");

              x[n] = f2.Eval(run[n]);
              ex[n] = 0;

              gey->Fit(&f2, "MRQ");
              y[n] = f2.Eval(run[n]);
              ey[n] = 0;

            }
        }
      n++;

    }
  flog.close();

  TGraphErrors * gx = new TGraphErrors(n, run, x);
  TGraphErrors * gy = new TGraphErrors(n, run, y);
  gx->SetName("gx");
  gy->SetName("gy");

  gex->SetLineColor(kBlue);
  gey->SetLineColor(kBlue);
  gx->SetLineColor(kRed);
  gy->SetLineColor(kRed);
  gx->SetLineWidth(2);
  gy->SetLineWidth(2);

  cout << " n = " << n << ", run[0] = " << run[0] << endl;
//  gx -> Print();

// save

  cout << "Save to " << infile + ".Smooth.GetBeamPos" << endl;
  fstream flog(infile + ".Smooth.GetBeamPos", ios_base::out);
  for (int i = 0; i < n; i++)
    {

      flog << (gx->GetX())[i] << "\t" << (gx->GetY())[i] << "\t"
          << (gy->GetY())[i] << "\t" << 0 << "\t" << 0 << endl;

    }
  flog.close();

//  ///direct/phenix+subsys+fvtx/jinhuang/GoldenDimuon/AlignDST_GD_DST/result/ALL.AlignDST.root
////  Average VTX shift VS FVTX = (-0.222241 0.099704) cm
//  const double VTX_Shift_X = -0.222241;
//  const double VTX_Shift_Y = 0.099704;

  ///Elog 797, 510pp taxi
  const double VTX_Shift_X = -(7.68547e-02 + 7.73623e-02) / 2;
  const double VTX_Shift_Y = -(3.95420e-02 + 4.44859e-02) / 2;

  cout << "Save to " << infile + ".Smooth_Shifted.GetBeamPos" << endl;
  fstream flog(infile + ".Smooth_Shifted.GetBeamPos", ios_base::out);
  for (int i = 0; i < n; i++)
    {

      flog << (gx->GetX())[i] << "\t" << (gx->GetY())[i] - VTX_Shift_X << "\t"
          << (gy->GetY())[i] - VTX_Shift_Y << "\t" << 0 << "\t" << 0 << endl;

    }
  flog.close();

  TCanvas *c1 = new TCanvas("GetBeamPos_Sum", "GetBeamPos_Sum", 1800, 900 * .8);
  c1->Divide(3, 1);
  int idx = 1;

  c1->cd(idx++);
  c1->Update();
  gPad->DrawFrame(min_run - 1, -.3, max_run + 1, .3,
      "Beam X VS run;Run;VTX Beam X (cm)");
  gex->Draw("*");
  gx->Draw("l");

  c1->cd(idx++);
  c1->Update();
  gPad->DrawFrame(min_run - 1, -.3, max_run + 1, .3,
      "Beam Y VS run;Run;VTX Beam Y (cm)");
  gey->Draw("*");
  gy->Draw("l");

  c1->cd(idx++);
  c1->Update();
  grun->Draw("a*l");

  SaveCanvas(c1, TString(infile) + TString("_") + TString(c1->GetName()),
      kFALSE);
}

void
GetBeamPos_FitOffSet(TString infile = "data/taxi_1264.filelist_hist.root")
{

  TFile *_file0 = TFile::Open(infile);

  DCA_Phi_Arm_SUM->GetZaxis()->SetRange(1, 1);
  TH2D * DCA_Phi_Arm_SUM_xy_arm0 = (TH2D *) DCA_Phi_Arm_SUM->Project3D(
      "xy_arm0");
  DCA_Phi_Arm_SUM->GetZaxis()->SetRange(2, 2);
  TH2D * DCA_Phi_Arm_SUM_xy_arm1 = (TH2D *) DCA_Phi_Arm_SUM->Project3D(
      "xy_arm1");

  DCA_Phi_Arm_SUM_xy_arm0->FitSlicesY();
  DCA_Phi_Arm_SUM_xy_arm1->FitSlicesY();

  TCanvas *c1 = new TCanvas("GetBeamPos_FitOffSet", "GetBeamPos_FitOffSet",
      1800, 900 * .8);
  c1->Divide(2, 1);
  int idx = 1;

  c1->cd(idx++);
  c1->Update();

  DCA_Phi_Arm_SUM_xy_arm0->Draw("COLZ");
  DCA_Phi_Arm_SUM_xy_arm0_1->Draw("same");

  TF1 * f = new TF1("f1", "[0] * cos(x) + [1] * sin(x)", -10, 10);
  f->SetLineColor(kMagenta);
  DCA_Phi_Arm_SUM_xy_arm0_1->Fit(f, "M0");
  f->Draw("same");

  c1->cd(idx++);
  c1->Update();

  DCA_Phi_Arm_SUM_xy_arm1->Draw("COLZ");
  DCA_Phi_Arm_SUM_xy_arm1_1->Draw("same");
  TF1 * f = new TF1("f2", "[0] * cos(x) + [1] * sin(x)", -10, 10);
  f->SetLineColor(kMagenta);
  DCA_Phi_Arm_SUM_xy_arm1_1->Fit(f, "M0");
  f->Draw("same");

  SaveCanvas(c1, TString(infile) + TString("_") + TString(c1->GetName()),
      kFALSE);
}

void
GetBeamPos_FitOffSet_FVTX(TString infile = "data/taxi_1264.filelist_hist.root")
{

  TFile *_file0 = TFile::Open(infile);

  FVTX_DCA_Phi_Arm_SUM->GetZaxis()->SetRange(1, 1);
  TH2D * FVTX_DCA_Phi_Arm_SUM_xy_arm0 =
      (TH2D *) FVTX_DCA_Phi_Arm_SUM->Project3D("xy_arm0");
  FVTX_DCA_Phi_Arm_SUM->GetZaxis()->SetRange(2, 2);
  TH2D * FVTX_DCA_Phi_Arm_SUM_xy_arm1 =
      (TH2D *) FVTX_DCA_Phi_Arm_SUM->Project3D("xy_arm1");

  FVTX_DCA_Phi_Arm_SUM_xy_arm0->FitSlicesY();
  FVTX_DCA_Phi_Arm_SUM_xy_arm1->FitSlicesY();

  TCanvas *c1 = new TCanvas("GetBeamPos_FitOffSet_FVTX",
      "GetBeamPos_FitOffSet_FVTX", 1800, 900 * .8);
  c1->Divide(2, 1);
  int idx = 1;

  c1->cd(idx++);
  c1->Update();

  FVTX_DCA_Phi_Arm_SUM_xy_arm0->Draw("COLZ");
  FVTX_DCA_Phi_Arm_SUM_xy_arm0_1->Draw("same");

  TF1 * f = new TF1("f1", "[0] * cos(x) + [1] * sin(x)", -10, 10);
  f->SetLineColor(kMagenta);
  FVTX_DCA_Phi_Arm_SUM_xy_arm0_1->Fit(f, "M0");
  f->Draw("same");

  c1->cd(idx++);
  c1->Update();

  FVTX_DCA_Phi_Arm_SUM_xy_arm1->Draw("COLZ");
  FVTX_DCA_Phi_Arm_SUM_xy_arm1_1->Draw("same");
  TF1 * f = new TF1("f2", "[0] * cos(x) + [1] * sin(x)", -10, 10);
  f->SetLineColor(kMagenta);
  FVTX_DCA_Phi_Arm_SUM_xy_arm1_1->Fit(f, "M0");
  f->Draw("same");

  SaveCanvas(c1, TString(infile) + TString("_") + TString(c1->GetName()),
      kFALSE);
}

TProfile *
FitProfile(TH2D *h2)
{
  TProfile * h1 = h2->ProfileX();

  TF1 * f = new TF1("f", "gaus", -1, 1);

  for (int i = 1; i <= h2->GetNbinsX(); i++)
    {
      h2->FitSlicesY(f, i, i)

    }

  return h1;
}
