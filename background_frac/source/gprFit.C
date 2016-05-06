//#include "/phenix/spin/phnxsp01/yuhw/CVS_PHENIX/install/include/GausProc.h"

#include "TFitter.h"

#include "gprFit.h"

#define __ERROR__ cout<<"ERROR: "<<__FILE__<<" : "<<__LINE__<<endl
#define __DEBUG__ cout<<"DEBUG: "<<__FILE__<<" : "<<__LINE__<<endl
#define __WARNING__ cout<<"WARNING: "<<__FILE__<<" : "<<__LINE__<<endl

const double JpsiMass = 3.096916;
const double PsiMass  = 3.68609;
const double xmin=1.5;
const double xmax=6.0;
const double xbinw = 0.05;
const double NSIGMA = 2.;

//---------------------------------------------------------------------------------------------------
double CrystalBall(double* x, double* par){
	//http://en.wikipedia.org/wiki/Crystal_Ball_function
	double xcur = x[0];
	double alpha = par[0];
	double n = par[1];
	double mu = par[2];
	double sigma = par[3];
	double N = par[4];
	TF1* exp = new TF1("exp","exp(x)",1e-20,1e20);
	double A; double B;
	if (alpha < 0){
		A = pow((n/(-1*alpha)),n)*exp->Eval((-1)*alpha*alpha/2);
		B = n/(-1*alpha) + alpha;
	} else {
		A = pow((n/alpha),n)*exp->Eval((-1)*alpha*alpha/2);
		B = n/alpha - alpha;
	} 
	double f;
	if ((xcur-mu)/sigma > (-1)*alpha)
		f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/
				(2*sigma*sigma));
	else
		f = N*A*pow((B- (xcur-mu)/sigma),(-1*n));
	delete exp;

	return f;
} 

//---------------------------------------------------------------------------------------------------
double CB_plus_CB_same_Jpsi_tail(double *x, double *par)
{
	double Jpsi = CrystalBall(x, par);
	double par_psip[] = {
		par[0],
		par[1],
		par[2]*PsiMass/JpsiMass,
		par[3]*PsiMass/JpsiMass,
		par[5]
	};
	double Psip = CrystalBall(x, par_psip);
	return Jpsi + Psip;
}

double CB_plus_CB(double *x, double *par)
{
	double Jpsi = CrystalBall(x, par);
	double par_psip[] = {
		par[5],
		par[6],
		par[2]*PsiMass/JpsiMass,
		par[3]*PsiMass/JpsiMass,
		par[7]
	};
	double Psip = CrystalBall(x, par_psip);
	return Jpsi + Psip;
}

//---------------------------------------------------------------------------------------------------
double CB_plus_Gaussian(double *x, double *par)
{
	double Jpsi = CrystalBall(x, par);
	double Psip = par[5]*TMath::Gaus(x[0], PsiMass/JpsiMass*par[2], PsiMass/JpsiMass*par[3], true);
	return Jpsi + Psip;
}

//---------------------------------------------------------------------------------------------------
double CB_plus_Gaussian_plus_pol3(double *x, double *par)
{
	double Jpsi = CrystalBall(x, par);
	double Psip = par[5]*TMath::Gaus(x[0], PsiMass/JpsiMass*par[2], PsiMass/JpsiMass*par[3], true);
	double Bkg = par[6] + par[7]*x[0] + par[8]*pow(x[0],2) + par[9]*pow(x[0],3);
	return Jpsi + Psip + Bkg;
}

//---------------------------------------------------------------------------------------------------
double Gaussian_plus_Gaussian_plus_pol3(double *x, double *par)
{
	double Jpsi = par[0]*TMath::Gaus(x[0], par[1], par[2], true);
	double Psip = par[3]*TMath::Gaus(x[0], PsiMass/JpsiMass*par[1], PsiMass/JpsiMass*par[2], true);
	double Bkg = par[4] + par[5]*x[0] + par[6]*pow(x[0],2) + par[7]*pow(x[0],3);
	return Jpsi + Psip + Bkg;
}


double CB_plus_Gaussian_nopeak_correlated(double *x, double *par)
{
	double Jpsi = CrystalBall(x, par);
	double Psip = par[5]*TMath::Gaus(x[0],par[6],par[7], true);
	return Jpsi+Psip;
}

//---------------------------------------------------------------------------------------------------
TF1 * tf1_CB_plus_Gaussian()
{
	TF1 *f = new TF1("model",CB_plus_Gaussian, xmin, xmax, 6);
	f->SetParNames(
			"alpha",
			"n",
			"mu",
			"sigma",
			"N_Jpsi",
			"N_Psip"
			);

	f->SetParLimits(0,1,10);
	f->SetParLimits(1,0,100);
	//f->SetParLimits(2,3.,4.);
	//f->SetParLimits(3,0.1,0.3);
	//f->SetParLimits(4,0,10000);
	f->SetParLimits(5,0,10000);

	f->SetParameters(2,1,3.09,0.18,50000,1000);
	return f;
}

//---------------------------------------------------------------------------------------------------
GausProc *gprFit::gp(
		//double &N_bkg,
		//double &N_bkg_err
		)
{
	if(use_fvtx == 0) cout<<"use_fvtx = 0" <<endl;
	else if(use_fvtx == 1) cout<<"use_fvtx = 1" <<endl; 
	else if(use_fvtx == -1) cout<<"use_fvtx = -1" <<endl; 
	else cout<<"use_fvtx out of range!!!"<<endl;

	std::string mass = "mass";
	if(use_fvtx == 1) mass = "mass_fvtxmutr";

	unsigned nPredictions=TMath::Nint((xmax - xmin)/xbinw);

	double zone1_min = 1.5;
	double zone1_max = 2.2;
	double zone2_min = 4.3;
	double zone2_max = 6.0;

	//gSystem->Load("/phenix/spin/phnxsp01/yuhw/CVS_PHENIX/install/lib/libGausProc.so");

	int verbosity=0;

	const std::string name_core(Form("BinningMode_%d_Arm%d_Charge%d_bin%03.1f_%03.1f",binning_mode,arm,charge_req,bin_min,bin_max));

	char *outname=Form("%s_GPR.root",name_core.data());
	char *data_out_name = Form("%s_data.root",name_core.data());

	//-----------------------------------------------------------------------
	const char *inname = "/gpfs/mnt/gpfs02/phenix/spin/spin2/chenxu/pAu_AN_Run15/Run15_pAu.root";
//	const char *inname = "/gpfs/mnt/gpfs02/phenix/spin/spin2/chenxu/spinAnalyzer/Xcheck/Run15_pp_goodrun.root";
	TFile *infile = TFile::Open(inname,"READ");
	if(!infile)
	{
		__ERROR__;	
		return NULL;
	}

//	TTree *dataTree = (TTree *)infile->Get("dimuon");
	TTree *dataTree = (TTree *)infile->Get("T");

	TH1D *hpm = new TH1D("hpm","hpm",nPredictions,xmin,xmax);
	TH1D *hpp = new TH1D("hpp","hpp",nPredictions,xmin,xmax);
	TH1D *hmm = new TH1D("hmm","hmm",nPredictions,xmin,xmax);

	std::string cut_base = "";
	if(binning_mode == 1)
		cut_base = Form("%s>%f&&%s<%f&&pT>%f&&pT<%f",mass.data(),xmin,mass.data(),xmax,bin_min,bin_max);
	else if(binning_mode == 2)
		cut_base = Form("%s>%f&&%s<%f&&TMath::Abs(rapidity)>%f&&TMath::Abs(rapidity)<%f",mass.data(),xmin,mass.data(),xmax,bin_min,bin_max);
	else
		cut_base = "(0)";
	if(arm==0) cut_base += "&&Tr0_pz>0&&Tr0_DG0<25.&&Tr1_DG0<25.&&((lvl1_trigscaled & 0x00400000)||(lvl1_trigscaled &0x00100000))";
   else if(arm==1) cut_base += "&&Tr0_pz<0&&Tr0_DG0<30.&&Tr1_DG0<30.&&((lvl1_trigscaled & 0x00800000)||(lvl1_trigscaled &0x00200000))";

	std::string cut_fvtx = "(";
	cut_fvtx += "TMath::Abs(Tr0_dr_fvtx)<2.";
	cut_fvtx += "&&TMath::Abs(Tr0_dtheta_fvtx)<.5";
	cut_fvtx += "&&TMath::Abs(Tr0_dphi_fvtx)<.5";
	cut_fvtx += "&&TMath::Abs(Tr1_dr_fvtx)<2.";
	cut_fvtx += "&&TMath::Abs(Tr1_dtheta_fvtx)<.5";
	cut_fvtx += "&&TMath::Abs(Tr1_dphi_fvtx)<.5";
	cut_fvtx += ")";

	std::string cut_nofvtx = "(";	
	cut_nofvtx += "Evt_vtxchi2<5.";
	cut_nofvtx += "&&abs(rapidity)<2.2&&abs(rapidity)>1.2";
	cut_nofvtx += "&&Tr0_DDG0<10. && Tr1_DDG0<10.";
	cut_nofvtx += "&&Tr0_ntrhits>9 && Tr1_ntrhits>9";
	cut_nofvtx += "&&Tr0_nidhits>5 && Tr1_nidhits>5";
	cut_nofvtx += "&&Tr0_lastgap>2 && Tr1_lastgap>2";
	cut_nofvtx += "&&Tr0_dca_r<5. && Tr1_dca_r<5.";
	cut_nofvtx += "&&Tr0_pz*Tr1_pz>0.";
	cut_nofvtx += "&&pT<6.";
	cut_nofvtx += "&&p<100.";
	cut_nofvtx += ")";

	if(use_fvtx == 1)
	{
		cut_base += "&&";
		cut_base += cut_fvtx;
	}
	else if(use_fvtx == -1)
	{
		cut_base += "&&!";
		cut_base += cut_fvtx;
	}

	else if(use_fvtx == 0)
	{
		cut_base += "&&";
		cut_base += cut_nofvtx;
	}
	dataTree->Project("hpm",mass.data(),(cut_base+"&& charge== 0").data());
	dataTree->Project("hpp",mass.data(),(cut_base+"&& charge== 2").data());
	dataTree->Project("hmm",mass.data(),(cut_base+"&& charge==-2").data());

	double npp = hpp->Integral();
	double nmm = hmm->Integral();
	double g = 2*sqrt(npp*nmm)/(npp+nmm);

	TH1D *hdata = NULL;
	if(charge_req == 0){
		hdata = (TH1D *)hpm->Clone("hdata");
	}
	else if(charge_req == 1){
		hdata = (TH1D *)hpp->Clone("hdata");	
		hdata->Add(hmm);
	}
	else if(charge_req == 2){
		hdata = (TH1D *)hpm->Clone("hdata");
		hdata->Sumw2();
		hdata->Add(hpp,-1*g);
		hdata->Add(hmm,-1*g);
	}

	std::vector<double> x,y,sigma_y;

	for(int i=1;i<=hdata->GetNbinsX();i++)
	{
		double bin_mean = hdata->GetBinLowEdge(i)+0.5*hdata->GetBinWidth(i);
		if((bin_mean>zone1_min && bin_mean<zone1_max) ||
				(bin_mean>zone2_min && bin_mean<zone2_max))
		{
			double y_err = hdata->GetBinError(i);
			if(TMath::Abs(y_err)<1e-10){
				cout<<__FILE__<<" : "<<__LINE__<<" : ERROR!"<<endl;
				continue;
			}
			x.push_back(bin_mean);
			y.push_back(hdata->GetBinContent(i));
			sigma_y.push_back(y_err);
		}
	}

	//-----------------------------------------------------------------------
	double xmin_new = xmin;
	double xmax_new = xmax;
	if(do_shift_bin_gpr)
	{
		xmin_new -= 0.5*xbinw;
		xmax_new += 0.5*xbinw;
	}
	nPredictions=TMath::Nint((xmax_new - xmin_new)/xbinw);

	GausProc a(x, y, sigma_y, xmin_new, xmax_new, nPredictions, outname);
	cout<<"GausProc a"<<"("<<"x"<<","<<"y"<<","<<"sigma_y"<<","<< xmin_new<<","<< xmax_new<<","<< nPredictions<<","<< outname<<");"<<endl;
	//gSystem->Exec(Form("rm -f %s",outname));

	for(unsigned int i=0;i<x.size();i++)
		cout<<x[i]<<" "<<y[i]<<" "<<sigma_y[i]<<endl;

	a.SetVerbosity(verbosity);
	a.SetKernel(GausProc::RBF);
	//a.warp(0);//should be used only if you have data spanning several orders of magnitude

	GPOptimizer c(&a,2,10);
	c.GPoptimize(2,0);
	cout<<c.getPar(0)<<" "<<c.getPar(1)<<endl;  

	GausProc *d = new GausProc(x, y, sigma_y, xmin_new, xmax_new, nPredictions, outname);
	//d.warp(0);//see comment above
	d->SetPar(0,c.getPar(0));
	d->SetPar(1,c.getPar(1));
	d->process();
	d->Write(-1);
	//d->unwarp(0);//see comment above

	TFile *data_out_file = new TFile(data_out_name,"RECREATE");
	data_out_file->cd();
	hdata->Write();
	data_out_file->Close();

	infile->Close();

	return d;
	//d.Integral(lowerBound_int, upperBound_int+0.05, integral, dIntegral);

	cout<<"End of gp"<<endl;
}

//---------------------------------------------------------------------------------------------------
void gprFit::gpr_extract
(
)
{
	const std::string name_core(Form("BinningMode_%d_Arm%d_Charge%d_bin%03.1f_%03.1f",binning_mode,arm,charge_req,bin_min,bin_max));


	const std::string gpr_name(Form("%s_GPR.root",name_core.data()));
	const std::string data_name(Form("%s_data.root",name_core.data()));
	char *outname = Form("%s_GPR_EXTRACT.root",name_core.data());

	TFile *data_file = new TFile(data_name.data(),"READ");
	//TH1D *data_hist = (TH1D *)data_file->GetObjectChecked(name_core.data(),"TH1D");
	TH1D *data_hist = (TH1D *)data_file->GetObjectChecked("hdata","TH1D");

	TFile *gpr_file = new TFile(gpr_name.data(),"READ");
	//TH1F *gpr_hist = (TH1F *)gpr_file->GetObjectChecked("ho","TH1F"); 
	TH1F *gpr_hist = (TH1F *)gpr_file->Get("ho"); 

	TH1F *ho_shift_back = (TH1F*)data_hist->Clone("ho_shift_back"); 
	ho_shift_back->SetTitle("ho_shift_back");
	ho_shift_back->Reset();

	TH1D *extract_hist = (TH1D *)data_hist->Clone("hextract");

	if(do_shift_bin_gpr)
		if(data_hist->GetNbinsX()!=gpr_hist->GetNbinsX()-1)
		{
			cout<<"data_hist->GetNbinsX()!=gpr_hist->GetNbinsX()-1"<<endl;
			return;
		}

	for(int i=1; i<=data_hist->GetNbinsX();i++)
	{
		double gpr_mean = gpr_hist->GetBinContent(i);
		double gpr_error = gpr_hist->GetBinError(i);

		if(do_shift_bin_gpr)
		{
			gpr_mean = 0.5*(gpr_hist->GetBinContent(i) + gpr_hist->GetBinContent(i+1));
			gpr_error = 0.5*(gpr_hist->GetBinError(i) + gpr_hist->GetBinError(i+1));
		}

		ho_shift_back->SetBinContent(i,gpr_mean);
		ho_shift_back->SetBinError(i,gpr_error);

		double data_mean = data_hist->GetBinContent(i);
		double data_error = data_hist->GetBinError(i);

		double extract_mean = data_mean - gpr_mean;
		double extract_error = sqrt(data_error*data_error+gpr_error*gpr_error);

		extract_hist->SetBinContent(i,extract_mean);
		extract_hist->SetBinError(i,extract_error);
	}

	TCanvas *c0 = new TCanvas((name_core+"_extract").data(),(name_core+"_extract").data());

	double hist_max = 1.1*data_hist->GetMaximum();
	double hist_min = 1.5*extract_hist->GetMinimum();

	data_hist->SetTitle(name_core.data());
	data_hist->SetMaximum(hist_max);
	data_hist->SetMinimum(hist_min);

	data_hist->SetLineColor(kBlack);
	data_hist->SetMarkerStyle(4);
	data_hist->SetMarkerColor(kBlack);
	data_hist->Draw("e");

	if(!plot_ho_shift_back)
	{
	gpr_hist->SetLineColor(kBlue);
	gpr_hist->SetMarkerStyle(1);
	gpr_hist->SetMarkerColor(kBlue);
	gpr_hist->Draw("same");
	}
	else
	{
		ho_shift_back->SetLineColor(kBlue);
		ho_shift_back->SetMarkerStyle(1);
		ho_shift_back->SetMarkerColor(kBlue);
		ho_shift_back->Draw("same");
	}

	extract_hist->SetLineColor(kRed);
	extract_hist->SetMarkerStyle(4);
	extract_hist->SetMarkerColor(kRed);
	extract_hist->Draw("same");

	TBox *b0 = new TBox(1.5,hist_min,2.1,hist_max);
	b0->SetLineColor(0);
	b0->SetFillColor(kGreen);
	b0->SetFillStyle(3013);
	b0->Draw();

	TBox *b1 = new TBox(4.3,hist_min,5.0,hist_max);
	b1->SetLineColor(0);
	b1->SetFillColor(kGreen);
	b1->SetFillStyle(3013);
	b1->Draw();

	c0->Update();

	TFile *out_file = new TFile(outname,"RECREATE");
	out_file->cd();
	c0->Write();
	data_hist->Write();
	gpr_hist->Write();
	ho_shift_back->Write();
	extract_hist->Write();
	out_file->Close();

	data_file->Close();
	gpr_file->Close();

	cout<<"End of gpr_extract"<<endl;
}

//---------------------------------------------------------------------------------------------------
bool gprFit::roofit
(
 GausProc *gaus_proc
 //double &N_bkg,
 //double &N_bkg_err
 )
{
	using namespace RooFit;

	//double xmin = 1.5;
	//double xmax = 6.0;
	//double xbinw = 0.05;//0.025;
	Int_t xnbin= TMath::Nint((xmax - xmin)/xbinw);

	double zone1_min = 1.5;
	double zone1_max = 2.2;
	double zone2_min = 4.3;
	double zone2_max = 6.0;

	const std::string name_core(Form("BinningMode_%d_Arm%d_Charge%d_bin%03.1f_%03.1f",binning_mode,arm,charge_req,bin_min,bin_max));


	gROOT->SetStyle("Plain");
	RooMsgService::instance().setSilentMode(true);

	char *f0 = Form("%s_GPR_EXTRACT.root",name_core.data());
	TFile *dataFile = new TFile(f0,"READ");

	TFile *fout = TFile::Open(Form("%s_roofit.root",name_core.data()),"recreate");

	TH1D *hdata = (TH1D *) dataFile->Get("hdata");
	TH1D *hextract = (TH1D *) dataFile->Get("hextract");
	TH1F *ho = (TH1F *) dataFile->Get("ho");
	if(plot_ho_shift_back)
		ho = (TH1F *) dataFile->Get("ho_shift_back");

	RooRealVar x("x", "x", xmin, xmax);

	RooWorkspace *w0 = new RooWorkspace("w0",kTRUE);
	RooDataHist *roohextract = new RooDataHist("roohextract","roohextract",x,hextract);
	RooDataHist *rooho = new RooDataHist("rooho","rooho",x,ho);
	RooDataHist *roohdata = new RooDataHist("roohdata","roohdata",x,hdata);
	TH1D *training_zone_1 = hmanip::RangeClone(hdata,zone1_min,zone1_max);
	TH1D *training_zone_2 = hmanip::RangeClone(hdata,zone2_min,zone2_max);

	TCanvas* c0 = new TCanvas((name_core+"_roofit").data(),(name_core+"_roofit").data(),100,100,1000,800);
	c0->SetLogy();

	w0->import(*roohextract);

	switch(model_choice)
	{
		case 4://Jpsi:CB, Psi(2s): Gaussian
			w0->factory("Gaussian::Psip(x,expr('1.190*Jpsi_mean',Jpsi_mean[3.09,3.0,3.2]),expr('1.190*Jpsi_sigma',Jpsi_sigma[0.18,0.1,0.3]))");
			//w0->factory("CBShape::Jpsi(x,Jpsi_mean,Jpsi_sigma,Jpsi_alpha[2.,1.,4.],Jpsi_n[150,0,1000.])");
			w0->factory("CBShape::Jpsi(x,Jpsi_mean,Jpsi_sigma,Jpsi_alpha[2.,1.,10.],Jpsi_n[1,0.0,10.])");
			w0->factory("SUM::model(N_Jpsi[60000,5000,200000]*Jpsi,N_Psip[1000,0,10000]*Psip)");
			break;
		case 3://Jpsi, Psi(2s): CB
			//w0->factory("CBShape::Psip(x,expr('1.190*Jpsi_mean',Jpsi_mean[3.09,3.0,3.2]),expr('1.190*Jpsi_sigma',Jpsi_sigma[0.18,0.1,0.3]),Psip_alpha[2.,1.,3],Psip_n[1,0,3])");
			w0->factory("CBShape::Psip(x,Psip_mean[3.686,3.3,3.9],Psip_sigma[0.20,0.1,0.3],Psip_alpha[2.,1.,4.],Psip_n[5,0,10.])");
			w0->factory("CBShape::Jpsi(x,Jpsi_mean[3.097,3.0,3.2],Jpsi_sigma[0.18,0.1,0.3],Jpsi_alpha[2.,1.,4.],Jpsi_n[50,0,10000.])");
			w0->factory("SUM::model(N_Jpsi[30000,1000,200000]*Jpsi,N_Psip[100,0,10000]*Psip)");
			break;
		case 2://Jpsi, Psi(2s): Gaussian
			//w0->factory("Gaussian::Psip(x,Psip_mean[3.686,3.3,3.9],Psip_sigma[0.20,0.1,0.3])");
			//w0->factory("Gaussian::Jpsi(x,Jpsi_mean[3.097,3.0,3.2],Jpsi_sigma[0.18,0.1,0.3])");
			//w0->factory("SUM::model(N_Jpsi[30000,1000,200000]*Jpsi,N_Psip[100,0,10000]*Psip)");
			w0->factory("Gaussian::Psip(x,expr('1.190*Jpsi_mean',Jpsi_mean[3.09,3.0,3.2]),expr('1.190*Jpsi_sigma',Jpsi_sigma[0.18,0.1,0.3]))");
			w0->factory("Gaussian::Jpsi(x,Jpsi_mean,Jpsi_sigma)");
			w0->factory("SUM::model(N_Jpsi[60000,5000,200000]*Jpsi,N_Psip[1000,0,10000]*Psip)");
			break;
		case 1://Jpsi: Double Gaussian, Psi(2s): Gaussian
			//w0->factory("Gaussian::Psip(x,Psip_mean[3.686,3.3,3.9],Psip_sigma[0.20,0.1,0.3])");
			//w0->factory("Gaussian::Jpsi_n(x,Jpsi_n_mean[3.097,3.0,3.2],Jpsi_n_sigma[0.14,0.1,0.2])");
			//w0->factory("Gaussian::Jpsi_w(x,Jpsi_w_mean[3.097,3.0,3.2],Jpsi_w_sigma[0.27,0.2,0.4])");
			w0->factory("Gaussian::Psip(x,expr('Jpsi_mean+0.589174',Jpsi_mean[3.09,3.0,3.2]),Psip_sigma[0.20,0.1,0.3])");
			w0->factory("Gaussian::Jpsi_n(x,Jpsi_mean,Jpsi_n_sigma[0.14,0.1,0.3])");
			w0->factory("Gaussian::Jpsi_w(x,Jpsi_mean,Jpsi_w_sigma[0.27,0.2,0.6])");
			w0->factory("SUM::Jpsi(frac_n[0.5,0,1]*Jpsi_n,Jpsi_w)");
			w0->factory("SUM::model(N_Jpsi[60000,1000,400000]*Jpsi,N_Psip[100,0,10000]*Psip)");
			break;

		default:
			cout<<"Error, bad input, quitting\n";
			break;
	}

	//w0->Print("t");
	RooFitResult *fr = NULL;
	fr = (w0->pdf("model")->fitTo(*roohextract,Save(),Verbose(false),PrintLevel(-1),Warnings(false),NumCPU(8)));

	double Jpsi_mean = 0;
	double Jpsi_sigma = 0;
	if(model_choice==1)
	{
		double frac_n = w0->var("frac_n")->getVal();
		double Jpsi_n_sigma = w0->var("Jpsi_n_sigma")->getVal();
		double Jpsi_w_sigma = w0->var("Jpsi_w_sigma")->getVal();

		//double Jpsi_n_mean = w0->var("Jpsi_n_mean")->getVal();
		//double Jpsi_w_mean = w0->var("Jpsi_w_mean")->getVal();
		//Jpsi_mean = frac_n*Jpsi_n_mean + (1-frac_n)*Jpsi_w_mean;

		double Jpsi_mean = w0->var("Jpsi_mean")->getVal();

		Jpsi_sigma = frac_n*Jpsi_n_sigma + (1-frac_n)*Jpsi_w_sigma;
	}
	else
	{
		Jpsi_mean = w0->var("Jpsi_mean")->getVal();
		Jpsi_sigma = w0->var("Jpsi_sigma")->getVal();
	}

	double xmin_win_2sigma = Jpsi_mean - 2. * Jpsi_sigma;
	double xmax_win_2sigma = Jpsi_mean + 2. * Jpsi_sigma;
	double xmin_win_3sigma = Jpsi_mean - 3. * Jpsi_sigma;
	double xmax_win_3sigma = Jpsi_mean + 3. * Jpsi_sigma;

	RooPlot *frame = x.frame(Title(name_core.data()));
	//RooPlot *frame = x.frame("");

	if(is_for_paper == true)
		roohextract->plotOn(frame,LineColor(kRed),MarkerStyle(6),MarkerColor(kRed),XErrorSize(0));
	else
		roohextract->plotOn(frame,LineColor(kRed),MarkerStyle(6),MarkerColor(kRed));

	//hextract->SetLineColor(kRed);
	//hextract->SetMarkerStyle(4);
	//hextract->SetMarkerColor(kRed);
	//hextract->Draw("9e1x0p");

	//w0->pdf("model")->plotOn(frame,LineColor(kRed));
	w0->pdf("model")->plotOn(frame,Normalization(w0->var("N_Jpsi")->getVal()+w0->var("N_Psip")->getVal(),RooAbsReal::NumEvent),LineColor(kRed));
	//w0->pdf("model")->paramOn(frame);
	Int_t nFitParam = 0;
	if(model_choice==1) nFitParam = 9;
	if(model_choice==4) nFitParam = 6;
	double Chi2FitX = frame->chiSquare(nFitParam);
	w0->pdf("Jpsi")->plotOn(frame,Normalization(w0->var("N_Jpsi")->getVal(),RooAbsReal::NumEvent),LineStyle(kDashed),LineColor(kGreen));
	w0->pdf("Psip")->plotOn(frame,Normalization(w0->var("N_Psip")->getVal(),RooAbsReal::NumEvent),LineStyle(kDashed),LineColor(kBlue));
	if(is_for_paper)
		frame->SetTitle("");
	frame->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} / GeV");
	frame->GetYaxis()->SetTitle(Form("Events / %.0fMeV",xbinw*1000));

	roohdata->plotOn(frame,LineColor(kBlack),MarkerStyle(6),MarkerColor(kBlack),XErrorSize(0));
	rooho->plotOn(frame,LineColor(kBlue),MarkerStyle(6),MarkerColor(kBlue),XErrorSize(0));

	//frame->SetMinimum(1);
	//frame->SetMaximum(100000);
	//frame->GetYaxis()->SetRangeUser(10,10000);
	frame->GetYaxis()->SetRangeUser(10,
			2.5*hextract->GetMaximum());
	frame->Draw();

	double N_Jpsi = w0->var("N_Jpsi")->getVal();
	double N_Psip = w0->var("N_Psip")->getVal();

	x.setRange("window_2sigma",xmin_win_2sigma,xmax_win_2sigma);
	double N_Jpsi_2sigma = (w0->pdf("Jpsi"))->createIntegral(x,NormSet(x),Range("window_2sigma"))->getVal() * w0->var("N_Jpsi")->getVal();
	double eN_Jpsi_2sigma = (w0->pdf("Jpsi"))->createIntegral(x,NormSet(x),Range("window_2sigma"))->getVal() * w0->var("N_Jpsi")->getError();
	double N_Psip_2sigma = (w0->pdf("Psip"))->createIntegral(x,NormSet(x),Range("window_2sigma"))->getVal() * w0->var("N_Psip")->getVal();

	x.setRange("window_3sigma",xmin_win_3sigma,xmax_win_3sigma);
	double N_Jpsi_3sigma = (w0->pdf("Jpsi"))->createIntegral(x,NormSet(x),Range("window_3sigma"))->getVal() * w0->var("N_Jpsi")->getVal();
	double eN_Jpsi_3sigma = (w0->pdf("Jpsi"))->createIntegral(x,NormSet(x),Range("window_3sigma"))->getVal() * w0->var("N_Jpsi")->getError();
	double N_Psip_3sigma = (w0->pdf("Psip"))->createIntegral(x,NormSet(x),Range("window_3sigma"))->getVal() * w0->var("N_Psip")->getVal();

	int lbin_2sigma = TMath::Nint((xmin_win_2sigma - ho->GetXaxis()->GetXmin())/xbinw);
	int hbin_2sigma = TMath::Nint((xmax_win_2sigma - ho->GetXaxis()->GetXmin())/xbinw);
	double N_BKG_2sigma = ho->Integral(lbin_2sigma,hbin_2sigma); 
	double N_ALL_2sigma = hdata->Integral(lbin_2sigma,hbin_2sigma); 

	int lbin_3sigma = TMath::Nint((xmin_win_3sigma - ho->GetXaxis()->GetXmin())/xbinw);
	int hbin_3sigma = TMath::Nint((xmax_win_3sigma - ho->GetXaxis()->GetXmin())/xbinw);
	double N_BKG_3sigma = ho->Integral(lbin_3sigma,hbin_3sigma); 
	double N_ALL_3sigma = hdata->Integral(lbin_3sigma,hbin_3sigma); 

	double r_2sigma = 1.0 - 1.0*N_Jpsi_2sigma/N_ALL_2sigma;
	double r_3sigma = 1.0 - 1.0*N_Jpsi_3sigma/N_ALL_3sigma;
	double er_2sigma = hmanip::get_ratio_error(N_Jpsi_2sigma,N_ALL_2sigma,eN_Jpsi_2sigma,TMath::Sqrt(N_ALL_2sigma));
	double er_3sigma = hmanip::get_ratio_error(N_Jpsi_3sigma,N_ALL_3sigma,eN_Jpsi_3sigma,TMath::Sqrt(N_ALL_3sigma));

	//use TF1::IntegralError to calculate error
	int error_calculation_method = 0;
	if(error_calculation_method == 1)
	{
		RooAbsPdf *pdf_Jpsi = w0->pdf("Jpsi");

		//RooArgSet *prodSet_Jpsi = new RooArgSet(*pdf_Jpsi);
		//prodSet_Jpsi->add(*w0->var("N_Jpsi"));
		//RooProduct *ext_pdf_Jpsi = new RooProduct("ext_pdf_Jpsi", "extended pdf Jpsi", *prodSet_Jpsi);
		//TF1 *tf_Jpsi = ext_pdf_Jpsi->asTF(RooArgList(x), *(ext_pdf_Jpsi->getParameters(RooArgSet(x))));

		TF1 *tf_Jpsi = pdf_Jpsi->asTF(RooArgList(x), *(pdf_Jpsi->getParameters(RooArgSet(x))));

		RooAbsPdf *pdf_model = w0->pdf("model");
		TF1 *tf_model = pdf_model->asTF(RooArgList(x), *(pdf_model->getParameters(RooArgSet(x))));
		cout
			<<"DEBUG: "
			<<endl<<"----------"<<endl;

		RooRealVar test_mean("test_mean","test_mean",0,-1,1);
		RooRealVar test_sigma("test_sigma","test_sigma",1,-10,10);
		RooGaussian *test_pdf = new RooGaussian("test_pdf", "test_pdf", x, test_mean, test_sigma);
		TF1 *tf_test_pdf = test_pdf->asTF(RooArgList(x), *(test_pdf->getParameters(RooArgSet(x))));

		tf_test_pdf->Print();

		//TCanvas *c_test = new TCanvas("c_test","c_test",100,100,1400,800);
		//c->Divide("")

		cout
			<<" 1: "<<tf_test_pdf->Eval(1.)
			<<" 2: "<<tf_test_pdf->Eval(0)
			<<" 3: "<<tf_test_pdf->Eval(-1)
			<<endl;

		cout
			<<" 1: "<<tf_model->Eval(1.)
			<<" 2: "<<tf_model->Eval(2)
			<<" 3: "<<tf_model->Eval(3)
			<<endl;

		tf_Jpsi->Print();
		cout<<"ExpFormula: "<<tf_Jpsi->GetExpFormula().Data()<<endl;;
		TCanvas *c_temp = new TCanvas("c_temp","c_temp");
		tf_Jpsi->GetHistogram()->Draw();

		cout
			<<"DEBUG: "
			<<" 1: "<<tf_Jpsi->Eval(1.)
			<<" 2: "<<tf_Jpsi->Eval(2.)
			<<" 3: "<<tf_Jpsi->Eval(3.)
			<<endl;

		fr->Print();

		pdf_Jpsi->Print();


		cout
			<<endl<<"----------"<<endl
			<<endl;

		//RooArgList(x).Print();
		//(pdf_Jpsi->getParameters(RooArgSet(x)))->Print();

		TMatrixDSym cov_matrix_jpsi_psip = fr->covarianceMatrix();
		//cov_matrix_jpsi_psip.Print();

		if(model_choice == 4)
		{
			int nparjpsi = 5;
			TMatrixDSym cov_matrix_jpsi(nparjpsi);
			for(int i=0;i<nparjpsi;i++)
			{
				for(int j=0;j<nparjpsi;j++)
				{
					cov_matrix_jpsi[i][j] = cov_matrix_jpsi_psip[i][j];
				}
			}

			cov_matrix_jpsi.Print();

			N_Jpsi_2sigma = tf_Jpsi->Integral(xmin_win_2sigma,xmax_win_2sigma)/xbinw;
			eN_Jpsi_2sigma = tf_Jpsi->IntegralError(xmin_win_2sigma,xmax_win_2sigma,tf_Jpsi->GetParameters(),cov_matrix_jpsi.GetMatrixArray())/xbinw;

			cout
				<<"DEBUG: "
				<<"IntegralError: "<< N_Jpsi_2sigma <<" +- "<<eN_Jpsi_2sigma
				<<endl;
		}
	}

	//double r_2sigma = 1.0*(N_BKG_2sigma+N_Psip_2sigma)/(N_Jpsi_2sigma+N_BKG_2sigma+N_Psip_2sigma);
	//double r_3sigma = 1.0*(N_BKG_3sigma+N_Psip_3sigma)/(N_Jpsi_3sigma+N_BKG_3sigma+N_Psip_3sigma);
	//double er_2sigma = get_binomial_error(N_BKG_2sigma+N_Psip_2sigma,N_Jpsi_2sigma+N_BKG_2sigma+N_Psip_2sigma);
	//double er_3sigma = get_binomial_error(N_BKG_3sigma+N_Psip_3sigma,N_Jpsi_3sigma+N_BKG_3sigma+N_Psip_3sigma);

	bool draw_nsigma_line = false;
	if(draw_nsigma_line)
	{
		TLine *lLowSig_2sigma = new TLine(xmin_win_2sigma,0,xmin_win_2sigma,frame->GetMaximum());
		lLowSig_2sigma->SetLineColor(kRed);
		lLowSig_2sigma->Draw();
		TLine *lHigSig_2sigma = new TLine(xmax_win_2sigma,0,xmax_win_2sigma,frame->GetMaximum());
		lHigSig_2sigma->SetLineColor(kRed);
		lHigSig_2sigma->Draw();

		TLine *lLowSig_3sigma = new TLine(xmin_win_3sigma,0,xmin_win_3sigma,frame->GetMaximum());
		lLowSig_3sigma->SetLineColor(kRed);
		lLowSig_3sigma->SetLineStyle(2);
		lLowSig_3sigma->Draw();
		TLine *lHigSig_3sigma = new TLine(xmax_win_3sigma,0,xmax_win_3sigma,frame->GetMaximum());
		lHigSig_3sigma->SetLineColor(kRed);
		lHigSig_3sigma->SetLineStyle(2);
		lHigSig_3sigma->Draw();
	}

	double hist_min = frame->GetMinimum();
	double hist_max = frame->GetMaximum();

	training_zone_1->SetFillColor(kGreen);
	training_zone_1->SetFillStyle(3013);
	training_zone_1->Draw("same");

	training_zone_2->SetFillColor(kGreen);
	training_zone_2->SetFillStyle(3013);
	training_zone_2->Draw("same");


	//TBox *b0 = new TBox(zone1_min,hist_min,zone1_max,hist_max);
	//b0->SetLineColor(0);
	//b0->SetFillColor(kGreen);
	//b0->SetFillStyle(3013);
	//b0->Draw();

	//TBox *b1 = new TBox(zone2_min,hist_min,zone2_max,hist_max);
	//b1->SetLineColor(0);
	//b1->SetFillColor(kGreen);
	//b1->SetFillStyle(3013);
	//b1->Draw();

	//double x0 = 0.65; 
	//double y0 = 0.8;
	//double decrease = 0.065;
	//TLatex *t = new TLatex();
	//t->SetNDC();
	//t->SetTextFont(12);
	////t->SetTextFont(62);
	//t->SetTextColor(36);
	//t->SetTextSize(0.04);
	//t->SetTextAlign(12);
	//t->DrawLatex(x0,y0,Form("#chi^{2}/n.d.f = %.2f",Chi2FitX)); y0 -= decrease;
	////t->DrawLatex(x0,y0,Form("Entries:%4.3e",_hMass[arm]->GetEntries())); y0 -= decrease;
	//t->DrawLatex(x0,y0,Form("#mu_{J/#psi} = %4.3f",Jpsi_mean)); y0 -= decrease;
	//t->DrawLatex(x0,y0,Form("#sigma_{J/#psi} = %4.3f",Jpsi_sigma)); y0 -= decrease;
	//t->DrawLatex(x0,y0,Form("N_{J/#psi} = %4.3e",N_Jpsi)); y0 -= decrease;
	//t->DrawLatex(x0,y0,Form("N_{#psi'} = %4.3e",N_Psip)); y0 -= decrease;
	////t->DrawLatex(x0,y0,Form("N_{BKG.}^{2#sigma} = %4.3e",N_BKG+N_Psip)); y0 -= decrease;
	////t->DrawLatex(x0,y0,Form("#frac{N_{sig.}}{N_{BKG.}} = %4.3f",S2B)); y0 -= decrease;
	//t->DrawLatex(x0,y0,Form("r^{2#sigma} = %4.3f #pm %4.3f",r_2sigma,er_2sigma)); y0 -= decrease;
	//t->DrawLatex(x0,y0,Form("r^{3#sigma} = %4.3f #pm %4.3f",r_3sigma,er_3sigma)); y0 -= decrease;
	//t->DrawLatex(x0,y0,Form("2#sigma window: (%4.3f,%4.3f)",Jpsi_mean-2*Jpsi_sigma,Jpsi_mean+2*Jpsi_sigma)); y0 -= decrease;
	//t->DrawLatex(x0,y0,Form("3#sigma window: (%4.3f,%4.3f)",Jpsi_mean-3*Jpsi_sigma,Jpsi_mean+3*Jpsi_sigma)); y0 -= decrease;

	if(is_for_paper)
	{
		TLegend *leg = new TLegend(0.6,0.88,0.88,0.6);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->AddEntry(frame->findObject("h_roohdata"),"Data","lep");
		leg->AddEntry(training_zone_1,"Training Zone","f");
		leg->AddEntry(frame->findObject("h_rooho"),"GPR BKG. estimation","lep");
		leg->AddEntry(frame->findObject("h_roohextract"),"Data after BKG. extraction","lep");
		leg->AddEntry(frame->findObject("model_Norm[x]"),"J/#psi + #psi'","l");
		leg->AddEntry(frame->findObject("Jpsi_Norm[x]"),"J/#psi","l");
		leg->AddEntry(frame->findObject("Psip_Norm[x]"),"#psi'","l");
		c0->cd();
		leg->Draw();
	}
	else
	{
		TPaveText *pt = new TPaveText(0.6,0.55,0.9,0.9,"blNDC");
		pt->SetBorderSize(1);
		pt->SetFillColor(0);
		pt->SetTextFont(42);
		pt->AddText(Form("#chi^{2}/n.d.f = %.2f",Chi2FitX));
		pt->AddText(Form("#mu_{J/#psi} = %4.3f",Jpsi_mean));
		pt->AddText(Form("#sigma_{J/#psi} = %4.3f",Jpsi_sigma));
		pt->AddText(Form("N_{ALL}^{2#sigma} = %4.3e",N_ALL_2sigma));
		pt->AddText(Form("N_{ALL}^{3#sigma} = %4.3e",N_ALL_3sigma));
		pt->AddText(Form("N_{J/#psi} = %4.3e",N_Jpsi));
		pt->AddText(Form("N_{#psi'} = %4.3e",N_Psip));
		pt->AddText(Form("r^{2#sigma} = %4.3f #pm %4.3f",r_2sigma,er_2sigma));
		pt->AddText(Form("r^{3#sigma} = %4.3f #pm %4.3f",r_3sigma,er_3sigma));
		pt->AddText(Form("2#sigma window: (%4.3f,%4.3f)",Jpsi_mean-2*Jpsi_sigma,Jpsi_mean+2*Jpsi_sigma));
		pt->AddText(Form("3#sigma window: (%4.3f,%4.3f)",Jpsi_mean-3*Jpsi_sigma,Jpsi_mean+3*Jpsi_sigma));
		c0->cd();
		pt->Draw();
		pt->Print();
	}

	c0->Update();	

	cout<<"=========================================================================="<<endl;
	w0->Print("t");
	fr->Print();
	cout<<"Chi2FitX: "<<Chi2FitX<<endl;
	cout<<"minNLL: "<< fr->minNll() << endl;
	cout<<"=========================================================================="<<endl;

	std::ofstream fitparam;
	fitparam.open("fitparam_roofit.dat",std::ofstream::out | std::ofstream::app);
	fitparam
		<<arm <<" "
		<<binning_mode<<" "
		<<bin_min <<" "<<bin_max<<" "
		<<r_2sigma<<" "<<r_3sigma<<" "
		<<Jpsi_mean-2*Jpsi_sigma<<" "<<Jpsi_mean+2*Jpsi_sigma<<" "
		<<Jpsi_mean-3*Jpsi_sigma<<" "<<Jpsi_mean+3*Jpsi_sigma<<" "
		<<endl;
	fitparam.close();

	dataFile->Close();

	//w->writeToFile(Form("model_%d_%s.root",model_choice,name_core.data()));
	//TFile *fout = TFile::Open(Form("%s_fit.root",name_core.data()),"recreate");
	fout->cd();
	c0->Write();
	fout->Close();

	return true;
}

//---------------------------------------------------------------------------------------------------
bool gprFit::fit
(
 GausProc *gaus_proc
 )
{
	//const float JpsiMass = 3.096916;
	//const float PsiMass  = 3.68609;
	//const float NSIGMA = 2.;

	//double xmin = 1.5;
	//double xmax = 6.0;
	//double xbinw = 0.05;//0.025;
	Int_t xnbin= TMath::Nint((xmax - xmin)/xbinw);

	double zone1_min = 1.5;
	double zone1_max = 2.2;
	double zone2_min = 4.3;
	double zone2_max = 6.0;

	const std::string name_core(Form("BinningMode_%d_Arm%d_Charge%d_bin%03.1f_%03.1f",binning_mode,arm,charge_req,bin_min,bin_max));
	gROOT->SetStyle("Plain");

	char *f0 = Form("%s_GPR_EXTRACT.root",name_core.data());
	TFile *dataFile = TFile::Open(f0,"READ");

	TFile *fout = TFile::Open(Form("%s_fit.root",name_core.data()),"recreate");

	TH1D *hdata = (TH1D *) dataFile->Get("hdata");
	TH1D *hextract = (TH1D *) dataFile->Get("hextract");
	TH1F *ho = (TH1F *) dataFile->Get("ho");
	TH1F *ho_shift = (TH1F *) dataFile->Get("ho_shift_back");	

	if(plot_ho_shift_back)
		ho = (TH1F *) dataFile->Get("ho_shift_back");
	TH1D *training_zone_1 = hmanip::RangeClone(hdata,zone1_min,zone1_max);
	TH1D *training_zone_2 = hmanip::RangeClone(hdata,zone2_min,zone2_max);
	TH1D *side_band_zone = hmanip::RangeClone(hdata,side_band_range_min,side_band_range_max);

	TCanvas* c0 = new TCanvas((name_core+"_fit").data(),(name_core+"_fit").data(),100,100,1000,800);
	c0->SetLogy();

	TF1 *model = NULL;
	int n_par = 0;

	switch(model_choice)
	{
		case 1://Jpsi: Double Gaussian, Psi(2s): Gaussian
			model = new TF1("model",Form("[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[1],[4],1) +  [5]*TMath::Gaus(x,%f*[1],%f*[2],1)",PsiMass/JpsiMass,PsiMass/JpsiMass),xmin,xmax);
			model->SetParNames(
					"N_Jpsi",
					"mu",
					"sigma",
					"N_Jpsi_tail",
					"sigma_tail",
					"N_Psip"
					);
			model->SetParameters(50000,3.09,0.18,1000,0.25,1000);
			model->SetParLimits(3,0,10000);
			model->SetParLimits(4,2.,100);
			model->SetParLimits(5,0,10000);
			break;
		case 2://Jpsi, Psi(2s): Gaussian
			model = new TF1("model",Form("[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,%f*[1],%f*[2],1)",PsiMass/JpsiMass,PsiMass/JpsiMass),xmin,xmax);
			model->SetParNames(
					"N_Jpsi",
					"mu",
					"sigma",
					"N_Psip"
					);
			model->SetParameters(50000,3.09,0.18,1000);
			break;
		case 3://Jpsi, Psi(2s): CB
			n_par = 8;
			model = new TF1("model",CB_plus_CB, xmin, xmax, n_par);
			model->SetParNames(
					"alpha",
					"n",
					"mu",
					"sigma",
					"N_Jpsi",
					"alpha_Psip",
					"n_Psip",
					"N_Psip"
					);
			model->SetParameters(2,1,3.09,0.18,50000,
					2,
					1,
					200);

			model->SetParLimits(0,1,10);
			model->SetParLimits(1,0.001,1000);
			//model->SetParLimits(2,3.,4.);
			model->SetParLimits(3,0.1,0.3);
			//model->SetParLimits(4,0,10000);
			model->SetParLimits(5,1,10);
			model->SetParLimits(6,0,1000);
			model->SetParLimits(7,0,1000);
			break;
		case 5://Jpsi, Psi(2s): CB psip has same tail as Jpsi
			n_par = 6;
			model = new TF1("model",CB_plus_CB_same_Jpsi_tail, xmin, xmax, n_par);
			model->SetParNames(
					"alpha",
					"n",
					"mu",
					"sigma",
					"N_Jpsi",
					"N_Psip"
					);
			model->SetParameters(2,1,3.09,0.18,50000,
					200);

			model->SetParLimits(0,1,10);
			model->SetParLimits(1,0.001,1000);
			//model->SetParLimits(2,3.,4.);
			model->SetParLimits(3,0.1,0.3);
			//model->SetParLimits(4,0,10000);
			model->SetParLimits(5,0,1000);

			model->FixParameter(1,100);
			break;
		case 4://Jpsi:CB, Psi(2s): Gaussian
			n_par = 6;
			model = new TF1("model",CB_plus_Gaussian, xmin, xmax, n_par);
			model->SetParNames(
					"alpha", "n", "mu",	"sigma", "N_Jpsi",
					"N_Psip"
					);
			model->SetParameters(
					2,1,3.09,0.18,50000,
					200);

			model->SetParLimits(0,1,10);
			model->SetParLimits(1,0.01,100);
			//model->SetParLimits(2,3.,4.);
			model->SetParLimits(3,0.1,0.3);
			//model->SetParLimits(4,0,10000);
			model->SetParLimits(5,0,1000);
			//model = tf1_CB_plus_Gaussian();
			break;
		case 6://CB plus Gaussian without j/psi and psip peak correlation
		n_par = 8;
		model = new TF1("model",CB_plus_Gaussian_nopeak_correlated,xmin,xmax, n_par);
		model->SetParNames("alpha","n","mu","sigma","N_Jpsi","N_Psip","mup","sigmap");
		model->SetParameters(2,1,3.09,0.18,50000,100,3.80,0.22);
		
		model->SetParLimits(0,1,10);
		model->SetParLimits(1,0.01,100);
		model->SetParLimits(3,0.1,0.3);
		model->SetParLimits(6,3.5,3.9);
		model->SetParLimits(5,10,1000);
		model->SetParLimits(7,0.15,0.30);
		break;
		default:
			cout<<"Error, bad input, quitting\n";
			break;
	}

	TVirtualFitter* fitter = TVirtualFitter::Fitter(hextract);
	fitter->SetPrecision(0.1);

	TFitResultPtr fr = hextract->Fit(model, "S");
	//TFitResultPtr fr = hextract->Fit(model, "VSL");





	//For pull distribution
	int Nbins = hextract->GetNbinsX();
	double data_error, Ndata, Nfit;
	double pull_value;

	TH1F *pull_distri = new TH1F(Form("%s_Pull_distri",name_core.data()),Form("%s_Pull_distri",name_core.data()), Nbins, xmin, xmax);
 	TH1F *pull = new TH1F(Form("%s_Pull_distribution",name_core.data()),Form("%s_Pull_distribution",name_core.data()),60,-6,6);	
	TH1F *fit_curve = new TH1F(Form("%s_fit_curve",name_core.data()), Form("%s_fit_curve",name_core.data()), Nbins, xmin, xmax);
	for (int i = 1; i<=Nbins; i++)
	{
		Ndata = hdata->GetBinContent(i);
		data_error = sqrt(Ndata);
//		Nfit = ho_shift->GetBinContent(i) + model->Eval(1.5+0.05*(i-1));
		Nfit = ho->GetBinContent(i) + model->Eval(1.5+0.05*(i-1./2));
		fit_curve->SetBinContent(i,Nfit);
		pull_value = (Nfit - Ndata)/data_error;
		pull_distri->SetBinContent(i,pull_value);

		pull->GetXaxis()->SetTitle("(N_{data}-N_{fit})/#deltaN_{data}");
		pull->Fill(pull_value);
		
	}
	pull_distri->SetFillColor(4);
	pull_distri->SetBarWidth(1);
	pull_distri->GetYaxis()->SetTitle("(N_{data}-N_{fit})/#deltaN_{data}");
	pull_distri->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}}[GeV/c^{2}]");	
	TCanvas *c1 = new TCanvas(Form("%s_Pull_Canvas",name_core.data()),Form("%s_Pull_Canvas",name_core.data()),1000,400);
	c1->Divide(2,1);
	c1->cd(1);
	pull_distri->Draw("b");
	c1->cd(2);
	pull->Draw();
	fout->cd();
	c1->Write();










	TF1 *tf_Jpsi = NULL;
	TF1 *tf_Psip = NULL;
	double Jpsi_mean = 0;
	double Jpsi_sigma = 0;
	if(model_choice == 1){
		Jpsi_mean = model->GetParameter("mu");
		Jpsi_sigma = model->GetParameter("sigma");
		tf_Jpsi = new TF1("tf_Jpsi","[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[1],[4],1)",xmin,xmax);
		for(int i=0;i<5;i++)
			tf_Jpsi->SetParameter(i,model->GetParameter(i));
		tf_Psip = new TF1("tf_Psip","[0]*TMath::Gaus(x,[1],[2],1)",xmin,xmax);
		tf_Psip->SetParameter(0,model->GetParameter("N_Psip"));
		tf_Psip->SetParameter(1,PsiMass/JpsiMass*Jpsi_mean);
		tf_Psip->SetParameter(2,PsiMass/JpsiMass*Jpsi_sigma);
	} else	if(model_choice == 2){
		Jpsi_mean = model->GetParameter("mu");
		Jpsi_sigma = model->GetParameter("sigma");
		tf_Jpsi = new TF1("tf_Jpsi","[0]*TMath::Gaus(x,[1],[2],1)",xmin,xmax);
		tf_Jpsi->SetParameter(0,model->GetParameter("N_Psip"));
		tf_Jpsi->SetParameter(1,Jpsi_mean);
		tf_Jpsi->SetParameter(2,Jpsi_sigma);
		for(int i=0;i<3;i++)
			tf_Jpsi->SetParameter(i,model->GetParameter(i));
		tf_Psip = new TF1("tf_Psip","[0]*TMath::Gaus(x,[1],[2],1)",xmin,xmax);
		tf_Psip->SetParameter(0,model->GetParameter("N_Psip"));
		tf_Psip->SetParameter(1,PsiMass/JpsiMass*Jpsi_mean);
		tf_Psip->SetParameter(2,PsiMass/JpsiMass*Jpsi_sigma);
	}	else if(model_choice == 3){
		Jpsi_mean = model->GetParameter("mu");
		Jpsi_sigma = model->GetParameter("sigma");
		tf_Jpsi = new TF1("tf_Jpsi",CrystalBall,xmin,xmax,5);
		for(int i=0;i<5;i++)
			tf_Jpsi->SetParameter(i,model->GetParameter(i));
		tf_Psip = new TF1("tf_Psip",CrystalBall,xmin,xmax,5);
		tf_Psip->SetParameter(5,model->GetParameter(0));
		tf_Psip->SetParameter(6,model->GetParameter(1));
		tf_Psip->SetParameter(2,PsiMass/JpsiMass*Jpsi_mean);
		tf_Psip->SetParameter(3,PsiMass/JpsiMass*Jpsi_sigma);
		tf_Psip->SetParameter(7,model->GetParameter("N_Psip"));
	}	else if(model_choice == 5){
		Jpsi_mean = model->GetParameter("mu");
		Jpsi_sigma = model->GetParameter("sigma");
		tf_Jpsi = new TF1("tf_Jpsi",CrystalBall,xmin,xmax,5);
		for(int i=0;i<5;i++)
			tf_Jpsi->SetParameter(i,model->GetParameter(i));
		tf_Psip = new TF1("tf_Psip",CrystalBall,xmin,xmax,5);
		tf_Psip->SetParameter(0,model->GetParameter(0));
		tf_Psip->SetParameter(1,model->GetParameter(1));
		tf_Psip->SetParameter(2,PsiMass/JpsiMass*Jpsi_mean);
		tf_Psip->SetParameter(3,PsiMass/JpsiMass*Jpsi_sigma);
		tf_Psip->SetParameter(5,model->GetParameter("N_Psip"));
		//__ERROR__;
		//return false;
	}	else if(model_choice == 4){
		Jpsi_mean = model->GetParameter("mu");
		Jpsi_sigma = model->GetParameter("sigma");
		tf_Jpsi = new TF1("tf_Jpsi",CrystalBall,xmin,xmax,5);
		for(int i=0;i<5;i++)
			tf_Jpsi->SetParameter(i,model->GetParameter(i));
		tf_Psip = new TF1("tf_Psip","[0]*TMath::Gaus(x,[1],[2],1)",xmin,xmax);
		tf_Psip->SetParameter(0,model->GetParameter("N_Psip"));
		tf_Psip->SetParameter(1,PsiMass/JpsiMass*Jpsi_mean);
		tf_Psip->SetParameter(2,PsiMass/JpsiMass*Jpsi_sigma);
	}
		else if(model_choice == 6){
		Jpsi_mean = model->GetParameter("mu");
                Jpsi_sigma = model->GetParameter("sigma");
		tf_Jpsi = new TF1("tf_Jpsi",CrystalBall,xmin,xmax,5);
		for(int i=0;i<5;i++)
                        tf_Jpsi->SetParameter(i,model->GetParameter(i));
		tf_Psip = new TF1("tf_Psip","[0]*TMath::Gaus(x,[1],[2],1)",xmin,xmax);
                tf_Psip->SetParameter(0,model->GetParameter("N_Psip"));
		tf_Psip->SetParameter(1,model->GetParameter("mup"));
                tf_Psip->SetParameter(2,model->GetParameter("sigmap"));
	}
	double xmin_win_2sigma = Jpsi_mean - 2. * Jpsi_sigma;
	double xmax_win_2sigma = Jpsi_mean + 2. * Jpsi_sigma;
	double xmin_win_3sigma = Jpsi_mean - 3. * Jpsi_sigma;
	double xmax_win_3sigma = Jpsi_mean + 3. * Jpsi_sigma;

	if(useFixedMassWin)
	{
		xmin_win_2sigma = 2.76;
		xmax_win_2sigma = 3.50;
		xmin_win_3sigma = 2.58;
		xmax_win_3sigma = 3.68;
	}

	hextract->SetLineColor(kRed);
	hextract->SetMarkerStyle(25);
	hextract->SetMarkerColor(kRed);
	hextract->SetMinimum(0.1);
	//hextract->SetMaximum(100000);
	//hextract->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} / GeV");
	//hextract->GetYaxis()->SetTitle(Form("Events / %.0fMeV",xbinw*1000));
	hextract->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}}[GeV/c^{2}]");
	hextract->GetYaxis()->SetTitle(Form("dN/dM_{#mu^{+}#mu^{-}} [(%.0f MeV/c^{2})^{-1}]",xbinw*1000));
	hextract->GetYaxis()->SetRangeUser(
			10,
			2.5*hextract->GetMaximum());
	if(is_for_paper)
		hextract->SetTitle("");
	else
		hextract->SetTitle(name_core.data());
	hextract->SetStats(0);
	hextract->Draw("9e1x0p");
	//hextract->Print("all");

	model->SetLineColor(kRed);
	model->Draw("same");

	tf_Jpsi->SetLineColor(kGreen);
	tf_Jpsi->SetLineStyle(7);
	tf_Jpsi->Draw("same");

	c0->cd();
	tf_Psip->SetLineColor(kBlue);
	tf_Psip->SetLineStyle(5);
	//tf_Psip->Print();
	//cout
	//	<<"DEBUG: tf_Psip->Eval: "
	//	<<" "<<tf_Psip->Eval(3.7)
	//	<<" "<<tf_Psip->Eval(3.8)
	//	<<" "<<tf_Psip->Eval(3.9)
	//	<<endl;
	tf_Psip->Draw("same");
	//tf_Psip->Draw();
	c0->Update();

	//Int_t nFitParam = 0;
	//if(model_choice==1) nFitParam = 9;
	//if(model_choice==4) nFitParam = 6;

	hdata->SetLineColor(kBlack);
	hdata->SetMarkerStyle(4);
	hdata->SetMarkerColor(kBlack);
	hdata->SetMinimum(1);
	hdata->SetStats(0);
	hdata->Draw("9e1x0psame");

	ho->SetLineColor(kBlue);
	ho->SetMarkerStyle(26);
	ho->SetMarkerSize(1.);
	ho->SetMarkerColor(kBlue);
	ho->SetMinimum(1);
	ho->Draw("9e1x0psame");

	fit_curve->Draw("csame");

	double N_Jpsi = model->Integral(xmin,xmax)/xbinw;
	double N_Psip = model->GetParameter("N_Psip")/xbinw;

	double N_Jpsi_2sigma_histo_count = hextract->Integral(
			hextract->FindBin(xmin_win_2sigma),
			hextract->FindBin(xmax_win_2sigma)
			);
	double N_Jpsi_2sigma = model->Integral(xmin_win_2sigma,xmax_win_2sigma)/xbinw;
	double eN_Jpsi_2sigma = model->IntegralError(xmin_win_2sigma,xmax_win_2sigma,
			model->GetParameters(),
			fr.Get()->GetCovarianceMatrix().GetMatrixArray())/xbinw;

	gaus_proc->Integral(xmin_win_2sigma,xmax_win_2sigma, integral_bkg, integral_error_bkg);
	double N_BKG_2sigma = integral_bkg/xbinw;
	double eN_BKG_2sigma = integral_error_bkg/xbinw;

	double r_2sigma = 1./(1. + N_Jpsi_2sigma/N_BKG_2sigma);
	double er_2sigma = (r_2sigma*r_2sigma) * (N_Jpsi_2sigma/N_BKG_2sigma) *
		sqrt( pow(eN_Jpsi_2sigma/N_Jpsi_2sigma,2) + pow(eN_BKG_2sigma/N_BKG_2sigma,2));

	cout
		<<"DEBUG: method 1: -----------------------------"<<endl
		<<" N_Jpsi_2sigma_histo_count: "<<N_Jpsi_2sigma_histo_count
		<<" N_Jpsi_2sigma: "<<N_Jpsi_2sigma
		<<" eN_Jpsi_2sigma: "<<eN_Jpsi_2sigma
		<<" N_BKG_2sigma: "<<N_BKG_2sigma
		<<" eN_BKG_2sigma: "<<eN_BKG_2sigma
		<<endl;

	if(r_extraction_method == 2)
	{
		double N_Incl = hdata->Integral(
				hdata->FindBin(xmin_win_2sigma),
				hdata->FindBin(xmax_win_2sigma)
				);
		r_2sigma = N_BKG_2sigma/N_Incl;
		er_2sigma = r_2sigma*sqrt(
				pow(eN_BKG_2sigma/N_BKG_2sigma,2) +
				1.0/N_Incl
				);
	}

	cout
		<<"DEBUG: method 2: -----------------------------"<<endl
		<<" N_Jpsi_2sigma_histo_count: "<<N_Jpsi_2sigma_histo_count
		<<" N_Jpsi_2sigma: "<<N_Jpsi_2sigma
		<<" eN_Jpsi_2sigma: "<<eN_Jpsi_2sigma
		<<" N_BKG_2sigma: "<<N_BKG_2sigma
		<<" eN_BKG_2sigma: "<<eN_BKG_2sigma
		<<endl;

	model->Print();

	fr.Get()->GetCovarianceMatrix().Print();//DEBUG

	double N_Jpsi_3sigma = model->Integral(xmin_win_3sigma,xmax_win_3sigma)/xbinw;
	double eN_Jpsi_3sigma = model->IntegralError(xmin_win_3sigma,xmax_win_3sigma,
			model->GetParameters(),
			fr.Get()->GetCovarianceMatrix().GetMatrixArray())/xbinw;

	gaus_proc->Integral(xmin_win_3sigma,xmax_win_3sigma, integral_bkg, integral_error_bkg);
	double N_BKG_3sigma = integral_bkg/xbinw;
	double eN_BKG_3sigma = integral_error_bkg/xbinw;

	double r_3sigma = 1./(1. + N_Jpsi_3sigma/N_BKG_3sigma);
	double er_3sigma = (r_3sigma*r_3sigma) * (N_Jpsi_3sigma/N_BKG_3sigma) *
		sqrt( pow(eN_Jpsi_3sigma/N_Jpsi_3sigma,2) + pow(eN_BKG_3sigma/N_BKG_3sigma,2));

	if(r_extraction_method == 2)
	{
		double N_Incl = hdata->Integral(
				hdata->FindBin(xmin_win_3sigma),
				hdata->FindBin(xmax_win_3sigma)
				);
		r_3sigma = N_BKG_3sigma/N_Incl;
		er_3sigma = r_3sigma*sqrt(
				pow(eN_BKG_3sigma/N_BKG_3sigma,2) +
				1.0/N_Incl
				);
	}

	bool draw_nsigma_line = true;
	if(is_for_paper)
		draw_nsigma_line = false;

	TLine *lLowSig_2sigma = new TLine(xmin_win_2sigma,0,xmin_win_2sigma,hextract->GetMaximum());
	lLowSig_2sigma->SetLineStyle(5);
	lLowSig_2sigma->SetLineColor(kRed);
	lLowSig_2sigma->Draw();
	TLine *lHigSig_2sigma = new TLine(xmax_win_2sigma,0,xmax_win_2sigma,hextract->GetMaximum());
	lHigSig_2sigma->SetLineStyle(5);
	lHigSig_2sigma->SetLineColor(kRed);
	lHigSig_2sigma->Draw();

	if(draw_nsigma_line)
	{
		TLine *lLowSig_3sigma = new TLine(xmin_win_3sigma,0,xmin_win_3sigma,hextract->GetMaximum());
		lLowSig_3sigma->SetLineColor(kRed);
		lLowSig_3sigma->SetLineStyle(2);
		lLowSig_3sigma->Draw();
		TLine *lHigSig_3sigma = new TLine(xmax_win_3sigma,0,xmax_win_3sigma,hextract->GetMaximum());
		lHigSig_3sigma->SetLineColor(kRed);
		lHigSig_3sigma->SetLineStyle(2);
		lHigSig_3sigma->Draw();
	}

	double hist_min = hextract->GetMinimum();
	double hist_max = hextract->GetMaximum();

	if(is_for_paper)
	{
		side_band_zone->SetFillColor(kGreen);
		side_band_zone->SetFillStyle(3013);
		side_band_zone->Draw("same");
	} else {

		training_zone_1->SetFillColor(kGreen);
		training_zone_1->SetFillStyle(3013);
		training_zone_1->Draw("same");

		training_zone_2->SetFillColor(kGreen);
		training_zone_2->SetFillStyle(3013);
		training_zone_2->Draw("same");
	}

	if(is_for_paper)
	{
		TLegend *leg = new TLegend(0.6,0.88,0.88,0.6);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->AddEntry(hdata,"Data","lep");
		//leg->AddEntry(training_zone_1,"Training Zone","f");
		leg->AddEntry(side_band_zone,"Sideband region","f");
		leg->AddEntry(ho,"GPR BKG. estimation","lep");
		leg->AddEntry(hextract,"Data after BKG. extraction","lep");
		leg->AddEntry(model,"J/#psi + #psi'","l");
		leg->AddEntry(tf_Jpsi,"J/#psi","l");
		leg->AddEntry(tf_Psip,"#psi'","l");
		leg->Draw();
	}
	else
	{
		TPaveText *pt = new TPaveText(0.6,0.55,0.9,0.9,"blNDC");
		pt->SetBorderSize(1);
		pt->SetFillColor(0);
		pt->SetTextFont(42);
		pt->AddText(Form("#chi^{2}/n.d.f=%.2f/%d=%4.3f",fr.Get()->Chi2(), fr.Get()->Ndf(),fr.Get()->Chi2()/fr.Get()->Ndf()));
		pt->AddText(Form("#mu_{J/#psi} = %4.3f",Jpsi_mean));
		pt->AddText(Form("#sigma_{J/#psi} = %4.3f",Jpsi_sigma));
		pt->AddText(Form("N_{ALL}^{2#sigma} = %4.3e",N_Jpsi_2sigma + N_BKG_2sigma));
		pt->AddText(Form("N_{ALL}^{2#sigma} = %4.3e",N_Jpsi_2sigma + N_BKG_2sigma));
		pt->AddText(Form("N_{ALL}^{3#sigma} = %4.3e",N_Jpsi_3sigma + N_BKG_3sigma));
		pt->AddText(Form("N_{J/#psi}^{count} = %4.3e",N_Jpsi_2sigma_histo_count));
		pt->AddText(Form("N_{J/#psi} = %4.3e",N_Jpsi));
		pt->AddText(Form("N_{#psi'} = %4.3e",N_Psip));
		pt->AddText(Form("r^{2#sigma} = %4.3f #pm %4.3f",r_2sigma,er_2sigma));
		pt->AddText(Form("r^{3#sigma} = %4.3f #pm %4.3f",r_3sigma,er_3sigma));
		pt->AddText(Form("2#sigma window: (%4.3f,%4.3f)",Jpsi_mean-2*Jpsi_sigma,Jpsi_mean+2*Jpsi_sigma));
		pt->AddText(Form("3#sigma window: (%4.3f,%4.3f)",Jpsi_mean-3*Jpsi_sigma,Jpsi_mean+3*Jpsi_sigma));
		pt->Draw();
	}

	c0->Update();	

	cout<<"=========================================================================="<<endl;
	fr.Get()->Print();
	cout<<"Chi2/n.d.f. "<<fr.Get()->Chi2() / fr.Get()->Ndf()<<endl;
	//cout<<"minNLL: "<< fr->minNll() << endl;
	cout<<"=========================================================================="<<endl;

	std::ofstream fitparam;
	fitparam.open("fitparam_fit.dat",std::ofstream::out | std::ofstream::app);
	fitparam
		<<arm <<" "
		<<binning_mode<<" "
		<<bin_min <<" "<<bin_max<<" "
		<<r_2sigma<<" "<<r_3sigma<<" "
		//<<Jpsi_mean-2*Jpsi_sigma<<" "<<Jpsi_mean+2*Jpsi_sigma<<" "
		//<<Jpsi_mean-3*Jpsi_sigma<<" "<<Jpsi_mean+3*Jpsi_sigma<<" "
		<<xmin_win_2sigma<<" "<<xmax_win_2sigma<<" "
		<<xmin_win_3sigma<<" "<<xmax_win_3sigma<<" "
		<<endl;
	fitparam.close();


	std::ofstream fiterr;
	fiterr.open("fiterr.dat",std::ofstream::out | std::ofstream::app);
	if(binning_mode == 1)
		fiterr<<"pT"<<" \t ";
	else if(binning_mode == 2)
		fiterr<<"eta"<<" \t ";

	fiterr
		<<arm<<" "
		<<Form("%5.3f",bin_min)<<" "
		<<Form("%5.3f",bin_max)<<" "
		<<Form("%5.3f",r_2sigma)<<" "
		<<Form("%5.3f",er_2sigma)<<" "
		<<Form("%5.3f",er_2sigma_syst)<<" "
		<<endl;
	fiterr.close();

	//w->writeToFile(Form("model_%d_%s.root",model_choice,name_core.data()));
	//TFile *fout = TFile::Open(Form("%s_fit.root",name_core.data()),"recreate");
	c0->Write();
	fout->Close();

	dataFile->Close();

	return true;
}

bool gprFit::function_formed_fit(
		GausProc *gaus_proc
		)
{
	Int_t xnbin= TMath::Nint((xmax - xmin)/xbinw);

	double zone1_min = 1.5;
	double zone1_max = 2.2;
	double zone2_min = 4.3;
	double zone2_max = 6.0;

	const std::string name_core(Form("BinningMode_%d_Arm%d_Charge%d_bin%03.1f_%03.1f",binning_mode,arm,charge_req,bin_min,bin_max));
	gROOT->SetStyle("Plain");

	char *f0 = Form("%s_GPR_EXTRACT.root",name_core.data());
	TFile *dataFile = TFile::Open(f0,"READ");

	TFile *fout = TFile::Open(Form("%s_fit.root",name_core.data()),"recreate");

	TH1D *hdata = (TH1D *) dataFile->Get("hdata");
	TH1D *hextract = (TH1D *) dataFile->Get("hextract");
	TH1F *ho = (TH1F *) dataFile->Get("ho");
	if(plot_ho_shift_back)
		ho = (TH1F *) dataFile->Get("ho_shift_back");
	TH1D *training_zone_1 = hmanip::RangeClone(hdata,zone1_min,zone1_max);
	TH1D *training_zone_2 = hmanip::RangeClone(hdata,zone2_min,zone2_max);
	TH1D *side_band_zone = hmanip::RangeClone(hdata,side_band_range_min,side_band_range_max);

	TCanvas* c0 = new TCanvas((name_core+"_fit").data(),(name_core+"_fit").data(),100,100,1000,800);
	c0->SetLogy();

	TF1 *model = NULL;
	int n_par = 0;

	double xmin_form = 1.5;
	double xmax_form = 5.0;

	switch(model_choice)
	{
		case 1://Jpsi: Double Gaussian, Psi(2s): Gaussian
			break;
		case 2://Jpsi, Psi(2s): Gaussian; Bkg: pol3
			n_par = 8;
			model = new TF1("model",Gaussian_plus_Gaussian_plus_pol3,xmin_form,xmax_form,n_par);
			model->SetParNames(
					"N_Jpsi", "mu",	"sigma",
					"N_Psip",
					"p0","p1","p2","p3"
					);
			model->SetParameters(
					2000,3.09,0.18,
					200,
					0,0,0,0	
					);
			model->SetParLimits(3,10,1000);
			break;
		case 3://Jpsi, Psi(2s): CB
			break;
		case 5://Jpsi, Psi(2s): CB psip has same tail as Jpsi
			break;
		case 4://Jpsi:CB, Psi(2s): Gaussian; Bkg: pol3
			break;
		default:
			break;
	}

	//TVirtualFitter* fitter = TVirtualFitter::Fitter(hextract);
	//fitter->SetPrecision(0.1);

	//TFitResultPtr fr = hextract->Fit(model, "QSM");
	TFitResultPtr fr = hdata->Fit(model, "SR");

	TF1 *tf_Jpsi = NULL;
	TF1 *tf_Psip = NULL;
	TF1 *tf_Bkg = NULL;
	double Jpsi_mean = 0;
	double Jpsi_sigma = 0;

	switch(model_choice)
	{
		case 1:
			break;
		case 2:
			Jpsi_mean = model->GetParameter("mu");
			Jpsi_sigma = model->GetParameter("sigma");
			tf_Jpsi = new TF1("tf_Jpsi","[0]*TMath::Gaus(x,[1],[2],1)",xmin_form,xmax_form);
			tf_Jpsi->SetParameter(0,model->GetParameter("N_Psip"));
			tf_Jpsi->SetParameter(1,Jpsi_mean);
			tf_Jpsi->SetParameter(2,Jpsi_sigma);
			for(int i=0;i<3;i++)
				tf_Jpsi->SetParameter(i,model->GetParameter(i));
			tf_Psip = new TF1("tf_Psip","[0]*TMath::Gaus(x,[1],[2],1)",xmin_form,xmax_form);
			tf_Psip->SetParameter(0,model->GetParameter("N_Psip"));
			tf_Psip->SetParameter(1,PsiMass/JpsiMass*Jpsi_mean);
			tf_Psip->SetParameter(2,PsiMass/JpsiMass*Jpsi_sigma);
			tf_Bkg = new TF1("tf_Bkg","[0] + [1]*x + [2]*x*x + [3]*x*x*x",1.5,4.5);
			for(int i=0;i<4;i++)
				tf_Bkg->SetParameter(i,model->GetParameter(i+4));
			break;
		default:
			break;
	}

	double xmin_win_2sigma = Jpsi_mean - 2. * Jpsi_sigma;
	double xmax_win_2sigma = Jpsi_mean + 2. * Jpsi_sigma;
	double xmin_win_3sigma = Jpsi_mean - 3. * Jpsi_sigma;
	double xmax_win_3sigma = Jpsi_mean + 3. * Jpsi_sigma;

	c0->cd();

	hextract->SetLineColor(kRed);
	hextract->SetMarkerStyle(25);
	hextract->SetMarkerColor(kRed);
	hextract->SetMinimum(0.1);
	//hextract->SetMaximum(100000);
	//hextract->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} / GeV");
	//hextract->GetYaxis()->SetTitle(Form("Events / %.0fMeV",xbinw*1000));
	hextract->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}}[GeV/c^{2}]");
	hextract->GetYaxis()->SetTitle(Form("dN/dM_{#mu^{+}#mu^{-}} [(%.0f MeV/c^{2})^{-1}]",xbinw*1000));
	hextract->GetYaxis()->SetRangeUser(
			10,
			2.5*hextract->GetMaximum());
	if(is_for_paper)
		hextract->SetTitle("");
	else
		hextract->SetTitle(name_core.data());
	hextract->SetStats(0);
	hextract->Draw("9e1x0p");
	//hextract->Print("all");

	model->SetLineColor(kRed);
	model->Draw("same");

	tf_Jpsi->SetLineColor(kGreen);
	tf_Jpsi->SetLineStyle(7);
	tf_Jpsi->Draw("same");

	tf_Psip->SetLineColor(kBlue);
	tf_Psip->SetLineStyle(5);
	tf_Psip->Draw("same");

	tf_Bkg->SetLineColor(kRed);
	tf_Bkg->SetLineStyle(7);
	tf_Bkg->Draw("same");

	hdata->SetLineColor(kBlack);
	hdata->SetMarkerStyle(4);
	hdata->SetMarkerColor(kBlack);
	hdata->SetMinimum(1);
	hdata->SetStats(0);
	hdata->Draw("9e1x0psame");

	ho->SetLineColor(kBlue);
	ho->SetMarkerStyle(26);
	ho->SetMarkerSize(1.);
	ho->SetMarkerColor(kBlue);
	ho->SetMinimum(1);
	ho->Draw("9e1x0psame");

	c0->Update();

	double N_Jpsi = -999;
	double N_Psip = -999;

	double N_Jpsi_2sigma = -999;
	double eN_Jpsi_2sigma = -999;
	double N_BKG_2sigma = -999;
	double eN_BKG_2sigma = -999;
	double N_Incl_2sigma = -999;
	double eN_Incl_2sigma = -999;
	double r_2sigma = -999;
	double er_2sigma = -999;

	double N_Jpsi_3sigma = -999;
	double eN_Jpsi_3sigma = -999;
	double N_BKG_3sigma = -999;
	double eN_BKG_3sigma = -999;
	double N_Incl_3sigma = -999;
	double eN_Incl_3sigma = -999;
	double r_3sigma = -999;
	double er_3sigma = -999;

	N_Jpsi = tf_Jpsi->Integral(xmin,xmax)/xbinw;
	N_Psip = tf_Psip->Integral(xmin,xmax)/xbinw;

	//2sigma

	N_Jpsi_2sigma = tf_Jpsi->Integral(xmin_win_2sigma,xmax_win_2sigma)/xbinw;
	eN_Jpsi_2sigma = 1;

	N_Incl_2sigma = model->Integral(xmin_win_2sigma,xmax_win_2sigma)/xbinw;
	eN_Incl_2sigma = model->IntegralError(xmin_win_3sigma,xmax_win_3sigma,
			model->GetParameters(),
			fr.Get()->GetCovarianceMatrix().GetMatrixArray())/xbinw;

	N_BKG_2sigma = tf_Bkg->Integral(xmin_win_2sigma,xmax_win_2sigma)/xbinw;
	eN_BKG_2sigma = 1;

	r_2sigma = N_BKG_2sigma/N_Incl_2sigma;
	er_2sigma = (r_2sigma*r_2sigma) * (N_Jpsi_2sigma/N_BKG_2sigma) *
		sqrt( pow(eN_Jpsi_2sigma/N_Jpsi_2sigma,2) + pow(eN_BKG_2sigma/N_BKG_2sigma,2));

	//3sigma
	//...

	bool draw_nsigma_line = true;
	if(is_for_paper)
		draw_nsigma_line = false;

	TLine *lLowSig_2sigma = new TLine(xmin_win_2sigma,0,xmin_win_2sigma,hextract->GetMaximum());
	lLowSig_2sigma->SetLineStyle(5);
	lLowSig_2sigma->SetLineColor(kRed);
	lLowSig_2sigma->Draw();
	TLine *lHigSig_2sigma = new TLine(xmax_win_2sigma,0,xmax_win_2sigma,hextract->GetMaximum());
	lHigSig_2sigma->SetLineStyle(5);
	lHigSig_2sigma->SetLineColor(kRed);
	lHigSig_2sigma->Draw();

	if(draw_nsigma_line)
	{
		TLine *lLowSig_3sigma = new TLine(xmin_win_3sigma,0,xmin_win_3sigma,hextract->GetMaximum());
		lLowSig_3sigma->SetLineColor(kRed);
		lLowSig_3sigma->SetLineStyle(2);
		lLowSig_3sigma->Draw();
		TLine *lHigSig_3sigma = new TLine(xmax_win_3sigma,0,xmax_win_3sigma,hextract->GetMaximum());
		lHigSig_3sigma->SetLineColor(kRed);
		lHigSig_3sigma->SetLineStyle(2);
		lHigSig_3sigma->Draw();
	}

	double hist_min = hextract->GetMinimum();
	double hist_max = hextract->GetMaximum();

	if(is_for_paper)
	{
		side_band_zone->SetFillColor(kGreen);
		side_band_zone->SetFillStyle(3013);
		side_band_zone->Draw("same");
	} else {

		training_zone_1->SetFillColor(kGreen);
		training_zone_1->SetFillStyle(3013);
		training_zone_1->Draw("same");

		training_zone_2->SetFillColor(kGreen);
		training_zone_2->SetFillStyle(3013);
		training_zone_2->Draw("same");
	}

	if(is_for_paper)
	{
		TLegend *leg = new TLegend(0.6,0.88,0.88,0.6);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->AddEntry(hdata,"Data","lep");
		//leg->AddEntry(training_zone_1,"Training Zone","f");
		leg->AddEntry(side_band_zone,"Sideband region","f");
		leg->AddEntry(ho,"GPR BKG. estimation","lep");
		leg->AddEntry(hextract,"Data after BKG. extraction","lep");
		leg->AddEntry(model,"J/#psi + #psi'","l");
		leg->AddEntry(tf_Jpsi,"J/#psi","l");
		leg->AddEntry(tf_Psip,"#psi'","l");
		leg->Draw();
	}
	else
	{
		TPaveText *pt = new TPaveText(0.6,0.55,0.9,0.9,"blNDC");
		pt->SetBorderSize(1);
		pt->SetFillColor(0);
		pt->SetTextFont(42);
		pt->AddText(Form("#chi^{2}/n.d.f=%.2f/%d=%4.3f",fr.Get()->Chi2(), fr.Get()->Ndf(),fr.Get()->Chi2()/fr.Get()->Ndf()));
		pt->AddText(Form("#mu_{J/#psi} = %4.3f",Jpsi_mean));
		pt->AddText(Form("#sigma_{J/#psi} = %4.3f",Jpsi_sigma));
		pt->AddText(Form("N_{ALL}^{2#sigma} = %4.3e",N_Jpsi_2sigma + N_BKG_2sigma));
		pt->AddText(Form("N_{ALL}^{3#sigma} = %4.3e",N_Jpsi_3sigma + N_BKG_3sigma));
		pt->AddText(Form("N_{J/#psi} = %4.3e",N_Jpsi));
		pt->AddText(Form("N_{#psi'} = %4.3e",N_Psip));
		pt->AddText(Form("r^{2#sigma} = %4.3f #pm %4.3f",r_2sigma,er_2sigma));
		pt->AddText(Form("r^{3#sigma} = %4.3f #pm %4.3f",r_3sigma,er_3sigma));
		pt->AddText(Form("2#sigma window: (%4.3f,%4.3f)",Jpsi_mean-2*Jpsi_sigma,Jpsi_mean+2*Jpsi_sigma));
		pt->AddText(Form("3#sigma window: (%4.3f,%4.3f)",Jpsi_mean-3*Jpsi_sigma,Jpsi_mean+3*Jpsi_sigma));
		pt->Draw();
	}

	c0->Update();	

	cout<<"=========================================================================="<<endl;
	fr.Get()->Print();
	cout<<"Chi2/n.d.f. "<<fr.Get()->Chi2() / fr.Get()->Ndf()<<endl;
	//cout<<"minNLL: "<< fr->minNll() << endl;
	cout<<"=========================================================================="<<endl;

	std::ofstream fitparam;
	fitparam.open("fitparam_fit.dat",std::ofstream::out | std::ofstream::app);
	fitparam
		<<arm <<" "
		<<binning_mode<<" "
		<<bin_min <<" "<<bin_max<<" "
		<<r_2sigma<<" "<<r_3sigma<<" "
		<<Jpsi_mean-2*Jpsi_sigma<<" "<<Jpsi_mean+2*Jpsi_sigma<<" "
		<<Jpsi_mean-3*Jpsi_sigma<<" "<<Jpsi_mean+3*Jpsi_sigma<<" "
		<<endl;
	fitparam.close();


	std::ofstream fiterr;
	fiterr.open("fiterr.dat",std::ofstream::out | std::ofstream::app);
	if(binning_mode == 1)
		fiterr<<"pT"<<" \t ";
	else if(binning_mode == 2)
		fiterr<<"eta"<<" \t ";

	fiterr
		<<arm<<" "
		<<Form("%5.3f",bin_min)<<" "
		<<Form("%5.3f",bin_max)<<" "
		<<Form("%5.3f",r_2sigma)<<" "
		<<Form("%5.3f",er_2sigma)<<" "
		<<Form("%5.3f",er_2sigma_syst)<<" "
		<<endl;
	fiterr.close();

	//w->writeToFile(Form("model_%d_%s.root",model_choice,name_core.data()));
	//TFile *fout = TFile::Open(Form("%s_fit.root",name_core.data()),"recreate");
	fout->cd();
	c0->Write();
	fout->Close();

	dataFile->Close();

	return true;

	return true;
}

//---------------------------------------------------------------------------------------------------

//int main(int argc, char **argv)
//{
//	int arm = 0; //0: North, 1: South
//	int binning_mode = 1;//1: pT; 2: eta
//	double bin_min = 0;
//	double bin_max = 10.;
//	int charge_req = 0; // 0: charge_req, 1:like_sign
//	int model_choice = 4;
//	int use_fvtx = 0;
//	int is_for_paper = 0;
//
//	gprFit *myfit = new gprFit(
//			arm,
//			binning_mode,
//			bin_min,
//			bin_max,
//			charge_req,
//			model_choice,
//			use_fvtx,
//			is_for_paper
//			);
//
//	GausProc *d = myfit->gp();
//	myfit->gpr_extract();
//	myfit->roofit(d);
//	myfit->fit(d);
//
//	return 0;
//}

//---------------------------------------------------------------------------------------------------
