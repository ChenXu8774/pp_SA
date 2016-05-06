#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TF1.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooKeysPdf.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooExponential.h"
#include "RooProduct.h"

#include "hmanip.h"

#include "GausProc.h"
#include "GPOptimizer.h"

//---------------------------------------------------------------------------------------------------
class gprFit
{
	public:
		int arm; //0: North, 1: South
		int binning_mode;//1: pT; 2: eta
		double bin_min;
		double bin_max;
		int charge_req; // 0: charge_req, 1:like_sign
		int model_choice;
		int use_fvtx;
		int is_for_paper;

		bool plot_ho_shift_back;
		bool do_shift_bin_gpr;

		int r_extraction_method;
		bool useFixedMassWin;

		double er_2sigma_syst;
		//double xmin;
		//double xmax;
		//double xbinw;

		//float JpsiMass;
		//float PsiMass;
		//float NSIGMA;

		double integral_bkg;
		double integral_error_bkg;

		double side_band_range_min;
		double side_band_range_max;

		gprFit(
				int 		t_arm = 0, //0: North, 1: South
				int 		t_binning_mode = 1,//1: pT; 2: eta
				double 	t_bin_min = 0,
				double 	t_bin_max = 10.,
				int 		t_charge_req = 0, // 0: charge_req, 1:like_sign
				int 		t_model_choice = 4,
				int 		t_use_fvtx = 0,
				bool 		t_is_for_paper = false
				):
			arm(t_arm),
			binning_mode(t_binning_mode),
			bin_min(t_bin_min),
			bin_max(t_bin_max),
			charge_req(t_charge_req),
			model_choice(t_model_choice),
			use_fvtx(t_use_fvtx),
			is_for_paper(t_is_for_paper),
			plot_ho_shift_back(true),
			do_shift_bin_gpr(true),
			r_extraction_method(2), //1: fitting; 2: gpr + counts
			useFixedMassWin(false),
			side_band_range_min(1.5),
			side_band_range_max(2.4)
	{
		//const double JpsiMass = 3.096916;
		//const double PsiMass  = 3.68609;
		//NSIGMA = 2.;

		er_2sigma_syst = 0.05;

		integral_bkg = 0;
		integral_error_bkg = 0;
	}

	public:

		GausProc *gp();

		void gpr_extract();

		bool roofit(GausProc *gaus_proc);

		bool fit(GausProc *gaus_proc);

		bool function_formed_fit(GausProc *gaus_proc);
};

