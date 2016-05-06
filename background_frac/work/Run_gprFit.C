/*#include "gprFit.h"
#include "GausProc.h"
#include "GPFitNewton.h"
#include "GPOptimizer.h"
#include "GausProcLinkDef.h"
*/
void Run_gprFit(
		int binning_mode = 1,
		int arm = 0,
		double min = 0,
		double max = 10,
		bool is_for_paper = false
		)
{
	gSystem->Load("/gpfs/mnt/gpfs02/phenix/spin/spin2/chenxu/spinAnalyzer/background_frac/install/lib/libgprFit.so");

	gprFit * myfit = new gprFit();

	myfit->binning_mode = binning_mode;
	myfit->arm = arm;
	myfit->bin_min = min;
	myfit->bin_max = max;
	myfit->is_for_paper = is_for_paper;
	
	//myfit->er_2sigma_syst = 0.05;
	//myfit->useFixedMassWin = true;

	myfit->model_choice = 6;

	GausProc *d = myfit->gp();
	myfit->gpr_extract();

	myfit->fit(d);
	//myfit->function_formed_fit(d);
}
