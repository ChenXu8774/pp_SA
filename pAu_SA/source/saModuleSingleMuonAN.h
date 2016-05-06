// $Id: saModuleSingleMuonAN.h,v 1.4 2016/02/28 05:32:50 jinhuang Exp $

/*!
 * \file saModuleSingleMuonAN.h
 * \brief Single Muon AN analysis.
 * \author Haiwang Yu <yuhw@rcf.rhic.bnl.gov>
 * \version $Revision: 1.4 $
 * \date $Date: 2016/02/28 05:32:50 $
 */

#ifndef _SAMODULESINGLEMUONAN_H_
#define _SAMODULESINGLEMUONAN_H_

//STL
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "saModuleBase.h"

//class DiMuon;
class SingleMuon;

class saModuleSingleMuonAN  : public saModuleBase
{
	public:
		saModuleSingleMuonAN(const std::string &name = "DimuonJpsiHaiwang", bool doControlData = false);
		virtual
			~saModuleSingleMuonAN();

		int verbosity;

		bool _use_bbc_cut;

		void set_use_bbc_cut(const bool a){_use_bbc_cut = a;}

		std::ofstream fout;

		bool _doControlData;
		int _FillCounter;
		double _pT_bins[4];

		enum EVENT_TYPE {MuM_N_LG4,MuM_N_LG23,MuM_S_LG4,MuM_S_LG23, BAD_EVENT};
		enum ARM_TYPE {NORTH, SOUTH};

#ifndef __CINT__

		//! global initialization
		virtual int
			init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

		//! Run initialization
		virtual int
			init_run(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

		//! event method
		virtual int
			event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

		//! global termination
		virtual int
			end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

		//! dimuon cut
		virtual bool
			pass_sngmuon_cut(const SingleMuon *dimuon);

		TH1 *makeTHist_pT_phi(
				std::string name,
				std::string title
				);
		saHist *makesaHist_pT_phi(
				sa_hist_mangager_ptr hm,
				TH1 *h
				);

#endif

		// TH2D bined by pT phi
		saHist * _sah_pT_phi_MuM_phi_N_LG4;  // mu-, North, last_gap = 4
                saHist * _sah_pT_phi_MuM_phi_N_LG23; // mu-, North, last_gap = 2,3
                saHist * _sah_pT_phi_MuM_phi_S_LG4;  // mu-, South, last_gap = 4
                saHist * _sah_pT_phi_MuM_phi_S_LG23; // mu-, South, last_gap = 2,3

		std::vector<saHist*> _v_sahist;
};

#endif /* _SAMODULESINGLEMUONAN_H_ */
