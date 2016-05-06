/*!
 * \file merge_saHists_hist_by_hist.C
 * \brief merge files with saHists hist-by-hist, modified from Jin's merge_saHist.C, slower than it, comsumes less memory.
 * \author Haiwang Yu <yuhw.pku@gmail.com>
 * \version $Revision: 1.2 $
 * \date $Date: 2016/03/18 18:12:55 $
 */


#include <iostream>
#include <cassert>
#include <stdexcept>

using namespace std;

/*!
 * \param 
 */
bool combine_saHists(
		char *infilelist = "test.lst",//input
		char *outfile = ""//output
		)
{
	char *vname = "_m_";
	gSystem->Load("libspin_analyzer");

	if(outfile == "")
		outfile = Form("%s.root",infilelist);


	TFile outf(outfile, "RECREATE");

	ifstream inf(infilelist);
	if(!inf.good())
	{ 
		cout << "Input saHist root files are not available!" << endl;
		return false;
	}

	char filename[500];

	vector<TString> saHist_names;
	while(!inf.eof())
	{
		inf.getline(filename, 500);
		if(inf.eof()) continue;

		TFile *rf = TFile::Open(filename,"read");

		TList *list = rf->GetListOfKeys();
		int len = list->GetSize();
		for(int i=0; i<len; i++)
		{
			TKey *key = (TKey *)list->At(i);
			TString type = key->GetClassName();
			if(type.CompareTo("saHist")) continue;

			TString keyname = key->GetName();
			if(!keyname.Contains(vname) && !keyname.Contains("RelLumi")) continue;
			if(keyname.Contains("RUN") || keyname.Contains("FILL")) continue;

			//	  cout << keyname <<", " << type << endl;
			cout<<"NO. "<<saHist_names.size()<<" saHist : "<<keyname.Data()<<endl;
			saHist_names.push_back(keyname);
			saHist *sah = (saHist *) rf->GetObjectChecked(keyname, "saHist");
			if(!sah)
			{
				cout<<"ERROR:: saHist: "<<keyname.Data()<<" Not Found!! \n";
				continue;
			}
			if(!(sah->AutoLoad())){
				cout<<"DEBUG: AutoLoad failed"<<endl;
				continue;
			}
			if(sah)
			{
				outf.cd();
				sah->Write();
				delete sah;
				sah = NULL;
			}
		}

		rf->Close();
	}// while

	outf.Close();

	cout << "Combining files are done!!" << endl;

	return true;
}
