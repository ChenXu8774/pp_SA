/*!
 * \file merge_saHists_hist_by_hist.C
 * \brief merge files with saHists hist-by-hist, modified from Jin's merge_saHist.C, slower than it, comsumes less memory.
 * \author Haiwang Yu <yuhw.pku@gmail.com>
 * \version $Revision: 1.3 $
 * \date $Date: 2016/04/13 00:01:47 $
 */


#include <iostream>
#include <cassert>
#include <stdexcept>

using namespace std;

/*!
 * \param 
 */
bool merge_saHists_hist_by_hist(
		int i_saHist = -1, //condor ID
		char *infilelist = "test.lst",//input
		char *outfile = ""//output
		)
{
	char *vname = "_m_";
	cout<<"i_saHist: "<<i_saHist<<endl;
	gSystem->Load("libspin_analyzer");

	if(outfile == "")
		if(i_saHist >= 0)
			outfile = Form("%s_%d.root",infilelist,i_saHist);
		else
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

		TFile *rf = new TFile(filename);

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
		}
		break;
	}// while

	saHist *sah_combined = NULL;
	for(int iname=0;iname<saHist_names.size();iname++)
	{
		if(i_saHist >= 0)
		{
			if(iname != i_saHist) continue;//continue saHist names loop
		}
		TString keyname = saHist_names[iname];
		cout<<"Processing: NO."<<iname<<" saHist : "<<keyname.Data()<<endl;

		int count = 0;
		bool first = true;
		inf.clear() ;
		inf.seekg(0, ios::beg) ;
		while(!inf.eof())
		{
			inf.getline(filename, 500);
			if(inf.eof()) continue;

			count++;
			cout << filename << " is loaded : \t" << count << endl;     

			TFile *rf = new TFile(filename);

			saHist *sah = (saHist *) rf->GetObjectChecked(keyname, "saHist");
			assert(sah);	  
			sah->Verbosity(0, true);
			if(!(sah->AutoLoad())){
				cout<<"DEBUG: AutoLoad failed"<<endl;
				continue;
			}
			if(first) 
			{	      
				sah_combined = sah;
			}
			else
			{
				if(!sah_combined)
				{
					cout<<"ERROR: sah_combined not initialized!"<<endl;
					break; //break while
				}

				sah_combined->MergeHistCount(sah); // merge master saHist
				sah_combined->AdoptFillSubHist(sah); // merge fill sub saHist
				//sah_combined->AdoptRunSubHist(sah); // merge run sub saHist	      

				delete sah;
			}

			rf->Close();
			delete rf;

			if(first) first = false;

		}// Loop all files
		if(sah_combined)
		{
			outf.cd();
			sah_combined->Write();
			delete sah_combined;
			sah_combined = NULL;
		}
	}// Loop all saHist names


	outf.Close();

	cout << "Merging files are done!!" << endl;

	return true;
}
