#ifndef JPSIAN
#define JPSIAN
	
#include <set>
#include <map>
#include <vector>
	
class TFile;
class TTree;
class TH1D;
class TCanvas;
class TH2D;	
	
class JpsiAN
{
public:
  JpsiAN();
  ~JpsiAN();

  enum SPIN {UP=1, DN=-1};
  enum SIDE {Left, Right};
  enum DET  {North, South};
  enum BEAM {Blue, Yelo};
  enum PAR {P0, P1, P2, P3, Njpsi, Mjpsi, Sjpsi, Npsi, Mpsi, Spsi};

private:
  double  _Par[2][8];
  float  Evt_Nmu, Run_Number, Evt_Number;
  int    Fill_Number;
  float  charge;
  float  Evt_bbcZ, Evt_vtxchi2, Evt_vtxooz;
  float Tr0_DDG0, Tr0_DG0, Tr0_DS3, Tr0_DS3ctp;
  float  Tr0_chi2, Tr0_idchi2;
  float Tr0_idhits,  Tr0_trhits;
  float  Tr0_px, Tr0_py, Tr0_pz;
  float Tr1_DDG0, Tr1_DG0, Tr1_DS3, Tr1_DS3ctp;
  float  Tr1_chi2, Tr1_idchi2;
  float Tr1_idhits, Tr1_trhits;
  float  Tr1_px, Tr1_py, Tr1_pz;
  float mass, p, pT, rapidity, x1, x2, xF;
  unsigned int  SpinX_ID, Pol_Y, Pol_B, GL1X_ID;
  float  X0, Y0, R0;
  float Tr0_nidhits,Tr1_nidhits;
  float Tr0_ntrhits,Tr1_ntrhits;
  short Tr0_lastgap, Tr1_lastgap;
  float Tr0_dca_r, Tr1_dca_r;
  float px,py,pz;
  unsigned int beamclk0;
  bool same_event;
  UInt_t lvl1_trigscaled;

  int    _nRuns;         // number of good runs;
  std::set<int> _runList;
  std::map<int, int> _fillTable;
  std::map<int, int> _RunTable;

  TFile *_file;
  TTree *_t;

  TH1D  *_hMass[2];
  TH1D  *_hXF[2];
  TH1D  *_hPt[2][2];	//[arm][pt_bin]
  TH2D  *_hIcl_u_u_B[2];
  TH2D  *_hIcl_u_d_B[2];
  TH2D  *_hIcl_d_u_B[2];
  TH2D  *_hIcl_d_d_B[2];
  TH2D  *_hBgr_u_u_B[2];
  TH2D  *_hBgr_u_d_B[2];
  TH2D  *_hBgr_d_u_B[2];
  TH2D  *_hBgr_d_d_B[2];
  TH2D  *_hit_plots[2];


  double  _MinJpsiMs[2];     // mean - 2*sigma
  double  _MaxJpsiMs[2];     // mean + 2*sigma

  double  _SB[2];
  double  _xF[2];
  double  _Pt[2];

  int     _nFills;
  double *_Fill;             // this is used as x-axis for graphs
  double *_eFill;            // used to define graphs

	
  // relative luminosity
  double *_R1[2];         // [beam][fill]
  double *_R2[2];         // [beam][fill]
  double *_R3[2];         // [beam][fill]
  double *_eR1[2];        // [beam][fill]
  double *_eR2[2];        // [beam][fill]
  double *_eR3[2];        // [beam][fill]

  double _L_uuB, _L_udB, _L_duB, _L_ddB,_L_uuY, _L_udY, _L_duY, _L_ddY;	
  double R1B,R2B,R3B,R1Y,R2Y,R3Y;
  double eR1B,eR2B,eR3B,eR1Y,eR2Y,eR3Y;
  double _PB, _ePB, _PY, _ePY;  

  // spin information
  std::map<int, float>  _PolB;
  std::map<int, float>  _ePolB;
  std::map<int, float>  _PolY;
  std::map<int, float>  _ePolY;
  std::map<int, short*> _spinB;
  std::map<int, short*> _spinY;
  std::map<int, double*> _BBCinB;
  std::map<int, double*> _BBCinY;
	
  std::vector<TCanvas*> _vCanvas;
  TCanvas *c0;
  TCanvas *c_yeild[2];
  TCanvas *cPt;
	
  void ReadGoodRunList(); 
  void ReadFillTable();
  void GetSpinInfo();
  void SetTree();
  void SetHistos();
  void SetArrays();

  bool Cut();
  void RelLumi();
  void GetYeilds();
  void GetAN();
};
#endif
