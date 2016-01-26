#ifndef VgammaSkim_h
#define VgammaSkim_h

#include "TInputOutputTree.h"

#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>


class VgammaSkim{
 public :
  VgammaSkim(TString inputFileName="/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_data_ggNtuple_DoubleEG_Run2015D_PromptReco-v4_25ns_JSON_Golden_1560pb_miniAOD.root",
	     TString outDir="./",
	     TString nameDir="ggNtuplizer",
	     TString nameTree="EventTree",
	     bool isMC = false,
	     bool basic=false);
  virtual ~VgammaSkim();
  void LoopOverInputTree();
  float deltaR(float eta0, float eta1, float phi0, float phi1);
  //the main function which is called from outside
  //const static int numberOfTrees=5;
  //tree [0] is input tree;
  //trees [1]-[4] are output trees

 private :
  TInputOutputTree _TREE;
  /* int _sample; */
  /* bool _isDebugMode; */

  TFile    *_fileOut; //output Files
  TTree    *_outputTree; //output Trees
  TH1F     *_hskim;
  TString  _inputFileName;
  TString  _skimmedFileName;
  TString  _nameDir;
  TString  _nameTree;

  //some leptons...
  const float _eMass = 0.000510998910; // GeV/c2
  const float _muMass = 0.1056583715; // GeV/c2
  const float _tauMass = 1.77682; // GeV/c2
  //photon...for safety
  const float _gMass = 0.; 
  TLorentzVector lepton1;
  TLorentzVector lepton2;
  TLorentzVector ll; 
  TLorentzVector photon;
  TLorentzVector llg;
  bool gWrite ;
  // some variables for the Filling of events:

 float _ePt[100]; 
 float _ePhi[100];
 float _eEta[100];
 vector<float> _llgm;
 vector<float> _llm;
 vector<float> _gammaEt;
 vector<float> _deltaR;
 float _llpt[100];
 float _llphi[100];
 float _lleta[100];
 float _phoEt[100];
 float _phoEta[100];
 float _phoPhi[100];
 float dR1=0;
 float dR2=0;
};

#endif


