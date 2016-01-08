#ifndef VgammaSkim_h
#define VgammaSkim_h

#include "TInputOutputTree.h"

#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>


class VgammaSkim{
 public :
  VgammaSkim(TString inputFileName="/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_data_ggNtuple_DoubleEG_Run2015D_PromptReco-v4_25ns_JSON_Golden_1560pb_miniAOD.root",
	     TString outDir="./",
	     TString nameDir="ggNtuplizer",
	     TString nameTree="EventTree");
  virtual ~VgammaSkim();
  void LoopOverInputTree();
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

};

#endif


