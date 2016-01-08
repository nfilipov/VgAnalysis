//////////////////////////////////////////////////////////
// Input and Output tree 
//////////////////////////////////////////////////////////

#ifndef TInputOutputTree_h
#define TInputOutputTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
/* #include "../Include/TEventTree.h" */
class TInputOutputTree{
 public : 
  TInputOutputTree();
  virtual ~TInputOutputTree();
  TString nameFile_;
  TString nameDir_;
  TString nameTree_; 
  void    InitOutput  (TTree* outputTree);
  void    InitInput   (TTree *tree);
 

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  TObjArray*      listOfBranches;
  
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; } 
  struct InputTreeLeaves{
    Int_t           nVtx;
  };

  
  InputTreeLeaves treeLeaf;
  TBranch *b_nVtx; //!
};


#endif
