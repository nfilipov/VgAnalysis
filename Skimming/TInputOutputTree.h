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
    Int_t     nVtx;		
    Int_t     run;		
    Long64_t  event;		
    Int_t     lumis;		
    Bool_t    isData;		
    Int_t     nTrksPV;		
    float     vtx;		
    float     vty;		
    float     vtz;		
    float     rho;		
    float     rhoCentral;	
    ULong64_t HLTEleMuX;	
    ULong64_t HLTPho;		
    ULong64_t HLTJet;		
    ULong64_t HLTEleMuXIsPrescaled;
    ULong64_t HLTPhoIsPrescaled;
    ULong64_t HLTJetIsPrescaled;
  };

  
  InputTreeLeaves treeLeaf;
  TBranch  *b_nVtx; //!
  TBranch  *b_run;
  TBranch  *b_event;
  TBranch  *b_lumis;
  TBranch  *b_isData;
  TBranch  *b_nTrksPV;
  TBranch  *b_vtx;
  TBranch  *b_vty;
  TBranch  *b_vtz;
  TBranch  *b_rho;
  TBranch  *b_rhoCentral;
  TBranch  *b_HLTEleMuX;
  TBranch  *b_HLTPho;
  TBranch  *b_HLTJet;
  TBranch  *b_HLTEleMuXIsPrescaled;
  TBranch  *b_HLTPhoIsPrescaled;
  TBranch  *b_HLTJetIsPrescaled;
};


#endif
