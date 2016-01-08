#define TInputOutputTree_cxx
#include "TInputOutputTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TInputOutputTree::TInputOutputTree()
{
   
}

TInputOutputTree::~TInputOutputTree()
{

}

void TInputOutputTree::InitInput(TTree *tree){
  treeLeaf.nVtx =0;
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("nVtx", &treeLeaf.nVtx, &b_nVtx);
}

void TInputOutputTree::InitOutput(TTree* outputTree)
{
  outputTree->Branch("nVtx", &treeLeaf.nVtx, "nVtx/I");
}
