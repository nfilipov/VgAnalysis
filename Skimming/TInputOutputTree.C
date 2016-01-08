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
  treeLeaf.nVtx=0;		
  treeLeaf.run=0;		
  treeLeaf.event=0;		
  treeLeaf.lumis=0;		
  treeLeaf.isData=0;		
  treeLeaf.nTrksPV=0;		
  treeLeaf.vtx=0;		
  treeLeaf.vty=0;		
  treeLeaf.vtz=0;		
  treeLeaf.rho=0;		
  treeLeaf.rhoCentral=0;	
  treeLeaf.HLTEleMuX=0;	
  treeLeaf.HLTPho=0;		
  treeLeaf.HLTJet=0;		
  treeLeaf.HLTEleMuXIsPrescaled=0;
  treeLeaf.HLTPhoIsPrescaled=0;
  treeLeaf.HLTJetIsPrescaled=0;

  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("nVtx", &treeLeaf.nVtx, &b_nVtx);
  fChain->SetBranchAddress("run", &treeLeaf.run, &b_run);
  fChain->SetBranchAddress("event", &treeLeaf.event, &b_event);
  fChain->SetBranchAddress("lumis", &treeLeaf.lumis, &b_lumis);
  fChain->SetBranchAddress("isData", &treeLeaf.isData, &b_isData);
  fChain->SetBranchAddress("nTrksPV", &treeLeaf.nTrksPV, &b_nTrksPV);
  fChain->SetBranchAddress("vtx", &treeLeaf.vtx, &b_vtx);
  fChain->SetBranchAddress("vty", &treeLeaf.vty, &b_vty);
  fChain->SetBranchAddress("vtz", &treeLeaf.vtz, &b_vtz);
  fChain->SetBranchAddress("rho", &treeLeaf.rho, &b_rho);
  fChain->SetBranchAddress("rhoCentral", &treeLeaf.rhoCentral, &b_rhoCentral);
  fChain->SetBranchAddress("HLTEleMuX", &treeLeaf.HLTEleMuX, &b_HLTEleMuX);
  fChain->SetBranchAddress("HLTPho", &treeLeaf.HLTPho, &b_HLTPho);
  fChain->SetBranchAddress("HLTJet", &treeLeaf.HLTJet, &b_HLTJet);
  fChain->SetBranchAddress("HLTEleMuXIsPrescaled", &treeLeaf.HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
  fChain->SetBranchAddress("HLTPhoIsPrescaled", &treeLeaf.HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
  fChain->SetBranchAddress("HLTJetIsPrescaled", &treeLeaf.HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
}

void TInputOutputTree::InitOutput(TTree* outputTree)
{
  outputTree->Branch("nVtx",                 &treeLeaf.nVtx,       "nVtx/I");
  outputTree->Branch("run",                  &treeLeaf.run,        "run/I");
  outputTree->Branch("event",                &treeLeaf.event,      "event/L");
  outputTree->Branch("lumis",                &treeLeaf.lumis,      "lumis/I");
  outputTree->Branch("isData",               &treeLeaf.isData,     "isData/O");
  outputTree->Branch("nTrksPV",              &treeLeaf.nTrksPV,    "nVtxPV/I");
  outputTree->Branch("vtx",                  &treeLeaf.vtx,        "vtx/F"); 
  outputTree->Branch("vty",                  &treeLeaf.vty,        "vty/F"); 
  outputTree->Branch("vtz",                  &treeLeaf.vtz,        "vtz/F"); 
  outputTree->Branch("rho",                  &treeLeaf.rho,        "rho/F");
  outputTree->Branch("rhoCentral",           &treeLeaf.rhoCentral, "rhoCentral/F");
  outputTree->Branch("HLTEleMuX",            &treeLeaf.HLTEleMuX);
  outputTree->Branch("HLTPho",               &treeLeaf.HLTPho);
  outputTree->Branch("HLTJet",               &treeLeaf.HLTJet);
  outputTree->Branch("HLTEleMuXIsPrescaled", &treeLeaf.HLTEleMuXIsPrescaled);
  outputTree->Branch("HLTPhoIsPrescaled",    &treeLeaf.HLTPhoIsPrescaled);
  outputTree->Branch("HLTJetIsPrescaled",    &treeLeaf.HLTJetIsPrescaled);

}
