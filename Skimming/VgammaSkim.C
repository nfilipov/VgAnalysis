#include "VgammaSkim.h"
#include "datasets.h"
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>

VgammaSkim::VgammaSkim(TString inputFileName, TString outDir, TString nameDir, TString nameTree, bool isMC,int config)
{
 
  _nameDir=nameDir;
  _nameTree=nameTree;
  _inputFileName=inputFileName;

  TString inpTemp = _inputFileName;
  while (inpTemp.Contains('/')){
    inpTemp.Remove(0,1);
  }
  inpTemp.ReplaceAll(".root","");

  // TString skimPartOfName[2][2];
  // skimPartOfName[_config.MUON][_config.W_GAMMA]="_MuPhoSkim";// muon + photon
  // skimPartOfName[_config.MUON][_config.Z_GAMMA]="_MuMuPhoSkim";// two muons + photon
  // skimPartOfName[_config.ELECTRON][_config.W_GAMMA]="_ElePhoSkim";// electron + photon
  // skimPartOfName[_config.ELECTRON][_config.Z_GAMMA]="_EleElePhoSkim";// two electrons + photon

  TString skimPartOfName=skimName[config]; //expand this!

  _skimmedFileName=outDir+skimPartOfName+".root";
  _fileOut = new TFile(_skimmedFileName,"recreate");
  // _fileOut->mkdir(_nameDir); // removed this to lose the useless folder
  // _fileOut->cd(_nameDir);   // so that trees can be directly accessed.
  _outputTree = new TTree(_nameTree,_nameTree);
  _TREE.InitOutput(_outputTree,isMC);      //method of TInputOutputTree
  _hStats= new TH1F("hStats","hStats",20,0,20);    
  if (isMC) { std::cout<< " this is simulation " <<std::endl;}
}

VgammaSkim::~VgammaSkim()
{

}


void VgammaSkim::LoopOverInputTree(bool isMC)
{
  TFile f(_inputFileName);
  std::cout<<"processing "<<_inputFileName<<std::endl;
  f.cd(_nameDir);
  TTree* tree =(TTree*)gDirectory->Get(_nameTree);

  _TREE.InitInput(tree,isMC);

  Long64_t nentries = _TREE.fChain->GetEntries();
  std::cout<<"nentries "<<nentries<<std::endl;
  // if (_isDebugMode)
  //  nentries=1e5;

  //  for (Long64_t entry=17023372; entry<19037962; entry++) { //electrons
  // for (Long64_t entry=16108076; entry<17225997; entry++){  //muons in run 258706
  for (Long64_t entry=0; entry<nentries; entry++){  //general setup
    //bool simple=true;
    gWrite=false; // tells us what events to print to skim.

    // //not needed
    // _gammaEt.clear();//init every new custom branch at the beginning of the entry
    // _deltaR.clear();
    // _llgm.clear();
    // _llm.clear();
    //
    if (entry < 0) break;
    if ((entry%50000)==0) std::cout<<"entry="<<entry<<std::endl;
    _TREE.GetEntry(entry); 
    
    // if (_TREE.treeLeaf.run==258706){
      //basic pre-selection corrected. 
      //today, testing the minimum settings...
      if (_TREE.treeLeaf.nPho>0){ //at least one photon 
	if( _TREE.treeLeaf.nEle > 1 ){ //electrons
	  //if( _TREE.treeLeaf.nMu> 1 ){ // muons
	    gWrite =true;
	  }
	  // if( _mSize >1){
	  //   if ( (_TREE.treeLeaf.HLTEleMuX>>20)&1 && _TREE.treeLeaf.muPt->at(0) > 20 ) gWrite =true;
	  // }
	  //}
      
      
      if(gWrite) _outputTree->Fill();
    }
    //alternative skim (simpler)
    //  if(simple) _outputTree->Fill();


  }//end of entry loop
  

  _fileOut->cd();
  // _fileOut->cd(_nameDir);
  _outputTree->Write(_nameTree,TObject::kOverwrite);
  
  //close output files
  _TREE.fChain = 0;
  std::cout<<"file "<<_fileOut->GetName()<<std::endl<<" closed..."<<std::endl;
  f.Close();
  _hStats->SaveAs("hstats"+_nameTree+".root");
}


float VgammaSkim::deltaR(float eta0, float eta1, float phi0, float phi1){
  return sqrt((eta0-eta1)*(eta0-eta1)+(phi0-phi1)*(phi0-phi1));
}

