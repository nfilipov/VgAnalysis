#include "VgammaSkim.h"
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>

VgammaSkim::VgammaSkim(TString inputFileName, TString outDir, TString nameDir, TString nameTree)
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

  TString skimPartOfName="_EleElePhoSkim";

      _skimmedFileName=outDir+skimPartOfName+".root";
      _fileOut = new TFile(_skimmedFileName,"recreate");
      // _fileOut->mkdir(_nameDir);
      // _fileOut->cd(_nameDir);
      _outputTree = new TTree(_nameTree,_nameTree);
      _TREE.InitOutput(_outputTree);
        //method of TInputOutputTree
      _hskim= new TH1F("hskim","hskim",2,0,2);    
      
}

VgammaSkim::~VgammaSkim()
{

}


void VgammaSkim::LoopOverInputTree()
{
  TFile f(_inputFileName,"READ");
  std::cout<<"processing "<<_inputFileName<<std::endl;
  f.cd(_nameDir);
  TTree* tree =(TTree*)gDirectory->Get(_nameTree);
  std::cout<<"Got to the loop."<<std::endl;

  _TREE.InitInput(tree);

  Long64_t nentries = _TREE.fChain->GetEntries();
  // if (_isDebugMode)
  nentries=1000;
  std::cout<<"nentries "<<nentries<<std::endl;


  for (Long64_t entry=0; entry<nentries; entry++) {
    if (entry < 0) break;
    if ((entry%1000000)==0) std::cout<<"entry="<<entry<<std::endl;
    _TREE.GetEntry(entry); // copy of ROOT's classic GetEntry
    _outputTree->Fill();
  }

  _fileOut->cd();
  // _fileOut->cd(_nameDir);
  _outputTree->Write(_nameTree);
  
  //close output files
  _TREE.fChain = 0;
  std::cout<<"file "<<_fileOut->GetName()<<std::endl<<" closed..."<<std::endl;
}
