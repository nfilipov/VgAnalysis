#include "VgammaSkim.h"
#include "datasets.h"
#include "TString.h"
#include "TBenchmark.h"

#include <iostream>

void Skim()
{
  int nS=6; // 0,1,2 :dielectron. 3,4,5:dimuon. 6:DYJets, 7:ZZ. 
  TString fileToSkimName(ggFile[nS]);

  TString outDir("./");
  TString nameDir("ggNtuplizer");
  TString nameTree("EventTree");
  //  bool basic = false;
  bool isMC=false;
  if (nS>5) isMC= true;
  TBenchmark time;
  time.Start("time");
  std::cout<<"CPU time = "<<time.GetCpuTime("time")<<", Real time = "<<time.GetRealTime("time")<<std::endl;  
  // //  TString fileToSkimFullName = "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/" + fileToSkimName;
  std::cout<<"file "<<fileToSkimName<<std::endl<<" will be skimmed"<<std::endl;
  //  VgammaSkim skimmer(TConfiguration::MUON, TConfiguration::W_GAMMA, TConfiguration::DATA,fileToSkimName);
  VgammaSkim skimmer(fileToSkimName,outDir,nameDir,nameTree,isMC,nS);
  skimmer.LoopOverInputTree(isMC);
  std::cout<<"file "<<fileToSkimName<<std::endl<<" was skimmed"<<std::endl;
  
  time.Stop("time");
  std::cout<<"CPU time = "<<time.GetCpuTime("time")<<", Real time = "<<time.GetRealTime("time")<<std::endl; 
  
  std::cout<< "help I'm a rock!"<<std::endl;
  
}
