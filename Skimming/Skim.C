#include "VgammaSkim.h"

#include "TString.h"
#include "TBenchmark.h"

#include <iostream>

/// data: fileToSkimName = "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_data_ggNtuple_DoubleEG_Run2015D_PromptReco-v4_25ns_JSON_Golden_1560pb_miniAOD.root";
void Skim()
{
  TString fileToSkimName("/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_data_ggNtuple_DoubleEG_Run2015D_PromptReco-v4_25ns_JSON_Golden_1560pb_miniAOD.root");
  TString outDir("./");
  TString nameDir("ggNtuplizer");
  TString nameTree("EventTree");

  TBenchmark time;
  time.Start("time");
  std::cout<<"CPU time = "<<time.GetCpuTime("time")<<", Real time = "<<time.GetRealTime("time")<<std::endl;  

  // //  TString fileToSkimFullName = "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/" + fileToSkimName;

  std::cout<<"file "<<fileToSkimName<<std::endl<<" will be skimmed"<<std::endl;

  //   //  VgammaSkim skimmer(TConfiguration::MUON, TConfiguration::W_GAMMA, TConfiguration::DATA,fileToSkimName);
  VgammaSkim skimmer(fileToSkimName,outDir,nameDir,nameTree);
  skimmer.LoopOverInputTree();
   
  std::cout<<"file "<<fileToSkimName<<std::endl<<" was skimmed"<<std::endl;

  time.Stop("time");
  std::cout<<"CPU time = "<<time.GetCpuTime("time")<<", Real time = "<<time.GetRealTime("time")<<std::endl; 

  std::cout<< "help I'm a rock!"<<std::endl;
  
}
