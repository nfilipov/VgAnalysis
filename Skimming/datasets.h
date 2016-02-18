#ifndef datasets_h
#define datasets_h


using namespace std;

std::string ggFile[8]={"/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_DoubleEG_Run2015C_Oct05_miniAOD.root", //Run2015C DoubleEG Oct05
			    "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_DoubleEG_Run2015D_Oct05_miniAOD.root", //Run2015D DoubleEG Oct05
			    "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_DoubleEG_Run2015D_PR_v4_miniAOD.root", //Run2015D DoubleEG PR_v4
			    "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_DoubleMu_Run2015C_Oct05_miniAOD.root",//Run2015C DoubleMu Oct05
			    "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_DoubleMu_Run2015D_Oct05_miniAOD.root",//Run2015D DoubleMu Oct05
		       "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_DoubleMu_Run2015D_PR_v4_miniAOD.root", //Run2015D DoubleMu PR_v4; // mc questionable
		       "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/mc/job_spring15_ggNtuple_DYJetsToLL_M-50_amcatnlo_pythia8_25ns_miniAOD.root", //DYJetsToLL, 25ns
		       "/afs/cern.ch/user/n/nfilipov/eos/cms/store/group/phys_smp/ggNtuples/13TeV/mc/job_spring15_ggNtuple_ZZTo4L_powheg_pythia8_25ns_miniAOD.root"}; //ZZ pythia8, 25ns.

//TString skimPartOfName="mc_DYJetsToLL_M10-50_LLGSkim_big"; //expand this!
TString skimName[8]={"DoubleEG_Run2015C_Oct05","DoubleEG_Run2015D_Oct05","DoubleEG_Run2015D_PR_v4","DoubleMu_Run2015C_Oct05","DoubleMu_Run2015D_Oct05","DoubleMu_Run2015D_PR_v4","DYJetsToLL_M-50","ZZTo4L"};
#endif
