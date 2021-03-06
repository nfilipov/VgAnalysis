#include <TString.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
using namespace ROOT;
using namespace std;


float deltaR(float eta0, float eta1, float phi0, float phi1){
  return sqrt((eta0-eta1)*(eta0-eta1)+(phi0-phi1)*(phi0-phi1));
}

void ApplySelection()
{

  //  TString _inputFileName("DoubleEG_Run2015D_PR_v4_run258706.root"); // DoubleEG_Run2015D_PR_v4.root or others, like Run2015[C/D]_Oct05
  TString _inputFileName("DYJetsToLL_M-50.root");
  //TString _inputFileName("DoubleEG_Run2015D_PR_v4.root"); // 
  TString _inputFilePrefix("/afs/cern.ch/work/n/nfilipov/private/Vgamma/CMSSW_7_4_14/src/VgAnalysis/Skimming/");
  TFile *f = TFile::Open(_inputFilePrefix+_inputFileName,"READ");
  std::cout<<"processing "<<_inputFileName<<std::endl;
  TTree* tree =(TTree*)f->Get("EventTree");
  

  Long64_t nentries = tree->GetEntries();
  std::cout<<"nentries "<<nentries<<std::endl;
  std::ofstream myfile; //for writing event numbers
  myfile.open (_inputFileName+"_summary.txt");


  TH1F* hPhotons = new TH1F("hPhotons","number of photon candidates",20,0,20);
  TH1F* hPID = new TH1F("hPID","number of photon passing ID cuts",20,0,20);
  TH1F* hPhotonPt = new TH1F("hPhotonPt","p_{T}^{#gamma}",500,0,1000);
  TH1F* hDeltaR = new TH1F("hDeltaR","#DeltaR(l,#gamma)",28,0.7,6.1);
  TH1F* hMll = new TH1F("hMll","M_{ll} [GeV/c^{2}]",500,0,1000);
  TH1F* hMllg = new TH1F("hMllg","M_{ll#gamma} [GeV/c^{2}]",500,0,1000);

  //  TH1F* hElectrons = new TH1F("hElectrons","number of electron candidates",20,0,20); //...
 
  //some leptons...
  const float _eMass = 0.000510998910; // GeV/c2
  const float _muMass = 0.1056583715; // GeV/c2
  const float _tauMass = 1.77682; // GeV/c2
  //and a photon...for safety
  const float _gMass = 0.; 
  TLorentzVector lepton1[1000]={};
  TLorentzVector lepton2[1000]={};
  TLorentzVector leptonPlus[1000]={};
  TLorentzVector leptonMinus[1000]={};
  TLorentzVector photon[1000]={};

  TLorentzVector Zee[1000]={}; // Z candidates and its properties
  std::vector<float> *ZeePt=0;
  std::vector<float> *ZeeRapidity=0;
  std::vector<float> *ZeeMass=0;
  std::vector<float> *ZeePhi=0;
  std::vector<float> *ZeeMt=0;
  TLorentzVector Zeeg[200]={}; // Zg candidates and its properties
  std::vector<float> *ZeegPt=0;
  std::vector<float> *ZeegPhi=0;
  std::vector<float> *ZeegRapidity=0;
  std::vector<float> *ZeegMass=0;
  std::vector<float> *ZeegMt=0;

  //for later
  TLorentzVector Zmm[1000]= {}; // Z candidates and its properties
  std::vector<float> *ZmmPt=0;
  std::vector<float> *ZmmRapidity=0;
  std::vector<float> *ZmmMass=0;
  std::vector<float> *ZmmPhi=0;
  std::vector<float> *ZmmMt=0;
  TLorentzVector Zmmg[200]= {}; // Zg candidates and its properties
  std::vector<float> *ZmmgPt=0;
  std::vector<float> *ZmmgPhi=0;
  std::vector<float> *ZmmgRapidity=0;
  std::vector<float> *ZmmgMass=0;
  std::vector<float> *ZmmgMt=0;

  TLorentzVector ll[1000]={}; // Z candidates
  TLorentzVector llg[1000]={}; // Zg candidates

  float _mllg[100]={};
  float _mll[100]={};
  float _phoEt[100]={};

  float dr1;
  float dr2;
  //  float _eplus,_eminus;
  // TBranch *b_mllg;
  // TBranch *b_mll;
  // TBranch *b_phoEt;
  // TBranch *b_dr;
  nentries = 200;

  gROOT->ProcessLine("#include <vector>");
  //global event variables
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

  tree->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
  tree->SetBranchAddress("run", &run, &b_run);
  tree->SetBranchAddress("event", &event, &b_event);
  tree->SetBranchAddress("lumis", &lumis, &b_lumis);
  tree->SetBranchAddress("isData", &isData, &b_isData);
  tree->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
  tree->SetBranchAddress("vtx", &vtx, &b_vtx);
  tree->SetBranchAddress("vty", &vty, &b_vty);
  tree->SetBranchAddress("vtz", &vtz, &b_vtz);
  tree->SetBranchAddress("rho", &rho, &b_rho);
  tree->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
  tree->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
  tree->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
  tree->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
  tree->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
  tree->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
  tree->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled); 
  //electron variables
  int nEle;
  float _ePt[100]={};
  float _ePhi[100]={};
  float _eEta[100]={};
  float _eSCEta[100]={};
  int nLep=0;
  int nZg=0;
  std::vector<float> *lepPt=0;
  std::vector<int> *lepCh=0;
  std::vector<float> *lepEta=0;
  std::vector<float> *lepSCEta=0;
  std::vector<float> *lepPhi=0;
  std::vector<float> *elePt=0;
  std::vector<float> *elePhi=0;
  std::vector<float> *eleEta=0;
  std::vector<float> *eleCharge=0;
  std::vector<float> *eleSCEta=0;
  std::vector<float> *eleIDMVATrg=0;	     
  std::vector<float> *eledEtaAtVtx=0;	     
  std::vector<float> *eledPhiAtVtx=0;	     
  std::vector<float> *eleDr03TkSumPt=0;	     
  std::vector<float> *elePFClusEcalIso=0;    
  std::vector<float> *elePFClusHcalIso=0;    
  std::vector<float> *eleConvVeto=0;    
  std::vector<float> *eleMissHits=0;    
  std::vector<float> *eleHoverE=0;	     
  std::vector<float> *eleSigmaIEtaIEtaFull5x5=0;
  tree->SetBranchAddress("nEle",&nEle);
  tree->SetBranchAddress("elePt",&elePt);
  tree->SetBranchAddress("eleCharge",&eleCharge);
  tree->SetBranchAddress("elePhi",&elePhi);
  tree->SetBranchAddress("eleEta",&eleEta);		     	     
  tree->SetBranchAddress("eleSCEta",&eleSCEta);		     	     
  tree->SetBranchAddress("eleIDMVATrg"   ,&eleIDMVATrg);	     	     
  tree->SetBranchAddress("eledEtaAtVtx"   ,&eledEtaAtVtx);	     	     
  tree->SetBranchAddress("eledPhiAtVtx"   ,&eledPhiAtVtx);	     	     
  tree->SetBranchAddress("eleDr03TkSumPt" ,&eleDr03TkSumPt);	     	     
  tree->SetBranchAddress("elePFClusEcalIso",&elePFClusEcalIso);         
  tree->SetBranchAddress("elePFClusHcalIso",&elePFClusHcalIso);        
  tree->SetBranchAddress("eleConvVeto",&eleConvVeto);        
  tree->SetBranchAddress("eleMissHits",&eleMissHits);        
  tree->SetBranchAddress("eleHoverE"	    ,&eleHoverE);	     
  tree->SetBranchAddress("eleSigmaIEtaIEtaFull5x5",&eleSigmaIEtaIEtaFull5x5);
  //photon variables
  int nPho;
  float _pPt=0;
  float _pPhi=0;
  float _pEta=0;
  std::vector<float> *PhotonPt=0;
  std::vector<float> *PhotonPhi=0;
  std::vector<float> *PhotonEta=0;
  std::vector<float> *phoEt=0;
  std::vector<float> *phoPhi=0;
  std::vector<float> *phoEta=0;
  std::vector<float> *phoSCEta=0;
  std::vector<float> *phoIDMVA=0;
  std::vector<float> *phoSigmaIEtaIEtaFull5x5=0;
  std::vector<float> *phoHoverE=0;
  std::vector<float> *phoPFPhoIso=0;
  std::vector<float> *phoPFChWorstIso=0;
  std::vector<float> *phoR9=0;
  std::vector<float> *phoR9Full5x5=0;
  std::vector<int>   *phohasPixelSeed=0;
  std::vector<int>   *phoEleVeto=0;
  tree->SetBranchAddress("nPho",&nPho);
  tree->SetBranchAddress("phoEt",&phoEt);
  tree->SetBranchAddress("phohasPixelSeed",&phohasPixelSeed);
  tree->SetBranchAddress("phoEleVeto",&phoEleVeto);
  tree->SetBranchAddress("phoPhi",&phoPhi);
  tree->SetBranchAddress("phoEta",&phoEta);
  tree->SetBranchAddress("phoSCEta",&phoSCEta);
  tree->SetBranchAddress("phoIDMVA",&phoIDMVA);	
  tree->SetBranchAddress("phoSigmaIEtaIEtaFull5x5",&phoSigmaIEtaIEtaFull5x5);
  tree->SetBranchAddress("phoHoverE",&phoHoverE);
  tree->SetBranchAddress("phoPFPhoIso",&phoPFPhoIso);
  tree->SetBranchAddress("phoPFChWorstIso",&phoPFChWorstIso);
  tree->SetBranchAddress("phoR9",&phoR9);
  //  tree->SetBranchAddress("phoR9Full5x5",&phoR9Full5x5);
  //moun variables
  //muon
  Int_t          nMu;
  std::vector<float>  *muPt=0;
  std::vector<float>  *muEn=0;
  std::vector<float>  *muEta=0;
  std::vector<float>  *muPhi=0;
  std::vector<int>    *muCharge=0;
  std::vector<int>    *muType=0;
  std::vector<Bool_t> *muIsLooseID=0;
  std::vector<Bool_t> *muIsMediumID=0;
  std::vector<Bool_t> *muIsTightID=0;
  std::vector<Bool_t> *muIsSoftID=0;
  std::vector<Bool_t> *muIsHighPtID=0;
  //global muon? PF muon ?
  std::vector<float>  *muD0=0;
  std::vector<float>  *muDz=0;
  std::vector<float>  *muChi2NDF=0;
  std::vector<int>    *muTrkLayers=0;
  std::vector<int>    *muPixelHits=0;
  std::vector<int>    *muStations=0;
  std::vector<float>  *muPFChIso=0;
  std::vector<float>  *muPFPhoIso=0;
  std::vector<float>  *muPFNeuIso=0;
  std::vector<float>  *muPFPUIso=0;
  std::vector<int>    *muFiredTrgs=0;
  //muon channel
  tree->SetBranchAddress("nMu",&nMu);
  tree->SetBranchAddress("muPt",&muPt);
  tree->SetBranchAddress("muEn",&muEn);
  tree->SetBranchAddress("muEta",&muEta);
  tree->SetBranchAddress("muPhi",&muPhi);
  tree->SetBranchAddress("muCharge",&muCharge); 
  tree->SetBranchAddress("muType",&muType);
  tree->SetBranchAddress("muIsLooseID",&muIsLooseID);
  tree->SetBranchAddress("muIsMediumID",&muIsMediumID);
  tree->SetBranchAddress("muIsTightID",&muIsTightID);
  tree->SetBranchAddress("muIsSoftID",&muIsSoftID);
  tree->SetBranchAddress("muIsHighPtID",&muIsHighPtID);
  tree->SetBranchAddress("muD0",&muD0	   );
  tree->SetBranchAddress("muDz",&muDz	   );
  tree->SetBranchAddress("muChi2NDF",&muChi2NDF); 
  tree->SetBranchAddress("muTrkLayers",&muTrkLayers);
  tree->SetBranchAddress("muPixelHits",&muPixelHits);
  tree->SetBranchAddress("muStations",&muStations);
  tree->SetBranchAddress("muPFChIso",&muPFChIso );
  tree->SetBranchAddress("muPFPhoIso",&muPFPhoIso);
  tree->SetBranchAddress("muPFNeuIso",&muPFNeuIso);
  tree->SetBranchAddress("muPFPUIso",&muPFPUIso );
  tree->SetBranchAddress("muFiredTrgs",&muFiredTrgs);

  //make reduced tree (only passing the cuts)
  TFile *zgf = new TFile("red_eeg_AllCuts_"+_inputFileName,"RECREATE");
  TTree *_zgt = new TTree("zgtree","RECREATE");
  // _zgt->Branch("nVtx",                 &nVtx,       "nVtx/I");
  _zgt->Branch("run",                  &run,        "run/I");
  _zgt->Branch("event",                &event,      "event/L");
  _zgt->Branch("phoEt",                &phoEt	 );  
  _zgt->Branch("lumis",                &lumis,      "lumis/I");
  _zgt->Branch("isData",               &isData,     "isData/O");
  _zgt->Branch("nTrksPV",              &nTrksPV,    "nVtxPV/I");
  _zgt->Branch("vtx",                  &vtx,        "vtx/F"); 
   _zgt->Branch("vty",                  &vty,        "vty/F"); 
   _zgt->Branch("vtz",                  &vtz,        "vtz/F"); 
   _zgt->Branch("rho",                  &rho,        "rho/F");
  _zgt->Branch("rhoCentral",           &rhoCentral, "rhoCentral/F");
  _zgt->Branch("HLTEleMuX",            &HLTEleMuX, "HLTEleMuX/l");
  _zgt->Branch("HLTPho",               &HLTPho, "HLTPho/l");
  _zgt->Branch("HLTJet",               &HLTJet, "HLTJet/l");
  _zgt->Branch("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, "HLTEleMuXIsPrescaled/l");
  _zgt->Branch("HLTPhoIsPrescaled",    &HLTPhoIsPrescaled, "HLTPhoIsPrescaled/l");
  _zgt->Branch("HLTJetIsPrescaled",    &HLTJetIsPrescaled, "HLTJetIsPrescaled/l");
  //electrons
  _zgt->Branch("nEle",&nEle);
  _zgt->Branch("elePt",&elePt);
  _zgt->Branch("elePhi",&elePhi);
  _zgt->Branch("eleCharge",&eleCharge);
  _zgt->Branch("eleEta",&eleEta);		     	     
  _zgt->Branch("eleSCEta",&eleSCEta);		     	     
  _zgt->Branch("eleIDMVATrg"   ,&eleIDMVATrg);	     	     
  _zgt->Branch("eledEtaAtVtx"   ,&eledEtaAtVtx);	     	     
  _zgt->Branch("eledPhiAtVtx"   ,&eledPhiAtVtx);	     	     
  _zgt->Branch("eleDr03TkSumPt" ,&eleDr03TkSumPt);	     	     
  _zgt->Branch("elePFClusEcalIso",&elePFClusEcalIso);         
  _zgt->Branch("elePFClusHcalIso",&elePFClusHcalIso);        
  _zgt->Branch("eleConvVeto",&eleConvVeto);        
  _zgt->Branch("eleMissHits",&eleMissHits);        
  _zgt->Branch("eleHoverE"	    ,&eleHoverE);	     
  _zgt->Branch("eleSigmaIEtaIEtaFull5x5",&eleSigmaIEtaIEtaFull5x5);
  //photons
  _zgt->Branch("nPho",&nPho);
  _zgt->Branch("phoEt",&phoEt);
  //  _zgt->Branch("phoEt",&phoCalibEt)
  _zgt->Branch("phoPhi",&phoPhi);
  _zgt->Branch("phoEta",&phoEta);
  _zgt->Branch("phohasPixelSeed",&phohasPixelSeed);
  _zgt->Branch("phoEleVeto",&phoEleVeto);
  _zgt->Branch("phoSCEta",&phoSCEta);
  _zgt->Branch("phoIDMVA",&phoIDMVA);	
  _zgt->Branch("phoSigmaIEtaIEtaFull5x5",&phoSigmaIEtaIEtaFull5x5);
  _zgt->Branch("phoHoverE",&phoHoverE);
  _zgt->Branch("phoPFPhoIso",&phoPFPhoIso);
  _zgt->Branch("phoPFChWorstIso",&phoPFChWorstIso);
  _zgt->Branch("phoR9",&phoR9);
  //  _zgt->Branch("phoR9Full5x5",&phoR9Full5x5); /// test, not used further.
  //muons
  // //muon channel
  // _zgt->Branch("nMu",&nMu);
  // _zgt->Branch("muPt",&muPt);
  // _zgt->Branch("muEn",&muEn);
  // _zgt->Branch("muEta",&muEta);
  // _zgt->Branch("muPhi",&muPhi);
  // _zgt->Branch("muCharge",&muCharge); 
  // _zgt->Branch("muType",&muType);
  // _zgt->Branch("muIsLooseID",&muIsLooseID);
  // _zgt->Branch("muIsMediumID",&muIsMediumID);
  // _zgt->Branch("muIsTightID",&muIsTightID);
  // _zgt->Branch("muIsSoftID",&muIsSoftID);
  // _zgt->Branch("muIsHighPtID",&muIsHighPtID);
  // _zgt->Branch("muD0",&muD0	   );
  // _zgt->Branch("muDz",&muDz	   );
  // _zgt->Branch("muChi2NDF",&muChi2NDF); 
  // _zgt->Branch("muTrkLayers",&muTrkLayers);
  // _zgt->Branch("muPixelHits",&muPixelHits);
  // _zgt->Branch("muStations",&muStations);
  // _zgt->Branch("muPFChIso",&muPFChIso );
  // _zgt->Branch("muPFPhoIso",&muPFPhoIso);
  // _zgt->Branch("muPFNeuIso",&muPFNeuIso);
  // _zgt->Branch("muPFPUIso",&muPFPUIso );
  // _zgt->Branch("muFiredTrgs",&muFiredTrgs);

  //my development stuff
  _zgt->Branch("ZeePt",&ZeePt);
  _zgt->Branch("ZeeRapidity",&ZeeRapidity);
  _zgt->Branch("ZeePhi",&ZeePhi);
  _zgt->Branch("ZeeMass",&ZeeMass);
  _zgt->Branch("ZeeMt",&ZeeMt);

  _zgt->Branch("nZg",&nLep,"nZg/I");
  _zgt->Branch("ZeegPt",&ZeegPt);
  _zgt->Branch("ZeegRapidity",&ZeegRapidity);
  _zgt->Branch("ZeegPhi",&ZeegPhi);
  _zgt->Branch("ZeegMass",&ZeegMass);
  _zgt->Branch("ZeegMt",&ZeegMt);

  _zgt->Branch("PhotonPt",&PhotonPt);
  _zgt->Branch("PhotonEta",&PhotonEta);
  _zgt->Branch("PhotonPhi",&PhotonPhi);
 
  _zgt->Branch("nLep",&nLep,"nLep/I");
  _zgt->Branch("lepPt",&lepPt);
  _zgt->Branch("lepCh",&lepCh);
  _zgt->Branch("lepEta",&lepEta);
  _zgt->Branch("lepSCEta",&lepSCEta);
  _zgt->Branch("lepPhi",&lepPhi);

  //entry loop
  for (int entry=0; entry<nentries;entry++)
    {
      nZg=0;
      nLep=0;
      bool Zgamma_in=false;
      ZeegMass->clear();
      ZeegPt->clear();
      ZeegRapidity->clear();
      ZeegPhi->clear();
      ZeegMt->clear();
      ZeePt->clear();
      ZeeRapidity->clear();
      ZeePhi->clear();
      ZeeMass->clear();
      ZeeMt->clear();				   
      PhotonPt->clear();
      PhotonEta->clear();
      PhotonPhi->clear();
      // int nep=0;//number of plus e
      // int nem=0;//number of minus e
      tree->GetEntry(entry);
      if (nEle<2  ){ continue; } else {

	int ig=0; //number of Zgamma candidates
   
	//START Zgamma candidate
	//start photon loop
	for (int ip =0 ; ip<nPho;ip++)
	  {
	    _pPt=phoEt->at(ip); 
	    _pPhi=phoPhi->at(ip); 
	    _pEta=phoSCEta->at(ip);
	    if( 
	       _pPt > 10  && ((HLTEleMuX>>7)&1)==1
	       &&  (phoEleVeto->at(ip)==true 
		    )
	       && ((abs(phoSCEta->at(ip))<1.444 || abs(phoSCEta->at(ip))>1.566) && abs(phoSCEta->at(ip))<2.5)
		)
	      {
		if((abs(phoSCEta->at(ip))<1.444    &&  /// barrel  
		    phoPFChWorstIso->at(ip)<15 && /// pho PF charged hadron worst iso
		    phoPFPhoIso->at(ip)<15 && /// PF Photon ECal iso
		    phoHoverE->at(ip)<0.08 &&// photon H / E
		    phoSigmaIEtaIEtaFull5x5->at(ip)<0.012 &&  /// sig_ietaieta
		    phoIDMVA->at(ip) > 0.4  /// photon  MVA id
		    )||
		   (abs(phoSCEta->at(ip))>1.566  && /// endcap  
		    phoR9->at(ip) > 0.85    &&       //// removed for now
		    phoPFChWorstIso->at(ip)<15 &&  /// pho PF charged hadron worst iso
		    phoPFPhoIso->at(ip)<15 &&      /// PF Photon ECal iso
		    phoHoverE->at(ip)<0.05 &&      // photon H / E
		    phoSigmaIEtaIEtaFull5x5->at(ip)<0.027 && /// sig_ietaieta
		    phoIDMVA->at(ip) > 0.3    
		    /// photon  MVA and id
		    )){

		  //START ELECTRON LOOP
		  for (int ie=0; ie<nEle;ie++)
		    {
		      _ePt[ie] = elePt->at(ie);
		      _eEta[ie] = eleEta->at(ie);
		      _eSCEta[ie] = eleSCEta->at(ie);
		      _ePhi[ie] = elePhi->at(ie);
		      if (((  abs(eleSCEta->at(ie))<0.8 && ( eleIDMVATrg->at(ie)>0.972153 ))
			   ||((abs(eleSCEta->at(ie))>0.8 && abs(eleSCEta->at(ie))<1.479) && eleIDMVATrg->at(ie)>0.922126)
			   ||( (abs(eleSCEta->at(ie))>1.479 && abs(eleSCEta->at(ie))<2.5) && eleIDMVATrg->at(ie)>0.610764 
			       )
			   )// && (eleConvVeto->at(ie)==true)
			  )
			{
			  std::cout<<"run="<<run<<", lumi= "<<lumis<<" , event="<<event<<std::endl;	   
			  if((abs(eleSCEta->at(ie))<1.479   &&     /// barrel
			      abs(eledEtaAtVtx->at(ie))<0.0095 &&     /// dEtaIN
			      abs(eledPhiAtVtx->at(ie))<0.065 &&     /// dPhiIn
			      eleDr03TkSumPt->at(ie)/_ePt[ie]<0.18  &&    /// Trk iso
			      //elePFClusEcalIso->at(ie)/_ePt[ie]<0.37   &&   /// ecal iso
			      //elePFClusHcalIso->at(ie)/_ePt[ie]<0.25  &&   /// hcal iso
			      eleHoverE->at(ie)<0.09 &&           /// H/E
			      eleSigmaIEtaIEtaFull5x5->at(ie)<0.012
			      ) /// sig_ietaieta
			     ||
			     (abs(eleSCEta->at(ie))>1.479  &&    ////endcaps
			      eleDr03TkSumPt->at(ie)/_ePt[ie]<0.18 &&    /// trk iso
			      //elePFClusEcalIso->at(ie)/_ePt[ie]<0.45  &&  ///ecal iso
			      //elePFClusHcalIso->at(ie)/_ePt[ie]<0.28 &&   /// Hcal iso.
			      eleHoverE->at(ie)<0.09 &&          /// H/E
			      eleSigmaIEtaIEtaFull5x5->at(ie)<0.033
			      )) ////Sig_ietaieta
			    {
			      dr1=deltaR(_eSCEta[ie],phoSCEta->at(ip),_ePhi[ie],phoPhi->at(ip));
			      if (dr1>0.7)
				{
				  //start second electron loop
				  //START Second ELECTRON LOOP
				  for (int ie2=ie+1; ie2<nEle;ie2++)
				    {
				      _ePt[ie2] = elePt->at(ie2);
				      _eEta[ie2] = eleEta->at(ie2);
				      _eSCEta[ie2] = eleSCEta->at(ie2);
				      _ePhi[ie2] = elePhi->at(ie2);
				      if (((  abs(eleSCEta->at(ie2))<0.8 && ( eleIDMVATrg->at(ie2)>0.972153 ))
					  ||((abs(eleSCEta->at(ie2))>0.8 && abs(eleSCEta->at(ie2))<1.479) && eleIDMVATrg->at(ie2)>0.922126)
					  ||( (abs(eleSCEta->at(ie2))>1.479 && abs(eleSCEta->at(ie2))<2.5) && eleIDMVATrg->at(ie2)>0.610764 
					      )//risky business
					   ) //&& (eleConvVeto->at(ie2)==true)
					  )
					{	   
					  if((abs(eleSCEta->at(ie2))<1.479   &&     /// barrel
					      abs(eledEtaAtVtx->at(ie2))<0.0095 &&     /// dEtaIN
					      abs(eledPhiAtVtx->at(ie2))<0.065 &&     /// dPhiIn
					      eleDr03TkSumPt->at(ie2)/_ePt[ie2]<0.18  &&    /// Trk iso
					      //elePFClusEcalIso->at(ie2)/_ePt[ie2]<0.37   &&   /// ecal iso
					      //elePFClusHcalIso->at(ie2)/_ePt[ie2]<0.25  &&   /// hcal iso
					      eleHoverE->at(ie2)<0.09 &&           /// H/E
					      eleSigmaIEtaIEtaFull5x5->at(ie2)<0.012
					      ) /// sig_ietaieta
					     ||
					     (abs(eleSCEta->at(ie2))>1.479  &&    ////endcaps
					      eleDr03TkSumPt->at(ie2)/_ePt[ie2]<0.18 &&    /// trk iso
					      //elePFClusEcalIso->at(ie2)/_ePt[ie2]<0.45  &&  ///ecal iso
					      //elePFClusHcalIso->at(ie2)/_ePt[ie2]<0.28 &&   /// Hcal iso.
					      eleHoverE->at(ie2)<0.09 &&          /// H/E
					      eleSigmaIEtaIEtaFull5x5->at(ie2)<0.033
					      )) ////Sig_ietaieta
					    {
					      ///second electron passes ID cuts. compute DR2
					      dr2=deltaR(_eSCEta[ie2],phoSCEta->at(ip),_ePhi[ie2],phoPhi->at(ip));
					      std::cout<<"here"<<std::endl;	   
					      if (dr2>0.7)
						{
						  if ((_ePt[ie] > 20 && _ePt[ie2] > 15)||( _ePt[ie] > 15 && _ePt[ie2] > 20))
						    {
						      //fill both lepton's info after making sure we have a working pair.
						      lepCh->push_back(eleCharge->at(ie));
						      lepPt->push_back(_ePt[ie]);
						      lepEta->push_back(_eEta[ie]);
						      lepSCEta->push_back(_eSCEta[ie]);
						      lepPhi->push_back(_ePhi[ie]);
						      lepton1[ie].SetPtEtaPhiM(_ePt[ie],_eSCEta[ie],_ePhi[ie],_eMass);
						      //
						      lepCh->push_back(eleCharge->at(ie2));
						      lepPt->push_back(_ePt[ie2]);
						      lepEta->push_back(_eEta[ie2]);
						      lepSCEta->push_back(_eSCEta[ie2]);
						      lepPhi->push_back(_ePhi[ie2]);
						      lepton2[ie2].SetPtEtaPhiM(_ePt[ie2],_eSCEta[ie2],_ePhi[ie2],_eMass);
						      ///can build the tmp Z. then Zg
						      TLorentzVector tmpZ;
						      tmpZ = lepton1[ie]+lepton2[ie2];
						      if(tmpZ.M() > 50) // check the Z mass requirement
							{
							  nLep=nLep+2;
							  nZg++;
							  photon[ip].SetPtEtaPhiM(_pPt,_pEta,_pPhi,_gMass);
							  Zeeg[ig]=photon[ip]+tmpZ;
							  ZeegMass->push_back(Zeeg[ig].M());
							  ZeegPt->push_back(Zeeg[ig].Pt());
							  ZeegRapidity->push_back(Zeeg[ig].Rapidity());
							  ZeegPhi->push_back(Zeeg[ig].Phi());
							  ZeegMt->push_back(Zeeg[ig].Mt());
							  ZeePt->push_back(tmpZ.Pt());
							  ZeeRapidity->push_back(tmpZ.Rapidity());
							  ZeePhi->push_back(tmpZ.Phi());
							  ZeeMass->push_back(tmpZ.M());
							  ZeeMt->push_back(tmpZ.Mt());
							  std::cout<<"run="<<run<<", lumi= "<<lumis<<", pt("<<ie<<") = "<<lepton1[ie].Pt()<<", pt("<<ie2<<") = "<<lepton2[ie2].Pt()<<",  Mll="<<tmpZ.M()<<" GeV/c2, phoEt("<<ip<<")="<<_pPt<<" GeV/c, and "<<ig<<"th Zg candidate, mass="<<Zeeg[ig].M()<<" GeV/c2"<<std::endl;
							  PhotonPt->push_back(_pPt);
							  PhotonEta->push_back(_pEta);
							  PhotonPhi->push_back(_pPhi);
							  ig++;
							  Zgamma_in=true;
							  // Zgamma_in = true;
							  myfile<<event<<" "<<lumis<<" "<<run<<endl;
							}
						    }//lepton pt cut 
						}//dr2 cut
					    }//electron2 ID cuts
					}//electron2 MVA cuts
				    }///electron 2 loop end
				}//dr1 cut
			    }//electron1 ID cuts
			}//electron1 MVA cuts
		    }//electron 1 loop end
		}//Photon MVA cut
	      }//photon ID cut ?
	  }//photon loop
      }///if nEle>2 cut

      if (Zgamma_in){	  _zgt->Fill();}
    }//end of entry loop

  //  myfile.close();
  // hPhotonPt->SaveAs("hist_PhotonPt_"+_inputFileName);
  // hDeltaR->SaveAs("hist_DeltaR_"+_inputFileName);
  // hMll->SaveAs("hist_Mll_"+_inputFileName);
  // hMllg->SaveAs("hist_Mllg_"+_inputFileName);
  // hPhotons->SaveAs("test_cuts.root");

  _zgt->GetDirectory()->cd();
  TH1F *h2 = (TH1F*)f->Get("hPUTrue");
  h2->Write();
  zgf->Write();
}



