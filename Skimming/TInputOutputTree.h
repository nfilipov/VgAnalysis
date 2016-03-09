//////////////////////////////////////////////////////////
// Input and Output tree 
//////////////////////////////////////////////////////////

#ifndef TInputOutputTree_h
#define TInputOutputTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
/* #include "../Include/TEventTree.h" */
class TInputOutputTree{
 public : 
  TInputOutputTree();
  virtual ~TInputOutputTree();
  TString nameFile_;
  TString nameDir_;
  TString nameTree_; 
  bool isMC=true;
  void    InitOutput  (TTree* outputTree, bool isMC);
  void    InitInput   (TTree *tree, bool isMC);

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  TObjArray*      listOfBranches;
  
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; } 
  struct InputTreeLeaves{
    //global event
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
    //gen info (if isMC==true)
    // (local) variables associated with tree branches
    vector<float>    *pdf;
    float            pthat;
    float            processID;
    float            genWeight;
    float            genHT;

    Int_t            nPUInfo;
    vector<int>      *nPU;
    vector<int>      *puBX;
    vector<float>    *puTrue;

    Int_t            nMC;          
    vector<int>      *mcPID;
    vector<float>    *mcVtx;
    vector<float>    *mcVty;
    vector<float>    *mcVtz;
    vector<float>    *mcPt;
    vector<float>    *mcMass;
    vector<float>    *mcEta;
    vector<float>    *mcPhi;
    vector<float>    *mcE;
    vector<float>    *mcEt;
    vector<int>      *mcGMomPID;
    vector<int>      *mcMomPID;
    vector<float>    *mcMomPt;
    vector<float>    *mcMomMass;
    vector<float>    *mcMomEta;
    vector<float>    *mcMomPhi;
    vector<int>      *mcIndex;
    vector<UShort_t> *mcStatusFlag;
    vector<int>      *mcParentage;
    vector<int>      *mcStatus;
    vector<float>    *mcCalIsoDR03;
    vector<float>    *mcTrkIsoDR03;
    vector<float>    *mcCalIsoDR04;
    vector<float>    *mcTrkIsoDR04;
    //muon
    Int_t          nMu;
    vector<float>  *muPt;
    vector<float>  *muEn;
    vector<float>  *muEta;
    vector<float>  *muPhi;
    vector<int>    *muCharge;
    vector<int>    *muType;
    vector<Bool_t> *muIsLooseID;
    vector<Bool_t> *muIsMediumID;
    vector<Bool_t> *muIsTightID;
    vector<Bool_t> *muIsSoftID;
    vector<Bool_t> *muIsHighPtID;
    //global muon? PF muon ?
    vector<float>  *muD0;
    vector<float>  *muDz;
    vector<float>  *muChi2NDF;
    /* vector<float>  *muInnerD0; */
    /* vector<float>  *muInnerDz; */
    vector<int>    *muTrkLayers;
    /* vector<int>    *muPixelLayers; */
    vector<int>    *muPixelHits;
    /* vector<int>    *muMuonHits; */
    vector<int>    *muStations;
    /* vector<int>    *muMatches; */
    /* vector<int>    *muTrkQuality; */
    /* vector<float>  *muIsoTrk; */
    vector<float>  *muPFChIso;
    vector<float>  *muPFPhoIso;
    vector<float>  *muPFNeuIso;
    vector<float>  *muPFPUIso;
    vector<int>    *muFiredTrgs;
    /* vector<float>  *muInnervalidFraction; */
    /* vector<float>  *musegmentCompatibility; */
    /* vector<float>  *muchi2LocalPosition; */
    /* vector<float>  *mutrkKink; */
    /* vector<float>  *muBestTrkPtError; */
    /* vector<float>  *muBestTrkPt; */
    //electron
    Int_t          nEle;
    vector<int>    *eleCharge;
    vector<int>    *eleChargeConsistent;
    vector<float>  *eleEn;
    vector<float>  *eleSCEn;
    /* vector<float>  eleESEn; */
    /* vector<float>  eleESEnP1; */
    /* vector<float>  eleESEnP2; */
    /* vector<float>  eleESEnP1Raw; */
    /* vector<float>  eleESEnP2Raw; */
    vector<float>  *eleD0;
    vector<float>  *eleDz;
    vector<float>  *elePt;
    vector<float>  *eleEta;
    vector<float>  *elePhi;
    vector<float>  *eleR9;
    vector<float>  *eleSCEta;
    vector<float>  *eleSCPhi;
    vector<float>  *eleSCRawEn;
    /* vector<float>  eleSCEtaWidth; */
    /* vector<float>  eleSCPhiWidth; */
    vector<float>  *eleHoverE;
    vector<float>  *eleEoverP;
    /* vector<float>  eleEoverPout; */
    /* vector<float>  eleEoverPInv; */
    vector<float>  *eleBrem;
    vector<float>  *eledEtaAtVtx;
    vector<float>  *eledPhiAtVtx;
    /* vector<float>  eledEtaAtCalo; */
    vector<float>  *eleSigmaIEtaIEta;
    vector<float>  *eleSigmaIEtaIPhi;
    vector<float>  *eleSigmaIPhiIPhi;
    vector<float>  *eleSigmaIEtaIEtaFull5x5;
    vector<float>  *eleSigmaIPhiIPhiFull5x5;
    vector<int>    *eleConvVeto;
    vector<int>    *eleMissHits;
    vector<float>  *eleESEffSigmaRR;
    vector<float>  *elePFChIso;
    vector<float>  *elePFPhoIso;
    vector<float>  *elePFNeuIso;
    vector<float>  *elePFPUIso;
    vector<float>  *elePFClusEcalIso;
    vector<float>  *elePFClusHcalIso;
    /* vector<float>  eleIDMVANonTrg; */
    vector<float>  *eleIDMVATrg;
    /* vector<float>  eledEtaseedAtVtx; */
    /* vector<float>  eleE1x5; */
    /* vector<float>  eleE2x5; */
    /* vector<float>  eleE5x5; */
    /* vector<float>  eleE1x5Full5x5; */
    /* vector<float>  eleE2x5Full5x5; */
    /* vector<float>  eleE5x5Full5x5; */
    /* vector<float>  eleR9Full5x5; */
    /* vector<int>    eleEcalDrivenSeed; */
    /* vector<float>  eleDr03EcalRecHitSumEt; */
    /* vector<float>  eleDr03HcalDepth1TowerSumEt; */
    /* vector<float>  eleDr03HcalDepth2TowerSumEt; */
    /* vector<float>  eleDr03HcalTowerSumEt; */
    vector<float>  *eleDr03TkSumPt;
    /* vector<float>  elecaloEnergy; */
    /* vector<float>  eleTrkdxy; */
    /* vector<float>  eleKFHits; */
    /* vector<float>  eleKFChi2; */
    /* vector<float>  eleGSFChi2; */
    /* vector<vector<float> > eleGSFPt; */
    /* vector<vector<float> > eleGSFEta; */
    /* vector<vector<float> > eleGSFPhi; */
    /* vector<vector<float> > eleGSFCharge; */
    /* vector<vector<int> >   eleGSFHits; */
    /* vector<vector<int> >   eleGSFMissHits; */
    /* vector<vector<int> >   eleGSFNHitsMax; */
    /* vector<vector<float> > eleGSFVtxProb; */
    /* vector<vector<float> > eleGSFlxyPV; */
    /* vector<vector<float> > eleGSFlxyBS; */
    /* vector<vector<float> > eleBCEn; */
    /* vector<vector<float> > eleBCEta; */
    /* vector<vector<float> > eleBCPhi; */
    /* vector<vector<float> > eleBCS25; */
    /* vector<vector<float> > eleBCS15; */
    /* vector<vector<float> > eleBCSieie; */
    /* vector<vector<float> > eleBCSieip; */
    /* vector<vector<float> > eleBCSipip; */
    /* vector<vector<float> > eleESEnEta; */
    /* vector<vector<float> > eleESEnPhi; */
    /* vector<vector<int> >   eleESEnZ; */
    /* vector<vector<int> >   eleESEnP; */
    /* vector<vector<int> >   eleESEnX; */
    /* vector<vector<int> >   eleESEnY; */
    /* vector<vector<int> >   eleESEnS; */
    /* vector<vector<float> > eleESEnE; */
    /* Int_t nGSFTrk; */
    /* vector<float> gsfPt; */
    /* vector<float> gsfEta; */
    /* vector<float> gsfPhi; */
    vector<int>    *eleFiredTrgs;
    vector<UShort_t> *eleIDbit; 

    //Photon

    Int_t          nPho;
    vector<float>  *phoE;
    vector<float>  *phoEt;
    vector<float>  *phoEta;
    vector<float>  *phoPhi;
    vector<float>  *phoSCE;
    vector<float>  *phoSCRawE;
    /* vector<float>  phoESEn; */
    /* vector<float>  phoESEnP1; */
    /* vector<float>  phoESEnP2; */
    vector<float>  *phoSCEta;
    vector<float>  *phoSCPhi;
    vector<float>  *phoSCEtaWidth;
    vector<float>  *phoSCPhiWidth;
    vector<float>  *phoSCBrem;
    vector<int>    *phohasPixelSeed;
    vector<int>    *phoEleVeto;
    vector<float>  *phoR9;
    vector<float>  *phoHoverE;
    vector<float>  *phoSigmaIEtaIEta;
    vector<float>  *phoSigmaIEtaIPhi;
    vector<float>  *phoSigmaIPhiIPhi;
    /* vector<float>  phoE1x3; */
    /* vector<float>  phoE1x5; */
    /* vector<float>  phoE2x2; */
    /* vector<float>  phoE2x5Max; */
    /* vector<float>  phoE5x5; */
    vector<float>  *phoESEffSigmaRR;
    vector<float>  *phoSigmaIEtaIEtaFull5x5;
    vector<float>  *phoSigmaIEtaIPhiFull5x5;
    vector<float>  *phoSigmaIPhiIPhiFull5x5;
    /* vector<float>  phoE1x3Full5x5; */
    /* vector<float>  phoE1x5Full5x5; */
    /* vector<float>  phoE2x2Full5x5; */
    /* vector<float>  phoE2x5MaxFull5x5; */
    /* vector<float>  phoE5x5Full5x5; */
    vector<float>  *phoR9Full5x5;
    vector<float>  *phoPFChIso;
    vector<float>  *phoPFPhoIso;
    vector<float>  *phoPFNeuIso;
    vector<float>  *phoPFChWorstIso; 
    /* vector<float>  phoPFChIsoFrix1; */
    /* vector<float>  phoPFChIsoFrix2; */
    /* vector<float>  phoPFChIsoFrix3; */
    /* vector<float>  phoPFChIsoFrix4; */
    /* vector<float>  phoPFChIsoFrix5; */
    /* vector<float>  phoPFChIsoFrix6; */
    /* vector<float>  phoPFChIsoFrix7; */
    /* vector<float>  phoPFChIsoFrix8; */
    /* vector<float>  phoPFPhoIsoFrix1; */
    /* vector<float>  phoPFPhoIsoFrix2; */
    /* vector<float>  phoPFPhoIsoFrix3; */
    /* vector<float>  phoPFPhoIsoFrix4; */
    /* vector<float>  phoPFPhoIsoFrix5; */
    /* vector<float>  phoPFPhoIsoFrix6; */
    /* vector<float>  phoPFPhoIsoFrix7; */
    /* vector<float>  phoPFPhoIsoFrix8; */
    /* vector<float>  phoPFNeuIsoFrix1; */
    /* vector<float>  phoPFNeuIsoFrix2; */
    /* vector<float>  phoPFNeuIsoFrix3; */
    /* vector<float>  phoPFNeuIsoFrix4; */
    /* vector<float>  phoPFNeuIsoFrix5; */
    /* vector<float>  phoPFNeuIsoFrix6; */
    /* vector<float>  phoPFNeuIsoFrix7; */
    /* vector<float>  phoPFNeuIsoFrix8; */
    /* vector<float>  phoSeedBCE; */
    /* vector<float>  phoSeedBCEta; */
    vector<float>  *phoIDMVA;
    vector<Int_t>  *phoFiredSingleTrgs;
    vector<Int_t>  *phoFiredDoubleTrgs;
    /* vector<float>  phoEcalRecHitSumEtConeDR03; */
    /* vector<float>  phohcalDepth1TowerSumEtConeDR03; */
    /* vector<float>  phohcalDepth2TowerSumEtConeDR03; */
    /* vector<float>  phohcalTowerSumEtConeDR03; */
    /* vector<float>  photrkSumPtHollowConeDR03; */
    /* vector<float>  photrkSumPtSolidConeDR03; */
    vector<UShort_t> *phoIDbit;

    //jets (minimal stuff)
    Int_t         nJet;
    vector<float> *jetPt;
    vector<float> *jetEn;
    vector<float> *jetEta;
    vector<float> *jetPhi;
    vector<float> *jetRawPt;
    vector<float> *jetRawEn;
    vector<float> *jetMt;
    vector<float> *jetJECUnc;
    vector<int>   *jetFiredTrgs;
    //gen-info for ak4
    vector<int>   *jetGenJetIndex;
    vector<float> *jetGenJetEn;
    vector<float> *jetGenJetPt;
    vector<float> *jetGenJetEta;
    vector<float> *jetGenJetPhi;
    vector<int>   *jetGenPartonID;
    vector<float> *jetGenEn;
    vector<float> *jetGenPt;
    vector<float> *jetGenEta;
    vector<float> *jetGenPhi;
    vector<int>   *jetGenPartonMomID;

    //met (minimal stuff)
    Int_t metFilters;
    float genMET;
    float genMETPhi;
    float pfMET;
    float pfMETPhi;
    float pfMETsumEt;
    float pfMETmEtSig;
    float pfMETSig;

  };

  
  InputTreeLeaves treeLeaf;
  //Output tree additional leaves.
  
  struct OutputTreeLeaves{
    vector<float> eeMass;
    vector<float> eegMass;
    vector<float> gammaEt;
    vector<float> deltaR;
  };
  
  OutputTreeLeaves outLeaf;

  //Global branches
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
  TBranch *b_nMu;
  TBranch *b_muPt;
  TBranch *b_muEn;
  TBranch *b_muEta;
  TBranch *b_muPhi;
  TBranch *b_muCharge;
  TBranch *b_muType;
  TBranch *b_muIsLooseID;
  TBranch *b_muIsMediumID;
  TBranch *b_muIsTightID;
  TBranch *b_muIsSoftID;
  TBranch *b_muIsHighPtID;

  //gen info (if isMC==true)
    // (local) variables associated with tree branches
  TBranch *b_pdf;
  TBranch *b_pthat;
  TBranch *b_processID;
  TBranch *b_genWeight;
  TBranch *b_genHT;

  TBranch *b_nPUInfo;
  TBranch *b_nPU;
  TBranch *b_puBX;
  TBranch *b_puTrue;

  TBranch *b_nMC;          
  TBranch *b_mcPID;
  TBranch *b_mcVtx;
  TBranch *b_mcVty;
  TBranch *b_mcVtz;
  TBranch *b_mcPt;
  TBranch *b_mcMass;
  TBranch *b_mcEta;
  TBranch *b_mcPhi;
  TBranch *b_mcE;
  TBranch *b_mcEt;
  TBranch *b_mcGMomPID;
  TBranch *b_mcMomPID;
  TBranch *b_mcMomPt;
  TBranch *b_mcMomMass;
  TBranch *b_mcMomEta;
  TBranch *b_mcMomPhi;
  TBranch *b_mcIndex;
  TBranch *b_mcStatusFlag;
  TBranch *b_mcParentage;
  TBranch *b_mcStatus;
  TBranch *b_mcCalIsoDR03;
  TBranch *b_mcTrkIsoDR03;
  TBranch *b_mcCalIsoDR04;
  TBranch *b_mcTrkIsoDR04;


    //global muon? PF muon ?
  TBranch *b_muD0;
  TBranch *b_muDz;
  TBranch *b_muChi2NDF;
    /* vector<float>  *muInnerD0; */
    /* vector<float>  *muInnerDz; */
  TBranch *b_muTrkLayers;
  /* vector<int>    *muPixelLayers; */
  TBranch *b_muPixelHits;
  /* vector<int>    *muMuonHits; */
  TBranch *b_muStations;
    /* vector<int>    *muMatches; */
    /* vector<int>    *muTrkQuality; */
    /* vector<float>  *muIsoTrk; */
  TBranch *b_muPFChIso;
  TBranch *b_muPFPhoIso;
  TBranch *b_muPFNeuIso;
  TBranch *b_muPFPUIso;
  TBranch *b_muFiredTrgs;
    /* vector<float>  *muInnervalidFraction; */
    /* vector<float>  *musegmentCompatibility; */
    /* vector<float>  *muchi2LocalPosition; */
    /* vector<float>  *mutrkKink; */
    /* vector<float>  *muBestTrkPtError; */
    /* vector<float>  *muBestTrkPt; */
  //electron branches
  TBranch  *b_nEle;		 
  TBranch  *b_eleCharge;           // /I
  TBranch  *b_eleChargeConsistent; // /I
  TBranch  *b_eleEn;		 
  TBranch  *b_eleSCEn;           
  /* vector<float>  eleESEn; */
  /* vector<float>  eleESEnP1; */
  /* vector<float>  eleESEnP2; */
  /* vector<float>  eleESEnP1Raw; */
  /* vector<float>  eleESEnP2Raw; */
 TBranch  *b_eleD0;
 TBranch  *b_eleDz;
 TBranch  *b_elePt;
 TBranch  *b_eleEta;
 TBranch  *b_elePhi;
 TBranch  *b_eleR9;
 TBranch  *b_eleSCEta;
 TBranch  *b_eleSCPhi;
 TBranch  *b_eleSCRawEn;
  /* vector<float>  eleSCEtaWidth; */
  /* vector<float>  eleSCPhiWidth; */
 TBranch  *b_eleHoverE;
 TBranch  *b_eleEoverP;
  /* vector<float>  eleEoverPout; */
  /* vector<float>  eleEoverPInv; */
 TBranch  *b_eleBrem;
 TBranch *b_eledEtaAtVtx;
 TBranch *b_eledPhiAtVtx;
  /* vector<float>  eledEtaAtCalo; */
 TBranch  *b_eleSigmaIEtaIEta;
 TBranch  *b_eleSigmaIEtaIPhi;
 TBranch  *b_eleSigmaIPhiIPhi;
 TBranch  *b_eleSigmaIEtaIEtaFull5x5;
 TBranch  *b_eleSigmaIPhiIPhiFull5x5;
 TBranch  *b_eleConvVeto;  // /I
 TBranch  *b_eleMissHits; // /I
 TBranch  *b_eleESEffSigmaRR;
 TBranch  *b_elePFChIso;
 TBranch  *b_elePFPhoIso;
 TBranch  *b_elePFNeuIso;
 TBranch  *b_elePFPUIso;
 TBranch  *b_elePFClusEcalIso;
 TBranch  *b_elePFClusHcalIso;
  /* vector<float>  eleIDMVANonTrg; */
 TBranch *b_eleIDMVATrg;
  /* vector<float>  eledEtaseedAtVtx; */
  /* vector<float>  eleE1x5; */
  /* vector<float>  eleE2x5; */
  /* vector<float>  eleE5x5; */
  /* vector<float>  eleE1x5Full5x5; */
  /* vector<float>  eleE2x5Full5x5; */
  /* vector<float>  eleE5x5Full5x5; */
  /* vector<float>  eleR9Full5x5; */
  /* vector<int>    eleEcalDrivenSeed; */
  /* vector<float>  eleDr03EcalRecHitSumEt; */
  /* vector<float>  eleDr03HcalDepth1TowerSumEt; */
  /* vector<float>  eleDr03HcalDepth2TowerSumEt; */
  /* vector<float>  eleDr03HcalTowerSumEt; */
 TBranch *b_eleDr03TkSumPt;
  /* vector<float>  elecaloEnergy; */
  /* vector<float>  eleTrkdxy; */
  /* vector<float>  eleKFHits; */
  /* vector<float>  eleKFChi2; */
  /* vector<float>  eleGSFChi2; */
  /* vector<vector<float> > eleGSFPt; */
  /* vector<vector<float> > eleGSFEta; */
  /* vector<vector<float> > eleGSFPhi; */
  /* vector<vector<float> > eleGSFCharge; */
  /* vector<vector<int> >   eleGSFHits; */
  /* vector<vector<int> >   eleGSFMissHits; */
  /* vector<vector<int> >   eleGSFNHitsMax; */
  /* vector<vector<float> > eleGSFVtxProb; */
  /* vector<vector<float> > eleGSFlxyPV; */
  /* vector<vector<float> > eleGSFlxyBS; */
  /* vector<vector<float> > eleBCEn; */
  /* vector<vector<float> > eleBCEta; */
  /* vector<vector<float> > eleBCPhi; */
  /* vector<vector<float> > eleBCS25; */
  /* vector<vector<float> > eleBCS15; */
  /* vector<vector<float> > eleBCSieie; */
  /* vector<vector<float> > eleBCSieip; */
  /* vector<vector<float> > eleBCSipip; */
  /* vector<vector<float> > eleESEnEta; */
  /* vector<vector<float> > eleESEnPhi; */
  /* vector<vector<int> >   eleESEnZ; */
  /* vector<vector<int> >   eleESEnP; */
  /* vector<vector<int> >   eleESEnX; */
  /* vector<vector<int> >   eleESEnY; */
  /* vector<vector<int> >   eleESEnS; */
  /* vector<vector<float> > eleESEnE; */
  /* Int_t nGSFTrk; */
  /* vector<float> gsfPt; */
  /* vector<float> gsfEta; */
  /* vector<float> gsfPhi; */
 TBranch  *b_eleFiredTrgs; // /I
 TBranch  *b_eleIDbit;  // /U

 //photon
 TBranch *b_nPho;		     
 TBranch *b_phoE;		     
 TBranch *b_phoEt;		     
 TBranch *b_phoEta;		     
 TBranch *b_phoPhi;		     
 TBranch *b_phoSCE;		     
 TBranch *b_phoSCRawE;		     
 TBranch *b_phoSCEta;		     
 TBranch *b_phoSCPhi;		     
 TBranch *b_phoSCEtaWidth;	     
 TBranch *b_phoSCPhiWidth;	     
 TBranch *b_phoSCBrem;		     
 TBranch *b_phohasPixelSeed;	     
 TBranch *b_phoEleVeto;		     
 TBranch *b_phoR9;
 TBranch *b_phoPFChWorstIso;
 TBranch *b_phoHoverE;		     
 TBranch *b_phoSigmaIEtaIEta;	     
 TBranch *b_phoSigmaIEtaIPhi;	     
 TBranch *b_phoSigmaIPhiIPhi;	     
 TBranch *b_phoESEffSigmaRR;	     
 TBranch *b_phoSigmaIEtaIEtaFull5x5; 
 TBranch *b_phoSigmaIEtaIPhiFull5x5; 
 TBranch *b_phoSigmaIPhiIPhiFull5x5; 
 TBranch *b_phoPFChIso;		     
 TBranch *b_phoPFPhoIso;	     
 TBranch *b_phoPFNeuIso;	     
 TBranch *b_phoFiredSingleTrgs;	     
 TBranch *b_phoFiredDoubleTrgs;   
 TBranch *b_phoIDbit;   

 /* vector<float>  phoESEn; */
 /* vector<float>  phoESEnP1; */
 /* vector<float>  phoESEnP2; */

 /* vector<float>  phoE1x3; */
 /* vector<float>  phoE1x5; */
 /* vector<float>  phoE2x2; */
 /* vector<float>  phoE2x5Max; */
 /* vector<float>  phoE5x5; */

 /* vector<float>  phoE1x3Full5x5; */
 /* vector<float>  phoE1x5Full5x5; */
 /* vector<float>  phoE2x2Full5x5; */
 /* vector<float>  phoE2x5MaxFull5x5; */
 /* vector<float>  phoE5x5Full5x5; */
  TBranch *b_phoR9Full5x5;


 /* vector<float>  phoPFChIsoFrix1; */
 /* vector<float>  phoPFChIsoFrix2; */
 /* vector<float>  phoPFChIsoFrix3; */
 /* vector<float>  phoPFChIsoFrix4; */
 /* vector<float>  phoPFChIsoFrix5; */
 /* vector<float>  phoPFChIsoFrix6; */
 /* vector<float>  phoPFChIsoFrix7; */
 /* vector<float>  phoPFChIsoFrix8; */
 /* vector<float>  phoPFPhoIsoFrix1; */
 /* vector<float>  phoPFPhoIsoFrix2; */
 /* vector<float>  phoPFPhoIsoFrix3; */
 /* vector<float>  phoPFPhoIsoFrix4; */
 /* vector<float>  phoPFPhoIsoFrix5; */
 /* vector<float>  phoPFPhoIsoFrix6; */
 /* vector<float>  phoPFPhoIsoFrix7; */
 /* vector<float>  phoPFPhoIsoFrix8; */
 /* vector<float>  phoPFNeuIsoFrix1; */
 /* vector<float>  phoPFNeuIsoFrix2; */
 /* vector<float>  phoPFNeuIsoFrix3; */
 /* vector<float>  phoPFNeuIsoFrix4; */
 /* vector<float>  phoPFNeuIsoFrix5; */
 /* vector<float>  phoPFNeuIsoFrix6; */
 /* vector<float>  phoPFNeuIsoFrix7; */
 /* vector<float>  phoPFNeuIsoFrix8; */
 /* vector<float>  phoSeedBCE; */
 /* vector<float>  phoSeedBCEta; */
 TBranch *b_phoIDMVA;

 /* vector<float>  phoEcalRecHitSumEtConeDR03; */
 /* vector<float>  phohcalDepth1TowerSumEtConeDR03; */
 /* vector<float>  phohcalDepth2TowerSumEtConeDR03; */
 /* vector<float>  phohcalTowerSumEtConeDR03; */
 /* vector<float>  photrkSumPtHollowConeDR03; */
 /* vector<float>  photrkSumPtSolidConeDR03; */

 //jets (minimal stuff)
 TBranch *b_nJet;
 TBranch *b_jetPt;
 TBranch *b_jetEn;
 TBranch *b_jetEta;
 TBranch *b_jetPhi;
 TBranch *b_jetRawPt;
 TBranch *b_jetRawEn;
 TBranch *b_jetMt;
 TBranch *b_jetJECUnc;
 TBranch *b_jetFiredTrgs;
 TBranch *b_jetGenJetIndex;
 TBranch *b_jetGenJetEn;
 TBranch *b_jetGenJetPt;
 TBranch *b_jetGenJetEta;
 TBranch *b_jetGenJetPhi;
 TBranch *b_jetGenPartonID;
 TBranch *b_jetGenEn;
 TBranch *b_jetGenPt;
 TBranch *b_jetGenEta;
 TBranch *b_jetGenPhi;
 TBranch *b_jetGenPartonMomID;

 //met (minimal stuff)
 TBranch *b_metFilters;
 TBranch *b_genMET;
 TBranch *b_genMETPhi;
 TBranch *b_pfMET;
 TBranch *b_pfMETPhi;
 TBranch *b_pfMETsumEt;
 TBranch *b_pfMETmEtSig;
 TBranch *b_pfMETSig;

 //extra
 TBranch *b_eeMass;
 TBranch *b_eegMass;
 TBranch *b_gammaEt;
 TBranch *b_deltaR;
};


#endif
