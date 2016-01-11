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
    /* vector<float>  eledEtaAtVtx; */
    /* vector<float>  eledPhiAtVtx; */
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
    /* vector<float>  eleIDMVATrg; */
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
    /* vector<float>  eleDr03TkSumPt; */
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
  };

  
  InputTreeLeaves treeLeaf;

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
  /* vector<float>  eledEtaAtVtx; */
  /* vector<float>  eledPhiAtVtx; */
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
  /* vector<float>  eleIDMVATrg; */
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
  /* vector<float>  eleDr03TkSumPt; */
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
};


#endif
