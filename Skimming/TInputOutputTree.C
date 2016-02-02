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

void TInputOutputTree::InitInput(TTree *tree, bool isMC){
  //global
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
  //gen info

  //gen info (if isMC==true)
  // (local) variables associated with tree branches
  if(isMC){
    treeLeaf.pdf=0;
    treeLeaf.pthat=0;
    treeLeaf.processID=0;
    treeLeaf.genWeight=0;
    treeLeaf.genHT=0;

    treeLeaf.nPUInfo=0;
    treeLeaf.nPU=0;
    treeLeaf.puBX=0;
    treeLeaf.puTrue=0;

    treeLeaf.nMC=0;          
    treeLeaf.mcPID=0;
    treeLeaf.mcVtx=0;
    treeLeaf.mcVty=0;
    treeLeaf.mcVtz=0;
    treeLeaf.mcPt=0;
    treeLeaf.mcMass=0;
    treeLeaf.mcEta=0;
    treeLeaf.mcPhi=0;
    treeLeaf.mcE=0;
    treeLeaf.mcEt=0;
    treeLeaf.mcGMomPID=0;
    treeLeaf.mcMomPID=0;
    treeLeaf.mcMomPt=0;
    treeLeaf.mcMomMass=0;
    treeLeaf.mcMomEta=0;
    treeLeaf.mcMomPhi=0;
    treeLeaf.mcIndex=0;
    treeLeaf.mcStatusFlag=0;
    treeLeaf.mcParentage=0;
    treeLeaf.mcStatus=0;
    treeLeaf.mcCalIsoDR03=0;
    treeLeaf.mcTrkIsoDR03=0;
    treeLeaf.mcCalIsoDR04=0;
    treeLeaf.mcTrkIsoDR04=0;
  }
  //muon   
  treeLeaf.nMu=0;
  treeLeaf.muPt=0;
  treeLeaf.muEn=0;
  treeLeaf.muEta=0;
  treeLeaf.muPhi=0;
  treeLeaf.muCharge=0;
  treeLeaf.muType=0;
  treeLeaf.muIsLooseID=0;
  treeLeaf.muIsMediumID=0;
  treeLeaf.muIsTightID=0;
  treeLeaf.muIsSoftID=0;
  treeLeaf.muIsHighPtID=0;
  treeLeaf.muD0=0;
  treeLeaf.muDz=0;
  treeLeaf.muChi2NDF=0;
  /* treeLeaf.muInnerD0=0; */
  /* treeLeaf.muInnerDz=0; */
  treeLeaf.muTrkLayers=0;
  /* treeLeaf.muPixelLayers=0; */
  treeLeaf.muPixelHits=0;
  /* treeLeaf.muTreeLeaf.MuonHits=0; */
  treeLeaf.muStations=0;
  /* treeLeaf.muMatches=0; */
  /* treeLeaf.muTrkQuality=0; */
  /* treeLeaf.muIsoTrk=0; */
  treeLeaf.muPFChIso=0;
  treeLeaf.muPFPhoIso=0;
  treeLeaf.muPFNeuIso=0;
  treeLeaf.muPFPUIso=0;
  treeLeaf.muFiredTrgs=0;
  /* treeLeaf.muInnervalidFraction=0; */
  /* treeLeaf.musegmentCompatibility=0; */
  /* treeLeaf.muchi2LocalPosition=0; */
  /* treeLeaf.mutrkKink=0; */
  /* treeLeaf.muBestTrkPtError=0; */
  /* treeLeaf.muBestTrkPt=0; */
  //electron
  treeLeaf.nEle=0;
  treeLeaf.eleCharge=0;
  treeLeaf.eleChargeConsistent=0;
  treeLeaf.eleEn=0;
  treeLeaf.eleSCEn=0;
  /*  treeLeaf.eleESEn=0; */
  /*  treeLeaf.eleESEnP1=0; */
  /*  treeLeaf.eleESEnP2=0; */
  /*  treeLeaf.eleESEnP1Raw=0; */
  /*  treeLeaf.eleESEnP2Raw=0; */
  treeLeaf.eleD0=0;
  treeLeaf.eleDz=0;
  treeLeaf.elePt=0;
  treeLeaf.eleEta=0;
  treeLeaf.elePhi=0;
  treeLeaf.eleR9=0;
  treeLeaf.eleSCEta=0;
  treeLeaf.eleSCPhi=0;
  treeLeaf.eleSCRawEn=0;
  /*   treeLeaf.eleSCEtaWidth=0; */
  /*   treeLeaf.eleSCPhiWidth=0; */
  treeLeaf.eleHoverE=0;
  treeLeaf.eleEoverP=0;
  /*  treeLeaf.eleEoverPout=0; */
  /*  treeLeaf.eleEoverPInv=0; */
  treeLeaf.eleBrem=0;
   treeLeaf.eledEtaAtVtx=0;
   treeLeaf.eledPhiAtVtx=0;
  /*  treeLeaf.eledEtaAtCalo=0; */
  treeLeaf.eleSigmaIEtaIEta=0;
  treeLeaf.eleSigmaIEtaIPhi=0;
  treeLeaf.eleSigmaIPhiIPhi=0;
  treeLeaf.eleSigmaIEtaIEtaFull5x5=0;
  treeLeaf.eleSigmaIPhiIPhiFull5x5=0;
  treeLeaf.eleConvVeto=0;
  treeLeaf.eleMissHits=0;
  treeLeaf.eleESEffSigmaRR=0;
  treeLeaf.elePFChIso=0;
  treeLeaf.elePFPhoIso=0;
  treeLeaf.elePFNeuIso=0;
  treeLeaf.elePFPUIso=0;
  treeLeaf.elePFClusEcalIso=0;
  treeLeaf.elePFClusHcalIso=0;
  /*  treeLeaf.eleIDMVANonTrg=0; */
   treeLeaf.eleIDMVATrg=0;
  /*  treeLeaf.eledEtaseedAtVtx=0; */
  /*  treeLeaf.eleE1x5=0; */
  /*  treeLeaf.eleE2x5=0; */
  /*  treeLeaf.eleE5x5=0; */
  /*  treeLeaf.eleE1x5Full5x5=0; */
  /*  treeLeaf.eleE2x5Full5x5=0; */
  /*  treeLeaf.eleE5x5Full5x5=0; */
  /*  treeLeaf.eleR9Full5x5=0; */
  /*  treeLeaf.eleEcalDrivenSeed=0; */
  /*  treeLeaf.eleDr03EcalRecHitSumEt=0; */
  /*  treeLeaf.eleDr03HcalDepth1TowerSumEt=0; */
  /*  treeLeaf.eleDr03HcalDepth2TowerSumEt=0; */
  /*  treeLeaf.eleDr03HcalTowerSumEt=0; */
   treeLeaf.eleDr03TkSumPt=0;
  /*  treeLeaf.elecaloEnergy=0; */
  /*  treeLeaf.eleTrkdxy=0; */
  /*  treeLeaf.eleKFHits=0; */
  /*  treeLeaf.eleKFChi2=0; */
  /*  treeLeaf.eleGSFChi2=0; */
  /*  treeLeaf.eleGSFPt=0; */
  /*  treeLeaf.eleGSFEta=0; */
  /*  treeLeaf.eleGSFPhi=0; */
  /*  treeLeaf.eleGSFCharge=0; */
  /*  treeLeaf.eleGSFHits=0; */
  /*  treeLeaf.eleGSFMissHits=0; */
  /*  treeLeaf.eleGSFNHitsMax=0; */
  /*  treeLeaf.eleGSFVtxProb=0; */
  /*  treeLeaf.eleGSFlxyPV=0; */
  /*  treeLeaf.eleGSFlxyBS=0; */
  /*  treeLeaf.eleBCEn=0; */
  /*  treeLeaf.eleBCEta=0; */
  /*  treeLeaf.eleBCPhi=0; */
  /*  treeLeaf.eleBCS25=0; */
  /*  treeLeaf.eleBCS15=0; */
  /*  treeLeaf.eleBCSieie=0; */
  /*  treeLeaf.eleBCSieip=0; */
  /*  treeLeaf.eleBCSipip=0; */
  /*  treeLeaf.eleESEnEta=0; */
  /*  treeLeaf.eleESEnPhi=0; */
  /*  treeLeaf.eleESEnZ=0; */
  /*  treeLeaf.eleESEnP=0; */
  /*  treeLeaf.eleESEnX=0; */
  /*  treeLeaf.eleESEnY=0; */
  /*  treeLeaf.eleESEnS=0; */
  /*  treeLeaf.eleESEnE=0; */
  /*  treeLeaf.nGSFTrk=0; */
  /*  treeLeaf.gsfPt=0; */
  /*  treeLeaf.gsfEta=0; */
  /*  treeLeaf.gsfPhi=0; */
  treeLeaf.eleFiredTrgs=0;
  treeLeaf.eleIDbit=0; 


  //photon
  treeLeaf.nPho=0;
  treeLeaf.phoE=0;
  treeLeaf.phoEt=0;
  treeLeaf.phoEta=0;
  treeLeaf.phoPhi=0;
  treeLeaf.phoSCE=0;
  treeLeaf.phoSCRawE=0;
  //  phoESEn=0; 
  //  phoESEnP1=0;
  //  phoESEnP2=0; 
  treeLeaf.phoSCEta=0;
  treeLeaf.phoSCPhi=0;
  treeLeaf.phoSCEtaWidth=0;
  treeLeaf.phoSCPhiWidth=0;
  treeLeaf.phoSCBrem=0;
  treeLeaf.phohasPixelSeed=0;
  treeLeaf.phoEleVeto=0;
  treeLeaf.phoR9=0;
  treeLeaf.phoHoverE=0;
  treeLeaf.phoSigmaIEtaIEta=0;
  treeLeaf.phoSigmaIEtaIPhi=0;
  treeLeaf.phoSigmaIPhiIPhi=0;
  /* vector<float>  phoE1x3=0; */
  /* vector<float>  phoE1x5=0; */
  /* vector<float>  phoE2x2=0; */
  /* vector<float>  phoE2x5Max=0; */
  /* vector<float>  phoE5x5=0; */
  treeLeaf.phoESEffSigmaRR=0;
  treeLeaf.phoSigmaIEtaIEtaFull5x5=0;
  treeLeaf.phoSigmaIEtaIPhiFull5x5=0;
  treeLeaf.phoSigmaIPhiIPhiFull5x5=0;
  /* vector<float>  phoE1x3Full5x5=0; */
  /* vector<float>  phoE1x5Full5x5=0; */
  /* vector<float>  phoE2x2Full5x5=0; */
  /* vector<float>  phoE2x5MaxFull5x5=0; */
  /* vector<float>  phoE5x5Full5x5=0; */
  /* vector<float>  phoR9Full5x5=0; */
  treeLeaf.phoPFChIso=0;
  treeLeaf.phoPFChWorstIso=0;
  treeLeaf.phoPFPhoIso=0;
  treeLeaf.phoPFNeuIso=0;
  /* vector<float>  phoPFChIsoFrix1=0; */
  /* vector<float>  phoPFChIsoFrix2=0; */
  /* vector<float>  phoPFChIsoFrix3=0; */
  /* vector<float>  phoPFChIsoFrix4=0; */
  /* vector<float>  phoPFChIsoFrix5=0; */
  /* vector<float>  phoPFChIsoFrix6=0; */
  /* vector<float>  phoPFChIsoFrix7=0; */
  /* vector<float>  phoPFChIsoFrix8=0; */
  /* vector<float>  phoPFPhoIsoFrix1=0; */
  /* vector<float>  phoPFPhoIsoFrix2=0; */
  /* vector<float>  phoPFPhoIsoFrix3=0; */
  /* vector<float>  phoPFPhoIsoFrix4=0; */
  /* vector<float>  phoPFPhoIsoFrix5=0; */
  /* vector<float>  phoPFPhoIsoFrix6=0; */
  /* vector<float>  phoPFPhoIsoFrix7=0; */
  /* vector<float>  phoPFPhoIsoFrix8=0; */
  /* vector<float>  phoPFNeuIsoFrix1=0; */
  /* vector<float>  phoPFNeuIsoFrix2=0; */
  /* vector<float>  phoPFNeuIsoFrix3=0; */
  /* vector<float>  phoPFNeuIsoFrix4=0; */
  /* vector<float>  phoPFNeuIsoFrix5=0; */
  /* vector<float>  phoPFNeuIsoFrix6=0; */
  /* vector<float>  phoPFNeuIsoFrix7=0; */
  /* vector<float>  phoPFNeuIsoFrix8=0; */
  /* vector<float>  phoSeedBCE=0; */
  /* vector<float>  phoSeedBCEta=0; */
  treeLeaf.phoIDMVA=0;
  treeLeaf.phoFiredSingleTrgs=0;
  treeLeaf.phoFiredDoubleTrgs=0;
  /* vector<float>  phoEcalRecHitSumEtConeDR03=0; */
  /* vector<float>  phohcalDepth1TowerSumEtConeDR03=0; */
  /* vector<float>  phohcalDepth2TowerSumEtConeDR03=0; */
  /* vector<float>  phohcalTowerSumEtConeDR03=0; */
  /* vector<float>  photrkSumPtHollowConeDR03=0; */
  /* vector<float>  photrkSumPtSolidConeDR03=0; */
  treeLeaf.phoIDbit=0;  

  //jets (minimal stuff)
  treeLeaf.nJet=0;
  treeLeaf.jetPt=0;
  treeLeaf.jetEn=0;
  treeLeaf.jetEta=0;
  treeLeaf.jetPhi=0;
  treeLeaf.jetRawPt=0;
  treeLeaf.jetRawEn=0;
  treeLeaf.jetMt=0;
  treeLeaf.jetJECUnc=0;
  treeLeaf.jetFiredTrgs=0;
  treeLeaf.jetGenJetIndex=0;
  treeLeaf.jetGenJetEn=0;
  treeLeaf.jetGenJetPt=0;
  treeLeaf.jetGenJetEta=0;
  treeLeaf.jetGenJetPhi=0;
  treeLeaf.jetGenPartonID=0;
  treeLeaf.jetGenEn=0;
  treeLeaf.jetGenPt=0;
  treeLeaf.jetGenEta=0;
  treeLeaf.jetGenPhi=0;
  treeLeaf.jetGenPartonMomID=0;

  //met (minimal stuff)
  treeLeaf.metFilters=0;
  treeLeaf.genMET=0;
  treeLeaf.genMETPhi=0;
  treeLeaf.pfMET=0;
  treeLeaf.pfMETPhi=0;
  treeLeaf.pfMETsumEt=0;
  treeLeaf.pfMETmEtSig=0;
  treeLeaf.pfMETSig=0;


  //extra
  //to initialise the members of (struct) TInputOutputTree::OutputTreeLeaves outLeaf
  //  memset(&outLeaf, 0, sizeof(OutputTreeLeaves));


  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  //global event
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
  //gen info
  // (local) variables associated with tree branches
  if(isMC){
    fChain->SetBranchAddress("pdf",&treeLeaf.pdf,&b_pdf);
    fChain->SetBranchAddress("pthat",&treeLeaf.pthat,&b_pthat);
    fChain->SetBranchAddress("processID",&treeLeaf.processID,&b_processID);
    fChain->SetBranchAddress("genWeight",&treeLeaf.genWeight,&b_genWeight);
    fChain->SetBranchAddress("genHT",&treeLeaf.genHT,&b_genHT);

    fChain->SetBranchAddress("nPUInfo",&treeLeaf.nPUInfo,&b_nPUInfo);
    fChain->SetBranchAddress("nPU",&treeLeaf.nPU,&b_nPU);
    fChain->SetBranchAddress("puBX",&treeLeaf.puBX,&b_puBX);
    fChain->SetBranchAddress("puTrue",&treeLeaf.puTrue,&b_puTrue);

    fChain->SetBranchAddress("nMC",         &treeLeaf.nMC,          &b_nMC);
    fChain->SetBranchAddress("mcPID",&treeLeaf.mcPID,  &b_mcPID);
    fChain->SetBranchAddress("mcVtx",	   &treeLeaf.mcVtx,	    &b_mcVtx);	    
    fChain->SetBranchAddress("mcVty",	   &treeLeaf.mcVty,	    &b_mcVty);	    
    fChain->SetBranchAddress("mcVtz",	   &treeLeaf.mcVtz,	    &b_mcVtz);	    
    fChain->SetBranchAddress("mcPt",	   &treeLeaf.mcPt,	    &b_mcPt);	    
    fChain->SetBranchAddress("mcMass",	   &treeLeaf.mcMass,	    &b_mcMass);	    
    fChain->SetBranchAddress("mcEta",	   &treeLeaf.mcEta,	    &b_mcEta);	    
    fChain->SetBranchAddress("mcPhi",	   &treeLeaf.mcPhi,	    &b_mcPhi);	    
    fChain->SetBranchAddress("mcE",	   &treeLeaf.mcE,	    &b_mcE);	    
    fChain->SetBranchAddress("mcEt",	   &treeLeaf.mcEt,	    &b_mcEt);	    
    fChain->SetBranchAddress("mcGMomPID",  &treeLeaf.mcGMomPID,   &b_mcGMomPID);   
    fChain->SetBranchAddress("mcMomPID",   &treeLeaf.mcMomPID,    &b_mcMomPID);    
    fChain->SetBranchAddress("mcMomPt",	   &treeLeaf.mcMomPt,	    &b_mcMomPt);	    
    fChain->SetBranchAddress("mcMomMass",  &treeLeaf.mcMomMass,   &b_mcMomMass);   
    fChain->SetBranchAddress("mcMomEta",   &treeLeaf.mcMomEta,    &b_mcMomEta);    
    fChain->SetBranchAddress("mcMomPhi",   &treeLeaf.mcMomPhi,    &b_mcMomPhi);    
    fChain->SetBranchAddress("mcIndex",	   &treeLeaf.mcIndex,	    &b_mcIndex);	    
    fChain->SetBranchAddress("mcStatusFlag",&treeLeaf.mcStatusFlag,&b_mcStatusFlag);
    fChain->SetBranchAddress("mcParentage",&treeLeaf.mcParentage, &b_mcParentage); 
    fChain->SetBranchAddress("mcStatus",   &treeLeaf.mcStatus,    &b_mcStatus);    
    fChain->SetBranchAddress("mcCalIsoDR03",&treeLeaf.mcCalIsoDR03,&b_mcCalIsoDR03);
    fChain->SetBranchAddress("mcTrkIsoDR03",&treeLeaf.mcTrkIsoDR03,&b_mcTrkIsoDR03);
    fChain->SetBranchAddress("mcCalIsoDR04",&treeLeaf.mcCalIsoDR04,&b_mcCalIsoDR04);
    fChain->SetBranchAddress("mcTrkIsoDR04",&treeLeaf.mcTrkIsoDR04,&b_mcTrkIsoDR04);
  }
  //muon channel
  fChain->SetBranchAddress("nMu",&treeLeaf.nMu,&b_nMu);
  fChain->SetBranchAddress("muPt",&treeLeaf.muPt	    ,&b_muPt);	     
  fChain->SetBranchAddress("muEn",&treeLeaf.muEn	    ,&b_muEn);	     
  fChain->SetBranchAddress("muEta",&treeLeaf.muEta	    ,&b_muEta);	     
  fChain->SetBranchAddress("muPhi",&treeLeaf.muPhi	    ,&b_muPhi);	     
  fChain->SetBranchAddress("muCharge",&treeLeaf.muCharge  ,&b_muCharge);   
  fChain->SetBranchAddress("muType",&treeLeaf.muType	    ,&b_muType);	     
  fChain->SetBranchAddress("muIsLooseID",&treeLeaf.muIsLooseID,&b_muIsLooseID); 
  fChain->SetBranchAddress("muIsMediumID",&treeLeaf.muIsMediumID,&b_muIsMediumID);
  fChain->SetBranchAddress("muIsTightID",&treeLeaf.muIsTightID,&b_muIsTightID); 
  fChain->SetBranchAddress("muIsSoftID",&treeLeaf.muIsSoftID,&b_muIsSoftID);  
  fChain->SetBranchAddress("muIsHighPtID",&treeLeaf.muIsHighPtID,&b_muIsHighPtID);
  fChain->SetBranchAddress("muD0",&treeLeaf.muD0	    ,&b_muD0);	     
  fChain->SetBranchAddress("muDz",&treeLeaf.muDz	    ,&b_muDz);	     
  fChain->SetBranchAddress("muChi2ND ",&treeLeaf.muChi2NDF ,&b_muChi2NDF);   
  fChain->SetBranchAddress("muTrkLayers",&treeLeaf.muTrkLayers,&b_muTrkLayers); 
  fChain->SetBranchAddress("muPixelHits",&treeLeaf.muPixelHits,&b_muPixelHits); 
  fChain->SetBranchAddress("muStations   ",&treeLeaf.muStations,&b_muStations);  
  fChain->SetBranchAddress("muPFChIso    ",&treeLeaf.muPFChIso ,&b_muPFChIso);   
  fChain->SetBranchAddress("muPFPhoIso   ",&treeLeaf.muPFPhoIso,&b_muPFPhoIso);  
  fChain->SetBranchAddress("muPFNeuIso   ",&treeLeaf.muPFNeuIso,&b_muPFNeuIso); 
  fChain->SetBranchAddress("muPFPUIso    ",&treeLeaf.muPFPUIso ,&b_muPFPUIso);   
  fChain->SetBranchAddress("muFiredTrgs  ",&treeLeaf.muFiredTrgs,&b_muFiredTrgs);  

  //electron channel
  fChain->SetBranchAddress("nEle",                    &treeLeaf.nEle,                &b_nEle	     );
  fChain->SetBranchAddress("eleCharge",               &treeLeaf.eleCharge,           &b_eleCharge            );
  fChain->SetBranchAddress("eleChargeConsistent",     &treeLeaf.eleChargeConsistent, &b_eleChargeConsistent    );
  fChain->SetBranchAddress("eleEn",                   &treeLeaf.eleEn,               &b_eleEn   );
  fChain->SetBranchAddress("eleSCEn",                 &treeLeaf.eleSCEn,             &b_eleSCEn );
  fChain->SetBranchAddress("eleD0",                   &treeLeaf.eleD0 , &b_eleD0);
  fChain->SetBranchAddress("eleDz",                   &treeLeaf.eleDz , &b_eleDz);
  fChain->SetBranchAddress("elePt", &treeLeaf.elePt, &b_elePt);
  fChain->SetBranchAddress("eleEta", &treeLeaf.eleEta, &b_eleEta);
  fChain->SetBranchAddress("elePhi", &treeLeaf.elePhi, &b_elePhi);
  fChain->SetBranchAddress("eleR9", &treeLeaf.eleR9, &b_eleR9);
  fChain->SetBranchAddress("eleSCEta", &treeLeaf.eleSCEta, &b_eleSCEta);
  fChain->SetBranchAddress("eleSCPhi", &treeLeaf.eleSCPhi, &b_eleSCPhi);
  fChain->SetBranchAddress("eleSCRawEn", &treeLeaf.eleSCRawEn, &b_eleSCRawEn);
  fChain->SetBranchAddress("eleHoverE",               &treeLeaf.eleHoverE, &b_eleHoverE);
  fChain->SetBranchAddress("eleEoverP",               &treeLeaf.eleEoverP, &b_eleEoverP);
  fChain->SetBranchAddress("eleBrem",                 &treeLeaf.eleBrem, &b_eleBrem);
  fChain->SetBranchAddress("eleSigmaIEtaIEta",        &treeLeaf.eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
  fChain->SetBranchAddress("eleSigmaIEtaIPhi",        &treeLeaf.eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
  fChain->SetBranchAddress("eleSigmaIPhiIPhi",        &treeLeaf.eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
  fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &treeLeaf.eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
  fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &treeLeaf.eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
  fChain->SetBranchAddress("eleConvVeto",             &treeLeaf.eleConvVeto, &b_eleConvVeto);
  fChain->SetBranchAddress("eleMissHits",             &treeLeaf.eleMissHits, &b_eleMissHits);
  fChain->SetBranchAddress("eleESEffSigmaRR",         &treeLeaf.eleESEffSigmaRR, &b_eleESEffSigmaRR);
  fChain->SetBranchAddress("elePFChIso",              &treeLeaf.elePFChIso, &b_elePFChIso);
  fChain->SetBranchAddress("elePFPhoIso",             &treeLeaf.elePFPhoIso, &b_elePFPhoIso);
  fChain->SetBranchAddress("elePFNeuIso",             &treeLeaf.elePFNeuIso, &b_elePFNeuIso);
  fChain->SetBranchAddress("elePFPUIso",              &treeLeaf.elePFPUIso, &b_elePFPUIso);
  fChain->SetBranchAddress("elePFClusEcalIso",        &treeLeaf.elePFClusEcalIso, &b_elePFClusEcalIso);
  fChain->SetBranchAddress("elePFClusHcalIso",        &treeLeaf.elePFClusHcalIso, &b_elePFClusHcalIso);
  fChain->SetBranchAddress("eleFiredTrgs",                &treeLeaf.eleFiredTrgs, &b_eleFiredTrgs);
  fChain->SetBranchAddress("eleIDbit",  &treeLeaf.eleIDbit, &b_eleIDbit);
  fChain->SetBranchAddress("eleIDMVATrg",  &treeLeaf.eleIDMVATrg, &b_eleIDMVATrg);
  fChain->SetBranchAddress("eledEtaAtVtx",  &treeLeaf.eledEtaAtVtx, &b_eledEtaAtVtx);
  fChain->SetBranchAddress("eledPhiAtVtx",  &treeLeaf.eledPhiAtVtx, &b_eledPhiAtVtx);
  fChain->SetBranchAddress("eleDr03TkSumPt",  &treeLeaf.eleDr03TkSumPt, &b_eleDr03TkSumPt);
  
  //photon channel 
  fChain->SetBranchAddress("nPho"			  ,&treeLeaf.nPho	       	   ,&b_nPho);		     
  fChain->SetBranchAddress("phoE"			  ,&treeLeaf.phoE      		   ,&b_phoE);		     
  fChain->SetBranchAddress("phoEt"		  ,&treeLeaf.phoEt		   ,&b_phoEt);		     
  fChain->SetBranchAddress("phoEta"		  ,&treeLeaf.phoEta		   ,&b_phoEta);		     
  fChain->SetBranchAddress("phoPhi"		  ,&treeLeaf.phoPhi		   ,&b_phoPhi);		     
  fChain->SetBranchAddress("phoSCE"		  ,&treeLeaf.phoSCE		   ,&b_phoSCE);		      
  fChain->SetBranchAddress("phoSCRawE"		  ,&treeLeaf.phoSCRawE		   ,&b_phoSCRawE);		      
  fChain->SetBranchAddress("phoSCEta"		  ,&treeLeaf.phoSCEta		   ,&b_phoSCEta);		      
  fChain->SetBranchAddress("phoSCPhi"		  ,&treeLeaf.phoSCPhi		   ,&b_phoSCPhi);		      
  fChain->SetBranchAddress("phoSCEtaWidth"	  ,&treeLeaf.phoSCEtaWidth	   ,&b_phoSCEtaWidth);	      
  fChain->SetBranchAddress("phoSCPhiWidth"	  ,&treeLeaf.phoSCPhiWidth	   ,&b_phoSCPhiWidth);	      
  fChain->SetBranchAddress("phoSCBrem"		  ,&treeLeaf.phoSCBrem		   ,&b_phoSCBrem);		      
  fChain->SetBranchAddress("phohasPixelSeed"	  ,&treeLeaf.phohasPixelSeed	   ,&b_phohasPixelSeed);	      
  fChain->SetBranchAddress("phoEleVeto"		  ,&treeLeaf.phoEleVeto		   ,&b_phoEleVeto);		      
  fChain->SetBranchAddress("phoR9"		  ,&treeLeaf.phoR9		   ,&b_phoR9);		      
  fChain->SetBranchAddress("phoHoverE"		  ,&treeLeaf.phoHoverE		   ,&b_phoHoverE);		      
  fChain->SetBranchAddress("phoSigmaIEtaIEta"	  ,&treeLeaf.phoSigmaIEtaIEta	   ,&b_phoSigmaIEtaIEta);	      
  fChain->SetBranchAddress("phoSigmaIEtaIPhi"	  ,&treeLeaf.phoSigmaIEtaIPhi	   ,&b_phoSigmaIEtaIPhi);	      
  fChain->SetBranchAddress("phoSigmaIPhiIPhi"	  ,&treeLeaf.phoSigmaIPhiIPhi	   ,&b_phoSigmaIPhiIPhi);	      
  fChain->SetBranchAddress("phoESEffSigmaRR"	  ,&treeLeaf.phoESEffSigmaRR	   ,&b_phoESEffSigmaRR);	      
  fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5",&treeLeaf.phoSigmaIEtaIEtaFull5x5 ,&b_phoSigmaIEtaIEtaFull5x5);   
  fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5",&treeLeaf.phoSigmaIEtaIPhiFull5x5 ,&b_phoSigmaIEtaIPhiFull5x5);   
  fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5",&treeLeaf.phoSigmaIPhiIPhiFull5x5 ,&b_phoSigmaIPhiIPhiFull5x5);   
  fChain->SetBranchAddress("phoPFChIso"		  ,&treeLeaf.phoPFChIso		   ,&b_phoPFChIso);		      
  fChain->SetBranchAddress("phoPFChWorstIso"		  ,&treeLeaf.phoPFChWorstIso		   ,&b_phoPFChWorstIso);		      
  fChain->SetBranchAddress("phoPFPhoIso"		  ,&treeLeaf.phoPFPhoIso		   ,&b_phoPFPhoIso);	              
  fChain->SetBranchAddress("phoPFNeuIso"		  ,&treeLeaf.phoPFNeuIso		   ,&b_phoPFNeuIso);	              
  fChain->SetBranchAddress("phoFiredSingleTrgs"	  ,&treeLeaf.phoFiredSingleTrgs	   ,&b_phoFiredSingleTrgs);	      
  fChain->SetBranchAddress("phoFiredDoubleTrgs"	  ,&treeLeaf.phoFiredDoubleTrgs	   ,&b_phoFiredDoubleTrgs);  
  fChain->SetBranchAddress("phoIDbit"               ,&treeLeaf.phoIDbit                ,&b_phoIDbit);              
  fChain->SetBranchAddress("phoIDMVA"               ,&treeLeaf.phoIDMVA                ,&b_phoIDMVA);    
//jets (minimal stuff)
  fChain->SetBranchAddress("nJet",&treeLeaf.nJet,&b_nJet);
  fChain->SetBranchAddress("jetPt",&treeLeaf.jetPt,&b_jetPt);
  fChain->SetBranchAddress("jetEn",&treeLeaf.jetEn,&b_jetEn);
  fChain->SetBranchAddress("jetEta",&treeLeaf.jetEta,&b_jetEta);
  fChain->SetBranchAddress("jetPhi",	     &treeLeaf.jetPhi,	         &b_jetPhi);	      
  fChain->SetBranchAddress("jetRawPt",	     &treeLeaf.jetRawPt,	 &b_jetRawPt);	          
  fChain->SetBranchAddress("jetRawEn",	     &treeLeaf.jetRawEn,	 &b_jetRawEn);	          
  fChain->SetBranchAddress("jetMt",	     &treeLeaf.jetMt,	      	 &b_jetMt);	      
 fChain->SetBranchAddress("jetJECUnc",	     &treeLeaf.jetJECUnc,	 &b_jetJECUnc);	          
 fChain->SetBranchAddress("jetFiredTrgs",    &treeLeaf.jetFiredTrgs,    &b_jetFiredTrgs);     
 fChain->SetBranchAddress("jetGenJetIndex",  &treeLeaf.jetGenJetIndex,  &b_jetGenJetIndex);   
 fChain->SetBranchAddress("jetGenJetEn",     &treeLeaf.jetGenJetEn,     &b_jetGenJetEn);      
 fChain->SetBranchAddress("jetGenJetPt",     &treeLeaf.jetGenJetPt,     &b_jetGenJetPt);      
 fChain->SetBranchAddress("jetGenJetEta",    &treeLeaf.jetGenJetEta,    &b_jetGenJetEta);     
 fChain->SetBranchAddress("jetGenJetPhi",    &treeLeaf.jetGenJetPhi,    &b_jetGenJetPhi);     
 fChain->SetBranchAddress("jetGenPartonID",  &treeLeaf.jetGenPartonID,  &b_jetGenPartonID);   
 fChain->SetBranchAddress("jetGenEn",	     &treeLeaf.jetGenEn,	 &b_jetGenEn);	          
 fChain->SetBranchAddress("jetGenPt",	     &treeLeaf.jetGenPt,	 &b_jetGenPt);	          
 fChain->SetBranchAddress("jetGenEta",	     &treeLeaf.jetGenEta,	 &b_jetGenEta);	          
 fChain->SetBranchAddress("jetGenPhi",	     &treeLeaf.jetGenPhi,	 &b_jetGenPhi);	          
 fChain->SetBranchAddress("jetGenPartonMomID",&treeLeaf.jetGenPartonMomID,&b_jetGenPartonMomID);
        
 //met (minimal stuff)
 fChain->SetBranchAddress("metFilters",&treeLeaf.metFilters,&b_metFilters);
 fChain->SetBranchAddress("genMET",&treeLeaf.genMET,&b_genMET);
 fChain->SetBranchAddress("genMETPhi",&treeLeaf.genMETPhi,&b_genMETPhi);
 fChain->SetBranchAddress("pfMET",&treeLeaf.pfMET,&b_pfMET);
 fChain->SetBranchAddress("pfMETPhi",&treeLeaf.pfMETPhi,&b_pfMETPhi);
 fChain->SetBranchAddress("pfMETsumEt",&treeLeaf.pfMETsumEt,&b_pfMETsumEt);
 fChain->SetBranchAddress("pfMETmEtSig",&treeLeaf.pfMETmEtSig,&b_pfMETmEtSig);
 fChain->SetBranchAddress("pfMETSig",&treeLeaf.pfMETSig,&b_pfMETSig);

  // fChain->SetBranchAddress("eeMass"               ,&outLeaf.eeMass                ,&b_eeMass);                  
  // fChain->SetBranchAddress("eegMass"               ,&outLeaf.eegMass                ,&b_eegMass);               
}

void TInputOutputTree::InitOutput(TTree* outputTree,bool isMC)
{
  gROOT->ProcessLine("#include <vector>");
  //global event (keep all, as general rule)
  outputTree->Branch("nVtx",                 &treeLeaf.nVtx,       "nVtx/I");
  outputTree->Branch("run",                  &treeLeaf.run,        "run/I");
  outputTree->Branch("event",                &treeLeaf.event,      "event/L");
  outputTree->Branch("phoEt"		  ,&treeLeaf.phoEt	 );	
  outputTree->Branch("lumis",                &treeLeaf.lumis,      "lumis/I");
  outputTree->Branch("isData",               &treeLeaf.isData,     "isData/O");
  outputTree->Branch("nTrksPV",              &treeLeaf.nTrksPV,    "nVtxPV/I");
  outputTree->Branch("vtx",                  &treeLeaf.vtx,        "vtx/F"); 
  outputTree->Branch("vty",                  &treeLeaf.vty,        "vty/F"); 
  outputTree->Branch("vtz",                  &treeLeaf.vtz,        "vtz/F"); 
  outputTree->Branch("rho",                  &treeLeaf.rho,        "rho/F");
  outputTree->Branch("rhoCentral",           &treeLeaf.rhoCentral, "rhoCentral/F");
  outputTree->Branch("HLTEleMuX",            &treeLeaf.HLTEleMuX, "HLTEleMuX/l");
  outputTree->Branch("HLTPho",               &treeLeaf.HLTPho, "HLTPho/l");
  outputTree->Branch("HLTJet",               &treeLeaf.HLTJet, "HLTJet/l");
  outputTree->Branch("HLTEleMuXIsPrescaled", &treeLeaf.HLTEleMuXIsPrescaled, "HLTEleMuXIsPrescaled/l");
  outputTree->Branch("HLTPhoIsPrescaled",    &treeLeaf.HLTPhoIsPrescaled, "HLTPhoIsPrescaled/l");
  outputTree->Branch("HLTJetIsPrescaled",    &treeLeaf.HLTJetIsPrescaled, "HLTJetIsPrescaled/l");
  //gen info
  // (local) variables associated with tree branches
  if(isMC){
    outputTree->Branch("pdf",&treeLeaf.pdf);
    outputTree->Branch("pthat",&treeLeaf.pthat);
    outputTree->Branch("processID",&treeLeaf.processID);
    outputTree->Branch("genWeight",&treeLeaf.genWeight);
    outputTree->Branch("genHT",&treeLeaf.genHT);

    outputTree->Branch("nPUInfo",&treeLeaf.nPUInfo);
    outputTree->Branch("nPU",&treeLeaf.nPU);
    outputTree->Branch("puBX",&treeLeaf.puBX);
    outputTree->Branch("puTrue",&treeLeaf.puTrue);

    outputTree->Branch("nMC",         &treeLeaf.nMC);
    outputTree->Branch("mcPID",       &treeLeaf.mcPID);
    outputTree->Branch("mcVtx",	   &treeLeaf.mcVtx);	    
    outputTree->Branch("mcVty",	   &treeLeaf.mcVty);	    
    outputTree->Branch("mcVtz",	   &treeLeaf.mcVtz);	    
    outputTree->Branch("mcPt",	   &treeLeaf.mcPt);	    
    outputTree->Branch("mcMass",	   &treeLeaf.mcMass);	    
    outputTree->Branch("mcEta",	   &treeLeaf.mcEta);	    
    outputTree->Branch("mcPhi",	   &treeLeaf.mcPhi);	    
    outputTree->Branch("mcE",	   &treeLeaf.mcE);	    
    outputTree->Branch("mcEt",	   &treeLeaf.mcEt);	    
    outputTree->Branch("mcGMomPID",  &treeLeaf.mcGMomPID);    
    outputTree->Branch("mcMomPID",   &treeLeaf.mcMomPID);     
    outputTree->Branch("mcMomPt",	   &treeLeaf.mcMomPt);		    
    outputTree->Branch("mcMomMass",  &treeLeaf.mcMomMass);    
    outputTree->Branch("mcMomEta",   &treeLeaf.mcMomEta);     
    outputTree->Branch("mcMomPhi",   &treeLeaf.mcMomPhi);     
    outputTree->Branch("mcIndex",	   &treeLeaf.mcIndex);		    
    outputTree->Branch("mcStatusFlag",&treeLeaf.mcStatusFlag);
    outputTree->Branch("mcParentage",&treeLeaf.mcParentage);
    outputTree->Branch("mcStatus",   &treeLeaf.mcStatus);
    outputTree->Branch("mcCalIsoDR03",&treeLeaf.mcCalIsoDR03);
    outputTree->Branch("mcTrkIsoDR03",&treeLeaf.mcTrkIsoDR03);
    outputTree->Branch("mcCalIsoDR04",&treeLeaf.mcCalIsoDR04);
    outputTree->Branch("mcTrkIsoDR04",&treeLeaf.mcTrkIsoDR04);
  }
  //muons
  outputTree->Branch("nMu",&treeLeaf.nMu,"nMu/I");
  outputTree->Branch("muPt",&treeLeaf.muPt);	     
  outputTree->Branch("muEn",&treeLeaf.muEn);	     
  outputTree->Branch("muEta",&treeLeaf.muEta);	     
  outputTree->Branch("muPhi",&treeLeaf.muPhi);	     
  outputTree->Branch("muCharge",&treeLeaf.muCharge);   
  outputTree->Branch("muType",&treeLeaf.muType   );	     
  outputTree->Branch("muIsLooseID",&treeLeaf.muIsLooseID); 
  outputTree->Branch("muIsMediumID",&treeLeaf.muIsMediumID);
  outputTree->Branch("muIsTightID",&treeLeaf.muIsTightID);
  outputTree->Branch("muIsSoftID",&treeLeaf.muIsSoftID);
  outputTree->Branch("muIsHighPtID",&treeLeaf.muIsHighPtID);
  outputTree->Branch("muD0",&treeLeaf.muD0	   );	     
  outputTree->Branch("muDz",&treeLeaf.muDz	   );	     
  outputTree->Branch("muChi2ND ",&treeLeaf.muChi2NDF );   
  outputTree->Branch("muTrkLayers",&treeLeaf.muTrkLayers); 
  outputTree->Branch("muPixelHits",&treeLeaf.muPixelHits); 
  outputTree->Branch("muStations   ",&treeLeaf.muStations);  
  outputTree->Branch("muPFChIso    ",&treeLeaf.muPFChIso );   
  outputTree->Branch("muPFPhoIso   ",&treeLeaf.muPFPhoIso);  
  outputTree->Branch("muPFNeuIso   ",&treeLeaf.muPFNeuIso); 
  outputTree->Branch("muPFPUIso    ",&treeLeaf.muPFPUIso );   
  outputTree->Branch("muFiredTrgs  ",&treeLeaf.muFiredTrgs);
  // electrons
  outputTree->Branch("nEle",                    &treeLeaf.nEle,                    "nEle/I");
  outputTree->Branch("eleCharge",               &treeLeaf.eleCharge);
  outputTree->Branch("eleChargeConsistent",     &treeLeaf.eleChargeConsistent);
  outputTree->Branch("eleEn",                   &treeLeaf.eleEn);
  outputTree->Branch("eleSCEn",                 &treeLeaf.eleSCEn);
  // outputTree->Branch("eleESEn",                 &treeLeaf.eleESEn);
  // outputTree->Branch("eleESEnP1",               &treeLeaf.eleESEnP1);
  // outputTree->Branch("eleESEnP2",               &treeLeaf.eleESEnP2);
  outputTree->Branch("eleD0",                   &treeLeaf.eleD0);
  outputTree->Branch("eleDz",                   &treeLeaf.eleDz);
  outputTree->Branch("elePt",                   &treeLeaf.elePt);
  outputTree->Branch("eleEta",                  &treeLeaf.eleEta);
  outputTree->Branch("elePhi",                  &treeLeaf.elePhi);
  outputTree->Branch("eleR9",                   &treeLeaf.eleR9);
  outputTree->Branch("eleSCEta",                &treeLeaf.eleSCEta);
  outputTree->Branch("eleSCPhi",                &treeLeaf.eleSCPhi);
  outputTree->Branch("eleSCRawEn",              &treeLeaf.eleSCRawEn);
  // outputTree->Branch("eleSCEtaWidth",           &treeLeaf.eleSCEtaWidth);
  // outputTree->Branch("eleSCPhiWidth",           &treeLeaf.eleSCPhiWidth);
  outputTree->Branch("eleHoverE",               &treeLeaf.eleHoverE);
  outputTree->Branch("eleEoverP",               &treeLeaf.eleEoverP);
  // outputTree->Branch("eleEoverPout",            &treeLeaf.eleEoverPout);
  // outputTree->Branch("eleEoverPInv",            &treeLeaf.eleEoverPInv);
  outputTree->Branch("eleBrem",                 &treeLeaf.eleBrem);
  outputTree->Branch("eledEtaAtVtx",            &treeLeaf.eledEtaAtVtx);
  outputTree->Branch("eledPhiAtVtx",            &treeLeaf.eledPhiAtVtx);
  // outputTree->Branch("eledEtaAtCalo",           &treeLeaf.eledEtaAtCalo);
  outputTree->Branch("eleSigmaIEtaIEta",        &treeLeaf.eleSigmaIEtaIEta);
  outputTree->Branch("eleSigmaIEtaIPhi",        &treeLeaf.eleSigmaIEtaIPhi);
  outputTree->Branch("eleSigmaIPhiIPhi",        &treeLeaf.eleSigmaIPhiIPhi);
  outputTree->Branch("eleSigmaIEtaIEtaFull5x5", &treeLeaf.eleSigmaIEtaIEtaFull5x5);
  outputTree->Branch("eleSigmaIPhiIPhiFull5x5", &treeLeaf.eleSigmaIPhiIPhiFull5x5);
  outputTree->Branch("eleConvVeto",             &treeLeaf.eleConvVeto);
  outputTree->Branch("eleMissHits",             &treeLeaf.eleMissHits);
  outputTree->Branch("eleESEffSigmaRR",         &treeLeaf.eleESEffSigmaRR);
  outputTree->Branch("elePFChIso",              &treeLeaf.elePFChIso);
  outputTree->Branch("elePFPhoIso",             &treeLeaf.elePFPhoIso);
  outputTree->Branch("elePFNeuIso",             &treeLeaf.elePFNeuIso);
  outputTree->Branch("elePFPUIso",              &treeLeaf.elePFPUIso);
  outputTree->Branch("elePFClusEcalIso",        &treeLeaf.elePFClusEcalIso);
  outputTree->Branch("elePFClusHcalIso",        &treeLeaf.elePFClusHcalIso);
  // outputTree->Branch("eleIDMVANonTrg",          &treeLeaf.eleIDMVANonTrg);
  outputTree->Branch("eleIDMVATrg",             &treeLeaf.eleIDMVATrg);
  // outputTree->Branch("eledEtaseedAtVtx",        &treeLeaf.eledEtaseedAtVtx);
  // outputTree->Branch("eleE1x5",                 &treeLeaf.eleE1x5);
  // outputTree->Branch("eleE2x5",                 &treeLeaf.eleE2x5);
  // outputTree->Branch("eleE5x5",                 &treeLeaf.eleE5x5);
  // outputTree->Branch("eleE1x5Full5x5",          &treeLeaf.eleE1x5Full5x5);
  // outputTree->Branch("eleE2x5Full5x5",          &treeLeaf.eleE2x5Full5x5);
  // outputTree->Branch("eleE5x5Full5x5",          &treeLeaf.eleE5x5Full5x5);
  // outputTree->Branch("eleR9Full5x5",                &treeLeaf.eleR9Full5x5);
  // outputTree->Branch("eleEcalDrivenSeed",           &treeLeaf.eleEcalDrivenSeed);
  // outputTree->Branch("eleDr03EcalRecHitSumEt",      &treeLeaf.eleDr03EcalRecHitSumEt);
  // outputTree->Branch("eleDr03HcalDepth1TowerSumEt", &treeLeaf.eleDr03HcalDepth1TowerSumEt);
  // outputTree->Branch("eleDr03HcalDepth2TowerSumEt", &treeLeaf.eleDr03HcalDepth2TowerSumEt);
  // outputTree->Branch("eleDr03HcalTowerSumEt",       &treeLeaf.eleDr03HcalTowerSumEt);
  outputTree->Branch("eleDr03TkSumPt",              &treeLeaf.eleDr03TkSumPt);
  // outputTree->Branch("elecaloEnergy",               &treeLeaf.elecaloEnergy);
  // outputTree->Branch("eleTrkdxy",                   &treeLeaf.eleTrkdxy);
  // outputTree->Branch("eleKFHits",                   &treeLeaf.eleKFHits);
  // outputTree->Branch("eleKFChi2",                   &treeLeaf.eleKFChi2);
  // outputTree->Branch("eleGSFPt",                    &treeLeaf.eleGSFPt);
  // outputTree->Branch("eleGSFEta",                   &treeLeaf.eleGSFEta);
  // outputTree->Branch("eleGSFPhi",                   &treeLeaf.eleGSFPhi);
  // outputTree->Branch("eleGSFCharge",                &treeLeaf.eleGSFCharge);
  // outputTree->Branch("eleGSFHits",                  &treeLeaf.eleGSFHits);
  // outputTree->Branch("eleGSFMissHits",              &treeLeaf.eleGSFMissHits);
  // outputTree->Branch("eleGSFNHitsMax",              &treeLeaf.eleGSFNHitsMax);
  // outputTree->Branch("eleGSFVtxProb",               &treeLeaf.eleGSFVtxProb);
  // outputTree->Branch("eleGSFlxyPV",                 &treeLeaf.eleGSFlxyPV);
  // outputTree->Branch("eleGSFlxyBS",                 &treeLeaf.eleGSFlxyBS);
  // outputTree->Branch("eleBCEn",                     &treeLeaf.eleBCEn);
  // outputTree->Branch("eleBCEta",                    &treeLeaf.eleBCEta);
  // outputTree->Branch("eleBCPhi",                    &treeLeaf.eleBCPhi);
  // outputTree->Branch("eleBCS25",                    &treeLeaf.eleBCS25);
  // outputTree->Branch("eleBCS15",                    &treeLeaf.eleBCS15);
  // outputTree->Branch("eleBCSieie",                  &treeLeaf.eleBCSieie);
  // outputTree->Branch("eleBCSieip",                  &treeLeaf.eleBCSieip);
  // outputTree->Branch("eleBCSipip",                  &treeLeaf.eleBCSipip);
  outputTree->Branch("eleFiredTrgs",                &treeLeaf.eleFiredTrgs);
  //  if (runeleIDVID)
  outputTree->Branch("eleIDbit",  &treeLeaf.eleIDbit);
  //devel only
  // outputTree->Branch("eleESEnP1Raw",              &treeLeaf.eleESEnP1Raw);
  // outputTree->Branch("eleESEnP2Raw",              &treeLeaf.eleESEnP2Raw);
  // outputTree->Branch("eleESEnEta",                &treeLeaf.eleESEnEta);
  // outputTree->Branch("eleESEnPhi",                &treeLeaf.eleESEnPhi);
  // outputTree->Branch("eleESEnZ",                  &treeLeaf.eleESEnZ);
  // outputTree->Branch("eleESEnP",                  &treeLeaf.eleESEnP);
  // outputTree->Branch("eleESEnX",                  &treeLeaf.eleESEnX);
  // outputTree->Branch("eleESEnY",                  &treeLeaf.eleESEnY);
  // outputTree->Branch("eleESEnS",                  &treeLeaf.eleESEnS);
  // outputTree->Branch("eleESEnE",                  &treeLeaf.eleESEnE);
  // outputTree->Branch("nGSFTrk",                   &treeLeaf.nGSFTrk);
  // outputTree->Branch("gsfPt",                     &treeLeaf.gsfPt);
  // outputTree->Branch("gsfEta",                    &treeLeaf.gsfEta);
  // outputTree->Branch("gsfPhi",                    &treeLeaf.gsfPhi);

  //photon
  outputTree->Branch("nPho"		  ,&treeLeaf.nPho	 );		     
  outputTree->Branch("phoE"		  ,&treeLeaf.phoE      	);		     
  outputTree->Branch("phoEta"		  ,&treeLeaf.phoEta	 );		     
  outputTree->Branch("phoPhi"		  ,&treeLeaf.phoPhi	 );		     
  outputTree->Branch("phoSCE"		  ,&treeLeaf.phoSCE	);		      
  outputTree->Branch("phoSCRawE"		  ,&treeLeaf.phoSCRawE	);		      
  outputTree->Branch("phoSCEta"		  ,&treeLeaf.phoSCEta		   );		      
  outputTree->Branch("phoSCPhi"		  ,&treeLeaf.phoSCPhi		   );		      
  outputTree->Branch("phoSCEtaWidth"	  ,&treeLeaf.phoSCEtaWidth	   );	      
  outputTree->Branch("phoSCPhiWidth"	  ,&treeLeaf.phoSCPhiWidth	   );	      
  outputTree->Branch("phoSCBrem"		  ,&treeLeaf.phoSCBrem	   );		      
  outputTree->Branch("phohasPixelSeed"	  ,&treeLeaf.phohasPixelSeed	  );	      
  outputTree->Branch("phoEleVeto"		  ,&treeLeaf.phoEleVeto	);		      
  outputTree->Branch("phoR9"		  ,&treeLeaf.phoR9		   );		      
  outputTree->Branch("phoHoverE"		  ,&treeLeaf.phoHoverE	);		      
  outputTree->Branch("phoSigmaIEtaIEta"	  ,&treeLeaf.phoSigmaIEtaIEta);	      
  outputTree->Branch("phoSigmaIEtaIPhi"	  ,&treeLeaf.phoSigmaIEtaIPhi);	      
  outputTree->Branch("phoSigmaIPhiIPhi"	  ,&treeLeaf.phoSigmaIPhiIPhi);	      
  outputTree->Branch("phoESEffSigmaRR"	  ,&treeLeaf.phoESEffSigmaRR	   );	      
  outputTree->Branch("phoSigmaIEtaIEtaFull5x5",&treeLeaf.phoSigmaIEtaIEtaFull5x5 );   
  outputTree->Branch("phoSigmaIEtaIPhiFull5x5",&treeLeaf.phoSigmaIEtaIPhiFull5x5 );   
  outputTree->Branch("phoSigmaIPhiIPhiFull5x5",&treeLeaf.phoSigmaIPhiIPhiFull5x5 );   
  outputTree->Branch("phoPFChIso"		  ,&treeLeaf.phoPFChIso		 );		      
  outputTree->Branch("phoPFChWorstIso"		  ,&treeLeaf.phoPFChWorstIso		 );		      
  outputTree->Branch("phoPFPhoIso"	  ,&treeLeaf.phoPFPhoIso		 );	              
  outputTree->Branch("phoPFNeuIso"	  ,&treeLeaf.phoPFNeuIso		 );	              
  outputTree->Branch("phoFiredSingleTrgs"	  ,&treeLeaf.phoFiredSingleTrgs	 );	      
  outputTree->Branch("phoFiredDoubleTrgs"	  ,&treeLeaf.phoFiredDoubleTrgs	 );  
  outputTree->Branch("phoIDbit"             ,&treeLeaf.phoIDbit              );    
  outputTree->Branch("phoIDMVA"             ,&treeLeaf.phoIDMVA              );    
  //jets (minimal stuff)
  outputTree->Branch("nJet",&treeLeaf.nJet, "nJet/I");
  outputTree->Branch("jetPt",&treeLeaf.jetPt);
  outputTree->Branch("jetEn",&treeLeaf.jetEn);
  outputTree->Branch("jetEta",&treeLeaf.jetEta);
  outputTree->Branch("jetPhi",	     &treeLeaf.jetPhi);
  outputTree->Branch("jetRawPt",	     &treeLeaf.jetRawPt);	          
  outputTree->Branch("jetRawEn",	     &treeLeaf.jetRawEn); 
  outputTree->Branch("jetMt",	     &treeLeaf.jetMt);	      	      
  outputTree->Branch("jetJECUnc",	     &treeLeaf.jetJECUnc);	          
  outputTree->Branch("jetFiredTrgs",    &treeLeaf.jetFiredTrgs);    
  outputTree->Branch("jetGenJetIndex",  &treeLeaf.jetGenJetIndex);  
  outputTree->Branch("jetGenJetEn",     &treeLeaf.jetGenJetEn);     
  outputTree->Branch("jetGenJetPt",     &treeLeaf.jetGenJetPt);     
  outputTree->Branch("jetGenJetEta",    &treeLeaf.jetGenJetEta);    
  outputTree->Branch("jetGenJetPhi",    &treeLeaf.jetGenJetPhi);    
  outputTree->Branch("jetGenPartonID",  &treeLeaf.jetGenPartonID);  
  outputTree->Branch("jetGenEn",	     &treeLeaf.jetGenEn);            
  outputTree->Branch("jetGenPt",	     &treeLeaf.jetGenPt);            
  outputTree->Branch("jetGenEta",	     &treeLeaf.jetGenEta);           
  outputTree->Branch("jetGenPhi",	     &treeLeaf.jetGenPhi);           
  outputTree->Branch("jetGenPartonMomID",&treeLeaf.jetGenPartonMomID);

  //met (minimal stuff)
  outputTree->Branch("metFilters",&treeLeaf.metFilters,"metFilters/I");
  outputTree->Branch("genMET",&treeLeaf.genMET,"genMET/F");
  outputTree->Branch("genMETPhi",&treeLeaf.genMETPhi,"genMETPhi/F");
  outputTree->Branch("pfMET",&treeLeaf.pfMET,"pfMET/F");
  outputTree->Branch("pfMETPhi",&treeLeaf.pfMETPhi,"pfMETPhi/F");
  outputTree->Branch("pfMETsumEt",&treeLeaf.pfMETsumEt,"pfMETsumEt/F");
  outputTree->Branch("pfMETmEtSig",&treeLeaf.pfMETmEtSig,"pfMETmEtSig");
  outputTree->Branch("pfMETSig",&treeLeaf.pfMETSig,"pfMETSig");
        
  //  outputTree->Branch("eeMass"      ,&outLeaf.eeMass);
  // outputTree->Branch("eegMass"     ,&outLeaf.eegMass);
  // outputTree->Branch("gammaEt"     ,&outLeaf.gammaEt);
  // outputTree->Branch("deltaR"     ,&outLeaf.deltaR);
}
