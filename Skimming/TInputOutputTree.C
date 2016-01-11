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
  /*  treeLeaf.eledEtaAtVtx=0; */
  /*  treeLeaf.eledPhiAtVtx=0; */
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
  /*  treeLeaf.eleIDMVATrg=0; */
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
  /*  treeLeaf.eleDr03TkSumPt=0; */
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
}

void TInputOutputTree::InitOutput(TTree* outputTree)
{
  //global event (keep all, as general rule)
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
  // electrons (only kinematics, trigger, for now)
  outputTree->Branch("nEle",                    &treeLeaf.nEle,                    "nEle/I");
  outputTree->Branch("eleCharge",               &treeLeaf.eleCharge,               "eleCharge/I");
  outputTree->Branch("eleChargeConsistent",     &treeLeaf.eleChargeConsistent,     "eleChargeConsistent/I");
  outputTree->Branch("eleEn",                   &treeLeaf.eleEn,                   "eleEn/F");
  outputTree->Branch("eleSCEn",                 &treeLeaf.eleSCEn,                 "eleSCEn/F");
  // outputTree->Branch("eleESEn",                 &treeLeaf.eleESEn);
  // outputTree->Branch("eleESEnP1",               &treeLeaf.eleESEnP1);
  // outputTree->Branch("eleESEnP2",               &treeLeaf.eleESEnP2);
  outputTree->Branch("eleD0",                   &treeLeaf.eleD0, "eleD0/F");
  outputTree->Branch("eleDz",                   &treeLeaf.eleDz, "eleDz/F");
  outputTree->Branch("elePt",                   &treeLeaf.elePt, "elePt/F");
  outputTree->Branch("eleEta",                  &treeLeaf.eleEta, "eleEta/F");
  outputTree->Branch("elePhi",                  &treeLeaf.elePhi, "elePhi/F");
  outputTree->Branch("eleR9",                   &treeLeaf.eleR9, "eleR9/F");
  outputTree->Branch("eleSCEta",                &treeLeaf.eleSCEta, "eleSCEta/F");
  outputTree->Branch("eleSCPhi",                &treeLeaf.eleSCPhi, "eleSCPhi/F");
  outputTree->Branch("eleSCRawEn",              &treeLeaf.eleSCRawEn, "eleSCRawEn/F");
  // outputTree->Branch("eleSCEtaWidth",           &treeLeaf.eleSCEtaWidth);
  // outputTree->Branch("eleSCPhiWidth",           &treeLeaf.eleSCPhiWidth);
  outputTree->Branch("eleHoverE",               &treeLeaf.eleHoverE, "eleHoverE/F");
  outputTree->Branch("eleEoverP",               &treeLeaf.eleEoverP, "eleEoverP/F");
  // outputTree->Branch("eleEoverPout",            &treeLeaf.eleEoverPout);
  // outputTree->Branch("eleEoverPInv",            &treeLeaf.eleEoverPInv);
  outputTree->Branch("eleBrem",                 &treeLeaf.eleBrem, "eleBrem/F");
  // outputTree->Branch("eledEtaAtVtx",            &treeLeaf.eledEtaAtVtx);
  // outputTree->Branch("eledPhiAtVtx",            &treeLeaf.eledPhiAtVtx);
  // outputTree->Branch("eledEtaAtCalo",           &treeLeaf.eledEtaAtCalo);
  outputTree->Branch("eleSigmaIEtaIEta",        &treeLeaf.eleSigmaIEtaIEta, "eleSigmaIEtaIEta/F");
  outputTree->Branch("eleSigmaIEtaIPhi",        &treeLeaf.eleSigmaIEtaIPhi, "eleSigmaIEtaIPhi/F");
  outputTree->Branch("eleSigmaIPhiIPhi",        &treeLeaf.eleSigmaIPhiIPhi, "eleSigmaIPhiIPhi/F");
  outputTree->Branch("eleSigmaIEtaIEtaFull5x5", &treeLeaf.eleSigmaIEtaIEtaFull5x5, "eleSigmaIEtaIEtaFull5x5/F");
  outputTree->Branch("eleSigmaIPhiIPhiFull5x5", &treeLeaf.eleSigmaIPhiIPhiFull5x5, "eleSigmaIPhiIPhi5x5/F");
  outputTree->Branch("eleConvVeto",             &treeLeaf.eleConvVeto, "eleConvVeto/I");
  outputTree->Branch("eleMissHits",             &treeLeaf.eleMissHits, "eleMissHits/I");
  outputTree->Branch("eleESEffSigmaRR",         &treeLeaf.eleESEffSigmaRR, "eleESEffSigmaRR/F");
  outputTree->Branch("elePFChIso",              &treeLeaf.elePFChIso, "elePFChIso/F");
  outputTree->Branch("elePFPhoIso",             &treeLeaf.elePFPhoIso, "elePFPhoIso/F");
  outputTree->Branch("elePFNeuIso",             &treeLeaf.elePFNeuIso, "elePFNeuIso/F");
  outputTree->Branch("elePFPUIso",              &treeLeaf.elePFPUIso, "elePFPUIso/F");
  outputTree->Branch("elePFClusEcalIso",        &treeLeaf.elePFClusEcalIso, "elePFClusEcalIso/F");
  outputTree->Branch("elePFClusHcalIso",        &treeLeaf.elePFClusHcalIso, "elePFClusHCalIso/F");
  // outputTree->Branch("eleIDMVANonTrg",          &treeLeaf.eleIDMVANonTrg);
  // outputTree->Branch("eleIDMVATrg",             &treeLeaf.eleIDMVATrg);
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
  // outputTree->Branch("eleDr03TkSumPt",              &treeLeaf.eleDr03TkSumPt);
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
  outputTree->Branch("eleFiredTrgs",                &treeLeaf.eleFiredTrgs, "eleFiredTrgs/I");

  //  if (runeleIDVID)
  outputTree->Branch("eleIDbit",  &treeLeaf.eleIDbit, "eleIDbit/s");

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
}
