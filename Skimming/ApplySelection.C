#include <TString.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
using namespace ROOT;
using namespace std;


float deltaR(float eta0, float eta1, float phi0, float phi1){
  return sqrt((eta0-eta1)*(eta0-eta1)+(phi0-phi1)*(phi0-phi1));
}

void ApplySelection()
{
 //this time, apply it for real
   TString _inputFileName("RD_EleElePhoSkim_big.root");
   // TString _inputFileName("mc_DYJetsToLL_M10-50_LLGSkim_big.root");
  TFile *f = TFile::Open(_inputFileName,"READ");
  std::cout<<"processing "<<_inputFileName<<std::endl;
  TTree* tree =(TTree*)f->Get("EventTree");
  

  Long64_t nentries = tree->GetEntries();
  std::cout<<"nentries "<<nentries<<std::endl;
  
  TH1F* hPhotons = new TH1F("hPhotons","number of photon candidates",20,0,20);
  TH1F* hPID = new TH1F("hPID","number of photon passing ID cuts",20,0,20);
  TH1F* hPhotonPt = new TH1F("hPhotonPt","p_{T}^{#gamma}",500,0,1000);
  TH1F* hDeltaR = new TH1F("hDeltaR","#DeltaR(l,#gamma)",28,0.7,6.1);
  TH1F* hMll = new TH1F("hMll","M_{ll} [GeV/c^{2}]",500,0,1000);
  TH1F* hMllg = new TH1F("hMllg","M_{ll#gamma} [GeV/c^{2}]",500,0,1000);

  TH1F* hElectrons = new TH1F("hElectrons","number of electron candidates",20,0,20); //...
 
  //some leptons...
  const float _eMass = 0.000510998910; // GeV/c2
  const float _muMass = 0.1056583715; // GeV/c2
  const float _tauMass = 1.77682; // GeV/c2
  //and a photon...for safety
  const float _gMass = 0.; 
  TLorentzVector lepton[1000]={};
  TLorentzVector photon[1000]={};
  TLorentzVector ll[1000]={}; // Z candidates
  TLorentzVector llg[1000]={}; // Zg candidates

  float _mllg[100]={};
  float _mll[100]={};
  float _phoEt[100]={};
  float _dr1;
  float _dr2;
  //  float _eplus,_eminus;
  // TBranch *b_mllg;
  // TBranch *b_mll;
  // TBranch *b_phoEt;
  // TBranch *b_dr;
  nentries = 1570000;

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
  bool goodElectron[100]={};
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
  tree->SetBranchAddress("eleHoverE"	    ,&eleHoverE);	     
  tree->SetBranchAddress("eleSigmaIEtaIEtaFull5x5",&eleSigmaIEtaIEtaFull5x5);
  //photon variables
  int nPho;
  float _pPt=0;
  float _pPhi=0;
  float _pEta=0;
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
  tree->SetBranchAddress("nPho",&nPho);
  tree->SetBranchAddress("phoEt",&phoEt);
  tree->SetBranchAddress("phoPhi",&phoPhi);
  tree->SetBranchAddress("phoEta",&phoEta);
  tree->SetBranchAddress("phoSCEta",&phoSCEta);
  tree->SetBranchAddress("phoIDMVA",&phoIDMVA);	
  tree->SetBranchAddress("phoSigmaIEtaIEtaFull5x5",&phoSigmaIEtaIEtaFull5x5);
  tree->SetBranchAddress("phoHoverE",&phoHoverE);
  tree->SetBranchAddress("phoPFPhoIso",&phoPFPhoIso);
  tree->SetBranchAddress("phoPFChWorstIso",&phoPFChWorstIso);
  tree->SetBranchAddress("phoR9",&phoR9);

  //make reduced tree (only passing the cuts)
  TFile *zgf = new TFile("red_eeg_"+_inputFileName,"RECREATE");
  TTree *_zgt = new TTree("zgtree","RECREATE");
  _zgt->Branch("nVtx",                 &nVtx,       "nVtx/I");
  _zgt->Branch("run",                  &run,        "run/I");
  _zgt->Branch("event",                &event,      "event/L");
  //  _zgt->Branch("phoEt",                &phoEt	 );  
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
  _zgt->Branch("eleHoverE"	    ,&eleHoverE);	     
  _zgt->Branch("eleSigmaIEtaIEtaFull5x5",&eleSigmaIEtaIEtaFull5x5);

  _zgt->Branch("_mllg",&_mllg,"_mllg[100]/F");
  //  _zgt->Branch("_eplus",&_mllg);
  _zgt->Branch("_mll",&_mll,"_mll[100]/F");
  _zgt->Branch("_phoEt",&_phoEt,"_phoEt[nPho]/F");
  _zgt->Branch("_dr1",&_dr1,"_dr1[nEle]/F");   
  _zgt->Branch("_dr2",&_dr2,"_dr2[nEle]/F");   





//   _zgt->Branch("nPho",&nPho);
//   _zgt->Branch("phoEt",&phoEt);
//   _zgt->Branch("phoPhi",&phoPhi);
//   _zgt->Branch("phoEta",&phoEta);
//   _zgt->Branch("phoSCEta",&phoSCEta);
//   _zgt->Branch("phoIDMVA",&phoIDMVA);	
//   _zgt->Branch("phoSigmaIEtaIEtaFull5x5",&phoSigmaIEtaIEtaFull5x5);
//   _zgt->Branch("phoHoverE",&phoHoverE);
//   _zgt->Branch("phoPFPhoIso",&phoPFPhoIso);
//   _zgt->Branch("phoPFChWorstIso",&phoPFChWorstIso);
//   _zgt->Branch("phoR9",&phoR9);
  //  _zgt->Branch();

  for (int entry=0; entry<nentries;entry++)
    {
      //      std::cout<<"hey"<<std::endl;
      int iz=0; //Z counter.
      int izg=0; //Zg counter.
      tree->GetEntry(entry);
      //      std::cout<<nEle<<std::endl;
      hElectrons->Fill(nEle);
      if(elePt->at(0)<20){ continue;}
      if(elePt->at(0) > 20.0 && elePt->at(1) > 15.0)
	{	
	  //make electron pairs
	  //      std::cout<< elePt<<std::endl;
	  for (int ie=0; ie<nEle;ie++)
	    {
	      goodElectron[ie]=false;
	      if((HLTEleMuX>>7&1)==1 && ((ie==0 && elePt->at(ie) > 20.0) || (ie>0 && elePt->at(ie) > 15.0)))
	  	{ 
	  	  // MVA Selection
	  	  if ((abs(eleSCEta->at(ie))<0.8 && ( eleIDMVATrg->at(ie)>0.972153 ))
	  	      ||((abs(eleSCEta->at(ie))>0.8 && abs(eleSCEta->at(ie))<1.479) && eleIDMVATrg->at(ie)>0.922126)
	  	      ||( (abs(eleSCEta->at(ie))>1.479 && abs(eleSCEta->at(ie))<2.5) && eleIDMVATrg->at(ie)>0.610764 ))
	  	    {
		      // std::cout<<"ok!"<<std::endl;
	  	      //apply everything else once (electron cuts)
	  	      //eta //(barrel:track matching) //relative trk iso //relative ECAL cluster iso //relative HCAL cluster iso //H over E //Sigma-Ieta-Ieta
	  	      if((abs(eleSCEta->at(ie))<1.442 && eledEtaAtVtx->at(ie)<0.0095 && eledPhiAtVtx->at(ie)<0.065 && eleDr03TkSumPt->at(ie)<0.18  && elePFClusEcalIso->at(ie)<0.37  && elePFClusHcalIso->at(ie)<0.25  && eleHoverE->at(ie)<0.09 && eleSigmaIEtaIEtaFull5x5->at(ie)<0.012)
	  	  	 ||(abs(eleSCEta->at(ie))>1.556 && eleDr03TkSumPt->at(ie)<0.18  && elePFClusEcalIso->at(ie)<0.45  && elePFClusHcalIso->at(ie)<0.28 && eleHoverE->at(ie)<0.09 && eleSigmaIEtaIEtaFull5x5->at(ie)<0.033
		  	    ))
	  	  	{
	  	  	  _ePt[ie] = elePt->at(ie);
	  	  	  _eEta[ie] = eleSCEta->at(ie);
	  	  	  _ePhi[ie] = elePhi->at(ie);
	  	  	  lepton[ie].SetPtEtaPhiM(_ePt[ie],_eEta[ie],_ePhi[ie],_eMass); //make a lepton TLorentzVector everytime there's a well-ID'd lepton
	  	  	  goodElectron[ie]=true;
	  	  	  //	  if(goodElectron[ie]==true) 		  std::cout<<entry<< "th entry, "<<" pt("<<ie<<") = "<<elePt->at(ie) << std::endl;
	  	  	}
	  	    }
	  	}
	    }
	  iz =0; //init Z and Zg counter now (before making the candidates)
	  izg =0; 
	  int pid=0; 
	  bool goodZg;
	  float llgm[100]={};
	  // for filling the ID cut histo
	  for(int ie1=0; ie1<nEle;ie1++){
	    for(int ie2 =ie1+1 ; ie2<nEle;ie2++)
	      {
		//		ll[iz].Clear();
		if(((_ePt[ie2]>20 && _ePt[ie1]>15) || (_ePt[ie1]>20 && _ePt[ie2]>15)) && (HLTEleMuX>>7&1)==1 && goodElectron[ie1]==true && goodElectron[ie2]==true)
		  {
		    ll[iz]= lepton[ie1]+lepton[ie2];
		    if(ll[iz].M() > 50.0 ){
		      int passingPhotons=0;
		      //start the photon loop here
		      for (int ip =0; ip<nPho; ip++)
			{
			  //kinematics (supercluster eta<2.5, avoiding gap region, and pt>15)
			  if(phoEt->at(ip) < 15 || ((abs(phoSCEta->at(ip))>1.442 && abs(phoSCEta->at(ip))<1.556) || abs(phoSCEta->at(ip))>2.5)) { pid=1; hPID->Fill(pid); continue;}
			  pid=2; hPID->Fill(pid);
			  float dr1=0, dr2=0;
			  photon[ip].Clear();
			  dr1=deltaR(_eEta[ie1],phoSCEta->at(ip),_ePhi[ie1],phoPhi->at(ip));
			  dr2=deltaR(_eEta[ie2],phoSCEta->at(ip),_ePhi[ie2],phoPhi->at(ip));
			  //			  std::cout<<"nentries "<<nentries<<std::endl;
			  //DeltaR cut
			  if (dr1 > 0.7 && dr2 > 0.7){
			    pid=3; hPID->Fill(pid);
			    //R9 (only cutting on endcaps)
			      if((abs(phoSCEta->at(ip))>1.556 &&  phoR9->at(ip) > 0.85  )|| abs(phoSCEta->at(ip))<1.442){
				pid=4; hPID->Fill(pid);
				//Worst PF charged iso
				if((abs(phoSCEta->at(ip))<1.442 || abs(phoSCEta->at(ip))>1.556) && phoPFChWorstIso->at(ip)<15 ){
				  pid=5;  hPID->Fill(pid);
				  //ECAL Photon ISO
				  if((abs(phoSCEta->at(ip))<1.442 || abs(phoSCEta->at(ip))>1.556) && phoPFPhoIso->at(ip)<15){
				    pid=6;  hPID->Fill(pid);
				    //HoverE
				    if((abs(phoSCEta->at(ip))<1.442 && phoHoverE->at(ip)<0.08) || (abs(phoSCEta->at(ip))>1.556 && phoHoverE->at(ip) <0.05 )){
				      pid=7;  hPID->Fill(pid);
				      //SigmaIetaIeta
				      if((abs(phoSCEta->at(ip))<1.442 && phoSigmaIEtaIEtaFull5x5->at(ip)<0.012 )||(abs(phoSCEta->at(ip))>1.556 && phoSigmaIEtaIEtaFull5x5->at(ip)<0.027 )){
					pid=8;  hPID->Fill(pid);
					//MVA photon
					if((phoIDMVA->at(ip) > 0.4 && abs(phoSCEta->at(ip))<1.442 ) || (phoIDMVA->at(ip) > 0.3 && abs(phoSCEta->at(ip))>1.556 )){
					  pid=9;  hPID->Fill(pid);
					  passingPhotons++;
					  _pPt=phoEt->at(ip); 
					  hPhotonPt->Fill(_pPt);
					  hDeltaR->Fill(dr1); 
					  hDeltaR->Fill(dr2);
					  hMll->Fill(ll[iz].M());
					  photon[ip].SetPtEtaPhiM(phoEt->at(ip),phoSCEta->at(ip),phoPhi->at(ip),_gMass);
					  llg[izg]= ll[iz]+photon[ip];
					  llgm[izg]=llg[izg].M(); hMllg->Fill(llgm[izg]);
					  std::cout<<entry<< "th entry, "<<" pt("<<ie1<<") = "<<elePt->at(ie1)<<", pt("<<ie2<<") = "<<elePt->at(ie2)<<", "<<iz<<"th Mll="<<ll[iz].M()<<" GeV/c2, phoEt("<<ip<<")="<<phoEt->at(ip)<<" GeV/c, and "<<izg<<"th Zg candidate, mass="<<llgm[izg]<<" GeV/c2"<<std::endl;
					  goodZg = true;
					  _phoEt[ip]=_pPt;
					  _mll[iz]=ll[iz].M();
					  _dr1 = dr1;
					  _dr2 = dr1;
					  _mllg[izg]=llgm[izg];
					  izg++;
					}//MVA
				      }//SigmaIetaIeta
				    }//HoverE
				 
				}//ECAL photon ISO
			      }//Worst Charged ISO
			    }//R9
			  } //DeltaR
			}//end of photon loop
		      hPhotons->Fill(passingPhotons);
		    }
		    iz++;
		  }
	      }//2nd lepton loop
	  }//first lepton loop
	  int mizg = izg; //max nb of zg cands in.
	  if (goodZg)
	    {
	      for(int iizg=0; iizg<mizg;iizg++){
		std::cout<<" izg = "<<iizg<<", Zg candidate mass="<<llgm[iizg]<<" GeV/c2"<<std::endl;
	       	
	       	_zgt->Fill();
	      }
	    }
	}
    }

  //Drawing some histos
  // TCanvas* c1 = new TCanvas ("c1","c1",600,600);
  // c1->Divide(2,2);
  // c1->cd(1);
  // hPhotonPt->Draw();
  // c1->cd(2);
  // hDeltaR->Draw();
  // c1->cd(3);
  // hMll->Draw();
  // c1->cd(4);
  // hMllg->Draw();
  // c1->Update();
  // c1->Draw();
  
  // c1->SaveAs("testOutputhistos.root");
  hPhotonPt->SaveAs("hist_PhotonPt_"+_inputFileName);
  hDeltaR->SaveAs("hist_DeltaR_"+_inputFileName);
  hMll->SaveAs("hist_Mll_"+_inputFileName);
  hMllg->SaveAs("hist_Mllg_"+_inputFileName);

  // c1->Close();
  // TCanvas* c2 = new TCanvas ("c2","c2",600,600);
  // c2->cd();
  // hPhotons->Draw();
  // c2->Update();
  // c2->Draw();
  hPhotons->SaveAs("test_cuts.root");
 
  hPID->Draw();
  hPID->SaveAs("hist_phoID_"+_inputFileName);
  _zgt->GetDirectory()->cd();
  //c2->Close();
  zgf->Write();
}

  //TCut fullcut("eleDr03TkSumPt[]<0.18 && elePt[0]>20 && elePt[]>15");
  // tree->Draw("phoEt>>hPhotonPt",fullcut,"");
  // hPhotonPt->SetFillColor(kYellow);
  // gStyle->SetHistLineColor(1);
  // gStyle->SetHistLineStyle(0);
  // gStyle->SetHistLineWidth(1);

    // if(!simple){
    //   // //first, loop over electrons to apply all ID selections.
    //   //      if(it>0){  j++;} //is it still necessary?
    //   _ePt[i] = elePt->at(i);
    //   if(i==0 && _ePt[0] > 20.0){ //be there only once for the first lepton.
	    
    //     /// MVA Selection
	    
    //     if ((abs(eleEta->at(0))<0.8 && ( eleIDMVATrg->at(0)>0.972153 ))
    // 	  ||((abs(eleEta->at(0))>0.8 && abs(eleEta->at(0))<1.479) && eleIDMVATrg->at(0)>0.922126)
    // 	  ||( (abs(eleEta->at(0))>1.479 && abs(eleEta->at(0))<2.5) && eleIDMVATrg->at(0)>0.610764 ))
	      
  
    //   else if( i>0 && _ePt[i] > 15){ //if it's not the first lepton, and it has pt>15
    //     iip=0; //reset photon counter.
    //     /// MVA Selection                                                                                                                                            
    //     if ((abs(eleEta->at(i))<0.8 && ( eleIDMVATrg->at(i)>0.972153 ))
    // 	  ||((abs(eleEta->at(i))>0.8 && abs(eleEta->at(i))<1.479) && eleIDMVATrg->at(i)>0.922126)
    // 	  ||( (abs(eleEta->at(i))>1.479 && abs(eleEta->at(i))<2.5) && eleIDMVATrg->at(i)>0.610764 ))
	      
    // 	//track-matching
    // 	//if(abs(eleEta->at(i))<1.442 && (  eledEtaAtVtx->at(i)<0.0095 &&  eledPhiAtVtx->at(i)<0.065  ))
	      
    // 	//relative trk iso
    // 	if((abs(eleEta->at(i))<1.442 && eleDr03TkSumPt->at(i)<0.18 )
    // 	   ||(abs(eleEta->at(i))>1.556 && eleDr03TkSumPt->at(i)<0.18 ))
		
    // 	  //ECAL cluster iso
    // 	  if((abs(eleEta->at(i))<1.442 && elePFClusEcalIso->at(i)<0.37 )
    // 	     ||(abs(eleEta->at(i))>1.556 && elePFClusEcalIso->at(i)<0.45 ))
		  
    // 	    //relative HCAL cluster iso                                                                                                                            
    // 	    if((abs(eleEta->at(i))<1.442 && elePFClusHcalIso->at(i)<0.25 )
    // 	       ||(abs(eleEta->at(i))>1.556 && elePFClusHcalIso->at(i)<0.28 ))
		    
    // 	      //H over E                                                                                                                                        
    // 	      if((abs(eleEta->at(i))<1.442 && eleHoverE->at(i)<0.09 )
    // 		 ||(abs(eleEta->at(i))>1.556 && eleHoverE->at(i)<0.09 ))
		      
    // 		//Full 5x5 sigma i eta i eta                                                                                                                           
    // 		if((abs(eleEta->at(i))<1.442 && eleSigmaIEtaIEtaFull5x5->at(i)<0.012 )
    // 		   ||(abs(eleEta->at(i))>1.556 && eleSigmaIEtaIEtaFull5x5->at(i)<0.033 ))
    // 		  ///if all these conditions are passed, build a Z.
    // 		  {
    // 		    for(int j=i; j>=1; j--)
    // 		      {
    // 			//the pt of leptons is at least > 15, in addition must make sure that one of them is > 20 
    // 			//(it can happen if the pair does not contain the zeroth lepton)
    // 			if (_ePt[j-1] > 20 || _ePt[i] > 20){ 
    // 			  _eEta[i] = eleEta->at(i);
    // 			  _ePhi[i] = elePhi->at(i);
				
    // 			  _eEta[j-1] = eleEta->at(j-1);
    // 			  _ePhi[j-1] = elePhi->at(j-1);
				
    // 			  _hStats->Fill(1);
    // 			  lepton2.SetPtEtaPhiM(_ePt[i],_eEta[i],_ePhi[i],_eMass); // latest lepton
    // 			  lepton1.SetPtEtaPhiM(_ePt[j-1],_eEta[j-1],_ePhi[j-1],_eMass); //previous.
    // 			  ll = lepton1 + lepton2; 
    // 			  Mll = (float) ll.M(); 
				
    // 			  if ( Mll > 50 ){ // if the two leptons pass the Z mass cut, look for a photon.
    // 			    {_hStats->Fill(2);
    // 			      //photon loop
    // 			      for(int ip =0;  ip< _pSize ;ip++){
    // 				_gammaEt.clear();//reinit for every new photon?
    // 				gWrite =false; // init for every new llg-pair.
    // 				//if scanning i-th photon, must make sure it passes the pt-eta cut
    // 				// out-of-gap-region check for all successive cuts (may be superfluous)
    // 				if(phoEt->at(ip)>15 && abs(phoSCEta->at(ip))<2.5)
    // 				  { _hStats->Fill(3);
    // 				    //MVA cut
    // 				    if( (phoIDMVA->at(ip) > 0.4 && abs(phoEta->at(ip))<1.442 ) ||
    // 					(phoIDMVA->at(ip) > 0.3 && abs(phoEta->at(ip))>1.556 ) )
    // 				      { _hStats->Fill(4);
    // 					//R9 (cut only on endcaps)
    // 					if((abs(phoEta->at(ip))>1.556 && phoR9->at(ip)>0.85) || abs(phoEta->at(ip))<1.442)
    // 					  { _hStats->Fill(5);
    // 					    //Worst PF charged ISO
    // 					    if(( abs(phoEta->at(ip))<1.442 &&  phoPFChWorstIso->at(ip) < 15 )
    // 					       ||( abs(phoEta->at(ip))>1.556 &&  phoPFChWorstIso->at(ip) < 15  ))
    // 					      { _hStats->Fill(6);
    // 						//PF photon ISO
    // 						if(( abs(phoEta->at(ip))<1.442 &&  phoPFChIso->at(ip) < 15 )
    // 						   ||( abs(phoEta->at(ip))>1.556 &&  phoPFChIso->at(ip) < 15  ))
    // 						  { _hStats->Fill(7);
    // 						    //H over E
    // 						    if((  abs(phoEta->at(ip))<1.442 &&  phoHoverE->at(ip) < 0.08 )
    // 						       ||( abs(phoEta->at(ip))>1.556 &&  phoHoverE->at(ip) < 0.05  ))
    // 						      { _hStats->Fill(8);
    // 							//Full 5x5 sigma i eta i eta
    // 							if((abs(phoEta->at(i))<1.442 && phoSigmaIEtaIEtaFull5x5->at(i)<0.012 )
    // 							   ||(abs(phoEta->at(i))>1.556 && phoSigmaIEtaIEtaFull5x5->at(i)<0.027 ))
    // 							  { _hStats->Fill(9);
    // 							    {
    // 							      _phoEt[ip] = phoEt->at(ip);
    // 							      _phoEta[ip] = phoEta->at(ip);
    // 							      _phoPhi[ip] = phoPhi->at(ip);
    // 							      dR1 = deltaR(_eEta[i],_phoEta[ip],_ePhi[i],_phoPhi[ip]);
    // 							      dR2 = deltaR(_eEta[j-1],_phoEta[ip],_ePhi[j-1],_phoPhi[ip]);
    // 							      //isolated photon!
								    
    // 							      if (dR1 > 0.7 && dR2 > 0.7){
    // 								//this photon is good, return the Z cand mass
    // 								std::cout <<entry<< "th evt | ("<<i<<","<<(j-1)<<"), pt("<<_ePt[i]<<","<<_ePt[j-1]<<"), m_ll = "<<Mll;
								      
    // 								//building the Zgamma candidate	
    // 								photon.SetPtEtaPhiM(_phoEt[ip],_phoEta[ip],_phoPhi[ip],_gMass);	
    // 								llg = ll + photon; 
    // 								float Mllg = (float) llg.M();		
								      
    // 								//if the photon was not already used, 
    // 								//keep it (in a separate collection from phoEt)
								      
    // 								if (_gammaEt.empty()){
    // 								  std::cout <<", photon("<<ip<<") Et= "<< _phoEt[ip];
    // 								  if(_goodPhotonsEt[ip]==_phoEt[ip]){
    // 								    std::cout <<" again, so skipped."<<std::endl; 
    // 								    continue;}
    // 								  else{
    // 								    _hStats->Fill(10);
    // 								    iip = ip; //for further photon double counting test
    // 								    _goodPhotonsEt[iip]=_phoEt[ip];			
    // 								    _gammaEt.push_back(_goodPhotonsEt[iip]);
    // 								    _llm.push_back(Mll); 
    // 								    _deltaR.push_back(dR1);
    // 								    _deltaR.push_back(dR2); // saving the two delta-R.
    // 								    _llgm.push_back(Mllg);
    // 								    // _TREE.outLeaf.gammaEt=_gammaEt;
    // 								    // _TREE.outLeaf.deltaR=_deltaR;
    // 								    // _TREE.outLeaf.eegMass=_llgm;
    // 								    // _TREE.outLeaf.eeMass=_llm;
									  
    // 								    _outputTree->Fill();
									  
    // 								  }
    // 								}
    // 								std::cout<<", llg mass = "<<Mllg<< std::endl;		      
    // 								//this event has passed all cuts, so we can safely say this:
    // 								gWrite =true;
    // 								// add selected gen photons later
    // 								// if 
    // 							      }
    // 							    }
    // 							  }
    // 						      }
    // 						  }
    // 					      }
    // 					  }
    // 				      }
    // 				  }
    // 			      }//end of the photon loop
    // 			    }
    // 			  }//if the Mll < 50, go to next lepton.
    // 			}
    // 		      }// end of the second lepton loop
    // 		  }// end of ith lepton test selection 
    //   }//end of ith lepton kinematic test 
    //   else continue; //if first lepton < 20, and any of the trailing leptons have pt<15, skip the lepton.
    // }
 
