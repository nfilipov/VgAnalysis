#include "VgammaSkim.h"
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>

VgammaSkim::VgammaSkim(TString inputFileName, TString outDir, TString nameDir, TString nameTree, bool isMC, bool basic)
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

  TString skimPartOfName="_EleElePhoSkim"; //expand this!

  _skimmedFileName=outDir+skimPartOfName+".root";
  _fileOut = new TFile(_skimmedFileName,"recreate");
  // _fileOut->mkdir(_nameDir); // removed this to lose the useless folder
  // _fileOut->cd(_nameDir);   // so that trees can be directly accessed.
  _outputTree = new TTree(_nameTree,_nameTree);
  _TREE.InitOutput(_outputTree, basic);      //method of TInputOutputTree
  _hskim= new TH1F("hskim","hskim",2,0,2);    
  if (isMC) { std::cout<< " this is simulation " <<std::endl;}
  
}

VgammaSkim::~VgammaSkim()
{

}


void VgammaSkim::LoopOverInputTree()
{

  _hGpt= new TH1F("hGpt","hGpt",2,0,2); 

  TFile f(_inputFileName,"READ");
  std::cout<<"processing "<<_inputFileName<<std::endl;
  f.cd(_nameDir);
  TTree* tree =(TTree*)gDirectory->Get(_nameTree);

  _TREE.InitInput(tree);

  Long64_t nentries = _TREE.fChain->GetEntries();
  std::cout<<"nentries "<<nentries<<std::endl;
  // if (_isDebugMode)
  nentries=50000;

  for (Long64_t entry=0; entry<nentries; entry++) {
    if (entry < 0) break;
    gWrite=false;
    //   if ((entry%1000000)==0) std::cout<<"entry="<<entry<<std::endl;
    _TREE.GetEntry(entry); 

    //basic pre-selection
    if (_TREE.treeLeaf.nPho<1) continue;
    if (_TREE.treeLeaf.phoEt->at(0)<15 ) continue;
    int _eSize = (int) _TREE.treeLeaf.elePt->size();
    int _pSize = (int) _TREE.treeLeaf.phoEt->size();
    if(_eSize < 2) continue ; // di-leptons' loop, skip single leptons

    //first, loop over electrons to apply all ID selections.

    for (int i=0; i< _eSize ;i++){
      gWrite =false;//init at every new lepton
      //      if(it>0){  j++;} //is it still necessary?
      _ePt[i] = _TREE.treeLeaf.elePt->at(i);
      if(i==0 && _ePt[0] > 20.0){ //be there only once for the first lepton.

	/// MVA Selection

	if ((abs(_TREE.treeLeaf.eleEta->at(0))<0.8 && ( _TREE.treeLeaf.eleIDMVATrg->at(0)>0.972153 ))
	    ||((abs(_TREE.treeLeaf.eleEta->at(0))>0.8 && abs(_TREE.treeLeaf.eleEta->at(0))<1.479) && _TREE.treeLeaf.eleIDMVATrg->at(0)>0.922126)
	    ||( (abs(_TREE.treeLeaf.eleEta->at(0))>1.479 && abs(_TREE.treeLeaf.eleEta->at(0))<2.5) && _TREE.treeLeaf.eleIDMVATrg->at(0)>0.610764 ))
	  
	  //track-matching
	  //	  if((abs(_TREE.treeLeaf.eleEta->at(0)<1.442) && (  _TREE.treeLeaf.eleGSFEta->at(0) )< )

	  //relative trk iso
	  //relative ECAL cluster iso
	  if((abs(_TREE.treeLeaf.eleEta->at(0))<1.442 && _TREE.treeLeaf.elePFClusEcalIso->at(0)<0.37 )
	      ||(abs(_TREE.treeLeaf.eleEta->at(0))>1.556 && _TREE.treeLeaf.elePFClusEcalIso->at(0)<0.45 ))
	    //relative HCAL cluster iso
	      if((abs(_TREE.treeLeaf.eleEta->at(0))<1.442 && _TREE.treeLeaf.elePFClusHcalIso->at(0)<0.25 )
		      ||(abs(_TREE.treeLeaf.eleEta->at(0))>1.556 && _TREE.treeLeaf.elePFClusHcalIso->at(0)<0.28 ))
		    //H over E
		      if((abs(_TREE.treeLeaf.eleEta->at(0))<1.442 && _TREE.treeLeaf.eleHoverE->at(0)<0.09 )
			||(abs(_TREE.treeLeaf.eleEta->at(0))>1.556 && _TREE.treeLeaf.eleHoverE->at(0)<0.09 ))
		      //Sigma-Ieta-Ieta
			if((abs(_TREE.treeLeaf.eleEta->at(0))<1.442 && _TREE.treeLeaf.eleSigmaIEtaIEtaFull5x5->at(0)<0.012 )
			  ||(abs(_TREE.treeLeaf.eleEta->at(0))>1.556 && _TREE.treeLeaf.eleSigmaIEtaIEtaFull5x5->at(0)<0.033 ))
			{
			  _eEta[0] = _TREE.treeLeaf.eleEta->at(i);
			  _ePhi[0] = _TREE.treeLeaf.elePhi->at(i);
			  lepton1.SetPtEtaPhiM(_ePt[0],_eEta[0],_ePhi[0],_eMass); //
			}
	continue; //go to next iteration
      }
      else if( i>0 && _ePt[i] > 15){ //if it's not the first lepton, and it has pt>15

	/// MVA Selection                                                                                                                                            
        if ((abs(_TREE.treeLeaf.eleEta->at(i))<0.8 && ( _TREE.treeLeaf.eleIDMVATrg->at(i)>0.972153 ))
            ||((abs(_TREE.treeLeaf.eleEta->at(i))>0.8 && abs(_TREE.treeLeaf.eleEta->at(i))<1.479) && _TREE.treeLeaf.eleIDMVATrg->at(i)>0.922126)
            ||( (abs(_TREE.treeLeaf.eleEta->at(i))>1.479 && abs(_TREE.treeLeaf.eleEta->at(i))<2.5) && _TREE.treeLeaf.eleIDMVATrg->at(i)>0.610764 ))
	  
	  //ECAL cluster iso
	  if((abs(_TREE.treeLeaf.eleEta->at(i))<1.442 && _TREE.treeLeaf.elePFClusEcalIso->at(i)<0.37 )
	     ||(abs(_TREE.treeLeaf.eleEta->at(i))>1.556 && _TREE.treeLeaf.elePFClusEcalIso->at(i)<0.45 ))
	  
	    //relative HCAL cluster iso                                                                                                                            
	    if((abs(_TREE.treeLeaf.eleEta->at(i))<1.442 && _TREE.treeLeaf.elePFClusHcalIso->at(i)<0.25 )
	       ||(abs(_TREE.treeLeaf.eleEta->at(i))>1.556 && _TREE.treeLeaf.elePFClusHcalIso->at(i)<0.28 ))
	     
	      //H over E                                                                                                                                        
	      if((abs(_TREE.treeLeaf.eleEta->at(i))<1.442 && _TREE.treeLeaf.eleHoverE->at(i)<0.09 )
		 ||(abs(_TREE.treeLeaf.eleEta->at(i))>1.556 && _TREE.treeLeaf.eleHoverE->at(i)<0.09 ))
		
		//Full 5x5 sigma i eta i eta                                                                                                                           
		if((abs(_TREE.treeLeaf.eleEta->at(i))<1.442 && _TREE.treeLeaf.eleSigmaIEtaIEtaFull5x5->at(i)<0.012 )
		   ||(abs(_TREE.treeLeaf.eleEta->at(i))>1.556 && _TREE.treeLeaf.eleSigmaIEtaIEtaFull5x5->at(i)<0.033 ))
		  ///if all these conditions are passed, build a Z.
		  {
		    for(int j=i; j>=1; j--)
		      {
			_eEta[i] = _TREE.treeLeaf.eleEta->at(i);
			_ePhi[i] = _TREE.treeLeaf.elePhi->at(i);
			lepton2.SetPtEtaPhiM(_ePt[i],_eEta[i],_ePhi[i],_eMass); // latest lepton
			_eEta[j-1] = _TREE.treeLeaf.eleEta->at(j-1);
			_ePhi[j-1] = _TREE.treeLeaf.elePhi->at(j-1);
			lepton1.SetPtEtaPhiM(_ePt[j-1],_eEta[j-1],_ePhi[j-1],_eMass); //previous.
			ll = lepton1 + lepton2; 
			float Mll = (float) ll.M(); 
			if ( Mll > 50 ){ // if the two leptons pass the Z mass cut, fill the Z branch.
			  _llm.push_back(Mll);  //and start the photon loop!
			  

			  //photon loop
			  for(int ip =0;  ip< _pSize ;ip++){
			    gWrite =false; // init for every new llg-pair.
			    //if scanning i-th photon, must make sure it passes the pt-eta cut
			    if(_TREE.treeLeaf.phoEt->at(ip)>15 && abs(_TREE.treeLeaf.phoSCEta->at(ip))<2.5)
			      
			      //MVA cut
			      if( (_TREE.treeLeaf.phoIDMVA->at(ip) > 0.4 && abs(_TREE.treeLeaf.phoEta->at(ip))<1.442 ) ||
				  (_TREE.treeLeaf.phoIDMVA->at(ip) > 0.3 && abs(_TREE.treeLeaf.phoEta->at(ip))>1.556 ) )
				
				//R9 (cut only on endcaps)
				if((abs(_TREE.treeLeaf.phoEta->at(ip))>1.556 && _TREE.treeLeaf.phoR9->at(ip)>0.85) || abs(_TREE.treeLeaf.phoEta->at(ip))<1.442)

				  //Worst PF charged ISO
				  if(( abs(_TREE.treeLeaf.phoEta->at(ip))<1.442 &&  _TREE.treeLeaf.phoPFChWorstIso->at(ip) < 15 )
				    ||( abs(_TREE.treeLeaf.phoEta->at(ip))>1.556 &&  _TREE.treeLeaf.phoPFChWorstIso->at(ip) < 15  ))
				    
				    //PF photon ISO
					if(( abs(_TREE.treeLeaf.phoEta->at(ip))<1.442 &&  _TREE.treeLeaf.phoPFChIso->at(ip) < 15 )
					   ||( abs(_TREE.treeLeaf.phoEta->at(ip))>1.556 &&  _TREE.treeLeaf.phoPFChIso->at(ip) < 15  ))

					  //H over E
					      if((  abs(_TREE.treeLeaf.phoEta->at(ip))<1.442 &&  _TREE.treeLeaf.phoHoverE->at(ip) < 0.08 )
						||( abs(_TREE.treeLeaf.phoEta->at(ip))>1.556 &&  _TREE.treeLeaf.phoHoverE->at(ip) < 0.05  ))
						
					  //Full 5x5 sigma i eta i eta
						if((abs(_TREE.treeLeaf.phoEta->at(i))<1.442 && _TREE.treeLeaf.phoSigmaIEtaIEtaFull5x5->at(i)<0.012 )
						   ||(abs(_TREE.treeLeaf.phoEta->at(i))>1.556 && _TREE.treeLeaf.phoSigmaIEtaIEtaFull5x5->at(i)<0.027 ))
					  
						  {
						    _phoEt[ip] = _TREE.treeLeaf.phoEt->at(ip);
						    _phoEta[ip] = _TREE.treeLeaf.phoEta->at(ip);
						    _phoPhi[ip] = _TREE.treeLeaf.phoPhi->at(ip);
						    dR1 = deltaR(_eEta[i],_phoEta[ip],_ePhi[i],_phoPhi[ip]);
						    dR2 = deltaR(_eEta[j-1],_phoEta[ip],_ePhi[j-1],_phoPhi[ip]);
						    //isolated photon!
						   
						    if (dR1 > 0.7 && dR2 > 0.7){
						      _deltaR.push_back(dR1);
						      _deltaR.push_back(dR2); // saving the two delta-R.
						      photon.SetPtEtaPhiM(_phoEt[ip],_phoEta[ip],_phoPhi[ip],_gMass);
						      //keeping the good photon in a separate collection
						      _gammaEt.push_back(_phoEt[ip]);
						      //building the Zgamma candidate
						      llg = ll + photon; 
						      float Mllg = (float) llg.M();
						      _llgm.push_back(Mllg);
						      std::cout <<entry<< "th event -- ll mass = "<<Mll << ", photonEt= "<< _phoEt[ip]<<", llg mass = "<<Mllg<< std::endl;		       	
						      //this event has passed all cuts, so we can safely say this:
						      gWrite =true;
						      // add selected gen photons later
						      _TREE.outLeaf.gammaEt=_gammaEt;
						      _TREE.outLeaf.deltaR=_deltaR;
						      _TREE.outLeaf.eegMass=_llgm;
						    }
						  }
			  }//end of the photon loop
			  _TREE.outLeaf.eeMass=_llm;
			}//if the Mll < 50, go to next lepton.
		      }// end of the second lepton loop
		  }// end of ith lepton test selection 
      }//end of ith lepton kinematic test 
      else continue; //if first lepton < 20, and any of the trailing leptons have pt<15, skip the lepton.
    }//end of lepton-photon loop!


    if(gWrite){ //writing part!
      _outputTree->Fill(); // test
      // add selected gen photons later
      _hskim->Fill(1);
    }
  }//end of entry loop
  

  _fileOut->cd();
  // _fileOut->cd(_nameDir);
  _outputTree->Write(_nameTree,TObject::kOverwrite);
  
  //close output files
  _TREE.fChain = 0;
  std::cout<<"file "<<_fileOut->GetName()<<std::endl<<" closed..."<<std::endl;
  f.Close();
  
}


float VgammaSkim::deltaR(float eta0, float eta1, float phi0, float phi1){
  return sqrt((eta0-eta1)*(eta0-eta1)+(phi0-phi1)*(phi0-phi1));
}
