#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>

#include "TCanvas.h"
#include "MuTrigNtExampleAna/ExampleLooper.h"
#include "MuTrigNtExampleAna/RpcAnalysis.h"
#include "TFile.h"

using namespace std;
using namespace TrigNtup;

/*--------------------------------------------------------------------------------*/
// ExampleLooper Constructor
/*--------------------------------------------------------------------------------*/
ExampleLooper::ExampleLooper(TChain *c, TString outputfilename, TString dataset) {
    m_input_chain = c;
    n_readin = 0;
    dataset_ = dataset;
    std::cout<< "dataset_____ = " <<dataset_<<std::endl;
    // dataset can be: zb or ex
    f = new TFile(outputfilename,"RECREATE");
}

/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
/*--------------------------------------------------------------------------------*/
void ExampleLooper::Begin(TTree* /*tree*/) {
    MuonTriggerNtAna::Begin(0);
    if (m_dbg) std::cout << "ExampleLooper::Begin\n";


    // for getting the trigger information
    //trigTools.init("root://atldtn1.slac.stanford.edu///atlas/rucio/user.asoffa:user.asoffa.12279820.EXT0._000009.muTrigNt.root"); //ZB
    //trigTools.init("root://atldtn1.slac.stanford.edu///atlas/rucio/user.asoffa:user.asoffa.12235124.EXT0._000397.muTrigNt.root"); // empty
    
    if (dataset_=="ex"){
      trigTools.init("express_stream_ex.root");
    }





    mdt_eta = -99;
    mdt_phi = -99;
    rho = -99;
     
    livetime = 900 * 1e-9; //in seconds 
    length = 100; //cm
    diameter = 3; //cm


    hm_lumi =  new TH1F("hm_lumi" ,"hm_lumi",13000,0.,13000.);    
    hm_time = new TH1F("hm_time","hm_time",1000,-1000,4000);  
    hm_etaphi = new TH2F("hm_etaphi","hm_etaphi",17,-8.5,8.5,9,0.5,8.5);
    hm_adc = new TH1F ("hm_adc","hm_adc",1000, -50,1500);    
    hm_rz = new TH2F("hm_rz","hm_rz", 500,-25,25,200,0,20);
    hm_summary_S = new TH2F ("hm_summary_S","hm_summary_S",17,-8.5,8.5,7,0.5,7.5); 
    hm_summary_L = new TH2F ("hm_summary_L","hm_summary_L",17,-8.5,8.5,7,0.5,7.5); 
    hm_event =new TH1F ("hm_events","hm_events",7,0,7);


    etasta_.clear();
    etasta_.push_back(-8);    etasta_.push_back(-7);    etasta_.push_back(-6);    etasta_.push_back(-5);    etasta_.push_back(-4);    etasta_.push_back(-3);    etasta_.push_back(-2);    etasta_.push_back(-1);    etasta_.push_back(1);    etasta_.push_back(2);    etasta_.push_back(3);    etasta_.push_back(4);    etasta_.push_back(5);    etasta_.push_back(6);    etasta_.push_back(7);    etasta_.push_back(8);

    dtype_.clear();
    dtype_.push_back("BIS");
    dtype_.push_back("BMS");
    dtype_.push_back("BOS");
    dtype_.push_back("EIS");
    dtype_.push_back("EMS");
    dtype_.push_back("EOS");
    dtype_.push_back("EES");
    dtype_.push_back("BIL");
    dtype_.push_back("BML");
    dtype_.push_back("BOL");
    dtype_.push_back("EIL");
    dtype_.push_back("EML");
    dtype_.push_back("EOL");
    dtype_.push_back("EEL");

    VarBin_.clear();
    VarBin_[dtype_[0]]  ={0.,200.,400.,600.,720.};
    VarBin_[dtype_[1]]  ={0.,200.,450.,850.,950.};
    VarBin_[dtype_[2]]  ={0.,400.,800.,1200.};
    VarBin_[dtype_[3]]  ={220.,320.,450.,600.};
    VarBin_[dtype_[4]]  ={180.,280.,500.,1000.};
    VarBin_[dtype_[5]]  ={300.,400.,600.,1200.};
    VarBin_[dtype_[6]]  ={220.,320.,450.,600.};          
    VarBin_[dtype_[7]]  ={0.,200.,400.,600.,720.};
    VarBin_[dtype_[8]]  ={0.,200.,450.,850.,950.};
    VarBin_[dtype_[9]]  ={0.,400.,800.,1200.};
    VarBin_[dtype_[10]]  ={220.,320.,450.,600.};
    VarBin_[dtype_[11]]  ={180.,280.,500.,1000.};
    VarBin_[dtype_[12]]  ={300.,400.,600.,1200.};
    VarBin_[dtype_[13]]  ={220.,320.,450.,600.};          

    std::vector<TH1F*> tmp1(etasta_.size(),0);   
    hitrateVsLumi_chamber_eta_.assign(dtype_.size(),tmp1);
    TString hname;


    

    for (int id=0; id<(int)dtype_.size() ; id++){
      hm_etaphi_[id] = new TH2F("hm_etaphi_"+dtype_[id], "hm_etaphi_"+dtype_[id], 17,-8.5,8.5,8,0.5,8.5);
      hm_eta_[id] = new TH1F("hm_eta_"+dtype_[id], "hm_eta_"+dtype_[id], 17,-8.5,8.5);
      const int nele = (int) VarBin_[dtype_[id]].size();
      Float_t a_arr[nele];
      for (int i=0; i<(int) nele; i++){
	a_arr[i] = VarBin_[dtype_[id]][i] ;}
      hitrateVsZ_[id] = new TH1F("hitrateVsZ_"+dtype_[id],"hitrateVsZ_"+dtype_[id],nele-1,a_arr);
      hitrateVsLumi_chamber_[id] = new TH1F("hitrateVsLumi_chamber_"+dtype_[id],"hitrateVsLumi_chamber_"+dtype_[id],13000,0.,13000);
      hitrateVsLumi_chamber_withoutSA_[id] = new TH1F("hitrateVsLumi_chamber_withoutSA_"+dtype_[id],"hitrateVsLumi_chamber_withoutSA_"+dtype_[id],13000,0.,13000);
      
      for(int iieta =0; iieta < (int)etasta_.size(); iieta++){
	TString etasta_string ;
	etasta_string.Form("%d", etasta_[iieta]);
	hname = "hitrateVsLumi_chamber_"+dtype_[id]+"_eta_"+etasta_string;
	hitrateVsLumi_chamber_eta_.at(id).at(iieta) = new TH1F(hname, hname,13000,0.,13000);
      }
    }
    obj.initialize();
    obj.SetArea();

}

Bool_t ExampleLooper::CheckLSQuality(unsigned int run_, unsigned int lumisection, bool ismc_){

  Bool_t isgood = false;

  // if running over MC sample run over whole sample
  if (ismc_ == true){
    isgood  =true;
  }


  //different runumber with good lumisection
  //------------------
  if (run_ == 311481) { 
    if ((lumisection >= 131 && lumisection <= 209 )
	|| (lumisection >= 211 && lumisection <= 234 )
	|| (lumisection >= 238 && lumisection <= 329 ) 
	|| (lumisection >= 331 && lumisection <= 885 )){
      isgood  =true;
    }
  }
  //--------------------


  //------------------
  if (run_ == 311402) { 
    if ((lumisection >= 88 && lumisection <= 267 )
	|| (lumisection >= 269 && lumisection <= 352 )
	|| (lumisection >= 355 && lumisection <= 359 ) 
	|| (lumisection >= 361 && lumisection <= 405 )
	|| (lumisection >= 407 && lumisection <= 816 ) ){
      isgood  =true;
    }
  }
  //--------------------



  //------------------
  if (run_ == 311287) { 
    if ((lumisection >= 70 && lumisection <= 259 )
	|| (lumisection >= 293 && lumisection <= 839 )
	|| (lumisection >= 841 && lumisection <= 843 ) 
	|| (lumisection >= 847 && lumisection <= 927 ) 
	|| (lumisection >= 929 && lumisection <= 944 )){
      isgood  =true;
    }
  }
  //--------------------
  //--------------------
  if (run_ == 311365) { 
    if ((lumisection >= 153 && lumisection <= 282 )
	|| (lumisection >= 284 && lumisection <= 406 )
	|| (lumisection >= 409 && lumisection <= 698 )){ 
      isgood  =true;
    }
  }



  if(run_ ==  309314 ){
    if ((lumisection >=267 && lumisection <= 288)
	|| (lumisection >= 290 && lumisection <= 291)
	|| (lumisection >= 293 && lumisection <= 311)
	|| (lumisection >= 313 && lumisection <= 318)
	|| (lumisection >= 320 && lumisection <= 353)
	|| (lumisection >= 355 && lumisection <= 356)
	|| (lumisection >= 358 && lumisection <= 371)
	|| (lumisection >= 373 && lumisection <= 388)
	|| (lumisection >= 390 && lumisection <= 391)
	|| (lumisection >= 393 && lumisection <= 395)
	|| (lumisection >= 397 && lumisection <= 399)
	|| (lumisection >= 401 && lumisection <= 404)
	|| (lumisection >= 406 && lumisection <= 411)
	|| (lumisection >= 413 && lumisection <= 433)){
      isgood  =true;
    }
  }



  if(run_ ==  310216){
    if (lumisection >=172  && lumisection <= 525){
      isgood  =true;
    }
  }

  if(run_ == 309375){
    if ((lumisection >=185 && lumisection <=205)
	|| (lumisection >= 207 && lumisection <= 231)
	|| (lumisection >= 233 && lumisection <= 250)
	|| (lumisection >= 252 && lumisection <= 299)
	|| (lumisection >= 301 && lumisection <= 314)
	|| (lumisection >= 406 && lumisection <= 453)){
      isgood  =true;
    }
  }
  if(run_ == 310574 ){
    isgood  =true;}
  //--------------------
  return isgood;
}

/*--------------------------------------------------------------------------------*/
// Main process loop function 
/*--------------------------------------------------------------------------------*/
Bool_t ExampleLooper::Process(Long64_t entry) {
    // Communicate tree entry number to MuonTriggerNtObject
    GetEntry(entry);
    MuonTriggerNtAna::clearObjects();

    

    // Chain entry not the same as tree entry
    //static Long64_t chainEntry = -1;
    m_chainEntry++;
    if (m_dbg)std::cout << typeid(nt).name()<<std::endl; 
    if (m_dbg || m_chainEntry % 1000 == 0) {
      std::cout << "\n**** Processing entry " << std::setw(6) << m_chainEntry
		<< " run " << std::setw(6)   << nt.evt()->run
		<< " event " << std::setw(7) << nt.evt()->eventNumber 
		<< " av lum "<< std::setw(7) << nt.evt()->lbAvgLumi
		<< " lb "<< std::setw(7) << nt.evt()->lb
		//<< "bcid " <<std::setw(7) << nt.evt()->lbBcidLumi
		<< "bunches "<<std::setw(7)  << nt.evt()->bunches
	        //<< "pile up in " <<std::setw(7)  << nt.evt()->avgMu
		//<< "pile up out" <<std::setw(7)  << nt.evt()->actMu
		<< " ****\n";
    }
    

    bool isMC_ = nt.evt()->isMC;
    unsigned int thisLS = nt.evt()->lb ;
    bool isGoodLS = CheckLSQuality(nt.evt()->run, thisLS,isMC_);
    
     
    //bool pass_HLT_j40_L1ZB = trigTools.passTrigger(nt.evt()->trigBits, "HLT_j40_L1ZB");
    //bool pass_HLT_noalg_zb_L1ZB = trigTools.passTrigger(nt.evt()->trigBits, "HLT_noalg_zb_L1ZB");
    bool pass_HLT_noalg_L1RD0_EMPTY ="false";
    if (dataset_=="ex"){pass_HLT_noalg_L1RD0_EMPTY =  trigTools.passTrigger(nt.evt()->trigBits, "HLT_noalg_L1RD0_EMPTY");}
    if (isGoodLS){
         
      ++n_readin;
      hm_lumi->Fill(nt.evt()->lbAvgLumi);


      //std::cout << " yes pass the triggers" << std::endl;
      if (dataset_=="ex"){
	if(pass_HLT_noalg_L1RD0_EMPTY){
	  dumpMdtHit();
	}
      }
      if(dataset_=="zb"){
	//if((pass_HLT_j40_L1ZB ) || (pass_HLT_noalg_zb_L1ZB)){ //no trigger on ZB
	dumpMdtHit();
      }
    }
    return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called
/*--------------------------------------------------------------------------------*/
void ExampleLooper::Terminate()
{
  MuonTriggerNtAna::Terminate();
  if(m_dbg) std::cout << "ExampleLooper::Terminate\n";
  hm_event->SetBinContent(1,n_readin);
  


    for (int id=0; id<(int)dtype_.size() ; id++){
      double denominator =  livetime*0.01; //0.01 for cm
      
      hm_etaphi_[id]->Scale(1./denominator);
      hm_eta_[id]->Scale(1./denominator);
      hitrateVsZ_[id]->Scale(1./denominator);      
      hitrateVsLumi_chamber_[id]->Scale(1./denominator);
      hitrateVsLumi_chamber_withoutSA_[id]->Scale(1./livetime);
      //hitrateVsLumi_chamber_[id]->Scale(1./livetime);



      for(int nb = 1; nb <= hm_eta_[id]->GetNbinsX(); nb++){
        int xlab = nb-9;
        int ylab = id+1;
	float rate = hm_eta_[id]->GetBinContent(nb);       
        if(id <= 6) { 
	  hm_summary_S->Fill(xlab,ylab,rate);
	}
        if(id > 6) { 
          ylab = id-6  ;
	  hm_summary_L->Fill(xlab,ylab,rate);
	}
      }

      for(int iieta =0; iieta< (int)etasta_.size(); iieta++){          
	hitrateVsLumi_chamber_eta_[id][iieta]->Scale(1./denominator);
      }
    }
    dumpEventCounters();
    //TFile *f = new TFile("out_mdt_adc.root","RECREATE");


    f->cd();



    for (int id=0; id<(int)dtype_.size() ; id++){
      hm_etaphi_[id]->Write();
      hm_eta_[id]->Write();
      hitrateVsZ_[id]->Write();
      hitrateVsLumi_chamber_[id]->Write();
      hitrateVsLumi_chamber_withoutSA_[id]->Write();
      for(int iieta =0; iieta<(int)etasta_.size(); iieta++){
       	hitrateVsLumi_chamber_eta_[id][iieta]->Write(); 
      }
    }

    hm_summary_S->Write();
    hm_summary_L->Write();
    hm_etaphi->Write();
    hm_adc->Write();
    hm_time->Write();
    hm_rz->Write();
    hm_lumi->Write();
    hm_event->Write();
    f->Close();
}

/*--------------------------------------------------------------------------------*/
// Event counters
/*--------------------------------------------------------------------------------*/
void ExampleLooper::dumpEventCounters()
{
    std::cout << "\n"
              << "ExampleLooper event counters\n"
              << "read in : " << n_readin << "\n";
}

/*--------------------------------------------------------------------------------*/
// Debug event
/*--------------------------------------------------------------------------------*/
bool ExampleLooper::debugEvent()
{
    return false;
}

/*--------------------------------------------------------------------------------*/
// Dump Offline muon
/*--------------------------------------------------------------------------------*/
void ExampleLooper::dumpMuon()
{
    for(uint iMu = 0; iMu < nt.muo()->size(); ++iMu) {
       TrigNtup::Muon* mu = &nt.muo()->at(iMu); 
       mu->print();
    }
             
}

/*--------------------------------------------------------------------------------*/
// Dump MDT hits info
/*--------------------------------------------------------------------------------*/
void ExampleLooper::dumpMdtHit()
{
    TVector3 v1; 
    for (uint iMdtHit = 0; iMdtHit < nt.mdtHit()->size(); ++iMdtHit) {
      MdtHit* mdtHit = &nt.mdtHit()->at(iMdtHit);

      if (mdtHit->adc < 80 ) continue;
      float zrange = mdtHit->globalPos.Z()/10.;
      rho = sqrt((mdtHit->globalPos.X()*mdtHit->globalPos.X())+(mdtHit->globalPos.Y()*mdtHit->globalPos.Y()));

      if (m_dbg)std::cout<< "no of MDT hit: "<<nt.mdtHit()->size()  <<std::endl;
      if (m_dbg) std::cout << "z :" <<mdtHit->globalPos.Z() << "rho: "  << rho <<std::endl;
      if (m_dbg) std::cout <<" Eta: "<<mdt_eta << "Phi:"<<mdt_phi<< "  On segment  :" <<mdtHit->isOnSegment  <<std::endl;

      v1.SetXYZ(mdtHit->globalPos.X(),mdtHit->globalPos.Y(),mdtHit->globalPos.Z());
      mdt_eta = v1.Eta();
      mdt_phi = v1.Phi();
     
      hm_time->Fill(mdtHit->driftTime);
      hm_etaphi->Fill(mdtHit->stationEta(),mdtHit->stationPhi());
      hm_adc->Fill(mdtHit->adc);
      hm_rz->Fill(mdtHit->globalPos.Z()/1000.,rho/1000.);
      
      TString stationname = mdtHit->stationNameString().substr(0,3);
      if (m_dbg) cout <<" stationame" <<stationname <<std::endl;


      float weight_area = 1./(obj.getarea(stationname,mdtHit->stationEta(),mdtHit->stationPhi()));
      if (m_dbg) cout << "stationname :"<<stationname<<" ieta : " << mdtHit->stationEta() << " iphi: " << mdtHit->stationPhi() << "   weight :" << weight_area << std::endl;
  	
      for (int id=0; id<(int)dtype_.size() ; id++){
	if(mdtHit->stationNameString().substr(0,3)==dtype_[id]){
	  hm_etaphi_[id]->Fill(mdtHit->stationEta(),mdtHit->stationPhi(),weight_area);
          hm_eta_[id]->Fill(mdtHit->stationEta(),weight_area);
	  hitrateVsLumi_chamber_[id]->Fill(nt.evt()->lbAvgLumi,weight_area);
	  hitrateVsLumi_chamber_withoutSA_[id]->Fill(nt.evt()->lbAvgLumi);
          if(mdtHit->stationNameString().substr(0,1)=="B"){
	    hitrateVsZ_[id]->Fill(fabs(zrange),weight_area);}
          if(mdtHit->stationNameString().substr(0,1)=="E"){
	    hitrateVsZ_[id]->Fill(rho/10.,weight_area);}
	  for(int iieta =0; iieta< (int)etasta_.size(); iieta++){          
	    if(mdtHit->stationEta() == (etasta_[iieta])){        
	      hitrateVsLumi_chamber_eta_[id][iieta]->Fill(nt.evt()->lbAvgLumi,weight_area);
	    }
	  }
	}
      }
    }
}

/*--------------------------------------------------------------------------------*/
// Dump MDT Segment info
/*--------------------------------------------------------------------------------*/
void ExampleLooper::dumpMdtSegment()
{
    for (uint iMdtSeg = 0; iMdtSeg < nt.mdtSegment()->size(); ++iMdtSeg) {
        Segment* mdtSeg = &nt.mdtSegment()->at(iMdtSeg);
        if(mdtSeg->nMu()==0) continue; //Not attached to offline muon
        if (m_dbg) std::cout << "Segment " << iMdtSeg << " nMu " << mdtSeg->nMu()
                  << " from muon " << mdtSeg->muIdx[0]
                  << " isMDT " << mdtSeg->isMDT 
                  << " isCSC " << mdtSeg->isCSC
                  << " nMdtHit " << mdtSeg->nHits() << std::endl;
        std::cout << " Hits Idx "; 
        for (uint iHit = 0; iHit < mdtSeg->hitsIdx.size(); ++iHit) {
            std::cout <<  mdtSeg->hitsIdx[iHit] << " ";
        }
        std::cout << std::endl;
    }
    
}








           
