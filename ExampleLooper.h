#pragma once
// Root Packages
#include "TTree.h"
#include "TChain.h"

#include "PhaseIIMuonTriggerNtuple/MuonTriggerNtAna.h"
#include "PhaseIIMuonTriggerNtuple/TriggerTools.h"
#include "PhaseIIMuonTriggerNtuple/Event.h"
#include "PhaseIIMuonTriggerNtuple/RpcHit.h"
#include "RpcRoi/RpcRoIManager.h"
#include "RpcRoi/RpcRoI.h"
#include "GetArea.h"
#include "RpcAnalysis.h"

#include <fstream>

using namespace std;
using namespace TrigNtup;


#include <string>
#include <sstream>
#include <vector>
#include <iterator>


class ExampleLooper : public MuonTriggerNtAna {
public:
  ExampleLooper(TChain *c, TString outputfilename="out.root",TString dataset="zb");
  TFile* f;
  virtual ~ExampleLooper(){};
    
  GetArea obj;
  TrigNtup::TriggerTools trigTools;
    // Begin is called before looping on entries
    virtual void Begin(TTree *tree);

    // Terminate is called after looping is finished
    virtual void Terminate();

    // Main event loop function
    virtual Bool_t  Process(Long64_t entry);

    // Dump cutflow - if derived class uses different cut ordering,
    // override this method
    virtual void dumpEventCounters();

    // debug check
    bool debugEvent();

    // Dump Offline muon info
    void dumpMuon();

    // Dump MDT hit info
    void dumpMdtHit();

    // Dump MDT segment info
    void dumpMdtSegment();


    // std::vector<string> split(const string& s, const string& s1, const bool keep_empty );
    //float GetArea(  TString s,int s1, int s2);
    //checking Data/MC and good lumi section
    bool CheckLSQuality(unsigned int,unsigned int, bool);


    TH2F * hm_etaphi;
    TH2F * hm_rz;
 
    TH1F * hm_lumi;

    TH2F * hm_etaphi_[14];
    TH1F * hitrateVsZ_[14];

    TH2F * hm_summary_L; 
    TH2F * hm_summary_S; 
    TH1F * hm_eta_[14];
    TH1F * hitrateVsLumi_chamber_[14];
    TH1F * hitrateVsLumi_chamber_withoutSA_[14];
    TH1F * hm_time;
    TH1F * hm_adc;
    TH1F *hm_event;
    //    TH1F * hit[14][17];A
    //TH1F * hitrateVsLumi_chamber_eta_[14];
    std::vector<std::vector<TH1F* > > hitrateVsLumi_chamber_eta_;
    double mdt_eta;
    double mdt_phi;
    double rho; 
    float livetime;
    float length;
    float diameter;



    std::map<TString,double> Area_;
    std::map<TString,int* > Bin_;
    int Bin;

    std::map<TString,std::vector<Float_t> >  VarBin_;   
    std::map<TString, std::pair<int,int> >  NTubes;

    std::vector<TString> dtype_;
    std::vector<TString> chamber_;
    std::vector<Int_t> etasta_;
    
    TString dataset_;
    
    ClassDef(ExampleLooper, 1);

protected:
    TChain *m_input_chain;  // input chain being processed
    uint    n_readin;       // Event counters

    
    //    std::string m_data_dir;

};

