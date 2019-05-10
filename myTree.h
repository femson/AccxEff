//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 12 12:29:44 2019 by ROOT version 6.06/00
// from TTree myTree/My TTree of dimuons
// found on file: /data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v1/HiForestAOD_ext_merged.root
//////////////////////////////////////////////////////////

#ifndef myTree_h
#define myTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <math.h>
#include <TPaveText.h>
#include <TUnfold.h>
#include <TLorentzVector.h>
#include <vector>
#include <TRandom.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>
#include <fstream>
#include <TLegend.h>
#include <TClonesArray.h>

using namespace std;

class myTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   bool isPbPb=false;
   bool isPr;
   bool isAcc;

//additional variables
   Float_t jpsi_m;
   Float_t jpsi_pt;
   Float_t jpsi_eta;
   Float_t jpsi_phi;
   Float_t jpsi_rap;
   Float_t weight;
   Float_t tnp_weight;
   int triggerIndex_PP = 0; //pp = 0; PbPb = 12;

   // Declaration of leaf types
   UInt_t          eventNb;
   UInt_t          runNb;
   UInt_t          LS;
   Float_t         zVtx;
   Float_t         nPV;
   Int_t           nTrig;
   Int_t           trigPrescale[22];   //[nTrig]
   ULong64_t       HLTriggers;
   Int_t           Reco_QQ_size;
   Int_t           Reco_QQ_type[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_sign[15];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_4mom;
   Int_t           Reco_QQ_mupl_idx[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_idx[15];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_trig[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_isCowboy[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_cosAlpha[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau3D[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr3D[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_cosAlpha3D[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_whichGen[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_VtxProb[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_dca[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_MassErr[15];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_vtx;
   Int_t           Reco_QQ_Ntrk[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxy_muonlessVtx[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxy_muonlessVtx[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxyErr_muonlessVtx[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxyErr_muonlessVtx[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dz_muonlessVtx[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dz_muonlessVtx[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dzErr_muonlessVtx[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dzErr_muonlessVtx[15];   //[Reco_QQ_size]
   Int_t           Reco_mu_size;
   Int_t           Reco_mu_type[6];   //[Reco_mu_size]
   Int_t           Reco_mu_whichGen[6];   //[Reco_mu_size]
   Int_t           Reco_mu_SelectionType[6];   //[Reco_mu_size]
   Int_t           Reco_mu_charge[6];   //[Reco_mu_size]
   TClonesArray    *Reco_mu_4mom;
   ULong64_t       Reco_mu_trig[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_highPurity[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_TrkMuArb[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_TMOneStaTight[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixValHits[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nMuValHits[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkHits[6];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_inner[6];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_global[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixWMea[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkWMea[6];   //[Reco_mu_size]
   Int_t           Reco_mu_StationsMatched[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dxy[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dxyErr[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dz[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dzErr[6];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_inner[6];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_global[6];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_inner[6];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_global[6];   //[Reco_mu_size]
   Int_t           Gen_QQ_size;
   Int_t           Gen_QQ_type[3];   //[Gen_QQ_size]
   TClonesArray    *Gen_QQ_4mom;
   Int_t           Gen_QQ_momId[3];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau[3];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau3D[3];   //[Gen_QQ_size]
   Int_t           Gen_QQ_mupl_idx[3];   //[Gen_QQ_size]
   Int_t           Gen_QQ_mumi_idx[3];   //[Gen_QQ_size]
   Int_t           Gen_QQ_whichRec[3];   //[Gen_QQ_size]
   Int_t           Gen_mu_size;
   Int_t           Gen_mu_type[9];   //[Gen_mu_size]
   Int_t           Gen_mu_charge[9];   //[Gen_mu_size]
   TClonesArray    *Gen_mu_4mom;
   Int_t           Gen_mu_whichRec[9];   //[Gen_mu_size]
   /*
   Int_t           pprimaryVertexFilter; // for PbPb
   Int_t           pPAprimaryVertexFilter; // for pp
   Int_t           pBeamScrapingFilter;
   Int_t           phfCoincFilter2Th4;
   */

   Int_t           Onia2MuMuPAT;
   Int_t           ana_step;
   Int_t           pHBHENoiseFilterResultProducer;
   Int_t           HBHENoiseFilterResult;
   Int_t           HBHENoiseFilterResultRun1;
   Int_t           HBHENoiseFilterResultRun2Loose;
   Int_t           HBHENoiseFilterResultRun2Tight;
   Int_t           HBHEIsoNoiseFilterResult;
   Int_t           pPAprimaryVertexFilter;
   Int_t           pBeamScrapingFilter;
   Int_t           pVertexFilterCutG;
   Int_t           pVertexFilterCutGloose;
   Int_t           pVertexFilterCutGtight;
   Int_t           pVertexFilterCutGplus;
   Int_t           pVertexFilterCutE;
   Int_t           pVertexFilterCutEandG;

   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_LS;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_nTrig;   //!
   TBranch        *b_trigPrescale;   //!
   TBranch        *b_HLTriggers;   //!
   TBranch        *b_Reco_QQ_size;   //!
   TBranch        *b_Reco_QQ_type;   //!
   TBranch        *b_Reco_QQ_sign;   //!
   TBranch        *b_Reco_QQ_4mom;   //!
   TBranch        *b_Reco_QQ_mupl_idx;   //!
   TBranch        *b_Reco_QQ_mumi_idx;   //!
   TBranch        *b_Reco_QQ_trig;   //!
   TBranch        *b_Reco_QQ_isCowboy;   //!
   TBranch        *b_Reco_QQ_ctau;   //!
   TBranch        *b_Reco_QQ_ctauErr;   //!
   TBranch        *b_Reco_QQ_cosAlpha;   //!
   TBranch        *b_Reco_QQ_ctau3D;   //!
   TBranch        *b_Reco_QQ_ctauErr3D;   //!
   TBranch        *b_Reco_QQ_cosAlpha3D;   //!
   TBranch        *b_Reco_QQ_whichGen;   //!
   TBranch        *b_Reco_QQ_VtxProb;   //!
   TBranch        *b_Reco_QQ_dca;   //!
   TBranch        *b_Reco_QQ_MassErr;   //!
   TBranch        *b_Reco_QQ_vtx;   //!
   TBranch        *b_Reco_QQ_Ntrk;   //!
   TBranch        *b_Reco_QQ_mupl_dxy_muonlessVtx;   //!
   TBranch        *b_Reco_QQ_mumi_dxy_muonlessVtx;   //!
   TBranch        *b_Reco_QQ_mupl_dxyErr_muonlessVtx;   //!
   TBranch        *b_Reco_QQ_mumi_dxyErr_muonlessVtx;   //!
   TBranch        *b_Reco_QQ_mupl_dz_muonlessVtx;   //!
   TBranch        *b_Reco_QQ_mumi_dz_muonlessVtx;   //!
   TBranch        *b_Reco_QQ_mupl_dzErr_muonlessVtx;   //!
   TBranch        *b_Reco_QQ_mumi_dzErr_muonlessVtx;   //!
   TBranch        *b_Reco_mu_size;   //!
   TBranch        *b_Reco_mu_type;   //!
   TBranch        *b_Reco_mu_whichGen;   //!
   TBranch        *b_Reco_mu_SelectionType;   //!
   TBranch        *b_Reco_mu_charge;   //!
   TBranch        *b_Reco_mu_4mom;   //!
   TBranch        *b_Reco_mu_trig;   //!
   TBranch        *b_Reco_mu_highPurity;   //!
   TBranch        *b_Reco_mu_TrkMuArb;   //!
   TBranch        *b_Reco_mu_TMOneStaTight;   //!
   TBranch        *b_Reco_mu_nPixValHits;   //!
   TBranch        *b_Reco_mu_nMuValHits;   //!
   TBranch        *b_Reco_mu_nTrkHits;   //!
   TBranch        *b_Reco_mu_normChi2_inner;   //!
   TBranch        *b_Reco_mu_normChi2_global;   //!
   TBranch        *b_Reco_mu_nPixWMea;   //!
   TBranch        *b_Reco_mu_nTrkWMea;   //!
   TBranch        *b_Reco_mu_StationsMatched;   //!
   TBranch        *b_Reco_mu_dxy;   //!
   TBranch        *b_Reco_mu_dxyErr;   //!
   TBranch        *b_Reco_mu_dz;   //!
   TBranch        *b_Reco_mu_dzErr;   //!
   TBranch        *b_Reco_mu_pt_inner;   //!
   TBranch        *b_Reco_mu_pt_global;   //!
   TBranch        *b_Reco_mu_ptErr_inner;   //!
   TBranch        *b_Reco_mu_ptErr_global;   //!
   TBranch        *b_Gen_QQ_size;   //!
   TBranch        *b_Gen_QQ_type;   //!
   TBranch        *b_Gen_QQ_4mom;   //!
   TBranch        *b_Gen_QQ_momId;   //!
   TBranch        *b_Gen_QQ_ctau;   //!
   TBranch        *b_Gen_QQ_ctau3D;   //!
   TBranch        *b_Gen_QQ_mupl_idx;   //!
   TBranch        *b_Gen_QQ_mumi_idx;   //!
   TBranch        *b_Gen_QQ_whichRec;   //!
   TBranch        *b_Gen_mu_size;   //!
   TBranch        *b_Gen_mu_type;   //!
   TBranch        *b_Gen_mu_charge;   //!
   TBranch        *b_Gen_mu_4mom;   //!
   TBranch        *b_Gen_mu_whichRec;   //!

   TBranch        *b_Onia2MuMuPAT;   //!
   TBranch        *b_ana_step;   //!
   TBranch        *b_pHBHENoiseFilterResultProducer;   //!
   TBranch        *b_HBHENoiseFilterResult;   //!
   TBranch        *b_HBHENoiseFilterResultRun1;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Loose;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Tight;   //!
   TBranch        *b_HBHEIsoNoiseFilterResult;   //!
   TBranch        *b_pPAprimaryVertexFilter;   //!
   TBranch        *b_pBeamScrapingFilter;   //!
   TBranch        *b_pVertexFilterCutG;   //!
   TBranch        *b_pVertexFilterCutGloose;   //!
   TBranch        *b_pVertexFilterCutGtight;   //!
   TBranch        *b_pVertexFilterCutGplus;   //!
   TBranch        *b_pVertexFilterCutE;   //!
   TBranch        *b_pVertexFilterCutEandG;   //!


   myTree(TTree *tree=0);
   virtual ~myTree();

   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Acc();
   virtual void     Eff();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Bool_t   isGlobalMuonInAccept2019(TLorentzVector* Muon);
   virtual Bool_t   areMuonsInAcceptance2019(Int_t iRecoQQ);
   virtual Bool_t   areGenMuonsInAcceptance2019(Int_t iGenQQ);
   virtual Bool_t   passQualityCuts2019 (Int_t iRecoQQ);
   virtual Bool_t   isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit);
};

#endif

#ifdef myTree_cxx
myTree::myTree(TTree *tree) : fChain(0) 
{
  TTree* skimTree = NULL;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v1/HiForestAOD_ext_merged.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v1/HiForestAOD_ext_merged.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v1/HiForestAOD_ext_merged.root:/hionia");
      dir->GetObject("myTree",tree);
      TDirectory * dir_skim = (TDirectory*)f->Get("/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v1/HiForestAOD_ext_merged.root:/skimanalysis");
      dir_skim->GetObject("HltTree",skimTree);
   }
   tree->AddFriend(skimTree);
   Init(tree);
}

myTree::~myTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Reco_QQ_4mom = 0;
   Reco_QQ_vtx = 0;
   Reco_mu_4mom = 0;
   Gen_QQ_4mom = 0;
   Gen_mu_4mom = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("LS", &LS, &b_LS);
   fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
   fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
   fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
   fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
   fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
   fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
   fChain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
   fChain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
   fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
   fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
   fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
   fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
   fChain->SetBranchAddress("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha, &b_Reco_QQ_cosAlpha);
   fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
   fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
   fChain->SetBranchAddress("Reco_QQ_cosAlpha3D", Reco_QQ_cosAlpha3D, &b_Reco_QQ_cosAlpha3D);
   fChain->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen, &b_Reco_QQ_whichGen);
   fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
   fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
   fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
   fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
   fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
   fChain->SetBranchAddress("Reco_QQ_mupl_dxy_muonlessVtx", Reco_QQ_mupl_dxy_muonlessVtx, &b_Reco_QQ_mupl_dxy_muonlessVtx);
   fChain->SetBranchAddress("Reco_QQ_mumi_dxy_muonlessVtx", Reco_QQ_mumi_dxy_muonlessVtx, &b_Reco_QQ_mumi_dxy_muonlessVtx);
   fChain->SetBranchAddress("Reco_QQ_mupl_dxyErr_muonlessVtx", Reco_QQ_mupl_dxyErr_muonlessVtx, &b_Reco_QQ_mupl_dxyErr_muonlessVtx);
   fChain->SetBranchAddress("Reco_QQ_mumi_dxyErr_muonlessVtx", Reco_QQ_mumi_dxyErr_muonlessVtx, &b_Reco_QQ_mumi_dxyErr_muonlessVtx);
   fChain->SetBranchAddress("Reco_QQ_mupl_dz_muonlessVtx", Reco_QQ_mupl_dz_muonlessVtx, &b_Reco_QQ_mupl_dz_muonlessVtx);
   fChain->SetBranchAddress("Reco_QQ_mumi_dz_muonlessVtx", Reco_QQ_mumi_dz_muonlessVtx, &b_Reco_QQ_mumi_dz_muonlessVtx);
   fChain->SetBranchAddress("Reco_QQ_mupl_dzErr_muonlessVtx", Reco_QQ_mupl_dzErr_muonlessVtx, &b_Reco_QQ_mupl_dzErr_muonlessVtx);
   fChain->SetBranchAddress("Reco_QQ_mumi_dzErr_muonlessVtx", Reco_QQ_mumi_dzErr_muonlessVtx, &b_Reco_QQ_mumi_dzErr_muonlessVtx);
   fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
   fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
   fChain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
   fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
   fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
   fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
   fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
   fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
   fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
   fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
   fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
   fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
   fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
   fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
   fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
   fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
   fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
   fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
   fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
   fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
   fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
   fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
   fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
   fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
   fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
   fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
   fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
   fChain->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
   fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   fChain->SetBranchAddress("Gen_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
   fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
   fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
   fChain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
   fChain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
   fChain->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec, &b_Gen_QQ_whichRec);
   fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
   fChain->SetBranchAddress("Gen_mu_type", Gen_mu_type, &b_Gen_mu_type);
   fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
   fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
   fChain->SetBranchAddress("Gen_mu_whichRec", Gen_mu_whichRec, &b_Gen_mu_whichRec);
   fChain->SetBranchAddress("Onia2MuMuPAT", &Onia2MuMuPAT, &b_Onia2MuMuPAT);
   fChain->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
   fChain->SetBranchAddress("pHBHENoiseFilterResultProducer", &pHBHENoiseFilterResultProducer, &b_pHBHENoiseFilterResultProducer);
   fChain->SetBranchAddress("HBHENoiseFilterResult", &HBHENoiseFilterResult, &b_HBHENoiseFilterResult);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun1", &HBHENoiseFilterResultRun1, &b_HBHENoiseFilterResultRun1);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose, &b_HBHENoiseFilterResultRun2Loose);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &HBHENoiseFilterResultRun2Tight, &b_HBHENoiseFilterResultRun2Tight);
   fChain->SetBranchAddress("HBHEIsoNoiseFilterResult", &HBHEIsoNoiseFilterResult, &b_HBHEIsoNoiseFilterResult);
   fChain->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter, &b_pPAprimaryVertexFilter);
   fChain->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter, &b_pBeamScrapingFilter);
   fChain->SetBranchAddress("pVertexFilterCutG", &pVertexFilterCutG, &b_pVertexFilterCutG);
   fChain->SetBranchAddress("pVertexFilterCutGloose", &pVertexFilterCutGloose, &b_pVertexFilterCutGloose);
   fChain->SetBranchAddress("pVertexFilterCutGtight", &pVertexFilterCutGtight, &b_pVertexFilterCutGtight);
   fChain->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus, &b_pVertexFilterCutGplus);
   fChain->SetBranchAddress("pVertexFilterCutE", &pVertexFilterCutE, &b_pVertexFilterCutE);
   fChain->SetBranchAddress("pVertexFilterCutEandG", &pVertexFilterCutEandG, &b_pVertexFilterCutEandG);

   Notify();
}

Bool_t myTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Bool_t myTree::isGlobalMuonInAccept2019 (TLorentzVector* Muon)
{
  return ( fabs(Muon->Eta()) < 2.4 &&
	   ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
	    (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.47-1.89*fabs(Muon->Eta())) ||
	    (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.5)));
};

Bool_t myTree::areMuonsInAcceptance2019 (Int_t iRecoQQ)
{
  TLorentzVector* RecoQQmupl = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[iRecoQQ]);
  TLorentzVector* RecoQQmumi = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[iRecoQQ]);
  
  return ( isGlobalMuonInAccept2019(RecoQQmupl) && isGlobalMuonInAccept2019(RecoQQmumi) );
};

Bool_t myTree::areGenMuonsInAcceptance2019 (Int_t iGenQQ)
{
  TLorentzVector* GenQQmupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGenQQ]);
  TLorentzVector* GenQQmumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGenQQ]);

  return ( isGlobalMuonInAccept2019(GenQQmupl) && isGlobalMuonInAccept2019(GenQQmumi) );
};


Bool_t myTree::passQualityCuts2019 (Int_t iRecoQQ)
{
  int iMupl = Reco_QQ_mupl_idx[iRecoQQ];
  int iMumi = Reco_QQ_mumi_idx[iRecoQQ];

  if ( ! (Reco_mu_SelectionType[iMumi] & ((ULong64_t) pow(2,1))) ) return false; //muons have to be global muons
  if ( ! (Reco_mu_SelectionType[iMumi] & ((ULong64_t) pow(2,3))) ) return false; //muons have to be tracker muons
  if ( ! (Reco_mu_TMOneStaTight[iMumi] == 1) ) return false;
  if ( ! (Reco_mu_nTrkWMea[iMumi] > 5) ) return false;
  if ( ! (Reco_mu_nPixWMea[iMumi] > 0) ) return false;
  if ( ! (fabs(Reco_mu_dxy[iMumi]) < 0.3) ) return false;
  if ( ! (fabs(Reco_mu_dz[iMumi]) < 20.0) ) return false;

  if ( ! (Reco_mu_SelectionType[iMupl] & ((ULong64_t) pow(2,1))) ) return false; //muons have to be global muons
  if ( ! (Reco_mu_SelectionType[iMupl] & ((ULong64_t) pow(2,3))) ) return false; //muons have to be tracker muons
  if ( ! (Reco_mu_TMOneStaTight[iMupl] == 1) ) return false;
  if ( ! (Reco_mu_nTrkWMea[iMupl] > 5) ) return false;
  if ( ! (Reco_mu_nPixWMea[iMupl] > 0) ) return false;
  if ( ! (fabs(Reco_mu_dxy[iMupl]) < 0.3) ) return false;
  if ( ! (fabs(Reco_mu_dz[iMupl]) < 20.0) ) return false;

  if ( ! (Reco_QQ_VtxProb[iRecoQQ] > 0.01) ) return false;

  return true;
}

Bool_t myTree::isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit)
{
  Bool_t cond = true;
  cond = cond && ( (HLTriggers & ((ULong64_t) pow(2, TriggerBit))) == ((ULong64_t) pow(2, TriggerBit)) );
  cond = cond && ( (Reco_QQ_trig[iRecoQQ] & ((ULong64_t) pow(2, TriggerBit))) == ((ULong64_t) pow(2, TriggerBit)) );
  return cond;
}

#endif // #ifdef myTree_cxx
