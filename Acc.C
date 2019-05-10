// This will become the acceptance macro

// rl /data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v1/HiForestAOD_ext_merged.root
/*
     Acceptance = Number of generated events(both muons pass ACC) / Number of generated events(total)
*/

#define myTree_cxx
#include "myTree.h"

Double_t ptbins [] = {3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 25., 30., 40., 50., 100.};
Double_t ybins []  = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

void myTree::Acc()
{
  int nptbins = sizeof(ptbins)/sizeof(double)-1;     //sizeof(ptbins) = 8 * #elements in array; sizeof(double) = 8
  int nybins  = sizeof(ybins)/sizeof(double)-1;

  TH1F* deno_pt = new TH1F("deno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbins, ptbins); deno_pt->Sumw2();
  TH1F* num_pt  = new TH1F("num_pt", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptbins, ptbins); num_pt->Sumw2();
  TH1F* deno_y  = new TH1F("deno_y", "N_{gen} vs y; y; N_{gen}", nybins, ybins); deno_y->Sumw2();
  TH1F* num_y   = new TH1F("num_y", "N_{reco} vs y; y; N_{reco}", nybins, ybins); num_y->Sumw2();

  TH2F* deno_pt_y = new TH2F("deno_pt_y", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); deno_pt_y->Sumw2();
  TH2F* num_pt_y = new TH2F("num_pt_y", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_pt_y->Sumw2();

  gStyle->SetOptStat(0);

  Long64_t nentries = fChain->GetEntries();
  //nentries = 100000;

  Long64_t nbytes = 0, nb = 0;
  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    
    if (jentry%10000==0) cout<<"[INFO] "<<jentry*100/nentries<<"% finished"<<'\r'<<flush;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    for(int iQQ = 0; iQQ < Gen_QQ_size; iQQ++){
      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
      jpsi_m = GenQQ4mom->M();
      jpsi_pt = GenQQ4mom->Pt();
      jpsi_rap = GenQQ4mom->Rapidity();

      if (jpsi_pt < 3 || jpsi_pt > 100) continue;
      if (fabs(jpsi_rap) > 2.4) continue;
      deno_pt_y->Fill(jpsi_rap, jpsi_pt);
      deno_pt->Fill(jpsi_pt);
      
      if (jpsi_pt >= 6.5) deno_y->Fill(fabs(jpsi_rap));

      if (!areGenMuonsInAcceptance2019(iQQ)) continue;

      num_pt_y->Fill(jpsi_rap, jpsi_pt);
      num_pt->Fill(jpsi_pt);

      if (jpsi_pt >= 6.5) num_y->Fill(fabs(jpsi_rap)); 
    }
  }

  gStyle->SetPalette(kBird);

  TCanvas* c = new TCanvas();
  c->Divide(2,1);
  c->cd(1);
  num_pt->Draw();
  c->cd(2);
  deno_pt->Draw();

  TCanvas* c_2D = new TCanvas();
  c_2D->Divide(2,1);
  c_2D->cd(1);
  num_pt_y->DrawClone("colz");
  // c_2D->cd(2);
  //num_pt_y->DrawClone("colz");
  c_2D->cd(2);
  deno_pt_y->DrawClone("colz");
  //c_2D->cd(4);
  //deno_pt_y->DrawClone("colz");


  auto AccQuot_pt = new TH1F(*num_pt);
  AccQuot_pt->Divide(deno_pt);
  auto AccQuot_y = new TH1F(*num_y);
  AccQuot_y->Divide(deno_y);
  auto AccQuot_2D = new TH2F(*num_pt_y);
  AccQuot_2D->Divide(deno_pt_y);

  TFile* fsave = new TFile("test.root", "recreate");
  deno_pt->Write("deno_pt");
  num_pt->Write("num_pt");
  deno_y->Write("deno_y");
  num_y->Write("num_y");
  deno_pt_y->Write("deno_pt_y");
  num_pt_y->Write("num_pt_y");

  AccQuot_pt->Write("AccQuot_pt");
  AccQuot_y->Write("AccQuot_y");
  AccQuot_2D->Write("AccQuot_2D");

  fsave->Close();

}
