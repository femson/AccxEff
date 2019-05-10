// This macro separately computes Acceptance and Efficiency
// Data: /data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v1/HiForestAOD_ext_merged.root

#define myTree_cxx
#include "myTree.h"

Double_t ptbins [] = {3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 25., 30., 40., 50., 100.};
Double_t ybins []  = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};


// ----------------------
// ----- ACCEPTANCE -----
// Acceptance = Number of generated events(both muons pass ACC) / Number of generated events(total)

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
  nentries = 100000;

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
    } // end of iQQ loop
  }   // end of jentry loop

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

  TFile* fsave = new TFile("AccHists.root", "recreate");
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
} //end of acceptance function
 

// ----------------------
// ----- EFFICIENCY -----
// ----------------------

void myTree::Eff()
{
  int nptbins = sizeof(ptbins)/sizeof(double)-1;
  int nybins  = sizeof(ybins)/sizeof(double)-1;

  TH1F* deno_pt = new TH1F("deno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbins, ptbins); deno_pt->Sumw2();
  TH1F* num_pt  = new TH1F("num_pt", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptbins, ptbins); num_pt->Sumw2();
  TH1F* deno_y  = new TH1F("deno_y", "N_{gen} vs y; y; N_{gen}", nybins, ybins); deno_y->Sumw2();
  TH1F* num_y   = new TH1F("num_y", "N_{reco} vs y; y; N_{reco}", nybins, ybins); num_y->Sumw2();

  TH2F* deno_pt_y       = new TH2F("deno_pt_y", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); deno_pt_y->Sumw2();
  TH2F* num_noweights   = new TH2F("num_noweights", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_noweights->Sumw2();
 
  Long64_t nentries = fChain->GetEntries();
  nentries = 1000000;

  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
      if (jentry%10000==0) cout<<"[INFO] "<<jentry*100/nentries<<"% finished"<<'\r'<<flush;

      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      for (int iQQ = 0; iQQ < Gen_QQ_size; iQQ++)
	{
	  /*
	  cout << "Event number: " << eventNb << endl;
	  cout << "Gen_QQ_size = " << Gen_QQ_size << "\t Gen_mu_size = " << Gen_mu_size << endl;
	  cout << "Reco_QQ_size = " << Reco_QQ_size << "\t Reco_mu_size = " << Reco_mu_size << endl;
	  */

	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  jpsi_m   = GenQQ4mom->M();
	  jpsi_pt  = GenQQ4mom->Pt();
	  jpsi_rap = GenQQ4mom->Rapidity();

	  if (jpsi_pt < 3 || jpsi_pt > 100) continue;
	  if (fabs(jpsi_rap) >= 2.4) continue;
	  if (!areGenMuonsInAcceptance2019(iQQ)) continue;

	  deno_pt_y->Fill(jpsi_rap, jpsi_pt);
	  deno_pt->Fill(jpsi_pt);

	  if (jpsi_pt > 6.5) deno_y->Fill(fabs(jpsi_rap));

	  int whichRec = Gen_QQ_whichRec[iQQ];
	  if (whichRec < 0) continue;
	  // cout << "Matched, Gen_QQ_whichRec[iQQ] = " << Gen_QQ_whichRec[iQQ] << endl;

	  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(whichRec);

	  if (!areMuonsInAcceptance2019(whichRec)) continue;
	  if (!passQualityCuts2019(whichRec)) continue;
	  if (!isTriggerMatch(whichRec, triggerIndex_PP)) continue;	
	  if (Reco_QQ_sign[whichRec] != 0) continue;
	  if (RecoQQ4mom->Pt() < 3 || RecoQQ4mom->Pt() > 100) continue;
	  if (fabs(RecoQQ4mom->Rapidity()) >= 2.4) continue;
	  if (RecoQQ4mom->M() < 2.6 || RecoQQ4mom->M() > 3.5) continue;


	  //if (isPbPb && !(pprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter2Th4)) continue;
	  if (!isPbPb && !(pPAprimaryVertexFilter && pBeamScrapingFilter)) continue;

	  //cout << "working" << endl;

	  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[whichRec]);
	  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[whichRec]);

	  tnp_weight = 1.0;
	  num_noweights->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt());

	  num_pt->Fill(RecoQQ4mom->Pt());
	  if (RecoQQ4mom->Pt() > 6.5) num_y->Fill(fabs(RecoQQ4mom->Rapidity()));
	}    //end of iQQ loop
    }        //end of jentry loop

  auto quotient_pt = new TH1F(*num_pt);
  quotient_pt->Divide(deno_pt);
  auto quotient_y = new TH1F(*num_y);
  quotient_y->Divide(deno_y);
  auto quotient_2D = new TH2F(*num_noweights);
  quotient_2D->Divide(deno_pt_y);

  TFile* fsave = new TFile ("EffHists.root", "recreate");
  deno_pt->Write("deno_pt");
  num_pt->Write("num_pt");
  deno_y->Write("deno_y");
  num_y->Write("num_y");
  deno_pt_y->Write("deno_2D");
  num_noweights->Write("num_2D_noweights");

  quotient_pt->Write("quotient_pt");
  quotient_y->Write("quotient_y");
  quotient_2D->Write("quotient_2D");

  fsave->Close();
} //end of Efficiency function
