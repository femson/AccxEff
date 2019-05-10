//This macro will create the efficiency histograms

#define myTree_cxx
#include "myTree.h"

Double_t ptbins [] = {3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 25., 30., 40., 50., 100.};
Double_t ybins []  = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

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
  /*
  TH2F* num_binned      = new TH2F("num_binned", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_binned->Sumw2();
  TH2F* num_pl1sig      = new TH2F("num_pl1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_pl1sig->Sumw2();
  TH2F* num_min1sig     = new TH2F("num_min1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_min1sig->Sumw2();
  TH2F* num_muid_sta    = new TH2F("num_muid_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_muid_sta->Sumw2();
  TH2F* num_muid        = new TH2F("num_muid", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_muid->Sumw2();
  TH2F* num_muid_pl1sig = new TH2F("num_muid_pl1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_muid_pl1sig->Sumw2();
  TH2F* num_muid_min1sig= new TH2F("num_muid_min1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_muid_min1sig->Sumw2();
  TH2F* num_sta         = new TH2F("num_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_sta->Sumw2();
  TH2F* num_sta_pl1sig  = new TH2F("num_sta_pl1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_sta_pl1sig->Sumw2();
  TH2F* num_sta_min1sig = new TH2F("num_sta_min1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_sta_min1sig->Sumw2();
  TH2F* num_trk_pl1sig  = new TH2F("num_trk_pl1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_trk_pl1sig->Sumw2();
  TH2F* num_trk_min1sig = new TH2F("num_trk_min1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); num_trk_min1sig->Sumw2();
  */

  Long64_t nentries = fChain->GetEntries();
  //nentries = 1000000;

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

	  /*
	    tnp = tag and probe
	    need a different tnp_weight function for each of the following:
	  */

	  /*
	  num_binned->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_pl1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_min1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_muid_sta->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_muid->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_muid_pl1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_muid_min1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_sta->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_sta_pl1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_sta_min1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_trk_pl1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  num_trk_min1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  */

	  num_pt->Fill(RecoQQ4mom->Pt());
	  if (RecoQQ4mom->Pt() > 6.5) num_y->Fill(fabs(RecoQQ4mom->Rapidity()));
	}
    }

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

  quotient_pt->Write("quotient_pt");
  quotient_y->Write("quotient_y");
  quotient_2D->Write("quotient_2D");

  deno_pt_y->Write("deno_2D");
  num_noweights->Write("num_2D_noweights");

  /*
  num_binned->Write("num_2D_binned");
  num_pl1sig->Write("num_2D_trg_pl1sig");
  num_min1sig->Write("num_2D_trg_min1sig");
  num_muid_sta->Write("num_2D_muid_sta");
  num_muid->Write("num_2D_muid");
  num_muid_pl1sig->Write("num_2D_muid_pl1sig");
  num_muid_min1sig->Write("num_2D_muid_min1sig");
  num_sta->Write("num_2D_sta");
  num_sta_pl1sig->Write("num_2D_sta_pl1sig");
  num_sta_min1sig->Write("num_2D_sta_min1sig");
  num_trk_pl1sig->Write("num_2D_trk_pl1sig");
  num_trk_min1sig->Write("num_2D_trk_min1sig");
  */

  fsave->Close();

}
