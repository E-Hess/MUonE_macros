#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TSystem.h"
#include "FairMCEventHeader.h"
#include "MUonETrack.h"
#include <cmath> // for std::abs and std::fabs
// #include "MUonETrackerStub.h"
#include <TAxis.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <fstream>
#include "TChain.h"
#include "FairRootManager.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
// #include "MUonERecoOutputBase.h"
#include "ElasticState.h"

using namespace std;
using namespace MuE;

// older code has reco + swapped angles

// CHANGED LINKING
// " new per-event variables which give you: https://indico.cern.ch/event/1505154/contributions/6335583/attachments/3001752/5290192/2025.01.23_GAbbiendi_FM-MCtruth-upd.pdf
// the interactionID of the produced particle (charged or neutral) with maximum energy (apart from the beam muon)
// the interactionID of the produced charged particle with maximum energy (apart from the beam muon)
// These two can be identical or not, depending on what happened. For example:
// -+ if there was an elastic scattering and the scattered electron is the most energetic particle both will be defined by it
// -+ if there was a radiated photon that at some point in the detector converts in a e+e- pair and those are the most energetic particles
// you will have the first variable saying ‘bremsstrahlung’ and the second ‘pair production’
// If nothing happens and there is just the beam muon the two variables are left at the value = 0"

// plots if reco but no swap
// plots if mc truth

void EventDisplay(TFile *file, vector<double> z_inmu, vector<double> y_inmu, vector<double> z_e, vector<double> y_e, vector<double> z_mu, vector<double> y_mu, int eventNumber)
{

    // cout << "Incoming Muon z positions: ";
    // for (const auto &z : z_inmu)
    // {
    //     cout << z << " ";
    // }
    // cout << endl;

    vector<double> xModules = {828.02, 899.92, 929.22, 1001.12};
    vector<double> yModules = {831.87, 903.77, 933.07, 1004, 97};
    vector<double> uvModules = {865.36, 866.62, 966.45, 967.82};

    double labelSize = 0.05;  // Adjust for larger labels
    double titleSize = 0.06;  // Adjust for larger titles
    double tickLength = 0.02; // Adjust for longer ticks

    double targetPosition = 912.7;
    double ll = 5.0;

    TCanvas *c1 = new TCanvas(Form("c1_event_%d", eventNumber), "Event Display", 800, 300);
    c1->Divide(1, 3);
    c1->cd(1);

    TH2D *h_inmu = new TH2D("h_inmu", Form("Event %d;Z;X plane", eventNumber), 100, 820, 1010, 100, -ll, ll);
    TH2D *h_outmu = new TH2D("h_outmu", Form("Event %d;Z;X plane", eventNumber), 100, 820, 1010, 100, -ll, ll);
    TH2D *h_oute = new TH2D("h_oute", Form("Event %d;Z;X plane", eventNumber), 100, 820, 1010, 100, -ll, ll);

    for (size_t i = 0; i < z_inmu.size(); ++i)
    {
        for (double ZofXmodules : xModules)
        {
            if (abs(z_inmu[i] - ZofXmodules) < 1.5)
            {
                h_inmu->Fill(z_inmu[i], y_inmu[i]);
            }
        }
    }

    h_inmu->SetStats(0);
    h_inmu->SetMarkerStyle(8);
    h_inmu->SetMarkerColor(4); // blue
    h_inmu->GetXaxis()->SetLabelSize(labelSize);
    h_inmu->GetXaxis()->SetTitleSize(titleSize);
    h_inmu->GetXaxis()->SetTickLength(tickLength);
    h_inmu->GetXaxis()->SetTitleOffset(1.2);
    h_inmu->GetYaxis()->SetLabelSize(labelSize);
    h_inmu->GetYaxis()->SetTitleSize(titleSize);
    h_inmu->GetYaxis()->SetTickLength(tickLength);
    h_inmu->Draw("SCAT");

    for (size_t i = 0; i < z_mu.size(); ++i)
    {
        for (double ZofXmodules : xModules)
        {
            if (abs(z_mu[i] - ZofXmodules) < 1.5)
            {
                h_outmu->Fill(z_mu[i], y_mu[i]);
            }
        }
    }

    h_outmu->SetStats(0);
    h_outmu->SetMarkerStyle(21);
    h_outmu->SetMarkerColor(2);
    h_outmu->Draw("SCAT SAME");

    for (size_t i = 0; i < z_e.size(); ++i)
    {
        for (double ZofXmodules : xModules)
        {
            if (abs(z_e[i] - ZofXmodules) < 1.5)
            {
                h_oute->Fill(z_e[i], y_e[i]);
            }
        }
    }

    h_oute->SetStats(0);
    h_oute->SetMarkerStyle(22);
    h_oute->SetMarkerColor(3);
    h_oute->Draw("SCAT SAME");

    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend1->AddEntry(h_inmu, "Incoming muon", "p");
    legend1->AddEntry(h_outmu, "Outgoing muon", "p");
    legend1->AddEntry(h_oute, "Outgoing electron", "p");
    legend1->Draw();

    // Add vertical lines for xModules
    for (double xModule : xModules)
    {
        TLine *line = new TLine(xModule, -5.0, xModule, 5.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(2); // Dashed line
        line->Draw();
    }

    TLine *targetLine = new TLine(targetPosition, -3.0, targetPosition, 3.0);
    targetLine->SetLineColor(kBlack);
    targetLine->SetLineStyle(1); // Solid line
    targetLine->SetLineWidth(3);
    targetLine->Draw();

    c1->Update();
    c1->Modified();

    c1->cd(2);

    TH2D *h_inmu2 = new TH2D("h_inmu2", Form("Event %d;Z;Y plane", eventNumber), 100, 820, 1010, 100, -ll, ll);
    TH2D *h_outmu2 = new TH2D("h_outmu2", Form("Event %d;Z;Y plane", eventNumber), 100, 820, 1010, 100, -ll, ll);
    TH2D *h_oute2 = new TH2D("h_oute2", Form("Event %d;Z;Y", eventNumber), 100, 820, 1010, 100, -ll, ll);

    for (size_t i = 0; i < z_inmu.size(); ++i)
    {
        for (double ZofYmodules : yModules)
        {
            if (abs(z_inmu[i] - ZofYmodules) < 1.5)
            {
                h_inmu2->Fill(z_inmu[i], y_inmu[i]);
            }
        }
    }

    h_inmu2->SetStats(0);
    h_inmu2->SetMarkerStyle(8);
    h_inmu2->SetMarkerColor(4); // blue
    h_inmu2->GetXaxis()->SetLabelSize(labelSize);
    h_inmu2->GetXaxis()->SetTitleSize(titleSize);
    h_inmu2->GetXaxis()->SetTickLength(tickLength);
    h_inmu2->GetXaxis()->SetTitleOffset(1.2);
    h_inmu2->GetYaxis()->SetLabelSize(labelSize);
    h_inmu2->GetYaxis()->SetTitleSize(titleSize);
    h_inmu2->GetYaxis()->SetTickLength(tickLength);
    h_inmu2->Draw("SCAT");

    for (size_t i = 0; i < z_mu.size(); ++i)
    {
        for (double ZofYmodules : yModules)
        {
            if (abs(z_mu[i] - ZofYmodules) < 1.5)
            {
                h_outmu2->Fill(z_mu[i], y_mu[i]);
            }
        }
    }

    h_outmu2->SetStats(0);
    h_outmu2->SetMarkerStyle(21);
    h_outmu2->SetMarkerColor(2);
    h_outmu2->Draw("SCAT SAME");

    for (size_t i = 0; i < z_e.size(); ++i)
    {
        for (double ZofYmodules : yModules)
        {
            if (abs(z_e[i] - ZofYmodules) < 1.5)
            {
                h_oute2->Fill(z_e[i], y_e[i]);
            }
        }
    }

    h_oute2->SetStats(0);
    h_oute2->SetMarkerStyle(22);
    h_oute2->SetMarkerColor(3);
    h_oute2->Draw("SCAT SAME");

    // c1->BuildLegend();
    TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2->AddEntry(h_inmu2, "Incoming muon", "p");
    legend2->AddEntry(h_outmu2, "Outgoing muon", "p");
    legend2->AddEntry(h_oute2, "Outgoing electron", "p");
    legend2->Draw();

    // Add vertical lines for yModules
    for (double yModule : yModules)
    {
        TLine *line = new TLine(yModule, -5.0, yModule, 5.0);
        line->SetLineColor(kGreen);
        line->SetLineStyle(2); // Dashed line
        line->Draw();
    }
    TLine *targetLine2 = new TLine(targetPosition, -3.0, targetPosition, 3.0);
    targetLine2->SetLineColor(kBlack);
    targetLine2->SetLineStyle(1); // Solid line
    targetLine2->SetLineWidth(3);
    targetLine2->Draw();

    c1->Update();
    c1->Modified();

    c1->cd(3);
    TH2D *h_inmu3 = new TH2D("h_inmu3", Form("Event %d;Z;UV plane", eventNumber), 100, 820, 1010, 100, -ll, ll);
    TH2D *h_outmu3 = new TH2D("h_outmu3", Form("Event %d;Z-axis;UV plane", eventNumber), 100, 820, 1010, 100, -ll, ll);
    TH2D *h_oute3 = new TH2D("h_oute3", Form("Event %d;Z;UV plane", eventNumber), 100, 820, 1010, 100, -ll, ll);

    for (size_t i = 0; i < z_inmu.size(); ++i)
    {
        for (double ZofUVmodules : uvModules)
        {
            if (abs(z_inmu[i] - ZofUVmodules) < 1.0)
            {
                h_inmu3->Fill(z_inmu[i], y_inmu[i]);
            }
        }
    }

    h_inmu3->SetStats(0);
    h_inmu3->SetMarkerStyle(8);
    h_inmu3->SetMarkerColor(4); // blue
    h_inmu3->GetXaxis()->SetLabelSize(labelSize);
    h_inmu3->GetXaxis()->SetTitleSize(titleSize);
    h_inmu3->GetXaxis()->SetTickLength(tickLength);
    h_inmu3->GetXaxis()->SetTitleOffset(1.2);
    h_inmu3->GetYaxis()->SetLabelSize(labelSize);
    h_inmu3->GetYaxis()->SetTitleSize(titleSize);
    h_inmu3->GetYaxis()->SetTickLength(tickLength);
    h_inmu3->Draw("SCAT");

    for (size_t i = 0; i < z_mu.size(); ++i)
    {
        for (double ZofUVmodules : uvModules)
        {
            if (abs(z_mu[i] - ZofUVmodules) < 1.0)
            {
                h_outmu3->Fill(z_mu[i], y_mu[i]);
            }
        }
    }

    h_outmu3->SetStats(0);
    h_outmu3->SetMarkerStyle(21);
    h_outmu3->SetMarkerColor(2);
    h_outmu3->Draw("SCAT SAME");

    for (size_t i = 0; i < z_e.size(); ++i)
    {
        for (double ZofUVmodules : uvModules)
        {
            if (abs(z_e[i] - ZofUVmodules) < 1.0)
            {
                h_oute3->Fill(z_e[i], y_e[i]);
            }
        }
    }

    h_oute3->SetStats(0);
    h_oute3->SetMarkerStyle(22);
    h_oute3->SetMarkerColor(3);
    h_oute3->Draw("SCAT SAME");

    // c1->BuildLegend();
    TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend3->AddEntry(h_inmu3, "Incoming muon", "p");
    legend3->AddEntry(h_outmu3, "Outgoing muon", "p");
    legend3->AddEntry(h_oute3, "Outgoing electron", "p");
    legend3->Draw();

    // Add vertical lines for uvModules
    for (double uvModule : uvModules)
    {
        TLine *line = new TLine(uvModule, -5.0, uvModule, 5.0);
        line->SetLineColor(kBlue);
        line->SetLineStyle(2); // Dashed line
        line->Draw();
    }
    TLine *targetLine3 = new TLine(targetPosition, -3.0, targetPosition, 3.0);
    targetLine3->SetLineColor(kBlack);
    targetLine3->SetLineStyle(1); // Solid line
    targetLine3->SetLineWidth(3);
    targetLine3->Draw();

    file->cd();
    c1->Write();
    delete h_inmu;
    delete h_outmu;
    delete h_oute;
    delete h_inmu2;
    delete h_outmu2;
    delete h_oute2;
    delete h_inmu3;
    delete h_outmu3;
    delete h_oute3;

    c1->Close();
}

void runMCAnalysis_NewLinking_allClass_swap()
{
    // Create and open a log file
    ofstream logFile("MCAnalisys_total_2AddCond_NewLink_allClass_swap_1.4.25.log");

    // Redirect std::cout to the log file
    streambuf *coutBuf = cout.rdbuf();
    cout.rdbuf(logFile.rdbuf());

    cout << "Starting the analysis..." << endl;

    // Disable interactive graphics updates
    gROOT->SetBatch(kTRUE);

    TString filepath = "/eos/user/e/ehess/fair_install_v16_RecoEmma_13.3.25/FairRootInst/share/MUonE/macros/OutputFiles/";
    TString fullpath1 = filepath + "MCMinimumBias_skimv2_SimDigiReco_NewLinking.root"; // cond:last modules

    TFile *f = TFile::Open(fullpath1);
    TTree *cbmsim = (TTree *)f->Get("cbmsim");
    cout << "File succesfully opened." << endl;

    TClonesArray *MCTrack = nullptr;
    TClonesArray *TrackerStubs = nullptr;
    MUonERecoOutputAnalysis *ReconstructionOutput = nullptr;

    cbmsim->SetBranchAddress("MCTrack", &MCTrack);
    cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
    cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

    Int_t Nreco = 0;                          // Number of reconstructed events
    Long64_t nEntries = cbmsim->GetEntries(); // Number of all the skimmed events

    TH1F *h_mctPdgCode = new TH1F("h_mctPdgCode", "PdgCode for all particles (tracks)", 60, -30, 30);
    TH1F *h_mctPdgCodeOutE = new TH1F("h_mctPdgCodeOutE", "PdgCode of the outgoing electron", 60, -30, 30);
    TH1F *h_mctInteractionID = new TH1F("h_mctInteractionID", "InteractionID for all particles (tracks)", 30, 0, 30);
    TH1F *h_mctInteractionIDoutE = new TH1F("h_mctInteractionIDoutE", "InteractionID of the outgoing electron", 30, 0, 30);
    TH1F *h_mctPdgCodeOther = new TH1F("h_mctPdgCodeOther", "PdgCode for other processes", 60, -300, 300);
    TH1F *h_mctInteractionIDOther = new TH1F("h_mctInteractionIDOther", "InteractionID for other processes", 40, 0, 40);
    TH1F *h_mctPdgCodeOther_afterCuts = new TH1F("h_mctPdgCodeOther_afterCuts", "PdgCode for other processes", 60, -300, 300);
    TH1F *h_mctInteractionIDOther_afterCuts = new TH1F("h_mctInteractionIDOther_afterCuts", "InteractionID for other processes", 40, 0, 40);
    TH1F *h_nStubs_1 = new TH1F("h_nStubs_1", "nStubs on S1", 30, 0, 30);

    TH1F *h_openingAngles = new TH1F("h_openingAngles", "Opening Angles", 100, 0, 0.1);
    TH1F *h_openingAngles5 = new TH1F("h_openingAngles5", "Opening Angles 5", 100, 0, 0.1);
    TH1F *h_openingAngles5phot = new TH1F("h_openingAngles5phot", "Opening Angles 5 phot", 100, 0, 0.1);
    TH1F *h_openingAngles9 = new TH1F("h_openingAngles9", "Opening Angles 9", 100, 0, 0.1);
    TH1F *h_openingAnglesOther = new TH1F("h_openingAnglesOther", "Opening Angles", 100, 0, 0.1);
    TH1F *h_openingAnglesNoMu = new TH1F("h_openingAnglesNoMu", "Opening Angles", 100, 0, 0.1);

    // 1. CUT: theta_max < 32 mrad && theta_min >=0.2 mrad
    // 2. CUT: 1. cut + N_stubs(S1)<=14
    // 3. CUT: 2. cut + |modified_Acoplanarity|<0.4 (or <1)
    // 4. CUT: 3. cut + chi2_vertex (KF) <= 20
    // 5. CUT: 4. cut + |z_vertex - z_target| <= 3 cm  (or 2cm depending on the target thickness)

    const int nCuts = 2;
    TH2F *h2_ThetaMuVsThetaE_badReco[nCuts];
    TH2F *h2_ThetaMuVsThetaE_notsameHardestReco[nCuts];
    TH2F *h2_ThetaMuVsThetaE_ppPhot[nCuts];
    TH2F *h2_ThetaMuVsThetaE_notSamePrimaryHardest[nCuts];
    TH2F *h2_ThetaMuVsThetaE_pp[nCuts];
    TH2F *h2_ThetaMuVsThetaE_signal[nCuts];
    TH2F *h2_ThetaMuVsThetaE_other[nCuts];
    TH2F *h2_ThetaMuVsThetaE_nuclear[nCuts];
    TH2F *h2_ThetaMuVsThetaE_decay[nCuts];
    TH2F *h2_ThetaMuVsThetaE_bremst[nCuts];

    TH1F *h_bvZposFit_badReco[nCuts];
    TH1F *h_bvZposFit_notsameHardestReco[nCuts];
    TH1F *h_bvZposFit_ppPhot[nCuts];
    TH1F *h_bvZposFit_notSamePrimaryHardest[nCuts];
    TH1F *h_bvZposFit_pp[nCuts];
    TH1F *h_bvZposFit_signal[nCuts];
    TH1F *h_bvZposFit_other[nCuts];

    for (int i = 0; i < nCuts; ++i)
    {
        h2_ThetaMuVsThetaE_badReco[i] = new TH2F(Form("h2_ThetaMuVsThetaE_badReco_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_badReco_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);
        h2_ThetaMuVsThetaE_notsameHardestReco[i] = new TH2F(Form("h2_ThetaMuVsThetaE_notsameHardestReco_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_notsameHardestReco_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);
        h2_ThetaMuVsThetaE_ppPhot[i] = new TH2F(Form("h2_ThetaMuVsThetaE_ppPhot_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_ppPhot_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);
        h2_ThetaMuVsThetaE_notSamePrimaryHardest[i] = new TH2F(Form("h2_ThetaMuVsThetaE_notSamePrimaryHardest_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_notSamePrimaryHardest_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);
        h2_ThetaMuVsThetaE_pp[i] = new TH2F(Form("h2_ThetaMuVsThetaE_pp_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_pp_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);
        h2_ThetaMuVsThetaE_signal[i] = new TH2F(Form("h2_ThetaMuVsThetaE_signal_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_signal_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);
        h2_ThetaMuVsThetaE_other[i] = new TH2F(Form("h2_ThetaMuVsThetaE_other_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_other_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);
        h2_ThetaMuVsThetaE_nuclear[i] = new TH2F(Form("h2_ThetaMuVsThetaE_nuclear_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_nuclear_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);
        h2_ThetaMuVsThetaE_decay[i] = new TH2F(Form("h2_ThetaMuVsThetaE_decay_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_decay_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);
        h2_ThetaMuVsThetaE_bremst[i] = new TH2F(Form("h2_ThetaMuVsThetaE_bremst_cut%d", i + 1), Form("h2_ThetaMuVsThetaE_bremst_cut%d", i + 1), 1000, 0, 0.032, 200, 0, 0.006);

        h_bvZposFit_badReco[i] = new TH1F(Form("h_bvZposFit_badReco_cut%d", i + 1), Form("h_bvZposFit_notsameHardestReco_cut%d", i + 1), 200, 850, 1050);
        h_bvZposFit_notsameHardestReco[i] = new TH1F(Form("h_bvZposFit_notsameHardestReco_cut%d", i + 1), Form("h_bvZposFit_notsameHardestReco_cut%d", i + 1), 200, 850, 1050);
        h_bvZposFit_ppPhot[i] = new TH1F(Form("h_bvZposFit_ppPhot_cut%d", i + 1), Form("h_bvZposFit_ppPhot_cut%d", i + 1), 200, 850, 1050);
        h_bvZposFit_notSamePrimaryHardest[i] = new TH1F(Form("h_bvZposFit_notSamePrimaryHardest_cut%d", i + 1), Form("h_bvZposFit_notSamePrimaryHardest_cut%d", i + 1), 200, 850, 1050);
        h_bvZposFit_pp[i] = new TH1F(Form("h_bvZposFit_pp_cut%d", i + 1), Form("h_bvZposFit_pp_cut%d", i + 1), 200, 850, 1050);
        h_bvZposFit_signal[i] = new TH1F(Form("h_bvZposFit_signal_cut%d", i + 1), Form("h_bvZposFit_signal_cut%d", i + 1), 200, 850, 1050);
        h_bvZposFit_other[i] = new TH1F(Form("h_bvZposFit_other_cut%d", i + 1), Form("h_bvZposFit_other_cut%d", i + 1), 200, 850, 1050);
    }

    vector<int> nMCTracks;
    int thetaMuminPPmu = 0;
    int thetaMuminPPphot = 0;
    int thetaMuminSignal = 0;
    int thetaMuminOther = 0;
    int thetaEmaxPPmu = 0;
    int thetaEmaxPPphot = 0;
    int thetaEmaxSignal = 0;
    int thetaEmaxOther = 0;

    int intE0 = 0;
    int bothMu = 0;
    int notInMu = 0;
    int other = 0;
    vector<int> pdgE0;
    int swapped = 0;
    int ppp = 0;
    int pair = 0;
    int signal = 0;
    int notSameHardestReco = 0;
    int notSameHardestPrimary = 0;
    int notTrack = 0, oneTrack = 0;
    int badRecoEvents = 0;
    int afterCuts = 0;
    vector<int> parentPDGpair;
    vector<double> xModules = {828.02, 899.92, 929.22, 1001.12};
    vector<double> yModules = {831.87, 903.77, 933.07, 1004, 97};
    vector<double> uvModules = {865.36, 866.62, 966.45, 967.82};
    vector<string> cuts = {"No cuts", "All cuts"};

    cout << "Number of events in file: " << nEntries << endl;

    // TFile *fEventDisplay = new TFile("EventDisplay_2AddCond_NewLinking_28.3.25.root", "RECREATE");

    for (Long64_t i = 0; i < nEntries; ++i) // nEntries
    {
        cbmsim->GetEntry(i);
        if (i % 50000 == 0)
            cout << "Event: " << i << endl;
        if (ReconstructionOutput->isReconstructed() == 0) // taking only reconstructed events (they have a vertex)
            continue;

        Nreco++;

        auto best_vrtx = ReconstructionOutput->bestVertex(); // MUonERecoOutputVertexAnalysis best_vrtx = ReconstructionOutput->bestVertex();
        auto oute = best_vrtx.outgoingElectron();            // MUonERecoOutputTrackAnalysis oute = best_vrtx.outgoingElectron();
        auto outmu = best_vrtx.outgoingMuon();               // MUonERecoOutputTrackAnalysis outmu = best_vrtx.outgoingMuon();
        auto inmu = best_vrtx.incomingMuon();                // MUonERecoOutputTrackAnalysis inmu = best_vrtx.incomingMuon();

        Int_t linkIDoutE = oute.linkedTrackID();
        Int_t linkIDoutMu = outmu.linkedTrackID();
        Int_t linkIDInMu = inmu.linkedTrackID();

        double chi2vert = best_vrtx.chi2perDegreeOfFreedom();
        int processIDE = oute.processIDofLinkedTrack();
        int processIDoutMu = outmu.processIDofLinkedTrack();
        double ZposFit = best_vrtx.zPositionFit();
        double modAcop = best_vrtx.modifiedAcoplanarity();
        double thetaE = best_vrtx.electronTheta();
        double thetamu = best_vrtx.muonTheta();

        auto hits = ReconstructionOutput->reconstructedHits();

        int recoHitsS1 = 0;
        for (int h = 0; h < hits.size(); h++)
        {
            if (hits.at(h).stationID() == 1)
                recoHitsS1++;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////// MOTHER ID, PDG, INTERACTION ID ///////////////////////

        int nTracks = MCTrack->GetEntries();
        nMCTracks.push_back(nTracks);

        Double_t Emax = 0.;
        Double_t pmax = 0.;
        const MUonETrack *primaryTrack = nullptr;
        const MUonETrack *chargedTrack = nullptr;

        Int_t ProcessID_primary = -1;
        Int_t ProcessID_hardestChargedParticle = -1;
        Int_t PDG_primary = 0;
        Int_t PDG_hardestChargedParticle = 0;

        if (nTracks < 2) // if reco event has only one track
        {
            cout << "----------------------------------------------------------" << endl;
            cout << "One track for reco event: " << i << endl;
            cout << "----------------------------------------------------------" << endl;
            oneTrack++;
            continue;
        }

        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
        {
            afterCuts++;
        }

        for (int n = 0; n < nTracks; n++) // Looping over all MC tracks
        {

            const auto MCTr = static_cast<const MUonETrack *>(MCTrack->At(n));

            int pdg = MCTr->pdgCode();
            int intID = MCTr->interactionID();
            int motherID = MCTr->motherID();
            // h_mctPdgCode->Fill(pdg);
            // h_mctInteractionID->Fill(intID);

            if (n < 1)
                continue; // skip incoming muon

            // New linking
            Double_t p = MCTr->p();
            Double_t E = MCTr->energy() - MCTr->mass(); // cout << "Impuls: " << p << endl; cout << "Energy: " << E << endl;
            Double_t charge = 0;

            auto pdgParticle = TDatabasePDG::Instance()->GetParticle(MCTr->pdgCode());
            if (pdgParticle)
                charge = pdgParticle->Charge(); // else charge = 0.0 because nucleai dont have charge so we initialized it to 0.

            if (E > Emax) // Update primary track (highest energy)
            {
                Emax = E;
                primaryTrack = MCTr;
            }

            if (charge != 0 && p > pmax) // Update hardest charged track (highest momentum among charged tracks)
            {
                pmax = p;
                chargedTrack = MCTr;
            }
        }

        if (primaryTrack && chargedTrack) // if primaryTrack && chargedTrack are found
        {
            // Use the found tracks' interaction IDs
            ProcessID_primary = primaryTrack->interactionID();
            ProcessID_hardestChargedParticle = chargedTrack->interactionID();
            PDG_primary = primaryTrack->pdgCode();
            PDG_hardestChargedParticle = chargedTrack->pdgCode();

            if (processIDE == 0) //processIDoutMu > 0 &&
            {
                swapped += 1;
                // Swap the particles
                swap(oute, outmu);
                // swap(linkIDoutE, linkIDoutMu);
                swap(thetaE, thetamu);
            }

            if (linkIDoutE < 0) // bad reconstruction
            {
                badRecoEvents++;

                if (ppp < 10) // for event display
                {
                    ppp++;
                    vector<double> z_inmu;
                    vector<double> y_inmu;
                    auto hits_inmu = inmu.hits();
                    for (const auto &hit : hits_inmu)
                    {
                        z_inmu.push_back(hit.z());
                        y_inmu.push_back(hit.positionPerpendicular());
                    }

                    vector<double> z_e;
                    vector<double> y_e;
                    auto hits_e = oute.hits();
                    for (const auto &hit : hits_e)
                    {
                        z_e.push_back(hit.z());
                        y_e.push_back(hit.positionPerpendicular());
                    }

                    vector<double> z_mu;
                    vector<double> y_mu;
                    auto hits_mu = outmu.hits();
                    for (const auto &hit : hits_mu)
                    {
                        z_mu.push_back(hit.z());
                        y_mu.push_back(hit.positionPerpendicular());
                    }

                    // cout << "Event " << i << ": Incoming Muon Hits: " << z_inmu.size() << endl;
                    // cout << "Event " << i << ": Outgoing Electron Hits: " << z_e.size() << endl;
                    // cout << "Event " << i << ": Outgoing Muon Hits: " << z_mu.size() << endl;

                    // EventDisplay(fEventDisplay, z_inmu, y_inmu, z_e, y_e, z_mu, y_mu, i);
                } // end of event display

                if (ProcessID_primary == ProcessID_hardestChargedParticle)
                {

                    if (ProcessID_primary == 5)
                    {
                        h2_ThetaMuVsThetaE_pp[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_pp[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_pp[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_pp[1]->Fill(ZposFit);
                        }
                    }

                    else if (ProcessID_primary == 9)
                    {
                        h2_ThetaMuVsThetaE_signal[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_signal[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_signal[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_signal[1]->Fill(ZposFit);
                        }
                    }
                    else if (ProcessID_primary == 25)
                    {

                        h2_ThetaMuVsThetaE_nuclear[0]->Fill(thetaE, thetamu);
                        // h_bvZposFit_nuclear[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_nuclear[1]->Fill(thetaE, thetamu);
                            // h_bvZposFit_nuclear[1]->Fill(ZposFit);
                        }
                    }
                    else if (ProcessID_primary == 4)
                    {

                        h2_ThetaMuVsThetaE_decay[0]->Fill(thetaE, thetamu);
                        // h_bvZposFit_[0]->Fill(ZposFit);
                        cout << "*Muon decay event: " << i << endl;
                        cout << "***ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                        cout << "***ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_decay[1]->Fill(thetaE, thetamu);
                            // h_bvZposFit_[1]->Fill(ZposFit);
                        }
                    }

                    else
                    {
                        h2_ThetaMuVsThetaE_other[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_other[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_other[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_other[1]->Fill(ZposFit);
                        }
                    }
                }
                else
                {

                    if (ProcessID_primary == 25)
                    {

                        h2_ThetaMuVsThetaE_nuclear[0]->Fill(thetaE, thetamu);
                        // h_bvZposFit_nuclear[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_nuclear[1]->Fill(thetaE, thetamu);
                            // h_bvZposFit_nuclear[1]->Fill(ZposFit);
                        }
                    }
                    else if (ProcessID_primary == 4)
                    {

                        h2_ThetaMuVsThetaE_decay[0]->Fill(thetaE, thetamu);
                        // h_bvZposFit_[0]->Fill(ZposFit);
                        cout << "*Muon decay event: " << i << endl;
                        cout << "***ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                        cout << "***ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_decay[1]->Fill(thetaE, thetamu);
                            // h_bvZposFit_[1]->Fill(ZposFit);
                        }
                    }
                    else if (ProcessID_primary == 8 && ProcessID_hardestChargedParticle == 5 && PDG_primary == 22)
                    {
                        h2_ThetaMuVsThetaE_ppPhot[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_ppPhot[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_ppPhot[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_ppPhot[1]->Fill(ZposFit);
                        }
                    }

                    else
                    {
                        // cout << "-Bad reconstruction: invalid linkIDoutE (" << linkIDoutE << ") at event " << i << endl;
                        // cout << "---ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                        // cout << "---ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                        h2_ThetaMuVsThetaE_badReco[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_badReco[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_badReco[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_badReco[1]->Fill(ZposFit);
                        }
                    }
                }

            } // end if reco bad

            else // if reconstructed good
            {
                // if (ProcessID_hardestChargedParticle != processIDE)
                // {
                //     notSameHardestReco++;

                //     if (ProcessID_hardestChargedParticle == 9 && processIDE == 0)
                //     {
                //         h2_ThetaMuVsThetaE_signal[0]->Fill(thetaE, thetamu);
                //         h_bvZposFit_signal[0]->Fill(ZposFit);

                //         if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 912.7) <= 3.0)
                //         {
                //             h2_ThetaMuVsThetaE_signal[1]->Fill(thetaE, thetamu);
                //             h_bvZposFit_signal[1]->Fill(ZposFit);
                //         }
                //     }
                //     else if (ProcessID_hardestChargedParticle == 5 && processIDE == 0)
                //     {
                //         h2_ThetaMuVsThetaE_pp[0]->Fill(thetaE, thetamu);
                //         h_bvZposFit_pp[0]->Fill(ZposFit);

                //         if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 912.7) <= 3.0)
                //         {
                //             h2_ThetaMuVsThetaE_pp[1]->Fill(thetaE, thetamu);
                //             h_bvZposFit_pp[1]->Fill(ZposFit);
                //         }
                //     }
                //     else if (ProcessID_hardestChargedParticle == 5 && processIDE == 9) // && ProcessID_primary == 5 && (PDG_hardestChargedParticle == -11 || PDG_hardestChargedParticle == 11)
                //     {                                                                  // it is a pair production but reconstructed as signal
                //         h2_ThetaMuVsThetaE_pp[0]->Fill(thetaE, thetamu);
                //         h_bvZposFit_pp[0]->Fill(ZposFit);

                //         if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 912.7) <= 3.0)
                //         {
                //             h2_ThetaMuVsThetaE_pp[1]->Fill(thetaE, thetamu);
                //             h_bvZposFit_pp[1]->Fill(ZposFit);
                //         }
                //     }

                //     else
                //     {
                //         cout << "#ProcessID_hardestChargedParticle and ProcessID of reconstructed outgoing electron not the same at event " << i << endl;
                //         cout << "###ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " ProcessID of reconstructed outgoing electron: " << processIDE << endl;
                //         cout << "###ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                //         cout << "###ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                //         h2_ThetaMuVsThetaE_notsameHardestReco[0]->Fill(thetaE, thetamu);
                //         h_bvZposFit_notsameHardestReco[0]->Fill(ZposFit);

                //         if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 912.7) <= 3.0)
                //         {
                //             h2_ThetaMuVsThetaE_notsameHardestReco[1]->Fill(thetaE, thetamu);
                //             h_bvZposFit_notsameHardestReco[1]->Fill(ZposFit);
                //             cout << "###ooooooo####" << endl;
                //             cout << "#Remaining event with different ProcessID_hardestChargedParticle and ProcessID of reconstructed outgoing electron" << endl;
                //             cout << "###ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " ProcessID of reconstructed outgoing electron: " << processIDE << endl;
                //             cout << "###ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                //             cout << "###ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;
                //         }
                //     }
                // }

                // else // if ProcessID_hardestChargedParticle == processIDE
                // {
                if (ProcessID_primary != ProcessID_hardestChargedParticle)
                {

                    if (ProcessID_primary == 8 && ProcessID_hardestChargedParticle == 5 && PDG_primary == 22)
                    {
                        h2_ThetaMuVsThetaE_ppPhot[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_ppPhot[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_ppPhot[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_ppPhot[1]->Fill(ZposFit);
                        }
                    }
                    else if (ProcessID_primary == 25)
                    {

                        h2_ThetaMuVsThetaE_nuclear[0]->Fill(thetaE, thetamu);
                        // h_bvZposFit_nuclear[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_nuclear[1]->Fill(thetaE, thetamu);
                            // h_bvZposFit_nuclear[1]->Fill(ZposFit);
                        }
                    }
                    else if ((ProcessID_primary == 4 || ProcessID_primary == 10) && PDG_primary == 22)
                    {
                        h2_ThetaMuVsThetaE_decay[0]->Fill(thetaE, thetamu);
                        cout << "*Muon decay event: " << i << endl;
                        cout << "***ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                        cout << "***ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_decay[1]->Fill(thetaE, thetamu);
                        }
                    }
                    else if (ProcessID_primary == 8 && PDG_primary == 22 && (ProcessID_hardestChargedParticle == 9 || ProcessID_hardestChargedParticle == 6))
                    {
                        h2_ThetaMuVsThetaE_bremst[0]->Fill(thetaE, thetamu);
                        cout << "--BremStralung event: " << i << endl;
                        cout << "-ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                        cout << "-ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_bremst[1]->Fill(thetaE, thetamu);
                        }
                    }
                    else
                    {
                        notSameHardestPrimary++;
                        // cout << "********************************************************************************************************************************************" << endl;
                        // cout << "*Primary particle and hardest charged particle not the same for event: " << i << endl;
                        // cout << "***ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                        // cout << "***ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                        h2_ThetaMuVsThetaE_notSamePrimaryHardest[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_notSamePrimaryHardest[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_notSamePrimaryHardest[1]->Fill(ZposFit);
                        }
                    }
                }
                else // if hardest and primary are the same
                {
                    if (ProcessID_primary == 5)
                    {
                        h2_ThetaMuVsThetaE_pp[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_pp[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_pp[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_pp[1]->Fill(ZposFit);
                        }
                    }

                    else if (ProcessID_primary == 9)
                    {
                        h2_ThetaMuVsThetaE_signal[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_signal[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_signal[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_signal[1]->Fill(ZposFit);
                        }
                    }
                    else if (ProcessID_primary == 25)
                    {

                        h2_ThetaMuVsThetaE_nuclear[0]->Fill(thetaE, thetamu);
                        // h_bvZposFit_nuclear[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_nuclear[1]->Fill(thetaE, thetamu);
                            // h_bvZposFit_nuclear[1]->Fill(ZposFit);
                        }
                    }
                    else if (ProcessID_primary == 4)
                    {

                        h2_ThetaMuVsThetaE_decay[0]->Fill(thetaE, thetamu);
                        cout << "*Muon decay event: " << i << endl;
                        cout << "***ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                        cout << "***ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                        // h_bvZposFit_[0]->Fill(ZposFit);

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_decay[1]->Fill(thetaE, thetamu);
                            // h_bvZposFit_[1]->Fill(ZposFit);
                        }
                    }

                    else
                    {
                        h2_ThetaMuVsThetaE_other[0]->Fill(thetaE, thetamu);
                        h_bvZposFit_other[0]->Fill(ZposFit);
                        // cout << "Other processes before the cuts event:" << i << endl;
                        // cout << "ProcessID of reconstructed outgoing electron: " << processIDE << endl;
                        // cout << "ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                        // cout << "ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                        if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                        {
                            h2_ThetaMuVsThetaE_other[1]->Fill(thetaE, thetamu);
                            h_bvZposFit_other[1]->Fill(ZposFit);
                        }
                    }
                }
                //}
            }
        }
        else
        {
            // cout << "Error: One or both of primaryTrack and chargedTrack are null at event " << i << "!" << endl;
            notTrack++;
        }

        // int pdgInMu = trackInMu->pdgCode();
        // int intIDInMu = trackInMu->interactionID();

        // int pdgOutE = trackOutE->pdgCode();
        // int intIDoutE = trackOutE->interactionID();
        // int motherIDoutE = trackOutE->motherID();
        // int startZoutE = trackOutE->startZ();

        // int pdgOutMu = trackOutMu->pdgCode();
        // int intIDoutMu = trackOutMu->interactionID();
        // int motherIDoutMu = trackOutMu->motherID();
        // int startZoutMu = trackOutMu->startZ();

        // double thetaE = best_vrtx.electronTheta();
        // double thetamu = best_vrtx.muonTheta();

        // // If particle lables are not correctly assigned (based on angles) we need to swap the containers to match particle type
        // // if ((pdgOutMu == 11 || pdgOutMu == -11) && pdgOutE == -13 && intIDoutE == 0)
        // // {
        // //     swapped += 1;
        // //     // Swap the particles
        // //     swap(oute, outmu);
        // //     swap(linkIDoutE, linkIDoutMu);
        // //     swap(pdgOutE, pdgOutMu);
        // //     swap(intIDoutE, intIDoutMu);
        // //     swap(motherIDoutE, motherIDoutMu);
        // //     swap(startZoutE, startZoutMu);
        // //     swap(thetaE, thetamu);
        // // }

        // // BOTH MU
        // // if (pdgOutMu == -13 && intIDoutMu == 0 && pdgOutE == -13 && intIDoutE == 0) // two outgoing muons
        // //{
        // //  here add primarycharged change from G
        // // bothMu++;
        // //  WeirdTree->Fill();
        // // h_bvZposFit2Mu->Fill(ZposFit);
        // // continue;
        // //}

        // // if (intIDoutE == 0) // only mu+
        // // {
        // //     intE0++;
        // //     pdgE0.push_back(pdgOutE);
        // // }

        // h_mctPdgCodeOutE->Fill(pdgOutE);
        // h_mctInteractionIDoutE->Fill(intIDoutE);
        // h_mctStartZ->Fill(startZoutE);

    } // end for loop for all events

    cout << "Number of all the skimmed events: " << nEntries << endl;
    cout << "Number of reconstructed events: " << Nreco << endl;
    cout << "Percentage of reconstructed events: " << double(Nreco) / double(nEntries) << endl;
    cout << "Number of events after the cuts: " << afterCuts << endl;
    // cout << "Number of tracks for one event: " << nMCTracks[0] << endl;
    cout << "Number of swapped events: " << swapped << endl;
    // cout << "Number of events which have reconstructed two muon tracks: " << bothMu << ", Percentage: " << 100 * double(bothMu) / double(Nreco) << endl;
    // cout << "Number of out e with interactionID = 0: " << intE0 << endl;
    // cout << "Number of events that don't have incoming muon: " << notInMu << endl;
    // // cout << pdgE0[0] << " " << pdgE0[1] << " " << pdgE0[2] << " " << pdgE0[3] << " " << pdgE0[4] << " " << pdgE0[5] << " " << pdgE0[6] << endl;

    // cout << "Percentage of pair production events: " << 100 * double(pair) / double(Nreco) << " %" << endl;
    // cout << "Percentage of signal events: " << 100 * double(signal) / double(Nreco) << " %" << endl;
    // cout << "Percentage of other events: " << 100 * double(other) / double(Nreco) << " %," << " Number: " << other << endl;

    // cout << "Percentage of pair producted events from muons (of pair production events): " << mu1 * 100.0 / totalParents << " %" << endl;
    // cout << "Percentage of pair producted events from photons (of pair production events): " << phot1 * 100.0 / totalParents << " %" << endl;
    // cout << "Percentage of pair producted events from other particles (of pair production events): " << other1 * 100.0 / totalParents << " %" << endl;

    // cout << "Percentage of pair producted events from muons (of total number of reco events): " << mu1 * 100.0 / double(Nreco) << " %" << endl;
    // cout << "Percentage of pair producted events from photons (of total number of reco events): " << phot1 * 100.0 / double(Nreco) << " %" << endl;
    // cout << "Percentage of pair producted events from other particles (of total number of reco events): " << other1 * 100.0 / double(Nreco) << " %" << endl;

    // cout << "Theta min cut" << endl;
    // cout << "Signal: " << thetaMuminSignal << endl;
    // cout << "Pair prod. from muons: " << thetaMuminPPmu << endl;
    // cout << "Pair prod. from photons: " << thetaMuminPPphot << endl;
    // cout << "Other: " << thetaMuminOther << endl;

    // cout << "Theta max cut" << endl;
    // cout << "Signal: " << thetaEmaxSignal << endl;
    // cout << "Pair prod. from muons: " << thetaEmaxPPmu << endl;
    // cout << "Pair prod. from photons: " << thetaEmaxPPphot << endl;
    // cout << "Other: " << thetaEmaxOther << endl;

    // cout << "Number of events that don't have same reco and mctruth hardest particle: " << notSameHardestReco* 100.0 / double(Nreco) << endl;
    // cout << "Number of events that don't have same hardest and primary particle: " << notSameHardestPrimary * 100.0 / double(Nreco)<< endl;

    //////////Elastic curve/////////////////////////////////////////
    /// using ElasticState.cc, ElasticState.h and MuE.cc (examle of usage)
    /// to make a library: g++ -O2 -fPIC -shared -o libElasticState.so ElasticState.cc $(root-config --cflags --libs)
    gSystem->Load("libElasticState.so"); // include library
    double energy_mu_in = 160.0;
    double mass_mu = 0.105658;
    double mass_e = 0.000511;
    double theta_e_out = 1.0;

    MuE::ElasticState elastic(energy_mu_in, mass_mu, mass_e, theta_e_out);
    // TCanvas *canvas = elastic.draw_elasticCurve();
    // canvas->Update();
    elastic.calculate_elasticCurve();
    // Retrieve the TGraph pointer
    TGraph *gr = elastic.GetGraph();

    // Convert graph points from mrad to rad (if needed)
    int nPoints = gr->GetN();
    TGraph *grMirrored = new TGraph(nPoints);
    TGraph *grShiftedPlus = new TGraph(nPoints);
    TGraph *grShiftedMinus = new TGraph(nPoints);
    TGraph *grMirroredPlus = new TGraph(nPoints);
    TGraph *grMirroredMinus = new TGraph(nPoints);

    Double_t *xvals = gr->GetX();
    Double_t *yvals = gr->GetY();
    double dx = 0.0002 * cos(M_PI / 4); // 0.2 * cos(45°)
    double dy = 0.0002 * sin(M_PI / 4); // 0.2 * sin(45°)

    for (int i = 0; i < nPoints; i++)
    {
        xvals[i] /= 1000.0; // converting from mrad to rad
        yvals[i] /= 1000.0;

        grMirrored->SetPoint(i, yvals[i], xvals[i]);
        grMirroredPlus->SetPoint(i, yvals[i] + 0.0002, xvals[i]);
        grMirroredMinus->SetPoint(i, yvals[i] - 0.0002, xvals[i]);

        // grShiftedPlus->SetPoint(i, xvals[i] + dx, yvals[i] + dy);
        // grShiftedMinus->SetPoint(i, xvals[i] - dx, yvals[i] - dy);
        grShiftedPlus->SetPoint(i, xvals[i], yvals[i] + 0.0002);
        grShiftedMinus->SetPoint(i, xvals[i], yvals[i] - 0.0002);
    }

    // Mirrored curve

    //////////////////////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////////////////////
    // /// Plot angles for signal and background
    TCanvas *cAngBS = new TCanvas("cAngBS", "cAngBS", 800, 800);
    //cAngBS->Divide(2);
    //int j = 0;
    //for (int i = 0; i < nCuts; ++i)
    //{
        //cAngBS->cd(i + 1);
        cAngBS->cd();
        h2_ThetaMuVsThetaE_signal[1]->SetStats(0);
        h2_ThetaMuVsThetaE_signal[1]->SetMarkerColor(4); // blue
        h2_ThetaMuVsThetaE_signal[1]->SetLineColor(4);
        h2_ThetaMuVsThetaE_signal[1]->SetMarkerSize(0.5);
        h2_ThetaMuVsThetaE_signal[1]->SetMarkerStyle(4);
        h2_ThetaMuVsThetaE_signal[1]->SetTitle(""); // Cut %d", i
        //j++;
        h2_ThetaMuVsThetaE_signal[1]->GetXaxis()->SetTitle("Electron angle [rad]");
        h2_ThetaMuVsThetaE_signal[1]->GetXaxis()->SetTitleSize(0.04);
        h2_ThetaMuVsThetaE_signal[1]->GetXaxis()->SetTitleOffset(1.1);
        h2_ThetaMuVsThetaE_signal[1]->GetYaxis()->SetTitle("Muon angle [rad]");
        h2_ThetaMuVsThetaE_signal[1]->GetYaxis()->SetTitleSize(0.04);
        h2_ThetaMuVsThetaE_signal[1]->Draw("SCAT");

        h2_ThetaMuVsThetaE_pp[1]->SetStats(0);
        h2_ThetaMuVsThetaE_pp[1]->SetMarkerColor(8); // 1 black 5 yellow
        h2_ThetaMuVsThetaE_pp[1]->SetLineColor(5);
        h2_ThetaMuVsThetaE_pp[1]->SetMarkerSize(1);
        h2_ThetaMuVsThetaE_pp[1]->SetMarkerStyle(22);
        //h2_ThetaMuVsThetaE_pp[1]->Draw("SAME SCAT");

        h2_ThetaMuVsThetaE_other[1]->SetStats(0);
        h2_ThetaMuVsThetaE_other[1]->SetMarkerColor(2); // red
        h2_ThetaMuVsThetaE_other[1]->SetLineColor(2);
        h2_ThetaMuVsThetaE_other[1]->SetMarkerSize(1);
        h2_ThetaMuVsThetaE_other[1]->SetMarkerStyle(25);
        h2_ThetaMuVsThetaE_other[1]->Draw("SAME SCAT");

        h2_ThetaMuVsThetaE_ppPhot[1]->SetStats(0);
        h2_ThetaMuVsThetaE_ppPhot[1]->SetMarkerStyle(33);
        h2_ThetaMuVsThetaE_ppPhot[1]->SetMarkerColor(6); // violet
        h2_ThetaMuVsThetaE_ppPhot[1]->SetLineColor(6);
        h2_ThetaMuVsThetaE_ppPhot[1]->SetMarkerSize(1.2);
        h2_ThetaMuVsThetaE_ppPhot[1]->Draw("SAME SCAT");

        h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetStats(0);
        h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetMarkerStyle(28);
        h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetMarkerColor(8); // green mat
        h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetLineColor(8);
        h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetMarkerSize(1);
        h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->Draw("SAME SCAT");

        h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetStats(0);
        h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetMarkerStyle(28);
        h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetMarkerColor(5); // yellow
        h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetLineColor(5);
        h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetMarkerSize(1);
        h2_ThetaMuVsThetaE_notsameHardestReco[1]->Draw("SAME SCAT");

        h2_ThetaMuVsThetaE_badReco[1]->SetStats(0);
        h2_ThetaMuVsThetaE_badReco[1]->SetMarkerStyle(28);
        h2_ThetaMuVsThetaE_badReco[1]->SetMarkerColor(46); // redish brown
        h2_ThetaMuVsThetaE_badReco[1]->SetLineColor(46);
        h2_ThetaMuVsThetaE_badReco[1]->SetMarkerSize(1);
        h2_ThetaMuVsThetaE_badReco[1]->Draw("SAME SCAT");

        h2_ThetaMuVsThetaE_nuclear[1]->SetStats(0);
        h2_ThetaMuVsThetaE_nuclear[1]->SetMarkerStyle(34);
        h2_ThetaMuVsThetaE_nuclear[1]->SetMarkerColor(7); //
        h2_ThetaMuVsThetaE_nuclear[1]->SetLineColor(7);
        h2_ThetaMuVsThetaE_nuclear[1]->SetMarkerSize(1.2);
        h2_ThetaMuVsThetaE_nuclear[1]->Draw("SAME SCAT");

        h2_ThetaMuVsThetaE_decay[1]->SetStats(0);
        h2_ThetaMuVsThetaE_decay[1]->SetMarkerStyle(28);
        h2_ThetaMuVsThetaE_decay[1]->SetMarkerColor(51); //
        h2_ThetaMuVsThetaE_decay[1]->SetLineColor(51);
        h2_ThetaMuVsThetaE_decay[1]->SetMarkerSize(1);
        h2_ThetaMuVsThetaE_decay[1]->Draw("SAME SCAT");

        h2_ThetaMuVsThetaE_bremst[1]->SetStats(0);
        h2_ThetaMuVsThetaE_bremst[1]->SetMarkerStyle(28);
        h2_ThetaMuVsThetaE_bremst[1]->SetMarkerColor(28); //
        h2_ThetaMuVsThetaE_bremst[1]->SetLineColor(28);
        h2_ThetaMuVsThetaE_bremst[1]->SetMarkerSize(1);
        h2_ThetaMuVsThetaE_bremst[1]->Draw("SAME SCAT");

        gr->SetLineColor(205);
        gr->SetLineWidth(2);
        gr->Draw("L same");

        // grMirrored->SetLineColor(17);
        // grMirrored->SetLineWidth(2);
        // grMirrored->Draw("L same");

        grShiftedPlus->SetLineColor(205);
        grShiftedPlus->SetLineWidth(2);
        grShiftedPlus->SetLineStyle(9);
        grShiftedPlus->Draw("L same");

        grShiftedMinus->SetLineColor(205);
        grShiftedMinus->SetLineWidth(2);
        grShiftedMinus->SetLineStyle(9);
        grShiftedMinus->Draw("L same");


        h2_ThetaMuVsThetaE_pp[1]->Draw("SAME SCAT");

        // grMirroredPlus->SetLineColor(17);
        // grMirroredPlus->SetLineStyle(9);
        // grMirroredPlus->SetLineWidth(2);
        // grMirroredPlus->Draw("L same");

        // grMirroredMinus->SetLineColor(17);
        // grMirroredMinus->SetLineWidth(2);
        // grMirroredMinus->SetLineStyle(9);
        // grMirroredMinus->Draw("L same");

        cAngBS->Update();
        cAngBS->Modified();

        // Legend
        int entries9 = h2_ThetaMuVsThetaE_signal[1]->GetEntries();
        int entries5 = h2_ThetaMuVsThetaE_pp[1]->GetEntries();
        int entries5phot = h2_ThetaMuVsThetaE_ppPhot[1]->GetEntries();
        int entriesOther = h2_ThetaMuVsThetaE_other[1]->GetEntries();
        int enotsamePH = h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->GetEntries();
        int enotsameHR = h2_ThetaMuVsThetaE_notsameHardestReco[1]->GetEntries();
        int ebadReco = h2_ThetaMuVsThetaE_badReco[1]->GetEntries();
        int enuclear = h2_ThetaMuVsThetaE_nuclear[1]->GetEntries();
        int edecay = h2_ThetaMuVsThetaE_decay[1]->GetEntries();
        int ebremst = h2_ThetaMuVsThetaE_bremst[1]->GetEntries();

        TLegend *legend = new TLegend(0.15, 0.73, 0.9, 0.9);
        legend->SetTextSize(0.04);
        legend->AddEntry(h2_ThetaMuVsThetaE_signal[1], Form("Signal: %d events", entries9), "p");
        legend->AddEntry(h2_ThetaMuVsThetaE_pp[1], Form("Pair prod. from muons: %d events", entries5), "p");
        legend->AddEntry(h2_ThetaMuVsThetaE_ppPhot[1], Form("Pair prod. from photons: %d events", entries5phot), "p");
        //legend->AddEntry(h2_ThetaMuVsThetaE_signal[1], "Signal", "p");
        //legend->AddEntry(h2_ThetaMuVsThetaE_pp[1], "Pair prod. from muons", "p");
        //legend->AddEntry(h2_ThetaMuVsThetaE_ppPhot[1], "Pair prod. from photons", "p");
        // legend->AddEntry(h2_ThetaMuVsThetaE_NoMu[i], Form("Events from not muons: %d events", entriesNoMu), "p");
        // legend->AddEntry(h2_ThetaMuVsThetaE_other[i], Form("Other background: %d events", entriesOther), "p");
        legend->AddEntry(h2_ThetaMuVsThetaE_nuclear[1], Form("Nuclear interaction: %d events", enuclear), "p");
        // legend->AddEntry(h2_ThetaMuVsThetaE_decay[i], Form("Muon decay: %d events", edecay), "p");
        // legend->AddEntry(h2_ThetaMuVsThetaE_bremst[i], Form("Bremsstrahlung: %d events", ebremst), "p");

        // legend->AddEntry(h2_ThetaMuVsThetaE_notSamePrimaryHardest[i], Form("Not same hardest and primary: %d events", enotsamePH), "p");
        // legend->AddEntry(h2_ThetaMuVsThetaE_notsameHardestReco[i], Form("Not same hardest and reco: %d events", enotsameHR), "p");
        // legend->AddEntry(h2_ThetaMuVsThetaE_badReco[i], Form("Bad reconstruction: %d events", ebadReco), "p");

        // cAngBS->BuildLegend();
        legend->Draw();
        cAngBS->Update();
        cAngBS->Modified();
        cAngBS->SaveAs("MCAnalisys_total_2AddCond_NewLink_allClass_swap_7.4.25.pdf");
    //}


    TCanvas *cfitZBS = new TCanvas("cfitZBS", "cfitZBS", 800, 600);

    cfitZBS->cd();
    h_bvZposFit_signal[0]->Draw("");
    cfitZBS->Update();
    cfitZBS->Modified();

    // // fEventDisplay->Close();

    TFile *outFile = new TFile("MCAnalisys_total_2AddCond_NewLink_allClass_swap_2.4.25.root", "RECREATE");
    cAngBS->Write();
    cfitZBS->Write();
    // outFile->Close();
    // // fEventDisplay->Close();

    // Restore the original cout buffer
    cout.rdbuf(coutBuf);
    logFile.close();

    // f->Close();
    // return 0;
}