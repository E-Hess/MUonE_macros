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

// CHANGED LINKING
// " new per-event variables which give you: https://indico.cern.ch/event/1505154/contributions/6335583/attachments/3001752/5290192/2025.01.23_GAbbiendi_FM-MCtruth-upd.pdf
// the interactionID of the produced particle (charged or neutral) with maximum energy (apart from the beam muon)
// the interactionID of the produced charged particle with maximum energy (apart from the beam muon)
// These two can be identical or not, depending on what happened. For example:
// -+ if there was an elastic scattering and the scattered electron is the most energetic particle both will be defined by it
// -+ if there was a radiated photon that at some point in the detector converts in a e+e- pair and those are the most energetic particles
// you will have the first variable saying ‘bremsstrahlung’ and the second ‘pair production’
// If nothing happens and there is just the beam muon the two variables are left at the value = 0"

// plots if reco with swap

void runMCAnalysis_NewLinking_allClass_swap()
{
    // Create and open a log file
    ofstream logFile("MCAnalysis_skimv2_NewLink_swap_2.5.25.log");

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

    int intE0 = 0;
    int bothMu = 0;
    int notInMu = 0;
    int swapped = 0;
    int ppp = 0;
    int notSameHardestReco = 0;
    int notSameHardestPrimary = 0;
    int notTrack = 0, oneTrack = 0;
    int badRecoEvents = 0;
    int afterCuts = 0;
    vector<string> cuts = {"No cuts", "All cuts"};

    cout << "Number of events in file: " << nEntries << endl;

    for (Long64_t i = 0; i < nEntries; ++i) // nEntries
    {
        cbmsim->GetEntry(i);
        if (i % 50000 == 0)
            cout << "Event: " << i << endl;
        if (ReconstructionOutput->isReconstructed() == 0) // Taking only reconstructed events (they have a vertex)
            continue;

        Nreco++;

        auto best_vrtx = ReconstructionOutput->bestVertex();
        auto oute = best_vrtx.outgoingElectron();
        auto outmu = best_vrtx.outgoingMuon();
        auto inmu = best_vrtx.incomingMuon();
        auto hits = ReconstructionOutput->reconstructedHits();

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

        int recoHitsS1 = 0;
        for (int h = 0; h < hits.size(); h++)
        {
            if (hits.at(h).stationID() == 1)
                recoHitsS1++;
        }

        bool passSelectionCuts =
            (recoHitsS1 <= 14) &&
            (std::fabs(modAcop) < 0.4) &&
            (thetamu >= 0.0002) &&
            (thetaE < 0.032) &&
            (chi2vert <= 20) &&
            (std::fabs(ZposFit - 911.2) <= 3.0); // 900?

        if (passSelectionCuts)
            afterCuts++;

        //////////////// Hardest charged and primary particles from the interaction///////////////////////

        Double_t Emax = 0.;
        Double_t pmax = 0.;
        Double_t Emax_ch = 0.;
        Double_t Zmin = 10000.;
        const MUonETrack *primaryTrack = nullptr;
        const MUonETrack *chargedTrack = nullptr;
        const MUonETrack *ZminTrack = nullptr;
        //const MUonETrack *primaryMotherTrack = nullptr;
        //const MUonETrack *chargedMotherTrack = nullptr;
        int primaryMother = -1;
        int chargedMother = -1;
        int ZminMother = -1;

        Int_t ProcessID_primary = -1;
        Int_t ProcessID_hardestChargedParticle = -1;
        Int_t ProcessID_Zmin = -1;
        Int_t PDG_primary = 0;
        Int_t PDG_hardestChargedParticle = 0;
        Int_t PDG_Zmin = 0;

        int nTracks = MCTrack->GetEntries();
        if (nTracks < 2) // If reco event has only one track (it is wrong)
        {
            cout << "----------------------------------------------------------" << endl;
            cout << "One track for reco event: " << i << endl;
            oneTrack++;
            continue;
        }

        for (int n = 0; n < nTracks; n++) // Looping over all MC tracks
        {
            const auto MCTr = static_cast<const MUonETrack *>(MCTrack->At(n));

            int pdg = MCTr->pdgCode();
            int intID = MCTr->interactionID();
            int motherID = MCTr->motherID();
            double startZ = MCTr->startZ();

            //const MUonETrack *motherTrack = static_cast<const MUonETrack *>(MCTrack->At(motherID));
            //int motherPDG = motherTrack->pdgCode();

            if (n < 1)
                continue; // Skip incoming muon

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
                //primaryMother = motherPDG;
            }

            if (charge != 0 && p > pmax) // Update hardest charged track (highest momentum among charged tracks)
            {
                pmax = p;
                chargedTrack = MCTr;
                //chargedMother = motherPDG;
            }

            if (startZ < Zmin)
            {
                Zmin = startZ;
                ZminTrack = MCTr;
                //ZminMother = motherPDG;
            }
        }

        if (primaryTrack && chargedTrack) // if primaryTrack && chargedTrack are found
        {
            // Use the found tracks' interaction IDs
            ProcessID_primary = primaryTrack->interactionID();
            ProcessID_hardestChargedParticle = chargedTrack->interactionID();
            ProcessID_Zmin = ZminTrack->interactionID();
            PDG_primary = primaryTrack->pdgCode();
            PDG_hardestChargedParticle = chargedTrack->pdgCode();
            PDG_Zmin = ZminTrack->pdgCode();

            if (processIDE == 0 && processIDoutMu > 0) //
            {
                swapped += 1;
                // Swap the particles
                swap(oute, outmu);
                // swap(linkIDoutE, linkIDoutMu);
                swap(thetaE, thetamu);
            }

            if (ProcessID_primary != ProcessID_hardestChargedParticle)
            {
                notSameHardestPrimary++;

                h2_ThetaMuVsThetaE_notSamePrimaryHardest[0]->Fill(thetaE, thetamu);
                // h_bvZposFit_notSamePrimaryHardest[0]->Fill(ZposFit);

                if (passSelectionCuts)
                {
                    h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->Fill(thetaE, thetamu);
                    // h_bvZposFit_notSamePrimaryHardest[1]->Fill(ZposFit);
                }

                if (ProcessID_primary == 8) // Bremsstrahlung
                {
                    if (ProcessID_hardestChargedParticle == 5 && PDG_primary == 22)
                    {
                        h2_ThetaMuVsThetaE_ppPhot[0]->Fill(thetaE, thetamu);
                        // h_bvZposFit_ppPhot[0]->Fill(ZposFit);

                        if (passSelectionCuts)
                        {
                            h2_ThetaMuVsThetaE_ppPhot[1]->Fill(thetaE, thetamu);
                            // h_bvZposFit_ppPhot[1]->Fill(ZposFit);
                        }
                    }

                    else // Other cases are bremsstrahlung + any other secondary process
                    {
                        h2_ThetaMuVsThetaE_bremst[0]->Fill(thetaE, thetamu);

                        if (passSelectionCuts)
                        {
                            h2_ThetaMuVsThetaE_bremst[1]->Fill(thetaE, thetamu);
                        }
                    }
                } // end Bremsstrahlung

                else if (ProcessID_primary == 25)
                {
                    h2_ThetaMuVsThetaE_nuclear[0]->Fill(thetaE, thetamu);
                    // h_bvZposFit_nuclear[0]->Fill(ZposFit);

                    if (passSelectionCuts)
                    {
                        h2_ThetaMuVsThetaE_nuclear[1]->Fill(thetaE, thetamu);
                        // h_bvZposFit_nuclear[1]->Fill(ZposFit);
                    }
                }

                else
                {
                    h2_ThetaMuVsThetaE_other[0]->Fill(thetaE, thetamu);
                    // h_bvZposFit_other[0]->Fill(ZposFit);

                    if (passSelectionCuts)
                    {
                        h2_ThetaMuVsThetaE_other[1]->Fill(thetaE, thetamu);
                        // h_bvZposFit_other[1]->Fill(ZposFit);
                    }
                }
            } // end if hardest and primary are not the same

            else // if hardest and primary are the same
            {
                if (ProcessID_primary == 5)
                {
                    h2_ThetaMuVsThetaE_pp[0]->Fill(thetaE, thetamu);
                    // h_bvZposFit_pp[0]->Fill(ZposFit);

                    if (passSelectionCuts)
                    {
                        h2_ThetaMuVsThetaE_pp[1]->Fill(thetaE, thetamu);
                        // h_bvZposFit_pp[1]->Fill(ZposFit);
                    }
                }

                else if (ProcessID_primary == 9)
                {
                    h2_ThetaMuVsThetaE_signal[0]->Fill(thetaE, thetamu);
                    // h_bvZposFit_signal[0]->Fill(ZposFit);

                    if (passSelectionCuts)
                    {
                        h2_ThetaMuVsThetaE_signal[1]->Fill(thetaE, thetamu);
                        // h_bvZposFit_signal[1]->Fill(ZposFit);
                    }
                }
                else if (ProcessID_primary == 25)
                {

                    h2_ThetaMuVsThetaE_nuclear[0]->Fill(thetaE, thetamu);
                    // h_bvZposFit_nuclear[0]->Fill(ZposFit);

                    if (passSelectionCuts)
                    {
                        h2_ThetaMuVsThetaE_nuclear[1]->Fill(thetaE, thetamu);
                        // h_bvZposFit_nuclear[1]->Fill(ZposFit);
                    }
                }
                else if (ProcessID_primary == 4)
                {

                    h2_ThetaMuVsThetaE_decay[0]->Fill(thetaE, thetamu);
                    // cout << "*Muon decay event: " << i << endl;
                    // cout << "***ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                    // cout << "***ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                    // h_bvZposFit_[0]->Fill(ZposFit);

                    if (passSelectionCuts)
                    {
                        h2_ThetaMuVsThetaE_decay[1]->Fill(thetaE, thetamu);
                        // h_bvZposFit_[1]->Fill(ZposFit);
                    }
                }

                else
                {
                    h2_ThetaMuVsThetaE_other[0]->Fill(thetaE, thetamu);
                    // h_bvZposFit_other[0]->Fill(ZposFit);
                    //  cout << "Other processes before the cuts event:" << i << endl;
                    //  cout << "ProcessID of reconstructed outgoing electron: " << processIDE << endl;
                    //  cout << "ProcessID_primary: " << ProcessID_primary << " PDG_primary: " << PDG_primary << endl;
                    //  cout << "ProcessID_hardestChargedParticle: " << ProcessID_hardestChargedParticle << " PDG_hardestChargedParticle: " << PDG_hardestChargedParticle << endl;

                    if (recoHitsS1 <= 14 && fabs(modAcop) < 0.4 && thetamu >= 0.0002 && thetaE < 0.032 && chi2vert <= 20 && fabs(ZposFit - 911.2) <= 3.0)
                    {
                        h2_ThetaMuVsThetaE_other[1]->Fill(thetaE, thetamu);
                        // h_bvZposFit_other[1]->Fill(ZposFit);
                    }
                }
            }

        } // end if primaryTrack && chargedTrack are found

        else // if primaryTrack && chargedTrack are NOT found
        {
            // cout << "Error: One or both of primaryTrack and chargedTrack are null at event " << i << "!" << endl;
            notTrack++;
        }

    } // end for loop for all events

    int signal[nCuts], pair[nCuts], pairPhot[nCuts], nuclear[nCuts], mudecay[nCuts], bremst[nCuts], other[nCuts];

    for (int i = 0; i < nCuts; ++i)
    {
        signal[i] = h2_ThetaMuVsThetaE_signal[i]->GetEntries();
        pair[i] = h2_ThetaMuVsThetaE_pp[i]->GetEntries();
        pairPhot[i] = h2_ThetaMuVsThetaE_ppPhot[i]->GetEntries();
        nuclear[i] = h2_ThetaMuVsThetaE_nuclear[i]->GetEntries();
        mudecay[i] = h2_ThetaMuVsThetaE_decay[i]->GetEntries();
        bremst[i] = h2_ThetaMuVsThetaE_bremst[i]->GetEntries();
        other[i] = h2_ThetaMuVsThetaE_other[i]->GetEntries();
    }
    // int enotsamePH = h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->GetEntries();
    // int enotsameHR = h2_ThetaMuVsThetaE_notsameHardestReco[1]->GetEntries();
    // int ebadReco = h2_ThetaMuVsThetaE_badReco[1]->GetEntries();

    cout << "Number of all the skimmed events: " << nEntries << endl;
    cout << "Number of reconstructed events: " << Nreco << endl;
    cout << "Percentage of reconstructed events: " << double(Nreco) / double(nEntries) << endl;
    cout << "Number of events after the cuts: " << afterCuts << endl;
    cout << "Number of swapped events: " << swapped << endl;

    cout << "Number of events where primary != hardest " << notSameHardestPrimary << endl;
    cout << "Number of events with null primary or hardest-charged track: " << notTrack << endl;

    cout << "Number of muIoni (SIGNAL) events: " << signal[0] << ", selected: " << signal[1] << endl;
    cout << "Number of muPair          events: " << pair[0] << ", selected: " << pair[1] << endl;
    cout << "Number of muPhotPair      events: " << pairPhot[0] << ", selected: " << pairPhot[1] << endl;
    cout << "Number of muNucl          events: " << nuclear[0] << ", selected: " << nuclear[1] << endl;
    cout << "Number of muDecay         events: " << mudecay[0] << ", selected: " << mudecay[1] << endl;
    cout << "Number of muBremsstrahlung event: " << bremst[0] << ", selected: " << bremst[1] << endl;
    cout << "Number of unclassified events   : " << other[0] << ", selected: " << other[1] << endl;

    // cout << "Number of events which have reconstructed two muon tracks: " << bothMu << ", Percentage: " << 100 * double(bothMu) / double(Nreco) << endl;
    // cout << "Number of out e with interactionID = 0: " << intE0 << endl;
    // cout << "Number of events that don't have incoming muon: " << notInMu << endl;
    // // cout << pdgE0[0] << " " << pdgE0[1] << " " << pdgE0[2] << " " << pdgE0[3] << " " << pdgE0[4] << " " << pdgE0[5] << " " << pdgE0[6] << endl;

    // cout << "Percentage of pair production events: " << 100 * double(pair) / double(Nreco) << " %" << endl;
    // cout << "Percentage of signal events: " << 100 * double(signal) / double(Nreco) << " %" << endl;
    // cout << "Percentage of other events: " << 100 * double(other) / double(Nreco) << " %," << " Number: " << other << endl;

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

    // //////////////////////////////////////////////////////////////////////////////////////////////////
    // /// Plot angles for signal and background
    TCanvas *cAngBS = new TCanvas("cAngBS", "cAngBS", 800, 800);
    // cAngBS->Divide(2);
    // int j = 0;
    // for (int i = 0; i < nCuts; ++i)
    //{
    // cAngBS->cd(i + 1);
    cAngBS->cd();
    h2_ThetaMuVsThetaE_signal[1]->SetStats(0);
    h2_ThetaMuVsThetaE_signal[1]->SetMarkerColor(4); // blue
    h2_ThetaMuVsThetaE_signal[1]->SetLineColor(4);
    h2_ThetaMuVsThetaE_signal[1]->SetMarkerSize(0.5);
    h2_ThetaMuVsThetaE_signal[1]->SetMarkerStyle(4);
    h2_ThetaMuVsThetaE_signal[1]->SetTitle(""); // Cut %d", i
    // j++;
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
    // h2_ThetaMuVsThetaE_pp[1]->Draw("SAME SCAT");

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

    // h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetStats(0);
    // h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetMarkerStyle(28);
    // h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetMarkerColor(8); // green mat
    // h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetLineColor(8);
    // h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->SetMarkerSize(1);
    // h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->Draw("SAME SCAT");

    // h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetStats(0);
    // h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetMarkerStyle(28);
    // h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetMarkerColor(5); // yellow
    // h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetLineColor(5);
    // h2_ThetaMuVsThetaE_notsameHardestReco[1]->SetMarkerSize(1);
    // h2_ThetaMuVsThetaE_notsameHardestReco[1]->Draw("SAME SCAT");

    // h2_ThetaMuVsThetaE_badReco[1]->SetStats(0);
    // h2_ThetaMuVsThetaE_badReco[1]->SetMarkerStyle(28);
    // h2_ThetaMuVsThetaE_badReco[1]->SetMarkerColor(46); // redish brown
    // h2_ThetaMuVsThetaE_badReco[1]->SetLineColor(46);
    // h2_ThetaMuVsThetaE_badReco[1]->SetMarkerSize(1);
    // h2_ThetaMuVsThetaE_badReco[1]->Draw("SAME SCAT");

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
    // int entries9 = h2_ThetaMuVsThetaE_signal[1]->GetEntries();
    // int entries5 = h2_ThetaMuVsThetaE_pp[1]->GetEntries();
    // int entries5phot = h2_ThetaMuVsThetaE_ppPhot[1]->GetEntries();
    // int entriesOther = h2_ThetaMuVsThetaE_other[1]->GetEntries();
    // int enotsamePH = h2_ThetaMuVsThetaE_notSamePrimaryHardest[1]->GetEntries();
    // int enotsameHR = h2_ThetaMuVsThetaE_notsameHardestReco[1]->GetEntries();
    // int ebadReco = h2_ThetaMuVsThetaE_badReco[1]->GetEntries();
    // int enuclear = h2_ThetaMuVsThetaE_nuclear[1]->GetEntries();
    // int edecay = h2_ThetaMuVsThetaE_decay[1]->GetEntries();
    // int ebremst = h2_ThetaMuVsThetaE_bremst[1]->GetEntries();

    TLegend *legend = new TLegend(0.15, 0.73, 0.9, 0.9);
    legend->SetTextSize(0.04);
    legend->AddEntry(h2_ThetaMuVsThetaE_signal[1], Form("Signal: %d events", signal[1]), "p");
    legend->AddEntry(h2_ThetaMuVsThetaE_pp[1], Form("Pair prod. from muons: %d events", pair[1]), "p");
    legend->AddEntry(h2_ThetaMuVsThetaE_ppPhot[1], Form("Pair prod. from photons: %d events", pairPhot[1]), "p");
    // legend->AddEntry(h2_ThetaMuVsThetaE_signal[1], "Signal", "p");
    // legend->AddEntry(h2_ThetaMuVsThetaE_pp[1], "Pair prod. from muons", "p");
    // legend->AddEntry(h2_ThetaMuVsThetaE_ppPhot[1], "Pair prod. from photons", "p");
    //  legend->AddEntry(h2_ThetaMuVsThetaE_NoMu[i], Form("Events from not muons: %d events", entriesNoMu), "p");
    legend->AddEntry(h2_ThetaMuVsThetaE_other[1], Form("Other background: %d events", other[1]), "p");
    legend->AddEntry(h2_ThetaMuVsThetaE_nuclear[1], Form("Nuclear interaction: %d events", nuclear[1]), "p");
    legend->AddEntry(h2_ThetaMuVsThetaE_decay[1], Form("Muon decay: %d events", mudecay[1]), "p");
    legend->AddEntry(h2_ThetaMuVsThetaE_bremst[1], Form("Bremsstrahlung: %d events", bremst[1]), "p");

    // legend->AddEntry(h2_ThetaMuVsThetaE_notSamePrimaryHardest[i], Form("Not same hardest and primary: %d events", enotsamePH), "p");
    // legend->AddEntry(h2_ThetaMuVsThetaE_notsameHardestReco[i], Form("Not same hardest and reco: %d events", enotsameHR), "p");
    // legend->AddEntry(h2_ThetaMuVsThetaE_badReco[i], Form("Bad reconstruction: %d events", ebadReco), "p");

    // cAngBS->BuildLegend();
    legend->Draw();
    cAngBS->Update();
    cAngBS->Modified();
    cAngBS->SaveAs("MCAnalysis_skimv2_NewLink_swap_2.5.25.pdf");
    //}

    TCanvas *cfitZBS = new TCanvas("cfitZBS", "cfitZBS", 800, 600);

    cfitZBS->cd();
    h_bvZposFit_signal[0]->Draw("");
    cfitZBS->Update();
    cfitZBS->Modified();

    // // fEventDisplay->Close();

    TFile *outFile = new TFile("MCAnalysis_skimv2_NewLink_swap_2.5.25.root", "RECREATE");
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