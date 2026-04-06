/*
 * delta_t_minimization.C
 *
 * Computes the ΔT discriminant for tt̄ → bb̄lτ_h event reconstruction
 * (the "tt1" topology — standard dileptonic tt̄ decay with 2 neutrinos).
 *
 * This implements the ψ_T component of the combined Ψ discriminator, first
 * introduced by Cookman & Harris (2019). For each event, the momenta of two
 * neutrinos are varied to minimise a penalty function ΔT that measures
 * kinematic consistency with the dileptonic tt̄ hypothesis. Events that
 * reconstruct well as tt̄ are assigned a high ψ_T and are therefore less
 * likely to be di-Higgs signal.
 *
 * 71% of tt̄ MC events achieve ΔT < 400 — hence the original filename.
 * The complementary file delta_t3_8combinations.C handles the more complex
 * tt3 topology (6 neutrinos), achieving 94% reconstruction efficiency.
 *
 * For each selected event (≥1e, ≥1μ, ≥2 b-jets, pT > 25 GeV, OS leptons),
 * minimization is run 100 times across 2 b-jet assignment combinations.
 * The minimum ΔT over all restarts and combinations is recorded.
 *
 * Additional histograms are filled for the reconstructed W and top masses,
 * useful for validating the quality of the reconstruction.
 *
 * Part of the ATLAS Di-Higgs Discrimination Project.
 * University of Manchester, MPhys Project, Sep 2019 – Jan 2020.
 * Supervisor: Prof. Terry Wyatt.
 *
 * Usage (in ROOT):
 *   root [0] .L delta_t_minimization.C
 *   root [1] mini t("your_data.root", "mini")
 *   root [2] t.Loop()
 */

#define mini_cxx
#include "mini.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TError.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>
#include <vector>

using namespace std;

// ============================================================
// Physical constants (all masses in MeV)
// ============================================================

static const double MW         = 80379.0;   // W boson mass [MeV]
static const double MT         = 173000.0;  // Top quark mass [MeV]
static const double M_ELECTRON = 0.511;     // Electron mass [MeV]
static const double M_MUON     = 105.66;    // Muon mass [MeV]
static const double PT_CUT     = 25000.0;   // Minimum pT cut [MeV]

// ============================================================
// Minimization settings
// ============================================================

static const int    N_REPEATS      = 100;   // Random restarts per event
static const int    N_COMBINATIONS = 2;     // b-jet assignment combinations
static const double STEP_MOMENTUM  = 10.0;  // Step size for neutrino momenta

// ============================================================
// Utility: random double in [lower, upper]
// ============================================================

double getRandomDouble(double lower, double upper)
{
    static std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(lower, upper);
    return distribution(generator);
}

// ============================================================
// Four-vector utilities
// ============================================================

using FourVec = array<double, 4>;  // {E, px, py, pz}

FourVec makeFourVec(double E, double px, double py, double pz)
{
    return {E, px, py, pz};
}

FourVec addFourVecs(std::initializer_list<FourVec> vecs)
{
    FourVec result = {0, 0, 0, 0};
    for (const auto& v : vecs)
        for (int i = 0; i < 4; ++i)
            result[i] += v[i];
    return result;
}

double invariantMass(const FourVec& A)
{
    double m2 = A[0]*A[0] - A[1]*A[1] - A[2]*A[2] - A[3]*A[3];
    return (m2 >= 0) ? sqrt(m2) : 0.0;
}

FourVec masslessFourVec(double px, double py, double pz)
{
    return {sqrt(px*px + py*py + pz*pz), px, py, pz};
}

FourVec massiveFourVec(double px, double py, double pz, double mass)
{
    return {sqrt(px*px + py*py + pz*pz + mass*mass), px, py, pz};
}

void polarToCartesian(double pt, double phi, double eta,
                      double& px, double& py, double& pz)
{
    px = pt * cos(phi);
    py = pt * sin(phi);
    pz = pt * sinh(eta);
}

// ============================================================
// ΔT penalty function (tt1 topology — 2 neutrinos)
//
// Free variables (xx[0..5]):
//   xx[0..2] : three-momentum of ν_l1 (neutrino from top decay)
//   xx[3..5] : three-momentum of ν_l2 (neutrino from antitop decay)
//
// Fixed inputs (xx[6..23]):
//   Kinematic info of leptons l1, l2, MET, and two b-jets.
// ============================================================

double Tfunction(const double* xx)
{
    // --- Free variables ---
    const double neut_px_l1 = xx[0];
    const double neut_py_l1 = xx[1];
    const double neut_pz_l1 = xx[2];
    const double neut_px_l2 = xx[3];
    const double neut_py_l2 = xx[4];
    const double neut_pz_l2 = xx[5];

    // --- Fixed inputs: lepton 1 (from top decay) ---
    const double lep_pt_l1   = xx[6];
    const double lep_phi_l1  = xx[7];
    const double lep_eta_l1  = xx[8];
    const double lep_mass_l1 = xx[9];

    // --- Fixed inputs: lepton 2 (from antitop decay) ---
    const double lep_pt_l2   = xx[10];
    const double lep_phi_l2  = xx[11];
    const double lep_eta_l2  = xx[12];
    const double lep_mass_l2 = xx[13];

    // --- Fixed inputs: MET ---
    const double et_miss  = xx[14];
    const double phi_miss = xx[15];

    // --- Fixed inputs: b-jets ---
    const double bjet_E_1   = xx[16];
    const double bjet_pt_1  = xx[17];
    const double bjet_phi_1 = xx[18];
    const double bjet_eta_1 = xx[19];
    const double bjet_E_2   = xx[20];
    const double bjet_pt_2  = xx[21];
    const double bjet_phi_2 = xx[22];
    const double bjet_eta_2 = xx[23];

    // --- Coordinate conversion ---
    double lep_px_l1, lep_py_l1, lep_pz_l1;
    double lep_px_l2, lep_py_l2, lep_pz_l2;
    double bjet_px_1, bjet_py_1, bjet_pz_1;
    double bjet_px_2, bjet_py_2, bjet_pz_2;

    polarToCartesian(lep_pt_l1, lep_phi_l1, lep_eta_l1, lep_px_l1, lep_py_l1, lep_pz_l1);
    polarToCartesian(lep_pt_l2, lep_phi_l2, lep_eta_l2, lep_px_l2, lep_py_l2, lep_pz_l2);
    polarToCartesian(bjet_pt_1, bjet_phi_1, bjet_eta_1, bjet_px_1, bjet_py_1, bjet_pz_1);
    polarToCartesian(bjet_pt_2, bjet_phi_2, bjet_eta_2, bjet_px_2, bjet_py_2, bjet_pz_2);

    const double et_miss_x = et_miss * cos(phi_miss);
    const double et_miss_y = et_miss * sin(phi_miss);

    // --- Build four-vectors ---
    const FourVec p_neut_l1 = masslessFourVec(neut_px_l1, neut_py_l1, neut_pz_l1);
    const FourVec p_neut_l2 = masslessFourVec(neut_px_l2, neut_py_l2, neut_pz_l2);
    const FourVec p_l1      = massiveFourVec(lep_px_l1, lep_py_l1, lep_pz_l1, lep_mass_l1);
    const FourVec p_l2      = massiveFourVec(lep_px_l2, lep_py_l2, lep_pz_l2, lep_mass_l2);
    const FourVec p_bjet_1  = makeFourVec(bjet_E_1, bjet_px_1, bjet_py_1, bjet_pz_1);
    const FourVec p_bjet_2  = makeFourVec(bjet_E_2, bjet_px_2, bjet_py_2, bjet_pz_2);

    // --- Intermediate particles ---
    const FourVec p_w1 = addFourVecs({p_l1, p_neut_l1});           // W+
    const FourVec p_w2 = addFourVecs({p_l2, p_neut_l2});           // W-
    const FourVec p_t1 = addFourVecs({p_l1, p_neut_l1, p_bjet_1}); // top
    const FourVec p_t2 = addFourVecs({p_l2, p_neut_l2, p_bjet_2}); // antitop

    const double m_w1 = invariantMass(p_w1);
    const double m_w2 = invariantMass(p_w2);
    const double m_t1 = invariantMass(p_t1);
    const double m_t2 = invariantMass(p_t2);

    // --- ΔT² discriminant ---
    return pow(m_w1 - MW, 2) + pow(m_w2 - MW, 2) +
           pow(m_t1 - MT, 2) + pow(m_t2 - MT, 2) +
           pow(neut_px_l1 + neut_px_l2 - et_miss_x, 2) +
           pow(neut_py_l1 + neut_py_l2 - et_miss_y, 2);
}

// ============================================================
// Result struct — bundles all outputs from one minimization
// ============================================================

struct MinResult {
    double delta_T = std::numeric_limits<double>::max();
    double m_w1    = 0;
    double m_w2    = 0;
    double m_t1    = 0;
    double m_t2    = 0;
};

// ============================================================
// Minimise Tfunction for a single b-jet combination
// Returns the best MinResult found across N_REPEATS restarts
// ============================================================

MinResult MinimizeTfunction(
    double lep_pt_l1, double lep_phi_l1, double lep_eta_l1, double lep_mass_l1,
    double lep_pt_l2, double lep_phi_l2, double lep_eta_l2, double lep_mass_l2,
    double et_miss, double phi_miss,
    double bjet_E_1, double bjet_pt_1, double bjet_phi_1, double bjet_eta_1,
    double bjet_E_2, double bjet_pt_2, double bjet_phi_2, double bjet_eta_2,
    const char* minName = "Minuit2", const char* algoName = "")
{
    double lep_px_l1, lep_py_l1, lep_pz_l1;
    double lep_px_l2, lep_py_l2, lep_pz_l2;
    polarToCartesian(lep_pt_l1, lep_phi_l1, lep_eta_l1, lep_px_l1, lep_py_l1, lep_pz_l1);
    polarToCartesian(lep_pt_l2, lep_phi_l2, lep_eta_l2, lep_px_l2, lep_py_l2, lep_pz_l2);

    ROOT::Math::Minimizer* minimizer =
        ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetMaxIterations(10000);
    minimizer->SetTolerance(0.1);
    minimizer->SetPrintLevel(0);

    ROOT::Math::Functor f(&Tfunction, 24);
    minimizer->SetFunction(f);

    const double step = STEP_MOMENTUM;
    MinResult best;

    for (int rep = 0; rep < N_REPEATS; ++rep)
    {
        // Randomise starting neutrino momenta around lepton momenta
        minimizer->SetVariable(0, "neut_px_l1", getRandomDouble(-2*lep_px_l1, 2*lep_px_l1), step);
        minimizer->SetVariable(1, "neut_py_l1", getRandomDouble(-2*lep_py_l1, 2*lep_py_l1), step);
        minimizer->SetVariable(2, "neut_pz_l1", getRandomDouble(-2*lep_pz_l1, 2*lep_pz_l1), step);
        minimizer->SetVariable(3, "neut_px_l2", getRandomDouble(-2*lep_px_l2, 2*lep_px_l2), step);
        minimizer->SetVariable(4, "neut_py_l2", getRandomDouble(-2*lep_py_l2, 2*lep_py_l2), step);
        minimizer->SetVariable(5, "neut_pz_l2", getRandomDouble(-2*lep_pz_l2, 2*lep_pz_l2), step);

        // Fixed inputs
        minimizer->SetFixedVariable(6,  "lep_pt_l1",   lep_pt_l1);
        minimizer->SetFixedVariable(7,  "lep_phi_l1",  lep_phi_l1);
        minimizer->SetFixedVariable(8,  "lep_eta_l1",  lep_eta_l1);
        minimizer->SetFixedVariable(9,  "lep_mass_l1", lep_mass_l1);
        minimizer->SetFixedVariable(10, "lep_pt_l2",   lep_pt_l2);
        minimizer->SetFixedVariable(11, "lep_phi_l2",  lep_phi_l2);
        minimizer->SetFixedVariable(12, "lep_eta_l2",  lep_eta_l2);
        minimizer->SetFixedVariable(13, "lep_mass_l2", lep_mass_l2);
        minimizer->SetFixedVariable(14, "et_miss",     et_miss);
        minimizer->SetFixedVariable(15, "phi_miss",    phi_miss);
        minimizer->SetFixedVariable(16, "bjet_E_1",    bjet_E_1);
        minimizer->SetFixedVariable(17, "bjet_pt_1",   bjet_pt_1);
        minimizer->SetFixedVariable(18, "bjet_phi_1",  bjet_phi_1);
        minimizer->SetFixedVariable(19, "bjet_eta_1",  bjet_eta_1);
        minimizer->SetFixedVariable(20, "bjet_E_2",    bjet_E_2);
        minimizer->SetFixedVariable(21, "bjet_pt_2",   bjet_pt_2);
        minimizer->SetFixedVariable(22, "bjet_phi_2",  bjet_phi_2);
        minimizer->SetFixedVariable(23, "bjet_eta_2",  bjet_eta_2);

        minimizer->Minimize();

        double val = sqrt(minimizer->MinValue()); // ΔT = sqrt(ΔT²)
        if (val < best.delta_T)
        {
            const double* xs = minimizer->X();
            best.delta_T = val;

            // Reconstruct intermediate masses at the best-fit neutrino momenta
            double lep_px_l1c, lep_py_l1c, lep_pz_l1c;
            double lep_px_l2c, lep_py_l2c, lep_pz_l2c;
            double bjet_px_1c, bjet_py_1c, bjet_pz_1c;
            double bjet_px_2c, bjet_py_2c, bjet_pz_2c;
            polarToCartesian(lep_pt_l1, lep_phi_l1, lep_eta_l1, lep_px_l1c, lep_py_l1c, lep_pz_l1c);
            polarToCartesian(lep_pt_l2, lep_phi_l2, lep_eta_l2, lep_px_l2c, lep_py_l2c, lep_pz_l2c);
            polarToCartesian(bjet_pt_1, bjet_phi_1, bjet_eta_1, bjet_px_1c, bjet_py_1c, bjet_pz_1c);
            polarToCartesian(bjet_pt_2, bjet_phi_2, bjet_eta_2, bjet_px_2c, bjet_py_2c, bjet_pz_2c);

            FourVec p_nl1 = masslessFourVec(xs[0], xs[1], xs[2]);
            FourVec p_nl2 = masslessFourVec(xs[3], xs[4], xs[5]);
            FourVec p_l1  = massiveFourVec(lep_px_l1c, lep_py_l1c, lep_pz_l1c, lep_mass_l1);
            FourVec p_l2  = massiveFourVec(lep_px_l2c, lep_py_l2c, lep_pz_l2c, lep_mass_l2);
            FourVec p_b1  = makeFourVec(bjet_E_1, bjet_px_1c, bjet_py_1c, bjet_pz_1c);
            FourVec p_b2  = makeFourVec(bjet_E_2, bjet_px_2c, bjet_py_2c, bjet_pz_2c);

            best.m_w1 = invariantMass(addFourVecs({p_l1, p_nl1}));
            best.m_w2 = invariantMass(addFourVecs({p_l2, p_nl2}));
            best.m_t1 = invariantMass(addFourVecs({p_l1, p_nl1, p_b1}));
            best.m_t2 = invariantMass(addFourVecs({p_l2, p_nl2, p_b2}));
        }
    }

    delete minimizer;
    return best;
}

// ============================================================
// Main event loop
// ============================================================

void mini::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    // Histograms
    TH1F* h_delta_T  = new TH1F("h_delta_T",  "#DeltaT discriminant | #geq1e, #geq1#mu, #geq2 b-jets", 1000, 0, 400);
    TH1F* h_w1_mass  = new TH1F("h_w1_mass",  "Reconstructed W^{+} mass",   1000, 20000,  140000);
    TH1F* h_w2_mass  = new TH1F("h_w2_mass",  "Reconstructed W^{-} mass",   1000, 20000,  140000);
    TH1F* h_t1_mass  = new TH1F("h_t1_mass",  "Reconstructed top mass",     1000, 100000, 340000);
    TH1F* h_t2_mass  = new TH1F("h_t2_mass",  "Reconstructed antitop mass", 1000, 100000, 340000);

    for (Long64_t jentry = 0; jentry < 1000; jentry++)
    {
        GetEntry(jentry);

        // Count particle types
        int eCount = 0, muCount = 0, bjetCount = 0;
        for (size_t i = 0; i < lep_n; i++) {
            if (lep_type->at(i) == 11) eCount++;
            if (lep_type->at(i) == 13) muCount++;
        }
        for (size_t i = 0; i < jet_trueflav->size(); i++) {
            if (jet_trueflav->at(i) == 5) bjetCount++;
        }

        if (eCount < 1 || muCount < 1 || bjetCount < 2) continue;

        // Collect pT values
        double e_pt[10]    = {};
        double mu_pt[10]   = {};
        double bjet_pt[20] = {};

        for (size_t i = 0; i < lep_n; i++) {
            if (lep_type->at(i) == 11) e_pt[i]  = lep_pt->at(i);
            if (lep_type->at(i) == 13) mu_pt[i] = lep_pt->at(i);
        }
        for (size_t i = 0; i < jet_trueflav->size(); i++) {
            if (jet_trueflav->at(i) == 5) bjet_pt[i] = jet_pt->at(i);
        }

        // Find highest-pT candidates
        int e_index  = std::distance(std::begin(e_pt),
                           std::max_element(std::begin(e_pt), std::end(e_pt)));
        int mu_index = std::distance(std::begin(mu_pt),
                           std::max_element(std::begin(mu_pt), std::end(mu_pt)));
        int bjet_index_1 = std::distance(std::begin(bjet_pt),
                               std::max_element(std::begin(bjet_pt), std::end(bjet_pt)));

        int bjet_index_2 = 0;
        double bjet_pt_2nd = 0;
        for (int i = 0; i < 20; i++) {
            if (bjet_pt[i] == jet_pt->at(bjet_index_1)) continue;
            if (bjet_pt[i] > bjet_pt_2nd) {
                bjet_pt_2nd  = bjet_pt[i];
                bjet_index_2 = i;
            }
        }

        // Apply pT > 25 GeV and opposite-sign cuts
        bool pass_pt = (lep_pt->at(e_index)       >= PT_CUT &&
                        lep_pt->at(mu_index)       >= PT_CUT &&
                        jet_pt->at(bjet_index_1)   >= PT_CUT &&
                        jet_pt->at(bjet_index_2)   >= PT_CUT);
        bool opposite_sign = ((lep_charge->at(e_index) ==  1 && lep_charge->at(mu_index) == -1) ||
                              (lep_charge->at(e_index) == -1 && lep_charge->at(mu_index) ==  1));

        if (!pass_pt || !opposite_sign) continue;

        bool e_positive = (lep_charge->at(e_index) == 1);

        // Run minimization for both b-jet combinations and keep the best
        MinResult combo1, combo2;

        if (e_positive) {
            // positron + muon
            combo1 = MinimizeTfunction(
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                met_et, met_phi,
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2));

            combo2 = MinimizeTfunction(
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                met_et, met_phi,
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2),
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1));
        } else {
            // electron + antimuon
            combo1 = MinimizeTfunction(
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                met_et, met_phi,
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2));

            combo2 = MinimizeTfunction(
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                met_et, met_phi,
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2),
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1));
        }

        // Select the combination with the smaller ΔT
        const MinResult& best = (combo1.delta_T <= combo2.delta_T) ? combo1 : combo2;

        h_delta_T->Fill(best.delta_T);
        h_w1_mass->Fill(best.m_w1);
        h_w2_mass->Fill(best.m_w2);
        h_t1_mass->Fill(best.m_t1);
        h_t2_mass->Fill(best.m_t2);

    } // end event loop

    h_delta_T->Draw();
}
