/*
 * delta_t3_8combinations.C
 *
 * Computes the ΔT3 discriminant for tt̄ → bbτ_h τ_l event reconstruction.
 *
 * The ΔT3 penalty function reconstructs the "tt3" topology in which both W
 * bosons from the top pair decay via tau leptons, resulting in up to 6
 * neutrinos in the final state. The collinear approximation is used to combine
 * the secondary neutrino pairs into effective momenta proportional to the
 * observed light leptons (scaled by free parameters β₁ and β₂).
 *
 * For each event passing selection (≥1e, ≥1μ, ≥2 b-jets, all pT > 25 GeV,
 * opposite-sign leptons), minimization is performed 100 times across 8
 * possible lepton–bjet assignment combinations. The minimum ΔT3 over all
 * combinations and restarts is recorded.
 *
 * Part of the ATLAS Di-Higgs Discrimination Project.
 * University of Manchester, MPhys Project, Sep 2019 – Jan 2020.
 * Supervisor: Prof. Terry Wyatt.
 *
 * Builds on the ψ_T discriminant of Cookman & Harris (2019), which used ΔT
 * for the simpler tt1 topology (2 neutrinos). ΔT3 extends this to 6 neutrinos.
 *
 * Usage (in ROOT):
 *   root [0] .L delta_t3_8combinations.C
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

static const double MW        = 80379.0;   // W boson mass [MeV]
static const double MT        = 173000.0;  // Top quark mass [MeV]
static const double M_ELECTRON = 0.511;    // Electron mass [MeV]
static const double M_MUON    = 105.66;    // Muon mass [MeV]
static const double PT_CUT    = 25000.0;   // Minimum pT cut [MeV]

// ============================================================
// Minimization settings
// ============================================================

static const int    N_REPEATS       = 100;   // Random restarts per event
static const int    N_COMBINATIONS  = 8;     // Lepton-bjet assignment combos
static const double STEP_MOMENTUM   = 10.0;  // Step size for neutrino momenta
static const double STEP_BETA       = 0.001; // Step size for beta scaling

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

// Add any number of four-vectors
FourVec addFourVecs(std::initializer_list<FourVec> vecs)
{
    FourVec result = {0, 0, 0, 0};
    for (const auto& v : vecs)
        for (int i = 0; i < 4; ++i)
            result[i] += v[i];
    return result;
}

// Minkowski dot product (metric signature +---)
double dot(const FourVec& A, const FourVec& B)
{
    return A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
}

// Invariant mass of a four-vector
double invariantMass(const FourVec& A)
{
    double m2 = dot(A, A);
    return (m2 >= 0) ? sqrt(m2) : 0.0;
}

// Build a massless four-vector from three-momentum components
FourVec masslessFourVec(double px, double py, double pz)
{
    double E = sqrt(px*px + py*py + pz*pz);
    return {E, px, py, pz};
}

// Build a massive four-vector from three-momentum and mass
FourVec massiveFourVec(double px, double py, double pz, double mass)
{
    double E = sqrt(px*px + py*py + pz*pz + mass*mass);
    return {E, px, py, pz};
}

// Coordinate conversion: polar (pT, phi, eta) -> Cartesian (px, py, pz)
void polarToCartesian(double pt, double phi, double eta,
                      double& px, double& py, double& pz)
{
    px = pt * cos(phi);
    py = pt * sin(phi);
    pz = pt * sinh(eta);
}

// ============================================================
// ΔT3 penalty function
//
// Free variables (xx[0..6]):
//   xx[0..2] : three-momentum of ν_l1  (neutrino from prompt W decay)
//   xx[3..5] : three-momentum of ν_τ1  (tau neutrino from prompt W decay)
//   xx[6]    : β (scaling factor for combined ν_τ2 + ν_l2 momentum)
//
// Fixed inputs (xx[7..24]):
//   Kinematic info of leptons l1, l2, MET, and two b-jets (all constants
//   passed as fixed variables through the Minuit2 interface).
// ============================================================

double T3Function(const double* xx)
{
    // --- Free variables ---
    const double neut_px_l1   = xx[0];
    const double neut_py_l1   = xx[1];
    const double neut_pz_l1   = xx[2];
    const double neut_px_tau1 = xx[3];
    const double neut_py_tau1 = xx[4];
    const double neut_pz_tau1 = xx[5];
    const double beta         = xx[6];  // must be > 0

    // --- Fixed inputs: lepton 1 (from prompt W) ---
    const double lep_pt_l1   = xx[7];
    const double lep_phi_l1  = xx[8];
    const double lep_eta_l1  = xx[9];
    const double lep_mass_l1 = xx[10];

    // --- Fixed inputs: lepton 2 (from secondary W) ---
    const double lep_pt_l2   = xx[11];
    const double lep_phi_l2  = xx[12];
    const double lep_eta_l2  = xx[13];
    const double lep_mass_l2 = xx[14];

    // --- Fixed inputs: MET ---
    const double et_miss  = xx[15];
    const double phi_miss = xx[16];

    // --- Fixed inputs: b-jet 1 (paired with prompt W top) ---
    const double bjet_E_1   = xx[17];
    const double bjet_pt_1  = xx[18];
    const double bjet_phi_1 = xx[19];
    const double bjet_eta_1 = xx[20];

    // --- Fixed inputs: b-jet 2 (paired with secondary W top) ---
    const double bjet_E_2   = xx[21];
    const double bjet_pt_2  = xx[22];
    const double bjet_phi_2 = xx[23];
    const double bjet_eta_2 = xx[24];

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

    // --- Build four-vectors for final state particles ---

    // Combined four-vector for l2 + ν_τ2 + ν_l2 (collinear approximation)
    // The secondary neutrinos carry momentum β * p_l2
    const double E_l2          = sqrt(lep_px_l2*lep_px_l2 + lep_py_l2*lep_py_l2 +
                                       lep_pz_l2*lep_pz_l2 + lep_mass_l2*lep_mass_l2);
    const double E_vtau2_vl2   = sqrt(beta*beta * (lep_px_l2*lep_px_l2 +
                                       lep_py_l2*lep_py_l2 + lep_pz_l2*lep_pz_l2));
    const FourVec p_l2_vtau2_vl2 = makeFourVec(E_l2 + E_vtau2_vl2,
                                                (1+beta)*lep_px_l2,
                                                (1+beta)*lep_py_l2,
                                                (1+beta)*lep_pz_l2);

    const FourVec p_neut_l1   = masslessFourVec(neut_px_l1, neut_py_l1, neut_pz_l1);
    const FourVec p_neut_tau1 = masslessFourVec(neut_px_tau1, neut_py_tau1, neut_pz_tau1);
    const FourVec p_l1        = massiveFourVec(lep_px_l1, lep_py_l1, lep_pz_l1, lep_mass_l1);
    const FourVec p_bjet_1    = makeFourVec(bjet_E_1, bjet_px_1, bjet_py_1, bjet_pz_1);
    const FourVec p_bjet_2    = makeFourVec(bjet_E_2, bjet_px_2, bjet_py_2, bjet_pz_2);

    // --- Intermediate particles ---
    const FourVec p_w1_1 = addFourVecs({p_l1, p_neut_l1});                    // W from top 1
    const FourVec p_w1_2 = addFourVecs({p_neut_tau1, p_l2_vtau2_vl2});        // W from top 2
    const FourVec p_t_1  = addFourVecs({p_l1, p_neut_l1, p_bjet_1});          // Top quark 1
    const FourVec p_t_2  = addFourVecs({p_bjet_2, p_neut_tau1, p_l2_vtau2_vl2}); // Top quark 2

    const double m_w1_1 = invariantMass(p_w1_1);
    const double m_w1_2 = invariantMass(p_w1_2);
    const double m_t_1  = invariantMass(p_t_1);
    const double m_t_2  = invariantMass(p_t_2);

    // --- ΔT3 discriminant ---
    // Penalises deviations of reconstructed masses from known W and top masses,
    // and deviations of reconstructed MET from measured MET.
    return sqrt(
        pow(m_w1_1 - MW, 2) +
        pow(m_w1_2 - MW, 2) +
        pow(m_t_1  - MT, 2) +
        pow(m_t_2  - MT, 2) +
        pow(neut_px_l1 + neut_px_tau1 + beta*lep_px_l2 - et_miss_x, 2) +
        pow(neut_py_l1 + neut_py_tau1 + beta*lep_py_l2 - et_miss_y, 2)
    );
}

// ============================================================
// Minimise T3Function for a single lepton-bjet combination
// Returns the minimum ΔT3 value found
// ============================================================

double MinimizeT3function(
    double lep_pt_l1, double lep_phi_l1, double lep_eta_l1, double lep_mass_l1,
    double lep_pt_l2, double lep_phi_l2, double lep_eta_l2, double lep_mass_l2,
    double et_miss, double phi_miss,
    double bjet_E_1, double bjet_pt_1, double bjet_phi_1, double bjet_eta_1,
    double bjet_E_2, double bjet_pt_2, double bjet_phi_2, double bjet_eta_2,
    const char* minName = "Minuit2", const char* algoName = "")
{
    // Coordinate conversion for initialising starting points
    double lep_px_l1, lep_py_l1, lep_pz_l1;
    double lep_px_l2, lep_py_l2, lep_pz_l2;
    polarToCartesian(lep_pt_l1, lep_phi_l1, lep_eta_l1, lep_px_l1, lep_py_l1, lep_pz_l1);
    polarToCartesian(lep_pt_l2, lep_phi_l2, lep_eta_l2, lep_px_l2, lep_py_l2, lep_pz_l2);

    ROOT::Math::Minimizer* minimizer =
        ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    minimizer->SetMaxFunctionCalls(100000000);
    minimizer->SetMaxIterations(1000000);
    minimizer->SetTolerance(0.1);
    minimizer->SetPrintLevel(0);  // suppress per-event output

    ROOT::Math::Functor f(&T3Function, 25);
    minimizer->SetFunction(f);

    const double step[7] = {STEP_MOMENTUM, STEP_MOMENTUM, STEP_MOMENTUM,
                             STEP_MOMENTUM, STEP_MOMENTUM, STEP_MOMENTUM,
                             STEP_BETA};

    double bestDeltaT3 = std::numeric_limits<double>::max();

    for (int rep = 0; rep < N_REPEATS; ++rep)
    {
        // Randomise starting neutrino momenta around lepton momenta
        double var[7] = {
            getRandomDouble(-2*lep_px_l1, 2*lep_px_l1),
            getRandomDouble(-2*lep_py_l1, 2*lep_py_l1),
            getRandomDouble(-2*lep_pz_l1, 2*lep_pz_l1),
            getRandomDouble(-2*lep_px_l2, 2*lep_px_l2),
            getRandomDouble(-2*lep_py_l2, 2*lep_py_l2),
            getRandomDouble(-2*lep_pz_l2, 2*lep_pz_l2),
            getRandomDouble(0, 1)   // beta > 0
        };

        // Free variables
        minimizer->SetVariable(0, "neut_px_l1",   var[0], step[0]);
        minimizer->SetVariable(1, "neut_py_l1",   var[1], step[1]);
        minimizer->SetVariable(2, "neut_pz_l1",   var[2], step[2]);
        minimizer->SetVariable(3, "neut_px_tau1", var[3], step[3]);
        minimizer->SetVariable(4, "neut_py_tau1", var[4], step[4]);
        minimizer->SetVariable(5, "neut_pz_tau1", var[5], step[5]);
        minimizer->SetLimitedVariable(6, "beta",  var[6], step[6], 0, 100);

        // Fixed inputs (passed through variable array as constants)
        minimizer->SetFixedVariable(7,  "lep_pt_l1",   lep_pt_l1);
        minimizer->SetFixedVariable(8,  "lep_phi_l1",  lep_phi_l1);
        minimizer->SetFixedVariable(9,  "lep_eta_l1",  lep_eta_l1);
        minimizer->SetFixedVariable(10, "lep_mass_l1", lep_mass_l1);
        minimizer->SetFixedVariable(11, "lep_pt_l2",   lep_pt_l2);
        minimizer->SetFixedVariable(12, "lep_phi_l2",  lep_phi_l2);
        minimizer->SetFixedVariable(13, "lep_eta_l2",  lep_eta_l2);
        minimizer->SetFixedVariable(14, "lep_mass_l2", lep_mass_l2);
        minimizer->SetFixedVariable(15, "et_miss",     et_miss);
        minimizer->SetFixedVariable(16, "phi_miss",    phi_miss);
        minimizer->SetFixedVariable(17, "bjet_E_1",    bjet_E_1);
        minimizer->SetFixedVariable(18, "bjet_pt_1",   bjet_pt_1);
        minimizer->SetFixedVariable(19, "bjet_phi_1",  bjet_phi_1);
        minimizer->SetFixedVariable(20, "bjet_eta_1",  bjet_eta_1);
        minimizer->SetFixedVariable(21, "bjet_E_2",    bjet_E_2);
        minimizer->SetFixedVariable(22, "bjet_pt_2",   bjet_pt_2);
        minimizer->SetFixedVariable(23, "bjet_phi_2",  bjet_phi_2);
        minimizer->SetFixedVariable(24, "bjet_eta_2",  bjet_eta_2);

        minimizer->Minimize();

        double delta_T3 = minimizer->MinValue();
        if (delta_T3 < bestDeltaT3)
            bestDeltaT3 = delta_T3;
    }

    delete minimizer;
    return bestDeltaT3;
}

// ============================================================
// Main event loop
// ============================================================

void mini::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    // Histogram: ΔT3 discriminant distribution
    TH1F* h_delta_T3 = new TH1F("h_delta_T3",
        "#DeltaT_{3} discriminant | #geq1e, #geq1#mu, #geq2 b-jets",
        1000, 0, 100000);

    // Loop over events
    for (Long64_t jentry = 0; jentry < 1000; jentry++)
    {
        GetEntry(jentry);

        // Count particle types in this event
        int eCount    = 0;
        int muCount   = 0;
        int bjetCount = 0;

        for (size_t i = 0; i < lep_n; i++) {
            if (lep_type->at(i) == 11) eCount++;
            if (lep_type->at(i) == 13) muCount++;
        }
        for (size_t i = 0; i < jet_trueflav->size(); i++) {
            if (jet_trueflav->at(i) == 5) bjetCount++;
        }

        // Require at least 1 electron, 1 muon, 2 b-jets
        if (eCount < 1 || muCount < 1 || bjetCount < 2) continue;

        // Collect pT values to find highest-pT candidates
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

        // Find indices of highest-pT electron and muon
        int e_index  = std::distance(std::begin(e_pt),
                           std::max_element(std::begin(e_pt), std::end(e_pt)));
        int mu_index = std::distance(std::begin(mu_pt),
                           std::max_element(std::begin(mu_pt), std::end(mu_pt)));

        // Find indices of two highest-pT b-jets
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

        // Apply pT > 25 GeV cut and require opposite-sign lepton pair
        bool pass_pt = (lep_pt->at(e_index)  >= PT_CUT &&
                        lep_pt->at(mu_index) >= PT_CUT &&
                        jet_pt->at(bjet_index_1) >= PT_CUT &&
                        jet_pt->at(bjet_index_2) >= PT_CUT);

        bool opposite_sign = ((lep_charge->at(e_index) ==  1 && lep_charge->at(mu_index) == -1) ||
                              (lep_charge->at(e_index) == -1 && lep_charge->at(mu_index) ==  1));

        if (!pass_pt || !opposite_sign) continue;

        // Determine which lepton is positron/antimuon to set combination ordering
        bool e_positive = (lep_charge->at(e_index) == 1);

        // Run minimization for all 8 lepton-bjet combinations and collect results.
        // Combinations differ by:
        //   - which lepton is assigned to the "prompt" W (l1) vs "secondary" W (l2)
        //   - which b-jet is paired with which top quark
        // See project report for the full combination table.
        double combo_results[N_COMBINATIONS];
        for (int i = 0; i < N_COMBINATIONS; i++)
            combo_results[i] = std::numeric_limits<double>::max();

        if (e_positive) {
            // Positron + muon: combos 1-4
            combo_results[0] = MinimizeT3function(
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                met_et, met_phi,
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2),
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1));

            combo_results[1] = MinimizeT3function(
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                met_et, met_phi,
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2));

            combo_results[2] = MinimizeT3function(
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                met_et, met_phi,
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2));

            combo_results[3] = MinimizeT3function(
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                met_et, met_phi,
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2),
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1));
        } else {
            // Electron + antimuon: combos 5-8
            combo_results[4] = MinimizeT3function(
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                met_et, met_phi,
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2),
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1));

            combo_results[5] = MinimizeT3function(
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                met_et, met_phi,
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2));

            combo_results[6] = MinimizeT3function(
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                met_et, met_phi,
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2));

            combo_results[7] = MinimizeT3function(
                lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), M_MUON,
                lep_pt->at(e_index),  lep_phi->at(e_index),  lep_eta->at(e_index),  M_ELECTRON,
                met_et, met_phi,
                jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2),
                jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2),
                jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1),
                jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1));
        }

        // Find the minimum ΔT3 across all combinations (fix: was incorrectly
        // using max_element in the original code)
        double delta_T3 = *std::min_element(std::begin(combo_results),
                                             std::end(combo_results));

        h_delta_T3->Fill(delta_T3);

    } // end event loop

    h_delta_T3->Draw();
}
