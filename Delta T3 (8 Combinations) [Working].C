#define mini_cxx
#include "mini.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <cstdlib>
#include <tuple>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>     /* abs */
#include <limits>
#include <random>
#include <array>
#include <vector>

using namespace std;


// Function to generate a random double number between two limits
double getRandomDouble(double lower_bound,  double upper_bound)
{
   static std::default_random_engine generator;
   std::uniform_real_distribution<double> distribution(lower_bound,upper_bound);

   return distribution(generator) ;
}

// Function to create an energy-momentum four vector for a particle
array<double, 4> createMomentumFourVector(double energy, double px, double py, double pz)
{
	array<double, 4> fourVector = {energy, px, py, pz}; 
	return fourVector; 
}

// Functions to add two/three/four/five four vectors (Overloaded)
array<double, 4> addFourVectors(array<double, 4> A, array<double, 4> B)
{
	array<double, 4> result = { A[0] + B[0], A[1] + B[1], A[2] + B[2], A[3] + B[3] };
	return result;
}
array<double, 4> addFourVectors(array<double, 4> A, array<double, 4> B, array<double, 4> C)
{
	array<double, 4> result = { A[0] + B[0] + C[0], A[1] + B[1]+C[1], A[2] + B[2]+C[2], A[3] + B[3]+C[3] };
	return result;
}
array<double, 4> addFourVectors(array<double, 4> A, array<double, 4> B, array<double, 4> C, array<double, 4> D)
{
	array<double, 4> result = { A[0] + B[0] + C[0]+D[0], A[1] + B[1]+C[1]+D[1], A[2] + B[2]+C[2]+D[2], A[3] + B[3]+C[3]+D[3] };
	return result;
}
array<double, 4> addFourVectors(array<double, 4> A, array<double, 4> B, array<double, 4> C, array<double, 4> D, array<double, 4> E)
{
	array<double, 4> result = { A[0]+B[0]+C[0]+D[0]+E[0], A[1]+B[1]+C[1]+D[1]+E[1], A[2]+B[2]+C[2]+D[2]+E[2], A[3]+B[3]+C[3]+D[3]+E[3] };
	return result;
}


// Function to take dot product of two four vectors
double dotProduct(array<double, 4> A, array<double, 4> B)
{
	double result = (A[0] * B[0]) - (A[1] * B[1]) - (A[2] * B[2]) - (A[3] * B[3]);
	return result; 
}

// Function to calculate invariant mass of an energy-momentum four vector
double getInvariantMass(array<double, 4> A)
{
	return sqrt ( (A[0]*A[0]) - (A[1]*A[1]) - (A[2]*A[2]) - (A[3]*A[3]) );
}

// Function to calculate invariant W mass (using Eqn 3.4)
// IF given mass of lepton lep_mass, three-momentum of lepton lep_px, lep_py, lep_pz, and three-momentum of neutrino
//    neut_px, neut_py, neut_pz
double getWMass(double lep_mass, double lep_px, double lep_py, double lep_pz, double neut_px, double neut_py,
	double neut_pz)
{
	double A = lep_px*lep_px + lep_py*lep_py + lep_pz*lep_pz + lep_mass*lep_mass ;
	double B = neut_px*neut_px + neut_py*neut_py + neut_pz*neut_pz ; 
	double w_mass = sqrt( pow( (sqrt(A) + sqrt(B)) ,2) - pow( (lep_px + neut_px) ,2) - pow( (lep_py + neut_py) ,2) 
		- pow( (lep_pz + neut_pz) ,2) );

	return w_mass;
}

// Function to calculate invariant T mass
// Using Eqn (5) Daniel Cookman Report
double getTMass(double bjet_E, double bjet_px, double bjet_py, double bjet_pz, double lep_mass, double lep_px,
  double lep_py, double lep_pz, double neut_px, double neut_py, double neut_pz)
{
  double lep_E = sqrt(lep_px*lep_px + lep_py*lep_py + lep_pz*lep_pz + lep_mass*lep_mass);
  double neut_E = sqrt(neut_px*neut_px + neut_py*neut_py + neut_pz*neut_pz);
  double TMassSquared = pow( bjet_E + neut_E + lep_E ,2) - pow( bjet_px + lep_px + neut_px ,2) 
                        - pow( bjet_py + lep_py + neut_py ,2) - pow( bjet_pz + lep_pz + neut_pz ,2);

  return sqrt(TMassSquared);
}

double T3Function(const double *xx )
{
  // Declare Variables (Momenta of the Four Neutrinos)
  const Double_t neut_px_l1 = xx[0]; // Three-Momentum of electron/muon neutrino from prompt W decay (v_l1) 
  const Double_t neut_py_l1 = xx[1]; // 
  const Double_t neut_pz_l1 = xx[2]; // (Variables)

  const Double_t neut_px_tau1 = xx[3]; // Three-Momentum of tau neutrino from prompt W decay (v_tau1)
  const Double_t neut_py_tau1 = xx[4]; //
  const Double_t neut_pz_tau1 = xx[5]; // (Variables)

  const Double_t beta = xx[6];        // Proportion of three momentum carried by v_tau2 and v_l2 compared to l2
                                      // beta > 0
  // End of Declare Variables

  // Declare Constants (Momenta of Electron, Muon, B Jets, Missing Energy)
  const Double_t lep_pt_l1 = xx[7];   // Kinematic Info of Electron/Muon from prompt W decay 
  const Double_t lep_phi_l1 = xx[8];  // 
  const Double_t lep_eta_l1 = xx[9];  // 
  const Double_t lep_mass_l1 = xx[10]; // (Constants)

  const Double_t lep_pt_l2 = xx[11];  // Kinematic Info of Electron/Muon from secondary W decay 
  const Double_t lep_phi_l2 = xx[12]; // 
  const Double_t lep_eta_l2 = xx[13]; // 
  const Double_t lep_mass_l2 = xx[14]; // (Constants)

  const Double_t et_miss = xx[15];  // Missing Transverse Energy (Constants)
  const Double_t phi_miss = xx[16]; // Missing Transverse Phi (Constants)

  const Double_t bjet_E_1 = xx[17]; // Kinematic Info of Bjet that accompanies only prompt W decay (Constants)
  const Double_t bjet_pt_1 = xx[18];
  const Double_t bjet_phi_1 = xx[19];
  const Double_t bjet_eta_1 = xx[20];

  const Double_t bjet_E_2 = xx[21]; // Kinematic Info of Bjet that accompanies secondary W decay (Constants)
  const Double_t bjet_pt_2 = xx[22];
  const Double_t bjet_phi_2 = xx[23];
  const Double_t bjet_eta_2 = xx[24];
  // End of Declare Constants 

  // Coordinate Conversion: Polar to Cartesian
  const Double_t lep_px_l1 = lep_pt_l1 * cos(lep_phi_l1) ; //  Three-Momentum of Electron in Cartesian (Constants)
  const Double_t lep_py_l1 = lep_pt_l1 * sin(lep_phi_l1) ; //
  const Double_t lep_pz_l1 = lep_pt_l1 * sinh(lep_eta_l1) ; //

  const Double_t lep_px_l2 = lep_pt_l2 * cos(lep_phi_l2) ; //  Three-Momentum of Muon in Cartesian (Constants)
  const Double_t lep_py_l2 = lep_pt_l2 * sin(lep_phi_l2) ; // 
  const Double_t lep_pz_l2 = lep_pt_l2 * sinh(lep_eta_l2) ; // 

  const Double_t et_miss_x = et_miss * cos(phi_miss); //  Two-Momentum of Missing Tranverse Energy in Cartesian (Constants)
  const Double_t et_miss_y = et_miss * sin(phi_miss); //

  const Double_t bjet_px_1 = bjet_pt_1 * cos(bjet_phi_1) ; // Three-Momentum of First Bjet in Cartesian 
  const Double_t bjet_py_1 = bjet_pt_1 * sin(bjet_phi_1) ;
  const Double_t bjet_pz_1 = bjet_pt_1 * sinh(bjet_eta_1) ; 

  const Double_t bjet_px_2 = bjet_pt_2 * cos(bjet_phi_2) ; // Three-Momentum of Second Bjet in Cartesian 
  const Double_t bjet_py_2 = bjet_pt_2 * sin(bjet_phi_2) ;
  const Double_t bjet_pz_2 = bjet_pt_2 * sinh(bjet_eta_2) ; 
  // End of Coordinate Conversion

  // Create combined four vector for particles l2, v_tau2, v_l2
  const Double_t E_l2 = sqrt(lep_px_l2*lep_px_l2 + lep_py_l2*lep_py_l2 + lep_pz_l2*lep_pz_l2 + lep_mass_l2*lep_mass_l2); // Energy of secondary lepton l2
  const Double_t E_vtau2_vl2 = sqrt(beta*beta*lep_px_l2*lep_px_l2 + beta*beta*lep_py_l2*lep_py_l2 + beta*beta*lep_pz_l2*lep_pz_l2); // Energy of secondary neutrinos v_tau2 and v_l2 
  const array<double, 4> p_l2_vtau2_vl2 = createMomentumFourVector( E_l2 + E_vtau2_vl2, 
    (1+beta) * lep_px_l2 , (1+beta) * lep_py_l2, (1+beta) * lep_pz_l2 ); // combined particle of l2, v_tau2, v_l2


  // Create four vectors for final products in decay (bbllvvvv)
  const array<double, 4> p_neut_l1 = createMomentumFourVector( sqrt(neut_px_l1*neut_px_l1 + neut_py_l1*neut_py_l1 + 
  	neut_pz_l1*neut_pz_l1) , neut_px_l1, neut_py_l1, neut_pz_l1 );  // electron/muon neutrino from prompt W decay (v_l1) 
  const array<double, 4> p_neut_tau1 = createMomentumFourVector( sqrt(neut_px_tau1*neut_px_tau1 + neut_py_tau1*neut_py_tau1 + 
  	neut_pz_tau1*neut_pz_tau1) , neut_px_tau1, neut_py_tau1, neut_pz_tau1 ); //  tau neutrino from prompt W decay (v_tau1)

  const array<double, 4> p_l1 = createMomentumFourVector( sqrt(lep_px_l1*lep_px_l1 + lep_py_l1*lep_py_l1 + 
  	lep_pz_l1*lep_pz_l1 + lep_mass_l1*lep_mass_l1) , lep_px_l1, lep_py_l1, lep_pz_l1 ); // Electron/Muon from prompt W decay
  const array<double, 4> p_l2 = createMomentumFourVector( sqrt(lep_px_l2*lep_px_l2 + lep_py_l2*lep_py_l2 + 
  	lep_pz_l2*lep_pz_l2 + lep_mass_l2*lep_mass_l2) , lep_px_l2, lep_py_l2, lep_pz_l2); // Electron/Muon from secondary W decay 
  
  const array<double, 4> p_bjet_1 = createMomentumFourVector( bjet_E_1, bjet_px_1, bjet_py_1, bjet_pz_1); // Bjet that accompanies only prompt W decay
  const array<double, 4> p_bjet_2 = createMomentumFourVector( bjet_E_2, bjet_px_2, bjet_py_2, bjet_pz_2); // Bjet that accompanies secondary W decay
  // End of Creation of four vectors for final products

  // Calculate four vectors for intermediate products in decay ( w1_1, w1_2, t_1, t_2)
  const array<double, 4> p_w1_1 = addFourVectors(p_l1, p_neut_l1);             // primary W boson that does not decay further
  const array<double, 4> p_w1_2 = addFourVectors(p_neut_tau1, p_l2_vtau2_vl2); // primary W boson that decays further to secondary W boson
  const array<double, 4> p_t_1 = addFourVectors(p_l1, p_neut_l1, p_bjet_1);    //  top quark that only decays to primay W 
  const array<double, 4> p_t_2 = addFourVectors(p_bjet_2, p_neut_tau1, p_l2_vtau2_vl2); // top quark that decays to primary and secondary W bosons
  // End of Caculate four vectors for intermediate products

  // Calculate invariant masses of intermediate products (tau, w1_1, w1_2, w2, t_1, t_2)
  const Double_t m_w1_1 = getInvariantMass(p_w1_1);
  const Double_t m_w1_2 = getInvariantMass(p_w1_2);
  const Double_t m_t_1 = getInvariantMass(p_t_1);
  const Double_t m_t_2 = getInvariantMass(p_t_2);
  // End of Cacluate invariant masses of intermediate products

  // Return discriminant delta_T3 
  return sqrt(pow(m_w1_1-80379, 2) + pow(m_w1_2-80379, 2) + pow(m_t_1-173000, 2) + 
  		 pow(m_t_2-173000, 2) + pow(neut_px_l1 + neut_px_tau1 + (beta*lep_px_l2) - et_miss_x, 2) + 
  		 pow(neut_py_l1 + neut_py_tau1 + (beta*lep_py_l2) - et_miss_y, 2));

}

std::tuple<double>  MinimizeT3function(
	double lep_pt_l1, double lep_phi_l1, double lep_eta_l1, double lep_mass_l1, 
	double lep_pt_l2, double lep_phi_l2, double lep_eta_l2, double lep_mass_l2, 
	double et_miss,  double phi_miss, 
	double bjet_E_1, double bjet_pt_1, double bjet_phi_1, double bjet_eta_1, 
	double bjet_E_2, double bjet_pt_2, double bjet_phi_2, double bjet_eta_2, 
	const char * minName = "Minuit2", const char *algoName = "")
{
   // create minimizer giving a name and a name (optionally) for the specific
   // algorithm
   // possible choices are:
   //     minName                  algoName
   // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
   //  Minuit2                     Fumili2
   //  Fumili
   //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
   //                              BFGS2, SteepestDescent
   //  GSLMultiFit
   //   GSLSimAn
   //   Genetic
   ROOT::Math::Minimizer* minimum =
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);

   // set tolerance , etc...
   minimum->SetMaxFunctionCalls(100000000); // for Minuit/Minuit2
   minimum->SetMaxIterations(1000000);  // for GSL
   minimum->SetTolerance(0.1);
   minimum->SetPrintLevel(1);

   // Transformation of Coordinates for Input Values (Polar to Cartesian)
  double lep_px_l1 = lep_pt_l1 * cos(lep_phi_l1) ; 
  double lep_py_l1 = lep_pt_l1 * sin(lep_phi_l1) ; 
  double lep_pz_l1 = lep_pt_l1 * sinh(lep_eta_l1) ; 

  double lep_px_l2 = lep_pt_l2 * cos(lep_phi_l2) ; 
  double lep_py_l2 = lep_pt_l2 * sin(lep_phi_l2) ;  
  double lep_pz_l2 = lep_pt_l2 * sinh(lep_eta_l2) ; 

  double et_miss_x = et_miss * cos(phi_miss); 
  double et_miss_y = et_miss * sin(phi_miss); 

  double bjet_px_1 = bjet_pt_1 * cos(bjet_phi_1) ; 
  double bjet_py_1 = bjet_pt_1 * sin(bjet_phi_1) ;
  double bjet_pz_1 = bjet_pt_1 * sinh(bjet_eta_1) ; 

  double bjet_px_2 = bjet_pt_2 * cos(bjet_phi_2) ; 
  double bjet_py_2 = bjet_pt_2 * sin(bjet_phi_2) ;
  double bjet_pz_2 = bjet_pt_2 * sinh(bjet_eta_2) ; 
  // End of Transformation of Coordinates 


   // create function wrapper for minimizer
   // a IMultiGenFunction type
   ROOT::Math::Functor f(&T3Function, 25);

   // Minimization Step Size 
   double step[7] = {10,10,10,10,10,10,0.001};

   // Minimization Starting Point
   double variable[7] = {}; 
   // Randomize the starting momenta of the four neutrinos 
   // ***NOTE: MIGHT NOT BE THE BESY WAY***
    variable[0] = getRandomDouble(-2*lep_px_l1, 2*lep_px_l1);
    variable[1] =getRandomDouble(-2*lep_py_l1,2*lep_py_l1); 
    variable[2] = getRandomDouble(-2*lep_pz_l1,2*lep_pz_l1);
    variable[3] = getRandomDouble(-2*lep_px_l2,2*lep_px_l2);
    variable[4] = getRandomDouble(-2*lep_py_l2,2*lep_py_l2);
    variable[5] = getRandomDouble(-2*lep_pz_l2,2*lep_pz_l2);
    variable[6] = getRandomDouble(0,1);
 
  //
   minimum->SetFunction(f);

   // Set the free variables to be minimized !
   minimum->SetVariable(0, "neut_px_l1", variable[0], step[0]);
   minimum->SetVariable(1, "neut_py_l1", variable[1], step[1] );
   minimum->SetVariable(2, "neut_pz_l1", variable[2], step[2] );

   minimum->SetVariable(3, "neut_px_tau1", variable[3], step[3]);
   minimum->SetVariable(4, "neut_py_tau1", variable[4], step[4] );
   minimum->SetVariable(5, "neut_pz_tau1", variable[5], step[5] );

   minimum->SetLimitedVariable(6, "beta", variable[6], step[6], 0, 100);
   // Set the fixed variables !

   minimum->SetFixedVariable(7, "lep_pt_l1", lep_pt_l1);
   minimum->SetFixedVariable(8, "lep_phi_l1",lep_phi_l1);
   minimum->SetFixedVariable(9, "lep_eta_l1",lep_eta_l1);
   minimum->SetFixedVariable(10, "lep_mass_l1",lep_mass_l1);

   minimum->SetFixedVariable(11, "lep_pt_l2",lep_pt_l2);
   minimum->SetFixedVariable(12, "lep_phi_l2",lep_phi_l2);
   minimum->SetFixedVariable(13, "lep_eta_l2",lep_eta_l2);
   minimum->SetFixedVariable(14, "lep_mass_l2",lep_mass_l2);

   minimum->SetFixedVariable(15, "et_miss",et_miss);
   minimum->SetFixedVariable(16, "phi_miss",phi_miss);

   minimum->SetFixedVariable(17, "bjet_E_1",bjet_E_1);
   minimum->SetFixedVariable(18, "bjet_pt_1",bjet_pt_1);
   minimum->SetFixedVariable(19, "bjet_phi_1",bjet_phi_1);
   minimum->SetFixedVariable(20, "bjet_eta_1",bjet_eta_1);

   minimum->SetFixedVariable(21, "bjet_E_2",bjet_E_2);
   minimum->SetFixedVariable(22, "bjet_pt_2",bjet_pt_2);
   minimum->SetFixedVariable(23, "bjet_phi_2",bjet_phi_2);
   minimum->SetFixedVariable(24, "bjet_eta_2",bjet_eta_2);



   // Do the minimization !
   minimum->Minimize();
   const double *xs = minimum->X();

   // Determine four vectors for final products  after minimization 
    double E_l2 = sqrt(lep_px_l2*lep_px_l2 + lep_py_l2*lep_py_l2 + lep_pz_l2*lep_pz_l2 + lep_mass_l2*lep_mass_l2); // Energy of secondary lepton l2
    double E_vtau2_vl2 = sqrt(xs[6]*xs[6]*lep_px_l2*lep_px_l2 + xs[6]*xs[6]*lep_py_l2*lep_py_l2 + xs[6]*xs[6]*lep_pz_l2*lep_pz_l2); // Energy of secondary neutrinos v_tau2 and v_l2 
    array<double, 4> p_l2_vtau2_vl2 = createMomentumFourVector( E_l2 + E_vtau2_vl2, 
    (1+xs[6]) * lep_px_l2 , (1+xs[6]) * lep_py_l2, (1+xs[6]) * lep_pz_l2 ); // combined particle of l2, v_tau2, v_l2

    array<double, 4> p_neut_l1 = createMomentumFourVector( sqrt(xs[0]*xs[0] + xs[1]*xs[1] + 
    xs[2]*xs[2]) , xs[0], xs[1], xs[2] );  // electron/muon neutrino from prompt W decay (v_l1) 
    array<double, 4> p_neut_tau1 = createMomentumFourVector( sqrt(xs[3]*xs[3] + xs[4]*xs[4] + 
    xs[5]*xs[5]) , xs[3], xs[4], xs[5] ); //  tau neutrino from prompt W decay (v_tau1)

    array<double, 4> p_l1 = createMomentumFourVector( sqrt(lep_px_l1*lep_px_l1 + lep_py_l1*lep_py_l1 + 
    lep_pz_l1*lep_pz_l1 + lep_mass_l1*lep_mass_l1) , lep_px_l1, lep_py_l1, lep_pz_l1 ); // Electron/Muon from prompt W decay
    array<double, 4> p_l2 = createMomentumFourVector( sqrt(lep_px_l2*lep_px_l2 + lep_py_l2*lep_py_l2 + 
    lep_pz_l2*lep_pz_l2 + lep_mass_l2*lep_mass_l2) , lep_px_l2, lep_py_l2, lep_pz_l2); // Electron/Muon from secondary W decay 
  
    array<double, 4> p_bjet_1 = createMomentumFourVector( bjet_E_1, bjet_px_1, bjet_py_1, bjet_pz_1); // Bjet that accompanies only prompt W decay
    array<double, 4> p_bjet_2 = createMomentumFourVector( bjet_E_2, bjet_px_2, bjet_py_2, bjet_pz_2); // Bjet that accompanies secondary W decay  
  // End of Creation of four vectors for final products

  // Calculate four vectors for intermediate products in decay ( w1_1, w1_2, t_1, t_2)
  array<double, 4> p_w1_1 = addFourVectors(p_l1, p_neut_l1);             // primary W boson that does not decay further
  array<double, 4> p_w1_2 = addFourVectors(p_neut_tau1, p_l2_vtau2_vl2); // primary W boson that decays further to secondary W boson
  array<double, 4> p_t_1 = addFourVectors(p_l1, p_neut_l1, p_bjet_1);    //  top quark that only decays to primay W 
  array<double, 4> p_t_2 = addFourVectors(p_bjet_2, p_neut_tau1, p_l2_vtau2_vl2); // top quark that decays to primary and secondary W bosons
  // End of Caculate four vectors for intermediate products

  // Calculate invariant masses of intermediate products (tau, w1_1, w1_2, w2, t_1, t_2)
  double m_w1_1 = getInvariantMass(p_w1_1);
  double m_w1_2 = getInvariantMass(p_w1_2);
  double m_t_1 = getInvariantMass(p_t_1);
  double m_t_2 = getInvariantMass(p_t_2);
  // End of Cacluate invariant masses of intermediate products

  // Caculate discriminant delta T3
  double delta_T3 = sqrt(pow(m_w1_1-80379, 2) + pow(m_w1_2-80379, 2) + pow(m_t_1-173000, 2) + 
       pow(m_t_2-173000, 2) + pow(xs[0] + xs[3] + (xs[6]*lep_px_l2) - et_miss_x, 2) + 
       pow(xs[1] + xs[4] + (xs[6]*lep_py_l2) - et_miss_y, 2));


   return std::make_tuple(delta_T3);
}


void mini::Loop()
{

//   In a ROOT session, you can do:
//      root> .L mini.C
//      root> mini t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   if (fChain == 0) return;

   // Get Total No. of Events Recorded
   Long64_t nentries = fChain->GetEntriesFast();

   // =========Create Histograms =========
   // Syntsetpax: TH1F("histogram_name", "Title", number_of_bins, x_min, x_max);
   TH1F *h_delta_T3 = new TH1F("h_delta_T3", "delta_T3 discriminant | events with at least 1 electron 1 muon 2 bjets ", 1000, 0,100000);
	
   // ========Loop Over All Events=========
   for (Long64_t jentry=0; jentry<1000;jentry++)
   {

    // Read the Information of the No. jentry Event
   	GetEntry(jentry);
    // Count Electrons, Muons, BJets in Event
    int muCount = 0; //  number of muons in event
    int eCount = 0; //  number of electrons in event
    int bjetCount = 0; //  number of bjets in event

    for (size_t lep_i=0; lep_i<lep_n; lep_i++) // loop over leptons
    {
    	if (lep_type->at(lep_i) == 11)  eCount++;  
    	if (lep_type->at(lep_i) == 13) muCount++;
    }

    for(size_t jet_i=0; jet_i<jet_trueflav->size(); jet_i++) // loop over jets
    {
    	if (jet_trueflav->at(jet_i) ==5) bjetCount++; 
    }

      // Only select events with more than 1 electron, 1 muon, 2 bjets
      // Select higher pt electron and muon
		if (eCount>=1 && muCount >=1 && bjetCount >= 2)
    {
      double e_pt[10] = {}; double mu_pt[10] = {}; double bjet_pt[20] = {};

      for (size_t lep_i=0; lep_i<lep_n; lep_i++) // loop over leptons
      {
        if (lep_type->at(lep_i) == 11)   // record kinematic info of electrons into array
        {
          e_pt[lep_i] = lep_pt->at(lep_i); 
        }
        if (lep_type->at(lep_i) == 13)  // record kinematic info of muons into array
        {
          mu_pt[lep_i] =  lep_pt->at(lep_i);
        }
      } // end of loop over leptons

      for(size_t jet_i=0; jet_i<jet_trueflav->size(); jet_i++) // loop over jets
      {
        if (jet_trueflav->at(jet_i) ==5) 
          {
            bjet_pt[jet_i] = jet_pt->at(jet_i); // record kinematic info of bjets into array
          } 
      }// end of loop over jets


      // Select highest pt electron, muon and b jets
      // Find index of largest pt value in vectors e_pt, mu_pt
      // Note: e_index is same as highest pt lep_i electron
      int e_index = std::distance(std::begin(e_pt), std::max_element(std::begin(e_pt), std::end(e_pt))); 
      int mu_index = std::distance(std::begin(mu_pt), std::max_element(std::begin(mu_pt), std::end(mu_pt)));
      int bjet_index_1 = std::distance(std::begin(bjet_pt), std::max_element(std::begin(bjet_pt), std::end(bjet_pt)));

      // find value of 2nd largest pt in jet_pt array
      int bjet_index_2 = 0;  
      double bjet_pt_2 = 0;

      for(int i=0; i<20; i++)
      {
        if(bjet_pt[i] == jet_pt->at(bjet_index_1)) {continue;} // do not consider largest pt in array
        if(bjet_pt[i] > bjet_pt_2) 
        {
          bjet_pt_2 = bjet_pt[i];
          bjet_index_2 = i;
        }
      }
      // End of selecting higehest pt 

	   // Additional requirement that event has to have e_pt,mu_pt,b_pt> 25Gev and opposite charged e/u 
    if( (lep_pt->at(e_index)>=25000 && lep_pt->at(mu_index)>=25000 && jet_pt->at(bjet_index_1)>=25000
      && jet_pt->at(bjet_index_2) >=25000 && lep_charge->at(e_index)==1 
      && lep_charge->at(mu_index)==-1) || (lep_pt->at(e_index)>=25000 && lep_pt->at(mu_index)>=25000
      && jet_pt->at(bjet_index_1)>=25000 && jet_pt->at(bjet_index_2) >=25000 
      && lep_charge->at(e_index)==-1 && lep_charge->at(mu_index)==1) )
	   {
        /* Minimize T3Function (possibly many times) */
        double delta_T3; 

        const int repeat = 100; // number of times to repeat the minimization
        const int outputs = 1;  // number of outputs from MinimizaTfunction()

        // Initialise 2D Matrix that stores output from each minimization attempt
        // E.g. { {delta_T3},   <- Minimization Attempt 1
        //        {delta_T3},   <- Minimization Attempt 2
        //        ...                            }
        double minimized_combo1[repeat][outputs];  // Combination 1 (refer to sketch on lab book)
        double minimized_combo2[repeat][outputs];  // Combination 2 
        double minimized_combo3[repeat][outputs];  // ...
        double minimized_combo4[repeat][outputs];  // ...
        double minimized_combo5[repeat][outputs];  //
        double minimized_combo6[repeat][outputs];  //
        double minimized_combo7[repeat][outputs];  //
        double minimized_combo8[repeat][outputs];  // 

        // Repeat minimization n times
        for (int i=0; i<repeat; i++)
        {
          // If final products are positron and muon
          if (lep_charge->at(e_index)==1 && lep_charge->at(mu_index)==-1)
          {
            // Do the minimization for combo 1
            tie(minimized_combo1[i][0]) = MinimizeT3function( 
            lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), 105.66, 
            lep_pt->at(e_index), lep_phi->at(e_index), lep_eta->at(e_index), 0.511, 
            met_et,  met_phi, 
            jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2), jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2), 
            jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1), jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1) 
            );

            // Do the minimization for combo 2
            tie(minimized_combo2[i][0]) = MinimizeT3function( 
            lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), 105.66, 
            lep_pt->at(e_index), lep_phi->at(e_index), lep_eta->at(e_index), 0.511, 
            met_et,  met_phi, 
            jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1), jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
            jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2), jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2)
            );

            // Do the minimization for combo 3
            tie(minimized_combo3[i][0]) = MinimizeT3function( 
            lep_pt->at(e_index), lep_phi->at(e_index), lep_eta->at(e_index), 0.511, 
            lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), 105.66, 
            met_et,  met_phi, 
            jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1), jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
            jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2), jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2)
            );

            // Do the minimization for combo 4
            tie(minimized_combo4[i][0]) = MinimizeT3function( 
            lep_pt->at(e_index), lep_phi->at(e_index), lep_eta->at(e_index), 0.511, 
            lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), 105.66, 
            met_et,  met_phi, 
            jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2), jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2),
            jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1), jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1)
            );

          }

          // If final products are electron and antimuon
          if (lep_charge->at(e_index)==-1 && lep_charge->at(mu_index)==1)
          {
            // Do the minimization for combo 5
            tie(minimized_combo5[i][0]) = MinimizeT3function( 
            lep_pt->at(e_index), lep_phi->at(e_index), lep_eta->at(e_index), 0.511, 
            lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), 105.66, 
            met_et,  met_phi, 
            jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2), jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2), 
            jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1), jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1) 
            );

            // Do the minimization for combo 6
            tie(minimized_combo6[i][0]) = MinimizeT3function( 
            lep_pt->at(e_index), lep_phi->at(e_index), lep_eta->at(e_index), 0.511, 
            lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), 105.66, 
            met_et,  met_phi, 
            jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1), jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
            jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2), jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2) 
            );

            // Do the minimization for combo 7
            tie(minimized_combo7[i][0]) = MinimizeT3function( 
            lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), 105.66, 
            lep_pt->at(e_index), lep_phi->at(e_index), lep_eta->at(e_index), 0.511, 
            met_et,  met_phi, 
            jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1), jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1),
            jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2), jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2)
            );

            // Do the minimization for combo 8
            tie(minimized_combo8[i][0]) = MinimizeT3function( 
            lep_pt->at(mu_index), lep_phi->at(mu_index), lep_eta->at(mu_index), 105.66, 
            lep_pt->at(e_index), lep_phi->at(e_index), lep_eta->at(e_index), 0.511, 
            met_et,  met_phi, 
            jet_E->at(bjet_index_2), jet_pt->at(bjet_index_2), jet_phi->at(bjet_index_2), jet_eta->at(bjet_index_2),
            jet_E->at(bjet_index_1), jet_pt->at(bjet_index_1), jet_phi->at(bjet_index_1), jet_eta->at(bjet_index_1)
            );
          }


        } /* End of Minimizing T3Function n times */

        // Find smallest value of delta T (a stupid but straightforward way)
        double delta_T3_combo1[repeat]; double delta_T3_combo2[repeat]; // create delta_T arrays
        double delta_T3_combo3[repeat]; double delta_T3_combo4[repeat];
        double delta_T3_combo5[repeat]; double delta_T3_combo6[repeat];
        double delta_T3_combo7[repeat]; double delta_T3_combo8[repeat];
        for(int i=0; i<repeat; i++)  // fill delta_T arrays
        {
          delta_T3_combo1[i] = minimized_combo1[i][0]; delta_T3_combo2[i] = minimized_combo2[i][0];
          delta_T3_combo3[i] = minimized_combo3[i][0]; delta_T3_combo4[i] = minimized_combo4[i][0];
          delta_T3_combo5[i] = minimized_combo5[i][0]; delta_T3_combo6[i] = minimized_combo6[i][0];
          delta_T3_combo7[i] = minimized_combo7[i][0]; delta_T3_combo8[i] = minimized_combo8[i][0];
        }
        // Find index of smallest delta_T in delta_T arrays
        int delta_T3_combo1_index = std::distance(std::begin(delta_T3_combo1), std::max_element(std::begin(delta_T3_combo1), std::end(delta_T3_combo1))); 
        int delta_T3_combo2_index = std::distance(std::begin(delta_T3_combo2), std::max_element(std::begin(delta_T3_combo2), std::end(delta_T3_combo2))); 
        int delta_T3_combo3_index = std::distance(std::begin(delta_T3_combo3), std::max_element(std::begin(delta_T3_combo3), std::end(delta_T3_combo3))); 
        int delta_T3_combo4_index = std::distance(std::begin(delta_T3_combo4), std::max_element(std::begin(delta_T3_combo4), std::end(delta_T3_combo4))); 
        int delta_T3_combo5_index = std::distance(std::begin(delta_T3_combo5), std::max_element(std::begin(delta_T3_combo5), std::end(delta_T3_combo5))); 
        int delta_T3_combo6_index = std::distance(std::begin(delta_T3_combo6), std::max_element(std::begin(delta_T3_combo6), std::end(delta_T3_combo6))); 
        int delta_T3_combo7_index = std::distance(std::begin(delta_T3_combo7), std::max_element(std::begin(delta_T3_combo7), std::end(delta_T3_combo7))); 
        int delta_T3_combo8_index = std::distance(std::begin(delta_T3_combo8), std::max_element(std::begin(delta_T3_combo8), std::end(delta_T3_combo8))); 

        // Create array that contain the smallest delta T from each comination's n minimizations
        double delta_T3_final[8] = {delta_T3_combo1[delta_T3_combo1_index],
          delta_T3_combo2[delta_T3_combo2_index], delta_T3_combo3[delta_T3_combo3_index],
          delta_T3_combo4[delta_T3_combo4_index],delta_T3_combo5[delta_T3_combo5_index],
          delta_T3_combo6[delta_T3_combo6_index],delta_T3_combo7[delta_T3_combo7_index],
          delta_T3_combo8[delta_T3_combo8_index]};
        // Find index of smallest delta_T in delta_T_final array
        int delta_T3_final_index = std::distance(std::begin(delta_T3_final), std::max_element(std::begin(delta_T3_final), std::end(delta_T3_final))); 
        // Assign minimized values with this smallest delta_T found from all minimizations of the 8 combinations
        if (delta_T3_final_index == 0)
        {
          delta_T3 = minimized_combo1[delta_T3_combo1_index][0];
        }
        if (delta_T3_final_index == 1)
        {
          delta_T3 = minimized_combo2[delta_T3_combo2_index][0];
        }
        if (delta_T3_final_index == 2)
        {
          delta_T3 = minimized_combo3[delta_T3_combo3_index][0];
        }
        if (delta_T3_final_index == 3)
        {
          delta_T3 = minimized_combo4[delta_T3_combo4_index][0];
        }
        if (delta_T3_final_index == 4)
        {
          delta_T3 = minimized_combo5[delta_T3_combo5_index][0];
        }
        if (delta_T3_final_index == 5)
        {
          delta_T3 = minimized_combo6[delta_T3_combo6_index][0];
        }
        if (delta_T3_final_index == 6)
        {
          delta_T3 = minimized_combo7[delta_T3_combo7_index][0];
        }
        if (delta_T3_final_index == 7)
        {
          delta_T3 = minimized_combo8[delta_T3_combo8_index][0];
        }

	      // Fill histograms
        h_delta_T3 -> Fill(delta_T3);

      } // End of requirement pt>25GeV
	  } // End of requirement at least 1e 1u 2b in event
  } // End of Loop Over All Events



   // ======= Draw Histograms ===========

	h_delta_T3 -> Draw();




   

   


}
