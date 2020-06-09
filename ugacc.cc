/*
 * @BEGIN LICENSE
 *
 * ugacc by T. Daniel Crawford, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"

#include "psi4/libmints/mintshelper.h"
#include "psi4/libtrans/integraltransform.h"
#include <map>

#include "hamiltonian.h"
#include "ccwfn.h"
#include "hbar.h"
#include "cclambda.h"
#include "ccdensity.h"
#include "perturbation.h"
#include "ccpert.h"
#include "ccresp.h"

#include "array.h"

using namespace std;

namespace psi { namespace ugacc {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
  if(name == "UGACC" || options.read_globals()) {
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF");
    options.add_str("WFN", "CCSD");
    options.add_str("DERTYPE", "NONE");
    options.add_str("MYHAND", "RIGHT");
    options.add_int("MAXITER", 100);
    options.add_bool("DIIS", true);
    options.add_double("R_CONVERGENCE", 1e-7);
    options.add_double("MY_OMEGA", 0.00);
    options.add_bool("OOC", false);
    options.add_str("PROPERTY", "POLARIZABILITY", "POLARIZABILITY ROTATION ROA ROA_TENSOR ALL");

  }

  return true;
}

extern "C" PSI_API
SharedWavefunction ugacc(SharedWavefunction ref, Options& options)
{
  outfile->Printf("\t**************************\n");
  outfile->Printf("\t*                        *\n");
  outfile->Printf("\t*         UGA-CC         *\n");
  outfile->Printf("\t*                        *\n");
  outfile->Printf("\t**************************\n");
  outfile->Printf("\n");

  outfile->Printf("\tWave function  = %s\n", options.get_str("WFN").c_str());
  outfile->Printf("\tMaxiter        = %d\n", options.get_int("MAXITER"));
  outfile->Printf("\tConvergence    = %3.1e\n", options.get_double("R_CONVERGENCE"));
  outfile->Printf("\tDIIS           = %s\n", options.get_bool("DIIS") ? "Yes" : "No");
  outfile->Printf("\tOut-of-core    = %s\n", options.get_bool("OOC") ? "Yes" : "No");
  outfile->Printf("\tDertype        = %s\n", options.get_str("DERTYPE").c_str());
  outfile->Printf("\tOMEGA          = %3.1e\n", options.get_double("MY_OMEGA"));
  outfile->Printf("\tHAND           = %s\n", options.get_str("MYHAND").c_str()); 

  // Error trapping – need closed-shell SCF in place
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");
  if(options.get_str("REFERENCE") != "RHF")
    throw PSIEXCEPTION("Only for use with RHF references.");
  for(int h=0; h < ref->nirrep(); h++)
    if(ref->soccpi()[h]) throw PSIEXCEPTION("UGACC is for closed-shell systems only.");
  
  // Set up I/O object
  shared_ptr<PSIO> psio(_default_psio_lib_);

  // Prepare MO space vector that runs over all orbitals
  std::vector<shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);

  // Prepare Hamiltonian (transform the integrals and sort them into member arrays)
  shared_ptr<Hamiltonian> H(new Hamiltonian(psio, ref, spaces));

  // Prepare to build a CC wave function (preps denominators and storage for amps and DIIS vectors)
  shared_ptr<CCWfn> cc(new CCWfn(ref, H, options));

  // Solve the T-amplitude equations and compute the CC energy (including triples, if needed)
  double ecc = cc->compute_energy();

  if(options.get_str("DERTYPE") == "NONE") return ref;
 
  // Build the similarity-transformed Hamiltonian
  shared_ptr<HBAR> hbar(new HBAR(H, cc));

  // Solve for Lambda amplitude equations (Lagrange multipliers)
  shared_ptr<CCLambda> cclambda(new CCLambda(cc, hbar));
  cclambda->compute_lambda();

  // Build the one- and two-electron densities
  shared_ptr<CCDensity> ccdensity(new CCDensity(cc, cclambda));

  double eone = ccdensity->onepdm();
  double etwo = ccdensity->twopdm();
  double eref = ref->energy();

  outfile->Printf("\tOne-Electron Energy        = %20.14f\n", eone);
  outfile->Printf("\tTwo-Electron Energy        = %20.14f\n", etwo);
  std::string wfn = options.get_str("WFN") == "CCSD_T" ? "CCSD(T)" : options.get_str("WFN");
  outfile->Printf("\t%s Correlation Energy    = %20.14f (from density)\n", wfn.c_str(), eone+etwo);
  outfile->Printf("\t%s Correlation Energy    = %20.14f (from ccwfn)\n", wfn.c_str(), ecc);
  outfile->Printf("\t%s Total Energy          = %20.14f (from density)\n", wfn.c_str(), eone+etwo+eref);
  outfile->Printf("\t%s Total Energy          = %20.14f (from ccwfn)\n", wfn.c_str(), ecc + eref);

  Process::environment.globals["CURRENT ENERGY"] = ecc + eref;

  // Prepare property integrals for perturbed wave functions

  shared_ptr<MintsHelper> mints(new MintsHelper(ref->basisset(), options, 0));
  shared_ptr<Perturbation> mu(new Perturbation("Mu", ref, mints, false));
  shared_ptr<Perturbation> am(new Perturbation("L", ref, mints, false));


  // Solve perturbed wave function equations for given perturbation and +/- field frequency

  map<string, shared_ptr<CCPert> > cc_perts; 
  map<string, double > polars; 
  map<string, double > rots; 
  double polar;
  double rotation;
  double omega = options.get_double("MY_OMEGA");
  vector<string> cart(3); cart[0] = "X"; cart[1] = "Y"; cart[2] = "Z";

  // const char * rol = options.get_str("MYHAND").c_str() ;
  hand my_hand ;
  if (!strcmp(options.get_str("MYHAND").c_str(),"RIGHT")) 
     my_hand = right;
  else  
     my_hand = left;

  /* Below is the recipe for calculating length gauge optical rotation and polarizability*/

    outfile->Printf("\n\tSolving right hand perturbed CC amplitudes\n");
    for(vector<string>::size_type iter = 0; iter != cart.size(); iter++) {
    string entry = "Mu" + cart[iter] + std::to_string(omega);
    string entry_1 = "L" + cart[iter] + std::to_string(omega);
    outfile->Printf("\n\tCC RH Perturbed Wavefunction: %s\n", entry.c_str());
    cc_perts[entry] = shared_ptr<CCPert>(new CCPert(mu->prop_p((int) iter), omega, cc, hbar, cclambda));
    cc_perts[entry]->solve(right);
    outfile->Printf("\n\tCC RH Perturbed Wavefunction: %s\n", entry_1.c_str());
    cc_perts[entry_1] = shared_ptr<CCPert>(new CCPert(am->prop_p((int) iter), omega, cc, hbar, cclambda));
    cc_perts[entry_1]->solve(right);
    if(omega != 0.0 && my_hand == right) {
      entry = "Mu" + cart[iter] + std::to_string(-omega);
      outfile->Printf("\n\tCC RH Perturbed Wavefunction: %s\n", entry.c_str());
      cc_perts[entry] = shared_ptr<CCPert>(new CCPert(mu->prop_p((int) iter), -omega, cc, hbar, cclambda));
      cc_perts[entry]->solve(right);
     }

    if (my_hand == left){
       outfile->Printf("\n\tSolving left hand perturbed CC amplitudes\n");
       outfile->Printf("\n\tCC LH Perturbed Wavefunction: %s\n", entry.c_str());
       cc_perts[entry]->solve(left);
       outfile->Printf("\n\tCC LH Perturbed Wavefunction: %s\n", entry_1.c_str());
       cc_perts[entry_1]->solve(left);
    }
   }
 
  /* Dipole polarizabilities */ 

  for(vector<string>::size_type p = 0; p != cart.size(); p++){
    for(vector<string>::size_type q = 0 ; q <= p; q++) {
      string pert_p = "Mu" + cart[p] + std::to_string(omega);
      string pert_q = "Mu" + cart[q] + std::to_string(omega);
      shared_ptr<CCResp> ccpolar(new CCResp(cc_perts[pert_p], cc_perts[pert_q]));
      if (p == q)
         polar = ccpolar->linresp(cc_perts[pert_p], cc_perts[pert_q]);
      else {
         polar = 0.5 * ccpolar->linresp(cc_perts[pert_p], cc_perts[pert_q]);
         polar += 0.5 * ccpolar->linresp(cc_perts[pert_q], cc_perts[pert_p]);
      }     
      string label = "<<Mu_" + cart[p] + ";" "Mu_" + cart[q] + ">>";
      polars[label] = polar;
      label = "<<Mu_" + cart[q] + ";" "Mu_" + cart[p] + ">>";
      polars[label] = polar;
    }
  }

  /* Length gauge optical rotation */ 

  for(vector<string>::size_type p = 0; p != cart.size(); p++){
    for(vector<string>::size_type q = 0 ; q != cart.size(); q++) {
      string pert_p = "Mu" + cart[p] + std::to_string(omega);
      string pert_q = "L" + cart[q] + std::to_string(omega);
      shared_ptr<CCResp> ccrot(new CCResp(cc_perts[pert_p], cc_perts[pert_q]));
      rotation = 0.5 * ccrot->linresp(cc_perts[pert_p], cc_perts[pert_q]);
      rotation -= 0.5 * ccrot->linresp(cc_perts[pert_q], cc_perts[pert_p]);
      string label = "<<Mu_" + cart[p] + ";" "L_" + cart[q] + ">>";
      rots[label] = -1.0 * rotation;
    }
  }



    outfile->Printf("\n\tDipole Polarizabilities (au)\n\n");
   for(auto elem : polars){
      outfile->Printf("\t%s : %20.14lf ", elem.first.c_str(), elem.second);
      outfile->Printf("\n\n");
      }

    outfile->Printf("\n\tOptical rotation (length gauge) (a.u)\n\n");
   for(auto elem : rots){
      outfile->Printf("\t%s : %20.14lf ", elem.first.c_str(), elem.second);
      outfile->Printf("\n\n");
      }

  return cc;
}

}} // End namespaces

