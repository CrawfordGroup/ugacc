#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <map>

#include "hamiltonian.h"
#include "ccwfn.h"
#include "hbar.h"
#include "cclambda.h"
#include "ccdensity.h"
#include "perturbation.h"
#include "ccpert.h"

#include "array.h"

INIT_PLUGIN

using namespace boost;
using namespace std;

namespace psi { namespace ugacc {

extern "C" 
int read_options(std::string name, Options& options)
{
  if(name == "UGACC" || options.read_globals()) {
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF");
    options.add_str("WFN", "CCSD");
    options.add_str("DERTYPE", "NONE");
    options.add_int("MAXITER", 100);
    options.add_bool("DIIS", true);
    options.add_double("R_CONVERGENCE", 1e-7);
    options.add_bool("OOC", false);
  }

  return true;
}

extern "C" 
SharedWavefunction ugacc(SharedWavefunction ref, Options& options)
{
  outfile->Printf("\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t*         UGA-CC         *\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\n");

  outfile->Printf("\tWave function  = %s\n", options.get_str("WFN").c_str());
  outfile->Printf("\tMaxiter        = %d\n", options.get_int("MAXITER"));
  outfile->Printf("\tConvergence    = %3.1e\n", options.get_double("R_CONVERGENCE"));
  outfile->Printf("\tDIIS           = %s\n", options.get_bool("DIIS") ? "Yes" : "No");
  outfile->Printf("\tOut-of-core    = %s\n", options.get_bool("OOC") ? "Yes" : "No");
  outfile->Printf("\tDertype        = %s\n", options.get_str("DERTYPE").c_str());

  // Error trapping – need closed-shell SCF in place
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");
  if(options.get_str("REFERENCE") != "RHF")
    throw PSIEXCEPTION("Only for use with RHF references.");
  for(int h=0; h < ref->nirrep(); h++)
    if(ref->soccpi()[h]) throw PSIEXCEPTION("UGACC is for closed-shell systems only.");

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);
  boost::shared_ptr<Hamiltonian> H(new Hamiltonian(psio, ref, spaces));
  boost::shared_ptr<CCWfn> cc(new CCWfn(ref, H, options));

  double ecc = cc->compute_energy();

  if(options.get_str("DERTYPE") == "NONE") return ref;
 
  boost::shared_ptr<HBAR> hbar(new HBAR(H, cc));
  boost::shared_ptr<CCLambda> cclambda(new CCLambda(cc, hbar));
  cclambda->compute_lambda();

  boost::shared_ptr<CCDensity> ccdensity(new CCDensity(cc, cclambda));

  double eone = ccdensity->onepdm();
  double etwo = ccdensity->twopdm();
  double eref = cc->reference_energy();

  outfile->Printf("\tOne-Electron Energy        = %20.14f\n", eone);
  outfile->Printf("\tTwo-Electron Energy        = %20.14f\n", etwo);
  std::string wfn = options.get_str("WFN") == "CCSD_T" ? "CCSD(T)" : options.get_str("WFN");
  outfile->Printf("\t%s Correlation Energy    = %20.14f (from density)\n", wfn.c_str(), eone+etwo);
  outfile->Printf("\t%s Correlation Energy    = %20.14f (from ccwfn)\n", wfn.c_str(), ecc);
  outfile->Printf("\t%s Total Energy          = %20.14f (from density)\n", wfn.c_str(), eone+etwo+eref);
  outfile->Printf("\t%s Total Energy          = %20.14f (from ccwfn)\n", wfn.c_str(), ecc + eref);

  Process::environment.globals["CURRENT ENERGY"] = ecc + eref;

  // Prepare property integrals for perturbed wave functions
  boost::shared_ptr<MintsHelper> mints(new MintsHelper(ref->basisset(), options, 0));
  boost::shared_ptr<Perturbation> mu(new Perturbation("Mu", ref, mints, false));

  // Solve perturbed wave function equations for given perturbation and +/- field frequency
  map<string, boost::shared_ptr<CCPert> > cc_perts; 
  double omega = 0.00;
  vector<string> cart(3); cart[0] = "X"; cart[1] = "Y"; cart[2] = "Z";

  for(vector<string>::size_type iter = 0; iter != cart.size(); iter++) {
    string entry = "Mu" + cart[iter] + std::to_string(omega);
    outfile->Printf("\n\tCC Perturbed Wavefunction: %s\n", entry.c_str());
    cc_perts[entry] = boost::shared_ptr<CCPert>(new CCPert(mu->prop_p((int) iter), omega, cc, hbar));
    cc_perts[entry]->solve(right);
    if(omega != 0.0) {
      entry = "Mu" + cart[iter] + std::to_string(-omega);
      outfile->Printf("\n\tCC Perturbed Wavefunction: %s\n", entry.c_str());
      cc_perts[entry] = boost::shared_ptr<CCPert>(new CCPert(mu->prop_p((int) iter), -omega, cc, hbar));
      cc_perts[entry]->solve(right);
    }
  }

  // Use perturbed wfns to construct the linear response function
  // Alternatively, build density-based linear response function
  std::string mu1 = "Mu" + cart[2] + std::to_string(omega);
  std::string mu2 = "Mu" + cart[2] + std::to_string(omega);
//  boost::shared_ptr<CCLinResp> ccpolar(new CCLinResp(cc_perts[mu1], cc_perts[mu2]));

  return cc;
}

}} // End namespaces

