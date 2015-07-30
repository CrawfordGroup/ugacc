#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>

#include "hamiltonian.h"
#include "ccrhwavefunction.h"
// #include "cclhwavefunction.h"
// #include "perturbation.h"
#include "globals.h"

INIT_PLUGIN

using namespace boost;

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
PsiReturnType ugacc(Options& options)
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

  boost::shared_ptr<Wavefunction> ref = Process::environment.wavefunction();

  // Error trapping – need closed-shell SCF in place
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");
  if(options.get_str("REFERENCE") != "RHF")
    throw PSIEXCEPTION("Only for use with RHF references determinants.");
  for(int h=0; h < ref->nirrep(); h++)
    if(ref->soccpi()[h]) throw PSIEXCEPTION("UGACC is for closed-shell systems only.");

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);
  boost::shared_ptr<Hamiltonian> H(new Hamiltonian(psio, ref, spaces));
  boost::shared_ptr<CCRHWavefunction> CC(new CCRHWavefunction(ref, H, options, psio));

  double ecc = CC->compute_energy();

  if(!CC->dertype()) return Success;

//  ccwfn->hbar();
//  ccwfn->compute_lambda();

//  double eone = ccwfn->onepdm();
//  double etwo = ccwfn->twopdm();
//  double eref = ccwfn->reference_energy();

//  outfile->Printf("\tOne-Electron Energy        = %20.14f\n", eone);
//  outfile->Printf("\tTwo-Electron Energy        = %20.14f\n", etwo);
  if(CC->wfn() == "CCSD") {
//    outfile->Printf("\tCCSD Correlation Energy    = %20.14f (from density)\n", eone+etwo);
//    outfile->Printf("\tCCSD Total Energy          = %20.14f (from density)\n", eone+etwo+eref);
    outfile->Printf("\tCCSD Total Energy          = %20.14f (from ccwfn)\n", ecc);
  }
  else if(CC->wfn() == "CCSD_T") {
//    outfile->Printf("\tCCSD(T) Correlation Energy = %20.14f (from density)\n", eone+etwo);
//    outfile->Printf("\tCCSD(T) Total Energy       = %20.14f (from density)\n", eone+etwo+eref);
    outfile->Printf("\tCCSD(T) Total Energy       = %20.14f (from ccwfn)\n", ecc);
  }

  // Prepare property integrals for perturbed wave functions
//  boost::shared_ptr<Perturbation> mu(new Perturbation("Mu", ref));

  // Create similarity transformed property integrals
//  boost::shared_ptr<Pertbar> mubar(new PertBar(mu, ccwfn));

  // Solve perturbed wave function equations for give perturbation and field frequency
//  boost::shared_ptr<PertCC> mucc(new PertCC(mu, ccwfn));



  // Use perturbed wfns to construct the linear response function
  // Alternatively, build density-based linear response function

  return Success;
}

}} // End namespaces

