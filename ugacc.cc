#include "ccwavefunction.h"

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
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
  boost::shared_ptr<CCWavefunction> ccwfn(new CCWavefunction(ref, H, options, psio));

  double ecc = ccwfn->compute_energy();

  if(!ccwfn->dertype()) return Success;

  ccwfn->hbar();
  ccwfn->compute_lambda();

  double Eone = ccwfn->onepdm();
  double Etwo = ccwfn->twopdm();
  double eref = ccwfn->reference_energy();

  outfile->Printf("\tOne-Electron Energy        = %20.14f\n", Eone);
  outfile->Printf("\tTwo-Electron Energy        = %20.14f\n", Etwo);
  if(ccwfn->wfn() == "CCSD") {
    outfile->Printf("\tCCSD Correlation Energy    = %20.14f (from density)\n", Eone+Etwo);
    outfile->Printf("\tCCSD Total Energy          = %20.14f (from density)\n", Eone+Etwo+eref);
    outfile->Printf("\tCCSD Total Energy          = %20.14f (from ccwfn)\n", ecc);
  }
  else if(ccwfn->wfn() == "CCSD_T") {
    outfile->Printf("\tCCSD(T) Correlation Energy = %20.14f (from density)\n", Eone+Etwo);
    outfile->Printf("\tCCSD(T) Total Energy       = %20.14f (from density)\n", Eone+Etwo+eref);
    outfile->Printf("\tCCSD(T) Total Energy       = %20.14f (from ccwfn)\n", ecc);
  }


  return Success;
}

}} // End namespaces

