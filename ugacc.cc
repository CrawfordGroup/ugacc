#include "ccwavefunction.h"

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
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
  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  boost::shared_ptr<Wavefunction> ref = Process::environment.wavefunction();
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");
  boost::shared_ptr<Hamiltonian> H(new Hamiltonian(ref));
  boost::shared_ptr<CCWavefunction> ccwfn(new CCWavefunction(ref, H, options, psio));

  outfile->Printf("\n\tThe Coupled-Cluster Iteration:\n");
  outfile->Printf(  "\t---------------------------------------------\n");
  outfile->Printf(  "\t Iter   Correlation Energy  T1 Norm    RMS   \n");
  outfile->Printf(  "\t---------------------------------------------\n");
  outfile->Printf(  "\t  %3d  %20.15f\n", 0, ccwfn->energy());

  return Success;
}

}} // End namespaces

