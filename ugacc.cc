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

void amp_write(int, double **, double ****, int, int, std::string);

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

  double rms = 0.0;
  for(int iter=1; iter <= ccwfn->maxiter(); iter++) {
    ccwfn->build_F();
    ccwfn->build_W();
    ccwfn->amp_save();
    ccwfn->build_t1();
    ccwfn->build_t2();
    rms = ccwfn->increment_amps();
    ccwfn->build_tau();

    outfile->Printf(  "\t  %3d  %20.15f  %8.6f  %8.6e\n",iter, ccwfn->energy(), ccwfn->t1norm(), rms);
    if(rms < ccwfn->convergence()) break;
//    if(ccwfn->do_diis()) ccwfn->diis(iter);
  }

  ccwfn->build_tau();

  if(rms >= ccwfn->convergence())
    throw PSIEXCEPTION("Computation has not converged.");

  amp_write(20, ccwfn->t1_p(), ccwfn->t2_p(), ccwfn->no(), ccwfn->nv(), "T"); 
  outfile->Printf("\n");

  return Success;
}

}} // End namespaces

