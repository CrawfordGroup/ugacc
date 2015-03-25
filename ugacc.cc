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
  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  boost::shared_ptr<Wavefunction> ref = Process::environment.wavefunction();
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");

  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);

  boost::shared_ptr<Hamiltonian> H(new Hamiltonian(ref, spaces));

  return Success;

  boost::shared_ptr<CCWavefunction> ccwfn(new CCWavefunction(ref, H, options, psio));

  double eref, emp2, eccsd, et;
  eref = ccwfn->reference_energy();

  outfile->Printf("\n\tThe Coupled-Cluster Iteration:\n");
  outfile->Printf(  "\t---------------------------------------------------\n");
  outfile->Printf(  "\t Iter   Correlation Energy   T1 Norm      RMS   \n");
  outfile->Printf(  "\t---------------------------------------------------\n");
  outfile->Printf(  "\t  %3d  %20.15f\n", 0, emp2 = ccwfn->energy());


  double rms = 0.0;
  for(int iter=1; iter <= ccwfn->maxiter(); iter++) {
    ccwfn->amp_save("T");
    ccwfn->build_F();
    ccwfn->build_W();
    ccwfn->build_t1();
    ccwfn->build_t2();
    rms = ccwfn->increment_amps("T");
    if(rms < ccwfn->convergence()) break;
    if(ccwfn->do_diis()) ccwfn->diis(iter, "T");
    ccwfn->build_tau();
    outfile->Printf(  "\t  %3d  %20.15f  %8.6f  %8.6e\n",iter, eccsd = ccwfn->energy(), ccwfn->t1norm(), rms);
  }

  ccwfn->build_tau();
  ccwfn->build_tstar();

  if(rms >= ccwfn->convergence())
    throw PSIEXCEPTION("Computation has not converged.");

  outfile->Printf(  "\n\tMP2 Correlation Energy     = %20.14f\n", emp2);
  outfile->Printf(  "\tMP2 Total Energy           = %20.14f\n", emp2 + eref);
  outfile->Printf(  "\tCCSD Correlation Energy    = %20.14f\n", eccsd);
  outfile->Printf(  "\tCCSD Total Energy          = %20.14f\n", eccsd + eref);
  if(ccwfn->wfn() == "CCSD_T") {
    if(ccwfn->ooc()) {
      outfile->Printf("\t(T) Correction             = %20.14f (ooc)\n", et = ccwfn->tcorr_ooc());
      outfile->Printf("\t(T) Correction             = %20.14f (TJL)\n", ccwfn->tcorr_ooc_TJL());
    }
    else outfile->Printf("\t(T) Correction             = %20.14f\n", et = ccwfn->tcorr());
    outfile->Printf("\tCCSD(T) Correlation Energy = %20.14f\n", eccsd + et);
    outfile->Printf("\tCCSD(T) Total Energy       = %20.14f\n", eref + eccsd + et);
  }

  ccwfn->amp_write(20, "T"); 
  outfile->Printf("\n");

  if(!ccwfn->dertype()) return Success;

  ccwfn->hbar();
  ccwfn->init_lambda();
  ccwfn->init_density();

  if(ccwfn->wfn() == "CCSD_T") {
    if(ccwfn->ooc()) ccwfn->tgrad_ooc();
    else ccwfn->tgrad();
  }

  outfile->Printf("\n\tThe Coupled-Cluster Lambda Iteration:\n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t Iter   Correlation Energy  RMS   \n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t  %3d  %20.15f\n", 0, ccwfn->pseudoenergy());

  rms = 0.0;
  for(int iter=1; iter <= ccwfn->maxiter(); iter++) {
    ccwfn->amp_save("L");
    ccwfn->build_G();
    ccwfn->build_l1();
    ccwfn->build_l2();
    rms = ccwfn->increment_amps("L");
    if(rms < ccwfn->convergence()) break;
    if(ccwfn->do_diis()) ccwfn->diis(iter, "L");
    outfile->Printf(  "\t  %3d  %20.15f  %5.3e\n",iter, ccwfn->pseudoenergy(), rms);
  }
  if(rms >= ccwfn->convergence())
    throw PSIEXCEPTION("Computation has not converged.");

  ccwfn->amp_write(20, "L");
  outfile->Printf("\n");

  double Eone = ccwfn->onepdm();
  double Etwo = ccwfn->twopdm();
  outfile->Printf("\tOne-Electron Energy        = %20.14f\n", Eone);
  outfile->Printf("\tTwo-Electron Energy        = %20.14f\n", Etwo);
  if(ccwfn->wfn() == "CCSD") {
    outfile->Printf("\tCCSD Correlation Energy    = %20.14f (from density)\n", Eone+Etwo);
    outfile->Printf("\tCCSD Correlation Energy    = %20.14f (from ccwfn)\n", eccsd);
    outfile->Printf("\tCCSD Total Energy          = %20.14f (from density)\n", Eone+Etwo+eref);
    outfile->Printf("\tCCSD Total Energy          = %20.14f (from ccwfn)\n", eccsd+eref);
  }
  else if(ccwfn->wfn() == "CCSD_T") {
    outfile->Printf("\tCCSD(T) Correlation Energy = %20.14f (from density)\n", Eone+Etwo);
    outfile->Printf("\tCCSD(T) Correlation Energy = %20.14f (from ccwfn)\n", eccsd + et);
    outfile->Printf("\tCCSD(T) Total Energy       = %20.14f (from density)\n", Eone+Etwo+eref);
    outfile->Printf("\tCCSD(T) Total Energy       = %20.14f (from ccwfn)\n", eccsd+et+eref);
  }


  return Success;
}

}} // End namespaces

