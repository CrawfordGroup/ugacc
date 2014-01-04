#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

INIT_PLUGIN

using namespace boost;

namespace psi { namespace ugacc {

void title(void);
void get_moinfo(boost::shared_ptr<Wavefunction> wfn, boost::shared_ptr<Chkpt> chkpt);
void integrals(void);
void denom(void);
void init_T_amps(void);
double energy(void);
void amp_save(double **, double **, double ****, double ****);
void tau_build(int);
void F_build(void);
void W_build(void);
void t1_build(void);
void t2_build(void);
double increment_amps(double **, double **, double ****, double ****);
double t1norm(void);
void diis(int);

void init_L_amps(void);
void hbar(void);
void G_build(int);

void ccdump(void);
void amp_write(int, double **, double ****);
void cleanup(void);

extern "C" 
int read_options(std::string name, Options& options)
{
  if(name == "UGACC" || options.read_globals()) {
    /*- The amount of information printed to the output file -*/
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF");
    options.add_str("WFN", "CCSD");
    options.add_int("MAXITER", 100);
    options.add_bool("DIIS", true);
    options.add_double("CONVERGENCE", 1e-7);
  }

  return true;
}

extern "C" 
PsiReturnType ugacc(Options& options)
{
  title();
  params.ref = options.get_str("REFERENCE");
  params.wfn = options.get_str("WFN");
  params.convergence = options.get_double("CONVERGENCE");
  params.do_diis = options.get_bool("DIIS");
  params.maxiter = options.get_int("MAXITER");

  fprintf(outfile, "\tWave function  = %s\n", params.wfn.c_str());
  fprintf(outfile, "\tReference      = %s\n", params.ref.c_str());
  fprintf(outfile, "\tMaxiter        = %d\n", params.maxiter);
  fprintf(outfile, "\tConvergence    = %3.1e\n", params.convergence);
  fprintf(outfile, "\tDIIS           = %s\n", params.do_diis ? "Yes" : "No");
  fflush(outfile);

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
  if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");
  boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  get_moinfo(wfn, chkpt);
  integrals();
  denom();
  init_T_amps();

  fprintf(outfile, "\n\tThe Coupled-Cluster Iteration:\n");
  fprintf(outfile,   "\t---------------------------------------------\n");
  fprintf(outfile,   "\t Iter   Correlation Energy  T1 Norm    RMS   \n");
  fprintf(outfile,   "\t---------------------------------------------\n");
  fprintf(outfile,   "\t  %3d  %20.15f\n", 0,moinfo.eccsd = energy());
  fflush(outfile);

  double rms = 0.0;
  for(int iter=1; iter <= params.maxiter; iter++) {
    amp_save(moinfo.t1, moinfo.t1old, moinfo.t2, moinfo.t2old);
    tau_build(iter);
    F_build(); 
    W_build();  
    t1_build();
    t2_build();
    rms = increment_amps(moinfo.t1, moinfo.t1old, moinfo.t2, moinfo.t2old);

    fprintf(outfile,   "\t  %3d  %20.15f  %5.3f  %5.3e\n",iter, moinfo.eccsd = energy(), t1norm(), rms);
    fflush(outfile);
    if(rms < params.convergence) {
      fprintf(outfile, "\n\tComputation has converged!\n");
      fflush(outfile);
      break;
    }
    if(params.do_diis) diis(iter, 90, 91, moinfo.t1, moinfo.t1old,
                            moinfo.t2, moinfo.t2old);
  }

  if(rms >= params.convergence)
    throw PSIEXCEPTION("Computation has not converged.");

  amp_write(20, moinfo.t1, moinfo.t2); fprintf(outfile, "\n");
  fprintf(outfile,  "\tCCSD Energy    = %20.14f\n",moinfo.eccsd+moinfo.escf);

  // ****** Lambda equations

  hbar();
  init_L_amps();

  fprintf(outfile, "\n\tThe Coupled-Cluster Lambda Iteration:\n");
  fprintf(outfile,   "\t-------------------------------------\n");
  fprintf(outfile,   "\t Iter   Correlation Energy  RMS   \n");
  fprintf(outfile,   "\t-------------------------------------\n");
  fprintf(outfile,   "\t  %3d  %20.15f\n", 0, pseudoenergy());
  fflush(outfile);

  rms = 0.0;
  for(iter=1; iter <= params.maxiter; iter++) {
    amp_save(moinfo.l1, moinfo.l1old, moinfo.l2, moinfo.l2old);
    G_build(iter);
/*
    l1_build();
    l2_build();
*/
    rms = increment_amps(moinfo.l1, moinfo.l1old, moinfo.l2, moinfo.l2old);
    fprintf(outfile,   "\t  %3d  %20.15f  %5.3e\n",iter, pseudoenergy(), rms);
    fflush(outfile);
    if(rms < params.convergence) {
      fprintf(outfile, "\n\tComputation has converged!\n");
      fflush(outfile);
      break;
    }
    if(params.do_diis) diis(iter, 92, 93, moinfo.l1, moinfo.l1old,
                            moinfo.l2, moinfo.l2old);
  }
  if(rms >= params.convergence)
    throw PSIEXCEPTION("Computation has not converged.");

  amp_write(20, moinfo.l1, moinfo.l2); fprintf(outfile, "\n");

  ccdump();
  cleanup();

  return Success;
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*         UGA-CC         *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

}} // End namespaces

