#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>
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
void amp_save(double ***, double ***, double *****, double *****);
void tau_build(int, double **, double ****);
void F_build(void);
void W_build(void);
void t1_build(void);
void t2_build(void);
double increment_amps(double **, double **, double ****, double ****);
double t1norm(void);
void diis(int error_file, int amp_file, int iter, double **t1,
          double **t1old, double ****t2, double ****t2old);
double triples(void);

void init_L_amps(void);
void hbar(void);
void G_build(int);
double pseudoenergy(void);
void l1_build(void);
void l2_build(void);
void make_Z_amps(double **l1, double ****l2);

double onepdm(void);
double twopdm(void);
void dipole(boost::shared_ptr<Chkpt> chkpt);

void ccdump(void);
void amp_write(int, double **, double ****, std::string);
void cleanup(void);

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
  }

  return true;
}

extern "C" 
PsiReturnType ugacc(Options& options)
{
  title();
  params.ref = options.get_str("REFERENCE");
  params.wfn = options.get_str("WFN");
  if(options.get_str("DERTYPE") == "NONE") params.dertype = 0;
  else if(options.get_str("DERTYPE") == "FIRST") params.dertype = 1;
  params.convergence = options.get_double("R_CONVERGENCE");
  params.do_diis = options.get_bool("DIIS");
  params.maxiter = options.get_int("MAXITER");

  fprintf(outfile, "\tWave function  = %s\n", params.wfn.c_str());
  fprintf(outfile, "\tReference      = %s\n", params.ref.c_str());
  fprintf(outfile, "\tComputation    = %s\n", params.dertype ? "Gradient" : "Energy");
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

  // ****** T-amplitude equations

  init_T_amps();

  fprintf(outfile, "\n\tThe Coupled-Cluster Iteration:\n");
  fprintf(outfile,   "\t---------------------------------------------\n");
  fprintf(outfile,   "\t Iter   Correlation Energy  T1 Norm    RMS   \n");
  fprintf(outfile,   "\t---------------------------------------------\n");
  fprintf(outfile,   "\t  %3d  %20.15f\n", 0,moinfo.eccsd = energy());
  fflush(outfile);

  double rms = 0.0;
  for(int iter=1; iter <= params.maxiter; iter++) {
    amp_save(&moinfo.t1, &moinfo.t1old, &moinfo.t2, &moinfo.t2old);
    tau_build(iter, moinfo.t1old, moinfo.t2old);
    F_build(); 
    W_build();  
    t1_build();
    t2_build();
    rms = increment_amps(moinfo.t1, moinfo.t1old, moinfo.t2, moinfo.t2old);

    fprintf(outfile,   "\t  %3d  %20.15f  %5.3f  %5.3e\n",iter, moinfo.eccsd = energy(), t1norm(), rms);
    fflush(outfile);
    if(rms < params.convergence) break;
    if(params.do_diis) diis(iter, 90, 91, moinfo.t1, moinfo.t1old,
                            moinfo.t2, moinfo.t2old);
  }

  if(rms >= params.convergence)
    throw PSIEXCEPTION("Computation has not converged.");

  fprintf(outfile,   "\n\tCCSD Energy    = %20.14f\n",moinfo.eccsd+moinfo.escf);
  if(params.wfn == "CCSD_T") {
    fprintf(outfile, "\t(T) Correction = %20.14f\n",moinfo.e_t = triples());
    fprintf(outfile, "\tCCSD(T) Energy = %20.14f\n",moinfo.escf+moinfo.eccsd+moinfo.e_t);
  }

  tau_build(2, moinfo.t1, moinfo.t2);
  amp_write(20, moinfo.t1, moinfo.t2, "T"); fprintf(outfile, "\n");

  if(!params.dertype) return Success;

  // ****** Lambda-amplitude equations

  hbar();
  init_L_amps();

  fprintf(outfile, "\n\tThe Coupled-Cluster Lambda Iteration:\n");
  fprintf(outfile,   "\t-------------------------------------\n");
  fprintf(outfile,   "\t Iter   Correlation Energy  RMS   \n");
  fprintf(outfile,   "\t-------------------------------------\n");
  fprintf(outfile,   "\t  %3d  %20.15f\n", 0, pseudoenergy());
  fflush(outfile);

  rms = 0.0;
  for(int iter=1; iter <= params.maxiter; iter++) {
    amp_save(&moinfo.l1, &moinfo.l1old, &moinfo.l2, &moinfo.l2old);
    G_build(iter);
    l1_build();
    l2_build();
    rms = increment_amps(moinfo.l1, moinfo.l1old, moinfo.l2, moinfo.l2old);

    fprintf(outfile,   "\t  %3d  %20.15f  %5.3e\n",iter, pseudoenergy(), rms);
    fflush(outfile);
    if(rms < params.convergence) break;
    if(params.do_diis) diis(iter, 92, 93, moinfo.l1, moinfo.l1old,
                            moinfo.l2, moinfo.l2old);
  }
  if(rms >= params.convergence)
    throw PSIEXCEPTION("Computation has not converged.");

  amp_write(20, moinfo.l1, moinfo.l2, "L"); fprintf(outfile, "\n");

  // Also print non-UGA version of lambda amps for comparison to PSI4 UHF-CCSD(T) code
  make_Z_amps(moinfo.l1, moinfo.l2);

  double Eone = onepdm();
  double Etwo = twopdm();
  fprintf(outfile, "One-electron energy        = %20.14f\n", Eone);
  fprintf(outfile, "Two-electron energy        = %20.14f\n", Etwo);
  if(params.wfn == "CCSD")
    fprintf(outfile, "CCSD correlation energy    = %20.14f\n", Eone+Etwo);
  else if(params.wfn == "CCSD_T")
    fprintf(outfile, "CCSD(T) correlation energy = %20.14f\n", Eone+Etwo);

  if(chkpt->rd_nirreps() == 1 && moinfo.nact == moinfo.nmo) dipole(chkpt);

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

void make_Z_amps(double **l1, double ****l2)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **Z1 = block_matrix(no, nv);
  double ****Z2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      Z1[i][a] = 0.5 * l1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          Z2[i][j][a][b] = (1./3.)*l2[i][j][a][b] + (1./6.)*l2[i][j][b][a];
   }
  amp_write(20, Z1, Z2, "Z"); fprintf(outfile, "\n");
  free_block(Z1);
  free_4d_array(Z2, no, no, nv);
}

}} // End namespaces

