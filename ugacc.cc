#include "ccwavefunction.h"

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
#include "libparallel/ParallelPrinter.h"

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
void tstar_build(double **, double ****);
void F_build(void);
void W_build(void);
void t1_build(void);
void t2_build(void);
double increment_amps(double **, double **, double ****, double ****);
double t1norm(void);
void diis(int error_file, int amp_file, int iter, double **t1,
          double **t1old, double ****t2, double ****t2old);

double tcorr(void);
double tcorr_ooc(void);
double tcorr_ooc_TJL(void);

void init_L_amps(void);
void hbar(void);
void G_build(int);
double pseudoenergy(void);
void l1_build(void);
void l2_build(void);
void make_Z_amps(double **l1, double ****l2);
void tgrad(void);

void tgrad(void);
void tgrad_ooc(void);
void init_density(void);
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
  boost::shared_ptr<CCWavefunction> ccwfn(new CCWavefunction(ref, options, psio));

  boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  get_moinfo(ref, chkpt);
  integrals();
  denom();

  return Success;

  // ****** T-amplitude equations

  init_T_amps();

  outfile->Printf("\n\tThe Coupled-Cluster Iteration:\n");
  outfile->Printf(  "\t---------------------------------------------\n");
  outfile->Printf(  "\t Iter   Correlation Energy  T1 Norm    RMS   \n");
  outfile->Printf(  "\t---------------------------------------------\n");
  outfile->Printf(  "\t  %3d  %20.15f\n", 0,moinfo.eccsd = energy());

  double rms = 0.0;
  for(unsigned int iter=1; iter <= ccwfn->maxiter(); iter++) {
    amp_save(&moinfo.t1, &moinfo.t1old, &moinfo.t2, &moinfo.t2old);
    tau_build(iter, moinfo.t1old, moinfo.t2old);
    F_build(); 
    W_build();  
    t1_build();
    t2_build();
    rms = increment_amps(moinfo.t1, moinfo.t1old, moinfo.t2, moinfo.t2old);

    outfile->Printf(  "\t  %3d  %20.15f  %5.3f  %5.3e\n",iter, moinfo.eccsd = energy(), t1norm(), rms);
    if(rms < ccwfn->convergence()) break;
    if(ccwfn->do_diis()) diis(iter, 90, 91, moinfo.t1, moinfo.t1old, moinfo.t2, moinfo.t2old);
  }

  tau_build(2, moinfo.t1, moinfo.t2);
  tstar_build(moinfo.t1old, moinfo.t2old);

  if(rms >= ccwfn->convergence())
    throw PSIEXCEPTION("Computation has not converged.");

  outfile->Printf(  "\n\tCCSD Energy    = %20.14f\n",moinfo.eccsd+moinfo.escf);
  if(ccwfn->wfn() == "CCSD_T") {
    if(ccwfn->ooc()) {
      outfile->Printf("\t(T) Correction = %20.14f (occ)\n", moinfo.e_t = tcorr_ooc());
      outfile->Printf("\t(T) Correction = %20.14f (TJL)\n", tcorr_ooc_TJL());
    }
    else outfile->Printf("\t(T) Correction = %20.14f\n",moinfo.e_t = tcorr());
    outfile->Printf("\tCCSD(T) Energy = %20.14f\n",moinfo.escf+moinfo.eccsd+moinfo.e_t);
  }

  amp_write(20, moinfo.t1, moinfo.t2, "T"); outfile->Printf("\n");

  if(options.get_str("DERTYPE") == "NONE") {
    ccdump();
    cleanup();
    return Success;
  }

  // ****** Lambda-amplitude equations

  hbar();
  init_L_amps();
  init_density();

  if(ccwfn->wfn() == "CCSD_T") {
    if(ccwfn->ooc()) tgrad_ooc();
    else tgrad();
  }

  outfile->Printf("\n\tThe Coupled-Cluster Lambda Iteration:\n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t Iter   Correlation Energy  RMS   \n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t  %3d  %20.15f\n", 0, pseudoenergy());

  rms = 0.0;
  for(unsigned int iter=1; iter <= ccwfn->maxiter(); iter++) {
    amp_save(&moinfo.l1, &moinfo.l1old, &moinfo.l2, &moinfo.l2old);
    G_build(iter);
    l1_build();
    l2_build();
    rms = increment_amps(moinfo.l1, moinfo.l1old, moinfo.l2, moinfo.l2old);

    outfile->Printf(  "\t  %3d  %20.15f  %5.3e\n",iter, pseudoenergy(), rms);
    if(rms < ccwfn->convergence()) break;
    if(ccwfn->do_diis()) diis(iter, 92, 93, moinfo.l1, moinfo.l1old, moinfo.l2, moinfo.l2old);
  }
  if(rms >= ccwfn->convergence())
    throw PSIEXCEPTION("Computation has not converged.");

  amp_write(20, moinfo.l1, moinfo.l2, "L"); outfile->Printf("\n");

  // Also print non-UGA version of lambda amps for comparison to PSI4 UHF-CCSD(T) code
  make_Z_amps(moinfo.l1, moinfo.l2);

  double Eone = onepdm();
  double Etwo = twopdm();
  outfile->Printf("\tOne-electron energy        = %20.14f\n", Eone);
  outfile->Printf("\tTwo-electron energy        = %20.14f\n", Etwo);
  if(ccwfn->wfn() == "CCSD") {
    outfile->Printf("\tCCSD correlation energy    = %20.14f (from density)\n", Eone+Etwo);
    outfile->Printf("\tCCSD correlation energy    = %20.14f (from moinfo)\n", moinfo.eccsd);
    outfile->Printf("\tCCSD total energy          = %20.14f (from density)\n", Eone+Etwo+moinfo.escf);
  }
  else if(ccwfn->wfn() == "CCSD_T") {
    outfile->Printf("\tCCSD(T) correlation energy = %20.14f (from density)\n", Eone+Etwo);
    outfile->Printf("\tCCSD(T) correlation energy = %20.14f (from moinfo)\n", moinfo.eccsd+moinfo.e_t);
    outfile->Printf("\tCCSD(T) total energy       = %20.14f (from density)\n", Eone+Etwo+moinfo.escf);
  }

  if(chkpt->rd_nirreps() == 1 && moinfo.nact == moinfo.nmo) dipole(chkpt);

  ccdump();
  cleanup();

  return Success;
}

void title(void)
{
  outfile->Printf("\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t*         UGA-CC         *\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\n");
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
  amp_write(20, Z1, Z2, "Z"); outfile->Printf("\n");
  free_block(Z1);
  free_4d_array(Z2, no, no, nv);
}

}} // End namespaces

