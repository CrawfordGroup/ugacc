#include "ccwavefunction.h"
#include "globals.h"
#include "hamiltonian.h"
#include <boost/shared_ptr.hpp>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>

namespace psi {

CCWavefunction::CCWavefunction(boost::shared_ptr<Wavefunction> reference, boost::shared_ptr<Hamiltonian> H, Options &options, boost::shared_ptr<PSIO> psio) : Wavefunction(options, psio)
{
  outfile->Printf("\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t*         UGA-CC         *\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\n");

  if(options.get_str("REFERENCE") != "RHF")
    throw PSIEXCEPTION("Only for use with RHF references determinants.");

  wfn_ = options.get_str("WFN");
  convergence_ = options.get_double("R_CONVERGENCE");
  maxiter_ = options.get_int("MAXITER");
  do_diis_ = options.get_bool("DIIS");
  ooc_ = options.get_bool("OOC");

  outfile->Printf("\tWave function  = %s\n", wfn().c_str());
  outfile->Printf("\tMaxiter        = %d\n", maxiter());
  outfile->Printf("\tConvergence    = %3.1e\n", convergence());
  outfile->Printf("\tDIIS           = %s\n", do_diis() ? "Yes" : "No");
  outfile->Printf("\tOut-of-core    = %s\n", ooc() ? "Yes" : "No");

  set_reference_wavefunction(reference);
  copy(reference);

  int nfrzv = 0;
  no_ = nv_ = 0;
  for(int i=0; i < nirrep_; i++) {
    no_ += doccpi_[i] - frzcpi_[i];
    nv_ += nmopi_[i] - doccpi_[i] - frzvpi_[i];
    nfrzv += frzvpi_[i];
  }
  char ** labels = molecule_->irrep_labels();

  outfile->Printf("\n\tReference Wfn Parameters:\n");
  outfile->Printf("\t---------------------------\n");
  outfile->Printf("\tNumber of irreps        = %d\n", nirrep_);
  outfile->Printf("\tNumber of MOs           = %d\n", nmo_);
  outfile->Printf("\tNumber of active MOs    = %d\n", no_+nv_);
  outfile->Printf("\tNumber of active occ    = %d\n", no_);
  outfile->Printf("\tNumber of active vir    = %d\n", nv_);
  outfile->Printf("\tNumber of frozen occ    = %d\n", nfrzc_);
  outfile->Printf("\tNumber of frozen vir    = %d\n\n", nfrzv);
  outfile->Printf("\tLabel\t# MOs\t# FZDC\t# DOCC\t# VIRT\t# FZVR\n");
  outfile->Printf("\t-----\t-----\t------\t------\t------\t------\n");
  for(int i=0; i < nirrep_; i++) {
      outfile->Printf("\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\n",
              labels[i],nmopi_[i],frzcpi_[i],doccpi_[i],nmopi_[i]-doccpi_[i],frzvpi_[i]);
    }
  outfile->Printf("\n\tNuclear Repulsion Energy    = %20.15f\n", molecule_->nuclear_repulsion_energy());
  outfile->Printf( "\tFrozen Core Energy          = %20.15f\n", efzc_);
  outfile->Printf( "\tTotal SCF Energy (chkpt)    = %20.15f\n", reference_wavefunction_->reference_energy());

  for(int i=0; i < nirrep_; i++) free(labels[i]);
  free(labels);

  H_ = H; // does this copy properly?

  // Prepare energy denominators
  int no = no_;
  int nv = nv_;
  double **fock = H_->fock_p();
  double ****ints = H_->ints_p();

  D1_ = block_matrix(no,nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      D1_[i][a] = fock[i][i] - fock[a+no][a+no];

  D2_ = init_4d_array(no,no,nv,nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          D2_[i][j][a][b] = fock[i][i] + fock[j][j] - fock[a+no][a+no] - fock[b+no][b+no];

  t1_ = block_matrix(no,nv);
  t1old_ = block_matrix(no,nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      t1_[i][a] = fock[i][a+no]/D1_[i][a];

  t2_ = init_4d_array(no,no,nv,nv);
  t2old_ = init_4d_array(no,no,nv,nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          t2_[i][j][a][b] = ints[i][j][a+no][b+no]/D2_[i][j][a][b];

  tau_ = init_4d_array(no,no,nv,nv);
  ttau_ = init_4d_array(no,no,nv,nv);

  build_tau();
}

CCWavefunction::~CCWavefunction()
{
  free_block(D1_);
  free_4d_array(D2_, no_, no_, nv_);
  free_block(t1_);
  free_block(t1old_);
  free_4d_array(t2_, no_, no_, nv_);
  free_4d_array(t2old_, no_, no_, nv_);
  free_4d_array(tau_, no_, no_, nv_);
  free_4d_array(ttau_, no_, no_, nv_);
}

double CCWavefunction::compute_energy() { return 0.0; }

double CCWavefunction::energy()
{
  int no = no_;
  int nv = nv_;
  double **fock = H_->fock_p();
  double ****L = H_->L_p();
  double **t1 = t1_;
  double ****tau = tau_;

  double one_energy=0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      one_energy = 2 * fock[i][a+no] * t1[i][a];

  double two_energy=0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          two_energy += tau[i][j][a][b] * L[i][j][a+no][b+no];

  return one_energy + two_energy;
}

void CCWavefunction::build_tau()
{
  int no = no_;
  int nv = nv_;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          tau_[i][j][a][b] = t2_[i][j][a][b] + t1_[i][a] * t1_[j][b];
          ttau_[i][j][a][b] = t2_[i][j][a][b] + 0.5 * t1_[i][a] * t1_[j][b];
        }
}

} // psi
