#include "ccwavefunction.h"
#include "globals.h"
#include "hamiltonian.h"
#include <boost/shared_ptr.hpp>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <cmath>
#include <libpsio/psio.h>

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
  if(options.get_str("DERTYPE") == "NONE") dertype_ = 0;
  else if(options.get_str("DERTYPE") == "FIRST") dertype_ = 1;

  outfile->Printf("\tWave function  = %s\n", wfn().c_str());
  outfile->Printf("\tMaxiter        = %d\n", maxiter());
  outfile->Printf("\tConvergence    = %3.1e\n", convergence());
  outfile->Printf("\tDIIS           = %s\n", do_diis() ? "Yes" : "No");
  outfile->Printf("\tOut-of-core    = %s\n", ooc() ? "Yes" : "No");
  outfile->Printf("\tDertype        = %d\n", dertype());

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

  Fvv_ = block_matrix(nv,nv);
  Foo_ = block_matrix(no,no);
  Fov_ = block_matrix(no,nv);

  Woooo_ = init_4d_array(no,no,no,no);
  Wovov_ = init_4d_array(no,nv,no,nv);
  Wovvo_ = init_4d_array(no,nv,nv,no);
}

CCWavefunction::~CCWavefunction()
{
  int no = no_;
  int nv = nv_;

  free_block(D1_);
  free_4d_array(D2_, no, no, nv);
  free_block(t1_);
  free_block(t1old_);
  free_4d_array(t2_, no, no, nv);
  free_4d_array(t2old_, no, no, nv);
  free_4d_array(tau_, no, no, nv);
  free_4d_array(ttau_, no, no, nv);
  free_block(Fvv_);
  free_block(Foo_);
  free_block(Fov_);
  free_4d_array(Woooo_, no, no, no);
  free_4d_array(Wovov_, no, nv, no);
  free_4d_array(Wovvo_, no, nv, nv);

  if(dertype()) {
    free_block(l1_);
    free_block(l1old_);
    free_4d_array(l2_, no, no, nv);
    free_4d_array(l2old_, no, no, nv);

    free_block(Hoo_);
    free_block(Hvv_);
    free_block(Hov_);
    free_4d_array(Hoooo_, no, no, no);
    free_4d_array(Hvvvv_, nv, nv, nv);
    free_4d_array(Hovov_, no, nv, no);
    free_4d_array(Hovvo_, no, nv, nv);
    free_4d_array(Hvovv_, nv, no, nv);
    free_4d_array(Hooov_, no, no, no);
    free_4d_array(Hovoo_, no, nv, no);
    free_4d_array(Hvvvo_, nv, nv, nv);

    free_block(Gvv_);
    free_block(Goo_);

    free_block(Doo_);
    free_block(Dvv_);
    free_block(Dov_);
    free_4d_array(Goooo_, no, no, no);
    free_4d_array(Gvvvv_, nv, nv, nv);
    free_4d_array(Gooov_, no, no, no);
    free_4d_array(Gvvvo_, nv, nv, nv);
    free_4d_array(Goovv_, no, no, nv);
    free_4d_array(Govov_, no, nv, no);
  }
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
          two_energy += tau[i][j][a][b]*L[i][j][a+no][b+no];
//          two_energy += (t2[i][j][a][b]+t1[i][a]*t1[j][b])*L[i][j][a+no][b+no];

  return one_energy + two_energy;
}

void CCWavefunction::build_tau()
{
  int no = no_;
  int nv = nv_;
  double **t1 = t1_;
  double ****t2 = t2_;

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          tau_[i][j][a][b] = t2[i][j][a][b] + t1[i][a] * t1[j][b];
          ttau_[i][j][a][b] = t2[i][j][a][b] + 0.5 * t1[i][a] * t1[j][b];
        }
}

void CCWavefunction::amp_save(std::string wfn)
{
  if(wfn == "T") {
    double ****t2tmp = t2_;
    t2_ = t2old_;
    t2old_ = t2tmp;

    double **t1tmp = t1_;
    t1_ = t1old_;
    t1old_ = t1tmp;
  }
  else if(wfn == "L") {
    double ****t2tmp = l2_;
    l2_ = l2old_;
    l2old_ = t2tmp;

    double **t1tmp = l1_;
    l1_ = l1old_;
    l1old_ = t1tmp;
  }
}

void CCWavefunction::build_F()
{
  int no = no_;
  int nv = nv_;
  double **fock = H_->fock_p();
  double ****L = H_->L_p();
  double **t1 = t1old_;
  double ****ttau = ttau_;

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = fock[a+no][e+no];
      for(int m=0; m < no; m++) {
        value -= 0.5*fock[m][e+no]*t1[m][a];
        for(int f=0; f < nv; f++) {
          value += t1[m][f]*L[m][a+no][f+no][e+no];
          for(int n=0; n < no; n++)
            value -= ttau[m][n][a][f]*L[m][n][e+no][f+no];
        }
      }
      Fvv_[a][e] = value;
    }

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = fock[m][i];
      for(int e=0; e < nv; e++) {
        value += 0.5*t1[i][e]*fock[m][e+no];
        for(int n=0; n < no; n++) {
          value += t1[n][e]*L[m][n][i][e+no];
          for(int f=0; f < nv; f++)
            value += ttau[i][n][e][f]*L[m][n][e+no][f+no];
        }
      }
      Foo_[m][i] = value;
    }

  for(int m=0; m < no; m++)
    for(int e=0; e < nv; e++) {
      double value = fock[m][e+no];
      for(int n=0; n < no; n++)
        for(int f=0; f < nv; f++)
          value += t1[n][f]*L[m][n][e+no][f+no];
      Fov_[m][e] = value;
    }
}

void CCWavefunction::build_W()
{
  int no = no_;
  int nv = nv_;
  double ****tau = tau_;
  double **t1 = t1old_;
  double ****t2 = t2old_;
  double ****ints = H_->ints_p();
  double ****L = H_->L_p();

  for(int m=0; m < no; m++)
    for(int n=0; n < no; n++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          double value = ints[m][n][i][j];
          for(int e=0; e < nv; e++) {
            value += t1[j][e]*ints[m][n][i][e+no] +
                     t1[i][e]*ints[m][n][e+no][j];
            for(int f=0; f < nv; f++)
              value += tau[i][j][e][f]*ints[m][n][e+no][f+no];
          }
          Woooo_[m][n][i][j] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int j=0; j < no; j++)
        for(int e=0; e < nv; e++) {
          double value = -ints[m][b+no][j][e+no];
          for(int f=0; f < nv; f++)
            value -= t1[j][f]*ints[m][b+no][f+no][e+no];
          for(int n=0; n < no; n++)
            value += t1[n][b]*ints[m][n][j][e+no];
          for(int n=0; n < no; n++) {
            for(int f=0; f < nv; f++)
              value += ints[m][n][f+no][e+no]*
                (0.5*t2[j][n][f][b] + t1[j][f]*t1[n][b]);
          }
          Wovov_[m][b][j][e] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int j=0; j < no; j++) {
          double value = ints[m][b+no][e+no][j];
          for(int f=0; f < nv; f++)
            value += t1[j][f]*ints[m][b+no][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][e+no][j];
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++)
              value -= ints[m][n][e+no][f+no]*
                (0.5*t2[j][n][f][b] + t1[j][f]*t1[n][b]);
          for(int n=0; n < no; n++) {
            for(int f=0; f < nv; f++)
              value += 0.5*L[m][n][e+no][f+no]*t2[n][j][f][b];
          }
          Wovvo_[m][b][e][j] = value;
        }
}

void CCWavefunction::build_t1()
{
  int no = no_;
  int nv = nv_;
  double **t1new = t1_;
  double **t1 = t1old_;
  double **fock = H_->fock_p();
  double ****ints = H_->ints_p();
  double ****L = H_->L_p();
  double **Fae = Fvv_;
  double **Fmi = Foo_;
  double **Fme = Fov_;
  double ****t2 = t2old_;

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      double value = fock[a+no][i];
      for(int e=0; e < nv; e++)
        value += t1[i][e]*Fae[a][e];
      for(int m=0; m < no; m++)
        value -= t1[m][a]*Fmi[m][i];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          value += (2*t2[i][m][a][e]-t2[i][m][e][a])*Fme[m][e];
      for(int n=0; n < no; n++)
        for(int f=0; f < nv; f++)
          value += t1[n][f]*L[n][a+no][f+no][i];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            value += (2*t2[m][i][e][f]-t2[m][i][f][e])*ints[m][a+no][e+no][f+no];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int n=0; n < no; n++)
            value -= t2[m][n][a][e]*L[n][m][e+no][i];
      t1new[i][a] = value;
   }
}

void CCWavefunction::build_t2()
{
  int no = no_;
  int nv = nv_;
  double ****t2new = t2_;
  double ****t2 = t2old_;
  double ****tau = tau_;
  double ****Wmnij = Woooo_;
  double ****Wmbej = Wovvo_;
  double ****Wmbje = Wovov_;
  double **t1 = t1old_;
  double **Fme = Fov_;
  double **Fmi = Foo_;
  double **Fae = Fvv_;
  double ****ints = H_->ints_p();

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          double value = ints[i][j][a+no][b+no];
          for(int e=0; e < nv; e++)
            value += t2[i][j][a][e]*Fae[b][e] + t2[j][i][b][e]*Fae[a][e];
          for(int e=0; e < nv; e++)
            for(int m=0; m < no; m++)
              value -= 0.5*(t2[i][j][a][e]*t1[m][b]*Fme[m][e] +
                      t2[i][j][e][b]*t1[m][a]*Fme[m][e]);
          for(int m=0; m < no; m++)
            value -=  t2[i][m][a][b]*Fmi[m][j] + t2[m][j][a][b]*Fmi[m][i];
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              value -= 0.5*(t2[i][m][a][b]*t1[j][e]*Fme[m][e] +
                      t2[m][j][a][b]*t1[i][e]*Fme[m][e]);
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              value += tau[m][n][a][b]*Wmnij[m][n][i][j];
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              value += tau[i][j][e][f]*ints[a+no][b+no][e+no][f+no];
          for(int e=0; e < nv; e++)
            value += t1[i][e]*ints[a+no][b+no][e+no][j] +
                     t1[j][e]*ints[b+no][a+no][e+no][i];
          for(int m=0; m < no; m++)
            value -= t1[m][a]*ints[m][b+no][i][j]+t1[m][b]*ints[m][a+no][j][i];
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++) {
              value += (t2[i][m][a][e] - t2[i][m][e][a]) * Wmbej[m][b][e][j];
              value += t2[i][m][a][e] * (Wmbej[m][b][e][j] + Wmbje[m][b][j][e]);
              value += t2[m][j][a][e] * Wmbje[m][b][i][e];
              value += t2[i][m][e][b] * Wmbje[m][a][j][e];
              value += t2[j][m][b][e] * (Wmbej[m][a][e][i] + Wmbje[m][a][i][e]);
              value += (t2[j][m][b][e] - t2[j][m][e][b]) * Wmbej[m][a][e][i];
            }
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++) {
              value -= t1[i][e]*t1[m][a]*ints[m][b+no][e+no][j];
              value -= t1[i][e]*t1[m][b]*ints[m][a+no][j][e+no];
              value -= t1[j][e]*t1[m][a]*ints[m][b+no][i][e+no];
              value -= t1[j][e]*t1[m][b]*ints[m][a+no][e+no][i];
            }
          t2new[i][j][a][b] = value;
        }

  double ****Zmbij = init_4d_array(no, nv, no, no);
  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          Zmbij[m][b][i][j] = 0.0;
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              Zmbij[m][b][i][j] += ints[m][b+no][e+no][f+no] * tau[i][j][e][f];
        }

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          double value = 0.0;
          for(int m=0; m < no; m++)
            value -= t1[m][a]*Zmbij[m][b][i][j];
          t2new[i][j][a][b] += value;
          t2new[j][i][b][a] += value;
        }
  free_4d_array(Zmbij,no,nv,no);
}

double CCWavefunction::t1norm()
{
  int no = no_;
  int nv = nv_;
  double **t1 = t1_;

  double diag = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      diag += t1[i][a] * t1[i][a];

  return sqrt(diag/(2*no));
}

void CCWavefunction::diis(int iter, std::string wfn)
{
  int nvector=8;  /* Number of error vectors to keep */
  int word;
  int diis_cycle;
  int vector_length;
  double *error;
  div_t fraction;
  double **B, *C, **vector;
  double product, determinant, maximum;
  psio_address start, end;
  int error_file, amp_file;
  if(wfn == "T") {
    error_file = 90;
    amp_file = 91;
  } 
  else if(wfn == "L") {
    error_file = 92;
    amp_file = 93;
  }

  int no = no_; 
  int nv = nv_;

  double **t1, **t1old, ****t2, ****t2old;
  if(wfn == "T") {
    t1 = t1_;
    t1old = t1old_;
    t2 = t2_;
    t2old = t2old_;
  }
  else if(wfn == "L") {
    t1 = l1_;
    t1old = l1old_;
    t2 = l2_;
    t2old = l2old_;
  }

  /* Calculate the length of a single error vector */
  vector_length = no*nv + no*no*nv*nv;

  /* If we haven't already, open the vector files for reading/writing */
  if(iter == 1) { 
    psio_open(error_file, PSIO_OPEN_NEW);
    psio_open(amp_file, PSIO_OPEN_NEW);
  }

  /* Set the diis cycle value */
  fraction = div((iter-1),nvector);
  diis_cycle = fraction.rem;

  /* Build the current error vector and dump it to disk */
  error = init_array(vector_length);
  word=0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      error[word++] = t1[i][a] - t1old[i][a];
    }

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          error[word++] = t2[i][j][a][b] - t2old[i][j][a][b];
        }

  start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
  psio_write(error_file, "DIIS Error Vectors", (char *) error, 
	     vector_length*sizeof(double), start, &end);

  /* Store the amplitudes, too */
  word=0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      error[word++] = t1[i][a];

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)  {
          error[word++] = t2[i][j][a][b];
        }

  start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
  psio_write(amp_file, "DIIS Amplitude Vectors", (char *) error, 
	     vector_length*sizeof(double), start, &end);
  
  free(error);
    
  /* If we haven't run through enough iterations, set the correct dimensions
     for the extrapolation */
  if(!(iter >= (nvector))) {
    if(iter < 2) return; /* Leave if we can't extrapolate at all */
    nvector = iter;
  }

  /* Now grab the full set of error vectors from the file */
  vector = init_matrix(nvector, vector_length);
  for(int p=0; p < nvector; p++) {
    start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));
    psio_read(error_file, "DIIS Error Vectors", (char *) vector[p], 
	      vector_length*sizeof(double), start, &end);
  }

  /* Build B matrix of error vector products */
  B = init_matrix(nvector+1,nvector+1);

  for(int p=0; p < nvector; p++)
    for(int q=0; q < nvector; q++) {
      dot_arr(vector[p], vector[q], vector_length, &product); 
      B[p][q] = product;
    }

  for(int p=0; p < nvector; p++) {
    B[p][nvector] = -1;
    B[nvector][p] = -1;
  }

  B[nvector][nvector] = 0;

  /* Find the maximum value in B and scale all its elements */
  maximum = fabs(B[0][0]);
  for(int p=0; p < nvector; p++)
    for(int q=0; q < nvector; q++)
      if(fabs(B[p][q]) > maximum) maximum = fabs(B[p][q]);

  for(int p=0; p < nvector; p++)
    for(int q=0; q < nvector; q++)
      B[p][q] /= maximum; 

  /* Build the constant vector */
  C = init_array(nvector+1);
  C[nvector] = -1;

  /* Solve the linear equations */
  flin(B, C, nvector+1, 1, &determinant);

  /* Grab the old amplitude vectors */
  for(int p=0; p < nvector; p++) {
    start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));
    psio_read(amp_file, "DIIS Amplitude Vectors", (char *) vector[p], 
	      vector_length*sizeof(double), start, &end);
  }
  
  /* Build the new amplitude vector from the old ones */
  word=0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      t1[i][a] = 0.0;
      for(int p=0; p < nvector; p++)
        t1[i][a] += C[p]*vector[p][word];
        word++;
    }

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          t2[i][j][a][b] = 0.0;
          for(int p=0; p < nvector; p++)
            t2[i][j][a][b] += C[p]*vector[p][word];
            word++;
        }

  /* Release memory and return */
  free_matrix(vector, nvector);
  free_matrix(B, nvector+1);
  free(C);

  return;
}

double CCWavefunction::increment_amps(std::string wfn)
{
  int no = no_;
  int nv = nv_;
  double **D1 = D1_;
  double ****D2 = D2_;
  double **t1, **t1old, ****t2, ****t2old;
  if(wfn == "T") {
    t1 = t1_;
    t1old = t1old_;
    t2 = t2_;
    t2old = t2old_;
  }
  else if(wfn == "L") {
    t1 = l1_;
    t1old = l1old_;
    t2 = l2_;
    t2old = l2old_;
  }

  double residual1 = 0.0;
  double residual2 = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      residual1 += t1[i][a] * t1[i][a];     
      t1[i][a] = t1old[i][a] + t1[i][a]/D1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          residual2 += t2[i][j][a][b] * t2[i][j][a][b];
          t2[i][j][a][b] = t2old[i][j][a][b] + t2[i][j][a][b]/D2[i][j][a][b];
        }
    }

  return sqrt(residual1 + residual2);
}

void CCWavefunction::hbar()
{
  int no = no_;
  int nv = nv_;

  Hvv_ = block_matrix(nv, nv);
  Hoo_ = block_matrix(no, no);
  Hov_ = block_matrix(no, nv);
  Hoooo_ = init_4d_array(no, no, no, no);
  Hvvvv_ = init_4d_array(nv, nv, nv, nv);
  Hovov_ = init_4d_array(no, nv, no, nv);
  Hovvo_ = init_4d_array(no, nv, nv, no);
  Hvovv_ = init_4d_array(nv, no, nv, nv);
  Hooov_ = init_4d_array(no, no, no, nv);
  Hovoo_ = init_4d_array(no, nv, no, no);
  Hvvvo_ = init_4d_array(nv, nv, nv, no);

  double **fock = H_->fock_p();
  double ****ints = H_->ints_p();
  double ****L = H_->L_p();
  double **t1 = t1_;
  double ****t2 = t2_;
  double ****tau = tau_;

  for(int m=0; m < no; m++)
    for(int e=0; e < nv; e++) {
      double value = fock[m][e+no];
      for(int n=0; n < no; n++)
        for(int f=0; f < nv; f++)
          value += t1[n][f]*L[m][n][e+no][f+no];
      Hov_[m][e] = value;
    }

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = fock[m][i];
      for(int e=0; e < nv; e++) {
        value += t1[i][e]*fock[m][e+no];
        for(int n=0; n < no; n++) {
          value += t1[n][e]*L[m][n][i][e+no];
          for(int f=0; f < nv; f++)
            value += tau[i][n][e][f]*L[m][n][e+no][f+no];
        }
      }
      Hoo_[m][i] = value;
    }

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = fock[a+no][e+no];
      for(int m=0; m < no; m++) {
        value -= fock[m][e+no]*t1[m][a];
        for(int f=0; f < nv; f++) {
          value += t1[m][f]*L[a+no][m][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= tau[m][n][f][a]*L[m][n][f+no][e+no];
        }
      }
      Hvv_[a][e] = value;
    }

  for(int m=0; m < no; m++)
    for(int n=0; n < no; n++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          double value = ints[m][n][i][j];
          for(int e=0; e < nv; e++) {
            value += 2.0*t1[j][e]*ints[m][n][i][e+no];
            for(int f=0; f < nv; f++)
              value += tau[i][j][e][f]*ints[m][n][e+no][f+no];
          }
          Hoooo_[m][n][i][j] = value;
        }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++) {
          double value = ints[a+no][b+no][e+no][f+no];
          for(int m=0; m < no; m++) {
            value -= 2.0*t1[m][b]*ints[a+no][m][e+no][f+no];
            for(int n=0; n < no; n++)
              value += tau[m][n][a][b]*ints[m][n][e+no][f+no];
          }
          Hvvvv_[a][b][e][f] = value;
        }

  for(int a=0; a < nv; a++)
    for(int m=0; m < no; m++)
      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++) {
          double value = ints[a+no][m][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][a] * ints[n][m][e+no][f+no];
          Hvovv_[a][m][e][f] = value;
        }

  for(int m=0; m < no; m++)
    for(int n=0; n < no; n++)
      for(int i=0; i < no; i++)
        for(int e=0; e < nv; e++) {
          double value = ints[m][n][i][e+no];
          for(int f=0; f < nv; f++)
            value += t1[i][f] * ints[n][m][e+no][f+no];
          Hooov_[m][n][i][e] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int j=0; j < no; j++) {
          double value = ints[m][b+no][e+no][j];
          for(int f=0; f < nv; f++)
            value += t1[j][f]*ints[m][b+no][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][e+no][j];
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++)
              value -= ints[m][n][e+no][f+no]*tau[n][j][b][f];
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++)
              value += t2[n][j][f][b]*L[m][n][e+no][f+no];
          Hovvo_[m][b][e][j] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int j=0; j < no; j++)
        for(int e=0; e < nv; e++) {
          double value = ints[m][b+no][j][e+no];
          for(int f=0; f < nv; f++)
            value += t1[j][f]*ints[b+no][m][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][j][e+no];
          for(int n=0; n < no; n++) {
            for(int f=0; f < nv; f++)
              value -= ints[n][m][e+no][f+no]*tau[j][n][f][b];
          }
          Hovov_[m][b][j][e] = value;
        }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int i=0; i < no; i++) {
          double value = ints[a+no][b+no][e+no][i];
          for(int f=0; f < nv; f++)
            value += t1[i][f]*ints[a+no][b+no][e+no][f+no];
          for(int m=0; m < no; m++)
            value -= t1[m][b]*ints[a+no][m][e+no][i] +
                     t1[m][a]*ints[b+no][m][i][e+no];
          for(int m=0; m < no; m++)
            for(int f=0; f < nv; f++)
              value -= tau[i][m][f][a]*ints[m][b+no][e+no][f+no] +
                       tau[i][m][f][b]*ints[a+no][m][e+no][f+no];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              value += tau[m][n][a][b]*ints[m][n][e+no][i];
          for(int m=0; m < no; m++)
            value -= fock[m][e+no]*t2[m][i][a][b];
          for(int m=0; m < no; m++)
            for(int f=0; f < nv; f++)
              value += t2[m][i][f][b]*L[a+no][m][e+no][f+no];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += ints[m][n][e+no][f+no]*(t1[i][f]*t2[m][n][a][b]+t1[m][a]*t2[i][n][f][b]+t1[n][b]*t2[m][i][a][f]);
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= (t1[m][f]*t2[n][i][a][b]+t1[n][a]*t2[m][i][f][b])*L[m][n][f+no][e+no];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += t1[i][f]*t1[m][a]*t1[n][b]*ints[m][n][e+no][f+no];

          Hvvvo_[a][b][e][i] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          double value = ints[m][b+no][i][j];
          for(int e=0; e < nv; e++)
            value += t1[j][e]*ints[m][b+no][i][e+no] +
                     t1[i][e]*ints[b+no][m][j][e+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][i][j];
          for(int n=0; n < no; n++)
            for(int e=0; e < nv; e++)
              value -= tau[i][n][e][b]*ints[n][m][j][e+no] +
                       tau[j][n][e][b]*ints[m][n][i][e+no];
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              value += ints[m][b+no][e+no][f+no]*tau[i][j][e][f];
          for(int e=0; e < nv; e++)
            value += fock[m][e+no]*t2[i][j][e][b];
          for(int n=0; n < no; n++)
            for(int e=0; e < nv; e++)
              value += L[m][n][i][e+no]*t2[n][j][e][b];
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= ints[m][n][e+no][f+no]*(t1[j][f]*t2[i][n][e][b]+t1[i][e]*t2[j][n][f][b]+t1[n][b]*t2[i][j][e][f]);
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += (t1[i][e]*t2[n][j][f][b]+t1[n][f]*t2[i][j][e][b])*L[m][n][e+no][f+no];
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= t1[j][f]*t1[i][e]*t1[n][b]*ints[m][n][e+no][f+no];

          Hovoo_[m][b][i][j] = value;
        }

}

void CCWavefunction::init_lambda()
{
  int no = no_;
  int nv = nv_;
  double **t1 = t1_;
  double ****t2 = t2_;

  l1_ = block_matrix(no,nv);
  l1old_ = block_matrix(no,nv);

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      l1_[i][a] = 2*t1[i][a];

  l2_ = init_4d_array(no,no,nv,nv);
  l2old_ = init_4d_array(no,no,nv,nv);

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          l2_[i][j][a][b] = 2*(2*t2[i][j][a][b]-t2[i][j][b][a]);

  Gvv_ = block_matrix(nv, nv);
  Goo_ = block_matrix(no, no);
}

void CCWavefunction::build_G()
{
  int no = no_;
  int nv = nv_;
  double ****t2 = t2_;
  double ****l2 = l2old_;

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = 0.0;
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          for(int j=0; j < no; j++)
            value += t2[m][j][a][b] * l2[i][j][a][b];
      Goo_[m][i] = value;
    }

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = 0.0;
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++)
          for(int b=0; b < nv; b++)
            value -= t2[i][j][e][b] * l2[i][j][a][b];
      Gvv_[a][e] = value;
    }
}

void CCWavefunction::build_l1()
{
  int no = no_;
  int nv = nv_;
  double **l1 = l1old_;
  double ****l2 = l2old_;
  double **l1new = l1_;
  double **Gvv = Gvv_;
  double **Goo = Goo_;
  double **Hov = Hov_;
  double **Hvv = Hvv_;
  double **Hoo = Hoo_;
  double ****Hvvvo = Hvvvo_;
  double ****Hovoo = Hovoo_;
  double ****Hovvo = Hovvo_;
  double ****Hovov = Hovov_;
  double ****Hvovv = Hvovv_;
  double ****Hooov = Hooov_;

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      double value = 2 * Hov[i][a];

      for(int e=0; e < nv; e++)
        value += l1[i][e] * Hvv[e][a];

      for(int m=0; m < no; m++)
        value -= l1[m][a] * Hoo[i][m];

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          value += l1[m][e] * (2*Hovvo[i][e][a][m] - Hovov[i][e][m][a]);

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            value += l2[i][m][e][f] * Hvvvo[e][f][a][m];

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            value -= l2[m][n][a][e] * Hovoo[i][e][m][n];

      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++)
          value -= Gvv[e][f] * (2*Hvovv[e][i][f][a] - Hvovv[e][i][a][f]);

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          value -= Goo[m][n] * (2*Hooov[m][i][n][a] - Hooov[i][m][n][a]);

      l1new[i][a] = value;
    }
}

void CCWavefunction::build_l2()
{
  int no = no_;
  int nv = nv_;
  double **l1 = l1old_;
  double ****l2 = l2old_;
  double ****l2new = l2_;
  double ****L = H_->L_p();
  double **Gvv = Gvv_;
  double **Goo = Goo_;
  double **Hov = Hov_;
  double **Hvv = Hvv_;
  double **Hoo = Hoo_;

  double ****Hovvo = Hovvo_;
  double ****Hovov = Hovov_;
  double ****Hoooo = Hoooo_;
  double ****Hvvvv = Hvvvv_;
  double ****Hvovv = Hvovv_;
  double ****Hooov = Hooov_;

  double ****Z = init_4d_array(no, no, nv, nv);

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          double value = L[i][j][a+no][b+no];

          value += 2.0*l1[i][a]*Hov[j][b] - l1[j][a]*Hov[i][b];

          for(int e=0; e < nv; e++)
            value += l2[i][j][e][b]*Hvv[e][a];

          for(int m=0; m < no; m++)
            value -= l2[m][j][a][b]*Hoo[i][m];

          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              value += 0.5 * l2[m][n][a][b] * Hoooo[i][j][m][n];

          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              value += 0.5 * l2[i][j][e][f] * Hvvvv[e][f][a][b];

          for(int e=0; e < nv; e++)
              value += l1[i][e]*(2*Hvovv[e][j][a][b] - Hvovv[e][j][b][a]);

          for(int m=0; m < no; m++)
              value -= l1[m][b]*(2*Hooov[j][i][m][a] - Hooov[i][j][m][a]);

          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++) {
              value += (2*Hovvo[i][e][a][m] - Hovov[i][e][m][a])*l2[m][j][e][b];
              value -= Hovov[j][e][m][a]*l2[m][i][b][e];
              value -= Hovvo[j][e][a][m]*l2[m][i][e][b];
            }

          for(int e=0; e < nv; e++)
            value += Gvv[a][e]*L[i][j][e+no][b+no];
          for(int m=0; m < no; m++)
            value -= Goo[m][i]*L[m][j][a+no][b+no];

          Z[i][j][a][b] = value;
    }

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          l2new[i][j][a][b] = Z[i][j][a][b] + Z[j][i][b][a];

   free_4d_array(Z, no, no, nv);
}

/*
** pseudoenergy(): Evaluates an energy-like expression for the Lambda
doubles
** amplitudes: 
**   E = <0|L2 H|0> = 1/2 <ab|ij> L(ij,ab)
** This expression is derived in the UGA formalism.
*/
double CCWavefunction::pseudoenergy()
{
  int no = no_;
  int nv = nv_;
  double ****ints = H_->ints_p();
  double ****l2 = l2_;

  double energy=0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          energy += 0.5*ints[i][j][a+no][b+no]*l2[i][j][a][b];

  return energy;
}

void CCWavefunction::init_density()
{
  int no = no_;
  int nv = nv_;

  Doo_ = block_matrix(no, no);
  Dvv_ = block_matrix(nv, nv);
  Dov_ = block_matrix(no, nv);
  Dvo_ = block_matrix(nv, no);
  Goooo_ = init_4d_array(no, no, no, no);
  Gvvvv_ = init_4d_array(nv, nv, nv, nv);
  Gooov_ = init_4d_array(no, no, no, nv);
  Gvvvo_ = init_4d_array(nv, nv, nv, no);
  Govov_ = init_4d_array(no, nv, no, nv);
  Goovv_ = init_4d_array(no, no, nv, nv);
}

double CCWavefunction::onepdm()
{
  int no = no_;
  int nv = nv_;
  double **t1 = t1_;
  double ****t2 = t2_;
  double ****tau = tau_;
  double **l1 = l1_;
  double ****l2 = l2_;
  double **fock = H_->fock_p();

  double **Doo = Doo_;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++) {
      for(int e=0; e < nv; e++)
        Doo[i][j] -= t1[i][e] * l1[j][e];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            Doo[i][j] -= t2[i][m][e][f] * l2[j][m][e][f];
    }

  double **Dvv = Dvv_;
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++) {
      for(int m=0; m < no; m++)
        Dvv[a][b] += t1[m][b] * l1[m][a];
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            Dvv[a][b] += t2[m][n][b][e] * l2[m][n][a][e];
    }

  double **Dvo = Dvo_;
  for(int a=0; a < nv; a++)
    for(int i=0; i < no; i++) {
      Dvo[a][i] = l1[i][a];
    }

  double **Dov = Dov_;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      Dov[i][a] += 2.0 * t1[i][a];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          Dov[i][a] += l1[m][e] * (2.0 * t2[i][m][a][e] - tau[m][i][a][e]);
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              Dov[i][a] -= l2[m][n][e][f] * (t1[m][a] * t2[i][n][e][f] + t1[i][e] * t2[m][n][a][f]);
    }

  double energy = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      energy += fock[i][j] * Doo[i][j];
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
    energy += fock[a+no][b+no] * Dvv[a][b];

  return energy;
}

double CCWavefunction::twopdm()
{
  int no = no_;
  int nv = nv_;
  double **t1 = t1_;
  double ****t2 = t2_;
  double ****tau = tau_;
  double **l1 = l1_;
  double ****l2 = l2_;
  double ****ints = H_->ints_p();

  double ****Goooo = Goooo_;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++) {
          Goooo[i][j][k][l] = 0.0;
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              Goooo[i][j][k][l] += tau[i][j][e][f] * l2[k][l][e][f];
        }

  double Eoooo = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++)
          Eoooo += 0.5 * ints[i][j][k][l] * Goooo[i][j][k][l];
  outfile->Printf("\tOOOO Energy = %20.14f\n", Eoooo);

  double ****Gvvvv = Gvvvv_;
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++)
        for(int d=0; d < nv; d++) {
          Gvvvv[a][b][c][d] = 0.0;
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              Gvvvv[a][b][c][d] += tau[m][n][a][b] * l2[m][n][c][d];
        }

  double Evvvv = 0.0;
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++)
        for(int d=0; d < nv; d++)
          Evvvv += 0.5 * ints[a+no][b+no][c+no][d+no] * Gvvvv[a][b][c][d];
  outfile->Printf("\tVVVV Energy = %20.14f\n", Evvvv);

  double ****Gooov = Gooov_;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int a=0; a < nv; a++) {
          for(int e=0; e < nv; e++)
            Gooov[i][j][k][a] -= l1[k][e] * (2.0 * tau[i][j][e][a] - tau[i][j][a][e]);
          for(int e=0; e < nv; e++)
            Gooov[i][j][k][a] -= t1[i][e] * l2[j][k][a][e];
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              for(int f=0; f < nv; f++)
                Gooov[i][j][k][a] -= l2[k][m][e][f] *
                           (2.0 * (t1[j][a] * t2[i][m][e][f] + t1[i][e] * t2[j][m][a][f])
                                - (t1[i][a] * t2[j][m][e][f] + t1[j][e] * t2[i][m][a][f])
                                - (t1[m][a] * t2[i][j][e][f] + t1[i][e] * t2[m][j][a][f] + t1[j][f] * t2[i][m][e][a]));
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              for(int f=0; f < nv; f++)
                Gooov[i][j][k][a] += l2[k][m][e][f] * t1[m][a] * t1[i][e] * t1[j][f];

        }

  double Eooov = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int a=0; a < nv; a++)
          Eooov += ints[i][j][k][a+no] * Gooov[i][j][k][a];
  outfile->Printf("\tOOOV Energy = %20.14f\n", Eooov);

  double ****Gvvvo = Gvvvo_;
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++)
        for(int i=0; i < no; i++) {
          for(int m=0; m < no; m++)
            Gvvvo[a][b][c][i] += l1[m][c] * (2.0 * tau[m][i][a][b] - tau[i][m][a][b]);
          for(int m=0; m < no; m++)
            Gvvvo[a][b][c][i] += t1[m][a] * l2[i][m][b][c];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                Gvvvo[a][b][c][i] += l2[n][m][c][e] * 
                           (2.0 * (t1[n][a] * t2[i][m][b][e] + t1[i][b] * t2[n][m][a][e])
                                - (t1[n][b] * t2[i][m][a][e] + t1[i][a] * t2[n][m][b][e])
                                - (t1[m][b] * t2[n][i][a][e] + t1[n][a] * t2[m][i][b][e] + t1[i][e] * t2[n][m][a][b]));
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                Gvvvo[a][b][c][i] -= l2[n][m][c][e] * t1[m][b] * t1[n][a] * t1[i][e];
        }

  double Evvvo = 0.0;
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++)
        for(int i=0; i < no; i++)
          Evvvo += ints[a+no][b+no][c+no][i] * Gvvvo[a][b][c][i];
  outfile->Printf("\tVVVO Energy = %20.14f\n", Evvvo);

  double ****Govov = Govov_;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          Govov[i][a][j][b] -= t1[i][a] * l1[j][b];
          for(int m=0; m < no; m++) 
            for(int e=0; e < nv; e++)
              Govov[i][a][j][b] -= (tau[m][i][b][e] * l2[j][m][e][a] + t2[i][m][b][e] * l2[m][j][e][a]);
        }

  double Eovov = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          Eovov += ints[i][a+no][j][b+no] * Govov[i][a][j][b];
  outfile->Printf("\tOVOV Energy = %20.14f\n", Eovov);

  double ****Goovv = Goovv_;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          Goovv[i][j][a][b] += 4.0 * t1[i][a] * l1[j][b];
          Goovv[i][j][a][b] += 4.0 * tau[i][j][a][b] - 2.0 * tau[i][j][b][a];
          Goovv[i][j][a][b] += l2[i][j][a][b];
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              Goovv[i][j][a][b] += 2.0 * l1[m][e] * 
                         (2.0 * t1[i][a] * (2.0 * t2[j][m][b][e] - t2[j][m][e][b])
                              - t1[i][b] * (2.0 * t2[j][m][a][e] - t2[j][m][e][a])
                              - t1[i][e] * (2.0 * tau[j][m][b][a] - tau[j][m][a][b])
                              - t1[m][a] * (2.0 * t2[i][j][e][b] - t2[i][j][b][e]));
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              Goovv[i][j][a][b] += 4.0 * t2[i][m][a][e] * l2[m][j][e][b] - 2.0 * tau[m][j][b][e] * l2[i][m][a][e];

          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                for(int f=0; f < nv; f++)
                  Goovv[i][j][a][b] += 2.0 * l2[m][n][e][f] * 
                      (0.5 * (t2[m][n][a][b] * t2[i][j][e][f] + t2[m][i][a][e] * t2[n][j][b][f] + t2[n][j][a][e] * t2[i][m][f][b])
                     + 2.0 * (- tau[i][j][a][e] * t2[m][n][b][f] - tau[i][m][a][b] * t2[j][n][e][f] 
                              - tau[m][j][b][e] * t2[i][n][a][f] + t2[i][m][a][e] * t2[j][n][b][f])
                           - (- tau[i][j][b][e] * t2[m][n][a][f] - tau[i][m][b][a] * t2[j][n][e][f] 
                              - tau[m][j][a][e] * t2[i][n][b][f] + t2[i][m][b][e] * t2[j][n][a][f]));

          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                for(int f=0; f < nv; f++)
                  Goovv[i][j][a][b] += l2[m][n][e][f] *
                         (t1[m][a] * t1[n][b] * t2[i][j][e][f] + t1[i][e] * t1[j][f] * t2[m][n][a][b]
                        + t1[m][a] * t1[i][e] * t2[n][j][b][f] + t1[n][b] * t1[j][f] * t2[m][i][a][e]
                        + t1[n][a] * t1[j][e] * t2[i][m][f][b] + t1[i][f] * t1[m][b] * t2[n][j][a][e]);

          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                for(int f=0; f < nv; f++)
                  Goovv[i][j][a][b] += l2[m][n][e][f] * t1[m][a] * t1[n][b] * t1[i][e] * t1[j][f];

        }

  double Eoovv = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          Eoovv += 0.5 * ints[i][j][a+no][b+no] * Goovv[i][j][a][b];
  outfile->Printf("\tOOVV Energy = %20.14f\n", Eoovv);
//  outfile->Printf("\tOOVV + OVOV = %20.14f\n", Eoovv+Eovov);

  return Eoooo+Evvvv+Eooov+Evvvo+Eovov+Eoovv;
}

} // psi
