#include "ccpert.h"

#include <boost/shared_ptr.hpp>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libdiis/diismanager.h>
#include "array.h"

namespace psi { namespace ugacc {

CCPert::CCPert(double **pert, double omega, boost::shared_ptr<CCWfn> CC, boost::shared_ptr<HBAR> HBAR)
{
  CC_ = CC;
  HBAR_ = HBAR;
  pert_ = pert;
  omega_ = omega;

  H_ = CC_->H_;
  no_ = CC_->no_;
  nv_ = CC_->nv_;

  int no = no_;
  int nv = nv_;

  Xov_ = block_matrix(no, nv);
  Xoo_ = block_matrix(no, no);
  Xvv_ = block_matrix(nv, nv);
  Xvo_ = block_matrix(nv, no);
  Xovoo_ = init_4d_array(no, nv, no, no);
  Xvvvo_ = init_4d_array(nv, nv, nv, no);
  Xvvoo_ = init_4d_array(nv, nv, no, no);

  pertbar();

  // Denominators
  D1_ = block_matrix(no, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      D1_[i][a] = HBAR_->Hoo_[i][i] - HBAR_->Hvv_[a][a];

  D2_ = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          D2_[i][j][a][b] = HBAR_->Hoo_[i][i] + HBAR_->Hoo_[j][j] - HBAR_->Hvv_[a][a] - HBAR_->Hvv_[b][b];

  // Initial guess perturbed amplitudes
  X1_ = block_matrix(no, nv);
  X1old_ = block_matrix(no, nv);
  Y1_ = block_matrix(no, nv);
  Y1old_ = block_matrix(no, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      X1_[i][a] = Xvo_[a][i]/(D1_[i][a] + omega_);
      Y1_[i][a] = Xvo_[a][i]/(D1_[i][a] + omega_);
    }

  X2_ = init_4d_array(no, no, nv, nv);
  X2old_ = init_4d_array(no, no, nv, nv);
  Y2_ = init_4d_array(no, no, nv, nv);
  Y2old_ = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          X2_[i][j][a][b] = (Xvvoo_[a][b][i][j]+Xvvoo_[b][a][j][i])/(D2_[i][j][a][b] + omega_);
          Y2_[i][j][a][b] = (Xvvoo_[a][b][i][j]+Xvvoo_[b][a][j][i])/(D2_[i][j][a][b] + omega_);
        }

  // DIIS Vectors
  X1diis_.resize(no_*nv_);
  X2diis_.resize(no_*no_*nv_*nv_);
  X1err_.resize(no_*nv_);
  X2err_.resize(no_*no_*nv_*nv_);
  Y1diis_.resize(no_*nv_);
  Y2diis_.resize(no_*no_*nv_*nv_);
  Y1err_.resize(no_*nv_);
  Y2err_.resize(no_*no_*nv_*nv_);
}

CCPert::~CCPert()
{
  int no = no_;
  int nv = nv_;

  free_block(Xov_);
  free_block(Xoo_);
  free_block(Xvv_);
  free_block(Xvo_);
  free_4d_array(Xovoo_, no, nv, no);
  free_4d_array(Xvvvo_, nv, nv, nv);
  free_4d_array(Xvvoo_, nv, nv, no);

  free_block(X1_);
  free_block(X1old_);
  free_4d_array(X2_, no, no, nv);
  free_4d_array(X2old_, no, no, nv);

  free_block(Y1_);
  free_block(Y1old_);
  free_4d_array(Y2_, no, no, nv);
  free_4d_array(Y2old_, no, no, nv);

  free_block(D1_);
  free_4d_array(D2_, no, no, nv);

  X1diis_.resize(0);
  X2diis_.resize(0);
  X1err_.resize(0);
  X2err_.resize(0);
  Y1diis_.resize(0);
  Y2diis_.resize(0);
  Y1err_.resize(0);
  Y2err_.resize(0);
}

void CCPert::pertbar()
{
  int no = no_;
  int nv = nv_;
  double **t1 = CC_->t1_;
  double ****t2 = CC_->t2_;

  for(int m=0; m < no; m++)
    for(int e=0; e < nv; e++) {
      Xov_[m][e] = pert_[m][e+no];
    }

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      Xoo_[m][i] = pert_[m][i];
      for(int e=0; e < nv; e++)
        Xoo_[m][i] += t1[i][e] * pert_[m][e+no];
    }

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      Xvv_[a][e] = pert_[a+no][e+no];
      for(int m=0; m < no; m++)
        Xvv_[a][e] -= t1[m][a] * pert_[m][e+no];
    }

  for(int a=0; a < nv; a++)
    for(int i=0; i < no; i++) {
      Xvo_[a][i] = pert_[a+no][i];
      for(int e=0; e < nv; e++)
        Xvo_[a][i] += t1[i][e] * pert_[a+no][e+no];
      for(int m=0; m < no; m++)
        Xvo_[a][i] -= t1[m][a] * pert_[m][i];
      for(int m=0; m < no; m++) 
        for(int e=0; e < nv; e++)
          Xvo_[a][i] += (2.0*t2[m][i][e][a] - t2[i][m][e][a] - t1[i][e] * t1[m][a]) * pert_[m][e+no];
    }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          Xovoo_[m][b][i][j] =0.0;
          for(int e=0; e < nv; e++)
            Xovoo_[m][b][i][j] += t2[i][j][e][b] * pert_[m][e+no];
        }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int i=0; i < no; i++) {
          Xvvvo_[a][b][e][i] = 0.0;
          for(int m=0; m < no; m++)
            Xvvvo_[a][b][e][i] -= t2[m][i][a][b] * pert_[m][e+no];
        }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          Xvvoo_[a][b][i][j] = 0.0;
          for(int e=0; e < nv; e++)
            Xvvoo_[a][b][i][j] += t2[i][j][e][b] * Xvv_[a][e];
          for(int m=0; m < no; m++)
            Xvvoo_[a][b][i][j] -= t2[m][j][a][b] * Xoo_[m][i];
        }
}

void CCPert::solve(enum hand myhand)
{
  outfile->Printf("\n\tCoupled-Cluster RH Perturbed Wfn Iteration:\n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t Iter   Pseudoresponse      RMS   \n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t  %3d  %20.15f\n", 0, pseudoresponse(myhand));

  double rms = 0.0;
  boost::shared_ptr<DIISManager> diis(new DIISManager(8, "CCPert DIIS",
    DIISManager::LargestError, DIISManager::InCore));
  diis->set_error_vector_size(2, DIISEntry::Pointer, no_*nv_, DIISEntry::Pointer, no_*no_*nv_*nv_);
  diis->set_vector_size(2, DIISEntry::Pointer, no_*nv_, DIISEntry::Pointer, no_*no_*nv_*nv_);
  for(int iter=1; iter <= CC_->maxiter_; iter++) {
    amp_save(myhand);

    if(myhand == right) { build_X1(); build_X2(); }
    else { build_Y1(); build_Y2(); }
    rms = increment_amps(myhand);

    if(rms < CC_->convergence_) break;
    if(CC_->do_diis_) {
      build_diis_error(myhand);
      if(myhand == right) {
        diis->add_entry(4, X1err_.data(), X2err_.data(), X1diis_.data(), X2diis_.data());
        if(diis->subspace_size() > 2) diis->extrapolate(2, X1diis_.data(), X2diis_.data());
      }
      else {
        diis->add_entry(4, Y1err_.data(), Y2err_.data(), Y1diis_.data(), Y2diis_.data());
        if(diis->subspace_size() > 2) diis->extrapolate(2, Y1diis_.data(), Y2diis_.data());
      }
      save_diis_vectors(myhand);
    }
    outfile->Printf(  "\t  %3d  %20.15f  %5.3e\n",iter, pseudoresponse(myhand), rms);
  }
  if(rms >= CC_->convergence_)
    throw PSIEXCEPTION("Computation has not converged.");
}

// Note that the doubles contribution expression assumes a lack of
// permutational symmetry of the Xvvoo_ quantity
double CCPert::pseudoresponse(enum hand myhand)
{
  int no = no_;
  int nv = nv_;
  double **X1, ****X2;

  if(myhand == right) { X1 = X1_;  X2 = X2_; }
  else {X1 = Y1_; X2 = Y2_; }

  double polar1 = 0.0;
  double polar2 = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      polar1 += 2.0 * X1[i][a] * Xvo_[a][i];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          polar2 += 2.0 * (2.0 * X2[i][j][a][b] - X2[i][j][b][a]) * Xvvoo_[a][b][i][j];
    }

  return -2.0*(polar1 + polar2);
}

void CCPert::build_diis_error(enum hand myhand)
{
  int no = no_;
  int nv = nv_;
  double **X1, ****X2, **X1old, ****X2old;
  double *X1diis, *X2diis, *X1err, *X2err;

  if(myhand == right) { 
    X1diis = X1diis_.data(); 
    X2diis = X2diis_.data(); 
    X1err = X1err_.data(); 
    X2err = X2err_.data();
    X1 = X1_; 
    X2 = X2_; 
    X1old = X1old_; 
    X2old = X2old_; 
  }
  else {
    X1diis = Y1diis_.data(); 
    X2diis = Y2diis_.data(); 
    X1err = Y1err_.data(); 
    X2err = Y2err_.data();
    X1 = Y1_; 
    X2 = Y2_; 
    X1old = Y1old_; 
    X2old = Y2old_; 
  }

  int X1len = 0;
  int X2len = 0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      X1diis[X1len] = X1[i][a];
      X1err[X1len++] = X1[i][a] - X1old[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          X2diis[X2len] = X2[i][j][a][b];
          X2err[X2len++] = X2[i][j][a][b] - X2old[i][j][a][b];
        }
    }
}

void CCPert::save_diis_vectors(enum hand myhand)
{
  int no = no_;
  int nv = nv_;

  double **X1, ****X2;
  double *X1diis, *X2diis;
 
  if(myhand == right) { 
    X1 = X1_; 
    X2 = X2_; 
    X1diis = X1diis_.data(); 
    X2diis = X2diis_.data(); 
  }
  else { 
    X1 = Y1_; 
    X2 = Y2_; 
    X1diis = Y1diis_.data(); 
    X2diis = Y2diis_.data(); 
  }

  int X1len = 0;
  int X2len = 0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      X1[i][a] = X1diis[X1len++];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          X2[i][j][a][b] = X2diis[X2len++];
        }
    }
}

void CCPert::amp_save(enum hand myhand)
{
  if(myhand == right) {
    double ****X2tmp = X2_;
    X2_ = X2old_;
    X2old_ = X2tmp;
    double **X1tmp = X1_;
    X1_ = X1old_;
    X1old_ = X1tmp;
  }
  else {
    double ****X2tmp = Y2_;
    Y2_ = Y2old_;
    Y2old_ = X2tmp;
    double **X1tmp = Y1_;
    Y1_ = Y1old_;
    Y1old_ = X1tmp;
  }
}

double CCPert::increment_amps(enum hand myhand)
{
  int no = no_;
  int nv = nv_;

  double **X1, **X1old;
  double ****X2, ****X2old;

  if(myhand == right) { X1 = X1_; X1old = X1old_; X2 = X2_; X2old = X2old_; }
  else { X1 = Y1_; X1old = Y1old_; X2 = Y2_; X2old = Y2old_; }

  double residual1 = 0.0;
  double residual2 = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      residual1 += X1[i][a] * X1[i][a];
      X1[i][a] = X1old[i][a] + X1[i][a]/(D1_[i][a] + omega_);
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          residual2 += X2[i][j][a][b] * X2[i][j][a][b];
          X2[i][j][a][b] = X2old[i][j][a][b] + X2[i][j][a][b]/(D2_[i][j][a][b] + omega_);
        }
    }

  return sqrt(residual1 + residual2);
}

void CCPert::build_X1()
{
  int no = no_;
  int nv = nv_;

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      X1_[i][a] = Xvo_[a][i] - omega_ * X1old_[i][a];
      for(int e=0; e < nv; e++) 
        X1_[i][a] += X1old_[i][e] * HBAR_->Hvv_[a][e];
      for(int m=0; m < no; m++)
        X1_[i][a] -= X1old_[m][a] * HBAR_->Hoo_[m][i];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          X1_[i][a] += X1old_[m][e] * (2.0*HBAR_->Hovvo_[m][a][e][i] - HBAR_->Hovov_[m][a][i][e]);
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          X1_[i][a] += HBAR_->Hov_[m][e] * (2.0*X2old_[m][i][e][a] - X2old_[i][m][e][a]);
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            X1_[i][a] += X2old_[i][m][e][f] * (2.0*HBAR_->Hvovv_[a][m][e][f] - HBAR_->Hvovv_[a][m][f][e]);
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            X1_[i][a] -= X2old_[m][n][a][e] * (2.0*HBAR_->Hooov_[m][n][i][e] - HBAR_->Hooov_[n][m][i][e]);
    }
}

void CCPert::build_X2()
{
  int no = no_;
  int nv = nv_;

  // Prep three-body terms
  double **Zvv = block_matrix(nv, nv);
  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      for(int m=0; m < no; m++)
        for(int f=0; f < nv; f++)
          Zvv[a][e] += (2.0 * HBAR_->Hvovv_[a][m][e][f] - HBAR_->Hvovv_[a][m][f][e]) * X1old_[m][f];
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int f=0; f < nv; f++)
            Zvv[a][e] -= H_->L_[m][n][e+no][f+no] * X2old_[m][n][a][f];
    }
  double **Zoo = block_matrix(no, no);
  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      for(int n=0; n < no; n++)
        for(int e=0; e < nv; e++)
          Zoo[m][i] -= (2.0 * HBAR_->Hooov_[m][n][i][e] - HBAR_->Hooov_[n][m][i][e]) * X1old_[n][e];
      for(int n=0; n < no; n++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            Zoo[m][i] -= H_->L_[m][n][e+no][f+no] * X2old_[i][n][e][f];
    }

  double ****R2 = init_4d_array(no, no, nv, nv); // residual
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          R2[i][j][a][b] = Xvvoo_[a][b][i][j] - 0.5 * omega_ * X2old_[i][j][a][b];
          for(int e=0; e < nv; e++)
            R2[i][j][a][b] += X1old_[i][e] * HBAR_->Hvvvo_[a][b][e][j];
          for(int m=0; m < no; m++)
            R2[i][j][a][b] -= X1old_[m][a] * HBAR_->Hovoo_[m][b][i][j];
          for(int m=0; m < no; m++)
            R2[i][j][a][b] += Zoo[m][i] * CC_->t2_[m][j][a][b];
          for(int e=0; e < nv; e++)
            R2[i][j][a][b] += Zvv[a][e] * CC_->t2_[i][j][e][b];
          for(int e=0; e < nv; e++)
            R2[i][j][a][b] += X2old_[i][j][e][b] * HBAR_->Hvv_[a][e];
          for(int m=0; m < no; m++)
            R2[i][j][a][b] -= X2old_[m][j][a][b] * HBAR_->Hoo_[m][i];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              R2[i][j][a][b] += 0.5 * X2old_[m][n][a][b] * HBAR_->Hoooo_[m][n][i][j];
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              R2[i][j][a][b] += 0.5 * X2old_[i][j][e][f] * HBAR_->Hvvvv_[a][b][e][f];
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              R2[i][j][a][b] -= X2old_[i][m][e][b] * HBAR_->Hovov_[m][a][j][e] + X2old_[i][m][e][a] * HBAR_->Hovvo_[m][b][e][j];
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              R2[i][j][a][b] += X2old_[m][i][e][a] * (2.0*HBAR_->Hovvo_[m][b][e][j] - HBAR_->Hovov_[m][b][j][e]);
        }

  // symmetrize and store
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          X2_[i][j][a][b] = (R2[i][j][a][b] + R2[j][i][b][a]);

  free_4d_array(R2, no, no, nv);
  free_block(Zoo);
  free_block(Zvv);
}

void CCPert::print_amps(enum hand myhand)
{
  int no = no_;
  int nv = nv_;

  double **X1, ****X2;
  if(myhand == right) { X1 = X1_; X2 = X2_; }
  else { X1 = Y1_; X2 = Y2_; }

  outfile->Printf("X1 Amplitudes:\n");
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      if(fabs(X1_[i][a]) > 1e-12)
        outfile->Printf("X1[%d][%d] = %20.15f\n", i, a, X1_[i][a]);

  outfile->Printf("X2 Amplitudes:\n");
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
 if(fabs(X2_[i][j][a][b]) > 1e-12)
 outfile->Printf("X2[%d][%d][%d][%d] = %20.15f\n", i, j, a, b, X2_[i][j][a][b]);

}

void CCPert::build_Y1()
{

}

void CCPert::build_Y2()
{

}

}} // psi::ugacc
