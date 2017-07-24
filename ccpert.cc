#include "ccpert.h"

#include <psi4/libciomr/libciomr.h>
#include <psi4/libqt/qt.h>
#include <psi4/libdiis/diismanager.h>
#include "array.h"

namespace psi { namespace ugacc {

CCPert::CCPert(double **pert, double omega, shared_ptr<CCWfn> CC, shared_ptr<HBAR> HBAR, shared_ptr<CCLambda> CCLambda)
{
  CC_ = CC;
  HBAR_ = HBAR;
  pert_ = pert;
  omega_ = omega;
  CCLambda_= CCLambda;


  H_ = CC_->H_;
  no_ = CC_->no_;
  nv_ = CC_->nv_;
  l1_ = CCLambda_->l1_;
  l2_ = CCLambda_->l2_;

  int no = no_;
  int nv = nv_;

  Aov_ = block_matrix(no, nv);
  Aoo_ = block_matrix(no, no);
  Avv_ = block_matrix(nv, nv);
  Avo_ = block_matrix(nv, no);
  Aovoo_ = init_4d_array(no, nv, no, no);
  Avvvo_ = init_4d_array(nv, nv, nv, no);
  Avvoo_ = init_4d_array(nv, nv, no, no);

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
      X1_[i][a] = Avo_[a][i]/(D1_[i][a] + omega_);
      Y1_[i][a] = Avo_[a][i]/(D1_[i][a] + omega_);
    }

  X2_ = init_4d_array(no, no, nv, nv);
  X2old_ = init_4d_array(no, no, nv, nv);
  Y2_ = init_4d_array(no, no, nv, nv);
  Y2old_ = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          X2_[i][j][a][b] = (Avvoo_[a][b][i][j]+Avvoo_[b][a][j][i])/(D2_[i][j][a][b] + omega_);
          Y2_[i][j][a][b] = (Avvoo_[a][b][i][j]+Avvoo_[b][a][j][i])/(D2_[i][j][a][b] + omega_);
        }

  Gvv_ = block_matrix(nv, nv);
  Goo_ = block_matrix(no, no);

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

  free_block(Aov_);
  free_block(Aoo_);
  free_block(Avv_);
  free_block(Avo_);
  free_4d_array(Aovoo_, no, nv, no);
  free_4d_array(Avvvo_, nv, nv, nv);
  free_4d_array(Avvoo_, nv, nv, no);

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

  free_block(Gvv_);
  free_block(Goo_);

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
      Aov_[m][e] = pert_[m][e+no];
    }

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      Aoo_[m][i] = pert_[m][i];
      for(int e=0; e < nv; e++)
        Aoo_[m][i] += t1[i][e] * pert_[m][e+no];
    }

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      Avv_[a][e] = pert_[a+no][e+no];
      for(int m=0; m < no; m++)
        Avv_[a][e] -= t1[m][a] * pert_[m][e+no];
    }

  for(int a=0; a < nv; a++)
    for(int i=0; i < no; i++) {
      Avo_[a][i] = pert_[a+no][i];
      for(int e=0; e < nv; e++)
        Avo_[a][i] += t1[i][e] * pert_[a+no][e+no];
      for(int m=0; m < no; m++)
        Avo_[a][i] -= t1[m][a] * pert_[m][i];
      for(int m=0; m < no; m++) 
        for(int e=0; e < nv; e++)
          Avo_[a][i] += (2.0*t2[m][i][e][a] - t2[i][m][e][a] - t1[i][e] * t1[m][a]) * pert_[m][e+no];
    }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          Aovoo_[m][b][i][j] =0.0;
          for(int e=0; e < nv; e++)
            Aovoo_[m][b][i][j] += t2[i][j][e][b] * pert_[m][e+no];
        }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int i=0; i < no; i++) {
          Avvvo_[a][b][e][i] = 0.0;
          for(int m=0; m < no; m++)
            Avvvo_[a][b][e][i] -= t2[m][i][a][b] * pert_[m][e+no];
        }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          Avvoo_[a][b][i][j] = 0.0;
          for(int e=0; e < nv; e++)
            Avvoo_[a][b][i][j] += t2[i][j][e][b] * Avv_[a][e];
          for(int m=0; m < no; m++)
            Avvoo_[a][b][i][j] -= t2[m][j][a][b] * Aoo_[m][i];
        }
}

void CCPert::solve(enum hand myhand)
{
  if (myhand == right)
  outfile->Printf("\n\tCoupled-Cluster RH Perturbed Wfn Iteration:\n");
  else
  outfile->Printf("\n\tCoupled-Cluster LH Perturbed Wfn Iteration:\n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t Iter   Pseudoresponse      RMS   \n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t  %3d  %20.15f\n", 0, pseudoresponse(myhand));

  double rms = 0.0;
  shared_ptr<DIISManager> diis(new DIISManager(8, "CCPert DIIS",
    DIISManager::LargestError, DIISManager::InCore));
  diis->set_error_vector_size(2, DIISEntry::Pointer, no_*nv_, DIISEntry::Pointer, no_*no_*nv_*nv_);
  diis->set_vector_size(2, DIISEntry::Pointer, no_*nv_, DIISEntry::Pointer, no_*no_*nv_*nv_);
  for(int iter=1; iter <= CC_->maxiter_; iter++) {
    amp_save(myhand);

    if(myhand == right) { build_X1(); build_X2(); }
    else { build_G(); build_Y1(); build_Y2(); }
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
      polar1 += 2.0 * X1[i][a] * Avo_[a][i];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          polar2 += 2.0 * (2.0 * X2[i][j][a][b] - X2[i][j][b][a]) * Avvoo_[a][b][i][j];
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
      X1_[i][a] = Avo_[a][i] - omega_ * X1old_[i][a];
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
          R2[i][j][a][b] = Avvoo_[a][b][i][j] - 0.5 * omega_ * X2old_[i][j][a][b];
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

void CCPert::build_G()
{
  int no = no_;
  int nv = nv_;

  double ****t2 = CC_->t2_;
  double ****Y2 = Y2old_;

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = 0.0;
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          for(int j=0; j < no; j++)
            value += t2[m][j][a][b] * Y2[i][j][a][b];
      Goo_[m][i] = value;
    }

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = 0.0;
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++)
          for(int b=0; b < nv; b++)
            value -= t2[i][j][e][b] * Y2[i][j][a][b];
      Gvv_[a][e] = value;
    }
}


void CCPert::build_Y1()
{
  int no = no_;
  int nv = nv_;
  double **Y1 = Y1old_;
  double ****Y2 = Y2old_;
  double **Y1new = Y1_;

  double **Gvv = Gvv_;
  double **Goo = Goo_;
  double **Hvv = HBAR_->Hvv_;
  double **Hoo = HBAR_->Hoo_;
  double **Hov = HBAR_->Hov_;
  double ****Hvvvo = HBAR_->Hvvvo_;
  double ****Hovoo = HBAR_->Hovoo_;
  double ****Hovvo = HBAR_->Hovvo_;
  double ****Hovov = HBAR_->Hovov_;
  double ****Hvovv = HBAR_->Hvovv_;
  double ****Hooov = HBAR_->Hooov_;
  double ****Hvvvv = HBAR_->Hvvvv_;
  double ****Hoooo = HBAR_->Hoooo_;

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      Y1_[i][a] = 2*pert_[i][a+no] + omega_ * Y1[i][a];

    for(int e=0; e < nv; e++)
      Y1_[i][a] += Y1[i][e] * Hvv[e][a];

    for(int m=0; m < no; m++)
      Y1_[i][a] -= Y1[m][a] * Hoo[i][m];

    for(int m=0; m < no; m++)
      for(int e=0; e < nv; e++)
        Y1_[i][a] += Y1[m][e] * (2 * Hovvo[i][e][a][m] - Hovov[i][e][m][a]);

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            Y1_[i][a] += Y2[i][m][e][f] * Hvvvo[e][f][a][m];

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            Y1_[i][a] -= Y2[m][n][a][e] * Hovoo[i][e][m][n];

      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++)
          Y1_[i][a] -= Gvv[e][f] * (2*Hvovv[e][i][f][a] - Hvovv[e][i][a][f]);

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          Y1_[i][a] -= Goo[m][n] * (2*Hooov[m][i][n][a] - Hooov[i][m][n][a]);
////
//Extra terms below 
////
     
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
	      Y1_[i][a] += 2.0 * X1_[m][e] * H_->L_[i][m][a+no][e+no];   // (i)

      for(int m=0; m < no; m++)
          Y1_[i][a] -= l1_[m][a] * Aoo_[i][m];  // (m)

      for(int e=0; e < nv; e++)
          Y1_[i][a] += l1_[i][e] * Avv_[e][a];  // (m)

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            Y1_[i][a] += l2_[i][m][f][e] * Avvvo_[f][e][a][m];     // (n)

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++){
            Y1_[i][a] -= 0.5 * l2_[m][n][e][a] * Aovoo_[i][e][n][m];    //  (n)
            Y1_[i][a] -= 0.5 * l2_[m][n][a][e] * Aovoo_[i][e][m][n];    //  (n)
	      }

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++){
	      Y1_[i][a] -= X1_[m][e] * Hov[m][a] * l1_[i][e];  // (q)	
	      Y1_[i][a] -= X1_[m][e] * Hov[i][e] * l1_[m][a];  // (q)	
          for(int n=0; n < no; n++){
            Y1_[i][a] -= X1_[m][e] * (2*Hooov[m][i][n][a] - Hooov[i][m][n][a]) * l1_[n][e];  // (q)  
            Y1_[i][a] -= X1_[m][e] * (2*Hooov[i][m][n][e] - Hooov[m][i][n][e]) * l1_[n][a];  // (q) 
		  }
          for(int f=0; f < nv; f++){
            Y1_[i][a] += X1_[m][e] * (2*Hvovv[f][m][a][e] - Hvovv[f][m][e][a]) * l1_[i][f];   // (q)
            Y1_[i][a] += X1_[m][e] * (2*Hvovv[f][i][e][a] - Hvovv[f][i][a][e]) * l1_[m][f];   // (q) 
          }
        }
  
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++){
       	      Y1_[i][a] += 2.0 * X2_[m][n][e][f] * H_->L_[i][m][a+no][e+no] * l1_[n][f];   // (r)
       	      Y1_[i][a] -= X2_[m][n][e][f] * H_->L_[i][m][a+no][f+no] * l1_[n][e];         // (r)
       	      Y1_[i][a] -= X2_[m][n][e][f] * H_->L_[m][i][e+no][f+no] * l1_[n][a];         // (r)
       	      Y1_[i][a] -= X2_[m][n][e][f] * H_->L_[n][m][f+no][a+no] * l1_[i][e];         // (r)
  			}

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++){
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++){
              Y1_[i][a] -= X1_[m][e] * Hovov[m][f][n][a] * l2_[n][i][e][f]; 	// (s)
	          Y1_[i][a] -= X1_[m][e] * Hovov[i][f][n][e] * l2_[n][m][a][f]; 	// (s)
	          Y1_[i][a] -= X1_[m][e] * Hovvo[m][f][a][n] * l2_[n][i][f][e];     // (s)	
	          Y1_[i][a] -= X1_[m][e] * Hovvo[i][f][e][n] * l2_[n][m][f][a];     // (s)
		    }

	    for(int f=0; f < nv; f++)
          for(int g=0; g < nv; g++){
            Y1_[i][a] += 0.5 * X1_[m][e] * Hvvvv[f][g][a][e] * l2_[i][m][f][g];   // (s) ---problem
            Y1_[i][a] += 0.5 * X1_[m][e] * Hvvvv[f][g][e][a] * l2_[m][i][f][g];   // (s) ---problem
          }
        for(int n=0; n < no; n++)
          for(int o=0; o < no; o++){
            Y1_[i][a] += 0.5 * X1_[m][e] * Hoooo[i][m][n][o] * l2_[n][o][a][e];   // (s)  ---problem
            Y1_[i][a] += 0.5 * X1_[m][e] * Hoooo[m][i][n][o] * l2_[n][o][e][a];   // (s)  ---problem
          }
        }    

// //double **Zvv = block_matrix(nv,nv);
 double ****t2 = CC_->t2_;
////
////       for(int f=0; f < nv; f++)
////          for(int b=0; b < nv; b++)
////    	     for(int j=0; j < no; j++)
////                 for(int m=0; m < no; m++)
////                     for(int e=0; e < nv; e++)
////                       Zvv[f][b] += t2[j][m][f][e] * l2_[j][m][b][e];
////
////        double **Zoo = block_matrix(no,no);
////
////       for(int n=0; n < no; n++)
////          for(int j=0; j < no; j++)
////             for(int b=0; b < nv; b++)
////               for(int f=0; f < nv; f++) {
////                          Zoo[n][i] += t2[j][n][b][f] * l2_[j][i][b][f];
////        }
////
////       for(int n=0; n < no; n++)
////          for(int b=0; b < nv; b++)
////             for(int f=0; f < nv; f++){
////              Y1_[i][a] -= 0.5 * H_->L_[i][n][a+no][f+no] * Zvv[f][b] * X1_[n][b];
////              Y1_[i][a] -= 0.5 * H_->L_[n][i][b+no][f+no] * Zvv[f][a] * X1_[n][b];
////        }
////
////        for(int m=0; m < no; m++)
////           for(int n=0; n < no; n++)
////              for(int e=0; e < nv; e++){
////                 Y1_[i][a] -= 0.5 * H_->L_[m][n][e+no][a+no] * Zoo[n][i] * X1_[m][e];
////                 Y1_[i][a] -= 0.5 * H_->L_[i][n][a+no][e+no] * Zoo[n][m] * X1_[m][e];
////        }
////

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              for(int j=0; j < no; j++)
                 for(int b=0; b < nv; b++){
                    Y1_[i][a] -=  H_->L_[i][n][a+no][f+no] * X1_[n][b] * t2[j][m][f][e] * l2_[j][m][b][e] ;
                    Y1_[i][a] -=  H_->L_[m][i][e+no][f+no] * X1_[m][e] * t2[j][n][f][b] * l2_[j][n][a][b] ;
		            Y1_[i][a] -=  H_->L_[m][n][e+no][a+no] * X1_[m][e] * t2[j][n][f][b] * l2_[j][i][f][b] ;
		            Y1_[i][a] -=  H_->L_[i][n][a+no][f+no] * X1_[j][f] * t2[n][m][e][b] * l2_[j][m][e][b] ;
                 }





      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++){
              Y1_[i][a] -=  X2_[m][n][e][f] * Hov[m][a] * l2_[n][i][f][e];   // (t)
	          Y1_[i][a] -=  X2_[m][n][e][f] * Hov[i][e] * l2_[n][m][f][a];   // (t)
              for(int g=0; g < nv; g++){
	            Y1_[i][a] -=  X2_[m][n][e][f] * Hvovv[g][n][e][a] * l2_[i][m][f][g] ; 	//(t)
	            Y1_[i][a] -=  X2_[m][n][e][f] * Hvovv[g][n][a][e] * l2_[m][i][f][g] ; 	//(t)
	            Y1_[i][a] -=  X2_[m][n][e][f] * Hvovv[g][i][e][f] * l2_[n][m][a][g] ;      //(t) 
	            Y1_[i][a] +=  X2_[m][n][e][f] * (2*Hvovv[g][m][a][e] - Hvovv[g][m][e][a]) * l2_[n][i][f][g] ; 	//(t)
	            Y1_[i][a] +=  X2_[m][n][e][f] * (2*Hvovv[g][i][e][a] - Hvovv[g][i][a][e]) * l2_[n][m][f][g] ;      //(t)
	          }
              for(int o=0; o < no; o++){
	            Y1_[i][a] +=  X2_[m][n][e][f] * Hooov[m][n][o][a] * l2_[o][i][e][f] ; // (t)	
	            Y1_[i][a] +=  X2_[m][n][e][f] * Hooov[i][n][o][e] * l2_[o][m][a][f] ; // (t) 	
	            Y1_[i][a] +=  X2_[m][n][e][f] * Hooov[m][i][o][f] * l2_[o][n][e][a] ; // (t)	
	            Y1_[i][a] -=  X2_[m][n][e][f] * (2*Hooov[m][i][o][a] - Hooov[i][m][o][a]) * l2_[n][o][f][e] ; 	// (t)
	            Y1_[i][a] -=  X2_[m][n][e][f] * (2*Hooov[i][m][o][e] - Hooov[m][i][o][e]) * l2_[n][o][f][a] ;      // (t)
              }
            }
     
  }
}

void CCPert::build_Y2()
{
  int no = no_;
  int nv = nv_;
  double **Y1 = Y1old_;
  double ****Y2 = Y2old_;
  double ****Y2new = Y2;

  double **Gvv = Gvv_;
  double **Goo = Goo_;
  double **Hvv = HBAR_->Hvv_;
  double **Hoo = HBAR_->Hoo_;
  double **Hov = HBAR_->Hov_;
  double ****Hvvvo = HBAR_->Hvvvo_;
  double ****Hovoo = HBAR_->Hovoo_;
  double ****Hovvo = HBAR_->Hovvo_;
  double ****Hovov = HBAR_->Hovov_;
  double ****Hvvvv = HBAR_->Hvvvv_;
  double ****Hoooo = HBAR_->Hoooo_;
  double ****Hvovv = HBAR_->Hvovv_;
  double ****Hooov = HBAR_->Hooov_;
  double ****ints =  H_->ints_;

  double ****R2 = init_4d_array(no, no, nv, nv); // residual


  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          R2[i][j][a][b] = 0.5 * omega_ * Y2[i][j][a][b]; // 1/2 because we'll be permuting (ia,jb)

          R2[i][j][a][b] += 2.0*Y1[i][a]*Hov[j][b] - Y1[j][a]*Hov[i][b];

          for(int e=0; e < nv; e++)
            R2[i][j][a][b] += Y2[i][j][e][b]*Hvv[e][a];

          for(int m=0; m < no; m++)
            R2[i][j][a][b] -= Y2[m][j][a][b]*Hoo[i][m];

          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              R2[i][j][a][b] += 0.5 * Y2[m][n][a][b] * Hoooo[i][j][m][n];

          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              R2[i][j][a][b] += 0.5 * Y2[i][j][e][f] * Hvvvv[e][f][a][b];

          for(int e=0; e < nv; e++)
              R2[i][j][a][b] += Y1[i][e]*(2*Hvovv[e][j][a][b] - Hvovv[e][j][b][a]);

          for(int m=0; m < no; m++)
              R2[i][j][a][b] -= Y1[m][b]*(2*Hooov[j][i][m][a] - Hooov[i][j][m][a]);

          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++) {
              R2[i][j][a][b] += Y2[m][j][e][b] * (2*Hovvo[i][e][a][m] - Hovov[i][e][m][a]);
              R2[i][j][a][b] -= Y2[m][i][b][e] * Hovov[j][e][m][a];
              R2[i][j][a][b] -= Y2[m][i][e][b] * Hovvo[j][e][a][m];
            }

          for(int e=0; e < nv; e++)
            R2[i][j][a][b] += Gvv[a][e]*H_->L_[i][j][e+no][b+no];
          for(int m=0; m < no; m++)
            R2[i][j][a][b] -= Goo[m][i]*H_->L_[m][j][a+no][b+no];

        // Additional terms : AK
	  
	      R2[i][j][a][b] += 2*l1_[j][b] * pert_[i][a+no] - l1_[i][b] * pert_[j][a+no];  // (o)

 	      for(int e=0; e < nv; e++)
            R2[i][j][a][b] += l2_[i][j][e][b]*Avv_[e][a];  // (p)

          for(int m=0; m < no; m++)
            R2[i][j][a][b] -= l2_[m][j][a][b]*Aoo_[i][m];  // (p)

          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++) {
	          R2[i][j][a][b] -=  X1_[m][e] * l1_[j][a] * H_->L_[m][i][e+no][b+no]; // (u)
              R2[i][j][a][b] -=  X1_[m][e] * l1_[m][b] * H_->L_[i][j][a+no][e+no]; // (u)
	          R2[i][j][a][b] -=  X1_[m][e] * l1_[i][e] * H_->L_[j][m][b+no][a+no]; // (u)
	          R2[i][j][a][b] += 2.0 * X1_[m][e] * l1_[j][b] * H_->L_[i][m][a+no][e+no]; // (u) 
            }

	      for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++) {
	          R2[i][j][a][b] -=  Hov[m][a] * l2_[j][i][b][e] * X1_[m][e]; // (w)
	          R2[i][j][a][b] -=  Hov[i][e] * l2_[j][m][b][a] * X1_[m][e]; // (w)
              for(int f=0; f < nv; f++){
                R2[i][j][a][b] -=  Hvovv[f][m][b][a] * l2_[j][i][f][e] * X1_[m][e];   // (w)
                R2[i][j][a][b] -=  Hvovv[f][j][e][a] * l2_[m][i][f][b] * X1_[m][e];   // (w)
                R2[i][j][a][b] -=  Hvovv[f][i][b][e] * l2_[j][m][f][a] * X1_[m][e];   // (w)
                R2[i][j][a][b] +=  (2 * Hvovv[f][m][a][e] - Hvovv[f][m][e][a]) * l2_[j][i][b][f] * X1_[m][e];  // (w)
                R2[i][j][a][b] +=  (2 * Hvovv[f][i][e][a] - Hvovv[f][i][a][e]) * l2_[j][m][b][f] * X1_[m][e];  // (w)
              }
	          for(int n=0; n < no; n++){
                R2[i][j][a][b] += Hooov[j][m][n][a] * l2_[n][i][b][e] * X1_[m][e];   // (w)
                R2[i][j][a][b] += Hooov[m][j][n][a] * l2_[n][i][e][b] * X1_[m][e];   // (w)
                R2[i][j][a][b] += Hooov[j][i][n][e] * l2_[n][m][b][a] * X1_[m][e];   // (w)
                R2[i][j][a][b] -= (2 * Hooov[m][i][n][a] - Hooov[i][m][n][a]) * l2_[j][n][b][e] * X1_[m][e];  // (w)
                R2[i][j][a][b] -= (2 * Hooov[i][m][n][e] - Hooov[m][i][n][e]) * l2_[j][n][b][a] * X1_[m][e];  // (w)
              }
	        }

	      for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                for(int f=0; f < nv; f++){					

	      	         /* x terms */

	              R2[i][j][a][b] += 0.25 * ints[n][m][a+no][b+no] * l2_[j][i][e][f] * X2_[m][n][e][f] ;
	              R2[i][j][a][b] += 0.25 * ints[m][n][a+no][b+no] * l2_[j][i][f][e] * X2_[m][n][e][f] ;

	              R2[i][j][a][b] += 0.25 * ints[j][n][a+no][e+no] * l2_[m][i][f][b] * X2_[m][n][e][f] ; 
	              R2[i][j][a][b] += 0.25 * ints[j][m][a+no][f+no] * l2_[n][i][e][b] * X2_[m][n][e][f] ; 

	              R2[i][j][a][b] += 0.25 * ints[n][j][a+no][e+no] * l2_[m][i][b][f] * X2_[m][n][e][f] ; 
	              R2[i][j][a][b] += 0.25 * ints[m][j][a+no][f+no] * l2_[n][i][b][e] * X2_[m][n][e][f] ; 

	              R2[i][j][a][b] += 0.25 * ints[j][n][e+no][a+no] * l2_[i][m][f][b] * X2_[m][n][e][f] ; 
	              R2[i][j][a][b] += 0.25 * ints[j][m][f+no][a+no] * l2_[i][n][e][b] * X2_[m][n][e][f] ; 

	              R2[i][j][a][b] += 0.25 * ints[n][j][e+no][a+no] * l2_[i][m][b][f] * X2_[m][n][e][f] ; 
	              R2[i][j][a][b] += 0.25 * ints[m][j][f+no][a+no] * l2_[i][n][b][e] * X2_[m][n][e][f] ; 

	              R2[i][j][a][b] += 0.25 * ints[j][i][e+no][f+no] * l2_[n][m][a][b] * X2_[m][n][e][f] ;
	              R2[i][j][a][b] += 0.25 * ints[j][i][f+no][e+no] * l2_[m][n][a][b] * X2_[m][n][e][f] ;

 	              R2[i][j][a][b] -=  0.5 * H_->L_[i][m][a+no][f+no] * l2_[j][n][b][e] * X2_[m][n][e][f] ; 
 	              R2[i][j][a][b] -=  0.5 * H_->L_[i][n][a+no][e+no] * l2_[j][m][b][f] * X2_[m][n][e][f] ; 
 	              R2[i][j][a][b] -=  0.5 * H_->L_[m][i][e+no][f+no] * l2_[j][n][b][a] * X2_[m][n][e][f] ; 
 	              R2[i][j][a][b] -=  0.5 * H_->L_[m][n][e+no][a+no] * l2_[j][i][b][f] * X2_[m][n][e][f] ; 
 	              R2[i][j][a][b] -=  0.5 * H_->L_[n][m][f+no][a+no] * l2_[j][i][b][e] * X2_[m][n][e][f] ; 
 	              R2[i][j][a][b] -=  0.5 * H_->L_[n][i][f+no][e+no] * l2_[j][m][b][a] * X2_[m][n][e][f] ; 

	              R2[i][j][a][b] -=   H_->L_[i][j][a+no][e+no] * l2_[n][m][f][b] * X2_[m][n][e][f] ;  
	              R2[i][j][a][b] -=   H_->L_[i][m][a+no][b+no] * l2_[n][j][f][e] * X2_[m][n][e][f] ; 
	              R2[i][j][a][b] -=   H_->L_[m][j][e+no][a+no] * l2_[n][i][f][b] * X2_[m][n][e][f] ;

	              R2[i][j][a][b] += 2.0 * H_->L_[i][m][a+no][e+no] * l2_[n][j][f][b] * X2_[m][n][e][f] ; 

                }
   }      

     for(int i=0; i < no; i++)
         for(int j=0; j < no; j++)
           for(int a=0; a < nv; a++)
             for(int b=0; b < nv; b++) 
                Y2_[i][j][a][b] = R2[i][j][a][b] + R2[j][i][b][a] ;

}


}} // psi::ugacc
