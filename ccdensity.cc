#include <boost/shared_ptr.hpp>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <cmath>
#include <libpsio/psio.h>

#include "ccdensity.h"
#include "ccrhwavefunction.h"
#include "cclhwavefunction.h"
#include "globals.h"

namespace psi { namespace ugacc {

CCDensity::CCDensity(boost::shared_ptr<CCRHWavefunction> CC, boost::shared_ptr<CCLHWavefunction> Lambda)
{
  no_ = CC->no_;
  nv_ = CC->nv_;

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

CCDensity::~CCDensity()
{
  int no = no_;
  int nv = nv_;

  free_block(Doo_);
  free_block(Dvv_);
  free_block(Dov_);
  free_block(Dvo_);
  free_4d_array(Goooo_, no, no, no);
  free_4d_array(Gvvvv_, nv, nv, nv);
  free_4d_array(Gooov_, no, no, no);
  free_4d_array(Gvvvo_, nv, nv, nv);
  free_4d_array(Govov_, no, nv, no);
  free_4d_array(Goovv_, no, no, nv);
}

double CCDensity::onepdm()
{
  int no = no_;
  int nv = nv_;
  double **t1 = t1_;
  double ****t2 = t2_;
  double ****tau = tau_;
  double **l1 = l1_;
  double ****l2 = l2_;
  double **fock = H_->fock_;

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

double CCDensity::twopdm()
{
  int no = no_;
  int nv = nv_;
  double **t1 = t1_;
  double ****t2 = t2_;
  double ****tau = tau_;
  double **l1 = l1_;
  double ****l2 = l2_;
  double ****ints = H_->ints_;

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
