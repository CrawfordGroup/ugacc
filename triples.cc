#include <string>
#include <cstdio>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void amp_write(int, double **, double ****, std::string);
void triples_gradient(int i, int j, int k, double ***W, double ***V);
void triples_gradient_viking(double ******);

double triples(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double ******t3;

  // Scandinavian expression for (T) correction
  t3 = init_6d_array(no, no, no, nv, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              double value = 0.0;
              for(int e=0; e < nv; e++) {
                value +=
                   + ints[i][e+no][a+no][b+no] * t2[k][j][c][e]
                   + ints[i][e+no][a+no][c+no] * t2[j][k][b][e]
                   + ints[k][e+no][c+no][a+no] * t2[j][i][b][e]
                   + ints[k][e+no][c+no][b+no] * t2[i][j][a][e]
                   + ints[j][e+no][b+no][c+no] * t2[i][k][a][e]
                   + ints[j][e+no][b+no][a+no] * t2[k][i][c][e];
              }
              for(int m=0; m < no; m++) {
                value -=
                   + ints[j][k][m][c+no] * t2[i][m][a][b]
                   + ints[k][j][m][b+no] * t2[i][m][a][c]
                   + ints[i][j][m][b+no] * t2[k][m][c][a]
                   + ints[j][i][m][a+no] * t2[k][m][c][b]
                   + ints[k][i][m][a+no] * t2[j][m][b][c]
                   + ints[i][k][m][c+no] * t2[j][m][b][a];
              }
              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              t3[i][j][k][a][b][c] = value/denom;
            }
      }

  double **X1 = block_matrix(no, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++)
          for(int c=0; c < nv; c++)
            for(int d=0; d < nv; d++)
              X1[i][a] += (t3[i][k][l][a][c][d] - t3[l][k][i][a][c][d]) * L[k][l][c+no][d+no];

  double ****X2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              X2[i][j][a][b] += (t3[i][j][k][a][b][c] - t3[k][j][i][a][b][c]) * fock[k][c+no];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int l=0; l < no; l++) {
                X2[i][j][a][b] -= t3[i][k][l][a][b][c] * L[k][l][j][c+no];
                X2[i][j][a][b] += t3[l][k][i][a][b][c] * ints[k][l][j][c+no];
              }

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int d=0; d < nv; d++) {
                X2[i][j][a][b] += t3[i][j][k][a][c][d] * L[b+no][k][c+no][d+no];
                X2[i][j][a][b] -= t3[k][j][i][a][c][d] * ints[b+no][k][c+no][d+no];
              }
        }
  

  double ET_UGA = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      ET_UGA += t1[i][a] * X1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          ET_UGA += (2.0*t2[i][j][a][b] - t2[i][j][b][a]) * X2[i][j][a][b];    
    }
  ET_UGA *= 2.0;
  fprintf(outfile, "\tViking's (T)   = %20.14f\n", ET_UGA);

  if(params.dertype) triples_gradient_viking(t3);

  free_block(X1);
  free_4d_array(X2, no, no, nv);
  free_6d_array(t3, no, no, no, nv, nv);

  double ***W = init_3d_array(nv, nv, nv);
  double ***V = init_3d_array(nv, nv, nv);
  double ET = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              double value = 0.0;
              for(int e=0; e < nv; e++) {
                value +=
                   + ints[i][e+no][a+no][b+no] * t2[k][j][c][e]
                   + ints[i][e+no][a+no][c+no] * t2[j][k][b][e]
                   + ints[k][e+no][c+no][a+no] * t2[j][i][b][e]
                   + ints[k][e+no][c+no][b+no] * t2[i][j][a][e]
                   + ints[j][e+no][b+no][c+no] * t2[i][k][a][e]
                   + ints[j][e+no][b+no][a+no] * t2[k][i][c][e];
              }
              for(int m=0; m < no; m++) {
                value -=
                   + ints[j][k][m][c+no] * t2[i][m][a][b]
                   + ints[k][j][m][b+no] * t2[i][m][a][c]
                   + ints[i][j][m][b+no] * t2[k][m][c][a]
                   + ints[j][i][m][a+no] * t2[k][m][c][b]
                   + ints[k][i][m][a+no] * t2[j][m][b][c]
                   + ints[i][k][m][c+no] * t2[j][m][b][a];
              }

              W[a][b][c] = value;


              V[a][b][c] = ints[j][k][b+no][c+no] * t1[i][a]
                         + ints[i][k][a+no][c+no] * t1[j][b]
                         + ints[i][j][a+no][b+no] * t1[k][c];

            } /* abc loop */

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              ET += (W[a][b][c] + V[a][b][c] - W[c][b][a] - V[c][b][a])*
                    (4.0*W[a][b][c] + W[b][c][a] + W[c][a][b])/(3.0*denom);
            } /* abc loop */

//        triples_gradient(i, j, k, W, V);

      } /* ijk loop */

  free_3d_array(W, nv, nv);
  free_3d_array(V, nv, nv);

  return ET;
}

void triples_gradient_viking(double ******t3)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double **s1 = moinfo.s1;
  double ****s2 = moinfo.s2;

  // T3 --> Lambda1
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++)
          for(int c=0; c < nv; c++)
            for(int d=0; d < nv; d++)
              s1[i][a] += 2.0 * (t3[i][k][l][a][c][d] - t3[k][i][l][a][c][d]) * L[k][l][c+no][d+no];

  // T3 --> Lambda2
  double ****Y2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          double value = 0.0;         

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              value += (t3[i][j][k][a][b][c] - t3[i][k][j][a][b][c]) * fock[k][c+no];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int l=0; l < no; l++)
                value -= (2.0 * t3[i][l][k][a][b][c] - t3[k][l][i][a][b][c] - t3[i][k][l][a][b][c]) * ints[l][k][j][c+no];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int d=0; d < nv; d++)
                value += (2.0 * t3[i][k][j][a][c][d] - t3[k][i][j][a][c][d] - t3[i][j][k][a][c][d]) * ints[b+no][k][d+no][c+no];

          Y2[i][j][a][b] = value;
        }

  double ****X2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          X2[i][j][a][b] = 2.0 * (Y2[i][j][a][b] + Y2[i][j][b][a]);

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          s2[i][j][a][b] += 2.0 * X2[i][j][a][b] - X2[i][j][b][a];

  // Special left-projection T amplitudes
  double **t1s = block_matrix(no, nv);
  double ****t2s = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      t1s[i][a] = 2.0 * t1[i][a];

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          t2s[i][j][a][b] = 4.0 * t2[i][j][a][b] - 2.0 * t2[i][j][b][a];

  // Lambda3 amplitudes
  double ******l3_tmp = init_6d_array(no, no, no, nv, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              double value = 0.0;
  
              value += fock[i][a+no] * t2s[j][k][b][c];
              value -= fock[j][a+no] * t2s[i][k][b][c];

              value += L[i][j][a+no][b+no] * t1s[k][c];
              value -= L[i][j][a+no][c+no] * t1s[k][b];

              for(int f=0; f < nv; f++) {
                value += L[f+no][j][a+no][b+no] * t2s[i][k][f][c];
                value -= ints[k][f+no][a+no][b+no] * t2s[i][j][c][f];
              }
              for(int n=0; n < no; n++) {
                value -= L[i][j][a+no][n] * t2s[n][k][b][c];
                value += ints[k][j][a+no][n] * t2s[n][i][b][c];
              }

              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              l3_tmp[i][j][k][a][b][c] = value/denom;
            }
      }

  double ******l3 = init_6d_array(no, no, no, nv, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++)
              l3[i][j][k][a][b][c] = l3_tmp[i][j][k][a][b][c]
                                   + l3_tmp[i][k][j][a][c][b]
                                   + l3_tmp[j][i][k][b][a][c]
                                   + l3_tmp[j][k][i][b][c][a]
                                   + l3_tmp[k][i][j][c][a][b]
                                   + l3_tmp[k][j][i][c][b][a];
  free_6d_array(l3_tmp, no, no, no, nv, nv);

  // Lambda3 --> Lambda2
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          double value = 0.0;

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int l=0; l < no; l++) 
                value -= l3[i][k][l][a][c][b] * ints[c+no][j][k][l];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int d=0; d < nv; d++) 
                value += l3[i][k][j][a][c][d] * ints[c+no][d+no][k][b+no];

          Y2[i][j][a][b] = value;
        }

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
//          s2[i][j][a][b] += Y2[i][j][a][b] + Y2[i][j][b][a];
        }

  free_block(t1s);
  free_4d_array(t2s, no, no, nv);
  free_4d_array(X2, no, no, nv);
  free_4d_array(Y2, no, no, nv);
  free_6d_array(l3, no, no, no, nv, nv);

  // Print non-UGA version of these amps for comparison to spin-orbital code
  double **Z1 = block_matrix(no, nv);
  double ****Z2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      Z1[i][a] = 0.5 * moinfo.s1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          Z2[i][j][a][b] = (1./3.)*moinfo.s2[i][j][a][b] + 
                           (1./6.)*moinfo.s2[i][j][b][a];
   }
  amp_write(20, Z1, Z2, "SZ"); fprintf(outfile, "\n");
  free_block(Z1);
  free_4d_array(Z2, no, no, nv);

  return;
}

void triples_gradient(int i, int j, int k, double ***W, double ***V)
{
  int no = moinfo.no;
  int nv = moinfo.nv;

  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double **s1 = moinfo.s1;
  double ****s2 = moinfo.s2;
  double ***M = init_3d_array(nv, nv, nv);

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++) {

        double denom = fock[i][i] + fock[j][j] + fock[k][k];
        denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

        s1[i][a] += ints[j][k][b+no][c+no] *
           (4.0 * W[a][b][c] + W[b][c][a] + W[c][a][b]
          - 2.0 * W[c][b][a] - 2.0 * W[a][c][b]
          - 2.0 * W[b][a][c])/denom;

        M[a][b][c] = 6.0 * (
            8.0 * W[a][b][c] + 2.0 * W[c][a][b] + 2.0 * W[b][c][a]
          - 4.0 * W[a][c][b] - 4.0 * W[b][a][c] - 4.0 * W[c][b][a]
          + 4.0 * V[a][b][c] + 1.0 * V[c][a][b] + 1.0 * V[b][c][a]
          - 2.0 * V[a][c][b] - 2.0 * V[b][a][c] - 2.0 * V[c][b][a]
          )/denom;

//        for(int d=0; d < nv; d++)
//          s2[k][j][c][d] += M[a][b][c] * ints[i][b+no][a+no][d+no];
//        for(int l=0; l < no; l++)
//          s2[k][l][c][b] -= M[a][b][c] * ints[i][j][a+no][l];

      }

  free_3d_array(M, nv, nv);

  amp_write(20, moinfo.s1, moinfo.s2, "S"); fprintf(outfile, "\n");

  // Also print non-UGA version of these amps for comparison to
  // spin-orbital code
  double **Z1 = block_matrix(no, nv);
  double ****Z2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      Z1[i][a] = 0.5 * moinfo.s1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          Z2[i][j][a][b] = (1./3.)*moinfo.s2[i][j][a][b] + 
                           (1./6.)*moinfo.s2[i][j][b][a];
   }
  amp_write(20, Z1, Z2, "SZ"); fprintf(outfile, "\n");
  free_block(Z1);
  free_4d_array(Z2, no, no, nv);

  return;
}

}} // namespace psi::ugacc
