#include <string>
#include <cstdio>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

double triples_ooc(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double ***t3;

  // Scandinavian expression for (T) correction
  t3 = init_3d_array(nv, nv, nv);
  double **X1 = block_matrix(no, nv);
  double ****X2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {
        t3_ijk(double ***t3, int i, int j, int k, double ****t2, double **fock, double ****ints);

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++)
              X1[i][a] += (t3[a][b][c] - t3[c][b][a]) * L[j][k][b+no][c+no];

        for(int a=0; a < nv; a++) 
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {

              for(int l=0; l < no; l++) {
                X2[i][l][a][b] -= t3[a][b][c] * L[j][k][l][c+no];
                X2[i][l][a][b] += t3[c][b][a] * ints[j][k][l][c+no];
              }
              for(int d=0; d < nv; d++) {
                X2[i][j][a][d] += t3[a][b][c] * L[d+no][k][b+no][c+no];
                X2[i][j][a][d] -= t3[c][b][a] * ints[d+no][k][b+no][c+no];
              }
            } // abc
      } // ijk

  double ET_UGA = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      ET_UGA += t1s[i][a] * X1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          ET_UGA += t2s[i][j][a][b] * X2[i][j][a][b];
    }

  free_block(X1);
  free_4d_array(X2, no, no, nv);
  free_3d_array(t3, nv, nv);

  return ET_UGA;
}

void t3_ijk(double ***t3, int i, int j, int k, double ****t2, double **fock, double ****ints)
{
  int no = moinfo.no;
  int nv = moinfo.nv;

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

        t3[a][b][c] = value/denom;
      } // abc

  return;
}

void t3_abc(double ***t3, int a, int b, int c, double ****t2, double **fock, double ****ints)
{
  int no = moinfo.no;
  int nv = moinfo.nv;

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {
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

        t3[i][j][k] = value/denom;
      } // ijk

  return;
}

void l3_ijk(double ***l3, int i, int, j, int k, double ****t2s, double **t1s, double
**fock, double ****L, double ****ints)
{
  int moinfo.no;
  int moinfo.nv;

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++) {
        double value = 0.0;
  
        value += fock[i][a+no] * t2s[j][k][b][c] - fock[j][a+no] * t2s[i][k][b][c];
        value += fock[i][a+no] * t2s[k][j][c][b] - fock[k][a+no] * t2s[i][j][c][b];
        value += fock[j][b+no] * t2s[i][k][a][c] - fock[i][b+no] * t2s[j][k][a][c];
        value += fock[j][b+no] * t2s[k][i][c][a] - fock[k][b+no] * t2s[j][i][c][a];
        value += fock[k][c+no] * t2s[i][j][a][b] - fock[i][c+no] * t2s[k][j][a][b];
        value += fock[k][c+no] * t2s[j][i][b][a] - fock[j][c+no] * t2s[k][i][b][a];

        value += L[i][j][a+no][b+no] * t1s[k][c] - L[i][j][a+no][c+no] * t1s[k][b];
        value += L[i][k][a+no][c+no] * t1s[j][b] - L[i][k][a+no][b+no] * t1s[j][c];
        value += L[j][i][b+no][a+no] * t1s[k][c] - L[j][i][b+no][c+no] * t1s[k][a];
        value += L[j][k][b+no][c+no] * t1s[i][a] - L[j][k][b+no][a+no] * t1s[i][c];
        value += L[k][j][c+no][a+no] * t1s[j][b] - L[k][i][c+no][b+no] * t1s[j][a];
        value += L[k][i][c+no][b+no] * t1s[i][a] - L[k][j][c+no][a+no] * t1s[i][b];

        for(int f=0; f < nv; f++) {
          value += L[f+no][j][a+no][b+no] * t2s[i][k][f][c];
          value += L[f+no][k][a+no][c+no] * t2s[i][j][f][b];
          value += L[f+no][i][b+no][a+no] * t2s[j][k][f][c];
          value += L[f+no][k][b+no][c+no] * t2s[j][i][f][a];
          value += L[f+no][i][c+no][a+no] * t2s[k][j][f][b];
          value += L[f+no][j][c+no][b+no] * t2s[k][i][f][a];

          value -= ints[k][f+no][a+no][b+no] * t2s[i][j][c][f];
          value -= ints[j][f+no][a+no][c+no] * t2s[i][k][b][f];
          value -= ints[k][f+no][b+no][a+no] * t2s[j][i][c][f];
          value -= ints[i][f+no][b+no][c+no] * t2s[j][k][a][f];
          value -= ints[j][f+no][c+no][a+no] * t2s[k][i][b][f];
          value -= ints[i][f+no][c+no][b+no] * t2s[k][j][a][f];
        }
        for(int n=0; n < no; n++) {
          value -= L[i][j][a+no][n] * t2s[n][k][b][c];
          value -= L[i][k][a+no][n] * t2s[n][j][c][b];
          value -= L[j][i][b+no][n] * t2s[n][k][a][c];
          value -= L[j][k][b+no][n] * t2s[n][i][c][a];
          value -= L[k][i][c+no][n] * t2s[n][j][a][b];
          value -= L[k][j][c+no][n] * t2s[n][i][b][a];

          value += ints[k][j][a+no][n] * t2s[n][i][b][c];
          value += ints[j][k][a+no][n] * t2s[n][i][c][b];
          value += ints[k][i][b+no][n] * t2s[n][j][a][c];
          value += ints[i][k][b+no][n] * t2s[n][j][c][a];
          value += ints[j][i][c+no][n] * t2s[n][k][a][b];
          value += ints[i][j][c+no][n] * t2s[n][k][b][a];
        }

        double denom = fock[i][i] + fock[j][j] + fock[k][k];
        denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

        l3[a][b][c] = value/denom;
      } // abc
  return;
}

void l3_abc(double ***l3, int a, int, b, int c, double ****t2s, double **t1s, double **fock, double ****L, double ****ints)
{
  int moinfo.no;
  int moinfo.nv;

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {
        double value = 0.0;

        value += fock[i][a+no] * t2s[j][k][b][c] - fock[j][a+no] * t2s[i][k][b][c];
        value += fock[i][a+no] * t2s[k][j][c][b] - fock[k][a+no] * t2s[i][j][c][b];
        value += fock[j][b+no] * t2s[i][k][a][c] - fock[i][b+no] * t2s[j][k][a][c];
        value += fock[j][b+no] * t2s[k][i][c][a] - fock[k][b+no] * t2s[j][i][c][a];
        value += fock[k][c+no] * t2s[i][j][a][b] - fock[i][c+no] * t2s[k][j][a][b];
        value += fock[k][c+no] * t2s[j][i][b][a] - fock[j][c+no] * t2s[k][i][b][a];

        value += L[i][j][a+no][b+no] * t1s[k][c] - L[i][j][a+no][c+no] * t1s[k][b];
        value += L[i][k][a+no][c+no] * t1s[j][b] - L[i][k][a+no][b+no] * t1s[j][c];
        value += L[j][i][b+no][a+no] * t1s[k][c] - L[j][i][b+no][c+no] * t1s[k][a];
        value += L[j][k][b+no][c+no] * t1s[i][a] - L[j][k][b+no][a+no] * t1s[i][c];
        value += L[k][j][c+no][a+no] * t1s[j][b] - L[k][i][c+no][b+no] * t1s[j][a];
        value += L[k][i][c+no][b+no] * t1s[i][a] - L[k][j][c+no][a+no] * t1s[i][b];

        for(int f=0; f < nv; f++) {
          value += L[f+no][j][a+no][b+no] * t2s[i][k][f][c];
          value += L[f+no][k][a+no][c+no] * t2s[i][j][f][b];
          value += L[f+no][i][b+no][a+no] * t2s[j][k][f][c];
          value += L[f+no][k][b+no][c+no] * t2s[j][i][f][a];
          value += L[f+no][i][c+no][a+no] * t2s[k][j][f][b];
          value += L[f+no][j][c+no][b+no] * t2s[k][i][f][a];

          value -= ints[k][f+no][a+no][b+no] * t2s[i][j][c][f];
          value -= ints[j][f+no][a+no][c+no] * t2s[i][k][b][f];
          value -= ints[k][f+no][b+no][a+no] * t2s[j][i][c][f];
          value -= ints[i][f+no][b+no][c+no] * t2s[j][k][a][f];
          value -= ints[j][f+no][c+no][a+no] * t2s[k][i][b][f];
          value -= ints[i][f+no][c+no][b+no] * t2s[k][j][a][f];
        }
        for(int n=0; n < no; n++) {
          value -= L[i][j][a+no][n] * t2s[n][k][b][c];
          value -= L[i][k][a+no][n] * t2s[n][j][c][b];
          value -= L[j][i][b+no][n] * t2s[n][k][a][c];
          value -= L[j][k][b+no][n] * t2s[n][i][c][a];
          value -= L[k][i][c+no][n] * t2s[n][j][a][b];
          value -= L[k][j][c+no][n] * t2s[n][i][b][a];

          value += ints[k][j][a+no][n] * t2s[n][i][b][c];
          value += ints[j][k][a+no][n] * t2s[n][i][c][b];
          value += ints[k][i][b+no][n] * t2s[n][j][a][c];
          value += ints[i][k][b+no][n] * t2s[n][j][c][a];
          value += ints[j][i][c+no][n] * t2s[n][k][a][b];
          value += ints[i][j][c+no][n] * t2s[n][k][b][a];
        }

        double denom = fock[i][i] + fock[j][j] + fock[k][k];
        denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

        l3[i][j][k] = value/denom;
      } // ijk
  return;
}
}} // namespace psi::ugacc
