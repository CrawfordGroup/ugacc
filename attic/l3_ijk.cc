#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void l3_ijk(double ***l3, int i, int j, int k, double ****t2s, double **t1s, double **fock, double ****L, double ****ints)
{
  int no = moinfo.no;
  int nv = moinfo.nv;

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
        value += L[k][i][c+no][a+no] * t1s[j][b] - L[k][i][c+no][b+no] * t1s[j][a];
        value += L[k][j][c+no][b+no] * t1s[i][a] - L[k][j][c+no][a+no] * t1s[i][b];

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

}} // namespace psi::ugacc
