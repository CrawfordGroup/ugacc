#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void M3_ijk(double ***M3, int i, int j, int k, double ****t2, double **fock, double ****ints)
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

        M3[a][b][c] = value/denom;
      } // abc

  return;
}

}} // namespace psi::ugacc
