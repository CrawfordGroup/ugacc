#include <string>
#include <cstdio>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void t3_ijk(double ***, int, int, int, double ****, double **, double ****);

double tcorr_ooc(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double **t1s = moinfo.t1s;
  double ****t2s = moinfo.t2s;

  double ***t3 = init_3d_array(nv, nv, nv);
  double **X1 = block_matrix(no, nv);
  double ****X2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {
        t3_ijk(t3, i, j, k, t2, fock, ints);

        for(int a=0; a < nv; a++) 
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              X1[i][a] += (t3[a][b][c] - t3[c][b][a]) * L[j][k][b+no][c+no];

              X2[i][j][a][b] += (t3[a][b][c] - t3[c][b][a]) * fock[k][c+no];

              for(int l=0; l < no; l++)
                X2[i][l][a][b] -= (2.0*t3[a][b][c] - t3[a][c][b] - t3[c][b][a]) * ints[j][k][l][c+no];

              for(int d=0; d < nv; d++)
                X2[i][j][a][d] += (2.0*t3[a][b][c] - t3[a][c][b] - t3[c][b][a]) * ints[d+no][k][b+no][c+no];
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

}} // namespace psi::ugacc
