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

void W3_ijk(double ***, int, int, int, double ****, double ****);

double tcorr_ooc_TJL(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;

  double ***W3 = init_3d_array(nv, nv, nv);
  double ***V3 = init_3d_array(nv, nv, nv);
  double ***X3 = init_3d_array(nv, nv, nv);
  double ***Y3 = init_3d_array(nv, nv, nv);
  double ***Z3 = init_3d_array(nv, nv, nv);
  double ET = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j <= i; j++)
      for(int k=0; k <= j; k++) {
        W3_ijk(W3, i, j, k, t2, ints);

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++)
              V3[a][b][c] = (W3[a][b][c] + ints[i][j][a+no][b+no] * t1[k][c]
                          + ints[i][k][a+no][c+no] * t1[j][b]
                          + ints[j][k][b+no][c+no] * t1[i][a])/(1.0+(a==b)+(a==c)+(b==c));

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              X3[a][b][c] = W3[a][b][c] * V3[a][b][c] + W3[a][c][b] * V3[a][c][b] 
                          + W3[b][a][c] * V3[b][a][c] + W3[b][c][a] * V3[b][c][a]
                          + W3[c][a][b] * V3[c][a][b] + W3[c][b][a] * V3[c][b][a];
              Y3[a][b][c] = V3[a][b][c] + V3[b][c][a] + V3[c][a][b];
              Z3[a][b][c] = V3[a][c][b] + V3[b][a][c] + V3[c][b][a];
            }

        for(int a=0; a < nv; a++)
          for(int b=0; b <= a; b++)
            for(int c=0; c <= b; c++) {
              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              ET += ((Y3[a][b][c] - 2.0 * Z3[a][b][c]) * (W3[a][b][c] + W3[b][c][a] + W3[c][a][b])
                  + (Z3[a][b][c] - 2.0 * Y3[a][b][c]) * (W3[a][c][b] + W3[b][a][c] + W3[c][b][a])
                  + 3.0 * X3[a][b][c]) * (2.0 - ((i==j) + (i==k) + (j==k)))/denom;
            }

      } // ijk

  free_3d_array(W3, nv, nv);
  free_3d_array(V3, nv, nv);
  free_3d_array(X3, nv, nv);
  free_3d_array(Y3, nv, nv);
  free_3d_array(Z3, nv, nv);

  return ET;
}

}} // namespace psi::ugacc
