#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

double triples(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;

  double ***W = init_3d_array(nv, nv, nv);
  double ***V = init_3d_array(nv, nv, nv);
  double ***X = init_3d_array(nv, nv, nv);
  double ***Y = init_3d_array(nv, nv, nv);
  double ***Z = init_3d_array(nv, nv, nv);

  double ET = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j <= i; j++)
      for(int k=0; k <= j; k++) {

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


              V[a][b][c] = value
                   + ints[j][k][b+no][c+no] * t1[i][a]
                   + ints[i][k][a+no][c+no] * t1[j][b]
                   + ints[i][j][a+no][b+no] * t1[k][c];

              V[a][b][c] /= (1 + (a==b) + (b==c) + (a==c));

            } /* abc loop */

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              X[a][b][c] = 
                 W[a][b][c] * V[a][b][c] + W[a][c][b] * V[a][c][b]
               + W[b][a][c] * V[b][a][c] + W[b][c][a] * V[b][c][a]
               + W[c][a][b] * V[c][a][b] + W[c][b][a] * V[c][b][a];

              Y[a][b][c] = V[a][b][c] + V[b][c][a] + V[c][a][b];
              Z[a][b][c] = V[a][c][b] + V[b][a][c] + V[c][b][a];
  
            } /* abc loop */

        for(int a=0; a < nv; a++)
          for(int b=0; b <= a; b++)
            for(int c=0; c <= b; c++) {
              double value1 = Y[a][b][c] - 2 * Z[a][b][c];
              double value2 = Z[a][b][c] - 2 * Y[a][b][c];
              double value3 = W[a][b][c] + W[b][c][a] + W[c][a][b];
              double value4 = W[a][c][b] + W[b][a][c] + W[c][b][a];
              double value5 = 3 * X[a][b][c];
              double value6 = 2 - ((i==j) + (j==k) + (i==k));

              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              ET += (value1*value3 + value2*value4 + value5)*value6/denom;
            } /* abc loop */
    
      } /* ijk loop */

  free_3d_array(W, nv, nv);
  free_3d_array(V, nv, nv);
  free_3d_array(X, nv, nv);
  free_3d_array(Y, nv, nv);
  free_3d_array(Z, nv, nv);

  return ET;
}

}} // namespace psi::ugacc
