#include <math.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

double increment_t(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **t1 = moinfo.t1;
  double **t1old = moinfo.t1old;
  double ****t2 = moinfo.t2;
  double ****t2old = moinfo.t2old;
  double **D1 = moinfo.D1;
  double ****D2 = moinfo.D2;

  double residual1 = 0.0;
  double residual2 = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      residual1 += t1[i][a] * t1[i][a];     
      t1[i][a] = t1old[i][a] + t1[i][a]/D1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          residual2 += t2[i][j][a][b] * t2[i][j][a][b];
          t2[i][j][a][b] = t2old[i][j][a][b] + t2[i][j][a][b]/D2[i][j][a][b];
        }
    }

  return sqrt(residual1 + residual2);
}

}} // namespace psi::ugacc
