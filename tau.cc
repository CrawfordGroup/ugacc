#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void tau_build(int iter)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double **t1 = moinfo.t1old;
  double ****t2 = moinfo.t2old;

  if(iter == 1) { // Only allocate on the first iteration
    moinfo.ttau = init_4d_array(no,no,nv,nv);
    moinfo.tau = init_4d_array(no,no,nv,nv);
  }

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          moinfo.tau[i][j][a][b] = t2[i][j][a][b] + t1[i][a] * t1[j][b];
          moinfo.ttau[i][j][a][b] = t2[i][j][a][b] + 0.5 * t1[i][a] * t1[j][b];
        }

  return;
}

}} // namespace psi::ugacc

