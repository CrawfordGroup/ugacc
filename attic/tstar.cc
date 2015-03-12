#include <libciomr/libciomr.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void tstar_build(double **t1, double ****t2)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;

  moinfo.t1s = block_matrix(no,nv);
  moinfo.t2s = init_4d_array(no,no,nv,nv);

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
          moinfo.t1s[i][a] = 2.0 * t1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          moinfo.t2s[i][j][a][b] = 4.0 * t2[i][j][a][b] - 2.0 * t2[i][j][b][a];
    }

  return;
}

}} // namespace psi::ugacc

