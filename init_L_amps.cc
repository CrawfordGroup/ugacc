#include <libciomr/libciomr.h>
#include <cstring>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void init_L_amps(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double **t1 = moinfo.t1; 
  double ****t2 = moinfo.t2;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double ****D2 = moinfo.D2;

  l1 = block_matrix(no,nv);
  l1old = block_matrix(no,nv);

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      l1[i][a] = t1[i][a];

  l2 = init_4d_array(no,no,nv,nv);
  l2old = init_4d_array(no,no,nv,nv);

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          l2[i][j][a][b] = t2[i][j][a][b];

  moinfo.l1 = l1; moinfo.l1old = l1old;
  moinfo.l2 = l2; moinfo.l2old = l2old;

  moinfo.Gae = block_matrix(nv, nv);
  moinfo.Gmi = block_matrix(no, no);

  return;
}

}} // namespace devel::rhfcclambda
