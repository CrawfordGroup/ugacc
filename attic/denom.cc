#include <libciomr/libciomr.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void denom(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;

  moinfo.D1 = block_matrix(no,nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      moinfo.D1[i][a] = fock[i][i] - fock[a+no][a+no];

  moinfo.D2 = init_4d_array(no,no,nv,nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          moinfo.D2[i][j][a][b] = fock[i][i] + fock[j][j] - fock[a+no][a+no] - fock[b+no][b+no];

  return;
}

}} // namespace psi::ugacc
