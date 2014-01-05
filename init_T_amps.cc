#include <libciomr/libciomr.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void init_T_amps(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double **D1 = moinfo.D1; 
  double ****D2 = moinfo.D2;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;

  double **t1 = block_matrix(no,nv);
  double **t1old = block_matrix(no,nv);

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      t1[i][a] = fock[i][a+no]/D1[i][a]; 

  double ****t2 = init_4d_array(no,no,nv,nv);
  double ****t2old = init_4d_array(no,no,nv,nv);

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          t2[i][j][a][b] = ints[i][j][a+no][b+no]/D2[i][j][a][b];

  moinfo.t1 = t1; moinfo.t1old = t1old;
  moinfo.t2 = t2; moinfo.t2old = t2old;
}

}} // namespace psi::ugacc

