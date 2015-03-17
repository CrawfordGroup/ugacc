#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

double energy(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;

  double one_energy=0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      one_energy = 2 * fock[i][a+no] * t1[i][a];

  double two_energy=0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          two_energy += (2*ints[i][j][a+no][b+no]-ints[i][j][b+no][a+no]) *
                        (t2[i][j][a][b] + t1[i][a] * t1[j][b]);

  return one_energy + two_energy;
}

}} // namespace psi::ugacc
