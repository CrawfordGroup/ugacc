#include <math.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

double t1norm(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double **t1 = moinfo.t1;

  double diag = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      diag += t1[i][a] * t1[i][a];

  return sqrt(diag/(2*no));
}

}} // namespace psi::ugacc
