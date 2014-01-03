#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void tsave(void)
{
  double ****t2tmp = moinfo.t2;
  moinfo.t2 = moinfo.t2old;
  moinfo.t2old = t2tmp;

  double **t1tmp = moinfo.t1;
  moinfo.t1 = moinfo.t1old;
  moinfo.t1old = t1tmp;
}

}} // namespace psi::ugacc

