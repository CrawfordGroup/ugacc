#include <cstdio>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void ccdump(void)
{
  psio_address next;
  int dumpfile = 89;
  int no = moinfo.no; 
  int nv = moinfo.nv;
  int nact = moinfo.nact;
  double **t1 = moinfo.t1; 
  double ****t2 = moinfo.t2;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **fock = moinfo.fock;

  psio_open(dumpfile, PSIO_OPEN_OLD);

  psio_write_entry(dumpfile, "nmo", (char *) &moinfo.nmo, sizeof(int));
  psio_write_entry(dumpfile, "nact", (char *) &moinfo.nact, sizeof(int));
  psio_write_entry(dumpfile, "no", (char *) &moinfo.no, sizeof(int));
  psio_write_entry(dumpfile, "nv", (char *) &moinfo.nv, sizeof(int));
  psio_write_entry(dumpfile, "nfzc", (char *) &moinfo.nfzc, sizeof(int));
  psio_write_entry(dumpfile, "nfzv", (char *) &moinfo.nfzv, sizeof(int));

  psio_write_entry(dumpfile, "escf", (char *) &moinfo.escf, sizeof(double));
  psio_write_entry(dumpfile, "efzc", (char *) &moinfo.efzc, sizeof(double));
  psio_write_entry(dumpfile, "eccsd", (char *) &moinfo.eccsd, sizeof(double));

  next = PSIO_ZERO;
  for(int i=0; i < no; i++) 
    psio_write(dumpfile, "T1 Amplitudes", (char *) t1[i], nv*sizeof(double),
	       next, &next);

  next = PSIO_ZERO;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
	psio_write(dumpfile, "T2 Amplitudes", (char *) t2[i][j][a], 
		   nv*sizeof(double), next, &next);

  next = PSIO_ZERO;
  for(int p=0; p < nact; p++)
    psio_write(dumpfile, "Fock Matrix", (char *) fock[p], nact*sizeof(double),
	       next, &next);

  next = PSIO_ZERO;
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++)
      for(int r=0; r < nact; r++)
	psio_write(dumpfile, "Two-Electron Integrals", (char *) ints[p][q][r],
		   nact*sizeof(double), next, &next);

  next = PSIO_ZERO;
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++)
      for(int r=0; r < nact; r++)
        psio_write(dumpfile, "L Integrals", (char *) L[p][q][r],
                   nact*sizeof(double), next, &next);

  psio_close(dumpfile, 1);
}
}} // namespace psi::ugacc

