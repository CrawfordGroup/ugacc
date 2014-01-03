#include <cstdio>
#include <psi4-dec.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void integrals(boost::shared_ptr<PSIO> psio)
{
  int *ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(int i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;

  int nmo = moinfo.nmo;
  int no = moinfo.no;
  int nact = moinfo.nact;
  int nfzc = moinfo.nfzc;
  int nfzv = moinfo.nfzv;

  /* One-electron integrals/frozen-core operator contribution */  
  int noei_all = nmo*(nmo+1)/2;
  int noei = nact*(nact+1)/2;
  double *scratch = init_array(noei_all); 
  double *oei = init_array(noei);
  iwl_rdone(PSIF_OEI, PSIF_MO_FZC, scratch, noei_all, 0, 0, outfile);
  filter(scratch, oei, ioff, nmo, nfzc, nfzv);

  int ntei = noei*(noei+1)/2;
  double *tei = init_array(ntei);
  struct iwlbuf Buf;
  iwl_buf_init(&Buf, PSIF_MO_TEI, 1e-14, 1, 1);
  iwl_buf_rd_all(&Buf, tei, ioff, ioff, 0, ioff, 0, outfile);
  iwl_buf_close(&Buf, 1); /* keep the integral file */

  // This does not appear to work -- don't know why, but would like to debug
  //  IWL::read_two(psio.get(), PSIF_MO_TEI, tei, ioff, nmo, nfzc, nfzv, 1, outfile);

  double ****ints = init_4d_array(nact, nact, nact, nact);
  for(int p=0; p < nact; p++)
    for(int r=0; r < nact; r++) {
      int pr = INDEX(p,r);
      for(int q=0; q < nact; q++)
        for(int s=0; s < nact; s++) {
          int qs = INDEX(q,s);
          int prqs = INDEX(pr,qs);
          ints[p][q][r][s] = tei[prqs];
        }
    }

  // L(pqrs) = 2<pq|rs> - <pq|sr>
  double ****L = init_4d_array(nact, nact, nact, nact);
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++)
      for(int r=0; r < nact; r++)
        for(int s=0; s < nact; s++)
          L[p][q][r][s] = 2*ints[p][q][r][s] - ints[p][q][s][r];

  // Build the Fock matrix
  double **fock = block_matrix(nact, nact);
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++) {
      fock[p][q] = oei[INDEX(p,q)];
      for(int m=0; m < no; m++)
        fock[p][q] += L[p][m][q][m];
    }

  // Recompute the SCF energy as a sanity check
  double escf = 0.0;
  for(int i=0; i < no; i++) {
    escf += 2 * oei[INDEX(i,i)];
    for(int j=0; j < no; j++)
      escf += 2 * ints[i][j][i][j] - ints[i][j][j][i];
  }
  fprintf(outfile, "\tSCF energy (recomputed)     = %20.15f\n", escf+moinfo.enuc+moinfo.efzc);

  free(ioff);
  free(oei); 
  free(tei);
  free(scratch);

  moinfo.fock = fock;
  moinfo.ints = ints;
  moinfo.L = L;
}

}} // namespace psi::ugacc

