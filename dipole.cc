#include <string>
#include <cstdio>
#include <cstdlib>
#include <psi4-dec.h>
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

using namespace std;

namespace psi { namespace ugacc {

void dipole(boost::shared_ptr<Chkpt> chkpt)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  int nso = moinfo.nso;
  int nmo = moinfo.nmo;

  // Prepare full-matrix density
  double **D = block_matrix(nmo, nmo);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      D[i][j] = moinfo.Doo[i][j];
//      D[i][j] = moinfo.Doo[i][j] + 2.0 * (i==j); // includes SCF contribution

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      D[a+no][b+no] = moinfo.Dvv[a][b];

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      D[i][a+no] = moinfo.Dov[i][a];
      D[a+no][i] = moinfo.Dvo[a][i];
    }

//  mat_print(D, nmo, nmo, outfile);

  // backtransform density to AO basis
  double **scf = chkpt->rd_scf();
  double **TMP = block_matrix(nso, nmo);
  double **DAO = block_matrix(nso, nmo);
  C_DGEMM('n', 'n', nso, nmo, nmo, 1, scf[0], nmo, D[0], nmo, 0, TMP[0], nmo);
  C_DGEMM('n', 't', nso, nso, nmo, 1, TMP[0], nmo, scf[0], nmo, 0, DAO[0], nmo);
  free_block(TMP);
  MintsHelper mints(Process::environment.options, 0);
  vector<SharedMatrix> dipole = mints.so_dipole();

  fprintf(outfile, "\n");
  for(int i=0; i < 3; i++) {
    double **TMP1 = dipole[i]->to_block_matrix();
    double mu = 0.0;
    for(int p=0; p < nso; p++)
      for(int q=0; q < nso; q++)
        mu += DAO[p][q] * TMP1[p][q];
    fprintf(outfile, "\tUnrelaxed Mu[%d] = %20.14f (AO density)\n", i, mu);
  }

  fprintf(outfile, "\n");
  // Try the MO basis
  double **TMP2 = block_matrix(nso,nso);
  for(int i=0; i < 3; i++) {
    double **TMP3 = dipole[i]->to_block_matrix();
    C_DGEMM('n','n',nso,nmo,nso,1,TMP3[0],nso,scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,scf[0],nmo,TMP2[0],nso,0,TMP3[0],nmo);
    double mu = 0.0;
    for(int p=0; p < nmo; p++)
      for(int q=0; q < nmo; q++)
        mu += D[p][q] * TMP3[p][q];
    fprintf(outfile, "\tUnrelaxed Mu[%d] = %20.14f (MO density)\n", i, mu);
  }

  free_block(TMP2);
  free_block(D);

  return;
}

}} // namespace psi::ugacc
