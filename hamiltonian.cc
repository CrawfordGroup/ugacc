#include "hamiltonian.h"
#include "globals.h"
#include <libiwl/iwl.h>
#include <libmints/wavefunction.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

namespace psi {

#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

Hamiltonian::Hamiltonian(boost::shared_ptr<Wavefunction> reference)
{

  nmo_ = reference->nmo();
  nfzc_ = reference->nfrzc();
  no_ = 0;
  nfzv_ = 0;
  for(int i=0; i < reference->nirrep(); i++) {
    no_ += reference->doccpi()[i] - reference->frzcpi()[i];
    nfzv_ += reference->frzvpi()[i];
  }
  nact_ = nmo_ - nfzc_ - nfzv_;

/*
  outfile->Printf("\n\tHamiltonian Parameters:\n");
  outfile->Printf("\t-------------------------\n");
  outfile->Printf("\tNumber of MOs           = %d\n", nmo_);
  outfile->Printf("\tNumber of active MOs    = %d\n", nact_);
  outfile->Printf("\tNumber of active occ    = %d\n", no_);
  outfile->Printf("\tNumber of frozen occ    = %d\n", nfzc_);
*/

  int *ioff = new int[IOFF_MAX];
  ioff[0] = 0;
  for(int i=0; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;

  int noei_all = nmo_*(nmo_+1)/2;
  int noei = nact_*(nact_+1)/2;
  double *scratch = new double[noei_all];
  std::memset(static_cast<void*>(scratch), '\0', noei_all*sizeof(double));
  double *oei = new double[noei];
  std::memset(static_cast<void*>(oei), '\0', noei*sizeof(double));
  iwl_rdone(PSIF_OEI, PSIF_MO_FZC, scratch, noei_all, 0, 0, "outfile");
  filter(scratch, oei, ioff, nmo_, nfzc_, nfzv_);

  int ntei = noei*(noei+1)/2;
  double *tei = new double[ntei];
  std::memset(static_cast<void*>(tei), '\0', ntei*sizeof(double));
  struct iwlbuf Buf;
  iwl_buf_init(&Buf, PSIF_MO_TEI, 1e-16, 1, 1);
  iwl_buf_rd_all(&Buf, tei, ioff, ioff, 0, ioff, 0, "outfile");
  iwl_buf_close(&Buf, 1);

  int no = no_;
  int nact = nact_;

  ints_ = init_4d_array(nact, nact, nact, nact);
  for(int p=0; p < nact; p++)
    for(int r=0; r < nact; r++) {
      int pr = INDEX(p,r);
      for(int q=0; q < nact; q++)
        for(int s=0; s < nact; s++) {
          int qs = INDEX(q,s);
          int prqs = INDEX(pr,qs);
          ints_[p][q][r][s] = tei[prqs];
        }
    }

  // L(pqrs) = 2<pq|rs> - <pq|sr>  
  L_ = init_4d_array(nact, nact, nact, nact);
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++)
      for(int r=0; r < nact; r++)
        for(int s=0; s < nact; s++)
          L_[p][q][r][s] = 2*ints_[p][q][r][s] - ints_[p][q][s][r];

  // Build the Fock matrix
  fock_ = block_matrix(nact, nact);
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++) {
      fock_[p][q] = oei[INDEX(p,q)];
      for(int m=0; m < no; m++)
        fock_[p][q] += L_[p][m][q][m];
    }

  delete[] ioff;
  delete[] oei;
  delete[] tei;
  delete[] scratch;
}

Hamiltonian::~Hamiltonian()
{
  free_4d_array(ints_, nact_, nact_, nact_);
  free_4d_array(L_, nact_, nact_, nact_);
  free_block(fock_); 
}

/*
Hamiltonian::Hamiltonian(const boost::shared_ptr<Hamiltonian> &H)
{
  nmo_ = H->nmo_;
  nact_ = H->nact_;
  no_ = H->nact_;
  nfzc_ = H->nfzc_;
  nfzv_ = H->nfzv_;

  int nact = nact_;

  fock_ = block_matrix(nact, nact);
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++)
      fock_[p][q] = H->fock_[p][q];

  ints_ = init_4d_array(nact, nact, nact, nact);
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++)
      for(int r=0; r < nact; r++)
        for(int s=0; s < nact; s++) {
          ints_[p][q][r][s] = H->ints_[p][q][r][s];
          L_[p][q][r][s] = H->L_[p][q][r][s];
        }
}
*/

} // namespace psi
