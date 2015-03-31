#include "hamiltonian.h"
#include "globals.h"
#include <libiwl/iwl.h>
#include <libmints/wavefunction.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libtrans/integraltransform.h>
#include <libdpd/dpd.h>

#define ID(x) ints.DPD_ID(x)

namespace psi {

#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

Hamiltonian::Hamiltonian(boost::shared_ptr<PSIO> psio, boost::shared_ptr<Wavefunction> reference, std::vector<boost::shared_ptr<MOSpace> > spaces)
{
  nmo_ = reference->nmo();
  nfzc_ = reference->nfrzc();
  nfzv_ = 0;
  for(int i=0; i < reference->nirrep(); i++) 
    nfzv_ += reference->frzvpi()[i];
  nact_ = nmo_ - nfzc_ - nfzv_;

/*
  outfile->Printf("\n\tHamiltonian Parameters:\n");
  outfile->Printf("\t-------------------------\n");
  outfile->Printf("\tNumber of MOs           = %d\n", nmo_);
  outfile->Printf("\tNumber of active MOs    = %d\n", nact_);
  outfile->Printf("\tNumber of active occ    = %d\n", no_);
  outfile->Printf("\tNumber of frozen occ    = %d\n", nfzc_);
*/

  int nact = nact_;
  SharedMatrix Fa = reference->Fa();
  SharedMatrix Ca = reference->Ca();
  Fa->transform(Ca);
  double **fock_ = block_matrix(nact, nact);
  for(int h=0; h < Fa->nirrep(); i++)
    for(int p=0; p < Fa->rowspi(h); p++)
      for(int q=0; q < Fa->colspi(h); q++)
        fock[p][q] = 
/*
  fock_ = block_matrix(nact, nact);
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++)
      fock_[p][q] = fock_p[p][q];
*/
  mat_print(fock_, nact, nact, "outfile");

  IntegralTransform ints(reference, spaces, IntegralTransform::Restricted);
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
  dpd_set_default(ints.get_dpd_id());

  dpdbuf4 K;
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
  global_dpd_->buf4_mat_irrep_init(&K, 0); // symmetry = c1
  global_dpd_->buf4_mat_irrep_rd(&K, 0);

  ints_ = init_4d_array(nact, nact, nact, nact);
  for(int pq=0; pq < K.params->rowtot[0]; pq++) {
    int p = K.params->roworb[0][pq][0];
    int q = K.params->roworb[0][pq][1];
    for(int rs=0; rs < K.params->coltot[0]; rs++) {
      int r = K.params->colorb[0][rs][0];
      int s = K.params->colorb[0][rs][1];
      ints_[p][r][q][s] = K.matrix[0][pq][rs];
    }
  }
  global_dpd_->buf4_mat_irrep_close(&K, 0);
  global_dpd_->buf4_close(&K);

  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);


  return;

/*

  double *scratch = new double[noei_all];
  std::memset(static_cast<void*>(scratch), '\0', noei_all*sizeof(double));
  double *oei = new double[noei];
  std::memset(static_cast<void*>(oei), '\0', noei*sizeof(double));
*/
//  iwl_rdone(PSIF_OEI, PSIF_MO_FZC, scratch, noei_all, 0, 0, "outfile");
//  filter(scratch, oei, ioff, nmo_, nfzc_, nfzv_);

//  struct iwlbuf Buf;
//  iwl_buf_init(&Buf, PSIF_MO_TEI, 1e-16, 1, 1);
//  iwl_buf_rd_all(&Buf, tei, ioff, ioff, 0, ioff, 0, "outfile");
//  iwl_buf_close(&Buf, 1);

/*
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
*/
}

Hamiltonian::~Hamiltonian()
{
  free_4d_array(ints_, nact_, nact_, nact_);
/*
  free_4d_array(L_, nact_, nact_, nact_);
  free_block(fock_); 
*/
}

} // namespace psi
