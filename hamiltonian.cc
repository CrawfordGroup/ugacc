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
#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

namespace psi {

Hamiltonian::Hamiltonian(boost::shared_ptr<PSIO> psio, boost::shared_ptr<Wavefunction> ref, std::vector<boost::shared_ptr<MOSpace> > spaces)
{
  nmo_ = ref->nmo();
  nfzc_ = ref->nfrzc();
  nfzv_ = 0;
  for(int i=0; i < ref->nirrep(); i++) 
    nfzv_ += ref->frzvpi()[i];
  nact_ = nmo_ - nfzc_ - nfzv_;

  int nact = nact_;

  // Prepare Fock matrix in MO basis
  SharedMatrix Fa = ref->Fa();
  SharedMatrix Ca = ref->Ca();
  Fa->transform(Ca);

  // Translate from Pitzer to QT
  int *map = init_int_array(nmo_);
  reorder_qt((int *) ref->doccpi(), (int *) ref->soccpi(), (int *) ref->frzcpi(), (int *) ref->frzvpi(), 
             map, (int *) ref->nmopi(), ref->nirrep());

  fock_ = block_matrix(nact, nact);
  int mo_offset=0;
  for(int h=0; h < ref->nirrep(); h++) {
    int nmo = ref->nmopi()[h]; int nfv = ref->frzvpi()[h]; int nfc = ref->frzcpi()[h];
    for(int p=nfc; p < nmo-nfv; p++) {
      for(int q=nfc; q < nmo-nfv; q++) {
      int P = map[p+mo_offset]; int Q = map[q+mo_offset];
      fock_[P-nfzc_][Q-nfzc_] = Fa->get(h,p,q);
      }
    }
    mo_offset += nmo;
  }


  IntegralTransform ints(ref, spaces, IntegralTransform::Restricted, IntegralTransform::IWLAndDPD);
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);

  int *ioff = new int[IOFF_MAX];
  ioff[0] = 0;
  for(int i=0; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;

  int noei = nact_*(nact_+1)/2;
  int ntei = noei*(noei+1)/2;
  double *tei = new double[ntei];
  std::memset(static_cast<void*>(tei), '\0', ntei*sizeof(double));
  struct iwlbuf Buf;
  iwl_buf_init(&Buf, PSIF_MO_TEI, 1e-16, 1, 1);
  iwl_buf_rd_all(&Buf, tei, ioff, ioff, 0, ioff, 0, "outfile");
  iwl_buf_close(&Buf, 1);

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
  delete[] tei;
  delete[] ioff;

/*
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  dpd_set_default(ints.get_dpd_id());
  dpdbuf4 K;
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
  for(int h=0; h < ref->nirrep(); h++) {
    global_dpd_->buf4_mat_irrep_init(&K, h);
    global_dpd_->buf4_mat_irrep_rd(&K, h);
    for(int pq=0; pq < K.params->rowtot[h]; pq++) {
      int p = map[ K.params->roworb[h][pq][0] ];
      int q = map[ K.params->roworb[h][pq][1] ];
      for(int rs=0; rs < K.params->coltot[h]; rs++) {
        int r = map[ K.params->colorb[h][rs][0] ];
        int s = map[ K.params->colorb[h][rs][1] ];
        ints_[p][r][q][s] = K.matrix[h][pq][rs];
        if(fabs(ints_[p][q][r][s]) > 1e-6) outfile->Printf("%d %d %d %d %20.14f\n", p, q, r, s, ints_[p][q][r][s]);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_LIBTRANS_DPD, 1);
*/

  // L(pqrs) = 2<pq|rs> - <pq|sr>  
  L_ = init_4d_array(nact, nact, nact, nact);
  for(int p=0; p < nact; p++)
    for(int q=0; q < nact; q++)
      for(int r=0; r < nact; r++)
        for(int s=0; s < nact; s++)
          L_[p][q][r][s] = 2*ints_[p][q][r][s] - ints_[p][q][s][r];
}

Hamiltonian::~Hamiltonian()
{
  free_4d_array(ints_, nact_, nact_, nact_);
  free_4d_array(L_, nact_, nact_, nact_);
  free_block(fock_); 
}

} // namespace psi
