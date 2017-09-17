#include "hamiltonian.h"
#include "array.h"
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/psifiles.h>
#include <psi4/libciomr/libciomr.h>
#include <psi4/libqt/qt.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libdpd/dpd.h>

#define ID(x) ints.DPD_ID(x)

namespace psi { namespace ugacc {

Hamiltonian::Hamiltonian(shared_ptr<PSIO> psio, shared_ptr<Wavefunction> ref, std::vector<shared_ptr<MOSpace> > spaces)
{
  nmo_ = ref->nmo();
  nso_ = ref->nso();
  nfzc_ = ref->nfrzc();
  nfzv_ = 0;
  for(int i=0; i < ref->nirrep(); i++)
    nfzv_ += ref->frzvpi()[i];
  nact_ = nmo_ - nfzc_ - nfzv_;

  int nact = nact_;

  int *mo_offset = init_int_array(ref->nirrep()); // Pitzer offsets
  for(int h=1; h < ref->nirrep(); h++) mo_offset[h] = mo_offset[h-1] + ref->nmopi()[h-1];

  int *map = init_int_array(nmo_); // Translates from Pitzer (including frozen docc) to QT
  Dimension doccpi = ref->doccpi();
  Dimension soccpi = ref->soccpi();
  Dimension frzcpi = ref->frzcpi();
  Dimension frzvpi = ref->frzvpi();
  Dimension nmopi = ref->nmopi();
  reorder_qt((int *) doccpi, (int *) soccpi, (int *) frzcpi, (int *) frzvpi,
             map, (int *) nmopi, ref->nirrep());

  // Prepare Fock matrix in MO basis in QT ordering
  SharedMatrix Fa = ref->Fa();
  SharedMatrix Ca = ref->Ca();
  Fa->transform(Ca);
  fock_ = block_matrix(nact, nact);
  for(int h=0; h < ref->nirrep(); h++) {
    int nmo = ref->nmopi()[h]; int nfv = ref->frzvpi()[h]; int nfc = ref->frzcpi()[h];
    for(int p=nfc; p < nmo-nfv; p++) {
      for(int q=nfc; q < nmo-nfv; q++) {
      int P = map[p+mo_offset[h]]; int Q = map[q+mo_offset[h]];
      fock_[P-nfzc_][Q-nfzc_] = Fa->get(h,p,q);
      }
    }
  }

  free(mo_offset);
  free(map);

  // Use reorder_qt() to generate a new mapping array w/o frozen core or virtual orbitals
  int *null = init_int_array(ref->nirrep());
  for(int h=0; h < ref->nirrep(); h++) {
    doccpi[h] = ref->doccpi()[h] - ref->frzcpi()[h];
    nmopi[h] = ref->nmopi()[h] - ref->frzcpi()[h] - ref->frzvpi()[h];
  }
  int *map2 = init_int_array(nact); // Translates from Pitzer (w/o frozen MOs) to QT
  reorder_qt((int*) doccpi, (int*) soccpi, null, null, map2, (int*) nmopi, ref->nirrep());

  IntegralTransform ints(ref, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
  efzc_ = ints.get_frozen_core_energy();

  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  dpd_set_default(ints.get_dpd_id());
  dpdbuf4 K;
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
  ints_ = init_4d_array(nact, nact, nact, nact);
  for(int h=0; h < ref->nirrep(); h++) {
    global_dpd_->buf4_mat_irrep_init(&K, h);
    global_dpd_->buf4_mat_irrep_rd(&K, h);
    for(int pq=0; pq < K.params->rowtot[h]; pq++) {
      int p = map2[ K.params->roworb[h][pq][0] ];
      int q = map2[ K.params->roworb[h][pq][1] ];
      for(int rs=0; rs < K.params->coltot[h]; rs++) {
        int r = map2[ K.params->colorb[h][rs][0] ];
        int s = map2[ K.params->colorb[h][rs][1] ];
        ints_[p][r][q][s] = K.matrix[h][pq][rs];
      }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_LIBTRANS_DPD, 1);

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

}} // namespace psi::ugacc
