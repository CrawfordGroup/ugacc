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

Hamiltonian::Hamiltonian(boost::shared_ptr<PSIO> psio, boost::shared_ptr<Wavefunction> reference, std::vector<boost::shared_ptr<MOSpace> > spaces)
{
  nmo_ = reference->nmo();
  nfzc_ = reference->nfrzc();
  nfzv_ = 0;
  for(int i=0; i < reference->nirrep(); i++) 
    nfzv_ += reference->frzvpi()[i];
  nact_ = nmo_ - nfzc_ - nfzv_;

  int nact = nact_;

  SharedMatrix Fa = reference->Fa();
  SharedMatrix Ca = reference->Ca();
  Fa->transform(Ca);
  Fa->print();
  fock_ = block_matrix(nact, nact);
  for(int h=0; h < Fa->nirrep(); h++)
    for(int p=0; p < nact; p++)
      for(int q=0; q < nact; q++)
        fock_[p][q] = Fa->get(h, p+nfzc_, q+nfzc_);

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
