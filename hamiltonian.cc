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

Hamiltonian::Hamiltonian(boost::shared_ptr<PSIO> psio, boost::shared_ptr<Wavefunction> ref, std::vector<boost::shared_ptr<MOSpace> > spaces)
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
  reorder_qt((int *) ref->doccpi(), (int *) ref->soccpi(), (int *) ref->frzcpi(), (int *) ref->frzvpi(), 
             map, (int *) ref->nmopi(), ref->nirrep());

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

  // Prepare MO-basis property integrals in QT ordering
  boost::shared_ptr<Molecule> mol = ref->molecule();
  boost::shared_ptr<IntegralFactory> fact = ref->integral();
  OperatorSymmetry dipsym(1, mol, fact);
  int *mu_irreps = new int[3];
  mu_irreps[0] = dipsym.component_symmetry(0);
  mu_irreps[1] = dipsym.component_symmetry(1);
  mu_irreps[2] = dipsym.component_symmetry(2);

  double **scf = Ca->to_block_matrix();
  MintsHelper mints(Process::environment.options, 0);
  std::vector<SharedMatrix> dipole = mints.so_dipole();
  mu_ = new double** [3];
  for(int i=0; i < 3; i++) {
    double **A = dipole[i]->to_block_matrix();
    double **B = block_matrix(nso_, nmo_);
    double **C = block_matrix(nmo_, nmo_);
    C_DGEMM('n','n',nso_,nmo_,nso_,1,A[0],nso_,scf[0],nmo_,0,B[0],nmo_);
    C_DGEMM('t','n',nmo_,nmo_,nso_,1,scf[0],nmo_,B[0],nmo_,0,C[0],nmo_);
    mu_[i] = block_matrix(nact, nact);
    for(int hl=0; hl < ref->nirrep(); hl++) {
      int hr = hl ^ mu_irreps[i];
      for(int p=ref->frzcpi()[hl]; p < (ref->nmopi()[hl]-ref->frzvpi()[hl]); p++) {
        for(int q=ref->frzcpi()[hr]; q < (ref->nmopi()[hr]-ref->frzvpi()[hr]); q++) {
          int P = map[p+mo_offset[hl]]; int Q = map[q+mo_offset[hr]];
          mu_[i][P-nfzc_][Q-nfzc_] = C[p+mo_offset[hl]][q+mo_offset[hr]];
        }
      }
    }
    free_block(A); free_block(B); free_block(C);
  }
  free_block(scf);
  free(mo_offset);
  free(mu_irreps);
  free(map);


  // Use reorder_qt() to generate a new mapping array w/o frozen core or virtual orbitals
  int *doccpi = init_int_array(ref->nirrep());
  int *nmopi = init_int_array(ref->nirrep());
  int *null = init_int_array(ref->nirrep());
  for(int h=0; h < ref->nirrep(); h++) {
    doccpi[h] = ref->doccpi()[h] - ref->frzcpi()[h];
    nmopi[h] = ref->nmopi()[h] - ref->frzcpi()[h] - ref->frzvpi()[h];
  }
  int *map2 = init_int_array(nact); // Translates from Pitzer (w/o frozen MOs) to QT
  reorder_qt(doccpi, (int *) ref->soccpi(), null, null, map2, nmopi, ref->nirrep());
  free(null); free(nmopi); free(doccpi);

  IntegralTransform ints(ref, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);

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

  for(int i=0; i < 3; i++) free_block(mu_[i]);
}

} // namespace psi
