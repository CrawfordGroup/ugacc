Perturbation::Perturbation(std::string operator, boost::shared_ptr<Wavefunction> ref)
{
  int *map = init_int_array(nmo_); // Translates from Pitzer (including frozen docc) to QT
  reorder_qt((int *) ref->doccpi(), (int *) ref->soccpi(), (int *) ref->frzcpi(), 
             (int *) ref->frzvpi(), map, (int *) ref->nmopi(), ref->nirrep());

  int *mo_offset = init_int_array(ref->nirrep()); // Pitzer offsets
  for(int h=1; h < ref->nirrep(); h++) mo_offset[h] = mo_offset[h-1] + ref->nmopi()[h-1];

  // Prepare MO-basis property integrals in QT ordering
  boost::shared_ptr<Molecule> mol = ref->molecule();
  boost::shared_ptr<IntegralFactory> fact = ref->integral();
  OperatorSymmetry dipsym(1, mol, fact);
  mu_irreps = new int[3];
  mu_irreps[0] = dipsym.component_symmetry(0);
  mu_irreps[1] = dipsym.component_symmetry(1);
  mu_irreps[2] = dipsym.component_symmetry(2);

  SharedMatrix Ca = ref->Ca();
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

  delete [] mu_irreps;
  free(mo_offset);
  free(map);
}

