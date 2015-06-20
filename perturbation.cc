#include "perturbation.h"

#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include "libparallel/ParallelPrinter.h"

namespace psi {

Perturbation::Perturbation(std::string op, boost::shared_ptr<Wavefunction> ref)
{
  operator_ = op;

  // Basis-set parameters
  nmo_ = ref->nmo();
  nso_ = ref->nso();
  nfzc_ = ref->nfrzc();
  nfzv_ = 0;
  for(int i=0; i < ref->nirrep(); i++)
    nfzv_ += ref->frzvpi()[i];
  nact_ = nmo_ - nfzc_ - nfzv_;

  int *mo_offset = init_int_array(ref->nirrep()); // Pitzer offsets
  for(int h=1; h < ref->nirrep(); h++) mo_offset[h] = mo_offset[h-1] + ref->nmopi()[h-1];

  int *map = init_int_array(nmo_); // Translates from Pitzer (including frozen docc) to QT
  reorder_qt((int *) ref->doccpi(), (int *) ref->soccpi(), (int *) ref->frzcpi(), (int *) ref->frzvpi(),
             map, (int *) ref->nmopi(), ref->nirrep());

  // Symmetry info
  boost::shared_ptr<Molecule> mol = ref->molecule();
  boost::shared_ptr<IntegralFactory> fact = ref->integral();
  OperatorSymmetry dipsym(1, mol, fact);
  int *prop_irreps = new int[3];
  if(operator_ == "Mu" || operator_ == "P" || operator_ == "P*" || operator_ == "Q") {
    prop_irreps[0] = dipsym.component_symmetry(0);
    prop_irreps[1] = dipsym.component_symmetry(1);
    prop_irreps[2] = dipsym.component_symmetry(2);
  }
  else if(operator_ == "L" || operator_ == "L*") {
    prop_irreps[0] = dipsym.component_symmetry(1) ^ dipsym.component_symmetry(2);
    prop_irreps[1] = dipsym.component_symmetry(2) ^ dipsym.component_symmetry(0);
    prop_irreps[2] = dipsym.component_symmetry(0) ^ dipsym.component_symmetry(1);
  }

  // Grab the raw SO integrals
  MintsHelper mints(Process::environment.options, 0);
  std::vector<SharedMatrix> prop;
  if(operator_ == "Mu") prop = mints.so_dipole();
  else if(operator_ == "P" || operator_ == "P*") prop = mints.so_nabla();
  else if(operator_ == "L" || operator_ == "L*") prop = mints.so_angular_momentum();
  else if(operator_ == "Q") prop = mints.so_traceless_quadrupole();

  // Transform and sort to QT ordering
  SharedMatrix Ca = ref->Ca();
  double **scf = Ca->to_block_matrix();

  if(onebody(operator_)) {
    prop2_ = new double** [3];
    if(operator_ == "P*" || operator_ == "L*") // take complex conjugate
      for(int i=0; i < 3; i++) prop[i]->scale(-1.0);
    if(operator_ == "L" || operator_ == "L*") // -1/2 in definition of magnetic dipole
      for(int i=0; i < 3; i++) prop[i]->scale(-0.5);
    for(int i=0; i < 3; i++) {
      double **A = prop[i]->to_block_matrix();
      double **B = block_matrix(nso_, nmo_);
      double **C = block_matrix(nmo_, nmo_);
      C_DGEMM('n','n',nso_,nmo_,nso_,1,A[0],nso_,scf[0],nmo_,0,B[0],nmo_);
      C_DGEMM('t','n',nmo_,nmo_,nso_,1,scf[0],nmo_,B[0],nmo_,0,C[0],nmo_);
      prop2_[i] = block_matrix(nact_, nact_);
      for(int hl=0; hl < ref->nirrep(); hl++) {
        int hr = hl ^ prop_irreps[i];
        for(int p=ref->frzcpi()[hl]; p < (ref->nmopi()[hl]-ref->frzvpi()[hl]); p++) {
          for(int q=ref->frzcpi()[hr]; q < (ref->nmopi()[hr]-ref->frzvpi()[hr]); q++) {
            int P = map[p+mo_offset[hl]]; int Q = map[q+mo_offset[hr]];
            prop2_[i][P-nfzc_][Q-nfzc_] = C[p+mo_offset[hl]][q+mo_offset[hr]];
          }
        }
      }
      free_block(A); free_block(B); free_block(C);
    }
  }
  else if(twobody(operator_)) {
    prop3_ = new double*** [3];
    for(int i=0; i < 3; i++) prop3_[i] = new double** [3];
    for(int i=0,ij=0; i < 3; i++) {
      for(int j=i; j < 3; j++,ij++) {
        double **A = prop[ij]->to_block_matrix();
        double **B = block_matrix(nso_, nmo_);
        double **C = block_matrix(nmo_, nmo_);
        C_DGEMM('n','n',nso_,nmo_,nso_,1,A[0],nso_,scf[0],nmo_,0,B[0],nmo_);
        C_DGEMM('t','n',nmo_,nmo_,nso_,1,scf[0],nmo_,B[0],nmo_,0,C[0],nmo_);
        prop3_[i][j] = block_matrix(nact_, nact_);  
        for(int hl=0; hl < ref->nirrep(); hl++) {
          int hr = hl ^ prop_irreps[i];
          for(int p=ref->frzcpi()[hl]; p < (ref->nmopi()[hl]-ref->frzvpi()[hl]); p++) {
            for(int q=ref->frzcpi()[hr]; q < (ref->nmopi()[hr]-ref->frzvpi()[hr]); q++) {
              int P = map[p+mo_offset[hl]]; int Q = map[q+mo_offset[hr]];
              prop3_[i][j][P-nfzc_][Q-nfzc_] = C[p+mo_offset[hl]][q+mo_offset[hr]];
            }
          }
        }
        if(i!=j) prop3_[j][i] = prop3_[i][j];
        free_block(A); free_block(B); free_block(C);
      }
    }
  }

  free_block(scf);
  delete [] prop_irreps;
  free(mo_offset);
  free(map);
}

Perturbation::~Perturbation()
{
  if(onebody(operator_)) {
    for(int i=0; i < 3; i++) free_block(prop2_[i]);
    delete [] prop2_;
  } 
  else if(twobody(operator_)) {
    for(int i=0; i < 3; i++) {
      for(int j=i; j < 3; j++) {
        free_block(prop3_[i][j]);
      }
      delete [] prop3_[i];
    }
    delete [] prop3_;
  }
}

void Perturbation::print(std::string out)
{
  boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:boost::shared_ptr<OutFile>(new OutFile(out)));

  std::string cart = "XYZ";

  printer->Printf("\n");
  if(onebody(operator_)) {
    for(int i=0; i < 3; i++) {
      printer->Printf("%s(%c)\n", operator_.c_str(), cart[i]);
      mat_print(prop2_[i], nact_, nact_, out);
    }
  }
  else if(twobody(operator_)) {
    for(int i=0; i < 3; i++)
      for(int j=0; j <= i; j++) {
        printer->Printf("%s(%c,%c)\n", operator_.c_str(), cart[i], cart[j]);
        mat_print(prop3_[i][j], nact_, nact_, out);
      }
  }
}

void Perturbation::print(int i, std::string out)
{
  boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:boost::shared_ptr<OutFile>(new OutFile(out)));

  std::string cart = "XYZ";

  printer->Printf("\n");
  if(onebody(operator_)) {
    printer->Printf("%s(%c)\n", operator_.c_str(), cart[i]);
    mat_print(prop2_[i], nact_, nact_, out);
  }
  else throw PSIEXCEPTION("Single Cartesian index given for multipolar property?");
}

void Perturbation::print(int i, int j, std::string out)
{
  boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:boost::shared_ptr<OutFile>(new OutFile(out)));

  std::string cart = "XYZ";

  printer->Printf("\n");
  if(onebody(operator_))
   throw PSIEXCEPTION("Two Cartesian indices given for dipolar property?");
  else if(twobody(operator_)) {
    printer->Printf("%s(%c)\n", operator_.c_str(), cart[i]);
    mat_print(prop3_[i][j], nact_, nact_, out);
  }
}

void Perturbation::print() { print("outfile"); }
void Perturbation::print(int i) { print(i, "outfile"); }
void Perturbation::print(int i, int j) { print(i, j, "outfile"); }

bool Perturbation::allowed(std::string op)
{
  if(op == "Mu" || op == "P" || op == "P*" || op == "L" || op == "L*" || op == "Q") return true;
  else return false;
}

bool Perturbation::onebody(std::string op)
{
  if(op == "Mu" || op == "P" || op == "P*" || op == "L" || op == "L*") return true;
  else return false;
}

bool Perturbation::twobody(std::string op)
{
  if(op == "Q") return true;
  else return false;
}

} // psi
