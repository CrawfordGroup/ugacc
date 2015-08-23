#include "perturbation.h"

#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include "libparallel/ParallelPrinter.h"

namespace psi { namespace ugacc {

Perturbation::Perturbation(std::string op, boost::shared_ptr<Wavefunction> ref, bool full_virtual_space)
{
  operator_ = op;

  // Basis-set parameters
  int nmo = ref->nmo();
  int nso = ref->nso();
  int nfzc = ref->nfrzc();
  int nfzv = 0;
  std::vector<int> frzvpi(ref->nirrep());
  if(!full_virtual_space)
    for(int i=0; i < ref->nirrep(); i++) { frzvpi[i] = ref->frzvpi()[i]; nfzv += ref->frzvpi()[i]; }
  else 
    for(int i=0; i < ref->nirrep(); i++) frzvpi[i] = 0;
  nact_ = nmo - nfzc - nfzv;

  outfile->Printf("\n\t==> Perturbation = %s <==\n", op.c_str());

  outfile->Printf("\tNMO    = %d; NSO = %d; NFZC = %d; NFZV = %d; NACT = %d\n", nmo, nso, nfzc, nfzv, nact_);

  int *mo_offset = init_int_array(ref->nirrep()); // Pitzer offsets
  for(int h=1; h < ref->nirrep(); h++) mo_offset[h] = mo_offset[h-1] + ref->nmopi()[h-1];

  int *map = init_int_array(nmo); // Translates from Pitzer (including frozen docc) to QT
  reorder_qt((int *) ref->doccpi(), (int *) ref->soccpi(), (int *) ref->frzcpi(), (int *) &frzvpi[0],
             map, (int *) ref->nmopi(), ref->nirrep());

  // Symmetry info
  boost::shared_ptr<Molecule> mol = ref->molecule();
  boost::shared_ptr<IntegralFactory> fact = ref->integral();
  OperatorSymmetry dipsym(1, mol, fact);
  int *prop_irreps;
  if(operator_ == "Mu" || operator_ == "P" || operator_ == "P*") {
    prop_irreps = new int[3];
    prop_irreps[0] = dipsym.component_symmetry(0);
    prop_irreps[1] = dipsym.component_symmetry(1);
    prop_irreps[2] = dipsym.component_symmetry(2);
  }
  else if(operator_ == "L" || operator_ == "L*") {
    prop_irreps = new int[3];
    prop_irreps[0] = dipsym.component_symmetry(1) ^ dipsym.component_symmetry(2);
    prop_irreps[1] = dipsym.component_symmetry(2) ^ dipsym.component_symmetry(0);
    prop_irreps[2] = dipsym.component_symmetry(0) ^ dipsym.component_symmetry(1);
  }
  else if(operator_ == "Q" || operator_ == "RR") {
    prop_irreps = new int[6];
    prop_irreps[0] = dipsym.component_symmetry(0) ^ dipsym.component_symmetry(0);
    prop_irreps[1] = dipsym.component_symmetry(0) ^ dipsym.component_symmetry(1);
    prop_irreps[2] = dipsym.component_symmetry(0) ^ dipsym.component_symmetry(2);
    prop_irreps[3] = dipsym.component_symmetry(1) ^ dipsym.component_symmetry(1);
    prop_irreps[4] = dipsym.component_symmetry(1) ^ dipsym.component_symmetry(2);
    prop_irreps[5] = dipsym.component_symmetry(2) ^ dipsym.component_symmetry(2);
  }

  // Grab the raw SO integrals
  MintsHelper mints(Process::environment.options, 0);
  std::vector<SharedMatrix> prop;
  if(operator_ == "Mu") prop = mints.so_dipole();
  else if(operator_ == "P" || operator_ == "P*") prop = mints.so_nabla();
  else if(operator_ == "L" || operator_ == "L*") prop = mints.so_angular_momentum();
  else if(operator_ == "Q") prop = mints.so_traceless_quadrupole();
  else if(operator_ == "RR") prop = mints.so_quadrupole();

  // Transform and sort to QT ordering
  SharedMatrix Ca = ref->Ca();
  double **scf = Ca->to_block_matrix();

  if(dipole(operator_)) prop_ = new double** [3];
  else if(quadrupole(operator_)) prop_ = new double** [6];

  if(operator_ == "P*" || operator_ == "L*") // take complex conjugate
    for(int i=0; i < 3; i++) prop[i]->scale(-1.0);
  if(operator_ == "L" || operator_ == "L*") // -1/2 in definition of magnetic dipole
    for(int i=0; i < 3; i++) prop[i]->scale(-0.5);

  for(int i=0; i < (dipole(operator_) ? 3 : 6); i++) {
    double **A = prop[i]->to_block_matrix();
    double **B = block_matrix(nso, nmo);
    double **C = block_matrix(nmo, nmo);
    C_DGEMM('n','n',nso,nmo,nso,1,A[0],nso,scf[0],nmo,0,B[0],nmo);
    C_DGEMM('t','n',nmo,nmo,nso,1,scf[0],nmo,B[0],nmo,0,C[0],nmo);
    prop_[i] = block_matrix(nact_, nact_);
    for(int hl=0; hl < ref->nirrep(); hl++) {
      int hr = hl ^ prop_irreps[i];
      for(int p=ref->frzcpi()[hl]; p < (ref->nmopi()[hl]-frzvpi[hl]); p++) {
        for(int q=ref->frzcpi()[hr]; q < (ref->nmopi()[hr]-frzvpi[hr]); q++) {
          int P = map[p+mo_offset[hl]]; int Q = map[q+mo_offset[hr]];
          prop_[i][P-nfzc][Q-nfzc] = C[p+mo_offset[hl]][q+mo_offset[hr]];
        }
      }
    }
    free_block(A); free_block(B); free_block(C);
  }

  free_block(scf);
  delete [] prop_irreps;
  free(mo_offset);
  free(map);
}

Perturbation::~Perturbation()
{
  for(int i=0; i < (dipole(operator_) ? 3 : 6); i++) free_block(prop_[i]);
  delete [] prop_;
}

void Perturbation::print(std::string out)
{
  boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:boost::shared_ptr<OutFile>(new OutFile(out)));

  std::string cart = "XYZ";

  printer->Printf("\n");
  if(dipole(operator_)) {
    for(int i=0; i < 3; i++) {
      printer->Printf("%s(%c)\n", operator_.c_str(), cart[i]);
      mat_print(prop_[i], nact_, nact_, out);
    }
  }
  else if(quadrupole(operator_)) {
    for(int i=0,ij=0; i < 3; i++) {
      for(int j=i; j < 3; j++,ij++) {
        printer->Printf("%s(%c,%c)\n", operator_.c_str(), cart[i], cart[j]);
        mat_print(prop_[ij], nact_, nact_, out);
      }
    }
  }
}

void Perturbation::print(int i, std::string out)
{
  boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:boost::shared_ptr<OutFile>(new OutFile(out)));

  std::string cart = "XYZ";

  printer->Printf("\n");
  if(dipole(operator_)) {
    printer->Printf("%s(%c)\n", operator_.c_str(), cart[i]);
    mat_print(prop_[i], nact_, nact_, out);
  }
  else throw PSIEXCEPTION("Single Cartesian index given for multipolar property?");
}

void Perturbation::print(int i, int j, std::string out)
{
  boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:boost::shared_ptr<OutFile>(new OutFile(out)));

  std::string cart = "XYZ";

  printer->Printf("\n");
  if(dipole(operator_))
   throw PSIEXCEPTION("Two Cartesian indices given for dipolar property?");
  else if(quadrupole(operator_)) {
    int ij = ((i) > (j) ? (i)*((i)+1)/2 + (j) : (j)*((j)+1)/2 + (i));
    printer->Printf("%s(%c,%c)\n", operator_.c_str(), cart[i], cart[j]);
    mat_print(prop_[ij], nact_, nact_, out);
  }
}

void Perturbation::print() { print("outfile"); }
void Perturbation::print(int i) { print(i, "outfile"); }
void Perturbation::print(int i, int j) { print(i, j, "outfile"); }

bool Perturbation::allowed(std::string op)
{
  if(op == "Mu" || op == "P" || op == "P*" || op == "L" || op == "L*" || op == "Q" || op == "RR") return true;
  else return false;
}

bool Perturbation::dipole(std::string op)
{
  if(op == "Mu" || op == "P" || op == "P*" || op == "L" || op == "L*") return true;
  else return false;
}

bool Perturbation::quadrupole(std::string op)
{
  if(op == "Q" || op == "RR") return true;
  else return false;
}

}} // psi::ugacc
