#include "ccwavefunction.h"
#include "hamiltonian.h"
#include <boost/shared_ptr.hpp>

namespace psi {

CCWavefunction::CCWavefunction(boost::shared_ptr<Wavefunction> reference, 
                               boost::shared_ptr<Hamiltonian> H,
                               Options &options, boost::shared_ptr<PSIO>
psio) : Wavefunction(options, psio)
{
  outfile->Printf("\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t*         UGA-CC         *\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\n");

  if(options.get_str("REFERENCE") != "RHF")
    throw PSIEXCEPTION("Only for use with RHF references determinants.");

  wfn_ = options.get_str("WFN");
  convergence_ = options.get_double("R_CONVERGENCE");
  maxiter_ = options.get_int("MAXITER");
  do_diis_ = options.get_bool("DIIS");
  ooc_ = options.get_bool("OOC");

  outfile->Printf("\tWave function  = %s\n", wfn().c_str());
  outfile->Printf("\tMaxiter        = %d\n", maxiter());
  outfile->Printf("\tConvergence    = %3.1e\n", convergence());
  outfile->Printf("\tDIIS           = %s\n", do_diis() ? "Yes" : "No");
  outfile->Printf("\tOut-of-core    = %s\n", ooc() ? "Yes" : "No");
}

CCWavefunction::~CCWavefunction() { }

double CCWavefunction::compute_energy() { return 0.0; }

} // psi
