#include "ccwavefunction.h"
#include "perturbation.h"

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include "globals.h"

INIT_PLUGIN

using namespace boost;

namespace psi { namespace ugacc {

extern "C" 
int read_options(std::string name, Options& options)
{
  if(name == "UGACC" || options.read_globals()) {
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF");
    options.add_str("WFN", "CCSD");
    options.add_str("DERTYPE", "NONE");
    options.add_int("MAXITER", 100);
    options.add_bool("DIIS", true);
    options.add_double("R_CONVERGENCE", 1e-7);
    options.add_bool("OOC", false);
    options.add("OMEGA",new ArrayType());
  }

  return true;
}

extern "C" 
PsiReturnType ugacc(Options& options)
{
  outfile->Printf("\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t*         UGA-CC         *\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\n");

  boost::shared_ptr<Wavefunction> ref = Process::environment.wavefunction();

  // Error trapping – need closed-shell SCF in place
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");
  if(options.get_str("REFERENCE") != "RHF")
    throw PSIEXCEPTION("Only for use with RHF references determinants.");
  for(int h=0; h < ref->nirrep(); h++)
    if(ref->soccpi()[h]) throw PSIEXCEPTION("UGACC is for closed-shell systems only.");

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);
  boost::shared_ptr<Hamiltonian> H(new Hamiltonian(psio, ref, spaces));
  boost::shared_ptr<CCWavefunction> ccwfn(new CCWavefunction(ref, H, options, psio));

  double ecc = ccwfn->compute_energy();

  if(!ccwfn->dertype()) return Success;

  ccwfn->hbar();
  ccwfn->compute_lambda();

  double eone = ccwfn->onepdm();
  double etwo = ccwfn->twopdm();
  double eref = ccwfn->reference_energy();

  outfile->Printf("\tOne-Electron Energy        = %20.14f\n", eone);
  outfile->Printf("\tTwo-Electron Energy        = %20.14f\n", etwo);
  if(ccwfn->wfn() == "CCSD") {
    outfile->Printf("\tCCSD Correlation Energy    = %20.14f (from density)\n", eone+etwo);
    outfile->Printf("\tCCSD Total Energy          = %20.14f (from density)\n", eone+etwo+eref);
    outfile->Printf("\tCCSD Total Energy          = %20.14f (from ccwfn)\n", ecc);
  }
  else if(ccwfn->wfn() == "CCSD_T") {
    outfile->Printf("\tCCSD(T) Correlation Energy = %20.14f (from density)\n", eone+etwo);
    outfile->Printf("\tCCSD(T) Total Energy       = %20.14f (from density)\n", eone+etwo+eref);
    outfile->Printf("\tCCSD(T) Total Energy       = %20.14f (from ccwfn)\n", ecc);
  }

  // Handle external field frequencies
  int nomega;
  std::vector<double> omega;
  if(options["OMEGA"].size() == 0) { // Assume 0.0 E_h for field energy
    nomega = 1;
    omega.push_back(0.0);
  }
  else if(options.["OMEGA"].size() == 1) { // Assume E_h for field energy and read value
    nomega = 1;
    omega.push_back(options["OMEGA"][0].to_double());
  }
  else if(options["OMEGA"].size() >= 2) {
    nomega = count-1;
    omega.resize(nomega);
    std::string units = options["OMEGA"][count-1].to_string();
    for(int i=0; i < nomega; i++) {
      omega[i] = options["OMEGA"][i].to_double();
      if(units == "HZ" || units == "Hz" || units == "hz")
        omega[i] *= pc_h / pc_hartree2J;
      else if(units == "AU" || units == "Au" || units == "au") continue;
      else if(units == "NM" || units == "nm")
        omega[i] = (pc_c*pc_h*1e9)/(omega[i]*pc_hartree2J);
      else if(units == "EV" || units == "ev" || units == "eV")
        omega[i] /= pc_hartree2ev;
      else
        throw PsiException("Error in unit for input field frequencies, should be au, Hz, nm, or eV", __FILE__,__LINE__);
    }
  }

  // Prepare property integrals for perturbed wave functions
  boost::shared_ptr<Perturbation> mu(new Perturbation("Mu", ref));

  // Solve perturbed wave function equations for give perturbation and field frequency
//  std::map<double, boost::shared_ptr<PertCC> > pertcc; // map for all perturbed CC wfns
  boost::shared_ptr<PertCC> pertcc(new PertCC(ccwfn, ));



  // Use perturbed wfns to construct the linear response function
  // Alternatively, build density-based linear response function

  return Success;
}

}} // End namespaces

