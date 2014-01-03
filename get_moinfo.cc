#include <cstdio>
#include <psi4-dec.h>
#include <libmints/mints.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void get_moinfo(boost::shared_ptr<Wavefunction> wfn, boost::shared_ptr<Chkpt> chkpt)
{
  moinfo.nmo = chkpt->rd_nmo();
  int nirreps = chkpt->rd_nirreps();
  int iopen = chkpt->rd_iopen();
  char **labels = chkpt->rd_irr_labs();
  int *orbspi = chkpt->rd_orbspi();
  int *clsdpi = chkpt->rd_clsdpi();
  moinfo.enuc = chkpt->rd_enuc();
  moinfo.escf = chkpt->rd_escf();
  moinfo.efzc = chkpt->rd_efzc();

  if(iopen)
    throw PSIEXCEPTION("UGA-CC code for closed-shell RHF references only.");

  int *frdocc = wfn->frzcpi();
  int *fruocc = wfn->frzvpi();

  int *virtpi = new int[nirreps];
  for(int i=0; i < nirreps; i++)
    virtpi[i] = orbspi[i]-clsdpi[i];

  moinfo.no = 0;
  moinfo.nv = 0;
  for(int i=0; i < nirreps; i++) {
    moinfo.no += clsdpi[i] - frdocc[i];
    moinfo.nv += virtpi[i] - fruocc[i];
  }

  moinfo.nfzc = 0; moinfo.nfzv = 0;
  for(int i=0; i < nirreps; i++) {
      moinfo.nfzc += frdocc[i];
      moinfo.nfzv += fruocc[i];
  }

  moinfo.nact = moinfo.nmo - moinfo.nfzc - moinfo.nfzv;
  moinfo.nact = moinfo.nmo - moinfo.nfzc - moinfo.nfzv;

  fprintf(outfile,"\n\tCheckpoint Parameters:\n");
  fprintf(outfile,"\t------------------------\n");
  fprintf(outfile,"\tNumber of irreps        = %d\n", nirreps);
  fprintf(outfile,"\tNumber of MOs           = %d\n", moinfo.nmo);
  fprintf(outfile,"\tNumber of active MOs    = %d\n", moinfo.nact);
  fprintf(outfile,"\tNumber of active occ    = %d\n", moinfo.no);
  fprintf(outfile,"\tNumber of active vir    = %d\n", moinfo.nv);
  fprintf(outfile,"\tNumber of frozen occ    = %d\n", moinfo.nfzc);
  fprintf(outfile,"\tNumber of frozen vir    = %d\n\n", moinfo.nfzv);
  fprintf(outfile, "\tLabel\t# MOs\t# FZDC\t# DOCC\t# VIRT\t# FZVR\n");
  fprintf(outfile, "\t-----\t-----\t------\t------\t------\t------\n");
  for(int i=0; i < nirreps; i++) {
      fprintf(outfile, "\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\n",
              labels[i],orbspi[i],frdocc[i],clsdpi[i],virtpi[i],fruocc[i]);
    }
  fprintf(outfile,"\n\tNuclear Repulsion Energy    = %20.15f\n", moinfo.enuc);
  fprintf(outfile,  "\tFrozen Core Energy          = %20.15f\n", moinfo.efzc);
  fprintf(outfile,  "\tTotal SCF Energy (chkpt)    = %20.15f\n", moinfo.escf);
  fflush(outfile);

  for(int i=0; i < nirreps; i++) delete [] labels[i];
  delete [] labels;
  delete [] orbspi;
  delete [] clsdpi;
  delete [] virtpi;
}

}} // namespace psi::ugacc

