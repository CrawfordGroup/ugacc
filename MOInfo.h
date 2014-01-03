#ifndef _psi_psi_ugacc_moinfo_h
#define _psi_psi_ugacc_moinfo_h

namespace psi { namespace ugacc {

struct MOInfo {
  int nmo;          /* # molecular orbitals */
  int nact;         /* # active orbitals */
  int no;           /* # occupied orbitals */
  int nv;           /* # unoccupied orbitals */
  int nfzc;         /* # frozen core orbitals */
  int nfzv;         /* # frozen virtual orbitals */
  double enuc;      /* nuclear repulsion energy */
  double escf;      /* SCF energy */
  double efzc;      /* frozen core energy */
  double eccsd;     /* CCSD energy */
  double **fock;    /* f(p,q) */
  double ****ints;  /* <pq|rs> */
  double ****L;     /* 2<pq|rs> - <pq|sr> */
  double **D1;      /* one-electron denominators */
  double ****D2;    /* two-electron denominators */
  double **t1;      /* current t1 amplitudes */
  double **t1old;   /* previous t1 amplitudes */
  double ****t2;    /* current t2 amplitudes */
  double ****t2old; /* previous t2 amplitudes */

  double **Fae;     /* CC intermediate */
  double **Fmi;     /* CC intermediate */
  double **Fme;     /* CC intermediate */
  double ****Wmnij; /* CC intermediate */
  double ****Wmbej; /* CC intermediate */
  double ****Wmbje; /* CC intermediate */
  double ****ttau;  /* tau-tilde array */
  double ****tau;   /* tau array */
};

}} // namespace psi::ugacc

#endif // _psi_psi_ugacc_moinfo_h
