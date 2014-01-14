#ifndef _psi_psi_ugacc_moinfo_h
#define _psi_psi_ugacc_moinfo_h

namespace psi { namespace ugacc {

struct MOInfo {
  int nmo;          /* # molecular orbitals */
  int nso;          /* # symmetry orbitals */
  int nact;         /* # active orbitals */
  int no;           /* # occupied orbitals */
  int nv;           /* # unoccupied orbitals */
  int nfzc;         /* # frozen core orbitals */
  int nfzv;         /* # frozen virtual orbitals */
  double enuc;      /* nuclear repulsion energy */
  double escf;      /* SCF energy */
  double efzc;      /* frozen core energy */
  double eccsd;     /* CCSD energy */
  double e_t;       /* (T) energy */
  double **fock;    /* f(p,q) */
  double ****ints;  /* <pq|rs> */
  double ****L;     /* 2<pq|rs> - <pq|sr> */

  // T-amplitude quantities
  double **D1;      /* one-electron denominators */
  double ****D2;    /* two-electron denominators */
  double **t1;      /* current t1 amplitudes */
  double **t1old;   /* previous t1 amplitudes */
  double ****t2;    /* current t2 amplitudes */
  double ****t2old; /* previous t2 amplitudes */
  double ****ttau;  /* tau-tilde effective doubles */
  double ****tau;   /* tau effective doubles */

  // L-amplitude quantities
  double **l1;      /* current l1 amplitudes */
  double **l1old;   /* previous l1 amplitudes */
  double ****l2;    /* current l2 amplitudes */
  double ****l2old; /* previous l2 amplitudes */

  // CCSD intermediates for amplitude equations
  double **Fae;     /* CC intermediate */
  double **Fmi;     /* CC intermediate */
  double **Fme;     /* CC intermediate */
  double ****Wmnij; /* CC intermediate */
  double ****Wmbej; /* CC intermediate */
  double ****Wmbje; /* CC intermediate */

  // CCSD HBAR
  double **Hoo;
  double **Hvv;
  double **Hov;
  double ****Hoooo;
  double ****Hvvvv;
  double ****Hovov;
  double ****Hovvo;
  double ****Hvovv;
  double ****Hooov;
  double ****Hovoo;
  double ****Hvvvo;

  // Three-body intermediates
  double **Gvv;
  double **Goo;

  // Additional contributions to Lambda equations from (T) correction
  double **s1;
  double ****s2;

  // One-electron density components
  double **Doo;
  double **Dvv;
  double **Dov;
  double **Dvo;

  // Triples
  double ******t3;
  double ******l3;
};

}} // namespace psi::ugacc

#endif // _psi_psi_ugacc_moinfo_h
