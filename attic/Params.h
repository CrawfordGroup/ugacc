#ifndef _psi_psi_ugacc_params_h
#define _psi_psi_ugacc_params_h

#include <string>

namespace psi { namespace ugacc {

struct Params {
    std::string ref;
    std::string wfn;
    int dertype;
    int maxiter;
    double convergence;
    int do_diis;        // DIIS boolean
    int ooc;            // out-of-core boolean
};

}} // namespace devel::ugacc

#endif // _psi_psi_ugacc_params_h
