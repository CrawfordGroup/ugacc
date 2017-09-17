#
# @BEGIN LICENSE
#
# ugacc by T. Daniel Crawford, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util

def run_ugacc(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    ugacc can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('ugacc')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    if ('wfn' in kwargs):
        if (kwargs['wfn'] == 'ccsd'):
            psi4.core.set_global_option('WFN', 'CCSD')
        elif (kwargs['wfn'] == 'ccsd(t)'):
            psi4.core.set_global_option('WFN', 'CCSD_T')
    scf_wfn = kwargs.get('ref_wfn', None)
    if scf_wfn is None:
        scf_wfn = psi4.driver.scf_helper(name, **kwargs)
    proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), scf_wfn)

    return psi4.core.plugin('ugacc.so', scf_wfn)

def run_ugacc_gradient(name, **kwargs):
    psi4.core.set_global_option('DERTYPE', 'FIRST')
    return run_ugacc(name, **kwargs)

# Integration with driver routines
psi4.driver.procedures['energy']['ugacc'] = run_ugacc
psi4.driver.procedures['gradient']['ugacc'] = run_ugacc_gradient


def exampleFN():
    # Your Python code goes here
    pass
