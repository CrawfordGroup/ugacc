import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from p4xcpt import *


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
            psi4.set_global_option('WFN', 'CCSD')
        elif (kwargs['wfn'] == 'ccsd(t)'):
            psi4.set_global_option('WFN', 'CCSD_T')
    scf_helper(name, **kwargs)
    psi4.transqt2()
    returnvalue = psi4.plugin('ugacc.so')
#    psi4.set_variable('CURRENT ENERGY', returnvalue)

def run_ugacc_gradient(name, **kwargs):
    psi4.set_global_option('DERTYPE', 'FIRST')
    run_ugacc(name, **kwargs)

# Integration with driver routines
procedures['energy']['ugacc'] = run_ugacc
procedures['gradient']['ugacc'] = run_ugacc_gradient


def exampleFN():
    # Your Python code goes here
    pass
