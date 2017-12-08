# -*- coding: utf-8 -*-

try:
    __SIGGEN_SETUP__
except NameError:
    __SIGGEN_SETUP__ = False

if not __SIGGEN_SETUP__:
    __all__ = ["Detector"]

    # from ._siggen import PyICPC,PyGEM, PyPPC
    # from .PyDetector import PyDetector
    from .Detector import Detector
    from .ICPC import ICPC
    from .PPC import PPC
    from .GEM import GEM
    # from .make_fields import read_fields
