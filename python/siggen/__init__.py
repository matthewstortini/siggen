# -*- coding: utf-8 -*-

try:
    __SIGGEN_SETUP__
except NameError:
    __SIGGEN_SETUP__ = False

if not __SIGGEN_SETUP__:
    __all__ = ["Detector"]

    from .Detector import Detector
    from .PPC import PPC
    from .ICPC import ICPC
    from .GEM import GEM
