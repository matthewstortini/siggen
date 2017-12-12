# distutils: language = c++

from _siggen cimport *

include "geometries/TemplateWrapper_PPC.pxi"
include "geometries/TemplateWrapper_GEM.pxi"
include "geometries/TemplateWrapper_ICPC.pxi"

include "geometries/PPCGeometry.pxi"
include "geometries/ICPCGeometry.pxi"
include "geometries/GEMGeometry.pxi"

include "PySiggen.pxi"

cdef class PySetup:
  cdef Setup* cobj
  def __cinit__(self, conf_file):
    self.cobj = new Setup(conf_file)

  def __dealloc__(self):
    del self.cobj

# include "fields.pxi"
