from cython.operator cimport dereference
import cython

cdef class GEMGeometry:
  cdef GEM* cobj

  def __cinit__(self, PySetup setup):
    self.cobj = new GEM(dereference(setup.cobj))

  def __dealloc__(self):
    del self.cobj

  property xtal_radius:
    def __get__(self):
      return self.cobj.get_xtal_radius()
    # def __set__(self, float rad):
    #     self.cobj.some_var = var
  property xtal_length:
    def __get__(self):
      return self.cobj.get_xtal_length()
    # def __set__(self, float rad):
    #     self.cobj.some_var = var
