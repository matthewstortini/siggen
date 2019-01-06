from cython.operator cimport dereference
import cython
import sys

cdef extern from "PPC.h":
  cdef cppclass PPC:
    PPC(Setup setup);
    float get_xtal_radius()
    float get_xtal_length()
    float get_pc_length()
    float get_pc_radius()
    bool get_bulletize_PC()
    float get_top_bullet_radius()
    float get_bottom_bullet_radius()
    float get_wrap_around_radius()
    float get_ditch_depth()
    float get_ditch_thickness()
    float get_taper_length()


cdef class PPCGeometry:
  cdef PPC* cobj

  def __cinit__(self, PySetup setup):
    self.cobj = new PPC(dereference(setup.cobj))

  def __dealloc__(self):
    del self.cobj

  property xtal_radius:
    def __get__(self):
      return self.cobj.get_xtal_radius()
    def __set__(self, float val):
      print("set not implemented for xtal_radius!")
      sys.exit()

  property xtal_length:
    def __get__(self):
      return self.cobj.get_xtal_length()
    def __set__(self, float val):
      print("set not implemented for xtal_length!")
      sys.exit()
  property pc_length:
    def __get__(self):
      return self.cobj.get_pc_length()
    def __set__(self, float val):
      print("set not implemented for pc_length!")
      sys.exit()

  property pc_radius:
    def __get__(self):
      return self.cobj.get_pc_radius()
    def __set__(self, float val):
      print("set not implemented for pc_radius!")
      sys.exit()
  property bulletize_PC:
    def __get__(self):
      return self.cobj.get_bulletize_PC()
    def __set__(self, float val):
      print("set not implemented for bulletize_PC!")
      sys.exit()

  property top_bullet_radius:
    def __get__(self):
      return self.cobj.get_top_bullet_radius()
    def __set__(self, float val):
      print("set not implemented for top_bullet_radius!")
      sys.exit()

  property bottom_bullet_radius:
    def __get__(self):
      return self.cobj.get_bottom_bullet_radius()
    def __set__(self, float val):
      print("set not implemented for bottom_bullet_radius!")
      sys.exit()

  property wrap_around_radius:
    def __get__(self):
      return self.cobj.get_wrap_around_radius()
    def __set__(self, float val):
      print("set not implemented for wrap_around_radius!")
      sys.exit()

  property ditch_depth:
    def __get__(self):
      return self.cobj.get_ditch_depth()
    def __set__(self, float val):
      print("set not implemented for get_ditch_depth!")
      sys.exit()

  property ditch_thickness:
    def __get__(self):
      return self.cobj.get_ditch_thickness()
    def __set__(self, float val):
      print("set not implemented for ditch_thickness!")
      sys.exit()

  property taper_length:
    def __get__(self):
      return self.cobj.get_taper_length()
    def __set__(self, float val):
      print("set not implemented for taper_length!")
      sys.exit()
