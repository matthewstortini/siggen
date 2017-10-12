from cython.operator cimport dereference
import numpy as np
cimport numpy as np
import cython

cdef class PySiggen_GEM:
  cdef Setup* setup
  cdef Detector[GEM]* detector
  cdef SignalGenerator[GEM]* siggen

  def __cinit__(self, conf_file):
    self.setup = new Setup(conf_file.encode('utf-8'))
    self.detector = new Detector[GEM](dereference(self.setup))
    self.siggen = new SignalGenerator[GEM](self.detector, dereference(self.setup))

  def __dealloc__(self):
    del self.setup
    del self.detector
    del self.siggen

  @cython.boundscheck(False)
  def MakeSignal(self, float x, float y, float z, float q, np.ndarray[np.float32_t, ndim=1] signal_array not None):
    cdef point pt
    cdef int val
    pt.x=x
    pt.y=y
    pt.z=z

    signal_array.fill(0.)
    val = self.siggen.make_signal(pt, &signal_array[0], q)

    for j in range(1, self.GetNumStepsCalc()):
      signal_array[j] += signal_array[j-1]

    return val


  def GetSignal(self, float x, float y, float z, np.ndarray[np.float32_t, ndim=1] signal_array not None):
    cdef point pt
    pt.x=x
    pt.y=y
    pt.z=z

    signal_array.fill(0.)
    val =  self.siggen.get_signal(pt, &signal_array[0])

    return val

  def InitializeFields(self):
    self.detector.field_setup()

  def InCrystal(self,  float x, float y, float z):
    cdef point pt
    cdef int result
    pt.x=x
    pt.y=y
    pt.z=z
    result = self.detector.outside_detector(pt)
    return not result

  def SetHoles(self,float mu0_100, float beta_100, float E_0_100,
                 float mu0_111, float beta_111, float E_0_111):
       self.detector.set_holes(mu0_100,beta_100,E_0_100,mu0_111,beta_111,E_0_111)
  def SetElectrons(self,float mu0_100, float beta_100, float E_0_100,float mu_n_100,
                float mu0_111, float beta_111, float E_0_111, float mu_n_111):
      self.detector.set_electrons(mu0_100,beta_100,E_0_100,mu_n_100,mu0_111,beta_111,E_0_111, mu_n_111)

  def GetMaxRadius(self):
    return self.setup.xtal_radius
  def GetMaxZ(self):
    return self.setup.xtal_length
  def GetNumSegments(self):
    return self.detector.get_nsegments()
  def GetNumSteps(self):
    return self.siggen.get_output_length()
  def GetNumStepsCalc(self):
    return self.siggen.get_calc_length()
  def GetLastDriftTime(self, float q):
    return self.siggen.get_last_drifttime(q)

  def GetEfield(self, float r, float z):
    cdef cyl_pt e
    cdef cyl_pt pt

    e.r=0.
    e.phi=0.
    e.z=0.

    pt.r = r;
    pt.phi=0.;
    pt.z=z;

    self.detector.efield(pt,e)
    return (e.r,e.z)

  def GetWpot(self, float r, float z):
    cdef vector[float] wp
    wp.resize(self.GetNumSegments())
    cdef point pt
    pt.x = r;
    pt.y = 0;
    pt.z = z;

    wp_numpy = np.zeros(self.GetNumSegments())

    self.detector.wpotential(pt,wp)
    for i in range(wp.size()):
      wp_numpy[i] = wp.at(i)
    return wp_numpy

  def GetPath(self, float q):
    cdef vector[point] dp
    dp = self.siggen.get_driftpath(q)
    dt = self.siggen.get_last_drifttime(q)
    if dt > dp.size():
      print("dt %d, size %d, how did this happen?" % (dt,dp.size()))

    #turn it into a np array
    dp_np = np.ones((3,dt))*np.nan
    for i in range(dt):
      if dp.at(i).x == 0 and dp.at(i).y==0 and dp.at(i).z==0:
        continue
      else: dp_np[:,i] = dp.at(i).x, dp.at(i).y, dp.at(i).z

    return dp_np
