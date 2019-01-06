import sys
import numpy as np
import cython, os

class PySiggen:

  def __init__(self, geometry_type, conf_file):
    self.conf_file = conf_file
    self.conf_path = os.path.dirname(self.conf_file)

    self.setup = PySetup(self.conf_file.encode('utf-8'))

    if geometry_type == "PPC":
      self.geometry = PPCGeometry(self.setup)
      self.detector = DetectorWrapper_PPC(self.geometry, self.setup)
      self.siggen = SignalGeneratorWrapper_PPC(self.detector, self.setup)
    elif geometry_type == "GEM":
      self.geometry = GEMGeometry(self.setup)
      self.detector = DetectorWrapper_GEM(self.geometry, self.setup)
      self.siggen = SignalGeneratorWrapper_GEM(self.detector, self.setup)
    elif geometry_type == "ICPC":
      self.geometry = ICPCGeometry(self.setup)
      self.detector = DetectorWrapper_ICPC(self.geometry, self.setup)
      self.siggen = SignalGeneratorWrapper_ICPC(self.detector, self.setup)
    else:
      print("Detector geometry {} is not implemented!".format(geometry_type))
      sys.exit()

  @cython.boundscheck(False)
  def MakeSignal(self, float x, float y, float z, double q, np.ndarray[np.float32_t, ndim=1] signal_array not None):
    cdef point pt
    cdef int val
    pt.x=x
    pt.y=y
    pt.z=z

    signal_array.fill(0.)
    val = self.siggen.make_signal(pt, signal_array, q)

    for j in range(1, self.GetNumStepsCalc()):
      signal_array[j] += signal_array[j-1]
    # signal_array = np.cumsum(signal_array)

    return val

  def GetSignal(self, float x, float y, float z, np.ndarray[np.float32_t, ndim=1] signal_array not None):
    cdef point pt
    pt.x=x
    pt.y=y
    pt.z=z

    signal_array.fill(0.)
    val =  self.siggen.get_signal(pt, signal_array)

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

  #Transparent setters/getters

  def UseVeloModel(self, bool usemodel):
    self.detector.set_use_velo_model(usemodel)
  def SetTemp(self, float temp):
    self.detector.set_temp(temp)

  def SetHoles(self,float mu0_100, float beta_100, float E_0_100,
                 float mu0_111, float beta_111, float E_0_111):
       self.detector.set_holes(mu0_100,beta_100,E_0_100,mu0_111,beta_111,E_0_111)
  def SetElectrons(self,float mu0_100, float beta_100, float E_0_100,float mu_n_100,
                float mu0_111, float beta_111, float E_0_111, float mu_n_111):
      self.detector.set_electrons(mu0_100,beta_100,E_0_100,mu_n_100,mu0_111,beta_111,E_0_111, mu_n_111)

  def GetMaxRadius(self):
    return self.geometry.xtal_radius
  def GetMaxZ(self):
    return self.geometry.xtal_length

  def GetNumSegments(self):
    return self.detector.nsegments
  def GetNumSteps(self):
    return self.siggen.output_length

  def GetNumStepsCalc(self):
    return self.siggen.calc_length

  def SetCalcLength(self, int nt):
    self.siggen.calc_length = nt

  def GetCalcTimeStep(self):
    return self.siggen.calc_timestep
  def SetCalcTimestep(self, float dt):
    self.siggen.calc_timestep = dt

  def SetTrapping(self, double trap_const):
    self.detector.trapping_constant = trap_const

  def GetFieldName(self):
    return os.path.join(self.conf_path, self.detector.get_field_name())
  def GetWpotName(self):
    return os.path.join(self.conf_path, self.detector.get_wpot_name())

  def GetDeadLayer(self):
    return self.detector.get_dead_layer()

  #Stuff that actually requires work

  def SetImpurityAvg(self, float imp, float grad):

    if imp + grad * self.GetMaxZ()/2/10 > 0:
      raise ValueError("Values would create a positive impurity at z_max!")

    self.detector.set_impurity_avg( imp,  grad)
  def SetImpurityEnds(self, float imp_0, float imp_max):
    avg = 0.5*(imp_0+imp_max)
    grad = (imp_max - imp_0)/(self.geometry.xtal_length/10)

    self.detector.set_impurity_avg( avg,  grad)

  def SetImpurityZ0(self, float imp, float grad):
    self.detector.set_impurity_z0( imp,  grad)
  def GetImpurity(self):
    return self.detector.get_impurity()
  def GetImpurityGradient(self):
    return self.detector.get_impurity_gradient()
  def GetXtalHV(self):
    return self.detector.get_xtal_HV()

  def SetImpurityZ0(self, float imp, float grad):
    self.detector.set_impurity_z0( imp,  grad)


  def GetLastDriftTime(self, float q):
    return self.siggen.get_last_drifttime(q)

  def GetEfield(self, float r, float z):

    (e_r, e_phi, e_z, val) = self.detector.efield(r,z)
    return (e_r,e_z)

  def GetWpot(self, float r, float z):
    cdef point pt
    pt.x = r;
    pt.y = 0;
    pt.z = z;

    wp_numpy = np.zeros(self.GetNumSegments(), dtype='f4')

    self.detector.wpotential(pt,wp_numpy)

    return wp_numpy

  def GetDWpot(self, float q):
    cdef vector[float] dwpot
    dwpot = self.siggen.get_dwpot()
    nt = self.GetNumStepsCalc()
    dt = self.siggen.get_last_drifttime(q)
    nsegs = self.GetNumSegments()
    dwp_np = np.ones((nsegs,dt))*np.nan

    for i in range(dt):
      for j in range(nsegs):
        dwp_np[j,i] = dwpot[j*nt+i]
    return dwp_np

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

  def save_efield(self, mat_full, efld_name, header_bytes):
    cdef EFieldPoint e_pt
    cdef cyl_pt c_pt
    cdef ofstream* outputter
    cdef unsigned int header_size
    # use try ... finally to ensure destructor is called
    outputter = new ofstream(efld_name.encode(), binary)

    r_num = mat_full.shape[0]
    z_num = mat_full.shape[1]
    imp_num = mat_full.shape[2]
    grad_num = mat_full.shape[3]

    header_size = np.uint32(len(header_bytes))

    outputter.write(<char*> & header_size,4);
    outputter.write(header_bytes, header_size);

    # mat_full = solve_efield()

    for i in range(r_num):
      for j in range(z_num):
        for k in range(imp_num):
          for m in range(grad_num):
            voltage, e, e_r, e_z = mat_full[i,j,k,m,:]
            c_pt.r = e_r
            c_pt.z = e_z
            c_pt.phi = 0
            e_pt.set_field(c_pt)
            e_pt.set_voltage(voltage)
            e_pt.serialize(outputter)

    del outputter
