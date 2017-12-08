#!/usr/local/bin/python
import numpy as np

try:
    from scipy import  signal, ndimage
    from scipy.interpolate import interp1d
except ImportError:
    pass
try:
    from dolfin import *
except ImportError:
    pass

from ._siggen import PySiggen_GEM, PySiggen_PPC, PySiggen_ICPC

#Does all the interfacing with siggen for you, stores/loads lookup tables, and does electronics shaping

class Detector:
  def __init__(self,  detector_geometry, conf_file, num_steps_calc=None, wf_padding=0, maxWfOutputLength = 1000, verbose=False):

    if detector_geometry == "GEM":
        self.siggenInst = PySiggen_GEM(conf_file)#PyICPC(conf_file)
    elif detector_geometry == "ICPC":
        self.siggenInst = PySiggen_ICPC(conf_file)
    elif detector_geometry == "PPC":
        self.siggenInst = PySiggen_PPC(conf_file)
    else:
        print("Detector geometry %s is not implemented!" % detector_geometry)
        exit(0)

    if num_steps_calc is not None:
        self.siggenInst.SetCalcLength(int(num_steps_calc))

    self.siggenInst.InitializeFields()

    self.nsegments=self.siggenInst.GetNumSegments()
    self.detector_radius = self.siggenInst.GetMaxRadius()
    self.detector_length = self.siggenInst.GetMaxZ()

    self.num_steps_out = self.siggenInst.GetNumSteps()

    #this is here so that the gaussian convolution has room on both sides
    self.wf_padding = wf_padding

    self.signal_array = np.zeros( (self.nsegments, wf_padding+self.num_steps_out), dtype='f4', order="C")
    self.signal_array_flat = np.zeros(self.nsegments*self.num_steps_out , dtype='f4', order="C")

    self.num_steps_calc = self.siggenInst.GetNumStepsCalc()
    self.total_steps = self.num_steps_calc + 2*wf_padding

    self.calc_array = np.zeros( (self.nsegments, wf_padding + self.num_steps_calc), dtype='f4', order="C")
    self.calc_array_flat = np.zeros(self.nsegments*self.num_steps_calc , dtype='f4', order="C")

    self.verbose = verbose

    self.time_step_size = self.siggenInst.GetCalcTimeStep()
    digitized_time_step = 10.
    data_to_siggen_size_ratio = np.around( digitized_time_step / self.time_step_size,3)
    self.data_to_siggen_size_ratio = np.int(data_to_siggen_size_ratio)

    self.maxWfOutputLength = maxWfOutputLength
    self.output_wf = np.zeros( maxWfOutputLength, dtype='f4', order="C" )

    self.interpType = 'linear'

  def GetWP(self,nr=500,nz=500):
      wpot = np.zeros((nr,nz))
      for i,r in enumerate(np.linspace(0,self.detector_radius,nr)):
        for j,z in enumerate(np.linspace(0,self.detector_length,nz)):
            wpot[i,j] = self.siggenInst.GetWpot(r,z)[0]
      return wpot

  def GetEfld(self,nr=500,nz=500):
      efld_r = np.zeros((nr,nz))
      efld_z = np.zeros((nr,nz))
      for i,r in enumerate(np.linspace(0,self.detector_radius,nr)):
        for j,z in enumerate(np.linspace(0,self.detector_length,nz)):
            efld_r[i,j], efld_z[i,j] = self.siggenInst.GetEfield(r,z)
      return efld_r, efld_z
  def GetInDetector(self, nr=500,nz=500):
      in_det = np.zeros((nr,nz))
      for i,r in enumerate(np.linspace(0,self.detector_radius,nr)):
        for j,z in enumerate(np.linspace(0,self.detector_length,nz)):
            in_det[i,j] = self.IsInDetector(r,0,z)
      return in_det

  def IsInDetector(self, r, phi, z):
      x = r * np.cos(phi)
      y = r * np.sin(phi)
      return self.siggenInst.InCrystal(x,y,z)

  def SetTransferFunctionRC(self, RC1_in_us, RC2_in_us, rc1_frac, digPeriod  = 1E9):
    RC1= 1E-6 * (RC1_in_us)
    self.rc1_for_tf = np.exp(-1./digPeriod/RC1)

    RC2 = 1E-6 * (RC2_in_us)
    self.rc2_for_tf = np.exp(-1./digPeriod/RC2)

    self.rc1_frac = rc1_frac

    num_term_1 = -1*(self.rc1_for_tf*(1 - self.rc1_frac )  + self.rc2_for_tf*self.rc1_frac + 1)
    num_term_2 = self.rc1_for_tf*(1 - self.rc1_frac )  + self.rc2_for_tf*self.rc1_frac

    self.hp_num = [1., num_term_1, num_term_2]
    self.hp_den = [1, -(self.rc1_for_tf + self.rc2_for_tf), self.rc1_for_tf * self.rc2_for_tf ]

  #makes wf for only one charge type
  def MakeWaveform(self, r, phi,z, q):
    self.calc_array_flat.fill(0)
    self.calc_array.fill(0)

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    calcFlag = self.siggenInst.MakeSignal(x, y, z, q, self.calc_array_flat[:]);

    if calcFlag == -1:
        if self.verbose: print ("Holes out of crystal alert! (%0.3f,%0.3f,%0.3f).  IsInDetector: %d" % (r,phi,z, self.IsInDetector(r,phi,z)))
        return None
    self.calc_array[:,self.wf_padding:] = self.calc_array_flat.reshape((self.nsegments,self.num_steps_calc))
    return self.calc_array

  #summed waveform of both charge types
  def GetWaveform(self, r,phi,z, energy=1):
    self.signal_array_flat.fill(0)
    self.signal_array.fill(0)

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    calcFlag = self.siggenInst.GetSignal(x, y, z, self.signal_array_flat[:]);

    if calcFlag == -1:
      if self.verbose: print ("Holes out of crystal alert! (%0.3f,%0.3f,%0.3f).  IsInDetector: %d" % (r,phi,z, self.IsInDetector(r,phi,z)))
      return None

    self.signal_array[:,self.wf_padding:] = self.signal_array_flat.reshape((self.nsegments,self.num_steps_out))
    return self.signal_array

# ###########################################################################################################################
#   def GetSimWaveform(self, r,phi,z,scale, align_point,  numSamples, smoothing=None):
#     sig_wf = self.GetRawSiggenWaveform(r, phi, z)
#     if sig_wf is None:
#       return None
#
#     sig_wf *= scale
#
#     #only do gaussian filtering of PC signal
#     if smoothing is not None:
#       ndimage.filters.gaussian_filter1d(sig_wf[self.pc_idx,:], smoothing, output=sig_wf[self.pc_idx,:])
#
#     sim_wf = self.ProcessWaveform(sig_wf, align_point,  numSamples)
#
#     return sim_wf
#
# ########################################################################################################
#   def ProcessWaveform(self, siggen_wf,  align_point, outputLength):
#     #low-pass filter (non-ppc)
#     temp_wf_sig = np.zeros_like(siggen_wf)
#     for i in range(self.nsegments):
#         temp_wf_sig[i,:] = signal.lfilter(self.lp_num[i], self.lp_den[i], siggen_wf[i,:])
#         # print(temp_wf_sig[i,:])
#         # print(self.lp_num[i] ,self.lp_den[i], np.sum(self.lp_den[i]))
#         temp_wf_sig[i,:] /= (np.sum(self.lp_num[i])/np.sum(self.lp_den[i]))
#
#     #RC filter for core
#     if self.core_rc > 0:
#         temp_wf_sig[-1,:] = signal.lfilter(self.core_num, self.core_den, siggen_wf[-1,:])
#     #hi-pass filter
#     temp_wf_sig= signal.lfilter(self.hp_num, self.hp_den, temp_wf_sig, axis=-1)
#
#     #uh, downsample it
#     temp_wf = temp_wf_sig[:,::self.data_to_siggen_size_ratio]
#
#     #enforces t0 time-alignment to the align_point
#     siggen_offset = self.t0_padding/self.data_to_siggen_size_ratio
#     first_idx= 0
#
#     num_wfs = self.nsegments
#     num_samples = self.total_steps/self.data_to_siggen_size_ratio
#
#     #find the first index after the align point
#     align_point_ceil = np.int( np.ceil(align_point) )
#     start_idx = align_point_ceil - first_idx
#     #and make sure its above 0!
#     if start_idx <0:
#         print("bad start idx")
#         return None
#
#     siggen_interp_fn = interpolate.interp1d(np.arange(num_samples), temp_wf, kind=self.interpType, copy="False", assume_sorted="True")
#
#     num_samples_to_fill = outputLength - start_idx
#     offset = align_point_ceil - align_point
#     sampled_idxs = np.arange(num_samples_to_fill) + offset + siggen_offset
#
#     processed_siggen_data = np.zeros((num_wfs, outputLength))
#     coarse_vals =   siggen_interp_fn(sampled_idxs)
#
#     try:
#         processed_siggen_data[:, start_idx:start_idx+num_samples_to_fill] = coarse_vals
#     except ValueError:
#         print( len(processed_siggen_data) )
#         print( start_idx)
#         print( num_samples_to_fill)
#         print( sampled_idxs)
#         exit(0)
#
#     return processed_siggen_data

  def solve_fields(self, meshmult, xtal_HV, impAvgRange, gradientRange, wp_name = "wpot.field", ef_name="ev.field", num_cpu=1):

      from multiprocessing import Pool, cpu_count

      max_cpu = cpu_count()
      if (num_cpu < 1 or num_cpu >  max_cpu): cpu_count = max_cpu

      boundary_pc, boundary_n = self.GetBoundaryConditions()

      print("Solving WP...")
      wp_mat = self.solve_wp(meshmult,boundary_pc, boundary_n)
      # with open(wpot_name, 'wb') as out:
      wp_mat.tofile(wp_name)

      nr = int(self.detector_radius*meshmult+1)
      nz = int(self.detector_length*meshmult+1)
      ngrad = len(gradientRange)
      nimp = len(impAvgRange)

      efield_args = []

      efield = np.zeros((nr,nz,nimp,ngrad,4), dtype=np.object)
      for i,avg in enumerate(impAvgRange):
        for j,grad in enumerate(gradientRange):
            print("Solving EF {} of {}... imp {}, grad {}".format(i*ngrad +j + 1, ngrad*nimp, avg, grad))
            if num_cpu  == 1:
                efield[:,:,i,j,:] = self.solve_efield(meshmult, xtal_HV, avg, grad, boundary_pc, boundary_n)
            else:
                efield_args.append[ (meshmult, xtal_HV, avg, grad, boundary_pc, boundary_n)  ]

      self.siggenInst.save_efield(efield, ef_name)

  def avg_to_z0(self, impurity_avg, impurity_gradient):
    return impurity_avg - impurity_gradient * (self.detector_length/10) / 2
  def z0_to_avg(self, impurity_z0, impurity_gradient):
    return impurity_z0 + impurity_gradient * (self.detector_length/10) / 2

  def solve_wp(self, meshmult, boundary_pc, boundary_n, num_threads=1):
    # parameters.add('num_threads', num_threads)
    # parameters['num_threads'] = num_threads

    # Create mesh and define function space
    mesh = RectangleMesh(Point(0.0, 0.0), Point(self.detector_radius+0.1, self.detector_length+0.1), int((self.detector_radius+0.1)*meshmult+1), int((self.detector_length+0.1)*meshmult+1), "crossed")
    V = FunctionSpace(mesh, "Lagrange", 1)

    # Define boundary condition
    u0 = Constant(0.0)
    u1 = Constant(1.0)

    bc0 = DirichletBC(V, u0, boundary_n)
    bc1 = DirichletBC(V, u1, boundary_pc)

    # Define variational problem
    T = TrialFunction(V)
    q = TestFunction(V)

    r = Expression('x[0]', degree=1)
    a = (Dx(T,0)*Dx(q,0) + Dx(T,1)*Dx(q,1))*r*dx()#inner(grad(u), grad(v))*dx

    f = Constant(0.0)
    g = Constant(0.0)
    L = f*q*dx + g*q*ds

    # Compute solution
    u = Function(V)
    solve(a == L, u, [bc0,bc1])

    nr = int(self.detector_radius*meshmult+1)
    nz = int(self.detector_length*meshmult+1)
    mat = np.zeros((nr,nz), dtype=np.float32)
    for i,r in enumerate(np.linspace(0, self.detector_radius, self.detector_radius*meshmult+1)):
        for j,z in enumerate(np.linspace(0, self.detector_length, self.detector_length*meshmult+1)):
            wp = u(Point(r,z))
            mat[i,j] = wp
    return mat

  def solve_efield(self, meshmult, xtal_HV, impurity_avg, impurity_gradient, boundary_pc, boundary_n, impurity_surface=0, num_threads=1):
    #   print ("{}".format(parameters))
    #   parameters['num_threads'] = num_threads

      mesh = RectangleMesh(Point(0.0, 0.0), Point(self.detector_radius+0.1, self.detector_length+0.1), int((self.detector_radius+0.1)*meshmult+1), int((self.detector_length+0.1)*meshmult+1), "crossed")
      V = FunctionSpace(mesh, "Lagrange", 1)

      # Define boundary condition
      u0 = Constant(xtal_HV)
      u1 = Constant(0.)

      bc0 = DirichletBC(V, u0, boundary_n)
      bc1 = DirichletBC(V, u1, boundary_pc)

      # Define variational problem
      T = TrialFunction(V)
      q = TestFunction(V)

      r = Expression('x[0]', degree=1)
      a = (Dx(T,0)*Dx(q,0) + Dx(T,1)*Dx(q,1))*r*dx()#inner(grad(u), grad(v))*dx

      #Need to express rho/epsilon in units of V / mm^2

      #For germanium, dielectric constant is 16, so epsilon = 16*epsilon_0
      epsilon_0 = 8.854187817E-15 #F/mm or C/(V*mm)
      epsilon_ge = 16*epsilon_0

      #impurity numbers are expressed in eE10/cm^3, so the charge density is like
      e = 1.602176E-19 #C/charge
      e_per_vol = e * 1E10 * 1E-3 #charge density in units of C/mm^3

      e_per_surf = e * 1E10 * 1E-2

      e_over_eps = e_per_vol / epsilon_ge #units of V/mm^2 #

      #we can multiply this by the density number as it comes from the conf file now

      #have to go from a volume density to a surface density.  factor of

      impurity_z0 = self.avg_to_z0(impurity_avg, impurity_gradient)
      if impurity_surface == 0:
          expstring = "({} + {}*x[1])*{}  ".format(impurity_z0,0.1*impurity_gradient, e_over_eps)
      else:
          expstring = "(x[1] < 0.01) ? {}  :({} + {}*x[1])*{}  ".format(impurity_surface*e_per_surf/epsilon_ge, impurity_z0,0.1*impurity_gradient, e_over_eps)

      #Constant(-5)#
      f = Expression(expstring, degree=1)
      g = Constant(0)
      L = f*q*r*dx() + g*q*ds

      # Compute solution
      u = Function(V)
      solve(a == L, u, [bc0,bc1])
      print("...done solving")

      V_vec = VectorFunctionSpace(mesh, "CG", 1)
      gradu = project(grad(u),V_vec)

      nr = int(self.detector_radius*meshmult+1)
      nz =int(self.detector_length*meshmult+1)
      mat_full = np.zeros((nr,nz,4), dtype=np.float32)

      for i,r in enumerate(np.linspace(0, self.detector_radius, self.detector_radius*meshmult+1)):
          for j,z in enumerate(np.linspace(0, self.detector_length, self.detector_length*meshmult+1)):
              voltage = u(Point(r,z))
              # mat[i,j] = voltage
              gradv = gradu(Point(r,z)) #presumably in V/mm
              e_r = -gradv[0] * 10 #to V/cm
              e_z = -gradv[1] * 10 #to V/cm
              e_theta = 0.
              e = np.sqrt(e_r**2+e_z**2)
              mat_full[i,j,:] = voltage, e, e_r, e_z

      return mat_full
