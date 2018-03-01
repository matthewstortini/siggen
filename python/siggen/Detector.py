#!/usr/local/bin/python
import numpy as np
import xml.etree.ElementTree as ET
import struct, os
try:
    from scipy import  signal, ndimage
    from scipy.interpolate import interp1d
except ImportError:
    pass
try:
    from dolfin import *
except ImportError:
    pass

from ._siggen import PySiggen

#Does all the interfacing with siggen for you, stores/loads lookup tables, and does electronics shaping

class Detector:
  def __init__(self,  detector_geometry, conf_file, num_steps_calc=None, wf_padding=0, maxWfOutputLength = 1000, verbose=False, doInit=True):

    self.siggenInst = PySiggen(detector_geometry, conf_file)
    self.conf_file = conf_file

    if num_steps_calc is not None:
        print("setting calc length in det")
        self.siggenInst.SetCalcLength(int(num_steps_calc))

    if doInit==True:
        self.siggenInst.InitializeFields()
        self.read_impurity_range()

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

  def SetTransferFunctionRC(self, RC1_in_us, RC2_in_us, rc1_frac, digFrequency  = 1E9):
    RC1= 1E-6 * (RC1_in_us)
    self.rc1_for_tf = np.exp(-1./digFrequency/RC1)

    RC2 = 1E-6 * (RC2_in_us)
    self.rc2_for_tf = np.exp(-1./digFrequency/RC2)

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
    return self.signal_array*energy

  def read_impurity_range(self):
      efield_file = self.siggenInst.GetFieldName()
      efield_file_ext = efield_file.split(".")[-1]
      if efield_file_ext != "field":
          print("Error: can currently only read impurity information from binary field files")
          return

      conf_dir, __ = os.path.split(self.conf_file)
      efield_file_path = os.path.join(conf_dir, efield_file)

      with open(efield_file_path, "rb") as fin:
          header_len, = struct.unpack('i', fin.read(4))
          header = fin.read(header_len)
          root = ET.fromstring(header.decode("ascii"))

          for var in root.findall('variable'):
              name = var.get('name')
              if name == "impurity_avg":
                  self.imp_avg_lims = [ float(var.find("min").text), float(var.find("max").text) ]
              elif name == "impurity_grad":
                  self.imp_grad_lims = [ float(var.find("min").text), float(var.find("max").text) ]
      return

  def solve_fields(self, meshmult, impAvgRange, gradientRange, wp_name = "wpot.field", ef_name="ev.field"):

    nr = int(self.detector_radius*meshmult+1)
    nz = int(self.detector_length*meshmult+1)
    ngrad = len(gradientRange)
    nimp = len(impAvgRange)

    #create some self-descriptive xml
    header_tree = ET.Element('field')

    r = ET.SubElement(header_tree, 'variable')
    r.set("name", "radialDimension")
    ET.SubElement(r, 'min').text = "{}".format(0)
    ET.SubElement(r, 'max').text = "{}".format(self.detector_radius)
    ET.SubElement(r, 'num').text = "{}".format(nr)

    l = ET.SubElement(header_tree, 'variable')
    l.set("name", "axialDimension")
    ET.SubElement(l, 'min').text = "{}".format(0)
    ET.SubElement(l, 'max').text = "{}".format(self.detector_length)
    ET.SubElement(l, 'num').text = "{}".format(nz)

    #No impurity info for WP
    header_bytes = ET.tostring(header_tree, encoding='ascii', method='xml')

    print("Solving WP...")
    wp_mat = self.solve_field("wpot", nr, nz)

    with open(wp_name, 'wb') as out:
        out.write(np.int32(len(header_bytes)))
        out.write(header_bytes)
        out.write(wp_mat.tobytes())

    #add impurity info for E field
    i = ET.SubElement(header_tree, 'variable')
    i.set("name", "impurity_avg")
    ET.SubElement(i, 'min').text = "{}".format(impAvgRange[0])
    ET.SubElement(i, 'max').text = "{}".format(impAvgRange[-1])
    ET.SubElement(i, 'num').text = "{}".format(len(impAvgRange))
    ET.SubElement(i, 'num').text = "{}".format(len(impAvgRange))

    g = ET.SubElement(header_tree, 'variable')
    g.set("name", "impurity_grad")
    ET.SubElement(g, 'min').text = "{}".format(gradientRange[0])
    ET.SubElement(g, 'max').text = "{}".format(gradientRange[-1])
    ET.SubElement(g, 'num').text = "{}".format(len(gradientRange))
    ET.SubElement(g, 'num').text = "{}".format(len(gradientRange))

    header_bytes = ET.tostring(header_tree, encoding='ascii', method='xml')

    efield = np.zeros((nr,nz,nimp,ngrad,4), dtype=np.float32)
    for i,avg in enumerate(impAvgRange):
        for j,grad in enumerate(gradientRange):
            print("Solving EF {} of {}... imp {}, grad {}".format(i*ngrad +j + 1, ngrad*nimp, avg, grad))
            efield[:,:,i,j,:] = self.solve_field("efield", nr, nz, impurity_gradient=grad, impurity_avg=avg)

    self.siggenInst.save_efield(efield, ef_name, header_bytes)
    return (wp_mat, efield)

  def avg_to_z0(self, impurity_avg, impurity_gradient):
    return impurity_avg - impurity_gradient * (self.detector_length/10) / 2
  def z0_to_avg(self, impurity_z0, impurity_gradient):
    return impurity_z0 + impurity_gradient * (self.detector_length/10) / 2
  #
  # def solve_wp(self, meshmult, boundary_pc, boundary_n, num_threads=1):
  #   # parameters.add('num_threads', num_threads)
  #   # parameters['num_threads'] = num_threads
  #
  #   # Create mesh and define function space
  #   mesh = RectangleMesh(Point(0.0, 0.0), Point(self.detector_radius+0.1, self.detector_length+0.1), int((self.detector_radius+0.1)*meshmult+1), int((self.detector_length+0.1)*meshmult+1), "crossed")
  #   V = FunctionSpace(mesh, "Lagrange", 1)
  #
  #   # Define boundary condition
  #   u0 = Constant(0.0)
  #   u1 = Constant(1.0)
  #
  #   bc0 = DirichletBC(V, u0, boundary_n)
  #   bc1 = DirichletBC(V, u1, boundary_pc)
  #
  #   # Define variational problem
  #   T = TrialFunction(V)
  #   q = TestFunction(V)
  #
  #   r = Expression('x[0]', degree=1)
  #   a = (Dx(T,0)*Dx(q,0) + Dx(T,1)*Dx(q,1))*r*dx()#inner(grad(u), grad(v))*dx
  #
  #   f = Constant(0.0)
  #   g = Constant(0.0)
  #   L = f*q*dx + g*q*ds
  #
  #   # Compute solution
  #   u = Function(V)
  #   solve(a == L, u, [bc0,bc1])
  #
  #   nr = int(self.detector_radius*meshmult+1)
  #   nz = int(self.detector_length*meshmult+1)
  #   mat = np.zeros((nr,nz), dtype=np.float32)
  #   for i,r in enumerate(np.linspace(0, self.detector_radius, self.detector_radius*meshmult+1)):
  #       for j,z in enumerate(np.linspace(0, self.detector_length, self.detector_length*meshmult+1)):
  #           wp = u(Point(r,z))
  #           mat[i,j] = wp
  #   return mat
  #
  # def solve_efield(self, meshmult, xtal_HV, impurity_avg, impurity_gradient, boundary_pc, boundary_n, impurity_surface=0, num_threads=1):
  #   #   print ("{}".format(parameters))
  #   #   parameters['num_threads'] = num_threads
  #
  #   #   class e_over_epsilon(Expression):
  #   #     def eval(self, value, x):
  #   #         "Set value[0] to value at point x"
  #   #         tol = 1E-14
  #   #         if x[1] <= 0.5 + tol:
  #   #             value[0] = self.k_0
  #   #         else:
  #   #             value[0] = self.k_1
  #
  #     mesh = RectangleMesh(Point(0.0, 0.0), Point(self.detector_radius+0.1, self.detector_length+0.1), int((self.detector_radius+0.1)*meshmult+1), int((self.detector_length+0.1)*meshmult+1), "crossed")
  #     V = FunctionSpace(mesh, "Lagrange", 1)
  #
  #     # Define boundary condition
  #     u0 = Constant(xtal_HV)
  #     u1 = Constant(0.)
  #
  #     bc0 = DirichletBC(V, u0, boundary_n)
  #     bc1 = DirichletBC(V, u1, boundary_pc)
  #
  #     # Define variational problem
  #     T = TrialFunction(V)
  #     q = TestFunction(V)
  #
  #     r = Expression('x[0]', degree=1)
  #     a = (Dx(T,0)*Dx(q,0) + Dx(T,1)*Dx(q,1))*r*dx()#inner(grad(u), grad(v))*dx
  #
  #     #Need to express rho/epsilon in units of V / mm^2
  #
  #     #For germanium, dielectric constant is 16, so epsilon = 16*epsilon_0
  #     epsilon_0 = 8.854187817E-15 #F/mm or C/(V*mm)
  #     epsilon_ge = 16*epsilon_0
  #
  #     #impurity numbers are expressed in eE10/cm^3, so the charge density is like
  #     e = 1.602176E-19 #C/charge
  #     e_per_vol = e * 1E10 * 1E-3 #charge density in units of C/mm^3
  #
  #     e_per_surf = e * 1E10 * 1E-2
  #
  #     e_over_eps = e_per_vol / epsilon_ge #units of V/mm^2 #
  #
  #     #we can multiply this by the density number as it comes from the conf file now
  #
  #     #have to go from a volume density to a surface density.  factor of
  #
  #     impurity_z0 = self.avg_to_z0(impurity_avg, impurity_gradient)
  #     if impurity_surface == 0:
  #         expstring = "({} + {}*x[1])*{}  ".format(impurity_z0,0.1*impurity_gradient, e_over_eps)
  #     else:
  #         expstring = "(x[1] < 0.01) ? {}  :({} + {}*x[1])*{}  ".format(impurity_surface*e_per_surf/epsilon_ge, impurity_z0,0.1*impurity_gradient, e_over_eps)
  #
  #     #Constant(-5)#
  #     f = Expression(expstring, degree=1)
  #     g = Constant(0)
  #     L = f*q*r*dx() + g*q*ds
  #
  #     # Compute solution
  #     u = Function(V)
  #     solve(a == L, u, [bc0,bc1])
  #     print("...done solving")
  #
  #     V_vec = VectorFunctionSpace(mesh, "CG", 1)
  #     gradu = project(grad(u),V_vec)
  #
  #     nr = int(self.detector_radius*meshmult+1)
  #     nz =int(self.detector_length*meshmult+1)
  #     mat_full = np.zeros((nr,nz,4), dtype=np.float32)
  #
  #     for i,r in enumerate(np.linspace(0, self.detector_radius, self.detector_radius*meshmult+1)):
  #         for j,z in enumerate(np.linspace(0, self.detector_length, self.detector_length*meshmult+1)):
  #             voltage = u(Point(r,z))
  #             # mat[i,j] = voltage
  #             gradv = gradu(Point(r,z)) #presumably in V/mm
  #             e_r = -gradv[0] * 10 #to V/cm
  #             e_z = -gradv[1] * 10 #to V/cm
  #             e_theta = 0.
  #             e = np.sqrt(e_r**2+e_z**2)
  #             mat_full[i,j,:] = voltage, e, e_r, e_z
  #
  #   #   import matplotlib.pyplot as plt
  #   #   plt.figure(figsize=(14,6))
  #   #   plt.subplot(131)
  #   #   plt.imshow(mat_full[:,:,1].T, origin="lower")
  #   #   plt.subplot(132)
  #   #   plt.imshow(mat_full[:,:,2].T, origin="lower")
  #   #   plt.subplot(133)
  #   #   plt.imshow(mat_full[:,:,3].T, origin="lower")
  #   #   plt.show()
  #   #   exit()
  #
  #     return mat_full
