#!/usr/local/bin/python
import sys
import numpy as np
try:
    from scipy import  signal, ndimage
    from scipy.interpolate import interp1d
    from scipy import special
except ImportError:
    pass
try:
    from dolfin import *
except ImportError:
    pass

from ..Detector import Detector

#Does all the interfacing with siggen for you, stores/loads lookup tables, and does electronics shaping

class PPC(Detector):
  def __init__(self, conf_file, **kwargs):
    super(PPC, self).__init__("PPC", conf_file, **kwargs)

    (self.detector_radius, self.detector_length) = np.floor( [self.detector_radius*10, self.detector_length*10] )/10.

    #TODO: figure out how to write a getter for this
    self.taper_length = 4.5
    self.top_bullet_radius = 1.2

    self.padded_siggen_data = np.zeros(self.total_steps, dtype='f4', order="C")

    self.smoothing_type = 0 #gaussian

    self.digital_filters = []

  def GetWaveform(self, r, phi, z, energy=1):
      #Overload Detector method to return a 1-D vector (since we only have one electrode)
      wf = super().GetWaveform(r,phi,z,energy)
      if wf is None: return wf
      return wf[0,:]

  def GetRadLims(self,theta):
      theta_eq = np.arctan(self.detector_length/self.detector_radius)
      theta_taper = np.arctan(self.taper_length/self.detector_radius)
      if theta <= theta_taper:
         z = np.tan(theta)*(self.detector_radius - self.taper_length) / (1-np.tan(theta))
         max_rad = z / np.sin(theta)
      elif theta <= theta_eq:
          max_rad = self.detector_radius / np.cos(theta)
      else:
          theta_comp = np.pi/2 - theta
          max_rad = self.detector_length / np.cos(theta_comp)

      #AND THE MINIMUM (from PC dimple)
      #min_rad  = 1./ ( np.cos(theta)**2/self.pcRad**2  +  np.sin(theta)**2/self.pcLen**2 )

      min_rad = np.amax([self.pcRad, self.pcLen])

      return (min_rad, max_rad)

  def GetThetaLims(self,rad):
    if rad < np.amin([self.detector_radius - self.taper_length, self.detector_length]):
        max_val = np.pi/2
        min_val = 0
    else:
        if rad < self.detector_radius - self.taper_length:
            #can't possibly hit the taper
            min_val = 0
        elif rad < np.sqrt(self.detector_radius**2 + self.taper_length**2):
            #low enough that it could hit the taper region
            a = self.detector_radius - self.taper_length
            z = 0.5 * (np.sqrt(2*rad**2-a**2) - a)
            min_val = np.arcsin(z/rad)
        else:
            #longer than could hit the taper
            min_val = np.arccos(self.detector_radius/rad)

        if rad < self.detector_length:
            max_val = np.pi/2
        else:
            max_val = np.pi/2 - np.arccos(self.detector_length/rad)
    return (min_val, max_val)


# ###########################################################################################################################
#   def GetSimWaveform(self, r,phi,z,scale, align_point, align_percent, numSamples, smoothing=None):
#     sig_wf = self.GetWaveform(r, phi, z)
#     if sig_wf is None:
#       return None
#     sig_wf *= scale
#
#     if smoothing is not None:
#       ndimage.filters.gaussian_filter1d(sig_wf[0,:], smoothing, output=sig_wf[0,:])
#
#     sim_wf = self.ProcessWaveform(sig_wf[0,:], align_point,align_percent, numSamples)

###########################################################################################################################
  def MakeSimWaveform(self, r,phi,z,scale, align_point, align_percent, numSamples, smoothing=None, e_smoothing=None, skew=None):
    hole_wf = self.MakeWaveform(r, phi, z, 1)
    if hole_wf is None:
      return None

    self.padded_siggen_data.fill(0.)
    self.padded_siggen_data[:len(hole_wf[0,:])] = hole_wf[0,:]
    self.padded_siggen_data[len(hole_wf[0,:]):] = hole_wf[0,-1]

    if smoothing is not None:
      if self.smoothing_type == 0:
          ndimage.filters.gaussian_filter1d(self.padded_siggen_data, smoothing, output=self.padded_siggen_data, mode="nearest")
      elif self.smoothing_type == 1:
          window = skew_norm(np.linspace(-200,100,301),  a=skew, w=smoothing )

          pad = len(window)
          wf_pad = np.pad(self.padded_siggen_data, (pad,pad), 'constant', constant_values=(0, self.padded_siggen_data[-1]))
          wf_pad= signal.convolve(wf_pad, window, 'same')
          self.padded_siggen_data = wf_pad[pad:-pad]

    electron_wf = self.MakeWaveform(r, phi, z, -1)

    # if e_smoothing is not None:
    #   ndimage.filters.gaussian_filter1d(electron_wf[0,:], e_smoothing, output=electron_wf[0,:], mode="nearest")

    if electron_wf is  None:
      return None
    self.padded_siggen_data[:len(electron_wf[0,:])] += electron_wf[0,:]
    self.padded_siggen_data[len(electron_wf[0,:]):] = self.padded_siggen_data[len(electron_wf[0,:])-1]

    self.padded_siggen_data *= scale

    # ret = self.padded_siggen_data[::10]
    # return ret[np.argmax(ret)-(numSamples-50): np.argmax(ret)+50]

    sim_wf = self.ProcessWaveform(self.padded_siggen_data, align_point,align_percent, numSamples)

    return sim_wf

########################################################################################################
  def AddDigitalFilter(self, filter):
      self.digital_filters.append(filter)

  def RemoveDigitalFilter(self, filter):
      self.digital_filters.remove(filter)

########################################################################################################
  def ProcessWaveform(self, siggen_wf,  align_point, align_percent, outputLength):
    interpType = "linear"

    data_to_siggen_size_ratio = self.data_to_siggen_size_ratio

    #we want to make sure we have at least enough from zero padding to end to give a "max wf"
    min_length = self.maxWfOutputLength * data_to_siggen_size_ratio + self.wf_padding
    temp_len = len(siggen_wf) if len(siggen_wf) > min_length else min_length
    #and need enough zero padding to give us the t0...
    extra_t0_pad = np.int(align_point*data_to_siggen_size_ratio)

    # if (self.maxWfOutputLength * data_to_siggen_size_ratio - (self.num_steps_calc + self.wf_padding) )

    temp_wf_sig = np.zeros( np.int(temp_len+extra_t0_pad))
    temp_wf_sig[extra_t0_pad:len(siggen_wf)+extra_t0_pad] = siggen_wf
    temp_wf_sig[extra_t0_pad+len(siggen_wf):] = siggen_wf[-1]

    for filter in self.digital_filters:
        temp_wf_sig = filter.apply_to_signal(temp_wf_sig)

    smax_idx = np.argmax(temp_wf_sig)
    smax = temp_wf_sig[smax_idx]
    if smax == 0:
    #   print("bad smax")
      return None

    #linear interpolation to find the alignPointIdx: find the align percentage in the simualted array
    alignarr = np.copy(temp_wf_sig)/smax
    first_idx = np.argmax(alignarr[:smax_idx] > align_percent) - 1

    if first_idx+1 == len(alignarr) or first_idx <0:
        # print("bad first idx")
        return None

    #linear interpolation as to where the true siggen align_percentage is (in the siggen wf)
    slope = (alignarr[first_idx+1] - alignarr[first_idx]) / 1.
    siggen_offset = ( align_percent -  alignarr[first_idx] ) / slope
    siggen_align = siggen_offset + first_idx

    #where (in the data wf) do we want to align?
    align_point_ceil = np.int( np.ceil(align_point) )

    start_idx_sig = siggen_align - align_point*data_to_siggen_size_ratio

    #align_point_ceil*data_to_siggen_size_ratio - (first_idx - self.t0_padding)
    #
    # if start_idx_sig <0:
    #     # print("bad start idx")
    #     return None

    #TODO: i only really _need_ to interp between like every data_to_siggen_size_ratio sample or something right?
    self.siggen_interp_fn = interp1d(np.arange(len(temp_wf_sig)), temp_wf_sig, kind=interpType, copy="False", assume_sorted="True")

    num_samples_to_fill = outputLength
    offset = (align_point_ceil - align_point)*data_to_siggen_size_ratio

    sampled_idxs = np.arange(num_samples_to_fill)*data_to_siggen_size_ratio + start_idx_sig #+ siggen_offset #+ offset +

    if sampled_idxs[0] < 0 or sampled_idxs[-1] > len(temp_wf_sig) - 1: return None
    # print( len(self.output_wf), outputLength )
    # print (offset , siggen_offset)
    # print (align_point_ceil, first_idx, data_to_siggen_size_ratio)
    # print (siggen_align, align_point, start_idx_sig, len(temp_wf_sig), )
    # print( num_samples_to_fill )
    # print( sampled_idxs[0], sampled_idxs[-1], len(temp_wf_sig))

    self.output_wf.fill(0.)

    try:
        coarse_vals =   self.siggen_interp_fn(sampled_idxs)
        self.output_wf[:num_samples_to_fill] = coarse_vals
    except ValueError:
        print( len(self.output_wf) )
        print( num_samples_to_fill)
        print( sampled_idxs)
        sys.exit()

    return self.output_wf[:outputLength]


# class OrtecPPC(Detector):
# ########################################################################################################
#   def GetBoundaryConditions(self):
#       xtal_radius           = self.siggenInst.geometry.xtal_radius
#       xtal_length           = self.siggenInst.geometry.xtal_length
#       pc_length             = self.siggenInst.geometry.pc_length
#       pc_radius             = self.siggenInst.geometry.pc_radius
#
#       bottom_bullet_radius  = self.siggenInst.geometry.bottom_bullet_radius
#       top_bullet_radius     = self.siggenInst.geometry.top_bullet_radius
#
#       wrap_around_radius    = self.siggenInst.geometry.wrap_around_radius
#       ditch_thickness       = self.siggenInst.geometry.ditch_thickness
#       ditch_depth           = self.siggenInst.geometry.ditch_depth
#
#       taper_length          = self.siggenInst.geometry.taper_length
#
#
#       def boundary_pc(x):
#           #Currently the PC is modeled as not-quite hemispherical: its circularly-rounded in whichever corner has the smaller radius
#           #TODO: matt busch thinks a hemisphere is a better idea
#           if x[0] > pc_radius or x[1] > pc_length: return False
#
#           if pc_radius < pc_length:
#               diff = pc_length - pc_radius
#               if x[0] < pc_radius + DOLFIN_EPS and x[1] < diff + DOLFIN_EPS:
#                   return True
#               elif x[0]*x[0] + (x[1] - diff)**2 < pc_radius**2: return True
#               else: return False
#           else:
#               diff = pc_radius - pc_length
#               if x[0] < diff + DOLFIN_EPS and x[1] < pc_length + DOLFIN_EPS:
#                   return True
#               elif (x[0] - diff)**2 + x[1]*x[1] < pc_length**2: return True
#               else:return False
#
#       def boundary_n(x):
#           #Rectangular boundary
#           if x[0] > xtal_radius - DOLFIN_EPS or x[1] > xtal_length - DOLFIN_EPS:
#               return True
#           #Taper
#           elif x[0] >= x[1] + xtal_radius - taper_length:
#               return True
#           #Top bullet radius
#           elif (x[0] > xtal_radius -top_bullet_radius - DOLFIN_EPS) and (x[1] > xtal_length -top_bullet_radius - DOLFIN_EPS):
#               if ((x[0] - xtal_radius +top_bullet_radius)  **2 + (x[1] - xtal_length +top_bullet_radius)**2) > top_bullet_radius**2 - DOLFIN_EPS:
#                   return True
#               else: return False
#           #Bottom bullet radius
#           elif taper_length == 0 and x[0] > (xtal_radius - bottom_bullet_radius - DOLFIN_EPS) and (x[1] < bottom_bullet_radius + DOLFIN_EPS):
#             if ((x[0] - xtal_radius + bottom_bullet_radius)  **2 + (x[1] - bottom_bullet_radius)**2) > bottom_bullet_radius**2 - DOLFIN_EPS:
#                 return True
#             else: return False
#           else:
#               return False
#       def boundary_ditch(x):
#         #Ditch
#         if wrap_around_radius > 0 and  ditch_thickness > 0 and ditch_depth > 0:
#             if x[0] < wrap_around_radius + DOLFIN_EPS and x[0] > wrap_around_radius - ditch_thickness - DOLFIN_EPS and x[1] < ditch_depth + DOLFIN_EPS:
#                 return False
#             elif x[0] > wrap_around_radius - DOLFIN_EPS and x[1] < DOLFIN_EPS:
#                 return True
#             else: return False
#         else:
#             return False
#
#       return (boundary_pc, boundary_n)

  def solve_field(self, field_type, n_r, n_z, impurity_gradient=None, impurity_avg=None, xtal_HV=None):

    if impurity_gradient is None: impurity_gradient = self.siggenInst.GetImpurityGradient()
    if impurity_avg is None:
        impurity_z0 = self.siggenInst.GetImpurity()
    else:
        impurity_z0 = self.avg_to_z0(impurity_avg, impurity_gradient)
    if xtal_HV is None:
        xtal_HV = self.siggenInst.GetXtalHV()

    xtal_radius           = self.siggenInst.geometry.xtal_radius
    xtal_length           = self.siggenInst.geometry.xtal_length
    pc_length             = self.siggenInst.geometry.pc_length
    pc_radius             = self.siggenInst.geometry.pc_radius
    bulletize_PC          = self.siggenInst.geometry.bulletize_PC

    bottom_bullet_radius  = self.siggenInst.geometry.bottom_bullet_radius
    top_bullet_radius     = self.siggenInst.geometry.top_bullet_radius

    taper_length           = self.siggenInst.geometry.taper_length
    wrap_around_radius    = self.siggenInst.geometry.wrap_around_radius
    ditch_thickness       = self.siggenInst.geometry.ditch_thickness
    ditch_depth           = self.siggenInst.geometry.ditch_depth

    is_BEGe = 0
    if taper_length > 0 and wrap_around_radius > 0:
        raise NotImplementedError("Field solver only implemented for wrap_around_radius>0 or taper_length>0: can't do both")
    elif wrap_around_radius >0:
        is_BEGe = 1

    class n_contact(SubDomain):
        def inside(self, x, on_boundary):
            if near(x[0], xtal_radius) or near(x[1], xtal_length):
                return True
            #Taper
            elif taper_length >0 and x[0] >= x[1] + xtal_radius - taper_length:
              return True
            #Wrap-around
            elif wrap_around_radius>0 and near(x[1], 0) and between(x[0], (wrap_around_radius, xtal_radius)):
                return True
            #Top bullet radius
            elif (x[0] > xtal_radius -top_bullet_radius - DOLFIN_EPS) and (x[1] > xtal_length -top_bullet_radius - DOLFIN_EPS):
              if ((x[0] - xtal_radius +top_bullet_radius)  **2 + (x[1] - xtal_length +top_bullet_radius)**2) > top_bullet_radius**2 - DOLFIN_EPS:
                  return True
              else: return False
            #Bottom bullet radius
            elif taper_length == 0 and x[0] > (xtal_radius - bottom_bullet_radius - DOLFIN_EPS) and (x[1] < bottom_bullet_radius + DOLFIN_EPS):
                if ((x[0] - xtal_radius + bottom_bullet_radius)  **2 + (x[1] - bottom_bullet_radius)**2) > bottom_bullet_radius**2 - DOLFIN_EPS:
                    return True
                else: return False
            else:
              return False

    class p_contact(SubDomain):
        def inside(self, x, on_boundary):
            if not (between(x[0], (0, pc_radius)) and between(x[1], (0, pc_length))):
                return False
            elif bulletize_PC:
                return between( (x[0]/pc_radius)**2 + (x[1]/pc_length)**2, (0,1))
            else:
                return True

    class passivated_surface(SubDomain):
        def inside(self, x, on_boundary):
            if is_BEGe:
                return near(x[0], 0) and between(x[0], (pc_radius, wrap_around_radius-ditch_thickness))
            else:
                return near(x[0], 0) and between(x[0], (pc_radius, xtal_radius-taper_length ))

    class Ditch(SubDomain):
        def inside(self, x, on_boundary):
            return (between(x[0], (wrap_around_radius-ditch_thickness, wrap_around_radius)) and between(x[0], (0, ditch_depth)))

    # Initialize sub-domain instances
    n_contact = n_contact()
    p_contact = p_contact()
    pass_surface = passivated_surface()

    if is_BEGe:
        ditch = Ditch()

    # Define mesh
    mesh = RectangleMesh(Point(0.0, 0.0), Point(xtal_radius, xtal_length), n_r, n_z)

    # Initialize mesh function for interior domains
    domains = CellFunction("size_t", mesh)
    domains.set_all(0)

    if is_BEGe:
        ditch.mark(domains, 1)

    # Initialize mesh function for boundary domains
    boundaries = FacetFunction("size_t", mesh)
    boundaries.set_all(0)
    n_contact.mark(boundaries, 1)
    p_contact.mark(boundaries, 2)
    pass_surface.mark(boundaries,3)

    # Define function space and basis functions
    V = FunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)

    # Define new measures associated with the interior domains and
    # exterior boundaries
    dx = Measure('dx', domain=mesh, subdomain_data=domains)
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    r = Expression('x[0]', degree=1)

    epsilon_0 = Constant(8.854187817E-15)
    epsilon_ge = 16*epsilon_0

    if field_type == "efield":
        ###############
        #Solve e-field
        ###############

        # Define charge gradient data (e-field)
        e = 1.602176E-19 #C/charge
        e_per_vol = e * 1E10 * 1E-3 #charge density in units of C/mm^3
        f=Expression("({} + {}*x[1])*{}  ".format(impurity_z0,0.1*impurity_gradient, e_per_vol), degree=1 )

        # Define Dirichlet boundary conditions at electrodes (e-field)
        u_n = Constant(xtal_HV) #n contact at bias voltage
        u_p = Constant(0)       #point contact at 0V
        bcs_efield = [DirichletBC(V, u_n, boundaries, 1), DirichletBC(V, u_p, boundaries, 2)]

        # Define variational form
        if is_BEGe:
            F = (inner(epsilon_ge*grad(u), grad(v))*r*dx(0) + inner(epsilon_0*grad(u), grad(v))*r*dx(1)
             - f*v*r*dx(0))
        else:
            F = (inner(epsilon_ge*grad(u), grad(v))*r*dx(0) - f*v*r*dx(0))

        # Separate left and right hand sides of equation
        a, L = lhs(F), rhs(F)

        # Solve problem
        voltage = Function(V)
        solve(a == L, voltage, bcs_efield)

        V_vec = VectorFunctionSpace(mesh, "CG", 1)
        gradu = project(grad(voltage),V_vec)
        mat_full = np.zeros((n_r,n_z,4), dtype=np.float32)

        for i,r in enumerate(np.linspace(0, xtal_radius, n_r)):
          for j,z in enumerate(np.linspace(0, xtal_length, n_z)):
              voltage_pt = voltage(Point(r,z))
              gradv = gradu(Point(r,z)) #in V/mm
              e_r = -gradv[0] * 10 #to V/cm
              e_z = -gradv[1] * 10 #to V/cm
              e = np.sqrt(e_r**2+e_z**2)
              mat_full[i,j,:] = voltage_pt, e, e_r, e_z

        return mat_full
    elif field_type == "wpot":
        ###############
        #Solve wpotential
        ###############

        # Define Dirichlet boundary conditions at electrodes (w-potential)
        wp_n = Constant(0)
        wp_p = Constant(1)
        bcs_wpot = [DirichletBC(V, wp_n, boundaries, 1), DirichletBC(V, wp_p, boundaries, 2)]

        # Define variational form
        if is_BEGe:
            F = (inner(epsilon_ge*grad(u), grad(v))*r*dx(0) + inner(epsilon_0*grad(u), grad(v))*r*dx(1))
        else:
            F = inner(epsilon_ge*grad(u), grad(v))*r*dx(0)

        # Separate left and right hand sides of equation
        a, L = lhs(F), rhs(F)

        # Solve problem
        wpot = Function(V)
        solve(a == L, wpot, bcs_wpot)

        mat_full = np.zeros((n_r,n_z), dtype=np.float32)
        for i,r in enumerate(np.linspace(0, xtal_radius, n_r)):
          for j,z in enumerate(np.linspace(0, xtal_length, n_z)):
              mat_full[i,j] = wpot(Point(r,z))

        return mat_full

def pdf(x):
    return 1/np.sqrt(2*np.pi) * np.exp(-x**2/2)

def cdf(x):
    return (1 + special.erf(x/np.sqrt(2))) / 2

def skew_norm(x,e=0,w=1,a=0):
    t = (x-e) / w
    return 2 / w * pdf(t) * cdf(a*t)
