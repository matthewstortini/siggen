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

from .Detector import Detector

#Does all the interfacing with siggen for you, stores/loads lookup tables, and does electronics shaping

class PPC(Detector):
  def __init__(self, conf_file, **kwargs):
    super().__init__("PPC", conf_file, **kwargs)

    (self.detector_radius, self.detector_length) = np.floor( [self.detector_radius*10, self.detector_length*10] )/10.

    #TODO: figure out how to write a getter for this
    self.taper_length = 4.5
    self.top_bullet_radius = 1.2

    self.padded_siggen_data = np.zeros(self.total_steps, dtype='f4', order="C")

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
  def MakeSimWaveform(self, r,phi,z,scale, align_point, align_percent, numSamples, smoothing=None, e_smoothing=None):
    hole_wf = self.MakeWaveform(r, phi, z, 1)
    if hole_wf is None:
      return None

    self.padded_siggen_data.fill(0.)
    self.padded_siggen_data[:len(hole_wf[0,:])] = hole_wf[0,:]
    self.padded_siggen_data[len(hole_wf[0,:]):] = hole_wf[0,-1]

    if smoothing is not None:
      ndimage.filters.gaussian_filter1d(self.padded_siggen_data, smoothing, output=self.padded_siggen_data)

    electron_wf = self.MakeWaveform(r, phi, z, -1)

    if e_smoothing is not None:
      ndimage.filters.gaussian_filter1d(electron_wf[0,:], e_smoothing, output=electron_wf[0,:])

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

    #low-pass filter
    temp_wf_sig = signal.lfilter(self.lp_num, self.lp_den, temp_wf_sig)
    temp_wf_sig /= (np.sum(self.lp_num)/np.sum(self.lp_den))

    #hi-pass filter
    temp_wf_sig= signal.lfilter(self.hp_num, self.hp_den, temp_wf_sig, axis=-1)

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
        exit(0)

    return self.output_wf[:outputLength]

########################################################################################################
  def GetBoundaryConditions(self):
      xtal_radius           = self.siggenInst.geometry.xtal_radius
      xtal_length           = self.siggenInst.geometry.xtal_length
      pc_length             = self.siggenInst.geometry.pc_length
      pc_radius             = self.siggenInst.geometry.pc_radius

      bottom_bullet_radius  = self.siggenInst.geometry.bottom_bullet_radius
      top_bullet_radius     = self.siggenInst.geometry.top_bullet_radius

      wrap_around_radius    = self.siggenInst.geometry.wrap_around_radius
      ditch_thickness       = self.siggenInst.geometry.ditch_thickness
      ditch_depth           = self.siggenInst.geometry.ditch_depth

      taper_length          = self.siggenInst.geometry.taper_length


      def boundary_pc(x):
          #Currently the PC is modeled as not-quite hemispherical: its circularly-rounded in whichever corner has the smaller radius
          #TODO: matt busch thinks a hemisphere is a better idea
          if x[0] > pc_radius or x[1] > pc_length: return False

          if pc_radius < pc_length:
              diff = pc_length - pc_radius
              if x[0] < pc_radius + DOLFIN_EPS and x[1] < diff + DOLFIN_EPS:
                  return True
              elif x[0]*x[0] + (x[1] - diff)**2 < pc_radius**2: return True
              else: return False
          else:
              diff = pc_radius - pc_length
              if x[0] < diff + DOLFIN_EPS and x[1] < pc_length + DOLFIN_EPS:
                  return True
              elif (x[0] - diff)**2 + x[1]*x[1] < pc_length**2: return True
              else:return False

      def boundary_n(x):
          #Rectangular boundary
          if x[0] > xtal_radius - DOLFIN_EPS or x[1] > xtal_length - DOLFIN_EPS:
              return True
          #Ditch
          elif wrap_around_radius > 0 and  ditch_thickness > 0 and ditch_depth > 0:
              if x[0] < wrap_around_radius + DOLFIN_EPS and x[0] > wrap_around_radius - ditch_thickness - DOLFIN_EPS and x[1] < ditch_depth + DOLFIN_EPS:
                  return False
              elif x[0] > wrap_around_radius - DOLFIN_EPS and x[1] < DOLFIN_EPS:
                  return True
              else: return False
          #Taper
          elif x[0] >= x[1] + xtal_radius - taper_length:
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

      return (boundary_pc, boundary_n)
