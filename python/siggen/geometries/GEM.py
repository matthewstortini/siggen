#!/usr/local/bin/python
import sys
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

from ..Detector import Detector

#Does all the interfacing with siggen for you, stores/loads lookup tables, and does electronics shaping

class GEM(Detector):
  def __init__(self, conf_file, maxWfOutputLength= 500, **kwargs):
    super(GEM, self).__init__().__init__("PPC", conf_file, **kwargs)

    (self.detector_radius, self.detector_length) = np.floor( [self.detector_radius*10, self.detector_length*10] )/10.

    #TODO: figure out how to write a getter for this
    self.taper_length = 4.5
    self.top_bullet_radius = 1.2

    self.maxWfOutputLength = maxWfOutputLength

    self.padded_siggen_data = np.zeros(self.total_steps, dtype='f4', order="C")


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
        print( start_idx)
        print( num_samples_to_fill)
        print( sampled_idxs)
        sys.exit()

    return self.output_wf[:outputLength]

########################################################################################################
  def GetBoundaryConditions(self):
      xtal_radius = self.detector_radius
      xtal_length = self.detector_length

      pc_length=            1.7  # point contact length
      pc_radius=            2.5  # point contact radius
      hole_length=          52   # length of hole, for inverted-coax style
      hole_radius=          5   # radius of hole, for inverted-coax style
      hole_bullet_radius=   5
      top_bullet_radius=    3   # bulletization radius at top of crystal
      bottom_taper_length = 3
      wrap_around_radius=   25   # wrap-around radius
      ditch_depth=          3   # depth of ditch next to wrap-around
      ditch_thickness=      3   # width of ditch next to wrap-around

          # Define Dirichlet boundary (x = 0 or x = 1)
      def boundary_pc(x):
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
          #Rectangular Boundary
          if x[0] > xtal_radius - DOLFIN_EPS or x[1] > xtal_length - DOLFIN_EPS:
              return True
          #Ditch
          elif x[0] < wrap_around_radius + DOLFIN_EPS and x[0] > wrap_around_radius - ditch_thickness - DOLFIN_EPS and x[1] < ditch_depth + DOLFIN_EPS:
              return True
          #Hole
          elif x[0] < hole_radius + DOLFIN_EPS and x[1] > xtal_length - hole_length - DOLFIN_EPS:
              rdiff = hole_radius - hole_bullet_radius
              zdiff = xtal_length - hole_length + hole_bullet_radius
              if x[1] > xtal_length - hole_length + hole_bullet_radius - DOLFIN_EPS:
                  return True
              elif x[0] > rdiff - DOLFIN_EPS and x[1] < zdiff + DOLFIN_EPS and (x[0] - rdiff)**2 + (x[1] - zdiff)**2 > hole_bullet_radius**2 -DOLFIN_EPS:
                  return False
              else:
                  return True
        #   #TODO: hole radius
          elif (x[0] > xtal_radius -top_bullet_radius - DOLFIN_EPS) and (x[1] > xtal_length -top_bullet_radius - DOLFIN_EPS):
              if ((x[0] - xtal_radius +top_bullet_radius)  **2 + (x[1] - xtal_length +top_bullet_radius)**2) > top_bullet_radius**2 - DOLFIN_EPS:
                  return True
              else: return False
          #Taper
          elif x[0] >= x[1] + xtal_radius - bottom_taper_length:
              return True
          else:
              return False

      return (boundary_pc, boundary_n)
