#!/usr/local/bin/python

import numpy as np
try:
    from scipy import  signal, interpolate, ndimage
except ImportError:
    pass

from ..Detector import Detector

#Does all the interfacing with siggen for you, stores/loads lookup tables, and does electronics shaping

class ICPC(Detector):
  def __init__(self, conf_file, **kwargs):
    super(ICPC, self).__init__("ICPC", conf_file, **kwargs)

    self.r_core = 5
    self.z_core = 25
    self.z_taper = 20
    self.r_taper_sub = 10.5
    self.pc_idx = 0
    self.core_rc = 0


  def GetRLims(self,z):
      if z < self.z_taper:
          return (0, self.detector_radius)
      else:
          r_min=0
          if z > self.z_core: r_min = self.r_core
          r_max = self.detector_radius - (self.r_taper_sub * (z - self.z_taper) / (self.detector_length - self.z_taper))
          return (r_min, r_max)

  def GetZLims(self,r):
      if r < self.r_core:
          return (0,self.z_core)
      elif r < (self.detector_radius - self.r_taper_sub):
          return (0,self.detector_length)
      else:
          zmax = self.z_taper + (self.detector_length - self.z_taper) * (self.detector_radius - r) / (self.r_taper_sub  )
          return (0,zmax)

  def SetCoreRC(self, RC_in_ns,  digPeriod  = 1E9):
      self.core_rc = RC_in_ns
      RC= 1E-9 * (RC_in_ns)

      self.core_num = [1 - np.exp(-1./digPeriod/RC) ,0]
      self.core_den = [1, -np.exp(-1./digPeriod/RC)]

###########################################################################################################################
  def GetSimWaveform(self, r,phi,z,scale, align_point,  numSamples, smoothing=None):
    sig_wf = self.GetRawSiggenWaveform(r, phi, z)
    if sig_wf is None:
      return None

    sig_wf *= scale

    #only do gaussian filtering of PC signal
    if smoothing is not None:
      ndimage.filters.gaussian_filter1d(sig_wf[self.pc_idx,:], smoothing, output=sig_wf[self.pc_idx,:])

    sim_wf = self.ProcessWaveform(sig_wf, align_point,  numSamples)

    return sim_wf

########################################################################################################
  def GetRawSiggenWaveform(self, r,phi,z, energy=1):
    self.signal_array_flat.fill(0)
    self.signal_array.fill(0)

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    calcFlag = self.siggenInst.GetSignal(x, y, z, self.signal_array_flat[:]);

    if calcFlag == -1:
      print ("Holes out of crystal alert! (%0.3f,%0.3f,%0.3f).  IsInDetector: %d" % (r,phi,z, self.IsInDetector(r,phi,z)))
      return None
    # if not np.any(self.raw_siggen_data):
    #   print( "found zero wf at r={0}, phi={1}, z={2} (calcflag is {3})".format(r, phi, z, calcFlag) )
    #   return None

    self.signal_array[:,self.t0_padding:] = self.signal_array_flat.reshape((self.nsegments,self.num_steps))

    return self.signal_array


########################################################################################################
  def ProcessWaveform(self, siggen_wf,  align_point, outputLength):
    #low-pass filter (non-ppc)
    temp_wf_sig = np.zeros_like(siggen_wf)
    for i in range(self.nsegments):
        temp_wf_sig[i,:] = signal.lfilter(self.lp_num[i], self.lp_den[i], siggen_wf[i,:])
        # print(temp_wf_sig[i,:])
        # print(self.lp_num[i] ,self.lp_den[i], np.sum(self.lp_den[i]))
        temp_wf_sig[i,:] /= (np.sum(self.lp_num[i])/np.sum(self.lp_den[i]))

    #RC filter for core
    if self.core_rc > 0:
        temp_wf_sig[-1,:] = signal.lfilter(self.core_num, self.core_den, siggen_wf[-1,:])
    #hi-pass filter
    temp_wf_sig= signal.lfilter(self.hp_num, self.hp_den, temp_wf_sig, axis=-1)

    #uh, downsample it
    temp_wf = temp_wf_sig[:,::self.data_to_siggen_size_ratio]

    #enforces t0 time-alignment to the align_point
    siggen_offset = self.t0_padding/self.data_to_siggen_size_ratio
    first_idx= 0

    num_wfs = self.nsegments
    num_samples = self.total_steps/self.data_to_siggen_size_ratio

    #find the first index after the align point
    align_point_ceil = np.int( np.ceil(align_point) )
    start_idx = align_point_ceil - first_idx
    #and make sure its above 0!
    if start_idx <0:
        print("bad start idx")
        return None

    siggen_interp_fn = interpolate.interp1d(np.arange(num_samples), temp_wf, kind=self.interpType, copy="False", assume_sorted="True")

    num_samples_to_fill = outputLength - start_idx
    offset = align_point_ceil - align_point
    sampled_idxs = np.arange(num_samples_to_fill) + offset + siggen_offset

    processed_siggen_data = np.zeros((num_wfs, outputLength))
    coarse_vals =   siggen_interp_fn(sampled_idxs)

    try:
        processed_siggen_data[:, start_idx:start_idx+num_samples_to_fill] = coarse_vals
    except ValueError:
        print( len(processed_siggen_data) )
        print( start_idx)
        print( num_samples_to_fill)
        print( sampled_idxs)
        exit(0)

    return processed_siggen_data
