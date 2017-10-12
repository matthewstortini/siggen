#!/usr/local/bin/python

import numpy as np
from scipy import  signal, interpolate, ndimage

from ._siggen import PySiggen_GEM, PySiggen_PPC, PySiggen_ICPC

#Does all the interfacing with siggen for you, stores/loads lookup tables, and does electronics shaping

class Detector:
  def __init__(self,  detector_geometry, conf_file, t0_padding=0, verbose=False):

    if detector_geometry == "GEM":
        self.siggenInst = PySiggen_GEM(conf_file)#PyICPC(conf_file)
    elif detector_geometry == "ICPC":
        self.siggenInst = PySiggen_ICPC(conf_file)
    elif detector_geometry == "PPC":
        self.siggenInst = PySiggen_PPC(conf_file)
    else:
        print("Detector geometry %s is not implemented!" % detector_geometry)
        exit(0)

    self.siggenInst.InitializeFields()

    self.nsegments=self.siggenInst.GetNumSegments()
    self.detector_radius = self.siggenInst.GetMaxRadius()
    self.detector_length = self.siggenInst.GetMaxZ()

    self.num_steps_out = self.siggenInst.GetNumSteps()
    self.signal_array = np.zeros( (self.nsegments, t0_padding+self.num_steps_out), dtype='f4', order="C")
    self.signal_array_flat = np.zeros(self.nsegments*self.num_steps_out , dtype='f4', order="C")

    self.num_steps_calc = self.siggenInst.GetNumStepsCalc()
    self.calc_array = np.zeros( (self.nsegments, t0_padding+self.num_steps_out), dtype='f4', order="C")
    self.calc_array_flat = np.zeros(self.nsegments*self.num_steps_out , dtype='f4', order="C")

    self.t0_padding = t0_padding
    self.verbose = verbose

    self.total_steps = self.num_steps_out + t0_padding
    self.time_step_size = 1
    data_to_siggen_size_ratio = np.around(10. / self.time_step_size,3)
    self.data_to_siggen_size_ratio = np.int(data_to_siggen_size_ratio)

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
    self.calc_array[:,self.t0_padding:] = self.calc_array_flat.reshape((self.nsegments,self.num_steps_out))
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

    self.signal_array[:,self.t0_padding:] = self.signal_array_flat.reshape((self.nsegments,self.num_steps_out))
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
