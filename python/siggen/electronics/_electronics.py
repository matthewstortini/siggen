# -*- coding: utf-8 -*-
import numpy as np
try:
    from scipy import  signal
except ImportError:
    pass

class DigitalFilter(object):

    def __init__(self, order):
        assert order == 1 or order == 2
        self.order = order

    def zpk_to_ba(self, mag, phi=0):
        if phi == 0:
            return [1, -mag]
        else:
            return [1, -2*mag*np.cos(phi), mag**2]

    def apply_to_signal(self, sig):
        sig = signal.lfilter(self.__num, self.__den, sig)
        if self.__num_sum != 0:
            sig /= (self.__num_sum/self.__den_sum)
        return sig

    def set_zeros(self, mag, phi=0):
        if self.order ==1:
            assert phi==0
            self.__zeros = [mag]
            self.num =  [1, -mag]
        else:
            if phi == 0:
                self.__zeros = [mag, mag]
            else:
                 self.__zeros = [mag * np.exp(1j*phi ), mag * np.exp(-1j*phi )]
            self.num =  self.zpk_to_ba(mag, phi)

    def set_poles(self, mag, phi=0):
        if self.order ==1:
            assert phi==0
            self.__poles = [mag]
            self.den =  [1, -mag]
        else:
            if phi == 0:
                self.__poles = [mag, mag]
            else:
                 self.__poles = [mag * np.exp(1j*phi ), mag * np.exp(-1j*phi )]
            self.den =  self.zpk_to_ba(mag, phi)

    @property
    def num(self):
        return self.__num

    @property
    def den(self):
        return self.__den

    @den.setter
    def den(self, den):
        self.__den =  den
        self.__den_sum = np.sum(self.__den)

    @num.setter
    def num(self, num):
        self.__num =  num
        self.__num_sum = np.sum(self.__num)


class GretinaOvershootFilter(DigitalFilter):
    def __init__(self, order, overshoot_frac=0):
        super(GretinaOvershootFilter, self).__init__(order)
        self.overshoot_frac = overshoot_frac

    def apply_to_signal(self, sig):
        sig += signal.lfilter(self.__num, self.__den, self.overshoot_frac*sig)
        return sig
