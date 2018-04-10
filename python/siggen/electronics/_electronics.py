# -*- coding: utf-8 -*-

class DigitalFilter(object):

    def __init__(order):
        assert order == 1 or order == 2
        self.order = order

    def zpk_to_ba(self, mag, phi=0):
        if phi == 0:
            return [1, -mag]
        else:
            return [1, -2*mag*np.cos(phi), mag**2]

    @property
    def zero(self):
        return self.__zero

    @property
    def pole(self):
        return self.__pole

    @property
    def num(self):
        return self.__num

    @property
    def den(self):
        return self.__den

    @zero.setter
    def zero(self, mag, phi=0):
        if order ==1: assert phi==0

        if phi == 0:
            self.__zero = mag
        else:
             self.__zero = mag * np.exp(1j*phi )
        self.__num =  self.zpk_to_ba(mag, phi)

    @pole.setter
    def pole(self, mag, phi=0):
        if order ==1: assert phi==0

        if phi == 0:
            self.__pole = mag
        else:
             self.__pole = mag * np.exp(1j*phi )
        self.__den =  self.zpk_to_ba(mag, phi)
