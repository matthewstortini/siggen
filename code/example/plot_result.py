#!/usr/local/bin/python
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal

# plt.ion()


def main():
    filename = "output.txt"

    data = np.loadtxt(filename, delimiter=",")

    plt.figure()

    if len(data.shape) > 1:
        # for i in range(1,9):
        for i in range(data.shape[1]):
            plt.plot(data[:,i])
    else:
        plt.plot(data[:])
        # plt.plot(sig)
    plt.show()

main()
