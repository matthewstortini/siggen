# distutils: language = c++

from libcpp.vector cimport vector

cdef extern from "Siggen.h" namespace "Siggen":

  cdef struct point:
    float x
    float y
    float z
  cdef struct cyl_pt:
    float r
    float phi
    float z

  cdef cppclass Setup:
    Setup (char* conf_file)

    float xtal_radius
    float xtal_length
    int nsegments

  cdef cppclass Detector[T]:
    Detector(Setup setup)

    int outside_detector(point pt)
    int field_setup()
    int efield(cyl_pt pt, cyl_pt e)
    int wpotential(point pt, vector[float] wp)
    int get_nsegments()

    void set_holes(float mu0_100, float beta_100, float E_0_100,
                   float mu0_111, float beta_111, float E_0_111)
    void set_electrons(float mu0_100, float beta_100, float E_0_100, float mu_n_100,
                       float mu0_111, float beta_111, float E_0_111, float mu_n_111)

  cdef cppclass SignalGenerator[T]:
    SignalGenerator(Detector[T]* detector, Setup setup)

    int get_signal(point pt, float* signal_out)
    int make_signal(point pt, float* signal, float q)
    int get_output_length()
    int get_calc_length()
    int get_last_drifttime(float q)
    vector[point] get_driftpath(float q)

cdef extern from "PPC.h":
  cdef cppclass PPC:
    PPC(Setup setup);
cdef extern from "ICPC.h":
  cdef cppclass ICPC:
    ICPC(Setup setup);
cdef extern from "GEM.h":
  cdef cppclass GEM:
    GEM(Setup setup);
