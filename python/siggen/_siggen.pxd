# distutils: language = c++

from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream& write(const char*, int) except +
    cdef cppclass istream:
        istream& read(char*, int) except +

cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary
    # you can define other constants as needed

cdef extern from "<fstream>" namespace "std":
    cdef cppclass ofstream(ostream):
        # constructors
        ofstream(const char*) except +
        ofstream(const char*, open_mode) except+
    cdef cppclass ifstream(istream):
        # constructors
        ifstream(const char*) except +
        ifstream(const char*, open_mode) except+
        bool eof() const

cdef extern from "Siggen.h" namespace "Siggen":

  cdef struct point:
    float x
    float y
    float z
  cdef struct cyl_pt:
    float r
    float phi
    float z

  cdef cppclass EFieldPoint:
    EFieldPoint()
    float r()
    float phi()
    float z()
    float get_voltage()

    cyl_pt get_field()
    void set_field(cyl_pt new_field)
    void set_voltage(float volt)
    void serialize(ofstream* stream)
    void deserialize(ifstream* stream)

  cdef cppclass Setup:
    Setup (char* conf_file)

    float xtal_radius
    float xtal_length
    int nsegments

  cdef cppclass Detector[T]:
    Detector(T geometry, Setup setup)

    int outside_detector(point pt)
    int field_setup()
    int efield(cyl_pt pt, cyl_pt e)
    int wpotential(point pt, vector[float] wp)
    int get_nsegments()
    float get_impurity()
    float get_impurity_gradient()
    float get_xtal_HV()
    float get_dead_layer()
    string get_field_name()
    string get_wpot_name()

    void set_trapping(double trap_c)
    void set_impurity_z0(float imp, float grad)
    void set_impurity_avg(float imp, float grad)
    void set_use_velo_model(bool usemodel)
    void set_temp(float temp)
    void set_holes(float mu0_100, float beta_100, float E_0_100,
                   float mu0_111, float beta_111, float E_0_111)
    void set_electrons(float mu0_100, float beta_100, float E_0_100, float mu_n_100,
                       float mu0_111, float beta_111, float E_0_111, float mu_n_111)

  cdef cppclass SignalGenerator[T]:
    SignalGenerator(Detector[T]* detector, Setup setup)

    int get_signal(point pt, float* signal_out)
    int make_signal(point pt, float* signal, double q)

    int get_output_length()
    int get_calc_length()
    float get_calc_timestep()

    int get_last_drifttime(float q)
    vector[point] get_driftpath(float q)
    vector[float] get_dwpot()

    void set_calc_timestep(float dt)
    void set_calc_length(int nt)


cdef extern from "ICPC.h":
  cdef cppclass ICPC:
    ICPC(Setup setup);
    float get_xtal_radius()
    float get_xtal_length()

cdef extern from "GEM.h":
  cdef cppclass GEM:
    GEM(Setup setup);
    float get_xtal_radius()
    float get_xtal_length()
