#ifndef Siggen_SignalGenerator
#define Siggen_SignalGenerator

// #include "point.h"
#include <cstdlib>
#include <vector>
#include "Detector.h"
#include "Utils.h"

namespace Siggen
{

template<class GeometryType>
class SignalGenerator
{
	private:

    // DetectorType detector;
		Detector<GeometryType>* detector;
    void parse_setup(std::map<std::string, std::string>& param_map);
    void initialize_arrays();

    //for make_signal
    std::vector<float> wpot, wpot_old, dwpot;
    std::vector<float> dwpot_hist;
		std::vector<float> signal_arr, sum, tmp;
    std::vector<point> dpath_e, dpath_h;
    // std::vector<Point> dpath_e, dpath_h;

    int nsegments; /*number of segments, including central contact*/
    // inline int get_segment_start_idx(int segment){return segment*tsteps;};

    int segment_max_wp(std::vector<float>& wp, float thresh);

    int last_hole_drift_time=0;
    int last_electron_drift_time=0;

    //Charge trapping
    double charge_trapping_per_step = 1.0;

    //Preamp
    float preamp_tau = 0.;

    //Individual signal params
    float charge_cloud_size = 0;
    float energy = 0;

    //TODO: don't use this
    float impurity_z0 = 0.;

    //Diffusion stuff
    int use_diffusion=0;
    float initial_vel, final_vel;  // initial and final drift velocities for charges collected to PC
    float dv_dE;     // derivative of drift velocity with field ((mm/ns) / (V/cm))
    float v_over_E;  // ratio of drift velocity to field ((mm/ns) / (V/cm))
    double final_charge_size;     // in mm

    //Signal calculation params
    int time_steps_calc=0;
    int ntsteps_out=0;
    float step_time_calc=0;
    float step_time_out=0;
    int max_iterations=0;

	public:
		// SignalGenerator () {};

		// Constructor: Pass in filename of configuration file
		SignalGenerator(Detector<GeometryType>* detector, Setup& setup_in);

    int get_signal(point pt, float* signal_out);
    int make_signal(point pt, float* signal, double q);
    int rc_integrate(std::vector<float>& s_in, std::vector<float>& s_out, float tau, int time_steps);

    inline int get_output_length(){return ntsteps_out;}
    inline int get_calc_length(){return time_steps_calc;}
    inline int get_nsegments(){return nsegments;}
    inline float get_calc_timestep(){return step_time_calc;}

    inline std::vector<float>& get_dwpot(){return dwpot_hist; }
    inline std::vector<point>& get_driftpath(float q){return (q > 0)?(dpath_h):(dpath_e); }
    inline int get_last_drifttime(float q){return (q > 0)?(last_hole_drift_time):(last_electron_drift_time); }

    inline void set_calc_timestep(float dt){step_time_calc = dt;}
    void set_calc_length(int nt);

};

} // namespace Siggen

#include "SignalGenerator.impl"
#endif
