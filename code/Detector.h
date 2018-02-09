#ifndef Siggen_Detector
#define Siggen_Detector

// #include "point.h"
#include <cstdlib>
#include <vector>
#include "Setup.h"

#include "VelocityModel.h"
#include "VelocityLookup.h"

namespace Siggen
{
template<class GeometryType>
class Detector
{
	private:

    void parse_setup(std::map<std::string, std::string>& detector_params);

    GeometryType geometry;

    //Sort of weird to have one of each, but this is the easiest way to be able to switch at runtime w/o using polymorphism
    VelocityLookup velocity_lookup;
    VelocityModel velocity_model;

    float xtal_radius=0;
    float xtal_length=0;

    //velocity stuff
    bool use_velo_model = true;
    std::string drift_name;
    float xtal_temp = 0.;

    //charge trapping -- exponential in time model (q = q_0 * e^(-t/tau))
    double trapping_constant = -1; // -1: no trapping

    float Li_thickness = 0; //Literally nothing is using this right now, I'll stick it here for now.

    //field-gen related...
    float xtal_grid=0;            // grid size in mm for field files (either 0.5 or 0.1 mm)
    float impurity_z0=0;          // net impurity concentration at Z=0, in 1e10 e/cm3
    float impurity_avg=0;          // net impurity concentration at Z=0.5*xtal_length, in 1e10 e/cm3
    float impurity_gradient=0;    // net impurity gradient, in 1e10 e/cm4
    float impurity_quadratic=0;   // net impurity difference from linear, at z=L/2, in 1e10 e/cm3
    float impurity_surface=0;     // surface impurity of passivation layer, in 1e10 e/cm2
    float impurity_radial_add=0;  // additive radial impurity at outside radius, in 1e10 e/cm3
    float impurity_radial_mult=1.0f; // multiplicative radial impurity at outside radius (neutral=1.0)
    float impurity_rpower=0;      // power for radial impurity increase with radius
    float xtal_HV=0;              // detector bias for fieldgen, in Volts
    int   max_iterations=0;       // maximum number of iterations to use in mjd_fieldgen
    int   write_field=0;          // set to 1 to write V and E to output file, 0 otherwise
    int   write_WP=0;             // set to 1 to calculate WP and write it to output file, 0 otherwise
    int   bulletize_PC=0;         // set to 1 for inside of point contact hemispherical, 0 for cylindrical

    float rmin, rmax, rstep;
    float zmin, zmax, zstep;
    float phimin, phimax, phistep;

    // file names
    std::string field_name;       // potential/efield file name
    std::string  wp_name;          // weighting potential file name

	public:
		// Detector () {};
    Detector (GeometryType& geometry, Setup& setup);

    int field_setup();

    //pass thru to the velocity calculator
    int drift_velocity(point pt, float q, float& v_over_E, float& dv_dE, vector *velo);
    inline void set_temp(float temp){xtal_temp=temp; velocity_lookup.set_xtal_temp(temp);}

    //velocity model
    void set_use_velo_model(bool usemodel);
    inline void set_holes(float mu0_100, float beta_100, float E_0_100,
                          float mu0_111, float beta_111, float E_0_111)
                          {velocity_model.set_holes(mu0_100,beta_100,E_0_100,mu0_111,beta_111,E_0_111);}
    inline void set_electrons(float mu0_100, float beta_100, float E_0_100, float mu_n_100,
                          float mu0_111, float beta_111, float E_0_111, float mu_n_111)
                          {velocity_model.set_electrons(mu0_100,beta_100,E_0_100,mu_n_100,mu0_111,beta_111,E_0_111, mu_n_111);}

    //Pass thru to the geometry info
    inline int wpotential(point pt, std::vector<float>& wp){return geometry.wpotential( pt,  wp);}
    inline int efield(cyl_pt pt, cyl_pt& ret_pt){return geometry.efield( pt, impurity_avg, impurity_gradient,  ret_pt);}
    inline int outside_detector(point pt){return geometry.outside_detector( pt);  }
    // inline int outside_detector_cyl(cyl_pt pt){return geometry.outside_detector_cyl( pt);  }

    inline float get_impurity(){return impurity_z0;}
    inline float get_nsegments(){return geometry.get_nsegments();}
    inline float get_trapping(){return trapping_constant;  }
    inline void set_trapping(double trap_c){ trapping_constant = trap_c;  }

    void set_impurity_avg(float imp, float impgrad);
    void set_impurity_z0(float imp, float impgrad);

};

} // namespace Siggen

#include "Detector.impl"
#endif
