#ifndef Siggen_ICPC
#define Siggen_ICPC

#include "Siggen.h"

using namespace Siggen;

class ICPC
{
	private:
    std::string field_name, wp_name;

    // int grid_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2][2]);

    float rmin, rmax, rstep;
    float zmin, zmax, zstep;
    float phimin, phimax, phistep;

    float bullet_radius; /*bulletization radius @ front of det (z=zmax)*/
    float contact_l, contact_r; /*length, radius of central contact*/
    float central_hole_l, central_hole_r; /*length and radius of central hole*/
    float taper_z;       /* z where taper begins */
    float taper_dr;      /* reduction in r at z=zlen due to tapering */

    Field<EFieldPoint> efld;
    std::vector<Field<float>> wpot;
    Field<float> az_wpot;

    int nsegs,rlen,zlen, alen;

	public:
		// Constructor
    // ICPC(){};
		ICPC(Setup& setup);
    void parse_setup(std::map<std::string,std::string>& geometry_params);

    /* wpotential
       gives (interpolated or extrapolated ) weighting potential
       at point pt. These values are stored in wp.
       returns 0 for success, 1 on failure.
    */
    int wpotential(point pt, std::vector<float>& wp);
    int efield(cyl_pt pt, cyl_pt& e);
    int in_crystal(point pt);

    int setup_efield();
    int setup_wp();


    inline int outside_detector(point pt){return !in_crystal(pt);};

    inline int segment_number(point pt){
      if (!in_crystal(pt)) return -1;
      return 0;
    }

    inline int get_rlen(){return rlen;}
    inline int get_zen(){return zlen;}
    inline int get_nsegments(){return nsegs;}
    inline int efield_exists(cyl_pt pt){return efld.field_exists(pt, *this);}

};



#endif
