#ifndef Siggen_GEM
#define Siggen_GEM

#include "Siggen.h"

using namespace Siggen;

class GEM
{
	private:
    std::string field_name, wp_name;

    // int grid_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2][2]);

    float rmin, rmax, rstep;
    float zmin, zmax, zstep;
    int rlen,zlen;

    float xtal_length,xtal_radius;
    float top_bullet_radius, bottom_bullet_radius;
    float pc_length, pc_radius;
    float bulletize_PC;
    float wrap_around_radius;
    float ditch_depth,ditch_thickness;
    float bottom_taper_length;
    float hole_length,hole_radius;
    float hole_bullet_radius;
    float outer_taper_length,inner_taper_length;
    float outer_taper_width,inner_taper_width;
    float taper_angle;

    Field<EFieldPoint> efld;
    Field<float> wpot;

    int nsegs;

	public:
		// Constructor
		GEM(Setup& setup);
    void parse_setup(std::map<std::string,std::string>& geometry_params);

    int wpotential(point pt, std::vector<float>& wp);
    int efield(cyl_pt pt, float imp_z0, float imp_grad,cyl_pt& e);
    int outside_detector(point pt);

    int setup_efield();
    int setup_wp();

    inline int get_rlen(){return rlen;}
    inline int get_zen(){return zlen;}
    inline int get_nsegments(){return nsegs;}
    inline int efield_exists(cyl_pt pt){return efld.field_exists(pt, *this);}

};



#endif
