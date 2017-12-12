#ifndef Siggen_PPC
#define Siggen_PPC

#include "Siggen.h"

using namespace Siggen;

class PPC
{
	private:
    std::string field_name, wp_name;

    // int grid_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2][2]);


    float xtal_length,xtal_radius;
    float top_bullet_radius, bottom_bullet_radius;
    float pc_length, pc_radius;
    bool bulletize_PC;
    float wrap_around_radius;
    float ditch_depth,ditch_thickness;
    float taper_length;

    float rmin, rmax, rstep;
    float zmin, zmax, zstep;
    int rlen,zlen;
    float gradmin, gradmax, gradstep, impmin, impmax, impstep;
    int gradnum, impnum;

    Field<EFieldPoint> efld;
    Field<float> wpot;

    int nsegs;

	public:
		// Constructor
		PPC(Setup& setup);
    void parse_setup(std::map<std::string,std::string>& geometry_params);

    int wpotential(point pt, std::vector<float>& wp);
    int efield(cyl_pt pt, float imp_avg, float imp_grad, cyl_pt& e);
    int outside_detector(point pt);

    int setup_efield();
    int setup_wp();

    inline int get_rlen(){return rlen;}
    inline int get_zen(){return zlen;}
    inline int get_nsegments(){return nsegs;}
    inline int efield_exists(cyl_pt pt){return efld.field_exists(pt, *this);}

    inline float get_xtal_radius(){return xtal_radius;}
    inline float get_xtal_length(){return xtal_length;}
    inline float get_pc_length(){return pc_length;}
    inline float get_pc_radius(){return pc_radius;}
    inline float get_wrap_around_radius(){return wrap_around_radius;}
    inline float get_top_bullet_radius(){return top_bullet_radius;}
    inline float get_bottom_bullet_radius(){return bottom_bullet_radius;}
    inline float get_ditch_depth(){return ditch_depth;}
    inline float get_ditch_thickness(){return ditch_thickness;}
    inline float get_taper_length(){return taper_length;}

    inline void set_xtal_radius(float val){xtal_radius = val;}
    inline void set_xtal_length(float val){xtal_length = val;}
    inline void set_pc_length(float val){pc_length = val;}
    inline void set_pc_radius(float val){pc_radius = val;}
    inline void set_wrap_around_radius(float val){wrap_around_radius = val;}
    inline void set_top_bullet_radius(float val){top_bullet_radius = val;}
    inline void set_bottom_bullet_radius(float val){bottom_bullet_radius = val;}
    inline void set_ditch_depth(float val){ditch_depth = val;}
    inline void set_ditch_thickness(float val){ditch_thickness = val;}
    inline void set_taper_length(float val){taper_length = val;}

};



#endif
