#ifndef Siggen_Setup
#define Siggen_Setup

#include "point.h"
#include "cyl_point.h"
#include <string>
#include <map>

namespace Siggen
{

class Setup
{
	// All SignalGenerator specialisations can access Setup members
	// template<class DetectorType>
	// friend class SignalGenerator;
  // friend class DetectorType;

	public:
    std::map<std::string, std::string> detector_map;
    std::map<std::string, std::string> siggen_map;
    std::map<std::string, std::string> geometry_map;

    //a few things everyone might like to know
    int verbosity=0;              // 0 = terse, 1 = normal, 2 = chatty/verbose
    float xtal_length=0;          // z length
    float xtal_radius=0;          // radius
    float xtal_grid=0;
    std::string field_name;
    std::string wp_name;

    float grad_min, grad_max, imp_min, imp_max;
    int grad_num =1;
    int imp_num=1;

    // float xtal_length=0;          // z length
    // float xtal_radius=0;          // radius
    // float top_bullet_radius=0;    // bulletization radius at top of crystal
    // float bottom_bullet_radius=0; // bulletization radius at bottom of BEGe crystal
    // float pc_length=0;            // point contact length
    // float pc_radius=0;            // point contact radius
    // float taper_length=0;         // size of 45-degree taper at bottom of ORTEC-type crystal
    // float wrap_around_radius=0;   // wrap-around radius for BEGes. Set to zero for ORTEC
    // float ditch_depth=0;          // depth of ditch next to wrap-around for BEGes. Set to zero for ORTEC
    // float ditch_thickness=0;      // width of ditch next to wrap-around for BEGes. Set to zero for ORTEC
    // // float Li_thickness=0;         // depth of full-charge-collection boundary for Li contact
    // float central_hole_l=0;
    // float central_hole_r=0;
    // float taper_z=0;
    // float taper_dr=0;
    // int nsegments=1;

		Setup() {};
		Setup( std::string filename);
		int read_config(std::string config_file_name);
    int set_value(std::string key, std::string value);

};

} // namespace Siggen


#endif
