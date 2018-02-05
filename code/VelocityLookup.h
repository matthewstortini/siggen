#ifndef Siggen_VelocityLookup
#define Siggen_VelocityLookup

#include "Setup.h"

namespace Siggen
{

  // from fields.c
  struct velocity_lookup{
    float e;
    float e100;
    float e110;
    float e111;
    float h100;
    float h110;
    float h111;
    float ea; //coefficients for anisotropic drift
    float eb;
    float ec;
    float ebp;
    float ecp;
    float ha;
    float hb;
    float hc;
    float hbp;
    float hcp;
    float hcorr;
    float ecorr;
  };

class VelocityLookup
{
	private:
    velocity_lookup* v_lookup;
    int v_lookup_len=0;
    int vlook_sz=0;
    float xtal_temp;
    std::string drift_name;
    bool is_setup;

	public:
    VelocityLookup(){};
    VelocityLookup (std::string drift_name);
    VelocityLookup (std::string drift_name, float xtal_temp);

    int setup_velo();

    int drift_velocity(point cart_en, float abse, float q, float& v_over_E, float& dv_dE, vector *velo);

    void set_drift_name(std::string drift){drift_name = drift;}

    inline bool get_is_setup(){return is_setup;}

    int set_xtal_temp(float temp);

};

} // namespace Siggen
#endif
