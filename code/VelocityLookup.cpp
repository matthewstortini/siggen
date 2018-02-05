#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "Utils.h"
#include "VelocityLookup.h"

namespace Siggen
{

/* Reference temperature for drift vel. corrections is 77K */
#define REF_TEMP 77.0
/* max, min temperatures for allowed range */
#define MIN_TEMP 40.0
#define MAX_TEMP 120.0

VelocityLookup::VelocityLookup(std::string drift_name):
xtal_temp(REF_TEMP),
drift_name(drift_name),
is_setup(0)
{
}

VelocityLookup::VelocityLookup(std::string drift_name, float xtal_temp):
xtal_temp(xtal_temp),
drift_name(drift_name)
{
  if (xtal_temp == 0) xtal_temp = REF_TEMP;
}

int VelocityLookup::setup_velo(){
  // static int vlook_sz = 0;
  // static struct velocity_lookup *v_lookup;

  if (xtal_temp < MIN_TEMP) xtal_temp = MIN_TEMP;
  if (xtal_temp > MAX_TEMP) xtal_temp = MAX_TEMP;

  char  line[MAX_LINE], *c;
  FILE  *fp;
  int   i;
  struct velocity_lookup *tmp, v, v0;
  float sumb_e, sumc_e, sumb_h, sumc_h;

  double be=1.3e7, bh=1.2e7, thetae=200.0, thetah=200.0;  // parameters for temperature correction
  double pwre=-1.680, pwrh=-2.398, mue=5.66e7, muh=1.63e9; //     adopted for Ge   DCR Feb 2015
  double mu_0_1, mu_0_2, v_s_1, v_s_2, E_c_1, E_c_2, e, f;

  if (vlook_sz == 0) {
    vlook_sz = 10;
    if ((v_lookup = (struct velocity_lookup *)
	 malloc(vlook_sz*sizeof(*v_lookup))) == NULL) {
      error("malloc failed in setup_velo\n");
      return -1;
    }
  }
  if ((fp = fopen(drift_name.c_str(), "r")) == NULL){
    error("failed to open velocity lookup table file: '%s'\n", drift_name.c_str());
    return -1;
  }
  line[0] = '#';
  c = line;
  while ((line[0] == '#' || line[0] == '\0') && c != NULL) c = fgets(line, MAX_LINE, fp);
  if (c == NULL) {
    error("Failed to read velocity lookup table from file: %s\n", drift_name.c_str());
    fclose(fp);
    return -1;
  }
  TELL_CHATTY("Drift velocity table:\n"
	      "  e          e100    e110    e111    h100    h110    h111\n");
  for (v_lookup_len = 0; ;v_lookup_len++){
    if (v_lookup_len == vlook_sz - 1){
      vlook_sz += 10;
      if ((tmp = (struct velocity_lookup *)
	   realloc(v_lookup, vlook_sz*sizeof(*v_lookup))) == NULL){
	error("realloc failed in setup_velo\n");
	fclose(fp);
	return -1;
      }
      v_lookup = tmp;
    }
    if (sscanf(line, "%f %f %f %f %f %f %f",
	       &v_lookup[v_lookup_len].e,
	       &v_lookup[v_lookup_len].e100,
	       &v_lookup[v_lookup_len].e110,
	       &v_lookup[v_lookup_len].e111,
	       &v_lookup[v_lookup_len].h100,
	       &v_lookup[v_lookup_len].h110,
	       &v_lookup[v_lookup_len].h111) != 7){
      break; //assume EOF
    }
    //v_lookup[v_lookup_len].e *= 100; /*V/m*/
    tmp = &v_lookup[v_lookup_len];
    TELL_CHATTY("%10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
		tmp->e, tmp->e100, tmp->e110, tmp->e111, tmp->h100, tmp->h110,tmp->h111);
    line[0] = '#';
    while ((line[0] == '#' || line[0] == '\0' ||
	    line[0] == '\n' || line[0] == '\r') && c != NULL) c = fgets(line, MAX_LINE, fp);
    if (c == NULL) break;
    if (line[0] == 'e' || line[0] == 'h') break; /* no more velocities data;
						    now reading temp correction data */
  }

  /* check for and decode temperature correction parameters */
  while (line[0] == 'e' || line[0] == 'h') {
    if (line[0] == 'e' &&
	sscanf(line+2, "%lf %lf %lf %lf",
	       &mue, &pwre, &be, &thetae) != 4) break;//asume EOF
    if (line[0] == 'h' &&
	sscanf(line+2, "%lf %lf %lf %lf",
	       &muh, &pwrh, &bh, &thetah) != 4) break;//asume EOF
    if (line[0] == 'e')
      TELL_CHATTY("electrons: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
		  mue, pwre, be, thetae);
    if (line[0] == 'h')
      TELL_CHATTY("    holes: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
		  muh, pwrh, bh, thetah);

    line[0] = '#';
    while ((line[0] == '#' || line[0] == '\0') && c != NULL) c = fgets(line, MAX_LINE, fp);
    if (c == NULL) break;
  }

  if (v_lookup_len == 0){
    error("Failed to read velocity lookup table from file: %s\n", drift_name.c_str());
    return -1;
  }
  v_lookup_len++;
  if (vlook_sz != v_lookup_len){
    if ((tmp = (struct velocity_lookup *)
	 realloc(v_lookup, v_lookup_len*sizeof(*v_lookup))) == NULL){
      error("realloc failed in setup_velo. This should not happen\n");
      fclose(fp);
      return -1;
    }
    v_lookup = tmp;
    vlook_sz = v_lookup_len;
  }
  TELL_NORMAL("Drift velocity table has %d rows of data\n", v_lookup_len);
  fclose(fp);

  /*
    apply temperature dependence to mobilities;
    see drift_velocities.doc and tempdep.c
    The drift velocity reduces at higher temperature due to the increasing of
    scattering with the lattice vibration. We used a model by M. Ali Omar and
    L. Reggiani (Solid-State Electronics Vol. 30, No. 12 (1987) 1351) to
    calculate the temperature dependence.
  */
  /* electrons */
  TELL_NORMAL("Adjusting mobilities for temperature, from %.1f to %.1f\n", REF_TEMP, xtal_temp);
  TELL_CHATTY("Index  field  vel_factor\n");
  mu_0_1 = mue * pow(REF_TEMP, pwre);
  v_s_1 = be * sqrt(tanh(0.5 * thetae / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = mue * pow(xtal_temp, pwre);
  v_s_2 = be * sqrt(tanh(0.5 * thetae / xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++){
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].e100 *= f;
    v_lookup[i].e110 *= f;
    v_lookup[i].e111 *= f;
    TELL_CHATTY("%2d %5.0f %f\n", i, e, f);
  }

  /* holes */
  mu_0_1 = muh * pow(REF_TEMP, pwrh);
  v_s_1 = bh * sqrt(tanh(0.5 * thetah / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = muh * pow(xtal_temp, pwrh);
  v_s_2 = bh * sqrt(tanh(0.5 * thetah / xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++){
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].h100 *= f;
    v_lookup[i].h110 *= f;
    v_lookup[i].h111 *= f;
    TELL_CHATTY("%2d %5.0f %f\n", i, e, f);
  }
  /* end of temperature correction */

  for (i = 0; i < vlook_sz; i++){
    v = v_lookup[i];
    v_lookup[i].ea =  0.5 * v.e100 -  4 * v.e110 +  4.5 * v.e111;
    v_lookup[i].eb = -2.5 * v.e100 + 16 * v.e110 - 13.5 * v.e111;
    v_lookup[i].ec =  3.0 * v.e100 - 12 * v.e110 +  9.0 * v.e111;
    v_lookup[i].ha =  0.5 * v.h100 -  4 * v.h110 +  4.5 * v.h111;
    v_lookup[i].hb = -2.5 * v.h100 + 16 * v.h110 - 13.5 * v.h111;
    v_lookup[i].hc =  3.0 * v.h100 - 12 * v.h110 +  9.0 * v.h111;
  }
  v_lookup[0].ebp = v_lookup[0].ecp = v_lookup[0].hbp = v_lookup[0].hcp = 0.0;
  sumb_e = sumc_e = sumb_h = sumc_h = 0.0;
  for (i = 1; i < vlook_sz; i++){
    v0 = v_lookup[i-1];
    v = v_lookup[i];
    sumb_e += (v.e - v0.e)*(v0.eb+v.eb)/2;
    sumc_e += (v.e - v0.e)*(v0.ec+v.ec)/2;
    sumb_h += (v.e - v0.e)*(v0.hb+v.hb)/2;
    sumc_h += (v.e - v0.e)*(v0.hc+v.hc)/2;
    v_lookup[i].ebp = sumb_e/v.e;
    v_lookup[i].ecp = sumc_e/v.e;
    v_lookup[i].hbp = sumb_h/v.e;
    v_lookup[i].hcp = sumc_h/v.e;
  }

  // setup.v_lookup = v_lookup;
  // setup->v_lookup_len = v_lookup_len;

  is_setup = 1;
  return 0;
}

int VelocityLookup::drift_velocity(point cart_en, float abse, float q, float& v_over_E, float& dv_dE, vector *velo)
{
  int   i, sign;
  float absv, f, a, b, c;
  float bp, cp, en4, en6;
  struct velocity_lookup *v_lookup1, *v_lookup2;

  /* find location in table to interpolate from */
  for (i = 0; i < v_lookup_len - 2 && abse > v_lookup[i+1].e; i++);
  v_lookup1 = v_lookup + i;
  v_lookup2 = v_lookup + i+1;
  f = (abse - v_lookup1->e)/(v_lookup2->e - v_lookup1->e);
  if (q > 0){
    a = (v_lookup2->ha - v_lookup1->ha)*f+v_lookup1->ha;
    b = (v_lookup2->hb- v_lookup1->hb)*f+v_lookup1->hb;
    c = (v_lookup2->hc - v_lookup1->hc)*f+v_lookup1->hc;
    bp = (v_lookup2->hbp- v_lookup1->hbp)*f+v_lookup1->hbp;
    cp = (v_lookup2->hcp - v_lookup1->hcp)*f+v_lookup1->hcp;
    dv_dE = (v_lookup2->h100 - v_lookup1->h100)/(v_lookup2->e - v_lookup1->e);
  }else{
    a = (v_lookup2->ea - v_lookup1->ea)*f+v_lookup1->ea;
    b = (v_lookup2->eb- v_lookup1->eb)*f+v_lookup1->eb;
    c = (v_lookup2->ec - v_lookup1->ec)*f+v_lookup1->ec;
    bp = (v_lookup2->ebp- v_lookup1->ebp)*f+v_lookup1->ebp;
    cp = (v_lookup2->ecp - v_lookup1->ecp)*f+v_lookup1->ecp;
    dv_dE = (v_lookup2->e100 - v_lookup1->e100)/(v_lookup2->e - v_lookup1->e);
  }
  /* velocity can vary from the direction of the el. field
     due to effect of crystal axes */
#define POW4(x) ((x)*(x)*(x)*(x))
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
  en4 = POW4(cart_en.x) + POW4(cart_en.y) + POW4(cart_en.z);
  en6 = POW6(cart_en.x) + POW6(cart_en.y) + POW6(cart_en.z);
  absv = a + b*en4 + c*en6;
  sign = (q < 0 ? -1 : 1);
  v_over_E = absv / abse;
  velo->x = sign*cart_en.x*(absv+bp*4*(cart_en.x*cart_en.x - en4)
			    + cp*6*(POW4(cart_en.x) - en6));
  velo->y = sign*cart_en.y*(absv+bp*4*(cart_en.y*cart_en.y - en4)
			    + cp*6*(POW4(cart_en.y) - en6));
  velo->z = sign*cart_en.z*(absv+bp*4*(cart_en.z*cart_en.z - en4)
			    + cp*6*(POW4(cart_en.z) - en6));
#undef POW4
#undef POW6

  return 0;

}

int VelocityLookup::set_xtal_temp(float temp){
  if (temp < MIN_TEMP || temp > MAX_TEMP){
    error("temperature out of range: %f\n", temp);
    return -1;
  }else{
    xtal_temp = temp;
    error("temperature set to %f\n", temp);
    /* re-read velocities and correct them to the new temperature value */
    return setup_velo();
  }
}


} // namespace Siggen
