#include "GEM.h"
#include "Siggen.h"

#include <iomanip>

using namespace std;
using namespace Siggen;

GEM::GEM(Setup& setup):
field_name(setup.field_name), wp_name(setup.wp_name),
rmin(0),rmax(setup.xtal_radius),rstep(setup.xtal_grid),
zmin(0),zmax(setup.xtal_length),zstep(setup.xtal_grid),
nsegs(1)
{
  parse_setup(setup.geometry_map);

  rmin  = 0;
  rmax  = xtal_radius;
  rstep = setup.xtal_grid;
  zmin  = 0;
  zmax  = xtal_length;
  zstep = setup.xtal_grid;

}

void GEM::parse_setup(std::map<std::string,std::string>& geometry_params){
  for (auto const& x : geometry_params)
  {
    std::stringstream valstream(x.second);
    std::string key =  x.first;

    if (key == "xtal_length"){
      valstream >> xtal_length;
    }else if (key == "xtal_radius"){
      valstream >> xtal_radius;
    }else if (key == "top_bullet_radius"){
      valstream >> top_bullet_radius;
    }else if (key == "bottom_bullet_radius"){
      valstream >> bottom_bullet_radius;
    }else if (key == "pc_length"){
      valstream >> pc_length;
    }else if (key == "pc_radius"){
      valstream >> pc_radius;
    }else if (key == "bulletize_PC"){
      valstream >> bulletize_PC;
    }else if (key == "wrap_around_radius"){
      valstream >> wrap_around_radius;
    }else if (key == "ditch_depth"){
      valstream >> ditch_depth;
    }else if (key == "ditch_thickness"){
      valstream >> ditch_thickness;
    }else if (key == "bottom_taper_length"){
      valstream >> bottom_taper_length;
    }else if (key == "hole_length"){
      valstream >> hole_length;
    }else if (key == "hole_radius"){
      valstream >> hole_radius;
    }else if (key == "hole_bullet_radius"){
      valstream >> hole_bullet_radius;
    }else if (key == "outer_taper_width"){
      valstream >> outer_taper_width;
    }else if (key == "inner_taper_width"){
      valstream >> inner_taper_width;
    }else if (key == "outer_taper_length"){
      valstream >> outer_taper_length;
    }else if (key == "inner_taper_length"){
      valstream >> inner_taper_length;
    }else if (key == "taper_angle"){
      valstream >> taper_angle;
    }else{
      std::cout << "WARNING: unrecognized geometry keyword " << key << "\n";
    }
  }
}


/* This may or may not break if we switch to a non-integer grid*/
/*setup_efield
  read electric field data from file, apply sanity checks
  returns 0 for success
*/
int GEM::setup_efield(){
  rlen = lrintf((rmax - rmin)/rstep) + 1;
  zlen = lrintf((zmax - zmin)/zstep) + 1;
  TELL_CHATTY("rlen, zlen: %d, %d\n", rlen, zlen);

  efld.set_grid(rlen, zlen, rstep,zstep, rmin, zmin);
  efld.read_data(field_name, *this);
  return 0;
}

/*setup_wp
  read weighting potential values from files. returns 0 on success*/
int GEM::setup_wp(){
  rlen = lrintf((rmax - rmin)/rstep) + 1;
  zlen = lrintf((zmax - zmin)/zstep) + 1;
  TELL_CHATTY("rlen, zlen: %d, %d\n", rlen, zlen);

  wpot.set_grid(rlen, zlen, rstep,zstep,rmin, zmin);
  wpot.read_data(wp_name, *this);
  return 0;
}

//0 for success, 1 for fail
int GEM::wpotential(point pt,std::vector<float>& wp){
  wp[0] = 0;
  int flag=0;
  cyl_pt cyl;
  cyl = cart_to_cyl(pt);
  flag = wpot.get_point_interp(cyl,  wp[0], *this);
  if (flag <0) return 1;

  return 0;
}

/* Find (interpolated or extrapolated) electric field for this point */
//0 for success, 1 for fail
int GEM::efield(cyl_pt pt, float imp_z0, float imp_grad,cyl_pt& e){
  e.r=0,e.phi=0,e.z=0;
  int flag=0;
  EFieldPoint efld_pt;
  flag = efld.get_point_interp(pt, imp_z0, imp_grad,efld_pt, *this);
  if (flag <0) return 1;

  e.r = efld_pt.r();
  e.z = efld_pt.z();
  e.phi = pt.phi;
  return 0;
}

int GEM::outside_detector(point pt){
  float r, z, br, a, b;

  z = pt.z;
  if (z > zmax || z < 0) return 1;

  r = sqrt(SQ(pt.x)+SQ(pt.y));
  if (r > rmax) return 1;
  br = top_bullet_radius;
  if (z > zmax - br &&
      r > (rmax - br) + sqrt(SQ(br)- SQ(z-(zmax - br)))) return 1;
  if (pc_radius > 0 &&
      z < pc_length && r < pc_radius) {
    if (!bulletize_PC) return 1;
    if (pc_length > pc_radius) {
      a = pc_length - pc_radius;
      if (z < a || SQ(z-a) + SQ(r) < SQ(pc_radius)) return 1;
    } else {
      a = pc_radius - pc_length;
      if (r < a || SQ(z) + SQ(r-a) < SQ(pc_length)) return 1;
    }
    return 0;
  }
  if (bottom_taper_length > 0 && z < bottom_taper_length &&
      r > zmax - bottom_taper_length + z) return 1;
  if (ditch_depth > 0 && z < ditch_depth  &&
      ditch_thickness > 0 && wrap_around_radius > 0 &&
      r < wrap_around_radius &&
      r > wrap_around_radius - ditch_thickness) return 1;

  /* check hole */
  if (r < hole_radius &&
      z > zmax - hole_length) {
    b = zmax - hole_length + hole_bullet_radius;
    if (z > b) return 1;
    a = hole_radius - hole_bullet_radius;
    if (r < a || SQ(b-z) + SQ(r-a) < SQ(hole_bullet_radius)) return 1;
    return 0;
  }

  /* check outer taper of crystal */
  if (outer_taper_length > 0 &&
      z > zmax - outer_taper_length &&
      r > rmax - ((z - zmax + outer_taper_length) *
                         outer_taper_width / outer_taper_length)) return 1;
  /* check inner taper of hole */
  if (inner_taper_length > 0 &&
      z > zmax - inner_taper_length &&
      r < hole_radius + ((z - zmax + inner_taper_length) *
                         inner_taper_width / inner_taper_length)) return 1;

  return 0;

}
