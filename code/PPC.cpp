#include "PPC.h"
#include "Siggen.h"

#include <iomanip>

using namespace std;
using namespace Siggen;

PPC::PPC(Setup& setup):
field_name(setup.field_name), wp_name(setup.wp_name),
rmin(0),rmax(setup.xtal_radius),rstep(setup.xtal_grid),
zmin(0),zmax(setup.xtal_length),zstep(setup.xtal_grid),
impmin(setup.imp_min),impmax(setup.imp_max),impnum(setup.imp_num),
gradmin(setup.grad_min),gradmax(setup.grad_max),gradnum(setup.grad_num),
nsegs(1)
{
  parse_setup(setup.geometry_map);

  if (gradnum > 1) gradstep = (gradmax - gradmin) / (gradnum-1);
  if (impnum > 1) impstep = (impmax - impmin) / (impnum-1);

  // std::cout <<"impmax: " << impmax <<"\n";
  // std::cout <<"impmin: " << impmin <<"\n";
  // std::cout <<"impstep: " << impstep <<"\n";
  // std::cout <<"impnum: " << impnum <<"\n";
  //
  // std::cout <<"gradmax: " << gradmax <<"\n";
  // std::cout <<"gradmin: " << gradmin <<"\n";
  // std::cout <<"gradstep: " << gradstep <<"\n";
  // std::cout <<"gradnum: " << gradnum <<"\n";

}

void PPC::parse_setup(std::map<std::string,std::string>& geometry_params){
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
    }else if (key == "taper_length"){
      valstream >> taper_length;
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
int PPC::setup_efield(){
  rlen = lrintf((rmax - rmin)/rstep) + 1;
  zlen = lrintf((zmax - zmin)/zstep) + 1;
  TELL_CHATTY("rlen, zlen: %d, %d\n", rlen, zlen);

  efld.set_grid(rlen, zlen, 1, impnum, gradnum,
                rstep,zstep, 1, impstep, gradstep,
                rmin, zmin,0, impmin, gradmin);
  efld.read_data(field_name, *this);
  return 0;
}

/*setup_wp
  read weighting potential values from files. returns 0 on success*/
int PPC::setup_wp(){
  rlen = lrintf((rmax - rmin)/rstep) + 1;
  zlen = lrintf((zmax - zmin)/zstep) + 1;
  TELL_CHATTY("rlen, zlen: %d, %d\n", rlen, zlen);

  wpot.set_grid(rlen, zlen, rstep,zstep,rmin, zmin);
  wpot.read_data(wp_name, *this);
  return 0;
}

//0 for success, 1 for fail
int PPC::wpotential(point pt, std::vector<float>& wp){
  wp[0] = 0;
  int flag=0;
  cyl_pt cyl;
  cyl = cart_to_cyl(pt);
  flag = wpot.get_point_interp(cyl, wp[0], *this);
  if (flag <0) return 1;

  return 0;
}

/* Find (interpolated or extrapolated) electric field for this point */
//0 for success, 1 for fail
int PPC::efield(cyl_pt pt, float imp_avg, float imp_grad, cyl_pt& e){
  e.r=0, e.z=0;
  int flag=0;
  EFieldPoint efld_pt;
  flag = efld.get_point_interp(pt, imp_avg, imp_grad, efld_pt, *this);
  if (flag <0) return 1;

  e.r = efld_pt.r();
  e.z = efld_pt.z();
  e.phi = pt.phi;
  return 0;
}

int PPC::outside_detector(point pt){
  float r, z, br, a;

  z = pt.z;
  if (z >= zmax || z < 0) return 1;

  r = sqrt(SQ(pt.x)+SQ(pt.y));
  if (r > rmax) return 1;
  br = top_bullet_radius;
  if (z > zmax - br &&
      r > (rmax - br) + sqrt(SQ(br)- SQ(z-(zmax - br)))) return 1;
  if (pc_radius > 0 &&
      z <= pc_length && r <= pc_radius) {
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
  if (taper_length > 0 && z < taper_length &&
      r > zmax - taper_length + z) return 1;
  if (ditch_depth > 0 && z < ditch_depth  &&
      ditch_thickness > 0 && wrap_around_radius > 0 &&
      r < wrap_around_radius &&
      r > wrap_around_radius - ditch_thickness) return 1;

  return 0;
}
