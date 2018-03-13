#include "ICPC.h"
#include "Siggen.h"

#include <iomanip>

using namespace std;
using namespace Siggen;

ICPC::ICPC(Setup& setup):
field_name(setup.field_name), wp_name(setup.wp_name),
xtal_radius(setup.xtal_radius), xtal_length(setup.xtal_length),
rmin(0),rmax(setup.xtal_radius),rstep(setup.xtal_grid),
zmin(0),zmax(setup.xtal_length),zstep(setup.xtal_grid),
phimin(0),phimax(359.),phistep(1),
wpot(20),
nsegs(20)
{
  parse_setup(setup.geometry_map);

  rmin  = 0;
  rmax  = setup.xtal_radius;
  rstep = setup.xtal_grid;
  zmin  = 0;
  zmax  = setup.xtal_length;
  zstep = setup.xtal_grid;


  phimin  /= 180.0/M_PI;  // degrees to radians
  phimax  /= 180.0/M_PI;
  phistep /= 180.0/M_PI;
}

void ICPC::parse_setup(std::map<std::string,std::string>& geometry_params){
  for (auto const& x : geometry_params)
  {
    std::stringstream valstream(x.second);
    std::string key =  x.first;

    if (key == "xtal_radius"){
      valstream >> rmax;
    }else if (key == "xtal_length"){
      valstream >> zmax;
    }else if (key == "top_bullet_radius"){
      valstream >> bullet_radius;
    }else if (key == "pc_length"){
      valstream >> contact_l;
    }else if (key == "pc_radius"){
      valstream >> contact_r;
    }else if (key == "central_hole_l"){
      valstream >> central_hole_l;
    }else if (key == "central_hole_r"){
      valstream >> central_hole_r;
    }else if (key == "taper_z"){
      valstream >> taper_z;
    }else if (key == "taper_dr"){
      valstream >> taper_dr;
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
int ICPC::setup_efield(){
  rlen = 1.5 + (rmax - rmin)/rstep; // round up and add 1
  zlen = 1.5 + (zmax - zmin)/zstep;
  efld.set_grid(rlen, zlen, 1, rstep,zstep,phistep, rmin, zmin, phimin);
  efld.read_data(field_name, *this);
  return 0;
}

/*setup_wp
  read weighting potential values from files. returns 0 on success*/
int ICPC::setup_wp(){
  int i;

  rlen = 1.5 + (rmax - rmin)/rstep;  // round up and add 1
  zlen = 1.5 + (zmax - zmin)/zstep;
  alen = 2.5 + (phimax - phimin)/phistep;  // round up and add 2

  wpot.resize(nsegs);

  for (i = 0; i < nsegs; i++){
    wpot.at(i).set_grid(rlen, zlen, 1, rstep,zstep,phistep, rmin, zmin, phimin);
  }
  az_wpot.set_grid(rlen, zlen, alen, rstep,zstep,phistep, rmin, zmin, phimin);

  /*now read the tables*/
  for (int segno = 0; segno < nsegs; segno++){
    std::ostringstream filename;
    filename << wp_name <<"_" << setfill('0') << setw(2) << segno << ".dat";

    if (segno == 1) {  // azimuthal segment
      az_wpot.read_data(filename.str(), *this);
      segno += 7;
    } else {  // rotationally symmetric segments
      wpot.at(segno).read_data(filename.str(), *this);
    }
  }
  return 0;
}

//0 for success, 1 for fail
int ICPC::wpotential(point pt, std::vector<float>& wp){
  float w[2][2][2];
  int i, j, k, a, aa;
  static int aa45deg=0;
  cyl_int_pt ipt;
  cyl_pt cyl;
  int res;

  if (aa45deg == 0) {
    aa45deg = (int) (1.5 + (phimax - phimin)/phistep) / 8;
    // printf(">>>>>>> aa45deg = %d <<<<<<<<<<\n", aa45deg);
    // printf("phimax %0.2f, phimin %0.2f, phistep %0.2f\n", phimax, phimin, phistep);
  }

  cyl = cart_to_cyl(pt);
  res = efld.nearest_field_grid_index(cyl, &ipt, *this);
  if (res < 0) return 1;

  az_wpot.position_weights(cyl, ipt, w);
  // printf("res: %d, %d, %d \n", ipt.r, ipt.phi, ipt.z);
  // printf("weights: %e, %e,%e, %e \n", w[0][0][0], w[0][1][0], w[1][0][1],w[1][1][1]);
  // exit(0);

  for (k=0; k<nsegs; k++) {
    wp[k] = 0.0;

    if (k > 0 && k < 9) {  // azimuthal segments, 45 degrees each; use 45-deg symm
      for (i = 0; i < 2; i++){
        for (j = 0; j < 2; j++){
          for (a = 0; a < 2; a++){
            aa = ipt.phi + a - aa45deg * (k-1);
            // printf("a %d, aa %d, iphi %d, phi %f, weight %e\n", a, aa, ipt.phi, cyl.phi, w[i][j][a]);
            if (aa < 0) aa += 8 * aa45deg;
            if (w[i][a][j] == 0.){
              wp[k] += 0.;
            }else{
              wp[k] += w[i][a][j]*az_wpot(ipt.r+i, aa, ipt.z+j);
            }
          }
        }
      }
    } else {  // rotationally symmetric segments
      for (i = 0; i < 2; i++){
        for (j = 0; j < 2; j++){
          // ignore angle, no dependence
            wp[k] += (w[i][0][j] + w[i][1][j])*wpot.at(k)( ipt.r+i, ipt.z+j);
        }
      }
    }
  }
  return 0;
}

/* Find (interpolated or extrapolated) electric field for this point */
//0 for success, 1 for fail
int ICPC::efield(cyl_pt pt,float imp_z0, float imp_grad, cyl_pt& e){
  int flag=0;
  EFieldPoint efld_pt;
  flag = efld.get_point_interp(pt, imp_z0, imp_grad,efld_pt, *this);
  if (flag <0) return 1;

  e.r = efld_pt.r();
  e.z = efld_pt.z();
  e.phi = pt.phi;
  return 0;
}

int ICPC::in_crystal(point pt){
  float r, z, br;
  z = pt.z;
  if (z >= zmax || z < 0){
    // printf("not in detector 5\n");
    return 0;
  }

  r = sqrt(SQ(pt.x)+SQ(pt.y));
  if (r > rmax){
    // printf("not in detector 4\n");
    return 0;
  }
  if (z > taper_z && r > rmax - (taper_dr * (z - taper_z) / (zmax - taper_z))){
    // printf("not in detector 3\n");
    return 0;
  }

  br = bullet_radius;
  if (z > zmax - br){
    if (r > (rmax - br) + sqrt(SQ(br)- SQ(z-(zmax - br))))
    // printf("not in detector 2\n");
      return 0;
  }
  if (contact_r > 0){
    if (z <= contact_l && r <= contact_r)
    // printf("not in detector 1\n");
      return 0;
  }

  if (central_hole_r > 0){
    if (z >= zmax - central_hole_l && r <= central_hole_r){
      // printf("not in detector 1\n");
      return 0;
    }
  }

  return 1;
}
