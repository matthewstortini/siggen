#include <math.h>
#include "VelocityModel.h"

#include <iostream>

namespace Siggen
{


AxisModel::AxisModel(float mu_0, float beta, float E_0):
mu_0(mu_0),
beta(beta),
E_0(E_0),
mu_n(0.)
{}

AxisModel::AxisModel(float mu_0, float beta, float E_0, float mu_n):
mu_0(mu_0),
beta(beta),
E_0(E_0),
mu_n(mu_n)
{}

float AxisModel::get_velocity(float E){
  float v;

  v = (mu_0 * E) / pow(1+pow(E/E_0,beta), 1./beta) - mu_n*E;
  return v * 10 * 1E-9;
}

void AxisModel::set_params(float mu_0_in, float beta_in, float E_0_in, float mu_n_in)
{
  mu_0 = mu_0_in;
  beta = beta_in;
  E_0  = E_0_in;
  mu_n = mu_n_in;
}

// void AxisModel::set_params(float mu_0_in, float beta_in, float E_0_in)
// {
//   mu_0 = mu_0_in;
//   beta = beta_in;
//   E_0  = E_0_in;
//   mu_n = 0;
// }



VelocityModel::VelocityModel()
{
  //Initialize gamma matrices for the electron velocity calculation
  int i;

  double m_l = 1.64;
  double m_t = 0.0819;
  double gamma_0[3][3] = { {1./m_t,0,0},{0,1./m_l,0},{0,0,1./m_t} };

  double beta=0, alpha=0;
  double R_j[3][3],R_j_T[3][3], gamma_Rj[3][3];

  double R_x[3][3] = { { 1., 0, 0}, {0, cos(beta), sin(beta)}, {0, -sin(beta), cos(beta)}};
  double R_z[3][3] = { {cos(alpha), sin(alpha), 0}, {-sin(alpha), cos(alpha), 0}, {0,0,1.} };

  for (i = 0; i<4;++i)
  {
    beta = acos(sqrt(2./3));
    alpha = (i) * M_PI / 2. + M_PI/4;

    R_x[1][1] = cos(beta), R_x[1][2] = sin(beta);
    R_x[2][1] = -sin(beta), R_x[2][2] = cos(beta);
    R_z[0][0] = cos(alpha), R_z[0][1] = sin(alpha);
    R_z[1][0] = -sin(alpha), R_z[1][1] = cos(alpha);

    mat_mult(R_x, R_z, R_j);
    transpose(R_j, R_j_T);

    mat_mult(gamma_0, R_j, gamma_Rj);
    mat_mult(R_j_T, gamma_Rj, gamma_j[i]);
  }

  //may as well also calculate Gamma_0 & n_j/n parameters
  double gamma_E[4][3];
  double EgammaE[4];

  double E_100[3] = {1,0,0};
  for (i = 0; i<4;++i){
    dot(gamma_j[i], E_100, gamma_E[i]);
    EgammaE[i] = dot(E_100, gamma_E[i]);

    //only need to look at first dimension of gamma_E, since E_0 is along <100>
    Gamma_0 += 0.25 * gamma_E[i][0] / sqrt(EgammaE[i]);
  }

  double theta_111 = M_PI/2. - acos(sqrt(2./3));
  double phi_111 = M_PI/4.;
  double E_111[3] = {sin(theta_111)*cos(phi_111),sin(theta_111)*sin(phi_111),cos(theta_111)};
  double nsum = 0;

  double coeffs[4]; // gamma_j*E/sqrt(E*gamma_j*E) parameter
  for (i = 0; i<4;++i){
    dot(gamma_j[i], E_111, gamma_E[i]);
    EgammaE[i] = dot(E_111, gamma_E[i]);
    coeffs[i] = gamma_E[i][0] / sqrt(EgammaE[i]); //only the 0 direction matters again
    nsum += sqrt(EgammaE[i]);
  }
  //sum up the coeffs from the (other) <111> directions
  coeff_111 = coeffs[0] + coeffs[1] + coeffs[2];
  coeff_1 = coeffs[3]; //and the one we're pointed in

  n1_ov_n = sqrt(EgammaE[3]) / nsum; //for direction we're pointed in
}

int VelocityModel::drift_velocity(point cart_en, float abse, float q, float& v_over_E, float& dv_dE, vector *velo)
{
  //TODO: doesn't do anything with diffusion parameters v_over_E and dv_dE for now

  //convert to spherical coordinates
  int flag;
  float phi = atan2(cart_en.y,cart_en.x);
  float theta = acos(cart_en.z  );

  if (q>0){
    flag = hole_velocity(abse, theta, phi, velo );
  }
  else{
    flag = electron_velocity(abse, theta, phi, velo );
  }

  if (flag == -1){
    velo->x = 0;
    velo->y = 0;
    velo->z = 0;
    return 0;
  }

  return 0;

}

int VelocityModel::hole_velocity(float field, float theta, float phi, vector *velo){

  float v_100 = h100.get_velocity(field);
  float v_111 = h111.get_velocity(field);

  if (v_100 == 0){
    return -1;
  }

  vector velo_local;

  float v_rel = v_111 / v_100;

  //Numbers from Bruyneel NIM A 569 (2006) p. 764 equations 24-26
  float k_0_0 =  9.2652;
  float k_0_1 =  -26.3467;
  float k_0_2 =  29.6137;
  float k_0_3 =  -12.3689;

  float k_0 = k_0_0 + k_0_1*v_rel + k_0_2*pow(v_rel,2)  + k_0_3* pow(v_rel,3);
  float lambda_k0 = -0.01322 * k_0 + 0.41145*pow(k_0,2) - 0.23657 * pow(k_0,3) + 0.04077*pow(k_0,4);
  float omega_k0 = 0.006550*k_0 - 0.19946*pow(k_0,2) + 0.09859*pow(k_0,3) - 0.01559*pow(k_0,4);

  //in rotated coordinates aligned with electric field
  velo_local.x = v_100 * (1.- lambda_k0*( pow(sin(theta),4) * pow(sin(2.*phi),2) + pow(sin(2.*theta),2) ) );
  velo_local.y = v_100 * omega_k0 * (2.*pow(sin(theta),3)*cos(theta)*pow(sin(2.*phi),2) + sin(4.*theta) );
  velo_local.z = v_100 * omega_k0 * pow(sin(theta),3)*sin(4.*phi);

  //rotate back to crystal coordates
  velo->x = cos(phi)*sin(theta)*velo_local.x + cos(phi)*cos(theta)*velo_local.y - sin(phi)*velo_local.z;
  velo->y = sin(phi)*sin(theta)*velo_local.x + sin(phi)*cos(theta)*velo_local.y + cos(phi)*velo_local.z;
  velo->z = cos(theta)*velo_local.x - sin(theta)*velo_local.y;

  return 0;
}

int VelocityModel::electron_velocity(float field, float theta, float phi,vector *velo){
  //Phenomenological approach from Nathan (Phys Rev 130 (6)), via Mihailescu (NIM A 447 p 350)
  //Best described in a GERDA thesis: https://mediatum.ub.tum.de/doc/701884/701884.pdf
  int i;
  double n1,n2,v_111_x, nsum=0, nj;
  double A_E, R_E;
  double gamma_E[4][3], EgammaE[4];

  double v_100 = e100.get_velocity(field);
  double v_111 = e111.get_velocity(field);
  if (v_100 == 0){
    return -1;
  }

  double E_0[3] = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};

  velo->x = 0, velo->y = 0, velo->z = 0;

  A_E = v_100 / Gamma_0;

  v_111_x = v_111*sin(M_PI/2 - acos(sqrt(2./3)))*cos(M_PI/4);
  n2 = (v_111_x/A_E - coeff_1)/(coeff_111 - 3*coeff_1);
  n1 = 1.-3.*n2;
  R_E = (n1-.25) / (  n1_ov_n -0.25);
  // printf("field: %f R_E: %f\n", field, R_E);

  nsum = 0;
  for (i = 0; i<4;++i){
    dot(gamma_j[i], E_0, gamma_E[i]);
    EgammaE[i] = dot(E_0, gamma_E[i]);
    nsum += sqrt(EgammaE[i]);
  }
  //-= to push electrons in correct direction
  for (i = 0; i<4;++i){
    nj = R_E* (sqrt(EgammaE[i]) / nsum  - 0.25) + 0.25;
    velo->x -= A_E *( nj ) *  gamma_E[i][0] / sqrt(EgammaE[i]);
    velo->y -= A_E *( nj ) *  gamma_E[i][1] / sqrt(EgammaE[i]);
    velo->z -= A_E *( nj ) *  gamma_E[i][2] / sqrt(EgammaE[i]);
  }

  return 0;
}

} // namespace Siggen
