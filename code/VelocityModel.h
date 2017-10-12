#ifndef Siggen_VelocityModel
#define Siggen_VelocityModel

#include "point.h"

namespace Siggen
{

class AxisModel
{
  private:
    float mu_0, beta, E_0, mu_n = 0;

  public:
    AxisModel(float mu_0, float beta, float E_0);
    AxisModel(float mu_0, float beta, float E_0, float mu_n);

    float get_velocity(float E);
    void set_params(float mu_0_in, float beta_in, float E_0_in, float mu_n_in=0);
    // void set_params(float mu_0, float beta, float E_0, float mu_n);
};

class VelocityModel
{
	private:
    // int hole_velocity(float field, float theta, float phi, vector& velo);
    // int electron_velocity(float field, float theta, float phi, vector& velo);

    //Axis velocity models.  Default parameters from Bruyneel NIM A 569 p 764

    //hole velocity based on Reggiani data (Phys Rev B 16 (6))
    AxisModel h100 = AxisModel(66333.,0.744,181.);
    AxisModel h111 = AxisModel(107270.,0.580,100.);

    //hole velocity based on Mihailescu data (NIM A 447 p 350)
    AxisModel e100 = AxisModel(40180.,0.72,493.,589.);
    AxisModel e111 = AxisModel(42420.,0.87,251.,62.);

    //gamma matrices for electron velocity calculation
    double gamma_j[4][3][3];
    double Gamma_0; //depends only on values of gamma_0 matrix, gives us A(E)
    double coeff_111, coeff_1, n1_ov_n; //for calculating R(E)

	public:
    VelocityModel ();
    int hole_velocity(float field, float theta, float phi, vector* velo);
    int electron_velocity(float field, float theta, float phi, vector* velo);
    int drift_velocity(point cart_en, float abse, float q, float& v_over_E, float& dv_dE, vector *velo);

    inline void set_holes(float mu0_100, float beta_100, float E_0_100,
                          float mu0_111, float beta_111, float E_0_111){
                          h100.set_params(mu0_100,beta_100,E_0_100);
                          h111.set_params(mu0_111,beta_111,E_0_111);
                         }
    inline void set_electrons(float mu0_100, float beta_100, float E_0_100, float mu_n_100,
                            float mu0_111, float beta_111, float E_0_111, float mu_n_111){
                            e100.set_params(mu0_100,beta_100,E_0_100,mu_n_100);
                            e111.set_params(mu0_111,beta_111,E_0_111,mu_n_111);
                          }


};

} // namespace Siggen
#endif
