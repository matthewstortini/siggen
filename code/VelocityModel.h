#ifndef Siggen_VelocityModel
#define Siggen_VelocityModel

#include "point.h"

namespace Siggen
{

class AxisModel
{
  private:
    float mu_0, beta, E_0, mu_n = 0.;

  public:
    AxisModel(float mu_0_in, float beta_in, float E_0_in, float mu_n_in=0.);
    void set_params(float mu_0_in, float beta_in, float E_0_in, float mu_n_in=0.);
    float get_velocity(float E);

};

class VelocityModel
{
	private:
    // int hole_velocity(float field, float theta, float phi, vector& velo);
    // int electron_velocity(float field, float theta, float phi, vector& velo);

    //Axis velocity models.
    AxisModel h100;
    AxisModel h111;
    AxisModel e100;
    AxisModel e111;

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
