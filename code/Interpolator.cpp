#include "Interpolator.h"


namespace Siggen
{
  void Interpolator::trilinear_interpolate(float a, float b, float c, float out[2][2][2]){
    out[0][0][0] = (1.0 - a) * (1.0 - b) * (1.0 - c);
    out[0][1][0] = (1.0 - a) *        b  * (1.0 - c);
    out[1][0][0] =        a  * (1.0 - b) * (1.0 - c);
    out[1][1][0] =        a  *        b  * (1.0 - c);

    out[0][0][1] = (1.0 - a) * (1.0 - b) * c;
    out[0][1][1] = (1.0 - a) *        b  * c;
    out[1][0][1] =        a  * (1.0 - b) * c;
    out[1][1][1] =        a  *        b  * c;
  }

  void Interpolator::bilinear_interpolate(float a, float b,  float out[2][2]){

    out[0][0] = (1.0 - a) * (1.0 - b);
    out[0][1] = (1.0 - a) *         b;
    out[1][0] =        a  * (1.0 - b);
    out[1][1] =        a  *         b;
  }

}//namespace Siggen
