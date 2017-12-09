#ifndef Siggen_Interpolator
#define Siggen_Interpolator

namespace Siggen
{

class Interpolator
{

	public:
		Interpolator() {};

    void trilinear_interpolate(float a, float b, float c, float out[2][2][2]);
    void bilinear_interpolate(float a, float b, float out[2][2]);

};

} // namespace Siggen


#endif
