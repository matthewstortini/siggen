#ifndef Siggen_Field
#define Siggen_Field

#include <string>

namespace Siggen
{

class EFieldPoint
{
  private:
    float voltage;
    cyl_pt field;
    float abs_field;

  public:
    EFieldPoint(){
      voltage=0;
      abs_field = 0;
      field.r = 0;
      field.phi=0;
      field.z=0;
    }

    inline cyl_pt get_field(){return field;}
    inline void set_field(cyl_pt new_field){field = new_field;}
    inline float r(){return field.r;}
    inline float phi(){return field.phi;}
    inline float z(){return field.z;}

    // friend ostream &operator<<( ostream &output, const EFieldPoint &F ) {
    //  output << "V : " << D.feet << " I : " << D.inches;
    //  return output;
    // }

    friend std::istream &operator>>( std::istream  &input, EFieldPoint &F ) {
     input >> F.voltage >> F.abs_field >> F.field.r >> F.field.z;
     F.field.phi=0;
     return input;
    }
    //overloaded operator for checking if efld is 0.
    bool operator==(const float& rhs)const{
      if (this->field.r==rhs && this->field.z==rhs && this->field.phi ==rhs ) return true;
      return false;
    }

    void operator= (const float& rhs){
      this->field.r = rhs;
      this->field.phi = rhs;
      this->field.z = rhs;
      this->voltage = rhs;
    }

    EFieldPoint& operator+= (const EFieldPoint& rhs){
      this->field.r += rhs.field.r ;
      this->field.phi += rhs.field.phi;
      this->field.z += rhs.field.z;
      this->voltage += rhs.voltage;
      return *this;
    }

    EFieldPoint& operator*= (float x){
      this->field.r *= x;
      this->field.phi *= x;
      this->field.z *= x;
      this->voltage *= x;
      return *this;
    }
};
inline EFieldPoint operator* (float x, EFieldPoint efld){return efld *= x;}
inline EFieldPoint operator* (EFieldPoint efld, float x){return efld *= x;}

template<class PointType>
class Field
{
	private:
    int rlen;
    int zlen;
    int philen;
    float rstep;
    float zstep;
    float phistep;

    float rmin;
    float zmin;
    float phimin;

    int last_ret =-99;
    cyl_pt  last_pt;
    cyl_int_pt last_ipt;

    std::vector<PointType> field_data;

    PointType& get_point(int ir, int iphi, int iz);

	public:
    Field(){}
    void set_grid(int rlen, int zlen, int philen,
          float rstep, float zstep, float phistep,
          float rmin, float zmin, float phimin);
    void set_grid(int rlen, int zlen,
          float rstep, float zstep,
          float rmin, float zmin);

    //Read in from text file
    template<class GeometryType>
    int read_data(std::string file_name,  GeometryType& geom);

    //Interpolators...
    void trilinear_interpolate(float a, float b, float c, float out[2][2][2]);
    void bilinear_interpolate(float a, float b, float out[2][2]);
    void position_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2]);
    void position_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2][2]);
    void interp_position2D(cyl_pt pt, cyl_int_pt ipt, PointType& d);

    template<class GeometryType>
    int get_point_interp(cyl_pt cyl, PointType& data, GeometryType& geometry);

    //for accessing PointType data
    PointType& operator()(int ir, int iz);
    PointType& operator()(int ir, int iphi, int iz);

    // int field_exists(cyl_pt)
    template<class GeometryType>
    int field_exists(cyl_pt pt, GeometryType& geometry);
    template<class GeometryType>
    int nearest_field_grid_index(cyl_pt pt, cyl_int_pt *ipt, GeometryType& geometry);
    cyl_int_pt field_grid_index(cyl_pt pt);

};

} // namespace Siggen

#include "Field.impl"
#endif
