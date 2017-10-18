#include "VectorXYZ.hpp"

vectorXYZ::vectorXYZ(void) {  // Void Constructor
	x=0.0;
	y=0.0;
	z=0.0;
}
vectorXYZ::vectorXYZ(double xi,double yi,double zi) {  // Constructor
	x=xi;
	y=yi;
	z=zi;
}
vectorXYZ::vectorXYZ(const vectorXYZ &a){   // Copy Vector Constructor
	x=a.x;
	y=a.y;
	z=a.z;
}

vectorXYZ &vectorXYZ::operator=(vectorXYZ a){  // Copy Vector
	x=a.x;
	y=a.y;
	z=a.z;
	return *this;
}
vectorXYZ &vectorXYZ::operator+=(vectorXYZ a) 
{
    x += a.x;
    y += a.y;
    z += a.z;
    return *this;
}
vectorXYZ &vectorXYZ::operator-=(vectorXYZ a) 
{
    x -= a.x;
    y -= a.y;
    z -= a.z;
    return *this;
}
vectorXYZ &vectorXYZ::operator*=(double scalar) 
{
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}
vectorXYZ &vectorXYZ::operator*=(vectorXYZ a) 
{
    x *= a.x;
    y *= a.y;
    z *= a.z;
    return *this;
}
vectorXYZ &vectorXYZ::operator/=(vectorXYZ a) 
{
    x /= a.x;
    y /= a.y;
    z /= a.z;
    return *this;
}


vectorXYZ  vectorXYZ::operator+(vectorXYZ a){  // Addition
  return vectorXYZ(*this)+=a;
}
vectorXYZ  vectorXYZ::operator-(vectorXYZ a){  // Subtraction
  return vectorXYZ(*this)-=a;
}
vectorXYZ  vectorXYZ::operator*(vectorXYZ a){  // Product
  return vectorXYZ(*this)*=a;
}
vectorXYZ  vectorXYZ::operator/(vectorXYZ a){  // Division
  return vectorXYZ(*this)/=a;
}


vectorXYZ::~vectorXYZ(){}    // Destructor


vectorXYZ operator*( double scalar,vectorXYZ a) {
  return vectorXYZ(a) *= scalar;
} 
vectorXYZ operator*(vectorXYZ a, double scalar) {
  return vectorXYZ(a) *= scalar;
} 

ostream &operator<<(ostream &out, vectorXYZ a)     //output
{
        out<<a.x<<" "<<a.y<<" "<<a.z;
        return out;
}

double scalar(vectorXYZ a,vectorXYZ b) {
  double result=a.x*b.x+a.y*b.y+a.z*b.z;
  
  return result;
}

vectorXYZ cross(vectorXYZ a,vectorXYZ b) {
  vectorXYZ result;

  result.x = a.y*b.z-a.z*b.y;
  result.y = a.z*b.x-a.x*b.z;
  result.z = a.x*b.y-a.y*b.x;

  return result;
}

istream &operator>>(istream &in, vectorXYZ &a)     //input
{
  in>>a.x;
  in>>a.y;
  in>>a.z;  
  return in;
}

bool operator==(vectorXYZ a,vectorXYZ b)
{
  return (a.x==b.x && a.y==b.y && a.z==b.z);
}
