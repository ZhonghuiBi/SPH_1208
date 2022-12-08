/*******************************************************************************
 *
 * Descritpion: 3D vector class
 *
 ******************************************************************************/
#pragma once
#include <iostream>
//using namespace std;
template<typename T>
class Vec3D
{
public:

	// Constructors
	Vec3D();
	Vec3D(const T& x, const T& y, const T& z);
	Vec3D(const T& xyz);	// x = y = z = xyz
	Vec3D(const T xyz[]);
	Vec3D(const Vec3D<T>& o);

	//Vec3D& operator=(const T& xyz);



	// Math operators (=, +, -, *, /, etc...) with self assignment
	Vec3D<T>& operator=(const Vec3D<T>& o);
	Vec3D<T>& operator+=(const Vec3D<T> &o);
	Vec3D<T>& operator-=(const Vec3D<T> &o);
	Vec3D<T>& operator*=(const T& scalar);
	Vec3D<T>& operator/=(const T& scalar);

	// Math operators without self-assignment (Create a temporary variable)
	Vec3D<T> operator+(const Vec3D<T>& o) const;
	Vec3D<T> operator-(const Vec3D<T>& o) const;
	Vec3D<T> operator*(const T& scalar) const;
	Vec3D<T> operator/(const T& scalar) const;
	template  <class ElemType>
	friend std::ostream& operator<<(std::ostream &os, const Vec3D<ElemType>& v); //向输出流输出

	T &operator()(int num) const;

	// Boolean operators
	bool operator==(const Vec3D<T>& o) const;

	// Vector operations
	T length2() const;
	T norm() const;

	Vec3D<T>& Normalize();	// Returns a reference to "this"
	Vec3D<T> getNormalizedVec() const;

	T dot(const Vec3D<T>& o) const;
	Vec3D<T> cross(const T& o) const;
	Vec3D<T>& applyCross(const T& o);	// a = cross(a,b). Returns a reference to "this"
public:
	// Member variables
	mutable T	x;
	mutable T	y;
	mutable T	z;
};

// Non-member operators
template<typename T> Vec3D<T> operator*(const T& scalar, const Vec3D<T>& v);
template<typename T> Vec3D<T> operator/(const T& scalar, const Vec3D<T>& v);

typedef Vec3D<double>			Vec3_t;
typedef Vec3D<float>			Vec3f;
typedef Vec3D<int>				Vec3i;
typedef Vec3D<unsigned int>		Vec3ui;
typedef Vec3D<long>				Vec3l;
typedef Vec3D<unsigned long>	Vec3ul;
typedef Vec3D<short>			Vec3s;
typedef Vec3D<unsigned short>	Vec3us;


//------------------------------------------------------------------------------
// Constructors(构造函数类外实现)
//------------------------------------------------------------------------------
template<typename T>
Vec3D<T>::Vec3D()
{}

template<typename T>
Vec3D<T>::Vec3D(const T& x, const T& y, const T& z)
	: x(x), y(y), z(z)
{}

template<typename T>
Vec3D<T>::Vec3D(const T& xyz)
	: x(xyz), y(xyz), z(xyz)
{}

template<typename T>
Vec3D<T>::Vec3D(const T xyz[])
	: x(xyz[0]), y(xyz[1]), z(xyz[2])
{}

template<typename T>
Vec3D<T>::Vec3D(const Vec3D<T>& o)
	: x(o.x), y(o.y), z(o.z)
{}

//------------------------------------------------------------------------------
// Math operators with self assignment
//------------------------------------------------------------------------------
template<typename T>
Vec3D<T>& Vec3D<T>::operator=(const Vec3D<T>& o)
{
	x = o.x;
	y = o.y;
	z = o.z;

	return *this;
}

template<typename T>
Vec3D<T>& Vec3D<T>::operator+=(const Vec3D<T> &o)
{
	x += o.x;
	y += o.y;
	z += o.z;

	return *this;
}

template<typename T>
Vec3D<T>& Vec3D<T>::operator-=(const Vec3D<T> &o)
{
	x -= o.x;
	y -= o.y;
	z -= o.z;

	return *this;
}

template<typename T>
Vec3D<T>& Vec3D<T>::operator*=(const T& scalar)
{
	x *= scalar;
	y *= scalar;
	z *= scalar;

	return *this;
}

template<typename T>
Vec3D<T>& Vec3D<T>::operator/=(const T& scalar)
{
	x /= scalar;
	y /= scalar;
	z /= scalar;

	return *this;
}

template <class ElemType>
std::ostream& operator<<(std::ostream& out, const Vec3D<ElemType>& v)
{
	out << v.x << "\t" << v.y << "\t" << v.z;
	return out;
}
//------------------------------------------------------------------------------
// Math operators without self assignment (Create a temporary variable)
//------------------------------------------------------------------------------
template<typename T>
Vec3D<T> Vec3D<T>::operator+(const Vec3D<T>& o) const
{
	Vec3D<T> result(*this);
	result += o;

	return result;
}

template<typename T>
Vec3D<T> Vec3D<T>::operator-(const Vec3D<T>& o) const
{
	Vec3D<T> result(*this);
	result -= o;

	return result;
}

template<typename T>
Vec3D<T> Vec3D<T>::operator*(const T& scalar) const
{
	Vec3D<T> result(*this);
	result *= scalar;

	return result;
}

template<typename T>
Vec3D<T> Vec3D<T>::operator/(const T& scalar) const
{
	Vec3D<T> result(*this);
	result /= scalar;

	return result;
}

template<typename T>
T& Vec3D<T>::operator()(int num) const
{
	if (num == 0) return x;
	if (num == 1) return y;
	if (num == 2) return z;
}

template<typename T>
Vec3D<T> operator*(const T& scalar, const Vec3D<T>& v)
{
	Vec3D<T> result(v);
	result *= scalar;

	return result;
}

template<typename T>
Vec3D<T> operator/(const T& scalar, const Vec3D<T>& v)
{
	Vec3D<T> result(v);
	result /= scalar;

	return result;
}

//------------------------------------------------------------------------------
// Boolean operators
//------------------------------------------------------------------------------
template<typename T>
bool Vec3D<T>::operator==(const Vec3D<T>& o) const
{
	return (x == o.x) && (y == o.y) && (z == o.z);
}

//------------------------------------------------------------------------------
// Vector operations
//------------------------------------------------------------------------------
template<typename T>
T Vec3D<T>::length2() const
{
	return x * x + y * y + z * z;
}

template<typename T>
T Vec3D<T>::norm() const
{
	return sqrt(length2());
}

template<typename T>
Vec3D<T>& Vec3D<T>::Normalize()
{
	T length = length();

	(*this) /= length;

	return *this;
}

template<typename T>
Vec3D<T> Vec3D<T>::getNormalizedVec() const
{
	Vec3D<T> result(*this);
	result.Normalize();

	return result;
}

template<typename T>
T Vec3D<T>::dot(const Vec3D<T>& o) const
{
	return x * o.x + y * o.y + z * o.z;
}

template<typename T>
Vec3D<T> Vec3D<T>::cross(const T& o) const
{
	Vec3D<T> result;
	result.x = y * o.z - z * o.y;
	result.y = z * o.x - x * o.z;
	result.z = x * o.y - y * o.x;

	return result;
}

template<typename T>
Vec3D<T>& Vec3D<T>::applyCross(const T& o)
{
	T x0 = x;
	T y0 = y;

	x = y * o.z - z * o.y;
	y = z * o.x - x0 - o.z;
	z = x0 * o.y - y0 * o.x;

	return *this;
}
