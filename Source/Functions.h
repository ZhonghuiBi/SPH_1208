#ifndef SPH_SPECIAL_FUNCTIONS_H
#define SPH_SPECIAL_FUNCTIONS_H
#pragma once
#include "../External/MatVec.h"
#include "../External/array.hpp"

double Kernel(size_t const & Dim, size_t const & KT, double const & q, double const & h);

double GradKernel(size_t const & Dim, size_t const & KT, double const & q, double const & h);

double LaplaceKernel(size_t const & Dim, size_t const & KT, double const & q, double const & h);

double SecDerivativeKernel(size_t const & Dim, size_t const & KT, double const & q, double const & h);

double EOS(size_t const & EQ, double const & Cs0, double const & P00, double const & Density, double const & Density0);

double SoundSpeed(size_t const & EQ, double const & Cs0, double const & Density, double const & Density0);

double DensitySolid(size_t const & EQ, double const & Cs0, double const & P00, double const & Pressure, double const & Density0);

void   Seepage(size_t const & ST, double const & k, double const & k2, double const & mu, double const & rho, double & SF1, double & SF2);

void   Viscous_Force(size_t const & VisEq, Vec3_t & VI, double const & Mu, double const & di, double const & dj, double const & GK, Vec3_t const & vab,
	size_t const & Dimension, double const & KernelType, double const & rij, double const & h, Vec3_t const & xij, Vec3_t const & vij);

void   Rotation(Mat3_t Input, Mat3_t & Vectors, Mat3_t & VectorsT, Mat3_t & Values);

Mat3_t abab(Mat3_t const & A, Mat3_t const & B);


//#include "Functions.cpp"

#endif // SPH_SPECIAL_FUNCTIONS_H
