#pragma once
#include "../External/Vec3D.hpp"
#include "../External/array.hpp"

class Boundary
{
public:
	// Data
	double	inDensity;	///< Apply a certain density to inflow particles
	double	outDensity;	///< Apply a certain density to outflow particles
	double	allDensity;	///< Apply a certain density to outflow particles

	Vec3_t 	inv;		///< Apply a certain velocity to inflow particles
	Vec3_t 	outv;		///< Apply a certain velocity to outflow particles
	Vec3_t	allv;		///< Apply a certain velocity to all particle

	bool 	Periodic[3];	///< Considering periodic in all directions => 0=X, 1=Y, 2=Z

	int 	InOutFlow;	///< Considering inflow in all directions  by adding and deleting particles=> [0]=X, [1]=Y, [2]=Z and 0=none, 1=-
	double	InFlowLoc1;
	double	InFlowLoc2;
	double	InFlowLoc3;
	double	OutFlowLoc;
	double	cellfac;
	int		inoutcounter;
	bool	MassConservation;

	Array <int>	OutPart;
	Array <int>	InPart;

	Boundary();
};
