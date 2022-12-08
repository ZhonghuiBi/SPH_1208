#include "Boundary_Condition.h"

Boundary::Boundary()
{
	allv = 0.0;
	inv = 0.0;
	outv = 0.0;
	Periodic[0] = Periodic[1] = Periodic[2] = false;
	inDensity = 0.0;
	outDensity = 0.0;
	allDensity = 0.0;
	InOutFlow = 0;
	InFlowLoc1 = 0.0;
	InFlowLoc2 = 0.0;
	InFlowLoc3 = 0.0;
	OutFlowLoc = 0.0;
	cellfac = 3.9;
	inoutcounter = 0;
	MassConservation = false;
}

