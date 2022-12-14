#pragma once
#ifndef SPH_DOMAIN_H
#define SPH_DOMAIN_H

#include <stdio.h>    // for NULL
#include <algorithm>  // for min,max
#include <iostream>
#include <string>




#include <fstream>
#include "Particle.h"
#include "Functions.h"
#include "Boundary_Condition.h"
//#include "../External/string.h"

//C++ Enum used for easiness of coding in the input files
enum Kernels_Type { Qubic_Spline = 0, Quintic = 1, Quintic_Spline = 2 };
enum Viscosity_Eq_Type { Morris = 0, Shao = 1, Incompressible_Full = 2, Takeda = 3 };
enum Gradient_Type { Squared_density = 0, Multiplied_density = 1 };


class Domain
{
public:
	typedef void(*PtVel) (Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry);
	typedef void(*PtOut) (Particle * Particles, double & Prop1, double & Prop2, double & Prop3);
	typedef void(*PtDom) (Domain & dom);
	// Constructor
	Domain();

	// Destructor
	~Domain();

	// Domain Part
	void AddSingleParticle(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed);		//Add one particle
	void AddBoxLength(int tag, Vec3_t const &V, double Lx, double Ly, double Lz, double r, double Density,
		double h, int type, int rotation, bool random, bool Fixed);									//Add a cube of particles with a defined dimensions
	void AddBoxNo(int tag, Vec3_t const &V, size_t nx, size_t ny, size_t nz, double r, double Density,
		double h, int type, int rotation, bool random, bool Fixed);									//Add a cube of particles with a defined numbers
	void DelParticles(int const & Tags);					//Delete particles by tag
	void CheckParticleLeave();													//Check if any particles leave the domain, they will be deleted

	void YZPlaneCellsNeighbourSearch(int q1);						//Create pairs of particles in cells of XZ plan
	void MainNeighbourSearch();									//Create pairs of particles in the whole domain
	void StartAcceleration(Vec3_t const & a = Vec3_t(0.0, 0.0, 0.0));	//Add a fixed acceleration such as the Gravity
	void PrimaryComputeAcceleration();									//Compute the solid boundary properties
	void LastComputeAcceleration();									//Compute the acceleration due to the other particles
	void CalcForce11(Particle * P1, Particle * P2);	//Calculates the contact force between fluid-fluid particles
	void CalcForce2233(Particle * P1, Particle * P2);	//Calculates the contact force between soil-soil/solid-solid particles
	void CalcForce12(Particle * P1, Particle * P2);	//Calculates the contact force between fluid-solid particles
	void CalcForce13(Particle * P1, Particle * P2);	//Calculates the contact force between fluid-soil particles
	void Move(double dt);										//Move particles

	void Solve(double tf, double dt, double dtOut, std::string const TheFileKey, size_t maxidx);		///< The solving function

	void CellInitiate();															//Find the size of the domain as a cube, make cells and HOCs
	void ListGenerate();															//Generate linked-list
	void CellReset();															//Reset HOCs and particles' LL to initial value of -1

	void WriteParaview(std::string const FileKey);					//Save a XDMF file for the visualization


	void InFlowBCLeave();
	void InFlowBCFresh();
	void WholeVelocity();

	void Kernel_Set(Kernels_Type const & KT);
	void Viscosity_Eq_Set(Viscosity_Eq_Type const & VQ);
	void Gradient_Approach_Set(Gradient_Type const & GT);
	// Data
	Array <Particle*>				Particles; 	///< Array of particles
	double					R;		///< Particle Radius in addrandombox

	double					sqrt_h_a;				//Coefficient for determining Time Step based on acceleration (can be defined by user)

	int 					Dimension;    	///< Dimension of the problem

	double					MuMax;		///< Max Dynamic viscosity for calculating the timestep
	double					CsMax;		///< Max speed of sound for calculating the timestep

	Vec3_t					Gravity;       	///< Gravity acceleration


	Vec3_t                 			TRPR;		///< Top right-hand point at rear of the domain as a cube
	Vec3_t                  			BLPF;           ///< Bottom left-hand point at front of the domain as a cube
	Vec3_t                  			CellSize;      	///< Calculated cell size according to (cell size >= 2h)
	int		                		CellNo[3];      ///< No. of cells for linked list
	double 					hmax;		///< Max of h for the cell size  determination
	Vec3_t                 	DomSize;	///< Each component of the vector is the domain size in that direction if periodic boundary condition is defined in that direction as well
	double					rhomax;

	int						*** HOC;	///< Array of "Head of Chain" for each cell

	size_t					SWIType;	///< Selecting variable to choose Soil-Water Interaction type
	bool					FSI;		///< Selecting variable to choose Fluid-Structure Interaction

	double 					XSPH;		///< Velocity correction factor
	double 					InitialDist;	///< Initial distance of particles for Inflow BC

	double					AvgVelocity;	///< Average velocity of the last two column for x periodic constant velocity

	size_t					Nproc;		///< No of threads which are going to use in parallel calculation
	omp_lock_t 				dom_lock;	///< Open MP lock to lock Interactions array
	Boundary				BC;
	PtOut					UserOutput;
	PtVel 					InCon;
	PtVel 					OutCon;
	PtVel 					AllCon;
	Vec3_t					DomMax;
	Vec3_t					DomMin;
	PtDom					GeneralBefore;	///< Pointer to a function: to modify particles properties before CalcForce function
	PtDom					GeneralAfter;	///< Pointer to a function: to modify particles properties after CalcForce function
	size_t					Scheme;		///< Integration scheme: 0 = Modified Verlet, 1 = Leapfrog

	Array<Array<std::pair<size_t, size_t> > >	SMPairs;
	Array<Array<std::pair<size_t, size_t> > >	NSMPairs;
	Array<Array<std::pair<size_t, size_t> > >	FSMPairs;
	Array< size_t > 				FixedParticles;
	Array< size_t >				FreeFSIParticles;

	Array<std::pair<size_t, size_t> >		Initial;
	Mat3_t I;
	std::string					OutputName[3];


private:
	void Periodic_X_Correction(Vec3_t & x, double const & h, Particle * P1, Particle * P2);		//Corrects xij for the periodic boundary condition
	void AdaptiveTimeStep();		//Uses the minimum time step to smoothly vary the time step

	void PrintInput(std::string const FileKey);		//Print out some initial parameters as a file
	void InitialChecks();		//Checks some parameter before proceeding to the solution
	void TimestepCheck();		//Checks the user time step with CFL approach

	size_t					VisEq;					//Choose viscosity Eq based on different SPH discretisation
	size_t					KernelType;			//Choose a kernel
	size_t					GradientType;		//Choose a Gradient approach 1/Rho i^2 + 1/Rho j^2 or 1/(Rho i * Rho j)
	double 					Cellfac;				//Define the compact support of a kernel

	double					Time;    				//Current time of simulation at each solving step
	double					deltat;					//Time Step
	double					deltatmin;			//Minimum Time Step
	double					deltatint;			//Initial Time Step

};

//#include "Interaction.cpp"
//#include "Domain.cpp"

#endif // SPH_DOMAIN_H
