#include "Particle.h"
#include "Functions.h"

Particle::Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0, bool Fixed)
{
	//Correction step for the Modified Verlet Algorithm
	ct = 0;
	a = 0.0;
	x = x0;
	//Prosity
	n = 0.0;
	//Initial Prosity
	n0 = 0.0;
	//Permeability
	k = 0.0;
	//Second Permeability for the Forchheimer Eq
	k2 = 0.0;
	//Speed of sound
	Cs = 0.0;
	//background pressure for equation of state
	P0 = 0.0;
	//Selecting variable to choose an equation of state
	PresEq = 0;
	//Dynamic viscosity coefficient of the fluid particle
	Alpha = 0.0;
	//Dynamic viscosity coefficient of the fluid particle
	Beta = 0.0;
	//Velocity of the particle n+1/2 (Leapfrog)
	va = 0.0;
	//Velocity of the particle n-1 (Modified Verlet)
	vb = 0.0;
	//Velocity of the fixed particle for no-slip BC
	NSv = 0.0;
	//Velocity of the fixed particle for no-slip BC in FSI
	FSINSv = 0.0;
	//Velocity of the particle n+1
	v = v0;
	//Mean Velocity of neighbor particles for updating the particle position (XSPH)
	VXSPH = 0.0;
	//Tensile instability factor
	TI = 0.0;
	//Tensile instability power
	TIn = 4.0;
	//Initial distance of particles for calculation of tensile instability
	TIInitDist = 0.0;
	//Density of the particle n+1/2 (Leapfrog)
	Densitya = 0.0;
	//Density of the particle n-1 (Modified Verlet)
	Densityb = 0.0;
	//Density of the particle n+1
	Density = Density0;
	//Reference Density of Particle
	RefDensity = Density0;
	//Mass of the particle
	Mass = Mass0;
	//Mass coefficient for fixed particles to avoid leaving particles
	FPMassC = 1.0;
	IsFree = !Fixed;
	//Smoothing length of the particle
	h = h0;
	//Pressure of the particle n+1
	Pressure = 0.0;
	//Pressure of the particle n+1 in FSI
	FSIPressure = 0.0;
	ID = Tag;
	//Current cell No for the particle (linked-list)
	CC[0] = CC[1] = CC[2] = 0;
	//Linked-List variable to show the next particle in the list of a cell
	LL = 0;
	//Summation of mb/db*Wab for neighbour particles of the particle a (for Shepard filter)
	ZWab = 0.0;
	//Summation of mb*Wab for neighbour particles of the particle a (for Shepard filter)
	SumDen = 0.0;
	//Rate of density change in time based on state equations n
	dDensity = 0.0;
	//Global shear rate for fluids
	ShearRate = 0.0;
	//Reference Dynamic viscosity coefficient
	//Dynamic viscosity coefficient of the fluid particle
	MuRef = Mu = 0.0;
	//Non-Newtonian viscosity method
	VisM = 0;
	//Yield stress for Bingham fluids
	T0 = 0.0;
	//Normalization value for Bingham fluids
	m = 300.0;
	//Summation of the kernel value for neighbour particles
	SumKernel = 0.0;
	//Summation of the kernel value for neighbour particles in FSI
	FSISumKernel = 0.0;
	//Shear modulus
	G = 0.0;
	// Bulk modulus
	K = 0.0;
	//an Integer value to identify the particle material type: 1 = Fluid, 2 = Solid, 3 = Soil
	Material = 0;
	//Failure criteria
	Fail = 0;
	//Cohesion
	c = 0.0;
	//Friction angel
	phi = 0.0;
	//Dilation angel
	psi = 0.0;
	//effective particle size
	d = 0.0;
	//Tensile yield stress
	Sigmay = 0.0;
	//No-Slip BC
	NoSlip = false;
	//Shepard Filter for the density
	Shepard = false;
	//Check the particle if it is in-flow or out-flow or not
	InOut = 0;
	//to initialize the integration scheme
	FirstStep = true;
	//Volume of a particle
	V = Mass / RefDensity;
	//Density of water or any other fluids
	RhoF = 0.0;
	//Check the particle if it is Saturated or not
	IsSat = false;
	//Check the particle Saturation at each time step
	SatCheck = false;
	//Cycle number for shepard correction
	ShepardStep = 40;
	//Count number of contributing particles
	ShepardCounter = 0;
	//Velocity derivative for surface erosion
	S = 0.0;
	//If yes, it will calculate porosity and permeability based on new calculated porosity
	VarPorosity = false;
	//Selecting variable to choose a Seepage method
	SeepageType = 0;
	//Large eddy simulation using sub-particle scale
	LES = false;
	//shear component for LES
	SBar = 0.0;
	//Coefficient of Smagorinsky-Lilly model
	CSmag = 0.17;


	//Total Strain n-1 (Modified Verlet)
	MatVec::set_to_zero(Strainb);
	//Total Strain n+1
	MatVec::set_to_zero(Strain);
	//Cauchy stress tensor (Total Stress) n-1 (Modified Verlet)
	MatVec::set_to_zero(Sigmab);
	//Cauchy stress tensor (Total Stress) n+1
	MatVec::set_to_zero(Sigma);
	//set_to_zero(FSISigma);
	//Cauchy stress tensor (Total Stress) n+1/2 (Leapfrog)
	MatVec::set_to_zero(Sigmaa);
	//Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n+1
	MatVec::set_to_zero(ShearStress);
	//Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n-1 (Modified Verlet)
	MatVec::set_to_zero(ShearStressb);
	//Tensile Instability stress tensor R
	MatVec::set_to_zero(TIR);
	//Global shear Strain rate tensor n
	MatVec::set_to_zero(StrainRate);
	//Global rotation tensor n
	MatVec::set_to_zero(RotationRate);
	omp_init_lock(&my_lock);

}

void Particle::Move(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, size_t Scheme, Mat3_t I)
{
	if (Scheme == 0)
		Move_MVerlet(I, dt);
	else
		Move_Leapfrog(I, dt);


	//Periodic BC particle position update
	if (Domainsize(0) > 0.0)
	{
		(x(0) > (domainmax(0))) ? x(0) -= Domainsize(0) : x(0);
		(x(0) < (domainmin(0))) ? x(0) += Domainsize(0) : x(0);
	}
	if (Domainsize(1) > 0.0)
	{
		(x(1) > (domainmax(1))) ? x(1) -= Domainsize(1) : x(1);
		(x(1) < (domainmin(1))) ? x(1) += Domainsize(1) : x(1);
	}
	if (Domainsize(2) > 0.0)
	{
		(x(2) > (domainmax(2))) ? x(2) -= Domainsize(2) : x(2);
		(x(2) < (domainmin(2))) ? x(2) += Domainsize(2) : x(2);
	}

}

void Particle::Mat1(double dt)
{
	Pressure = EOS(PresEq, Cs, P0, Density, RefDensity);
	double temp = (StrainRate(0, 0)*StrainRate(0, 0) + 2.0*StrainRate(0, 1)*StrainRate(1, 0) +
		2.0*StrainRate(0, 2)*StrainRate(2, 0) + StrainRate(1, 1)*StrainRate(1, 1) +
		2.0*StrainRate(1, 2)*StrainRate(2, 1) + StrainRate(2, 2)*StrainRate(2, 2));

	ShearRate = sqrt(0.5*temp);
	SBar = sqrt(2.0*temp);

	// LES model
	if (LES)
	{
		Mu = MuRef + RefDensity*pow((CSmag*h), 2.0)*SBar;
	}

	// Bingham viscosity calculation
	if (T0 > 0.0)
	{
		switch (VisM)
		{
		case 0:
			// Bingham
			if (ShearRate != 0.0)
				Mu = MuRef + T0*(1 - exp(-m*ShearRate)) / ShearRate;
			else
				Mu = MuRef + T0*m;
			break;
		case 1:
			// Cross
			Mu = (1000.0*MuRef + MuRef*MuRef*1000.0 / T0*ShearRate) / (1 + 1000.0*MuRef / T0*ShearRate);
			break;
		default:
			std::cout << "Non-Newtonian Viscosity Type No is out of range. Please correct it and run again" << std::endl;
			std::cout << "0 => Bingham" << std::endl;
			std::cout << "1 => Cross" << std::endl;
			abort();
			break;
		}
	}
}


void Particle::Move_MVerlet(Mat3_t I, double dt)
{
	if (FirstStep)
	{
		ct = 30;
		FirstStep = false;
	}

	x += dt*(v + VXSPH) + 0.5*dt*dt*a;

	if (ct == 30)
	{
		if (Shepard && ShepardCounter == ShepardStep)
		{
			if (ZWab > 0.6)
			{
				Densityb = SumDen / ZWab;
				//				Densityb	= Density;
				Density = SumDen / ZWab;
			}
			else
			{
				Densityb = Density;
				Density += dt*dDensity;
			}
		}
		else
		{
			Densityb = Density;
			Density += dt*dDensity;
		}

		vb = v;
		v += dt*a;
	}
	else
	{
		if (Shepard && ShepardCounter == ShepardStep)
		{
			if (ZWab > 0.6)
			{
				Densityb = SumDen / ZWab;
				//				Densityb	= Density;
				Density = SumDen / ZWab;
			}
			else
			{
				double dens = Density;
				Density = Densityb + 2.0*dt*dDensity;
				Densityb = dens;
			}
		}
		else
		{
			double dens = Density;
			Density = Densityb + 2.0*dt*dDensity;
			Densityb = dens;
		}

		Vec3_t temp;
		temp = v;
		v = vb + 2 * dt*a;
		vb = temp;
	}

	switch (Material)
	{
	case 1:
		Mat1(dt);
		break;
	case 2:
		Mat2MVerlet(dt);
		break;
	case 3:
		Mat3MVerlet(I, dt);
		break;
	default:
		std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "1 => Fluid" << std::endl;
		std::cout << "2 => Solid" << std::endl;
		std::cout << "3 => Soil" << std::endl;
		abort();
		break;
	}
	if (ct == 30) ct = 0; else ct++;
	if (ShepardCounter == ShepardStep) ShepardCounter = 0; else ShepardCounter++;
}

void Particle::Mat2MVerlet(double dt)
{
	Pressure = EOS(PresEq, Cs, P0, Density, RefDensity);

	// Jaumann rate terms
	Mat3_t RotationRateT, Stress, SRT, RS;
	MatVec::Trans(RotationRate, RotationRateT);
	MatVec::Mult(ShearStress, RotationRateT, SRT);
	MatVec::Mult(RotationRate, ShearStress, RS);

	// Elastic prediction step (ShearStress_e n+1)
	Stress = ShearStress;
	if (ct == 30)
		ShearStress = dt*(2.0*G*(StrainRate - 1.0 / 3.0*(StrainRate(0, 0) + StrainRate(1, 1) + StrainRate(2, 2))*OrthoSys::I) + SRT + RS) + ShearStress;
	else
		ShearStress = 2.0*dt*(2.0*G*(StrainRate - 1.0 / 3.0*(StrainRate(0, 0) + StrainRate(1, 1) + StrainRate(2, 2))*OrthoSys::I) + SRT + RS) + ShearStressb;
	ShearStressb = Stress;

	if (Fail == 1)
	{
		double J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
			2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
			2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));
		//Scale back
		ShearStress = std::min((Sigmay / sqrt(3.0*J2)), 1.0)*ShearStress;
	}

	Sigma = -Pressure * OrthoSys::I + ShearStress;

	Stress = Strain;
	if (ct == 30)
		Strain = dt*StrainRate + Strain;
	else
		Strain = 2.0*dt*StrainRate + Strainb;
	Strainb = Stress;


	if (Fail > 1)
	{
		std::cout << "Undefined failure criteria for solids" << std::endl;
		abort();
	}
}

void Particle::Mat3MVerlet(Mat3_t I, double dt)
{
	Mat3_t RotationRateT, Stress, SRT, RS;
	double I1, J2, alpha, kf, I1strain;

	// Jaumann rate terms
	MatVec::Trans(RotationRate, RotationRateT);
	MatVec::Mult(Sigma, RotationRateT, SRT);
	MatVec::Mult(RotationRate, Sigma, RS);

	// Volumetric strain
	I1strain = StrainRate(0, 0) + StrainRate(1, 1) + StrainRate(2, 2);

	// Elastic prediction step (Sigma_e n+1)
	Stress = Sigma;
	if (ct == 30)
		Sigma = dt*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate - 1.0 / 3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigma;
	else
		Sigma = 2.0*dt*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate - 1.0 / 3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigmab;
	Sigmab = Stress;

	if (Fail > 1)
	{
		if (I(2, 2) == 0.0)
		{
			// Drucker-Prager failure criterion for plane strain
			alpha = tan(phi) / sqrt(9.0 + 12.0*tan(phi)*tan(phi));
			kf = 3.0 * c / sqrt(9.0 + 12.0*tan(phi)*tan(phi));
		}
		else
		{
			// Drucker-Prager failure criterion for 3D
			alpha = (2.0*  sin(phi)) / (sqrt(3.0)*(3.0 - sin(phi)));
			kf = (6.0*c*cos(phi)) / (sqrt(3.0)*(3.0 - sin(phi)));
		}


		// Bring back stress to the apex of the failure criteria
		I1 = Sigma(0, 0) + Sigma(1, 1) + Sigma(2, 2);
		if ((kf - alpha*I1) < 0.0)
		{
			double Ratio;
			if (alpha == 0.0) Ratio = 0.0; else Ratio = kf / alpha;
			Sigma(0, 0) -= 1.0 / 3.0*(I1 - Ratio);
			Sigma(1, 1) -= 1.0 / 3.0*(I1 - Ratio);
			Sigma(2, 2) -= 1.0 / 3.0*(I1 - Ratio);
			I1 = Ratio;
		}

		// Shear stress based on the elastic assumption (S_e n+1)
		ShearStress = Sigma - 1.0 / 3.0* I1 *OrthoSys::I;
		J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
			2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
			2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));


		// Check the elastic prediction step by the failure criteria
		if ((sqrt(J2) + alpha*I1 - kf) > 0.0)
		{
			// Shear stress based on the existing stress (S n)
			ShearStress = Stress - 1.0 / 3.0*(Stress(0, 0) + Stress(1, 1) + Stress(2, 2))*OrthoSys::I;
			J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
				2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
				2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));

			if (sqrt(J2) > 0.0)
			{
				Mat3_t temp, Plastic;
				double sum, dLanda;

				// calculating the plastic term based on the existing shear stress and strain rate
				temp = abab(ShearStress, StrainRate);
				sum = temp(0, 0) + temp(0, 1) + temp(0, 2) + temp(1, 0) + temp(1, 1) + temp(1, 2) + temp(2, 0) + temp(2, 1) + temp(2, 2);
				switch (Fail)
				{
				case 2:
					dLanda = 1.0 / (9.0*alpha*alpha*K + G)*((3.0*alpha*K*I1strain) + (G / sqrt(J2))*sum);
					Plastic = 3.0*alpha*K*I + G / sqrt(J2)*ShearStress;
					break;
				case 3:
					dLanda = 1.0 / (9.0*alpha*K*3.0*sin(psi) + G)*((3.0*alpha*K*I1strain) + (G / sqrt(J2))*sum);
					Plastic = 3.0*3.0*sin(psi)*K*I + G / sqrt(J2)*ShearStress;
					break;
				default:
					std::cout << "Failure Type No is out of range. Please correct it and run again" << std::endl;
					std::cout << "2 => Associated flow rule" << std::endl;
					std::cout << "3 => non-associated flow rule" << std::endl;
					abort();
					break;
				}
				// Apply the plastic term
				if (ct == 30)
					Sigma = Sigma - dt*(dLanda*Plastic);
				else
					Sigma = Sigma - 2.0*dt*(dLanda*Plastic);
			}

			//Scale back
			I1 = Sigma(0, 0) + Sigma(1, 1) + Sigma(2, 2);
			if ((kf - alpha*I1) < 0.0)
			{
				double Ratio;
				if (alpha == 0.0) Ratio = 0.0; else Ratio = kf / alpha;
				Sigma(0, 0) -= 1.0 / 3.0*(I1 - Ratio);
				Sigma(1, 1) -= 1.0 / 3.0*(I1 - Ratio);
				Sigma(2, 2) -= 1.0 / 3.0*(I1 - Ratio);
				I1 = Ratio;
			}
			ShearStress = Sigma - 1.0 / 3.0* I1 *OrthoSys::I;
			J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
				2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
				2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));

			if ((sqrt(J2) + alpha*I1 - kf) > 0.0 && sqrt(J2) > 0.0) Sigma = I1 / 3.0*OrthoSys::I + (kf - alpha*I1) / sqrt(J2) * ShearStress;
		}
	}

	Stress = Strain;
	if (ct == 30)
		Strain = dt*StrainRate + Strain;
	else
		Strain = 2.0*dt*StrainRate + Strainb;
	Strainb = Stress;

	if (VarPorosity)
	{
		if (IsFree)
		{
			double ev = (Strain(0, 0) + Strain(1, 1) + Strain(2, 2));
			n = (n0 + ev) / (1.0 + ev);
			switch (SeepageType)
			{
			case 0:
				break;
			case 1:
				k = n*n*n*d*d / (180.0*(1.0 - n)*(1.0 - n));
				break;
			case 2:
				k = n*n*n*d*d / (150.0*(1.0 - n)*(1.0 - n));
				k2 = 1.75*(1.0 - n) / (n*n*n*d);
				break;
			case 3:
				k = n*n*n*d*d / (150.0*(1.0 - n)*(1.0 - n));
				k2 = 0.4 / (n*n*d);
				break;
			default:
				std::cout << "Seepage Type No is out of range. Please correct it and run again" << std::endl;
				std::cout << "0 => Darcy's Law" << std::endl;
				std::cout << "1 => Darcy's Law & Kozeny–Carman Eq" << std::endl;
				std::cout << "2 => The Forchheimer Eq & Ergun Coeffs" << std::endl;
				std::cout << "3 => The Forchheimer Eq & Den Adel Coeffs" << std::endl;
				abort();
				break;
			}
		}
		else
			n = n0;
	}
	else
		n = n0;


}

void Particle::ScalebackMat3(size_t Dimension, size_t Scheme)
{
	double I1, J2, alpha, kf;

	if (Dimension == 0.0)
	{
		// Drucker-Prager failure criterion for plane strain
		alpha = tan(phi) / sqrt(9.0 + 12.0*tan(phi)*tan(phi));
		kf = 3.0 * c / sqrt(9.0 + 12.0*tan(phi)*tan(phi));
	}
	else
	{
		// Drucker-Prager failure criterion for 3D
		alpha = (2.0*  sin(phi)) / (sqrt(3.0)*(3.0 - sin(phi)));
		kf = (6.0*c*cos(phi)) / (sqrt(3.0)*(3.0 - sin(phi)));
	}

	// Bring back stressb to the apex of the failure criteria
	I1 = Sigmab(0, 0) + Sigmab(1, 1) + Sigmab(2, 2);
	if ((kf - alpha*I1) < 0.0)
	{
		double Ratio;
		if (alpha == 0.0) Ratio = 0.0; else Ratio = kf / alpha;
		Sigmab(0, 0) -= 1.0 / 3.0*(I1 - Ratio);
		Sigmab(1, 1) -= 1.0 / 3.0*(I1 - Ratio);
		Sigmab(2, 2) -= 1.0 / 3.0*(I1 - Ratio);
		I1 = Ratio;
	}

	ShearStress = Sigmab - 1.0 / 3.0* I1 *OrthoSys::I;
	J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
		2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
		2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));

	if ((sqrt(J2) + alpha*I1 - kf) > 0.0 && sqrt(J2) > 0.0) Sigmab = I1 / 3.0*OrthoSys::I + (kf - alpha*I1) / sqrt(J2) * ShearStress;

	// Bring back stress to the apex of the failure criteria
	if (Scheme == 0)
	{
		I1 = Sigma(0, 0) + Sigma(1, 1) + Sigma(2, 2);
		if ((kf - alpha*I1) < 0.0)
		{
			double Ratio;
			if (alpha == 0.0) Ratio = 0.0; else Ratio = kf / alpha;
			Sigma(0, 0) -= 1.0 / 3.0*(I1 - Ratio);
			Sigma(1, 1) -= 1.0 / 3.0*(I1 - Ratio);
			Sigma(2, 2) -= 1.0 / 3.0*(I1 - Ratio);
			I1 = Ratio;
		}

		ShearStress = Sigma - 1.0 / 3.0* I1 *OrthoSys::I;
		J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
			2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
			2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));

		if ((sqrt(J2) + alpha*I1 - kf) > 0.0 && sqrt(J2) > 0.0) Sigma = I1 / 3.0*OrthoSys::I + (kf - alpha*I1) / sqrt(J2) * ShearStress;
	}
	else
	{
		I1 = Sigmaa(0, 0) + Sigmaa(1, 1) + Sigmaa(2, 2);
		if ((kf - alpha*I1) < 0.0)
		{
			double Ratio;
			if (alpha == 0.0) Ratio = 0.0; else Ratio = kf / alpha;
			Sigmaa(0, 0) -= 1.0 / 3.0*(I1 - Ratio);
			Sigmaa(1, 1) -= 1.0 / 3.0*(I1 - Ratio);
			Sigmaa(2, 2) -= 1.0 / 3.0*(I1 - Ratio);
			I1 = Ratio;
		}
		ShearStress = Sigmaa - 1.0 / 3.0* I1 *OrthoSys::I;
		J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
			2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
			2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));
		if ((sqrt(J2) + alpha*I1 - kf) > 0.0 && sqrt(J2) > 0.0) Sigmaa = I1 / 3.0*OrthoSys::I + (kf - alpha*I1) / sqrt(J2) * ShearStress;
	}
}

void Particle::Move_Leapfrog(Mat3_t I, double dt)
{
	if (FirstStep)
	{
		Densitya = Density - dt / 2.0*dDensity;
		va = v - dt / 2.0*a;
	}
	Densityb = Densitya;
	Densitya += dt*dDensity;
	Density = (Densitya + Densityb) / 2.0;
	vb = va;
	va += dt*a;
	v = (va + vb) / 2.0;
	x += dt*va;

	switch (Material)
	{
	case 1:
		Mat1(dt);
		break;
	case 2:
		Mat2Leapfrog(dt);
		break;
	case 3:
		Mat3Leapfrog(I, dt);
		break;
	default:
		std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "1 => Fluid" << std::endl;
		std::cout << "2 => Solid" << std::endl;
		std::cout << "3 => Soil" << std::endl;
		abort();
		break;
	}
	if (FirstStep) FirstStep = false;

}

void Particle::Mat2Leapfrog(double dt)
{
	Pressure = EOS(PresEq, Cs, P0, Density, RefDensity);

	// Jaumann rate terms
	Mat3_t RotationRateT, SRT, RS;
	MatVec::Trans(RotationRate, RotationRateT);
	MatVec::Mult(ShearStress, RotationRateT, SRT);
	MatVec::Mult(RotationRate, ShearStress, RS);

	// Elastic prediction step (ShearStress_e n+1)
	if (FirstStep)
		ShearStressa = -dt / 2.0*(2.0*G*(StrainRate - 1.0 / 3.0*(StrainRate(0, 0) + StrainRate(1, 1) + StrainRate(2, 2))*OrthoSys::I) + SRT + RS) + ShearStress;

	ShearStressb = ShearStressa;
	ShearStressa = dt*(2.0*G*(StrainRate - 1.0 / 3.0*(StrainRate(0, 0) + StrainRate(1, 1) + StrainRate(2, 2))*OrthoSys::I) + SRT + RS) + ShearStressa;

	if (Fail == 1)
	{
		double J2 = 0.5*(ShearStressa(0, 0)*ShearStressa(0, 0) + 2.0*ShearStressa(0, 1)*ShearStressa(1, 0) +
			2.0*ShearStressa(0, 2)*ShearStressa(2, 0) + ShearStressa(1, 1)*ShearStressa(1, 1) +
			2.0*ShearStressa(1, 2)*ShearStressa(2, 1) + ShearStressa(2, 2)*ShearStressa(2, 2));
		//Scale back
		ShearStressa = std::min((Sigmay / sqrt(3.0*J2)), 1.0)*ShearStressa;
	}
	ShearStress = 1.0 / 2.0*(ShearStressa + ShearStressb);

	Sigma = -Pressure * OrthoSys::I + ShearStress;

	if (FirstStep)
		Straina = -dt / 2.0*StrainRate + Strain;
	Strainb = Straina;
	Straina = dt*StrainRate + Straina;
	Strain = 1.0 / 2.0*(Straina + Strainb);


	if (Fail > 1)
	{
		std::cout << "Undefined failure criteria for solids" << std::endl;
		abort();
	}
}

void Particle::Mat3Leapfrog(Mat3_t I, double dt)
{
	Mat3_t RotationRateT, Stress, SRT, RS;
	double I1, J2, alpha, kf, I1strain;

	// Jaumann rate terms
	MatVec::Trans(RotationRate, RotationRateT);
	MatVec::Mult(Sigma, RotationRateT, SRT);
	MatVec::Mult(RotationRate, Sigma, RS);

	// Volumetric strain
	I1strain = StrainRate(0, 0) + StrainRate(1, 1) + StrainRate(2, 2);

	// Elastic prediction step (Sigma_e n+1)
	if (FirstStep)
		Sigmaa = -dt / 2.0*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate - 1.0 / 3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigma;

	Sigmab = Sigmaa;
	Sigmaa = dt*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate - 1.0 / 3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigmaa;

	if (Fail > 1)
	{
		if (I(2, 2) == 0.0)
		{
			// Drucker-Prager failure criterion for plane strain
			alpha = tan(phi) / sqrt(9.0 + 12.0*tan(phi)*tan(phi));
			kf = 3.0 * c / sqrt(9.0 + 12.0*tan(phi)*tan(phi));
		}
		else
		{
			// Drucker-Prager failure criterion for 3D
			alpha = (2.0*  sin(phi)) / (sqrt(3.0)*(3.0 - sin(phi)));
			kf = (6.0*c*cos(phi)) / (sqrt(3.0)*(3.0 - sin(phi)));
		}

		// Bring back stress to the apex of the failure criteria
		I1 = Sigmaa(0, 0) + Sigmaa(1, 1) + Sigmaa(2, 2);
		if ((kf - alpha*I1) < 0.0)
		{
			double Ratio;
			if (alpha == 0.0) Ratio = 0.0; else Ratio = kf / alpha;
			Sigmaa(0, 0) -= 1.0 / 3.0*(I1 - Ratio);
			Sigmaa(1, 1) -= 1.0 / 3.0*(I1 - Ratio);
			Sigmaa(2, 2) -= 1.0 / 3.0*(I1 - Ratio);
			I1 = Ratio;
		}

		// Shear stress based on the elastic assumption (S_e n+1)
		ShearStress = Sigmaa - 1.0 / 3.0* I1 *OrthoSys::I;
		J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
			2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
			2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));


		// Check the elastic prediction step by the failure criteria
		if ((sqrt(J2) + alpha*I1 - kf) > 0.0)
		{
			// Shear stress based on the existing stress (S n)
			ShearStress = Sigma - 1.0 / 3.0*(Sigma(0, 0) + Sigma(1, 1) + Sigma(2, 2))*OrthoSys::I;
			J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
				2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
				2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));

			if (sqrt(J2) > 0.0)
			{
				Mat3_t temp, Plastic;
				double sum, dLanda;

				// calculating the plastic term based on the existing shear stress and strain rate
				temp = abab(ShearStress, StrainRate);
				sum = temp(0, 0) + temp(0, 1) + temp(0, 2) + temp(1, 0) + temp(1, 1) + temp(1, 2) + temp(2, 0) + temp(2, 1) + temp(2, 2);
				switch (Fail)
				{
				case 2:
					dLanda = 1.0 / (9.0*alpha*alpha*K + G)*((3.0*alpha*K*I1strain) + (G / sqrt(J2))*sum);
					Plastic = 3.0*alpha*K*I + G / sqrt(J2)*ShearStress;
					break;
				case 3:
					dLanda = 1.0 / (9.0*alpha*K*3.0*sin(psi) + G)*((3.0*alpha*K*I1strain) + (G / sqrt(J2))*sum);
					Plastic = 3.0*3.0*sin(psi)*K*I + G / sqrt(J2)*ShearStress;
					break;
				default:
					std::cout << "Failure Type No is out of range. Please correct it and run again" << std::endl;
					std::cout << "2 => Associated flow rule" << std::endl;
					std::cout << "3 => non-associated flow rule" << std::endl;
					abort();
					break;
				}
				Sigmaa = Sigmaa - dt*(dLanda*Plastic);
			}

			I1 = Sigmaa(0, 0) + Sigmaa(1, 1) + Sigmaa(2, 2);
			if ((kf - alpha*I1) < 0.0)
			{
				double Ratio;
				if (alpha == 0.0) Ratio = 0.0; else Ratio = kf / alpha;
				Sigmaa(0, 0) -= 1.0 / 3.0*(I1 - Ratio);
				Sigmaa(1, 1) -= 1.0 / 3.0*(I1 - Ratio);
				Sigmaa(2, 2) -= 1.0 / 3.0*(I1 - Ratio);
				I1 = Ratio;
			}
			ShearStress = Sigmaa - 1.0 / 3.0* I1 *OrthoSys::I;
			J2 = 0.5*(ShearStress(0, 0)*ShearStress(0, 0) + 2.0*ShearStress(0, 1)*ShearStress(1, 0) +
				2.0*ShearStress(0, 2)*ShearStress(2, 0) + ShearStress(1, 1)*ShearStress(1, 1) +
				2.0*ShearStress(1, 2)*ShearStress(2, 1) + ShearStress(2, 2)*ShearStress(2, 2));
			if ((sqrt(J2) + alpha*I1 - kf) > 0.0 && sqrt(J2) > 0.0) Sigmaa = I1 / 3.0*OrthoSys::I + (kf - alpha*I1) / sqrt(J2) * ShearStress;
		}
	}
	Sigma = 1.0 / 2.0*(Sigmaa + Sigmab);

	if (FirstStep)
		Straina = -dt / 2.0*StrainRate + Strain;
	Strainb = Straina;
	Straina = dt*StrainRate + Straina;
	Strain = 1.0 / 2.0*(Straina + Strainb);

	if (VarPorosity)
	{
		if (IsFree)
		{
			double ev = (Strain(0, 0) + Strain(1, 1) + Strain(2, 2));
			n = (n0 + ev) / (1.0 + ev);
			switch (SeepageType)
			{
			case 0:
				break;
			case 1:
				k = n*n*n*d*d / (180.0*(1.0 - n)*(1.0 - n));
				break;
			case 2:
				k = n*n*n*d*d / (150.0*(1.0 - n)*(1.0 - n));
				k2 = 1.75*(1.0 - n) / (n*n*n*d);
				break;
			case 3:
				k = n*n*n*d*d / (150.0*(1.0 - n)*(1.0 - n));
				k2 = 0.4 / (n*n*d);
				break;
			default:
				std::cout << "Seepage Type No is out of range. Please correct it and run again" << std::endl;
				std::cout << "0 => Darcy's Law" << std::endl;
				std::cout << "1 => Darcy's Law & Kozeny–Carman Eq" << std::endl;
				std::cout << "2 => The Forchheimer Eq & Ergun Coeffs" << std::endl;
				std::cout << "3 => The Forchheimer Eq & Den Adel Coeffs" << std::endl;
				abort();
				break;
			}
		}
		else
			n = n0;
	}
	else
		n = n0;


}

void Particle::translate(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin)
{
	x = x + dt*v + 0.5*dt*dt*a;

	// Evolve velocity
	Vec3_t temp;
	temp = v;
	v = vb + 2 * dt*a;
	vb = temp;

	//Periodic BC particle position update
	if (Domainsize(0) > 0.0)
	{
		(x(0) > (domainmax(0))) ? x(0) -= Domainsize(0) : x(0);
		(x(0) < (domainmin(0))) ? x(0) += Domainsize(0) : x(0);
	}
	if (Domainsize(1) > 0.0)
	{
		(x(1) > (domainmax(1))) ? x(1) -= Domainsize(1) : x(1);
		(x(1) < (domainmin(1))) ? x(1) += Domainsize(1) : x(1);
	}
	if (Domainsize(2) > 0.0)
	{
		(x(2) > (domainmax(2))) ? x(2) -= Domainsize(2) : x(2);
		(x(2) < (domainmin(2))) ? x(2) += Domainsize(2) : x(2);
	}
}
