// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#undef max
#undef min
#include <fstream>
#include <time.h>

// Globals:
MACGrid target;


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 


#define FOR_EACH_FACE_X \
   for(int k = 0; k < theDim[MACGrid::Z]; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_FACE_Y \
   for(int k = 0; k < theDim[MACGrid::Z]; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE_Z \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 


MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);

   setUpAMatrix();
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    // TODO: Set initial values for density, temperature, and velocity.
	//std::cout<<mU(1,1,1)<<std::endl;
	mD(0,1,0) = 1.0;
	mT(0,1,0) = 30.0;
	//mU(0,1,0) = 0.0;
	mV(0,1,0) = 10.0;
	//mV(1,0,1) = 5;
	//mW(1,1,1) = 0;
}

void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target.
	target.mU = mU;
    target.mV = mV;
    target.mW = mW;

	FOR_EACH_FACE_X
	{
		double Xg, Yg, Zg;

		if(i< theDim[MACGrid::X])
		{
			Xg = getCenter(i,j,k).n[0] - 0.5*theCellSize;
			Yg = getCenter(i,j,k).n[1];
			Zg = getCenter(i,j,k).n[2];
		}
		else
		{
			Xg = getCenter(i-1,j,k).n[0] + 0.5*theCellSize;
			Yg = getCenter(i-1,j,k).n[1];
			Zg = getCenter(i-1,j,k).n[2];
		}
		vec3 posG = vec3(Xg,Yg,Zg);

		
		vec3 vel = vec3(target.mU(i,j,k),
						0.25 *(target.mV(i,j,k) + target.mV(i-1,j,k) + target.mV(i,j+1,k) + target.mV(i-1,j+1,k)),
						0.25 *(target.mW(i,j,k) + target.mW(i-1,j,k) + target.mW(i,j,k+1) + target.mW(i-1,j,k+1)));
		
		//vec3 vel = target.getVelocity(posG);
		vec3 posP = posG - dt * vel;
		target.mU(i,j,k) = mU.interpolate(posP);
	}

	FOR_EACH_FACE_Y
	{
		double Xg, Yg, Zg;

		if(j< theDim[MACGrid::Y])
		{
			Xg = getCenter(i,j,k).n[0];
			Yg = getCenter(i,j,k).n[1] - 0.5*theCellSize;
			Zg = getCenter(i,j,k).n[2];
		}
		else
		{
			Xg = getCenter(i,j-1,k).n[0];
			Yg = getCenter(i,j-1,k).n[1] + 0.5*theCellSize;
			Zg = getCenter(i,j-1,k).n[2];
		}
		vec3 posG = vec3(Xg,Yg,Zg);
		//vec3 vel = target.getVelocity(posG);
		
		vec3 vel = vec3(0.25 *(target.mU(i,j,k) + target.mU(i,j-1,k) + target.mU(i+1,j,k) + target.mU(i+1,j-1,k)),
						target.mV(i,j,k),
						0.25 *(target.mW(i,j,k) + target.mW(i,j-1,k) + target.mW(i,j,k+1) + target.mW(i,j-1,k+1)));
		
		vec3 posP = posG - dt * vel;
		//std::cout<<i<<j<<k<<' '<<posP<<std::endl;
		target.mV(i,j,k) = mV.interpolate(posP);
	}

	FOR_EACH_FACE_Z
	{
		double Xg, Yg, Zg;

		if(k< theDim[MACGrid::Z])
		{
			Xg = getCenter(i,j,k).n[0];
			Yg = getCenter(i,j,k).n[1];
			Zg = getCenter(i,j,k).n[2] - 0.5*theCellSize;
		}
		else
		{
			Xg = getCenter(i,j,k-1).n[0];
			Yg = getCenter(i,j,k-1).n[1];
			Zg = getCenter(i,j,k-1).n[2] + 0.5*theCellSize;
		}

		vec3 posG = vec3(Xg,Yg,Zg);
		//vec3 vel = target.getVelocity(posG);
		
		vec3 vel = vec3(0.25 *(target.mU(i,j,k) + target.mU(i,j,k-1) + target.mU(i+1,j,k) + target.mU(i+1,j,k-1)),
						0.25 *(target.mV(i,j,k) + target.mV(i,j,k-1) + target.mV(i,j+1,k) + target.mV(i,j+1,k-1)),
						mW(i,j,k));
		
		vec3 posP = posG - dt * vel;
		target.mW(i,j,k) = mW.interpolate(posP);
	}
	
    // Then save the result to our object.
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target.
	target.mT = mT;
	FOR_EACH_CELL
	{
		vec3 posG = target.getCenter(i,j,k);
		vec3 vel = target.getVelocity(posG);

		vec3 posP = posG - dt*vel;

		target.mT(i,j,k) = target.mT.interpolate(posP);
	}

    // Then save the result to our object.
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
	
    // TODO: Calculate new densitities and store in target.
	target.mD = mD;
	FOR_EACH_CELL
	{
		vec3 posG = getCenter(i,j,k);
		vec3 vel = target.getVelocity(posG);

		vec3 posP = posG - dt*vel;

		target.mD(i,j,k) = target.mD.interpolate(posP);
	}

    // Then save the result to our object.
    mD = target.mD;
}

void MACGrid::computeBouyancy(double dt)
{
	// TODO: Calculate bouyancy and store in target.
		target.mV = mV;
	
	double alpha = 0.1;
	double beta = 0.1;
	double T0 = 25.0;
	double T,S;
	double fbou;
	double fg = -9.8;

	FOR_EACH_FACE_Y
	{
		if(j == 0)
		{
			T = target.mT(i,j,k);
			S = target.mD(i,j,k);
			continue;
		}
		else if(j == theDim[MACGrid::Y])
		{
			T = target.mT(i,j-1,k);
			S = target.mD(i,j-1,k);
			continue;
		}
		else
		{
		T = (target.mT(i,j,k)+ target.mT(i,j-1,k))/2;
		S = (target.mD(i,j,k)+ target.mD(i,j-1,k))/2;
		}
		fbou = -alpha * S + beta * T;
		target.mV(i,j,k) = target.mV(i,j,k) + (fbou+fg)*dt;
	}

   // Then save the result to our object.
   mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	
	// Apply the forces to the current velocity and store the result in target.
	double epsilon = 1.0;
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	
	GridData U;
	GridData V;
	GridData W;
	
	U.initialize();
	V.initialize();
	W.initialize();

	GridData omega1;
	GridData omega2;
	GridData omega3;
	
	omega1.initialize();
	omega2.initialize();
	omega3.initialize();

	GridData grad1;
	GridData grad2;
	GridData grad3;
	//GridData grad;

	grad1.initialize();
	grad2.initialize();
	grad3.initialize();
	//grad.initialize();

	FOR_EACH_CELL
	{
		U(i,j,k) = 0.5*(mU(i,j,k)+mU(i+1,j,k));
		V(i,j,k) = 0.5*(mV(i,j,k)+mV(i,j+1,k));
		W(i,j,k) = 0.5*(mW(i,j,k)+mW(i,j,k+1));
		//std::cout<<U(i,j,k)<<' '<<V(i,j,k)<<' '<<W(i,j,k);
	}

	//Get omega
	double p[3][4];
	FOR_EACH_CELL
	{
		if(i == 0)
		{
			p[1][3] = 0;
			p[2][1] = 0;
		}
		else if(i == theDim[MACGrid::X]-1)
		{
			p[1][2] = 0;
			p[2][0] = 0;
		}
		else
		{
			p[1][3] = W(i-1,j,k);
			p[2][1] = V(i-1,j,k);
			p[1][2] = W(i+1,j,k);
			p[2][0] = V(i+1,j,k);
		}
		
		if(j == 0)
		{
			p[0][1] = 0;
			p[2][3] = 0;
		}
		else if(j == theDim[MACGrid::Y]-1)
		{
			p[0][0] = 0;
			p[2][2] = 0;
		}
		else
		{
			p[0][1] = W(i,j-1,k);
			p[2][3] = U(i,j-1,k);
			p[0][0] = W(i,j+1,k);
			p[2][2] = U(i,j+1,k);
		}

		if(k == 0)
		{
			p[0][3] = 0;
			p[1][1] = 0;
		}
		else if(k == theDim[MACGrid::Z]-1)
		{
			p[0][2] = 0;
			p[1][0] = 0;
		}
		else
		{
			p[0][3] = V(i,j,k-1);
			p[1][1] = U(i,j,k-1);
			p[0][2] = V(i,j,k+1);
			p[1][0] = U(i,j,k+1);
		}
		omega1(i,j,k) = 0.5/theCellSize*(p[0][0] - p[0][1] -p[0][2] + p[0][3]);
		omega2(i,j,k) = 0.5/theCellSize*(p[1][0] - p[1][1] -p[1][2] + p[1][3]);
		omega3(i,j,k) = 0.5/theCellSize*(p[2][0] - p[2][1] -p[2][2] + p[2][3]);
		//std::cout<<omega1(i,j,k)<<' '<<omega2(i,j,k)<<' '<<omega3(i,j,k)<<std::endl;
	}

	//Get the gradient of omega
	double q[3][2];

	FOR_EACH_CELL
	{
		if(i == 0)
		{
			q[0][1] = 0;
			//q[0][1] = 0.25*sqrt(omega1(i,j,k)*omega1(i,j,k)+omega2(i,j,k)*omega2(i,j,k)+omega3(i,j,k)*omega3(i,j,k));
		}
		else if(i == theDim[MACGrid::X]-1)
		{
			q[0][0] = 0;
			//q[0][0] = 0.25*sqrt(omega1(i,j,k)*omega1(i,j,k)+omega2(i,j,k)*omega2(i,j,k)+omega3(i,j,k)*omega3(i+1,j,k));
		}
		else
		{
			q[0][1] = sqrt(omega1(i-1,j,k)*omega1(i-1,j,k)+omega2(i-1,j,k)*omega2(i-1,j,k)+omega3(i-1,j,k)*omega3(i-1,j,k));
			q[0][0] = sqrt(omega1(i+1,j,k)*omega1(i+1,j,k)+omega2(i+1,j,k)*omega2(i+1,j,k)+omega3(i+1,j,k)*omega3(i+1,j,k));
		}
		
		if(j == 0)
		{
			//q[1][1] =  0.25*sqrt(omega1(i,j,k)*omega1(i,j,k)+omega2(i,j,k)*omega2(i,j,k)+omega3(i,j,k)*omega3(i,j,k));
			q[1][1] =  0;
		}
		else if(j == theDim[MACGrid::Y]-1)
		{
			//q[1][0] =  0.25*sqrt(omega1(i,j,k)*omega1(i,j,k)+omega2(i,j,k)*omega2(i,j,k)+omega3(i,j,k)*omega3(i,j,k));
			q[1][0] =  0;
		}
		else
		{
			q[1][1] = sqrt(omega1(i,j-1,k)*omega1(i,j-1,k)+omega2(i,j-1,k)*omega2(i,j-1,k)+omega3(i,j-1,k)*omega3(i,j-1,k));
			q[1][0] = sqrt(omega1(i,j+1,k)*omega1(i,j+1,k)+omega2(i,j+1,k)*omega2(i,j+1,k)+omega3(i,j+1,k)*omega3(i,j+1,k));
		}

		if(k == 0)
		{
			//q[2][1] =  0.25*sqrt(omega1(i,j,k)*omega1(i,j,k)+omega2(i,j,k)*omega2(i,j,k)+omega3(i,j,k)*omega3(i,j,k));
			q[2][1] =  0;
		}
		else if(k == theDim[MACGrid::Z]-1)
		{
			//q[2][0] =  0.25*sqrt(omega1(i,j,k)*omega1(i,j,k)+omega2(i,j,k)*omega2(i,j,k)+omega3(i,j,k)*omega3(i,j,k));
			q[2][0] =  0;
		}
		else
		{
			q[2][1] = sqrt(omega1(i,j,k-1)*omega1(i,j,k-1)+omega2(i,j,k-1)*omega2(i,j,k-1)+omega3(i,j,k-1)*omega3(i,j,k-1));
			q[1][0] = sqrt(omega1(i,j,k+1)*omega1(i,j,k+1)+omega2(i,j,k+1)*omega2(i,j,k+1)+omega3(i,j,k+1)*omega3(i,j,k+1));
		}
		grad1(i,j,k) = 0.5/theCellSize*(q[0][0] - q[0][1]);
		grad2(i,j,k) = 0.5/theCellSize*(q[1][0] - q[1][1]);
		grad3(i,j,k) = 0.5/theCellSize*(q[2][0] - q[2][1]);
		//std::cout<<grad1(i,j,k)<<' '<<grad2(i,j,k)<<' '<<grad3(i,j,k)<<std::endl;
	}

	FOR_EACH_CELL
	{
		double norm = sqrt(grad1(i,j,k)*grad1(i,j,k) + grad2(i,j,k)*grad2(i,j,k) + grad3(i,j,k)*grad3(i,j,k)) + 10e-20;
		//std::cout<<norm<<std::endl;
		
		double norm1 = grad1(i,j,k)/norm;
		double norm2 = grad2(i,j,k)/norm;
		double norm3 = grad3(i,j,k)/norm;
		//std::cout<<norm1<<' '<<norm2<<' '<<norm3<<std::endl;
		double fconf1 = epsilon * theCellSize * (omega3(i,j,k)*norm2 - omega2(i,j,k)*norm3);
		double fconf2 = epsilon * theCellSize * (omega1(i,j,k)*norm3 - omega3(i,j,k)*norm1);
		double fconf3 = epsilon * theCellSize * (omega2(i,j,k)*norm1 - omega1(i,j,k)*norm2);
		if(norm<10e-18)
			continue;
		target.mU(i,j,k) = target.mU(i,j,k) + fconf1*dt;
		target.mV(i,j,k) = target.mV(i,j,k) + fconf2*dt;
		target.mW(i,j,k) = target.mW(i,j,k) + fconf3*dt;
	}
	
	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
}

void MACGrid::project(double dt)
{
	// TODO: Solve Ap = d for pressure.
	target.mP = mP;
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	
	//GridData p;
	GridData d;
	double Vsolid = 0.0;

	//p.initialize();
	d.initialize();

	FOR_EACH_CELL 
	{
		if(i==0)
			target.mU(i,j,k) = 0;
		if(i == theDim[MACGrid::X]-1)
			target.mU(i+1,j,k) = 0;
		if(j==0)
			target.mV(i,j,k) = 0;
		if(j == theDim[MACGrid::Y]-1)
			target.mV(i,j+1,k) = 0;
		if(k==0)
			target.mW(i,j,k) = 0;
		if(k == theDim[MACGrid::Z]-1)
			target.mW(i,j,k+1) = 0;

		d(i,j,k) = -theCellSize/dt *(target.mU(i+1,j,k)-target.mU(i,j,k)+target.mV(i,j+1,k)-target.mV(i,j,k)+target.mW(i,j,k+1)-target.mW(i,j,k));

	}
	
	bool conj;
	clock_t startTime = clock();
	conj = conjugateGradient(AMatrix, target.mP, d, 1000, 0.00001);
	clock_t endTime = clock();
	std::cout<<endTime - startTime<<std::endl;
	FOR_EACH_FACE_X
	{
		if(i == 0)
		{
			target.mU(i,j,k) = 0;
			continue;
		}
		else if(i == theDim[MACGrid::X])
		{
			target.mU(i,j,k) = 0;
			continue;
		}
		target.mU(i,j,k) = target.mU(i,j,k) - dt/theCellSize *(target.mP(i,j,k) - target.mP(i-1,j,k));
	}

	FOR_EACH_FACE_Y
	{
		if(j == 0)
		{
			target.mV(i,j,k) = 0;
			continue;
		}
		else if(j == theDim[MACGrid::Y])
		{
			target.mV(i,j,k) = 0;
			continue;
		}
		target.mV(i,j,k) = target.mV(i,j,k) - dt/theCellSize *(target.mP(i,j,k) - target.mP(i,j-1,k));
	}

	FOR_EACH_FACE_Z
	{
		if(k == 0)
		{
			target.mW(i,j,k) = 0;
			continue;
		}
		else if(k == theDim[MACGrid::Z])
		{
			target.mW(i,j,k) = 0;
			continue;
		}
		target.mW(i,j,k) = target.mW(i,j,k) - dt/theCellSize *(target.mP(i,j,k) - target.mP(i,j,k-1));
	}

	// 1. Construct d
	// 2. Construct A
	// 3. Solve for p
	// Subtract pressure from our velocity and save in target.

	// Then save the result to our object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
	// IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
	
	//check the divergence
	FOR_EACH_CELL
	{
	//std::cout<<(target.mU(i+1,j,k)-target.mU(i,j,k)+target.mV(i,j+1,k)-target.mV(i,j,k)+target.mW(i,j,k+1)-target.mW(i,j,k))<<std::endl;
	}
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::setUpAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}





/////////////////////////////////////////////////////////////////////

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	// TODO: Apply a preconditioner here.
	
	GridData precon;
	precon.initialize();

	double tau = 0.97;
	FOR_EACH_CELL
	{
		double e = A.diag(i,j,k)-(A.plusI(i-1,j,k)*precon(i-1,j,k))*(A.plusI(i-1,j,k)*precon(i-1,j,k))-
								 (A.plusJ(i,j-1,k)*precon(i,j-1,k))*(A.plusJ(i,j-1,k)*precon(i,j-1,k))-
								 (A.plusK(i,j,k-1)*precon(i,j,k-1))*(A.plusK(i,j,k-1)*precon(i,j,k-1))-
								 tau*(A.plusI(i-1,j,k)*(A.plusJ(i-1,j,k)+A.plusK(i-1,j,k))*precon(i-1,j,k)*precon(i-1,j,k)+
									  A.plusJ(i,j-1,k)*(A.plusI(i,j-1,k)+A.plusK(i,j-1,k))*precon(i,j-1,k)*precon(i,j-1,k)+
									  A.plusK(i,j,k-1)*(A.plusI(i,j,k-1)+A.plusJ(i,j,k-1))*precon(i,j,k-1)*precon(i,j,k-1));
		precon(i,j,k) = 1/sqrt(e+10e-30);
	}

	
	// For now, just bypass the preconditioner:
	z = r;

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		// TODO: Apply a preconditioner here.
		
		GridData q;
		q.initialize();

		FOR_EACH_CELL
		{
			double t = r(i,j,k) - A.plusI(i-1,j,k)*precon(i-1,j,k)*q(i-1,j,k)-
								  A.plusJ(i,j-1,k)*precon(i,j-1,k)*q(i,j-1,k)-
								  A.plusK(i,j,k-1)*precon(i,j,k-1)*q(i,j,k-1);
			q(i,j,k) = t*precon(i,j,k);
		}

		FOR_EACH_CELL_REVERSE
		{
			double t = q(i,j,k) - A.plusI(i,j,k)*precon(i,j,k)*z(i+1,j,k)-
								  A.plusJ(i,j,k)*precon(i,j,k)*z(i,j+1,k)-
								  A.plusK(i,j,k)*precon(i,j,k)*z(i,j,k+1);
			z(i,j,k) = t*precon(i,j,k);
		}
		
		// For now, just bypass the preconditioner:
		//z = r;

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}




/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}




/////////////////////////////////////////////////////////////////////

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // Draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
		   // Un-comment the line below if you want all of the velocity lines to be the same length.
           //vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD(i, j, k); 
    return vec4(1.0, 1.0, 1.0, value);

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity(pt); 
    return vec4(1.0, 1.0, 1.0, value);

}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
