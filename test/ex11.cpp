#include <math.h>
#include "MG.hpp"
#include "lduAgglomeration.hpp"
#include "chebySmoother.hpp"
#include "lduGaussSeidelSmoother.hpp"
#include "matrixConversion.hpp"
#include "PBiCGStab.hpp"
#include "PCG.hpp"
#include "lduDiagPrecond.hpp"
#include "lduDICPrecond.hpp"
#include "lduDILUPrecond.hpp"
#include "readFromOpenFOAM.hpp"
#include <sstream>
#include <string.h>

//- test using data output from OpenFOAM

using namespace SMS;

#define LOCATEFILE(newName, fileName, dir) \
{ \
	std::ostringstream os; \
	os << dir << NPROCS << "/" << fileName << "_" << MYID << ".txt"; \
	strcpy(newName, os.str().c_str()); \
} \

int main(int argc, char *argv[])
{
	/* Initialize MPI */
   	MPI_Init(&argc, &argv);
   	smsMPI::init();

   	const char* dir = "../test/openfoam/cavity2.5w/cavity2.5wp";
   	char fileName[200];

   	lduMatrix lduA;
   	LOCATEFILE(fileName, "A_p", dir);
   	constructMatrixFromOpenFOAM(lduA, fileName);

   	if(PARRUN)
   	{
   		LOCATEFILE(fileName, "interfaces_p", dir);
   		constructInterfacesFromOpenFOAM(lduA, fileName);
   	}

   	label nCells = lduA.size();
   	scalarField b(nCells);

   	LOCATEFILE(fileName, "b_p", dir);
   	constructVectorFromOpenFOAM(b, fileName);

   	scalar tol = 0.0;
	scalar relTol = 1e-6;
	label  nFaces = lduA.upper().size();


	const bool useMG = true;
	const bool usePBiCGStab = false;
	scalarField x(nCells, 0.0);

	if(useMG)
	{
		scalarField weights(nFaces);
		forAll(i, nFaces)
		{
			weights[i] = mag(lduA.upper()[i]);
		}

		lduAgglomeration aggl(lduA);
		aggl.agglomerate(weights);
		PtrList<matrix::smoother> sm(aggl.size());

		forAll(i, aggl.size())
		{
			lduGaussSeidelSmoother* smLocPtr = new lduGaussSeidelSmoother;
			sm.setLevel(i, *smLocPtr);
		}

		// forAll(i, aggl.size())
		// {
		// 	chebySmoother* smLocPtr = new chebySmoother;
		// 	sm.setLevel(i, *smLocPtr);
		// }

		MGSolver MG(lduA, aggl, sm);

		MG.SET_tolerance(tol);
		MG.SET_relTol(relTol);
		MG.SET_nPreSweeps(0);

		matrix::solverPerformance solverPerf = MG.solve(x, lduA, b);
	}
	else if(usePBiCGStab)
	{
		// lduDiagPrecond precond(lduA);

		lduDICPrecond precond(lduA);

		// lduDILUPrecond precond(lduA);

		PBiCGStab PBiCGStabSolver(precond);

		// PBiCGStabSolver.SET_minIter(5);

		// PBiCGStabSolver.SET_maxIter(50);

		matrix::solverPerformance solverPerf = PBiCGStabSolver.solve(x, lduA, b);
	}
	else
	{
		// lduDiagPrecond precond(lduA);

		lduDICPrecond precond(lduA);

		PCG PCGSolver(precond);

		// PCGSolver.SET_minIter(5);

		// PCGSolver.SET_maxIter(50);

		matrix::solverPerformance solverPerf = PCGSolver.solve(x, lduA, b);

		COUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
		COUT << "Let me check now: " << ENDL;
	}
}
