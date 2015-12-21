/*********************************************
NAME: Henry Tran
CLASS: Chemistry 6510 -- Autumn 2015
INSTRUCTOR: M. Sotomayor
DESCRIPTION: Homework #2.6 -- Numerov's Method
SPECS: Program was compiled using Microsoft
Visual Studio 2013 Community Version in a .NET
framework on a system using an Intel Xeon E5-
1360 @ 3.7Ghz (4 Core) and 16 GB of RAM 
running 64-bit Windows 7 Professional Edition
*********************************************/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <cstdio>
#include <math.h>
#include <cmath>
//#include <tgmath.h>
#include <fstream>
#include <string>
//#include <vector>
#include <sstream>
#include <ctime>
#include <locale>
#include <Windows.h>

using namespace std;

double SquareWell(double x, double OutsidePotential, double InsidePotential, double LeftBound, double RightBound) // Potential for a square well or barrier.
{
	return x > LeftBound && x < RightBound ? InsidePotential : OutsidePotential;
}

double QuarticPotential(double x, double c)
{
	return c * x * x * x * x;
}

double LinearPotential(double x, double b)
{
	return b * x;
}

double BarrierWell(double x, double L, double V0)
{
	return fabs(x) < L / 4 ? V0 : 0;
}
double HOPlus(double x, double l, double k, double Mass)
{
	return 0.5 * k * x * x + l * (l + 1) / (2 * Mass * x * x);
}
double RadialG(double x, double l, double k, double Mass, double E)
{
	return Mass * (2 * E - k * x * x) - l * (l + 1) / (x * x);
}

double RadialEF(double MinX, double MaxX, double Energy, double Mass, double l, double k, int Parity, int Index, int NumberOfSteps)
{

	double EIncrement =  Energy / 100;

	double StepSize = (MaxX - MinX) / (double)NumberOfSteps; // h in HW statement

	/* Set up output file name */
	string IStr = to_string(Index);
	string OutputName = "RadialEF(Rr)_l=" + to_string(l) + "_" + IStr + ".txt"; // Name of output for n'th eigenfunction
	string CorrectOutputName = "RadialEF(R)_" + to_string(l) + "_" + IStr + ".txt";
	ofstream PrintEF(OutputName); // Output file which will store plot data for eigenvectors.
	ofstream PrintCorrectEF(CorrectOutputName);

	/* Y is the array that stores eigenfunction y values. I don't need to store all these values, but it makes normalization easier. */
	double* Y = new double[NumberOfSteps + 1];
	//Y[0] = 1; // Initializes the loop.

	double DerivStep1 = 1; // These will hold values to calculate the sign of the derivative. In this case, we only need to check their relative signs.
	double DerivStep2 = 1;

	int DerivStepCount = 0; // Tells us how many steps we have taken.

	while (true)
	{
		Energy += EIncrement; // Take energy step
		Y[NumberOfSteps] = 0; // Necessary initial condition.
		Y[NumberOfSteps - 1] = 0.1;

		DerivStepCount++;
		double X = MaxX - StepSize - StepSize; // Start from the second in the sequence.

		for (int i = NumberOfSteps - 2; i >= 0; i--)
		{
			Y[i] = ((2 - 5 * StepSize * StepSize * RadialG(X + StepSize, l, k, Mass, Energy) / 6) * Y[i + 1] - (1 + StepSize * StepSize * RadialG(X + StepSize + StepSize, l, k, Mass, Energy) / 12) * Y[i + 2]) / (1 + StepSize * StepSize * RadialG(X, l, k, Mass, Energy) / 12);
			X -= StepSize; // Increment X (down)
		}

		DerivStep1 = DerivStep2; // Shift up derivative variables.
		DerivStep2 = Y[0];
		cout << Energy << "\t" << Y[0] << endl;
		//Sleep(200);

		if (DerivStep1 * DerivStep2 < 0 && DerivStepCount > 1) // Means sign change (i.e. we missed the boundary condition at 0) and that we have taken two steps.
		{
			EIncrement /= -10; // Reverse the energy step and make it finer.
		}
		if (EIncrement < 10E-15 || fabs(Y[0]) < 10E-15)
		{
			break;
		}
	} // End while loop

	/* This normalizes the eigenfunction */
	double IntegralPP = 0; // Holds the integral
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		IntegralPP += (Y[i] * Y[i]) * fabs(StepSize); // Approximations integral by a finite sum.
	}
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		Y[i] /= sqrt(IntegralPP); // Divide by square root of integral for normalization.
	}

	/* This prints the eigenfunction */
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		PrintEF << MinX + i * StepSize << "\t" << Y[i] << endl;
		PrintCorrectEF << MinX + i * StepSize << "\t" << Y[i] / (MinX + i * StepSize) << endl;
	}
	delete[] Y; // Clear memory allocated to Y

	return Energy; // Eigenvalue for main function to print.
}
//double RadialEF(double MinX, double MaxX, double Energy, double Mass, double l, double k, int Parity, int Index, int NumberOfSteps)
//{
//
//	double StepSize = (MaxX - MinX) / (double)NumberOfSteps;
//	double EIncrement = fabs(Energy) / 1000; // Energy step size.
//
//	double* Y = new double[NumberOfSteps + 1]; // Array that holds the eigenfunction
//
//	int EStepCount = 0; // Counts steps in the loop
//	double RelativeHeight = 10;
//
//	while (true) // While the two eigenfunctions are discontinuous.
//	{
//		EStepCount++;
//		Energy += EIncrement; // Increments the energy.
//
//		double BarrierX = pow(Energy / k + pow(Mass * (-k * l - k * l * l + Energy * Energy * Mass), 0.5) / (k * Mass), 0.5); // This is when we hit the farthest barrier
//		int Barrier = 0;
//		for (int j = 0; j < NumberOfSteps; j++)
//		{
//			if (MinX + (double)j * StepSize >= BarrierX)
//			{
//				break;
//			}
//			Barrier++;
//		}
//		double* Y_Right = new double[NumberOfSteps + 1 - Barrier + 1];
//		Y = new double[NumberOfSteps + 1];
//
//		Y_Right[0] = 1; // Initialize to begin loop;
//		Y[Barrier - 1] = 2;
//
//		Y[0] = 0; // Necessary initial condition.
//		Y[1] = 0.1;
//
//		double X = MinX + StepSize + StepSize; // Start from the second in the sequence.
//
//		for (int i = 2; i < Barrier; i++)
//		{
//			//Y[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass * (Energy - HOPlus(X - StepSize, l, k, Mass)) / 6) * Y[i - 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - HOPlus(X - StepSize - StepSize, l, k, Mass)) / 12) * Y[i - 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - HOPlus(X, l, k, Mass)) / 12);
//			//X += StepSize; // Increment X.
//			Y[i] = ((2 - 5 * StepSize * StepSize * RadialG(X - StepSize, l, k, Mass, Energy) / 6) * Y[i - 1] - (1 + StepSize * StepSize * RadialG(X - StepSize - StepSize, l, k, Mass, Energy) / 12) * Y[i - 2]) / (1 + StepSize * StepSize * RadialG(X, l, k, Mass, Energy) / 12);
//			X += StepSize; // Increment X.
//		}
//		Y_Right[NumberOfSteps - Barrier + 1] = 0; // Necessary initial condition;
//		Y_Right[NumberOfSteps - Barrier] = 0.1; //* (double)Parity; // Positive for even, negative for odd.
//
//		double X_Right = MaxX - StepSize - StepSize; // Begin two steps down.
//
//		for (int i = NumberOfSteps - Barrier - 1; i >= 0; i--) // This builds from the right, but I keep the index in the same direction as Y, so we index down. Y_Right[0] is the leftmost point.
//		{
//			//Y_Right[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass * (Energy - HOPlus(X_Right + StepSize, l, k, Mass)) / 6) * Y_Right[i + 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - HOPlus(X_Right + StepSize + StepSize, l, k, Mass)) / 12) * Y_Right[i + 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - HOPlus(X_Right, l, k, Mass)) / 12);
//			//X_Right -= StepSize; // Increment X (down)
//			Y_Right[i] = ((2 - 5 * StepSize * StepSize * RadialG(X_Right + StepSize, l, k, Mass, Energy) / 6) * Y_Right[i + 1] - (1 + StepSize * StepSize * RadialG(X_Right + StepSize + StepSize, l, k , Mass, Energy) / 12) * Y_Right[i + 2]) / (1 + StepSize * StepSize * RadialG(X_Right, l, k, Mass, Energy) / 12);
//			X_Right -= StepSize; // Increment X (down)
//		}
//
//		if (RelativeHeight / (Y[Barrier - 1] - Y_Right[0]) < 0 && EStepCount > 1) // Means Y and Y_Right have flipped and this is the second E step.
//		{
//			EIncrement /= -10; // Reverse energy step and make it finer.
//		}
//		RelativeHeight = Y[Barrier - 1] - Y_Right[0];
//		cout << Energy << "\t" << RelativeHeight << "\t" << Y[Barrier - 1] << "\t" << Y_Right[0] << endl;
//		//Sleep(200);
//
//		if (fabs(RelativeHeight) < 10E-10 || fabs(EIncrement) < 10E-6)
//		{
//			for (int i = Barrier; i < NumberOfSteps + 1; i++)
//			{
//				Y[i] = Y_Right[i - Barrier + 1]; // Means Y_Right[1], ...
//			}
//			delete[] Y_Right; // Clear memory allocated to Y_Right.
//
//			int DStep = Barrier / 4; // I need a point to check the derivative. This number is ideally right after the last peak in the wave, the perfect spot to calculated the derivative.
//			if (DStep < 4) DStep = 4; // If the Step is too small then just make it 3.
//			/* Check if derivative is continous. If not, does not generate EF file. */ // This works but the program is too unstable for this to be reliable.
//			if ((Y[Barrier - 1] - Y[Barrier - DStep]))// * (double)Parity > 0) // Means derivative left of barrier has different sign than derivative to the right of the barrier.
//			{
//				delete[] Y;
//				return Energy; // Do not generate EF file. But I still need to return the found energy to increment the next one.
//			}
//			break;
//		}
//
//		delete[] Y_Right; // Clear memory allocated to Y_Right.
//	} // End while loop
//
//	//string EvenOrOdd; // Generate the output after checking the derivative.
//	//Parity == 1 ? EvenOrOdd = "Even" : EvenOrOdd = "Odd"; // Even or odd solutions.
//	string IStr = to_string(Index);
//	string OutputName = "RadialModEF_" + IStr + ".txt";
//	ofstream PrintEF(OutputName);
//
//	/* This normalizes the eigenfunction */
//	double IntegralPP = 0; // Holds the integral
//	for (int i = 0; i < NumberOfSteps + 1; i++)
//	{
//		IntegralPP += (Y[i] * Y[i]) * fabs(StepSize); // Approximations integral by a finite sum.
//	}
//	for (int i = 0; i < NumberOfSteps + 1; i++)
//	{
//		Y[i] /= sqrt(IntegralPP); // Divide by square root of integral for normalization.
//	}
//
//	/* This prints the eigenfunction */
//	for (int i = 0; i < NumberOfSteps + 1; i++)
//	{
//		PrintEF << MinX + i * StepSize << "\t" << Y[i] << endl;
//	}
//	delete[] Y; // Clear memory allocated to Y
//
//	return Energy; // Eigenvalue for main function to print.
//}

double BarrierWellEF(double Radius, double Energy, double Mass, double V, int Index, int NumberOfSteps)
{
	double MinX = -Radius;
	double MaxX = Radius;
	double L = 2 * Radius;

	double EIncrement = 0.0025 / (Mass * ((MaxX - MinX) * (MaxX - MinX))); // Increments energy. Using some prior knowledge, I know this is small enough.

	double StepSize = (MaxX - MinX) / (double)NumberOfSteps; 

	/* Set up output file name */
	string IStr = to_string(Index);
	string OutputName = "BarrierWell_" + IStr + ".txt"; // Name of output for n'th eigenfunction
	ofstream PrintEF(OutputName); // Output file which will store plot data for eigenvectors.

	/* Y is the array that stores eigenfunction y values. I don't need to store all these values, but it makes normalization easier. */
	double* Y = new double[NumberOfSteps + 1];
	Y[NumberOfSteps] = 1; // Initializes the loop.

	double DerivStep1 = 1; // These will hold values to calculate the sign of the derivative. In this case, we only need to check their relative signs.
	double DerivStep2 = 1;

	int DerivStepCount = 0; // Tells us how many steps we have taken.

	while (fabs(Y[NumberOfSteps]) > 1E-30) // Boundary Condition. 0 at barrier.
	{
		Energy += EIncrement; // Take energy step
		Y[0] = 0; // Boundary condition
		Y[1] = 0.01; // Sign affects the phase.

		DerivStepCount++;

		double X = MinX + StepSize + StepSize; // Start at step 2.

		for (int i = 2; i < NumberOfSteps + 1; i++) // y(x) is the eigenfunction. a(x) = (2m/h)(En - V(x))
		{
			Y[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass *(Energy - BarrierWell(X - StepSize, L, V)) / 6) * Y[i - 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - BarrierWell(X - 2 * StepSize, L, V)) / 12) * Y[i - 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - BarrierWell(X, L, V)) / 12);
			X += StepSize; // Increment X.
		}
		DerivStep1 = DerivStep2; // Shift up derivative variables.
		DerivStep2 = Y[NumberOfSteps];
		if (DerivStep1 * DerivStep2 < 0 && DerivStepCount > 1) // Means sign change (i.e. we missed the boundary condition at 0) and that we have taken two steps.
		{
			EIncrement /= -10; // Reverse the energy step and make it finer.
		}
		cout << Y[NumberOfSteps] << endl;
		if (fabs(EIncrement) < 10E-10)
		{
			break;
		}
	}

	/* This normalizes the eigenfunction */
	double IntegralPP = 0; // Holds the integral
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		IntegralPP += (Y[i] * Y[i]) * fabs(StepSize); // Approximations integral by a finite sum.
	}
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		Y[i] /= sqrt(IntegralPP); // Divide by square root of integral for normalization.
	}

	/* This prints the eigenfunction and eigenvalue */
	//PrintEF << "Eigenvalue = " << Energy << endl;
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		PrintEF << MinX + i * StepSize << "\t" << Y[i] << endl;
	}

	/* Clear memory allocated to Y */
	delete[] Y;

	return Energy;

}

double LinearPotentialEF(double MinX, double MaxX, double Energy, double Mass, double b, int Index, int NumberOfSteps)
{

	double EIncrement = Energy / 100;

	double StepSize = (MaxX - MinX) / (double)NumberOfSteps; // h in HW statement

	/* Set up output file name */
	string IStr = to_string(Index);
	string OutputName = "Linear_" + IStr + ".txt"; // Name of output for n'th eigenfunction
	ofstream PrintEF(OutputName); // Output file which will store plot data for eigenvectors.

	/* Y is the array that stores eigenfunction y values. I don't need to store all these values, but it makes normalization easier. */
	double* Y = new double[NumberOfSteps + 1];
	Y[0] = 1; // Initializes the loop.

	double DerivStep1 = 1; // These will hold values to calculate the sign of the derivative. In this case, we only need to check their relative signs.
	double DerivStep2 = 1;

	int DerivStepCount = 0; // Tells us how many steps we have taken.

	while (fabs(Y[0]) > 10E-10)
	{
		Energy += EIncrement; // Take energy step
		Y[NumberOfSteps] = 0; // Necessary initial condition.
		Y[NumberOfSteps - 1] = 0.1;

		DerivStepCount++;
		double X = MaxX - StepSize - StepSize; // Start from the second in the sequence.

		for (int i = NumberOfSteps - 2; i >= 0; i--)
		{
			Y[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass * (Energy - LinearPotential(X + StepSize, b)) / 6) * Y[i + 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - LinearPotential(X + StepSize + StepSize, b)) / 12) * Y[i + 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - LinearPotential(X, b)) / 12);
			X -= StepSize; // Increment X (down)
		}

		DerivStep1 = DerivStep2; // Shift up derivative variables.
		DerivStep2 = Y[0];
		cout << Y[0] << endl;
		if (DerivStep1 * DerivStep2 < 0 && DerivStepCount > 1) // Means sign change (i.e. we missed the boundary condition at 0) and that we have taken two steps.
		{
			EIncrement /= -10; // Reverse the energy step and make it finer.
		}
		if (EIncrement < 10E-6)
		{
			break;
		}
	} // End while loop

	/* This normalizes the eigenfunction */
	double IntegralPP = 0; // Holds the integral
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		IntegralPP += (Y[i] * Y[i]) * fabs(StepSize); // Approximations integral by a finite sum.
	}
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		Y[i] /= sqrt(IntegralPP); // Divide by square root of integral for normalization.
	}

	/* This prints the eigenfunction */
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		PrintEF << MinX + i * StepSize << "\t" << Y[i] << endl;
	}
	delete[] Y; // Clear memory allocated to Y

	return Energy; // Eigenvalue for main function to print.
}

double QuarticPotentialEF(double MinX, double MaxX, double Energy, double Mass, double c, int Parity, int Index, int NumberOfSteps)
{

	double StepSize = (MaxX - MinX) / (double)NumberOfSteps;
	double EIncrement = fabs(Energy) / 1000; // Energy step size.

	double* Y = new double[NumberOfSteps + 1]; // Array that holds the eigenfunction

	int EStepCount = 0; // Counts steps in the loop
	double RelativeHeight = 10;

	while (true) // While the two eigenfunctions are discontinuous.
	{
		EStepCount++;
		Energy += EIncrement; // Increments the energy.

		double BarrierX = pow(Energy, 0.25);
		int Barrier = 0;
		for (int j = 0; j < NumberOfSteps; j++)
		{
			if (MinX + (double)j * StepSize >= BarrierX)
			{
				break;
			}
			Barrier++;
		}
		double* Y_Right = new double[NumberOfSteps + 1 - Barrier + 1];
		Y = new double[NumberOfSteps + 1];

		Y_Right[0] = 1; // Initialize to begin loop;
		Y[Barrier - 1] = 2;

		Y[0] = 0; // Necessary initial condition.
		Y[1] = 0.1;

		double X = MinX + StepSize + StepSize; // Start from the second in the sequence.

		for (int i = 2; i < Barrier; i++)
		{
			Y[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass * (Energy - QuarticPotential(X - StepSize, c)) / 6) * Y[i - 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - QuarticPotential(X - StepSize - StepSize, c)) / 12) * Y[i - 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - QuarticPotential(X, c)) / 12);
			X += StepSize; // Increment X.
		}
		Y_Right[NumberOfSteps - Barrier + 1] = 0; // Necessary initial condition;
		Y_Right[NumberOfSteps - Barrier] = 0.1 * (double)Parity; // Positive for even, negative for odd.

		double X_Right = MaxX - StepSize - StepSize; // Begin two steps down.

		for (int i = NumberOfSteps - Barrier - 1; i >= 0; i--) // This builds from the right, but I keep the index in the same direction as Y, so we index down. Y_Right[0] is the leftmost point.
		{
			Y_Right[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass * (Energy - QuarticPotential(X_Right + StepSize, c)) / 6) * Y_Right[i + 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - QuarticPotential(X_Right + StepSize + StepSize, c)) / 12) * Y_Right[i + 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - QuarticPotential(X_Right, c)) / 12);
			X_Right -= StepSize; // Increment X (down)
		}

		if (RelativeHeight / (Y[Barrier - 1] - Y_Right[0]) < 0 && EStepCount > 1) // Means Y and Y_Right have flipped and this is the second E step.
		{
			EIncrement /= -10; // Reverse energy step and make it finer.
		}
		RelativeHeight = Y[Barrier - 1] - Y_Right[0];
		cout << Energy << "\t" << RelativeHeight << "\t" << Y[Barrier - 1] << "\t" << Y_Right[0] << endl;

		if (fabs(RelativeHeight) < 10E-10 || fabs(EIncrement) < 10E-6)
		{
			for (int i = Barrier; i < NumberOfSteps + 1; i++)
			{
				Y[i] = Y_Right[i - Barrier + 1]; // Means Y_Right[1], ...
			}
			delete[] Y_Right; // Clear memory allocated to Y_Right.

			int DStep = (NumberOfSteps - 2 * (NumberOfSteps - Barrier + 1)) / (Index * 2); // I need a point to check the derivative. This number is ideally right after the last peak in the wave, the perfect spot to calculated the derivative.
			if (DStep < 4) DStep = 4; // If the Step is too small then just make it 3.
			/* Check if derivative is continous. If not, does not generate EF file. */ // This works but the program is too unstable for this to be reliable.
			if ((Y[Barrier - 1] - Y[Barrier - DStep]) * (double)Parity > 0) // Means derivative left of barrier has different sign than derivative to the right of the barrier.
			{
				delete[] Y;
				return Energy; // Do not generate EF file. But I still need to return the found energy to increment the next one.
			}
			break;
		}

		delete[] Y_Right; // Clear memory allocated to Y_Right.
	} // End while loop

	string EvenOrOdd; // Generate the output after checking the derivative.
	Parity == 1 ? EvenOrOdd = "Even" : EvenOrOdd = "Odd"; // Even or odd solutions.
	string IStr = to_string(Index);
	string OutputName = "QuarticEF_" + IStr + "_" + EvenOrOdd + ".txt"; // SInce I generate even and odd solutions at the same time, I need to distinguish the output files.
	ofstream PrintEF(OutputName);

	/* This normalizes the eigenfunction */
	double IntegralPP = 0; // Holds the integral
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		IntegralPP += (Y[i] * Y[i]) * fabs(StepSize); // Approximations integral by a finite sum.
	}
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		Y[i] /= sqrt(IntegralPP); // Divide by square root of integral for normalization.
	}

	/* This prints the eigenfunction */
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		PrintEF << MinX + i * StepSize << "\t" << Y[i] << endl;
	}
	delete[] Y; // Clear memory allocated to Y

	return Energy; // Eigenvalue for main function to print.
}

double ParticleInBoxEF(double LeftBound, double RightBound, double WellHeight, double Mass, double Energy, int Index, int NumberOfSteps) // Prints the n'th eigenstate for Particle in a Box
{
	double MinX = LeftBound;
	double MaxX = RightBound;

	double EIncrement = 0.0025 / (Mass * ((RightBound - LeftBound) * (RightBound - LeftBound))); // Increments energy. Using some prior knowledge, I know this is small enough.

	double StepSize = (MaxX - MinX) / (double)NumberOfSteps; // h in HW statement

	/* Set up output file name */
	string IStr = to_string(Index);
	string OutputName = "ParticleInABoxEF_" + IStr + ".txt"; // Name of output for n'th eigenfunction
	ofstream PrintEF(OutputName); // Output file which will store plot data for eigenvectors.

	/* Y is the array that stores eigenfunction y values. I don't need to store all these values, but it makes normalization easier. */
	double* Y = new double[NumberOfSteps + 1];
	Y[NumberOfSteps] = 1; // Initializes the loop.

	double DerivStep1 = 1; // These will hold values to calculate the sign of the derivative. In this case, we only need to check their relative signs.
	double DerivStep2 = 1;

	int DerivStepCount = 0; // Tells us how many steps we have taken.

	while (fabs(Y[NumberOfSteps]) > 1E-10) // Boundary Condition. 0 at barrier.
	{
		Energy += EIncrement; // Take energy step
		Y[0] = 0; // Boundary condition
		Y[1] = 0.01; // Sign affects the phase.

		DerivStepCount++;

		double X = MinX + StepSize + StepSize; // Start at step 2.

		for (int i = 2; i < NumberOfSteps + 1; i++) // y(x) is the eigenfunction. a(x) = (2m/h)(En - V(x))
		{
			Y[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass *(Energy - SquareWell(X - StepSize, WellHeight, 0, LeftBound, RightBound)) / 6) * Y[i - 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X - 2 * StepSize, WellHeight, 0, LeftBound, RightBound)) / 12) * Y[i - 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X, WellHeight, 0, LeftBound, RightBound)) / 12);
			X += StepSize; // Increment X.
		}
		DerivStep1 = DerivStep2; // Shift up derivative variables.
		DerivStep2 = Y[NumberOfSteps];
		if (DerivStep1 * DerivStep2 < 0 && DerivStepCount > 1) // Means sign change (i.e. we missed the boundary condition at 0) and that we have taken two steps.
		{
			EIncrement /= -10; // Reverse the energy step and make it finer.
		}
	}

	/* This normalizes the eigenfunction */
	double IntegralPP = 0; // Holds the integral
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		IntegralPP += (Y[i] * Y[i]) * fabs(StepSize); // Approximations integral by a finite sum.
	}
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		Y[i] /= sqrt(IntegralPP); // Divide by square root of integral for normalization.
	}

	/* This prints the eigenfunction and eigenvalue */
	PrintEF << "Eigenvalue = " << Energy << endl;
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		PrintEF << MinX + i * StepSize << "\t" << Y[i] << endl;
	}

	/* Clear memory allocated to Y */
	delete[] Y;

	return Energy; // Eigenvalue to be printed in main fuction.
}

double FiniteWellEF(double MinX, double MaxX, double WellRadius, double V, double Energy, double Mass, int Parity, int Index, int NumberOfSteps) // Prints eigenfunction for a finite well.
{
	double StepSize = (MaxX - MinX) / (double)NumberOfSteps;
	double EIncrement = fabs(Energy - V) / 1000; // Energy step size.

	double* Y = new double[NumberOfSteps + 1]; // Array that holds the eigenfunction

	int Barrier = 0; // Index that marks the point right before the barrier.
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{	
		if ((MinX + (double)i * StepSize) * ((MaxX - MinX) / fabs(MaxX - MinX)) >= WellRadius) // I account for the case that we start from right and go left by generalizing for different signs.
		{
			break;
		}
		Barrier++;
	}
	double* Y_Right = new double[NumberOfSteps + 1 - Barrier + 1]; // This stores the eigenfunction where the algorithm starts from the right, to impose boundary conditions.

	Y_Right[0] = 1; // Initialize to begin loop;
	Y[Barrier - 1] = 2; 

	if (Energy > 0) // Means we are outside the well and don't have to worry about boundary conditions. The algorithm simple.
	{
		Y[0] = 0;
		Y[1] = 0.1;

		double X = MinX + StepSize + StepSize;

		for (int i = 2; i < NumberOfSteps + 1; i++) // y(x) is the eigenfunction. a(x) = (2m/h^2)(E - V(x))
		{
			Y[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X - StepSize, 0, V, -WellRadius, WellRadius)) / 6) * Y[i - 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X - StepSize - StepSize, 0, V, -WellRadius, WellRadius)) / 12) * Y[i - 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X, 0, V, -WellRadius, WellRadius)) / 12);
			X += StepSize; // Increment X.
		}
	}
	/* The case where we go into a barrier is more complicated. There are two solutions, e^(kx} and e^(-kx). The +kx solution is discarded by imposing normalizability. This is handled if we start
	inside a barrier by setting the first value to 0. However, we cannot do this if we end inside a barrier. Unfortunately, this algorithm picks the +kx solution unless we set the normalizability
	condition. In this case, we start inside and end inside a barrier, so it cannot be avoided that the algorithm will go into a barrier and pick the wrong solution. The way I have dealt with this
	is by running two algorithms. First runs through the well up to the barrier. Then the second algorithm starts from inside the barrier and stops at the well. The second algorithm allows me to
	impose the boundary condition and have a -kx solution. I then impose continuity to find energy and then check that derivatives agree. */
	else // Means we are inside the well and must treat the boundary conditions.
	{
		double RelativeHeight = Y[Barrier - 1] - Y_Right[0]; // Used to measure continuity.
		int EStepCount = 0; // Counts steps in the loop

		while (fabs(Y_Right[0] - Y[Barrier - 1]) > 0.001) // While the two eigenfunctions are discontinuous.
		{
			EStepCount++;
			Energy += EIncrement; // Increments the energy.

			Y[0] = 0; // Necessary initial condition.
			Y[1] = 0.1;

			bool RadiusIsTooSmall = false; // Make sure Y_Right will be large enough to use a two sides approach.
			if (MaxX - 2 * StepSize < WellRadius * ((MaxX - MinX) / fabs(MaxX - MinX)))
			{
				RadiusIsTooSmall = true;
				Barrier = NumberOfSteps + 1; // We'll just use Y in this case.
			}

			double X = MinX + StepSize + StepSize; // Start from the second in the sequence.

			for (int i = 2; i < Barrier; i++)
			{
				Y[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X - StepSize, 0, V, -WellRadius, WellRadius)) / 6) * Y[i - 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X - StepSize - StepSize, 0, V, -WellRadius, WellRadius)) / 12) * Y[i - 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X, 0, V, -WellRadius, WellRadius)) / 12);
				X += StepSize; // Increment X.
			}

			if (RadiusIsTooSmall == false)
			{
				Y_Right[NumberOfSteps - Barrier + 1] = 0; // Necessary initial condition;
				Y_Right[NumberOfSteps - Barrier] = 0.1 * (double)Parity; // Positive for even, negative for odd.

				double X_Right = MaxX - StepSize - StepSize; // Begin two steps down.

				for (int i = NumberOfSteps - Barrier - 1; i >= 0; i--) // This builds from the right, but I keep the index in the same direction as Y, so we index down. Y_Right[0] is the leftmost point.
				{
					Y_Right[i] = ((2 - 5 * StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X_Right + StepSize, 0, V, -WellRadius, WellRadius)) / 6) * Y_Right[i + 1] - (1 + StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X_Right + StepSize + StepSize, 0, V, -WellRadius, WellRadius)) / 12) * Y_Right[i + 2]) / (1 + StepSize * StepSize * 2 * Mass * (Energy - SquareWell(X_Right, 0, V, -WellRadius, WellRadius)) / 12);
					X_Right -= StepSize; // Increment X (down)
				}

				if (RelativeHeight / (Y[Barrier - 1] - Y_Right[0]) < 0 && EStepCount > 1) // Means Y and Y_Right have flipped and this is the second E step.
				{
					EIncrement /= -10; // Reverse energy step and make it finer.
				}
				RelativeHeight = Y[Barrier - 1] - Y_Right[0];
				cout << Energy << "\t" << RelativeHeight << endl;
				if (EIncrement < 1E-20) // Means energy is close enough.
				{
					break;
				}
			} // End if radius is too small.
			else // Generate Y and break loop.
			{
				break;
			}
		} // End while loop
	} // End else ( E > 0 );

	int DStep = (NumberOfSteps - 2 * (NumberOfSteps - Barrier + 1)) / (Index * 2); // I need a point to check the derivative. This number is ideally right after the last peak in the wave, the perfect spot to calculated the derivative.
	if (DStep < 4) DStep = 4; // If the Step is too small then just make it 3.
	/* Check if derivative is continous. If not, does not generate EF file. */ // This works but the program is too unstable for this to be reliable.
	if ((Y[Barrier - 1] - Y[Barrier - DStep]) * (double)Parity > 0) // Means derivative left of barrier has different sign than derivative to the right of the barrier.
	{
		delete[] Y;
		return Energy; // Do not generate EF file. But I still need to return the found energy to increment the next one.
	}

	string EvenOrOdd; // Generate the output after checking the derivative.
	Parity == 1 ? EvenOrOdd = "Even" : EvenOrOdd = "Odd"; // Even or odd solutions.
	string IStr = to_string(Index);
	string OutputName = "FiniteWell_" + IStr + "_" + EvenOrOdd + ".txt"; // SInce I generate even and odd solutions at the same time, I need to distinguish the output files.
	ofstream PrintEF(OutputName);

	for (int i = Barrier; i < NumberOfSteps + 1; i++)
	{
		Y[i] = Y_Right[i - Barrier + 1]; // Means Y_Right[1], ...
	}

	delete[] Y_Right; // Clear memory allocated to Y_Right.

	/* This normalizes the eigenfunction */
	double IntegralPP = 0; // Holds the integral
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		IntegralPP += (Y[i] * Y[i]) * fabs(StepSize); // Approximations integral by a finite sum.
	}
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		Y[i] /= sqrt(IntegralPP); // Divide by square root of integral for normalization.
	}

	/* This prints the eigenfunction */
	for (int i = 0; i < NumberOfSteps + 1; i++)
	{
		PrintEF << MinX + i * StepSize << "\t" << Y[i] << endl;
	}
	delete[] Y; // Clear memory allocated to Y

	return Energy; // Eigenvalue for main function to print.
}

int main()
{
	cout << "=========================================================================" << endl;
	cout << "=========== NUMEROV'S METHOD FOR VARIOUS TOY QUANTUM PROBLEMS ===========" << endl;
	cout << "=========================================================================" << endl;
	cout << "=========== AUTHOR: HENRY K. TRAN =======================================" << endl;
	cout << "=========== CLASS: CHEMISTRY 6510 AUTUMN 2015 ===========================" << endl;
	cout << "=========== HOMEWORK PROBLEM 2.6 ========================================" << endl;
	cout << "=========================================================================" << "\n" << endl;

	int Options;
	cout << "Choose problem:" << "\n" << "1. Print first n eigenstates of Particle in a Box" << "\n" << "2. Print eigenstate for specific E of Finite Potential Well" << endl;
	cin >> Options;

	clock_t StartTime;
	double Duration;

	if (Options == 1) // As a note, I could probably implement finite square well easy enough. I just have to break the algorithm into a left and right part.
	{
		double LeftBound; //Taken as range.
		double RightBound;
		int NumberOfSteps;// Steps used in Numerov Method.
		double Mass;
		double Energy;
		int Trials; // Number of eigenvalues to compute.
		//int MaxN;
		cout << "=========================================================================" << endl;
		cout << "=========== Problem 1 ===================================================" << endl;
		cout << "=========================================================================" << endl;
		cout << "                                                                         " << endl;
		cout << "            |                                               |            " << endl;
		cout << "            |                                               |            " << endl;
		cout << "            |                                               |            " << endl;
		cout << "            |              Infinte Square Well              |            " << endl;
		cout << "            |                                               |            " << endl;
		cout << "            |                                               |            " << endl;
		cout << "     V = 0  |_______________________________________________|            " << "\n" << endl;
		cout << "=========================================================================" << "\n" << endl;
		cout << "Enter left bound of well (au):" << endl;
		cin >> LeftBound;
		cout << "Enter right bound of well (au):" << endl;
		cin >> RightBound;
		cout << "Enter mass of particle (au):" << endl;
		cin >> Mass;
		cout << "Enter number of eigenstates to calculate:" << endl;
		cin >> Trials;
		cout << "Enter number of steps in algorithm:" << endl;
		cin >> NumberOfSteps;
		cout << "Initializing..." << endl;

		StartTime = clock();

		ofstream PrintEV("ParticleInABoxEV.txt"); // Output file which will store the eigenvalues.
		double EV; // Holds result of calculation. I will need this to increment energy.
		Energy = 2 * 2 / (8 * Mass * (RightBound - LeftBound) * (RightBound - LeftBound)); // Using some prior knowledge, I know this is less than the ground state energy and a decent starting point.

		for (int i = 0; i < Trials; i++)
		{
			EV = ParticleInBoxEF(LeftBound, RightBound, 1000, Mass, Energy, i + 1, NumberOfSteps);
			PrintEV << setprecision(10) << EV << endl;
			Energy = EV + 0.25 / (Mass * (RightBound - LeftBound) * (RightBound - LeftBound)); // This will start us away from the current EV but not over the next EV.
			cout << "Eigenvalue and Eigenvector " << i + 1 << " completed." << endl;
		}
	}

	if (Options == 2)
	{
		double V; // Potential of the middle part.
		double WellRadius; // Radius of the middle part.
		double Energy;
		double MinX; // Range
		double MaxX;
		double Mass;
		int Trials; // Number of times to find EV, hopefully the lowest.
		int NumberOfSteps; // Steps used in Numerov Method
		cout << "=========================================================================" << endl;
		cout << "=========== Problem 2 ===================================================" << endl;
		cout << "=========================================================================" << endl;
		cout << "_____________                                               _____________" << endl;
		cout << "            |  V = 0                                        |            " << endl;
		cout << "            |                                               |            " << endl;
		cout << "            |                                               |            " << endl;
		cout << "            |               Finite Square Well              |            " << endl;
		cout << "            |                                               |            " << endl;
		cout << "            |                                               |            " << endl;
		cout << "            |_______________________________________________|            " << "\n" << endl;
		cout << "=========================================================================" << "\n" << endl;
		cout << "Enter initial x value (au):" << endl;
		cin >> MinX;
		cout << "Enter final x value (au):" << endl;
		cin >> MaxX;
		cout << "Enter radius of well (au):" << endl;
		cin >> WellRadius; // a in the HW statement
		cout << "Enter mass of particle (au):" << endl;
		cin >> Mass;
		cout << "Enter potential (Note: Not the drop, the actual potential) (au):" << endl; // As in, I will take the sign as is.
		cin >> V; // Really should be negative.
		cout << "Enter number of trials:" << endl;
		cin >> Trials;
		cout << "Enter number of steps in algorithm:" << endl;
		cin >> NumberOfSteps;
		cout << "Initializing..." << endl;

		StartTime = clock();

		for (int j = -1; j <= 1; j += 2) // j = -1 corresponds to odd, j = 1 corresponds to evens.
		{
			Energy = V * 0.999; // First guess which assumes V is negative
			string EorO;
			j == 1 ? EorO = "Even" : EorO = "Odd";
			string OutputName = "FiniteWellEV_" + EorO + ".txt"; // Differentiate the outputs since I generate both at the same time.
			ofstream PrintEV(OutputName);
			double EV;

			for (int i = 0; i < Trials; i++)
			{
				if (Energy > 0) // Means we are over the barrier.
				{
					break;
				}
				EV = FiniteWellEF(MinX, MaxX, WellRadius, V, Energy, Mass, j, i + 1, NumberOfSteps);
				PrintEV << setprecision(5) << EV << endl;
				Energy = EV + 0.25 / (Mass * (2 * WellRadius) * (2 * WellRadius));
				cout << "Eigenvalue and Eigenvector " << i + 1 << " completed." << endl;
			}
		}
	}
	if (Options == 3)
	{
		// 4.32
		for (int j = 1; j >= -1; j -= 2) // j = -1 corresponds to odd, j = 1 corresponds to evens.
		{
			double MinX = -3;
			double MaxX = 3;
			double Mass = 1;
			double c = 1;
			int NumberOfSteps = 500;
			int Trials = 5;

			double Energy = 0.1;
			string EorO;
			j == 1 ? EorO = "Even" : EorO = "Odd";
			string OutputName = "QuarticEV_" + EorO + ".txt"; // Differentiate the outputs since I generate both at the same time.
			ofstream PrintEV(OutputName);
			double EV;

			for (int i = 0; i < Trials; i++)
			{
				EV = QuarticPotentialEF(MinX, MaxX, Energy, Mass, c, j, i + 1, NumberOfSteps);
				PrintEV << setprecision(5) << EV << endl;
				Energy = EV + 1;
				cout << "Eigenvalue and Eigenvector " << i + 1 << " completed." << endl;
			}
		}
	}
	if (Options == 4)
	{
		double MaxX = 10;
		double Energy = 0.1;
		double Mass = 1;
		double b = 1;
		int NumberOfSteps = 100;
		int Trials = 6;
		double EV;
		ofstream PrintEV("LinearPotentialEV.txt");
		for (int i = 0; i < Trials; i++)
		{
			EV = LinearPotentialEF(0, MaxX, Energy, Mass, b, i, NumberOfSteps);
			PrintEV << setprecision(5) << EV << endl;
			Energy = EV + 0.25;
			cout << "Eigenvalue and Eigenvector " << i + 1 << " completed." << endl;
		}
	}
	if (Options == 5)
	{
		double L = 10;
		
		double Mass = 1;
		double V = 100 / (L * L);
		double Radius = L / 2;
		int NumberOfSteps = 500;
		int Trials = 10;
		double EV;
		double Energy = 2 * 2 / (8 * Mass * L * L);
		ofstream PrintEV("BarrierWEllEV.txt");
		for (int i = 0; i < Trials; i++)
		{
			EV = BarrierWellEF(Radius, Energy, Mass, V, i, NumberOfSteps);
			PrintEV << setprecision(10) << EV << endl;
			Energy = EV + 0.0025 / (Mass * L * L);
			cout << "Eigenvalue and Eigenvector " << i + 1 << " completed." << endl;
		}
	}
	if (Options == 6)
	{
		double l = 1;
		double k = 1;
		double Mass = 1;
		double Energy = 0.6; // Min is k / 2 + 1 / 2 * l * (l + 1)
		int NumberOfSteps = 500;
		double MinX = 0.00000000001;
		double MaxX = 6.00000000001;
		int Trials = 4;
		double EV;
		string EVName = "RadialEV_l=" + to_string(l) + ".txt";
		ofstream PrintEV(EVName);
		for (int i = 0; i < Trials; i++)
		{
			EV = RadialEF(MinX, MaxX, Energy, Mass, l, k, 1, i, NumberOfSteps);
			PrintEV << setprecision(10) << EV << endl;
			Energy = EV + 0.1; // +2 * (i + 1) / (2 * (i + 1) * (i + 1) * (i + 2)*(i + 2));
			cout << "Eigenvalue and Eigenvector " << i + 1 << " completed." << endl;
		}

	}
	else
	{
		cout << "No problems selected: Terminating program." << endl;
		return 0;
	}

	//Duration = (clock() - StartTime) / (double)CLOCKS_PER_SEC; // Calculates the time of the process.
	//cout << "\n" << "Program took " << Duration << " seconds.";

	return 0;
}