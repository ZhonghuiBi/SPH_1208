#pragma once
#include <algorithm>
//#include<corecrt_math_defines.h>
#include "Vec3D.hpp"
#include "Matrix.hpp"
#include <map>
#include <vector>
#include <iomanip>
#include <math.h>

#define M_PI       3.14159265358979323846   // pi

class MatVec
{

public:

	//* Left multiplication. {B} = {A}^T*[M].  */
	static void Mult(Vec3_t const & A, Mat3_t const & M, Vec3_t & B)
	{
		B(0) = A(0)*M(0, 0) + A(1)*M(1, 0) + A(2)*M(2, 0);
		B(1) = A(0)*M(0, 1) + A(1)*M(1, 1) + A(2)*M(2, 1);
		B(2) = A(0)*M(0, 2) + A(1)*M(1, 2) + A(2)*M(2, 2);
	}
	//
	//** Right multiplication. {B} = [M]*{A}.  */
	static void Mult(Mat3_t const & M, Vec3_t const & A, Vec3_t & B)
	{
		B(0) = A(0)*M(0, 0) + A(1)*M(0, 1) + A(2)*M(0, 2);
		double c = A(0)*M(0, 0) + A(1)*M(0, 1) + A(2)*M(0, 2);

		B(1) = A(0)*M(1, 0) + A(1)*M(1, 1) + A(2)*M(1, 2);

		B(2) = A(0)*M(2, 0) + A(1)*M(2, 1) + A(2)*M(2, 2);

	}
	//
	//** Matrix multiplication. */
	static void Mult(Mat3_t const & A, Mat3_t const & B, Mat3_t & M)
	{
		M(0, 0) = A(0, 2)*B(2, 0) + A(0, 1)*B(1, 0) + A(0, 0)*B(0, 0);  M(0, 1) = A(0, 2)*B(2, 1) + A(0, 1)*B(1, 1) + A(0, 0)*B(0, 1);  M(0, 2) = A(0, 2)*B(2, 2) + A(0, 1)*B(1, 2) + A(0, 0)*B(0, 2);
		M(1, 0) = A(1, 2)*B(2, 0) + B(1, 0)*A(1, 1) + B(0, 0)*A(1, 0);  M(1, 1) = A(1, 2)*B(2, 1) + A(1, 1)*B(1, 1) + B(0, 1)*A(1, 0);  M(1, 2) = A(1, 2)*B(2, 2) + A(1, 1)*B(1, 2) + B(0, 2)*A(1, 0);
		M(2, 0) = A(2, 2)*B(2, 0) + B(1, 0)*A(2, 1) + B(0, 0)*A(2, 0);  M(2, 1) = B(2, 1)*A(2, 2) + B(1, 1)*A(2, 1) + B(0, 1)*A(2, 0);  M(2, 2) = A(2, 2)*B(2, 2) + B(1, 2)*A(2, 1) + B(0, 2)*A(2, 0);
	}

	static void set_to_zero(Mat3_t & M)
	{
		M(0, 0) = 0.0;
		M(0, 1) = 0.0;
		M(0, 2) = 0.0;
		M(1, 0) = 0.0;
		M(1, 1) = 0.0;
		M(1, 2) = 0.0;
		M(2, 0) = 0.0;
		M(2, 1) = 0.0;
		M(2, 2) = 0.0;
	}

	static void set_to_zero(Vec3_t & V)
	{
		V(0) = 0;
		V(1) = 0;
		V(2) = 0;
	}
	/** Transpose.*/
	static void Trans(Mat3_t const & M, Mat3_t & Mt)
	{
		Mt(0, 0) = M(0, 0);   Mt(0, 1) = M(1, 0);   Mt(0, 2) = M(2, 0);
		Mt(1, 0) = M(0, 1);   Mt(1, 1) = M(1, 1);   Mt(1, 2) = M(2, 1);
		Mt(2, 0) = M(0, 2);   Mt(2, 1) = M(1, 2);   Mt(2, 2) = M(2, 2);
	}
	/** Eigenvalues and eigenvectors. NOTE: This function changes the matrix M. */
	/**
	* @brief                        Jacobi eigenvalue algorithm
	* @param matrix				    n*n array
	* @param dim					dim represent n
	* @param eigenvectors			n*n array
	* @param eigenvalues			n*1 array
	* @param precision   			precision requirements
	* @param max					max number of iterations
	* @return
	*/


	static void Eig(Mat3_t & M, Vec3_t & L)
	{
		// precision requirements
		double precision = 1e-10;
		// max number of iterations
		int max = 999;

		Mat3_t matrix_M = M;
		//current iteration
		int dim = matrix_M.getCols();
		// Creat matrix E
		Mat3_t eigenvectors_M(dim, 1);

		Mat3_t eigenvalues_M(dim, 1);


		int nCount = 0; //current iteration
		while (1)
		{
			//find the largest element on the off-diagonal line of the matrix
			double dbMax = matrix_M(0, 1);
			int nRow = 0;
			int nCol = 1;
			for (int i = 0; i < dim; i++) {			//row
				for (int j = 0; j < dim; j++) {		//column
					double d = fabs(matrix_M(i, j));
					if ((i != j) && (d > dbMax))
					{
						dbMax = d;
						nRow = i;
						nCol = j;
					}
				}
			}

			if (dbMax < precision)     //precision check 
				break;
			if (nCount > max)       //iterations check
				break;
			nCount++;

			double dbAij = matrix_M(nRow, nCol);//ij
			double dbAjj = matrix_M(nCol, nCol);//jj
			double dbAii = matrix_M(nRow, nRow);//ii
												//compute rotate angle
			double dbAngle;
			if (dbAii != dbAjj)
				dbAngle = 0.5*atan2(-2 * dbAij, (-dbAii + dbAjj));
			else
				dbAngle = M_PI / 4.0;

			double dbSinTheta = sin(dbAngle);
			double dbCosTheta = cos(dbAngle);
			double dbSin2Theta = sin(2 * dbAngle);
			double dbCos2Theta = cos(2 * dbAngle);

			matrix_M(nRow, nRow) = dbAii * dbCosTheta*dbCosTheta +
				dbAjj * dbSinTheta*dbSinTheta + dbAij* dbSin2Theta;
			matrix_M(nCol, nCol) = dbAii * dbSinTheta*dbSinTheta +
				dbAjj * dbCosTheta*dbCosTheta - dbAij * dbSin2Theta;
			matrix_M(nCol, nRow) = 0.5*(dbAjj - dbAii)*dbSin2Theta + dbAij * dbCos2Theta;
			matrix_M(nRow, nCol) = matrix_M(nCol, nRow);

			for (int k = 0; k < dim; k++)
			{
				if ((k != nCol) && (k != nRow))
				{
					dbMax = matrix_M(k, nRow);
					matrix_M(k, nRow) = matrix_M(k, nCol) * dbSinTheta + dbMax * dbCosTheta;
					matrix_M(k, nCol) = matrix_M(k, nCol) * dbCosTheta - dbMax * dbSinTheta;

					matrix_M(nRow, k) = matrix_M(k, nRow);
					matrix_M(nCol, k) = matrix_M(k, nCol);
				}
			}
			//for (int k = 0; k < dim; k++)
			//{
			//	if ((k != nCol) && (k != nRow))
			//	{
			//		dbMax = matrix_M(nRow, k);
			//		matrix_M(nRow,k) = matrix_M(nCol, k) * dbSinTheta + dbMax * dbCosTheta;
			//		matrix_M(nCol,k) = matrix_M(nCol, k) * dbCosTheta - dbMax * dbSinTheta;
			//	}
			//}
			// compute eigenvector
			for (int k = 0; k < dim; k++)
			{
				dbMax = eigenvectors_M(k, nRow);
				eigenvectors_M(k, nRow) = eigenvectors_M(k, nCol) * dbSinTheta + dbMax * dbCosTheta;
				eigenvectors_M(k, nCol) = eigenvectors_M(k, nCol) * dbCosTheta - dbMax * dbSinTheta;
			}
		}
		// Output eigenvalues
		for (int i = 0; i < dim; i++)
		{
			L(i) = matrix_M(i, i);
		}
		std::cout << "==================No sort eigenvectors====================" << std::endl;
		for (int i = 0; i < matrix_M.getRows(); i++)
		{
			for (int j = 0; j < matrix_M.getCols(); j++)
			{
				std::cout << matrix_M(i, j) << "\t";
			}
			std::cout << std::endl;
		}

		std::cout << "==================No sort eigenvalues=======================" << std::endl;

		for (int i = 0; i < eigenvectors_M.getRows(); i++)
		{
			for (int j = 0; j < eigenvectors_M.getCols(); j++)
			{
				std::cout << eigenvectors_M(i, j) << "\t";
			}
			std::cout << std::endl;
		}
	}
	//void Eig(Mat3_t & M, Vec3_t & eigenvalues, Mat3_t * Q = NULL, bool Qtrans = false)
	static void Eig(Mat3_t & M, Vec3_t & L, Vec3_t & V0, Vec3_t & V1, Vec3_t & V2, bool SortAsc = false, bool SortDesc = false)
	{
		// precision requirements
		double precision = 1e-10;
		// max number of iterations
		int max = 999;

		Mat3_t matrix_M = M;
		//current iteration
		int dim = matrix_M.getCols();
		// Creat matrix E
		Mat3_t eigenvectors_M(dim, 1);

		Mat3_t eigenvalues_M(dim, 1);


		int nCount = 0;		//current iteration
		while (1)
		{
			//find the largest element on the off-diagonal line of the matrix
			double dbMax = matrix_M(0, 1);
			int nRow = 0;
			int nCol = 1;
			for (int i = 0; i < dim; i++) {			//row
				for (int j = 0; j < dim; j++) {		//column
					double d = fabs(matrix_M(i, j));
					if ((i != j) && (d > dbMax))
					{
						dbMax = d;
						nRow = i;
						nCol = j;
					}
				}
			}

			if (dbMax < precision)     //precision check 
				break;
			if (nCount > max)       //iterations check
				break;
			nCount++;

			double dbAij = matrix_M(nRow, nCol);//ij
			double dbAjj = matrix_M(nCol, nCol);//jj
			double dbAii = matrix_M(nRow, nRow);//ii
												//compute rotate angle
			double dbAngle = 0.5*atan2(-2 * dbAij, (-dbAii + dbAjj));
			double dbSinTheta = sin(dbAngle);
			double dbCosTheta = cos(dbAngle);
			double dbSin2Theta = sin(2 * dbAngle);
			double dbCos2Theta = cos(2 * dbAngle);

			matrix_M(nRow, nRow) = dbAii * dbCosTheta*dbCosTheta +
				dbAjj * dbSinTheta*dbSinTheta + dbAij * dbSin2Theta;
			matrix_M(nCol, nCol) = dbAii * dbSinTheta*dbSinTheta +
				dbAjj * dbCosTheta*dbCosTheta - dbAij * dbSin2Theta;
			matrix_M(nCol, nRow) = 0.5*(dbAjj - dbAii)*dbSin2Theta + dbAij * dbCos2Theta;
			matrix_M(nRow, nCol) = matrix_M(nCol, nRow);

			for (int k = 0; k < dim; k++)
			{
				if ((k != nCol) && (k != nRow))
				{
					dbMax = matrix_M(k, nRow);
					matrix_M(k, nRow) = matrix_M(k, nCol) * dbSinTheta + dbMax * dbCosTheta;
					matrix_M(k, nCol) = matrix_M(k, nCol) * dbCosTheta - dbMax * dbSinTheta;

					matrix_M(nRow, k) = matrix_M(k, nRow);
					matrix_M(nCol, k) = matrix_M(k, nCol);
				}
			}

			//for (int k = 0; k < dim; k++)
			//{
			//	if ((k != nCol) && (k != nRow))
			//	{
			//		dbMax = matrix_M(nRow, k);
			//		matrix_M(nRow, k) = matrix_M(nCol, k) * dbSinTheta + dbMax * dbCosTheta;
			//		matrix_M(nCol, k) = matrix_M(nCol, k) * dbCosTheta - dbMax * dbSinTheta;
			//	}
			//}

			// compute eigenvector
			for (int k = 0; k < dim; k++)
			{
				dbMax = eigenvectors_M(k, nRow);
				eigenvectors_M(k, nRow) = eigenvectors_M(k, nCol) * dbSinTheta + dbMax * dbCosTheta;
				eigenvectors_M(k, nCol) = eigenvectors_M(k, nCol) * dbCosTheta - dbMax * dbSinTheta;
			}
		}


		// Output eigenvalues
		for (int i = 0; i < dim; i++)
		{
			L(i) = matrix_M(i, i);
		}

		// sort
		if (SortAsc)
		{
			std::map<double, int> mapEigen;
			for (int i = 0; i < dim; i++)
			{
				L(i) = matrix_M(i, i);
				mapEigen.insert(std::make_pair(L(i), i));
			}

			double *pdbTmpVec = new double[dim*dim];
			std::map<double, int>::reverse_iterator iter = mapEigen.rbegin();
			for (int j = 0; iter != mapEigen.rend(), j < dim; ++iter, ++j) {
				for (int i = 0; i < dim; i++) {
					pdbTmpVec[i*dim + j] = eigenvectors_M(i, iter->second);
				}
				L(j) = iter->first;
			}
			for (int i = 0; i < dim; i++)
			{
				double dSumVec = 0;
				for (int j = 0; j < dim; j++)
					dSumVec += pdbTmpVec[j * dim + i];
				if (dSumVec < 0)
				{
					for (int j = 0; j < dim; j++)
						pdbTmpVec[j * dim + i] *= -1;
				}
			}
			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					eigenvectors_M(i, j) = pdbTmpVec[i * dim + j];
				}
			}
			delete[]pdbTmpVec;
		}

		// Output eigenvalues
		for (int i = 0; i < dim; i++)
		{
			L(i) = matrix_M(i, i);
			V0(i) = eigenvectors_M(i, 0);
			V1(i) = eigenvectors_M(i, 1);
			V2(i) = eigenvectors_M(i, 2);
		}
		//¾ØÕóÊä³ö

		//for (int i = 0; i < matrix_M.getRows(); i++)
		//{
		//	for (int j = 0; j < matrix_M.getCols(); j++)
		//	{
		//		std::cout << matrix_M(i, j) << "\t";
		//	}
		//	std::cout << std::endl;
		//}

		//std::cout << "=========================================" << std::endl;

		//for (int i = 0; i < eigenvectors_M.getRows(); i++)
		//{
		//	for (int j = 0; j < eigenvectors_M.getCols(); j++)
		//	{
		//		std::cout << eigenvectors_M(i, j) << "\t";
		//	}
		//	std::cout << std::endl;
		//}
	}

};

// Constants
namespace OrthoSys
{
	const Vec3_t O(0.0);        ///< Origin
	const Vec3_t e0(1.0, 0.0, 0.0), e1(0.0, 1.0, 0.0), e2(0.0, 0.0, 1.0); ///< Basis
	const Mat3_t I(3, 1.0);        ///< Identity
};

