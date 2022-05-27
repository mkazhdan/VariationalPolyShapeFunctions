#ifndef PRE_PROCESS_INCLUDED
#define PRE_PROCESS_INCLUDED

#define VERSION 2.5
#undef USE_UNORDERED_SET_MAP		// Surprisingly, this seems to make things slower

#undef EIGEN_USE_MKL_ALL			// Enabling this lets Eigen use the MKL library for efficiency
									// [NOTE] This requires setting the include and library directories accordingly

#ifdef EIGEN_USE_MKL_ALL
// Link to the necessary libraries 
// [NOTE] This may only work on Windows
#pragma comment( lib, "mkl_intel_lp64" )
#pragma comment( lib, "mkl_intel_thread" )
#pragma comment( lib, "mkl_core" )
#pragma comment( lib, "libiomp5md" )

#include "Eigen/PardisoSupport"
#endif // EIGEN_USE_MKL_ALL

#include "Eigen/Sparse"

struct SparseSolver
{
#ifdef EIGEN_USE_MKL_ALL
	typedef Eigen::PardisoLLT< Eigen::SparseMatrix< double > > LLT;
	typedef Eigen::PardisoLDLT< Eigen::SparseMatrix< double > > LDLT;
#else // !EIGEN_USE_MKL_ALL
	typedef Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > LLT;
	typedef Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > LDLT;
#endif // EIGEN_USE_MKL_ALL
};

#endif // PRE_PROCESS_INCLUDED
