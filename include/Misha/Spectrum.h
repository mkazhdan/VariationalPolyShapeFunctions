/*
Copyright (c) 2019
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef SPECTRUM_INCLUDED
#define SPECTRUM_INCLUDED

#include <Spectra/SymGEigsSolver.h>
#include <Eigen/Sparse>

template< class Real >
class EigenSolverCholeskyLDLt
{
	typedef Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > Eigen_Solver;
	typedef Eigen::VectorXd                                        Eigen_Vector;
	Eigen_Solver _solver;
	Eigen_Vector _eigenB;
	const Eigen::SparseMatrix< Real > &_M;
public:
	EigenSolverCholeskyLDLt( const Eigen::SparseMatrix< Real >& M ) : _M(M)
	{
		_solver.analyzePattern( _M );
		_solver.factorize( _M );
		if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCholeskyLDLt::EigenSolverCholeskyLDLt Failed to factorize matrix\n" ) , exit(0);
		_eigenB.resize( _M.rows() );
	}
	void multiply( const Real *x , Real *b )
	{
#pragma omp parallel for
		for( int i=0 ; i<_eigenB.size() ; i++ ) _eigenB[i] = x[i];
		Eigen_Vector eigenX = _M * _eigenB;
#pragma omp parallel for
		for( int i=0 ; i<eigenX.size() ; i++ ) b[i] = (Real)eigenX[i];
	}
	void solve( const Real *b , Real *x )
	{
#pragma omp parallel for
		for( int i=0 ; i<_eigenB.size() ; i++ ) _eigenB[i] = b[i];
		Eigen_Vector eigenX = _solver.solve( _eigenB );
#pragma omp parallel for
		for( int i=0 ; i<eigenX.size() ; i++ ) x[i] = (Real)eigenX[i];
	}
	size_t dimension( void ) const { return _eigenB.size(); }
};

template< typename Real >
struct Spectrum
{
	Spectrum( const Eigen::SparseMatrix< Real > &M , const Eigen::SparseMatrix< Real > &S , unsigned int dimension , Real offset );
	size_t size( void ) const { return _eigenvalues.size(); }
	const Real &eValue( unsigned int idx ) const { return _eigenvalues[idx]; }
	const std::vector< Real > &eVector( unsigned int idx ) const { return _eigenvectors[idx]; }
protected:
	static const unsigned long long _MAGIC_NUMBER;
	std::vector< Real > _eigenvalues;
	std::vector< std::vector< Real > > _eigenvectors;
};

//////////////
// Spectrum //
//////////////
template< typename Real > const unsigned long long Spectrum< Real >::_MAGIC_NUMBER = 0x2019ull;

template< typename Real >
Spectrum< Real >::Spectrum( const Eigen::SparseMatrix< Real > &M , const Eigen::SparseMatrix< Real > &S , unsigned int dimension , Real offset )
{
	// [Definition]
	//	We define the generalized eigensystem (A,B) to be the system A v = \lambda B v
	// [Fact 1]
	// If (v,\lambda) is a solution to the ges (A,B) then (v,\lambda+\epsilon) is a solution to the ges (A+\epsilon*B,B):
	//		(A + \epsilon B) v = (\lambda + \epsilon) B v
	//	<=>	A v + \epsilon B v = \lambda B v + \epsilon B V
	//	<=>                A v = \lambda B v
	// [Fact 2]
	// If (w,\delta) is a solution to the ges (A^{-1},B^{-1}) then (A^{-1}w,1/\delta) is a solution to the ges (A,B):
	//		A^{-1} w = \delta B^{-1} w
	//	<=> v = \delta B^{-1} A v
	//	<=> 1/\delta B v = A v
	// [Corollary]
	// If (w,\delta) is a solution to the ges ( (A+\epsilon*B)^{-1} , B^{-1} ) then:
	//	=> ( (A+\epsilon*B)^{-1} w , 1/\delta ) is a solution to the ges (A+\epsilon*B,B)
	//	=> ( (A+\epsilon*B)^{-1} w , 1\delta-\epsilon ) is a solution to the ges (A,B)

	typedef EigenSolverCholeskyLDLt< Real > Solver;
	struct InverseOperator
	{
		Solver solver;
		InverseOperator( const Eigen::SparseMatrix< Real > &M ) : solver( M ){}
		int rows( void ) const { return (int)solver.dimension(); }
		int cols( void ) const { return (int)solver.dimension(); }
		void perform_op( const Real *in , Real *out ) const { const_cast< Solver & >(solver).solve( in , out ); };
	};

	struct InverseBOperator
	{
		Solver solver;
		InverseBOperator( const Eigen::SparseMatrix< Real > &M ) : solver( M ){}
		int rows( void ) const { return (int)solver.dimension(); }
		int cols( void ) const { return (int)solver.dimension(); }
		void solve( const Real *in , Real *out ) const { const_cast< Solver & >(solver).multiply( in , out ); }
		void mat_prod( const Real *in , Real *out ) const { const_cast< Solver & >(solver).solve( in , out ); };
	};

	// Offset the stiffness matrix so that it becomes positive definite
	Eigen::SparseMatrix< double > _S = S + M * offset;

	InverseOperator op( _S );
	InverseBOperator Bop( M );

	Spectra::SymGEigsSolver< Real , Spectra::LARGEST_ALGE , InverseOperator , InverseBOperator , Spectra::GEIGS_REGULAR_INVERSE > geigs( &op , &Bop , dimension , 2*dimension );
	geigs.init();
	int nconv = geigs.compute();
	if( nconv!=dimension ) fprintf( stderr , "[WARNING] Number of converged is not equal to dimension: %d != %d\n" , nconv , dimension );
	Eigen::VectorXd evalues;
	Eigen::MatrixXd evecs;
	if( geigs.info()==Spectra::SUCCESSFUL )
	{
		evalues = geigs.eigenvalues();
		evecs = geigs.eigenvectors();
	}
	else if( geigs.info()==Spectra::NOT_COMPUTED )    fprintf( stderr , "[ERROR] Not computed\n" ) , exit(0);
	else if( geigs.info()==Spectra::NOT_CONVERGING 	) fprintf( stderr , "[ERROR] Not converging\n" ) , exit(0);
	else if( geigs.info()==Spectra::NUMERICAL_ISSUE ) fprintf( stderr , "[ERROR] Numerical issue\n" ) , exit(0);
	else                                              fprintf( stderr , "[ERROR] Failed\n" ) , exit(0);

	_eigenvalues.resize( evecs.cols() );
	_eigenvectors.resize( evecs.cols() );

	for( int i=0 ; i<evecs.cols() ; i++ )
	{
		_eigenvectors[i].resize( M.rows() );
		_eigenvalues[i] = (Real)1./evalues[i] - offset;
		std::vector< Real > w( M.rows() );
		for( int j=0 ; j<evecs.rows() ; j++ ) w[j] = evecs(j,i);
		op.perform_op( &w[0] , &_eigenvectors[i][0] );
#if 0
		Real l2Norm = 0;
#pragma omp parallel for reduction( + : l2Norm )
		for( int j=0 ; j<M.rows() ; j++ ) for( int k=0 ; k<M.rowSizes[j] ; k++ ) l2Norm += M[j][k].Value * _eigenvectors[i][j] * _eigenvectors[i][ M[j][k].N ];
		l2Norm = (Real)sqrt( l2Norm );
#pragma omp parallel for
		for( int j=0 ; j<_eigenvectors[i].size() ; j++ ) _eigenvectors[i][j] /= l2Norm;
#endif
	}
}
#endif // SPECTRUM_INCLUDED