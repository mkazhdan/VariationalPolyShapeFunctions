/*
Copyright (c) 2021
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

#ifndef LINEARLY_CONSTRAINED_QUADRATIC_OPTIMIZATION_INCLUDED
#define LINEARLY_CONSTRAINED_QUADRATIC_OPTIMIZATION_INCLUDED

#include <Eigen/Sparse>
#include <Eigen/Dense>

struct LCQO
{
	// Solve for the minimizer of the energy:
	//		E(x) = x^t * Q * x + 2 * x^t * q ( + c)
	// Subject to the linear constraint:
	//		C * x  = c
	// To solve for this, we need to find the value of x s.t.:
	// 	   \nabla E|_x 2 * Q * x + 2 * q \in ColSpan( c )
	// and:
	//		C * x = c
	// This can be formulated as solving for the vectors x and l s.t.:
	// 	   | Q C^t | | x | = | -q |
	// 	   | C  0  | | l | = |  c |
	// 

	// Running with reduce=true has the system try to remove as many (recursively) one-dimensional constraints as possible
	// so as to reduce the total size of the system being solved.
	LCQO( void ) : _variableCount(0) , _constraintCount(0){}
	LCQO( const Eigen::SparseMatrix< double > &Q , const Eigen::VectorXd &q , const Eigen::SparseMatrix< double , Eigen::RowMajor > &C , const Eigen::VectorXd &c , bool reduce=true );
	LCQO( const Eigen::SparseMatrix< double > &Q ,                            const Eigen::SparseMatrix< double , Eigen::RowMajor > &C , const Eigen::VectorXd &c , bool reduce=true );
	LCQO( const Eigen::SparseMatrix< double > &Q , const Eigen::VectorXd &q );
	LCQO( const Eigen::SparseMatrix< double > &Q                            );

	const Eigen::SparseMatrix< double > &reducedSystem( void ) const { return _M; }
	const Eigen::VectorXd &reducedConstraint( void ) const { return _b; }

	template< bool StableSolve >
	Eigen::VectorXd solve( void ) const;

	unsigned int reducedVariableCount( void ) const { return _variableCount; }
	unsigned int reducedConstraintCount( void ) const { return _constraintCount; }

	template< bool StableSolve >
	static Eigen::VectorXd Solve( const Eigen::SparseMatrix< double > &Q , const Eigen::SparseMatrix< double , Eigen::RowMajor > &C , const Eigen::VectorXd &c , bool reduce=true ){ return LCQO( Q , C , c , reduce ).solve< StableSolve >(); }

	template< bool StableSolve >
	static Eigen::VectorXd Solve( const Eigen::SparseMatrix< double > &Q , const Eigen::VectorXd &q , const Eigen::SparseMatrix< double , Eigen::RowMajor > &C , const Eigen::VectorXd &c , bool reduce=true ){ return LCQO( Q , q , C , c , reduce ).solve< StableSolve >(); }

	template< bool StableSolve >
	static Eigen::VectorXd Solve( const Eigen::SparseMatrix< double > &Q ){ return LCQO( Q ).solve< StableSolve >(); }

	template< bool StableSolve >
	static Eigen::VectorXd Solve( const Eigen::SparseMatrix< double > &Q , const Eigen::VectorXd &q ){ return LCQO( Q , q ).solve< StableSolve >(); }

protected:
	struct _VariableInfo
	{
		unsigned int idx;
		double lockedValue;

		_VariableInfo( unsigned int i=0 , double l=0 ) : idx(i) , lockedValue(l){}
	};

	std::vector< _VariableInfo > _variableInfo;
	unsigned int _variableCount , _constraintCount;

	Eigen::SparseMatrix< double > _M;
	Eigen::VectorXd _b;
};



LCQO::LCQO( const Eigen::SparseMatrix< double > &Q , const Eigen::SparseMatrix< double , Eigen::RowMajor > &C , const Eigen::VectorXd &c , bool reduce )
	: LCQO( Q , Eigen::VectorXd::Zero( Q.rows() ) , C , c , reduce ){}

LCQO::LCQO( const Eigen::SparseMatrix< double > &Q ) : LCQO( Q , Eigen::VectorXd::Zero( Q.rows() ) ){}

LCQO::LCQO( const Eigen::SparseMatrix< double > &Q , const Eigen::VectorXd &q , const Eigen::SparseMatrix< double , Eigen::RowMajor > &C , const Eigen::VectorXd &c , bool reduce )
{
	typedef Eigen::SparseMatrix< double > QMatrix;
	typedef Eigen::SparseMatrix< double , Eigen::RowMajor > CMatrix;
	if( Q.rows()!=Q.cols() ) THROW( "Quadratic energy term is not square: " , Q.rows() , " != " , Q.cols() );
	if( q.size()!=Q.rows() ) THROW( "Linear energy dimension does not match quadratic energy dimension: " , q.size() , " != " , Q.rows() );
	if( C.cols()!=Q.cols() ) THROW( "Constraint and energy dimensions do not match: " , C.cols() , " != " , Q.cols() );
	if( c.size()!=C.rows() ) THROW( "Constraint value dimensions do not match constraint term dimensions: " , c.size() , " != " , C.rows() );

	if( reduce )
	{
		_variableInfo.resize( Q.rows() );
		std::vector< unsigned int > constraintIndices( C.rows() , 0 );
		_variableCount = _constraintCount = 0;

		// First repeatedly identify the variables whose values should be locked.
		while( true )
		{
			bool foundLockedVariable = false;
			// Iterate over the rows of the constraint matrix
			for( unsigned int r=0 ; r<C.outerSize() ; r++ )
			{
				// Count the number of unlocked variables in the constraint row
				unsigned int count = 0;
				for( CMatrix::InnerIterator it( C , r ) ; it ; ++it ) if( _variableInfo[it.col()].idx!=-1 && it.value() ) count++;

				// If there is only one unlocked variable in the constraint row, it should be locked
				if( count==1 )
				{
					// First, get the value the variable should be locked to
					double value = c[r];
					for( CMatrix::InnerIterator it( C , r ) ; it ; ++it ) if( _variableInfo[it.col()].idx==-1 && it.value() )
						value -= it.value() * _variableInfo[ it.col() ].lockedValue;

					// Then, set the information for the locked variable
					for( CMatrix::InnerIterator it( C , r ) ; it ; ++it ) if( _variableInfo[it.col()].idx!=-1 && it.value() )
						_variableInfo[ it.col() ] = _VariableInfo( -1 , value/it.value() );

					foundLockedVariable = true;
				}
			}
			if( !foundLockedVariable ) break;
		}

		// Next identify those constraints that are redundant
		for( unsigned int r=0 ; r<C.outerSize() ; r++ )
		{
			unsigned int count = 0;
			double value = c[r];
			for( CMatrix::InnerIterator it( C , r ) ; it ; ++it ) if( _variableInfo[ it.col() ].idx==-1 )
				value -= it.value() * _variableInfo[ it.col() ].lockedValue;
			else count++;
			if( !count )
			{
				constraintIndices[r] = -1;
				if( fabs( value )>1e-8 ) THROW( "Inconsistent constraints: " , value );
			}
		}

		for( unsigned int i=0 ; i<_variableInfo.size() ; i++ ) if( _variableInfo[i].idx!=-1 ) _variableInfo[i].idx = _variableCount++;
		for( unsigned int i=0 ; i<constraintIndices.size() ; i++ ) if( constraintIndices[i]!=-1 ) constraintIndices[i] = _constraintCount++;

		_M.resize( _variableCount + _constraintCount , _variableCount + _constraintCount );
		_b.resize( _variableCount + _constraintCount );
		std::vector< Eigen::Triplet< double > > entries;

		entries.reserve( Q.nonZeros() + 2*C.nonZeros() );

		for( unsigned int i=0 ; i<_variableInfo.size() ; i++ ) if( _variableInfo[i].idx!=-1 ) _b[ _variableInfo[i].idx ] = -q[i];
		for( unsigned int i=0 ; i<constraintIndices.size() ; i++ ) if( constraintIndices[i]!=-1 ) _b[ constraintIndices[i]+_variableCount ] = c[i];

		// Add the coefficients from the quadratic energy (and update the RHS for locked variables);
		for( unsigned int o=0 ; o<Q.outerSize() ; o++ ) for( QMatrix::InnerIterator it( Q , o ) ; it ; ++it )
		{
			if( _variableInfo[ it.col() ].idx!=-1 && _variableInfo[ it.row() ].idx!=-1 )
				entries.push_back( Eigen::Triplet< double >( _variableInfo[ it.row() ].idx , _variableInfo[ it.col() ].idx , it.value() ) );
			else if( _variableInfo[ it.row() ].idx!=-1 )
				_b[ _variableInfo[ it.row() ].idx ] -= it.value() * _variableInfo[ it.col() ].lockedValue;
		}

		// Add the coefficients from the linear constraints (and update the RHS for locked variables);
		for( unsigned int r=0 ; r<C.outerSize() ; r++ ) if( constraintIndices[r]!=-1 ) 
			for( CMatrix::InnerIterator it( C , r ) ; it ; ++it )
			{
				if( _variableInfo[ it.col() ].idx!=-1 )
				{
					// Add C below
					entries.push_back( Eigen::Triplet< double >( constraintIndices[ it.row() ] + _variableCount , _variableInfo[ it.col() ].idx , it.value() ) );
					// Add C^t to the right
					entries.push_back( Eigen::Triplet< double >( _variableInfo[ it.col() ].idx , constraintIndices[ it.row() ] + _variableCount , it.value() ) );
				}
				else _b[ constraintIndices[ it.row() ]+_variableCount ] -= it.value() * _variableInfo[ it.col() ].lockedValue;
			}

		_M.setFromTriplets( entries.begin() , entries.end() );
	}
	else
	{
		_variableCount = (unsigned int)Q.rows();
		_constraintCount = (unsigned int)C.rows();

		_variableInfo.resize( Q.rows() );
		for( unsigned int i=0 ; i<_variableInfo.size() ; i++ ) _variableInfo[i].idx = i;

		_M.resize( Q.rows() + C.rows() , Q.cols() + C.rows() );
		_b.resize( Q.rows() + C.rows() );
		std::vector< Eigen::Triplet< double > > entries;

		entries.reserve( Q.nonZeros() + 2*C.nonZeros() );

		for( unsigned int o=0 ; o<Q.outerSize() ; o++ ) for( QMatrix::InnerIterator it( Q , o ) ; it ; ++it )
			entries.push_back( Eigen::Triplet< double >( (int)it.row() , (int)it.col() , it.value() ) );
		for( unsigned int o=0 ; o<C.outerSize() ; o++ ) for( CMatrix::InnerIterator it( C , o ) ; it ; ++it )
		{
			// Add C below
			entries.push_back( Eigen::Triplet< double >( (int)it.row() + (int)Q.rows() , (int)it.col() , it.value() ) );
			// Add C^t to the right
			entries.push_back( Eigen::Triplet< double >( (int)it.col() , (int)it.row() + (int)Q.cols() , it.value() ) );
		}
		for( unsigned int i=0 ; i<q.size() ; i++ ) _b[ i ]            = -q[i];
		for( unsigned int i=0 ; i<c.size() ; i++ ) _b[ i + Q.rows() ] =  c[i];

		_M.setFromTriplets( entries.begin() , entries.end() );
	}
}

LCQO::LCQO( const Eigen::SparseMatrix< double > &Q , const Eigen::VectorXd &q )
{
	typedef Eigen::SparseMatrix< double > QMatrix;
	typedef Eigen::SparseMatrix< double , Eigen::RowMajor > CMatrix;
	if( Q.rows()!=Q.cols() ) THROW( "Quadratic energy term is not square: " , Q.rows() , " != " , Q.cols() );
	if( q.size()!=Q.rows() ) THROW( "Linear energy dimension does not match quadratic energy dimension: " , q.size() , " != " , Q.rows() );

	_variableCount = (unsigned int)Q.rows();
	_constraintCount = 0;

	_variableInfo.resize( Q.rows() );
	for( unsigned int i=0 ; i<_variableInfo.size() ; i++ ) _variableInfo[i].idx = i;

	_M = Q;
	_b = -q;
}

template< bool StableSolve >
Eigen::VectorXd LCQO::solve( void ) const
{
	if( _M.rows()==0 || _M.cols()==0 )
	{
		Eigen::VectorXd solution( _variableInfo.size() );
		for( unsigned int i=0 ; i<_variableInfo.size() ; i++ )
			if( _variableInfo[i].idx==-1 ) solution[i] = _variableInfo[i].lockedValue;
			else                           solution[i] = 0;
		return solution;
	}
	else
	{
		typedef typename std::conditional< StableSolve , Eigen::SparseQR< Eigen::SparseMatrix< double > , Eigen::COLAMDOrdering< int > > , Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > >::type Solver;

		Solver solver( _M );

		switch( solver.info() )
		{
			case Eigen::NumericalIssue: THROW( "Eigen::Solver failed to factorize matrix -- numerical issue" );
			case Eigen::NoConvergence:  THROW( "Eigen::Solver failed to factorize matrix -- no convergence" );
			case Eigen::InvalidInput:   THROW( "Eigen::Solver failed to factorize matrix -- invalid input" );
			case Eigen::Success: ;
		}

		Eigen::VectorXd x = solver.solve( _b );
	
		Eigen::VectorXd solution( _variableInfo.size() );
		for( unsigned int i=0 ; i<_variableInfo.size() ; i++ )
			if( _variableInfo[i].idx==-1 ) solution[i] = _variableInfo[i].lockedValue;
			else                           solution[i] = x[ _variableInfo[i].idx ];
		return solution;
	}
}

#endif // LINEARLY_CONSTRAINED_QUADRATIC_OPTIMIZATION_INCLUDED