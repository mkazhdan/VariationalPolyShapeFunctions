/*
Copyright (c) 2022, Michael Kazhdan
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

#ifndef DEFORMABLE_SOLID_INCLUDED
#define DEFORMABLE_SOLID_INCLUDED
#include "Eigen/Sparse"
#include "Meshes.h"
#include "Misha/MGSolver.h"

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
struct DeformableSolid
{
	struct MaterialProperties
	{
		double density , youngsModulus , poissonRatio;
		MaterialProperties( void ) : density(1.) , youngsModulus(1e10) , poissonRatio(0.25) {}
		void lameCoefficients( double &fWeight , double &tWeight ) const { fWeight = youngsModulus / ( 2. + 2 * poissonRatio ) , tWeight = youngsModulus * poissonRatio / ( 1. + poissonRatio ) / ( 1. - 2. * poissonRatio ) ;}
	};

	struct AdvanceStats
	{
		double updateTime , solveTime;
		AdvanceStats( void ) : updateTime(0) , solveTime(0){}
	};

	MaterialProperties materialProperties( void ) const { return _materialProperties; }
	void setMaterialProperties( const MaterialProperties &materialProperties ){ _materialProperties = materialProperties ; _systemNeedsUpdate=true; }

	double timeStep( void ) const { return _timeStep; }
	void setTimeStep( double timeStep ){ _timeStep = timeStep ; _systemNeedsUpdate=true; }

	template< typename Real >
	DeformableSolid( const SolidMesh &solidMesh , unsigned int vNum , std::function< Point< Real , Dim > ( unsigned int ) > vFunction , Point< double , Dim > gravity , bool verbose=false );
	~DeformableSolid( void );

	template< typename Real >
	void setGeometricState( const std::vector< Point< Real , Dim > > &positions , const std::vector< Point< Real , Dim > > &velocities ){ setGeometricState( [&]( unsigned int idx ){ return Point< double , Dim >( positions[idx] ); } , [&]( unsigned int idx ){ return Point< double , Dim >( velocities[idx] ); } ); }

	template< typename PositionFunctor , typename VelocityFunctor >
	void setGeometricState( const PositionFunctor &pFunctor , const VelocityFunctor &vFunctor );

	void advance( void );
	void advance( AdvanceStats &aStats );

	void lock( const std::vector< bool > &lockedVertices );

	unsigned int vertices( void ) const { return _vNum; }

	void setVertex( unsigned int vIdx , Point< double , Dim > v );
	Point< double , Dim > vertex( unsigned int vIdx ) const;
	bool lockedVertex( unsigned int vIdx ) const;

	const Eigen::VectorXd &x( void ) const { return _x; }

	const SolidMesh &solidMesh( void ) const { return _solidMesh; }

	double energy( void ) const;

	double gravity;
	unsigned int vCycles , gsIters;
	bool verbose;
protected:
	typedef typename std::conditional< Hierarchical , MGSolver::Solver< MGSolver::ParallelGaussSeidelRelaxer< 20u > > , Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > >::type SolverType;

	const SolidMesh &_solidMesh;
	double _timeStep , _fWeight , _tWeight;
	MaterialProperties _materialProperties;
	Eigen::SparseMatrix< double > __A , _A;
	Eigen::SparseMatrix< double > _M , _frobeniusS , _traceS;
	Eigen::SparseMatrix< double > _Pt , _P;
	std::vector< Eigen::SparseMatrix< double > > __P;
	Eigen::VectorXd _b , _g , _x , _v;
	SolverType *_solver;
	std::vector< bool > _lockedNodes;
	Eigen::VectorXd _lockedOffsets;
	unsigned int _vNum;
	bool _systemNeedsUpdate;

	void _updateSystem();
	void _updateConstraintOffset();
};

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::~DeformableSolid( void ){ delete _solver; }

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
template< typename Real >
DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::DeformableSolid( const SolidMesh &solidMesh , unsigned int vNum , std::function< Point< Real , Dim > ( unsigned int ) > vFunction , Point< double , Dim > g , bool verbose )
	: _timeStep( 1e-6 ) , gravity(1) , _solidMesh( solidMesh ) , _solver(NULL)
{
	Miscellany::NestedTimer timer( "Got deformable solid system matrices" , verbose );
//	if( verbose ) std::cout << "Getting deformable solid system matrices" << std::endl;
//	Miscellany::Timer timer;
	_systemNeedsUpdate = true;
	this->verbose = false;

	_vNum = vNum;
	_lockedNodes.resize( _solidMesh.nodes() , false );

//	timer.reset();

#ifdef NEW_SOLID_SYSTEM_MATRIX
	if constexpr( Hierarchical )
	{
		_P = _solidMesh.P( _solidMesh.maxLevel() , _solidMesh.maxLevel()-1 );
		_Pt = _P.transpose();
		__P.resize( _solidMesh.maxLevel()-1 );
		for( unsigned int d=0 ; d<__P.size() ; d++ ) __P[d] = _solidMesh.P( d+1 , d );
		Eigen::SparseMatrix< double > M , frobeniusS , traceS;
		_solidMesh.solidSimplexMesh().setMassFrobeniusStiffnessAndTraceStiffnessMatrices( M , frobeniusS , traceS );
		_M = _Pt * M * _P;
		_frobeniusS = _Pt * frobeniusS * _P;
		_traceS = _Pt * traceS * _P;
	}
	else _solidMesh.setMassFrobeniusStiffnessAndTraceStiffnessMatrices( _M , _frobeniusS , _traceS );
#else // !NEW_SOLID_SYSTEM_MATRIX
	_M = _solidMesh.massMatrix();
	_frobeniusS = _solidMesh.frobeniusStiffnessMatrix();
	_traceS = _solidMesh.traceStiffnessMatrix();
#endif // NEW_SOLID_SYSTEM_MATRIX
//	if( verbose ) std::cout << "Got deformable solid system matrices: " << timer.elapsed() << std::endl;
	if constexpr( Hierarchical )
	{
		_b = _Pt * _solidMesh.solidSimplexMesh().stiffnessVector();
		_g = _Pt * _solidMesh.solidSimplexMesh().dcVector( g );
	}
	else
	{
		_b = _solidMesh.stiffnessVector();
		_g = _solidMesh.dcVector( g );
	}

	_x.resize( _solidMesh.nodes() * Dim );
	_v.resize( _solidMesh.nodes() * Dim );
	_lockedOffsets = Eigen::VectorXd::Zero( _solidMesh.nodes() * Dim );

	setGeometricState( vFunction , []( unsigned int idx ){ return Point< double , Dim >(); } );
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
template< typename PositionFunctor , typename VelocityFunctor >
void DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::setGeometricState( const PositionFunctor &pFunctor , const VelocityFunctor &vFunctor )
{
	auto &nodeMap = _solidMesh.nodeMap();
	for( auto iter=nodeMap.cbegin() ; iter!=nodeMap.cend() ; iter++ )
	{
		Point< double , Dim > p , v;
		for( unsigned int d=0 ; d<Degree ; d++ ) p += pFunctor( iter->first[d] ) , v += vFunctor( iter->first[d] );
		p /= Degree , v /= Degree;
		for( unsigned int d=0 ; d<Dim ; d++ ) _x[ iter->second*Dim+d ] = p[d] , _v[ iter->second*Dim+d ] = v[d];
	}
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
void DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::_updateConstraintOffset( void )
{
	for( unsigned int i=0 ; i<_lockedOffsets.size() ; i++ ) _lockedOffsets[i] = 0;
	for( unsigned int i=0 ; i<_A.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(_A,i) ; iter ; ++iter )
		if( !_lockedNodes[ iter.row()/Dim ] && _lockedNodes[ iter.col()/Dim ] ) _lockedOffsets[ iter.row() ] -= _x[ iter.col() ] * iter.value();
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
void DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::_updateSystem( void )
{
	_materialProperties.lameCoefficients( _fWeight , _tWeight );
	__A = _A = _M + 2 * _timeStep * _timeStep * _materialProperties.density * ( _frobeniusS * _fWeight + _traceS * _tWeight );

	_updateConstraintOffset();

	for( unsigned int i=0 ; i<__A.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(__A,i) ; iter ; ++iter )
		if( _lockedNodes[ iter.row()/Dim ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
		else if( _lockedNodes[ iter.col()/Dim ] ) iter.valueRef() = 0;

	if constexpr( Hierarchical )
	{
		if( _solver ) delete _solver;
		_solver = new SolverType( __A , __P , false );
		_solver->state.vCycles = vCycles;
		_solver->state.pIters = gsIters;
		_solver->state.verbose = verbose;
	}
	else
	{
		if( !_solver )
		{
			_solver = new SolverType();
			_solver->analyzePattern( __A );
		}
		_solver->factorize( __A );

		switch( _solver->info() )
		{
			case Eigen::NumericalIssue: ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- numerical issue" );
			case Eigen::NoConvergence:  ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- no convergence" );
			case Eigen::InvalidInput:   ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- invalid input" );
			case Eigen::Success: ;
		}
	}

	_systemNeedsUpdate = false;
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
void DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::advance( void )
{
	AdvanceStats aStats;
	advance( aStats );
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
void DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::advance( AdvanceStats &aStats )
{
	if( _systemNeedsUpdate )
	{
		Miscellany::Timer timer;
		_updateSystem();
		aStats.updateTime = timer.elapsed();
	}
	Eigen::VectorXd b = _M * _x + _M * _v * _timeStep + _b * 2 * _timeStep * _timeStep * _materialProperties.density * ( _fWeight + _tWeight * Dim ) + _timeStep * _timeStep * _materialProperties.density * _g * gravity;

	b += _lockedOffsets;
	for( unsigned int i=0 ; i<_lockedNodes.size() ; i++ ) if( _lockedNodes[i] ) for( unsigned int d=0 ; d<Dim ; d++ ) b[i*Dim+d] = _x[i*Dim+d];
	Miscellany::Timer timer;
	Eigen::VectorXd x;
	if constexpr( Hierarchical ) x = _solver->solve( b , _x );
	else                         x = _solver->solve( b );
	aStats.solveTime = timer.elapsed();
	if( _timeStep ) _v = ( x - _x ) / _timeStep;
	_x = x;
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
Point< double , Dim > DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::vertex( unsigned int vIdx ) const
{
	unsigned int idx[Degree];
	for( unsigned int d=0 ; d<Degree ; d++ ) idx[d] = vIdx;
	typename SolidSimplexMesh< Dim , Degree >::NodeMultiIndex ni( idx );
	unsigned int i = _solidMesh.nodeIndex(ni);
	Point< double , Dim > p;
	for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = _x[ i*Dim+d];
	return p;
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
void DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::setVertex( unsigned int vIdx , Point< double , Dim > v )
{
	unsigned int idx[Degree];
	for( unsigned int d=0 ; d<Degree ; d++ ) idx[d] = vIdx;
	typename SolidSimplexMesh< Dim , Degree >::NodeMultiIndex ni( idx );
	unsigned int i = _solidMesh.nodeIndex(ni);
	for( unsigned int d=0 ; d<Dim ; d++ ) _x[ i*Dim+d] = v[d];
	_updateConstraintOffset();
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
bool DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::lockedVertex( unsigned int vIdx ) const
{
	unsigned int idx[Degree];
	for( unsigned int d=0 ; d<Degree ; d++ ) idx[d] = vIdx;
	return _lockedNodes[ _solidMesh.nodeIndex( typename SolidSimplexMesh< Dim , Degree >::NodeMultiIndex( idx ) ) ];
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
void DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::lock( const std::vector< bool > &lockedVertices )
{
	for( unsigned int i=0 ; i<_lockedNodes.size() ; i++ ) _lockedNodes[i] = false;
	const auto &nodeMap = _solidMesh.nodeMap();
	for( auto iter=nodeMap.cbegin() ; iter!=nodeMap.cend() ; iter++ )
	{
		bool locked = true;
		for( unsigned int d=0 ; d<Degree ; d++ ) if( iter->first[d]>=lockedVertices.size() || !lockedVertices[ iter->first[d] ] ) locked = false;
		_lockedNodes[ iter->second ] = locked;
	}

	_systemNeedsUpdate = true;
}

template< unsigned int Dim , unsigned int Degree , typename SolidMesh , bool Hierarchical >
double DeformableSolid< Dim , Degree , SolidMesh , Hierarchical >::energy( void ) const
{
	double e = 0 , volume = 0;
	Eigen::VectorXd x = _frobeniusS*_x*_fWeight + _traceS*_x*_tWeight - 2. * ( _fWeight + _tWeight * Dim ) * _b;
	for( unsigned int i=0 ; i<_x.size() ; i++ ) e += _x[i] * x[i];
	for( unsigned int i=0 ; i<x.size() ; i++ ) x[i] = 1;
	x = _M * x;
	for( unsigned int i=0 ; i<x.size() ; i++ ) volume += x[i];
	volume /= Dim;
	e += volume * Dim * ( _fWeight + _tWeight * Dim );
	return e;
}

#endif // DEFORMABLE_SOLID_INCLUDED