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

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <set>
#include "Misha/PreProcess.h"
#include "Misha/Exceptions.h"
#include "Misha/Geometry.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Misha/Miscellany.h"
#include "Misha/SimplexRefinableMesh.h"
#include "Misha/Meshes.h"
#include "Misha/MGSolver.h"
#include "Misha/AutoDiff/AutoDiff.h"

static const unsigned int Dim = 3;

Misha::CmdLineParameter< std::string > Input( "in" );
Misha::CmdLineParameter< int > Degree( "degree" , 2 ) , CoarseNodeDimension( "coarseDim" , 1 ) , VCycles( "vCycles" , 3 ) , GSIterations( "gsIters" , 5 );
Misha::CmdLineParameter< double > PlanarityEpsilon( "pEps" , 0 );
Misha::CmdLineReadable Verbose( "verbose" ) , Multigrid( "mg" ) , FullVerbose( "fullVerbose" ) , LinearPrecision( "linear" );
Misha::CmdLineReadable* params[] = { &Input , &Degree , &Verbose , &CoarseNodeDimension , &Multigrid , &VCycles , &GSIterations , &FullVerbose , &LinearPrecision , &PlanarityEpsilon , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , Input.name.c_str() );
	printf( "\t[--%s <fem degree>=%d]\n" , Degree.name.c_str() , Degree.value );
	printf( "\t[--%s <coarse node dimensions>=%d]\n" , CoarseNodeDimension.name.c_str() , CoarseNodeDimension.value );
	printf( "\t[--%s <number of v-cycles>=%d]\n" , VCycles.name.c_str() , VCycles.value );
	printf( "\t[--%s <number of relaxation iterations per level>=%d]\n" , GSIterations.name.c_str() , GSIterations.value );
	printf( "\t[--%s <planarity epsilon>=%g]\n" , PlanarityEpsilon.name.c_str() , PlanarityEpsilon.value );
	printf( "\t[--%s]\n" , Multigrid.name.c_str() );
	printf( "\t[--%s]\n" , LinearPrecision.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
	printf( "\t[--%s]\n" , FullVerbose.name.c_str() );
}

template< unsigned int Degree , typename TestFunction >
void ExecuteDirect
(
	const SimplexRefinableCellMesh< Dim , Degree > &simplexRefinableCellMesh ,
	const std::function< Point< double , Dim > ( typename SimplexMesh< Dim , Degree >::NodeMultiIndex ) > &NodePosition ,
	const TestFunction &franke ,
	const std::function< bool ( typename SimplexMesh< Dim , Degree >::NodeMultiIndex ) > &IsBoundaryNode
)
{
	typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	Eigen::SparseMatrix< double > M , S;

	auto Franke = [&]( Point< double , Dim > p )
	{
		AutoDiff::Tensor< Dim > _p;
		for( unsigned int d=0 ; d<Dim ; d++ ) _p[d] = p[d];
		return (double)franke( _p );
	};

	auto franke_laplacian = franke.laplacian();
	auto FrankeLaplacian = [&]( Point< double , Dim > p )
	{
		AutoDiff::Tensor< Dim > _p;
		for( unsigned int d=0 ; d<Dim ; d++ ) _p[d] = p[d];
		return (double)franke_laplacian( _p );
	};

	// Get the system matrices
	{
		Miscellany::NestedTimer timer( "Got system matrices" , Verbose.set );
		M = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().mass() * simplexRefinableCellMesh.P();
		S = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().stiffness() * simplexRefinableCellMesh.P();
	}

	// Get the list of node multi-indices
	std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes() );
	for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap() ) nodeMultiIndices[idx] = nmi;

	// Identify the nodes that are locked
	std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes() );
	for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );

	Eigen::VectorXd b( simplexRefinableCellMesh.nodes() );

	// Initiailize the constraints to the evaluation of the Laplacian of the Franke function
	for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) b[i] = - FrankeLaplacian( NodePosition( nodeMultiIndices[i] ) );

	// Compute the dual representation of the constraints
	b = M * b;

	// Set the constraints at the locked vertices to the evluation of the Franke function
	for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) if( lockedNodes[i] ) b[i] = Franke( NodePosition( nodeMultiIndices[i] ) );

	// Adjust the right-hand-side to account for the locked nodes
	for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
		if( !lockedNodes[ iter.row() ] && lockedNodes[ iter.col() ] ) b[ iter.row() ] -= b[ iter.col() ] * iter.value();

	// Adjust the system matrix to account for the locked nodes
	for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
		if( lockedNodes[ iter.row() ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
		else if( lockedNodes[ iter.col() ] ) iter.valueRef() = 0;

	SparseSolver::LLT solver;

	// Compute the Cholesky factorization
	{
		Miscellany::NestedTimer timer( "Performed Cholesky (LLt) factorization" , Verbose.set );
		solver.compute( S );
		switch( solver.info() )
		{
			case Eigen::NumericalIssue: ERROR_OUT( "SparseLLT failed to factorize matrix -- numerical issue" );
			case Eigen::NoConvergence:  ERROR_OUT( "SparseLLT failed to factorize matrix -- no convergence" );
			case Eigen::InvalidInput:   ERROR_OUT( "SparseLLT failed to factorize matrix -- invalid input" );
			case Eigen::Success: ;
		}
		if( Verbose.set ) std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;
	}
	Eigen::VectorXd x;
	{
		Miscellany::NestedTimer timer( "Solved direct system" , Verbose.set );
		x = solver.solve( b );
		if( Verbose.set )
		{
			double eIn = 0 , eOut = 0;
			Eigen::VectorXd r = b - S* x;
			for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
			std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;
		}
	}

	double rms = 0;
	for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
	{
		Miscellany::StreamFloatPrecision sfp( std::cout , 3 , true );
		std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes() ) << std::endl;
	}
}

template< unsigned int Degree , typename RelaxerType , typename TestFunction >
void ExecuteMG
(
	const HierarchicalSimplexRefinableCellMesh< Dim , Degree > &simplexRefinableCellMesh ,
	const std::function< Point< double , Dim > ( typename SimplexMesh< Dim , Degree >::NodeMultiIndex ) > &NodePosition ,
	const TestFunction &franke ,
	const std::function< bool ( typename SimplexMesh< Dim , Degree >::NodeMultiIndex ) > &IsBoundaryNode ,
	unsigned int vCycles , unsigned int gsIters
)
{
	typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	Eigen::SparseMatrix< double > M , S;
	std::vector< Eigen::SparseMatrix< double > > P( CoarseNodeDimension.value );

	auto Franke = [&]( Point< double , Dim > p )
	{
		AutoDiff::Tensor< Dim > _p;
		for( unsigned int d=0 ; d<Dim ; d++ ) _p[d] = p[d];
		return (double)franke( _p );
	};

	auto franke_laplacian = franke.laplacian();
	auto FrankeLaplacian = [&]( Point< double , Dim > p )
	{
		AutoDiff::Tensor< Dim > _p;
		for( unsigned int d=0 ; d<Dim ; d++ ) _p[d] = p[d];
		return (double)franke_laplacian( _p );
	};

	// Get the system matrices
	{
		Miscellany::NestedTimer timer( "Got system matrices" , Verbose.set );
		for( unsigned int d=0 ; d<(unsigned int)CoarseNodeDimension.value ; d++ ) P[d] = simplexRefinableCellMesh.P( d+1 , d );
		{
			Eigen::SparseMatrix< double > _P , _Pt;
			_P = simplexRefinableCellMesh.P( simplexRefinableCellMesh.maxLevel() , CoarseNodeDimension.value );
			_Pt = _P.transpose();
			M = _Pt * simplexRefinableCellMesh.simplexMesh().mass() * _P;
			S = _Pt * simplexRefinableCellMesh.simplexMesh().stiffness() * _P;
		}
	}


	// Get the list of node multi-indices
	std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) );
	for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap( CoarseNodeDimension.value ) ) nodeMultiIndices[idx] = nmi;

	// Identify the nodes that are locked
	std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) );
	for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );

	Eigen::VectorXd b( simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) );

	// Initiailize the constraints to the evaluation of the Laplacian of the Franke function
	for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) ; i++ ) b[i] = - FrankeLaplacian( NodePosition( nodeMultiIndices[i] ) );

	// Compute the dual representation of the constraints
	b = M * b;

	// Set the constraints at the locked vertices to the evluation of the Franke function
	for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) if( lockedNodes[i] ) b[i] = Franke( NodePosition( nodeMultiIndices[i] ) );

	// Adjust the right-hand-side to account for the locked nodes
	for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
		if( !lockedNodes[ iter.row() ] && lockedNodes[ iter.col() ] ) b[ iter.row() ] -= b[ iter.col() ] * iter.value();

	// Adjust the system matrix to account for the locked nodes
	for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
		if( lockedNodes[ iter.row() ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
		else if( lockedNodes[ iter.col() ] ) iter.valueRef() = 0;

	// Compute the hierarchy of systems
	MGSolver::Solver< RelaxerType > *mgSolver = nullptr;
	{
		Miscellany::NestedTimer timer( "Constructed multigrid system" , Verbose.set );
		mgSolver = new MGSolver::Solver< RelaxerType >( S , P , FullVerbose.set );
	}

	if( Verbose.set && !FullVerbose.set ) std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;

	Eigen::VectorXd x;
	{
		Miscellany::NestedTimer timer( "Solved multigrid system" , Verbose.set );
		x = mgSolver->solve( b , vCycles , gsIters , gsIters , FullVerbose.set );
		if( Verbose.set )
		{
			double eIn = 0 , eOut = 0;
			Eigen::VectorXd r = b - S* x;
			for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
			std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;
		}
	}

	double rms = 0;
	for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
	{
		Miscellany::StreamFloatPrecision sfp( std::cout , 3 , true );
		std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) ) << std::endl;
	}

	delete mgSolver;
}

template< unsigned int Degree , typename TestFunction >
void Execute( const Meshes::PolyhedronMesh< unsigned int > &polyMesh , const std::vector< Point< double , Dim > > &vertices , const TestFunction &franke )
{
	typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;

	std::set< unsigned int > boundaryVertices;
	{
		std::vector< unsigned int > faceCount( polyMesh.polygons() , 0 );
		for( unsigned int i=0 ; i<polyMesh.polyhedra() ; i++ ) for( unsigned int j=0 ; j<polyMesh.polyhedronFaces(i) ; j++ ) 
			faceCount[ polyMesh.polyhedronFace(i,j).first ]++;
		for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ ) if( faceCount[i]==1 )
		{
			Meshes::Polygon< unsigned int > polygon = polyMesh.polygon( i );
			for( unsigned int j=0 ; j<polygon.size() ; j++ ) boundaryVertices.insert( polygon[j] );
			boundaryVertices.insert( (unsigned int)vertices.size() + i );
		}
	}

	auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
	{
		for( unsigned int d=0 ; d<Degree ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
		return true;
	};

	SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
	eWeights.kWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1;

	std::function< Point< double , Dim > ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
	{
		if( idx<vertices.size() ) return vertices[idx];
		ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
		return Point< double , Dim >();
	};
	std::function< Point< double , Dim > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

	auto NodePosition = [&]( NodeMultiIndex nmi )
	{
		Point< double , Dim > p;
		for( unsigned int d=0 ; d<Degree ; d++ ) p += fullVertexPositionFunction( nmi[d] );
		return p/(double)Degree;
	};

	if( Multigrid.set )
	{
		HierarchicalSimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh;
		simplexRefinableCellMesh = polyMesh.template hierarchicalSimplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , CoarseNodeDimension.value , false , LinearPrecision.set , PlanarityEpsilon.value , Verbose.set );
		return ExecuteMG< Degree , MGSolver::ParallelGaussSeidelRelaxer< 20 > >( simplexRefinableCellMesh , NodePosition , franke , IsBoundaryNode , VCycles.value , GSIterations.value );
	}
	else
	{
		SimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh;
		simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , CoarseNodeDimension.value , false , LinearPrecision.set , PlanarityEpsilon.value , Verbose.set );
		ExecuteDirect< Degree >( simplexRefinableCellMesh , NodePosition , franke , IsBoundaryNode );
	}
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );

	if( !Input.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	Verbose.set |= FullVerbose.set;

	std::stringstream sStream;
	sStream << "Running Time [Franke Test 3D (v. " << VERSION << ")]";
	Miscellany::NestedTimer timer( sStream.str() , Verbose.set );

	std::vector< std::vector< std::pair< unsigned int , bool > > > polyhedra;
	std::vector< Meshes::Polygon< unsigned int > > polygons;
	std::vector< Point< double , Dim > > vertices;

	Meshes::PolyhedronMesh< unsigned int >::ReadOVM( Input.value , polyhedra , polygons , vertices );
	Meshes::PolyhedronMesh< unsigned int > polyMesh = Meshes::PolyhedronMesh< unsigned int >( polyhedra , polygons );

	{
		using namespace AutoDiff;
		auto x = Linear< UIntPack<> , UIntPack< Dim > >( {} , {0} );
		auto y = Linear< UIntPack<> , UIntPack< Dim > >( {} , {1} );
		auto z = Linear< UIntPack<> , UIntPack< Dim > >( {} , {2} );

		auto cx2 = Pow( 9.*x-2. , 2. );
		auto cy2 = Pow( 9.*y-2. , 2. );
		auto cz2 = Pow( 9.*z-2. , 2. );

		auto cx1 = Pow( 9.*x+1. , 2. );
		auto cx7 = Pow( 9.*x+7. , 2. );

		auto cy3 = Pow( 9.*y-3. , 2. );
		auto cx4 = Pow( 9.*x-7. , 2. );

		auto cy7 = Pow( 9.*y-7. , 2. );
		auto cz5 = Pow( 9.*z-5. , 2. );

		auto e1 = (3./4.) * Exp( -(1./4.)*( cx2 + cy2 + cz2 ) );
		auto e2 = (3./4.) * Exp( -(1./49.)*cx1 - (9./10.)*y - (1./10.) - (9./10.)*z - (1./10.) );
		auto e3 = (1./2.) * Exp( -(1./4.)*( cx7 + cy3 + cz5 ) );
		auto e4 = (1./5.) * Exp( -( cx4 + cy7 + cz5 ) );

		auto franke = e1 + e2 + e3 - e4;

		switch( Degree.value )
		{
			case 1: Execute< 1 >( polyMesh , vertices , franke ) ; break;
			case 2: Execute< 2 >( polyMesh , vertices , franke ) ; break;
			case 3: Execute< 3 >( polyMesh , vertices , franke ) ; break;
			default: ERROR_OUT( "Only degrees 1, 2, or 3 supported" );
		}
	}
}