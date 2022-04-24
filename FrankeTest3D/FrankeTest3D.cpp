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
#include "Eigen/Sparse"
#include "Misha/Exceptions.h"
#include "Misha/Geometry.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Misha/Miscellany.h"
#include "Misha/SimplexRefinableMesh.h"
#include "Misha/Meshes.h"
#include "Misha/MGSolver.h"

static const unsigned int Dim = 3;

Misha::CmdLineParameter< std::string > Input( "in" );
Misha::CmdLineParameter< int > Degree( "degree" , 2 ) , CoarseNodeDimension( "coarseDim" , 1 ) , VCycles( "vCycles" , 1 ) , Iterations( "iters" , 0 ) , MGRelaxerType( "relaxer" , MGSolver::RELAXER_PARALLEL_GAUSS_SEIDEL );
Misha::CmdLineParameter< double > PlanarityEpsilon( "pEps" , 0 );
Misha::CmdLineReadable Verbose( "verbose" ) , Multigrid( "mg" ) , FullVerbose( "fullVerbose" ) , Debug( "debug" ) , LinearPrecision( "linear" );
Misha::CmdLineReadable* params[] = { &Input , &Degree , &Verbose , &CoarseNodeDimension , &Multigrid , &VCycles , &Iterations , &FullVerbose , &Debug , &MGRelaxerType , &LinearPrecision , &PlanarityEpsilon , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , Input.name.c_str() );
	printf( "\t[--%s <fem degree>=%d]\n" , Degree.name.c_str() , Degree.value );
	printf( "\t[--%s <coarse node dimensions>=%d]\n" , CoarseNodeDimension.name.c_str() , CoarseNodeDimension.value );
	printf( "\t[--%s <number of v-cycles>=%d]\n" , VCycles.name.c_str() , VCycles.value );
	printf( "\t[--%s <number of relaxation iterations per level>=%d]\n" , Iterations.name.c_str() , Iterations.value );
	printf( "\t[--%s <multigrid relaxer type>=%d]\n" , MGRelaxerType.name.c_str() , MGRelaxerType.value );
	for( unsigned int i=0 ; i<MGSolver::RELAXER_COUNT ; i++ ) std::cout << "\t\t" << i << "] " << MGSolver::RelaxerTypeNames[i] << std::endl;
	printf( "\t[--%s <planarity epsilon>=%g]\n" , PlanarityEpsilon.name.c_str() , PlanarityEpsilon.value );
	printf( "\t[--%s]\n" , Multigrid.name.c_str() );
	printf( "\t[--%s]\n" , LinearPrecision.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
	printf( "\t[--%s]\n" , FullVerbose.name.c_str() );
	printf( "\t[--%s]\n" , Debug.name.c_str() );
}

double Franke( Point< double , 2 > p )
{
	double x = p[0] , y = p[1];
	double cx2 = (9. * x - 2.) * (9. * x - 2.);
	double cy2 = (9. * y - 2.) * (9. * y - 2.);

	double cx1 = (9. * x + 1.) * (9. * x + 1.);
	double cx7 = (9. * x - 7.) * (9. * x - 7.);

	double cy3 = (9. * y - 3.) * (9. * y - 3.);
	double cx4 = (9. * x - 4.) * (9. * x - 4.);

	double cy7 = (9. * y - 7.) * (9. * y - 7.);

	return (3. / 4.) * exp(-(1. / 4.) * cx2 - (1. / 4.) * cy2) +
		(3. / 4.) * exp(-(1. / 49.) * cx1 - (9. / 10.) * y - 1. / 10.) +
		(1. / 2.) * exp(-(1. / 4.) * cx7 - (1. / 4.) * cy3) -
		(1. / 5.) * exp(-cx4 - cy7);
}

double FrankeLaplacian( Point< double , 2 > p )
{
	double x = p[0] , y = p[1];
	double mathematica =
		64.8 * exp(-pow(-4. + 9. * x, 2.0) - pow(-7. + 9. * y, 2.0)) -
		40.5 * exp(0.25 * (-pow(-7. + 9. * x, 2) - pow(-3. + 9. * y, 2))) -
		60.75 * exp(0.25 * (-pow(-2. + 9. * x, 2) - pow(-2. + 9. * y, 2))) -
		1.8720918367346937 * exp(-0.02040816326530612 * pow(1. + 9. * x, 2) -
			0.1 * (1. + 9. * y)) +
		10.125 * exp(0.25 * (-pow(-7. + 9. * x, 2) - pow(-3. + 9. * y, 2))) *
		pow(-7. + 9. * x, 2) -
		64.8 * exp(-pow(-4. + 9. * x, 2) - pow(-7. + 9. * y, 2)) *
		pow(-4. + 9. * x, 2) +
		15.1875 * exp(0.25 * (-pow(-2. + 9. * x, 2) - pow(-2. + 9. * y, 2))) *
		pow(-2. + 9. * x, 2) +
		0.1012078300708038 *
		exp(-0.02040816326530612 * pow(1. + 9. * x, 2) -
			0.1 * (1. + 9. * y)) *
		pow(1. + 9. * x, 2) -
		64.8 * exp(-pow(-4. + 9. * x, 2) - pow(-7. + 9. * y, 2)) *
		pow(-7. + 9. * y, 2) +
		10.125 * exp(0.25 * (-pow(-7. + 9. * x, 2) - pow(-3. + 9. * y, 2))) *
		pow(-3. + 9. * y, 2) +
		15.1875 * exp(0.25 * (-pow(-2. + 9. * x, 2) - pow(-2. + 9. * y, 2))) *
		pow(-2. + 9. * y, 2);

	return mathematica;
}

double Franke( Point< double , 3 > p )
{
	double x=p[0] , y=p[1] , z=p[2];

	double cx2 = (9. * x - 2.) * (9. * x - 2.);
	double cy2 = (9. * y - 2.) * (9. * y - 2.);
	double cz2 = (9. * z - 2.) * (9. * z - 2.);

	double cx1 = (9. * x + 1.) * (9. * x + 1.);
	double cx7 = (9. * x - 7.) * (9. * x - 7.);

	double cy3 = (9. * y - 3.) * (9. * y - 3.);
	double cx4 = (9. * x - 4.) * (9. * x - 4.);

	double cy7 = (9. * y - 7.) * (9. * y - 7.);
	double cz5 = (9. * z - 5.) * (9. * z - 5.);

	return (3. / 4.) *
		exp(-(1. / 4.) * cx2 - (1. / 4.) * cy2 - (1. / 4.) * cz2) +
		(3. / 4.) * exp(-(1. / 49.) * cx1 - (9. / 10.) * y - 1. / 10. -
			(9. / 10.) * z - 1. / 10.) +
		(1. / 2.) *
		exp(-(1. / 4.) * cx7 - (1. / 4.) * cy3 - (1. / 4.) * cz5) -
		(1. / 5.) * exp(-cx4 - cy7 - cz5);
}

double FrankeLaplacian( Point< double , 3 > p )
{
	double x=p[0] , y=p[1] , z=p[2];

	return (243. * (-2299. + 1800. * x * (2. + 9. * x))) /
		(480200. *
			exp((9. * (12. + 10. * x * (2. + 9. * x) + 49. * y + 49. * z)) /
				490.)) -
		(486. *
			exp(-pow(4. - 9. * x, 2) - pow(7. - 9. * y, 2) -
				pow(5. - 9. * z, 2)) *
			(59. + 6. * x * (-8. + 9. * x) + 6. * y * (-14. + 9. * y) +
				6. * z * (-10 + 9 * z))) /
		5. +
		(81. *
			exp((-pow(7. - 9. * x, 2) - 9 * pow(1. - 3. * y, 2) -
				pow(5. - 9. * z, 2)) /
				4.) *
			(77. + 9. * x * (-14. + 9. * x) + 27. * y * (-2. + 3. * y) +
				9 * z * (-10. + 9. * z))) /
		8. +
		(729. * (2. + 3. * x * (-4. + 9. * x) + 3. * y * (-4. + 9. * y) +
			3. * z * (-4. + 9. * z))) /
		(16. *
			exp((3. * (4. + 3. * x * (-4. + 9. * x) +
				3. * y * (-4. + 9. * y) + 3. * z * (-4. + 9. * z))) /
				4.0));
}

template< unsigned int Degree >
void ExecuteDirect
(
	const SimplexRefinableCellMesh< Dim , Degree > &simplexRefinableCellMesh ,
	const std::function< Point< double , Dim > ( typename SimplexMesh< Dim , Degree >::NodeMultiIndex ) > &NodePosition ,
	const std::function< bool ( typename SimplexMesh< Dim , Degree >::NodeMultiIndex ) > &IsBoundaryNode
)
{
	typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	Eigen::SparseMatrix< double > M , S;
	Timer timer;

	// Get the system matrices
	timer.reset();
	M = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().mass() * simplexRefinableCellMesh.P();

	S = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().stiffness() * simplexRefinableCellMesh.P();
	if( Verbose.set ) std::cout << "Got system matrices: " << timer.elapsed() << std::endl;

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

	Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > solver;

	// Compute the Cholesky factorization
	timer.reset();
	solver.compute( S );
	switch( solver.info() )
	{
		case Eigen::NumericalIssue: ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- numerical issue" );
		case Eigen::NoConvergence:  ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- no convergence" );
		case Eigen::InvalidInput:   ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- invalid input" );
		case Eigen::Success: ;
	}
	if( Verbose.set )
	{
		std::cout << "Performed Cholesky (LDLt) factorization: " << timer.elapsed() << std::endl;
		std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;
	}

	timer.reset();
	Eigen::VectorXd x = solver.solve( b );
	if( Verbose.set )
	{
		std::cout << "Solved direct system: " << timer.elapsed() << std::endl;
		double eIn = 0 , eOut = 0;
		Eigen::VectorXd r = b - S* x;
		for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
		std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;
	}

	double rms = 0;
	for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
	std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes() ) << std::endl;
}

template< unsigned int Degree , typename RelaxerType >
void ExecuteMG
(
	const HierarchicalSimplexRefinableCellMesh< Dim , Degree > &simplexRefinableCellMesh ,
	const std::function< Point< double , Dim > ( typename SimplexMesh< Dim , Degree >::NodeMultiIndex ) > &NodePosition ,
	const std::function< bool ( typename SimplexMesh< Dim , Degree >::NodeMultiIndex ) > &IsBoundaryNode ,
	unsigned int vCycles , unsigned int gsIters
)
{
	typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	Eigen::SparseMatrix< double > M , S;
	std::vector< Eigen::SparseMatrix< double > > P( CoarseNodeDimension.value );
	Timer timer;

	// Get the system matrices
	timer.reset();
	for( unsigned int d=0 ; d<(unsigned int)CoarseNodeDimension.value ; d++ ) P[d] = simplexRefinableCellMesh.P( d+1 , d );
	{
		Eigen::SparseMatrix< double > _P , _Pt;
		_P = simplexRefinableCellMesh.P( simplexRefinableCellMesh.maxLevel() , CoarseNodeDimension.value );
		_Pt = _P.transpose();
		M = _Pt * simplexRefinableCellMesh.simplexMesh().mass() * _P;
		S = _Pt * simplexRefinableCellMesh.simplexMesh().stiffness() * _P;
	}
	if( Verbose.set ) std::cout << "Got system matrices: " << timer.elapsed() << std::endl;


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
	timer.reset();
	MGSolver::Solver< RelaxerType > mgSolver( S , P , FullVerbose.set );
	if( Verbose.set ) std::cout << "Constructed multigrid system: " << timer.elapsed() << std::endl;

	if( Verbose.set && !FullVerbose.set ) std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;

	if( Debug.set )
	{
		for( unsigned int vCycles=1 ; vCycles<=32 ; vCycles*=2 )
		{
			std::cout << "V-Cycles[" << vCycles << "]" << std::endl;
			timer.reset();
			Eigen::VectorXd x = mgSolver.solve( b , vCycles , gsIters , gsIters , false );
			if( Verbose.set )
			{
				std::cout << "\tSolved multigrid system: " << timer.elapsed() << std::endl;
				double eIn = 0 , eOut = 0;
				Eigen::VectorXd r = b - S* x;
				for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
				std::cout << "\tError: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;
			}

			double rms = 0;
			for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
			std::cout << "\tRMS: " << sqrt( rms / simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) ) << std::endl;
		}
	}
	else
	{
		timer.reset();
		Eigen::VectorXd x = mgSolver.solve( b , vCycles , gsIters , gsIters , FullVerbose.set );
		if( Verbose.set )
		{
			std::cout << "Solved multigrid system: " << timer.elapsed() << std::endl;
			double eIn = 0 , eOut = 0;
			Eigen::VectorXd r = b - S* x;
			for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
			std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;
		}

		double rms = 0;
		for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
		std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes( CoarseNodeDimension.value ) ) << std::endl;
	}
}

template< unsigned int Degree >
void Execute( const Meshes::PolyhedronMesh< unsigned int > &polyMesh , const std::vector< Point< double , Dim > > &vertices )
{
	typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	Timer timer;

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
	eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

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
		switch( MGRelaxerType.value )
		{
			case MGSolver::RELAXER_JACOBI:
				return ExecuteMG< Degree , MGSolver::JacobiRelaxer >( simplexRefinableCellMesh , NodePosition , IsBoundaryNode , VCycles.value , Iterations.value );
			case MGSolver::RELAXER_GAUSS_SEIDEL:
				return ExecuteMG< Degree , MGSolver::GaussSeidelRelaxer >( simplexRefinableCellMesh , NodePosition , IsBoundaryNode , VCycles.value , Iterations.value );
			case MGSolver::RELAXER_PARALLEL_GAUSS_SEIDEL:
				return ExecuteMG< Degree , MGSolver::ParallelGaussSeidelRelaxer< 20 > >( simplexRefinableCellMesh , NodePosition , IsBoundaryNode , VCycles.value , Iterations.value );
			default: ERROR_OUT( "Unrecognized relaxer type: " , MGRelaxerType.value );
		}
	}
	else
	{
		SimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh;
		simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , CoarseNodeDimension.value , false , LinearPrecision.set , PlanarityEpsilon.value , Verbose.set );
		ExecuteDirect< Degree >( simplexRefinableCellMesh , NodePosition , IsBoundaryNode );
	}
}

int main( int argc , char* argv[] )
{
	if( argc<1 ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
	Misha::CmdLineParse( argc-1 , argv+1 , params );

	if( !Input.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	Verbose.set |= FullVerbose.set;

	std::vector< std::vector< std::pair< unsigned int , bool > > > polyhedra;
	std::vector< Meshes::Polygon< unsigned int > > polygons;
	std::vector< Point< double , Dim > > vertices;

	Meshes::PolyhedronMesh< unsigned int >::ReadOVM( Input.value , polyhedra , polygons , vertices );
	Meshes::PolyhedronMesh< unsigned int > polyMesh = Meshes::PolyhedronMesh< unsigned int >( polyhedra , polygons );


	switch( Degree.value )
	{
		case 1: Execute< 1 >( polyMesh , vertices ) ; break;
		case 2: Execute< 2 >( polyMesh , vertices ) ; break;
		case 3: Execute< 3 >( polyMesh , vertices ) ; break;
		case 4: Execute< 4 >( polyMesh , vertices ) ; break;
		default: ERROR_OUT( "Only degrees 1, 2, 3, or 4 supported" );
	}
}