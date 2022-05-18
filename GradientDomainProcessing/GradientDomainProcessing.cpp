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
#include <omp.h>
#include "Eigen/Sparse"
#include "Misha/Exceptions.h"
#include "Misha/Geometry.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Misha/Miscellany.h"
#include "Misha/MGSolver.h"
#include "Misha/SimplexRefinableMesh.h"
#include "Misha/Meshes.h"

const unsigned int ManifoldDimension = 2;

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< int > Degree( "degree" , 2 ) , CoarseNodeDimension( "coarseDim" , 1 ) , VCycles( "vCycles" , 3 ) , GSIters( "gsIters" , 5 );
Misha::CmdLineParameter< double > GradientWeight( "gWeight" , 1.e-5 ) , GradientScale( "gScale" , 1. );
Misha::CmdLineReadable Multigrid( "mg" ) , Verbose( "verbose" ) , Color( "color" );
Misha::CmdLineReadable* params[] = { &In , &Out , &GradientScale , &Degree , &GradientWeight , &Verbose , &CoarseNodeDimension , &Multigrid , &VCycles , &GSIters , &Color , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name.c_str() );
	printf( "\t[--%s <output mesh>]\n" , Out.name.c_str() );
	printf( "\t[--%s <element degree>=%d]\n" , Degree.name.c_str() , Degree.value );
	printf( "\t[--%s <coarse node dimensions>=%d]\n" , CoarseNodeDimension.name.c_str() , CoarseNodeDimension.value );
	printf( "\t[--%s <gradient weight>=%e]\n" , GradientWeight.name.c_str() , GradientWeight.value );
	printf( "\t[--%s <gradient scale>=%e]\n" , GradientScale.name.c_str() , GradientScale.value );
	printf( "\t[--%s <v-cycles>=%d]\n" , VCycles.name.c_str() , VCycles.value );
	printf( "\t[--%s <Gauss-Seidel iterations>=%d]\n" , GSIters.name.c_str() , GSIters.value );
	printf( "\t[--%s]\n" , Multigrid.name.c_str() );
	printf( "\t[--%s]\n" , Color.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

typedef VertexFactory::PositionFactory< double , 3 > Factory;
typedef typename Factory::VertexType Vertex;

typedef VertexFactory::Factory< double , VertexFactory::PositionFactory< double , 3 > , VertexFactory::RGBColorFactory< double > > RGBFactory;
typedef typename RGBFactory::VertexType RGBVertex;

template< unsigned int Degree , bool Hierarchical >
void Execute( const std::vector< RGBVertex > &vertices , const std::vector< std::vector< unsigned int > > &polygons , std::vector< Point< double , 3 > > &values )
{
	using MGSolver::PointVector;
	using MGSolver::operator +;
	using MGSolver::operator -;
	using MGSolver::operator *;
	using MGSolver::operator /;

	typedef typename std::conditional< Hierarchical , HierarchicalSimplexRefinableCellMesh< ManifoldDimension , Degree > , SimplexRefinableCellMesh< ManifoldDimension , Degree > >::type SimplexRefinableMesh;
	typedef typename SimplexMesh< ManifoldDimension , Degree >::NodeMultiIndex NodeMultiIndex;

	Meshes::VertexPositionFunction< 3 , unsigned int > vertexPositionFunction = [&]( unsigned int v ){ return Point3D< double >( vertices[v].get<0>() ); };

	SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
	eWeights.kWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1;

	Meshes::PolygonMesh< unsigned int > polygonMesh = Meshes::PolygonMesh< unsigned int >( polygons );
	Meshes::FullVertexPositionFunction< 3 > fullVertexPositionFunction = polygonMesh.fullVertexPositionFunction( vertexPositionFunction , true );

	SimplexRefinableMesh pMesh;
	if constexpr( Hierarchical )
	{
		pMesh = polygonMesh.template hierarchicalSimplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , CoarseNodeDimension.value , false , false , 0 , Verbose.set );
	}
	else
	{
		pMesh = polygonMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , CoarseNodeDimension.value , false , false , 0 , Verbose.set );
	}
	pMesh.makeUnitVolume();

	const SimplexMesh< ManifoldDimension , Degree > &sMesh = pMesh.simplexMesh();

	Miscellany::Timer timer;
	Eigen::SparseMatrix< double > S , Pt , P;
	std::vector< Eigen::SparseMatrix< double > > Ps;

	if constexpr( Hierarchical )
	{
		P = pMesh.P( pMesh.maxLevel() , CoarseNodeDimension.value );
		Ps.resize( CoarseNodeDimension.value );
		for( unsigned int d=0 ; d<(unsigned int)CoarseNodeDimension.value ; d++ ) Ps[d] = pMesh.P( d+1 , d );
	}
	else P = pMesh.P();
	Pt = P.transpose();
	S = Pt * pMesh.simplexMesh().stiffness() * P;

	Eigen::SparseMatrix< double > M = Pt * pMesh.simplexMesh().mass() * P;
	Eigen::SparseMatrix< double > L = M + S * GradientWeight.value;
	PointVector< 3 > nodeValues( pMesh.nodes() );

	// Use the vertex normals to set the node coefficients
	for( auto iter=pMesh.nodeMap().begin() ; iter!=pMesh.nodeMap().end() ; iter++ )
	{
		NodeMultiIndex nmi = iter->first;
		Point3D< double > v;
		for( unsigned int i=0 ; i<Degree ; i++ ) v += values[ nmi[i] ];
		nodeValues[iter->second] = v / Degree;
	}
	if( Verbose.set ) std::cout << "Got system: " << timer.elapsed() << " (s)" << std::endl;

	timer.reset();
	if( Ps.size() )
	{
		MGSolver::Solver< MGSolver::ParallelGaussSeidelRelaxer< 20u > > mgSolver( L , Ps , false );
		PointVector< 3 > b;
		b = M * nodeValues + S * nodeValues * GradientWeight.value * GradientScale.value;
		nodeValues = mgSolver.solve( b , nodeValues , VCycles.value , GSIters.value , GSIters.value , false );
	}
	else
	{
		Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > solver( L );
		Eigen::VectorXd x( pMesh.nodes() ) , b;
		for( unsigned int d=0 ; d<3 ; d++ )
		{
			for( unsigned int i=0 ; i<pMesh.nodes() ; i++ ) x[i] = nodeValues[i][d];
			b = M * x + S * x * GradientWeight.value * GradientScale.value;
			x = solver.solve( b );
			for( unsigned int i=0 ; i<pMesh.nodes() ; i++ ) nodeValues[i][d] = x[i];
		}
	}
	for( unsigned int i=0 ; i<values.size() ; i++ )
	{
		unsigned int idx[Degree];
		for( unsigned int d=0 ; d<Degree ; d++ ) idx[d] = i;
		values[i] = nodeValues[ pMesh.nodeIndex( NodeMultiIndex(idx) ) ];
	}

	if( Verbose.set ) std::cout << "Solved system: " << timer.elapsed() << " (s)" << std::endl;
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	auto Random = []( void ){ return 1.-2.*(double)rand()/RAND_MAX; };
	auto RandomBallPoint = [&]( void )
	{
		Point3D< double > r;
		while( true )
		{
			r = Point3D< double >( Random() , Random() , Random() );
			if( r.squareNorm()<1 ) return r;
		}
	};


	std::vector< std::vector< unsigned int > > polygons;
	std::vector< RGBVertex > vertices;
	std::vector< Point< double , 3 > > vectorField;

	RGBFactory factory;
	int file_type;

	bool *readFlags = new bool[ factory.plyReadNum() ];
	PLY::ReadPolygons< RGBFactory >( In.value , factory , vertices , polygons , readFlags , file_type );
	bool hasColor = factory.plyValidReadProperties<1>( readFlags );
	delete[] readFlags;
	if( !hasColor && Color.set ) ERROR_OUT( "Could not read in color" );

	if( Verbose.set ) std::cout << "Vertices / Polygons: " << vertices.size() << " / " << polygons.size() << std::endl;

	std::vector< Point< double , 3 > > values( vertices.size() );
	if( Color.set ) for( unsigned int i=0 ; i<vertices.size() ; i++ ) values[i] = vertices[i].get<1>();
	else            for( unsigned int i=0 ; i<vertices.size() ; i++ ) values[i] = vertices[i].get<0>();

	if( Multigrid.set )
	{
		switch( Degree.value )
		{
		case 1: Execute< 1 , true >( vertices , polygons , values ) ; break;
		case 2: Execute< 2 , true >( vertices , polygons , values ) ; break;
		case 3: Execute< 3 , true >( vertices , polygons , values ) ; break;
		default: ERROR_OUT( "Only degrees 1, 2, and 3 supported" );
		}
	}
	else
	{
		switch( Degree.value )
		{
		case 1: Execute< 1 , false >( vertices , polygons , values ) ; break;
		case 2: Execute< 2 , false >( vertices , polygons , values ) ; break;
		case 3: Execute< 3 , false >( vertices , polygons , values ) ; break;
		default: ERROR_OUT( "Only degrees 1, 2, and 3 supported" );
		}
	}

	auto ClampColor = []( Point< double , 3 > c )
	{
		for( unsigned int i=0 ; i<3 ; i++ ) c[i] = std::max< double >( 0. , std::min< double >( 255. , c[i] ) );
		return c;
	};


	if( Color.set ) for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].get<1>() = ClampColor( values[i] );
	else            for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].get<0>() = values[i];

	if( Out.set )
	{
		if( hasColor ) PLY::WritePolygons( Out.value , factory , vertices , polygons , file_type );
		else
		{
			Factory _factory;
			std::vector< Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i] = vertices[i].get<0>();
			PLY::WritePolygons( Out.value , _factory , _vertices , polygons , file_type );
		}
	}
	return EXIT_SUCCESS;
}