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
#include "Misha/PreProcess.h"
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

enum
{
	VALUE_POSITION=0 ,
	VALUE_NORMAL ,
	VALUE_COLOR ,
	VALUE_COUNT
};
const std::string ValueNames[] = { "position" , "normal" , "color" };

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< int > Degree( "degree" , 2 ) , CoarseNodeDimension( "coarseDim" , 1 ) , VCycles( "vCycles" , 3 ) , GSIters( "gsIters" , 5 ) , ValueType( "value" , VALUE_POSITION );
Misha::CmdLineParameter< double > GradientWeight( "gWeight" , 1.e-5 ) , GradientScale( "gScale" , 1. );
Misha::CmdLineReadable Multigrid( "mg" ) , Verbose( "verbose" );
Misha::CmdLineReadable* params[] = { &In , &Out , &GradientScale , &Degree , &GradientWeight , &Verbose , &CoarseNodeDimension , &Multigrid , &VCycles , &GSIters , &ValueType , NULL };

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
	printf( "\t[--%s <value type>=%d]\n" , ValueType.name.c_str() , ValueType.value );
	for( unsigned int i=0 ; i<VALUE_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i , ValueNames[i].c_str() );
	printf( "\t[--%s]\n" , Multigrid.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

typedef VertexFactory::Factory< double , VertexFactory::PositionFactory< double , 3 > , VertexFactory::NormalFactory< double , 3 > , VertexFactory::RGBColorFactory< double > > OrientedRGBFactory;
typedef typename OrientedRGBFactory::VertexType OrientedRGBVertex;

template< unsigned int Degree , bool Hierarchical >
void Execute( const std::vector< OrientedRGBVertex > &vertices , const std::vector< std::vector< unsigned int > > &polygons , std::vector< Point< double , 3 > > &values , bool normalize )
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

	Eigen::SparseMatrix< double > S , Pt , P;
	std::vector< Eigen::SparseMatrix< double > > Ps;
	PointVector< 3 > nodeValues;

	Eigen::SparseMatrix< double > M , L;
	{
		Miscellany::NestedTimer timer( "Got system" , Verbose.set );
		if constexpr( Hierarchical )
		{
			P = pMesh.P( pMesh.maxLevel() , CoarseNodeDimension.value );
			Ps.resize( CoarseNodeDimension.value );
			for( unsigned int d=0 ; d<(unsigned int)CoarseNodeDimension.value ; d++ ) Ps[d] = pMesh.P( d+1 , d );
		}
		else P = pMesh.P();
		Pt = P.transpose();
		S = Pt * pMesh.simplexMesh().stiffness() * P;

		M = Pt * pMesh.simplexMesh().mass() * P;
		L = M + S * GradientWeight.value;
		nodeValues.resize( pMesh.nodes() );

		// Use the vertex normals to set the node coefficients
		for( auto iter=pMesh.nodeMap().begin() ; iter!=pMesh.nodeMap().end() ; iter++ )
		{
			NodeMultiIndex nmi = iter->first;
			Point3D< double > v;
			for( unsigned int i=0 ; i<Degree ; i++ ) v += values[ nmi[i] ];
			if( normalize ) nodeValues[ iter->second ] = v / Point< double , 3 >::Length( v );
			else            nodeValues[ iter->second ] = v / Degree;
		}
	}
	{
		Miscellany::NestedTimer timer( "Solved system" , Verbose.set );
		if( Ps.size() )
		{
			MGSolver::Solver< MGSolver::ParallelGaussSeidelRelaxer< 20u > > mgSolver( L , Ps , false );
			PointVector< 3 > b;
			b = M * nodeValues + S * nodeValues * GradientWeight.value * GradientScale.value;
			nodeValues = mgSolver.solve( b , nodeValues , VCycles.value , GSIters.value , GSIters.value , false );
		}
		else
		{
			SparseSolver::LLT solver( L );
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
	}
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	std::stringstream sStream;
	sStream << "Running Time [Gradient Domain Process (V" << VERSION << ")]";
	Miscellany::NestedTimer timer( sStream.str() , Verbose.set );

	std::vector< std::vector< unsigned int > > polygons;
	std::vector< OrientedRGBVertex > vertices;

	OrientedRGBFactory factory;
	int file_type;

	bool *readFlags = new bool[ factory.plyReadNum() ];
	PLY::ReadPolygons< OrientedRGBFactory >( In.value , factory , vertices , polygons , readFlags , file_type );
	bool hasNormal = factory.plyValidReadProperties<1>( readFlags );
	bool hasColor = factory.plyValidReadProperties<2>( readFlags );
	delete[] readFlags;
	if( ValueType.value==VALUE_COLOR && !hasColor ) ERROR_OUT( "Could not read in color" );
	else if( ValueType.value==VALUE_NORMAL && !hasNormal )
	{
		WARN( "No normals found. Setting from face normals" );
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].get<1>() = Point< double , 3 >();
		for( unsigned int i=0 ; i<polygons.size() ; i++ )
		{
			Point< double , 3 > n;
			const std::vector< unsigned int > &p = polygons[i];
			Point< double , 3 > center;
			for( unsigned int j=0 ; j<p.size() ; j++ ) center += vertices[ p[j] ].get<0>();
			center /= (double)p.size();
			for( unsigned int j=0 ; j<p.size() ; j++ )
			{
				Point< double , 3 > v1 = vertices[ p[j] ].get<0>() , v2 = vertices[ p[(j+1)%p.size()] ].get<0>();
				n += Point< double , 3 >::CrossProduct( v1-center , v2-center );
			}
			for( unsigned int j=0 ; j<p.size() ; j++ ) vertices[ p[j] ].get<1>() += n;
		}
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].get<1>() = vertices[i].get<1>() / Point< double , 3 >::Length( vertices[i].get<1>() );
		hasNormal = true;
	}
	if( Verbose.set ) std::cout << "Vertices / Polygons: " << vertices.size() << " / " << polygons.size() << std::endl;

	std::vector< Point< double , 3 > > values( vertices.size() );

	switch( ValueType.value )
	{
	case VALUE_POSITION:
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) values[i] = vertices[i].get<0>();
		break;
	case VALUE_NORMAL:
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) values[i] = vertices[i].get<1>();
		break;
	case VALUE_COLOR:
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) values[i] = vertices[i].get<2>();
		break;
	default:
		ERROR_OUT( "Unrecognized value type: " , ValueType.value );
	}

	if( Multigrid.set )
	{
		switch( Degree.value )
		{
			case 1: Execute< 1 , true >( vertices , polygons , values , ValueType.value==VALUE_NORMAL ) ; break;
			case 2: Execute< 2 , true >( vertices , polygons , values , ValueType.value==VALUE_NORMAL ) ; break;
			case 3: Execute< 3 , true >( vertices , polygons , values , ValueType.value==VALUE_NORMAL ) ; break;
			default: ERROR_OUT( "Only degrees 1, 2, and 3 supported" );
		}
	}
	else
	{
		switch( Degree.value )
		{
			case 1: Execute< 1 , false >( vertices , polygons , values , ValueType.value==VALUE_NORMAL ) ; break;
			case 2: Execute< 2 , false >( vertices , polygons , values , ValueType.value==VALUE_NORMAL ) ; break;
			case 3: Execute< 3 , false >( vertices , polygons , values , ValueType.value==VALUE_NORMAL ) ; break;
			default: ERROR_OUT( "Only degrees 1, 2, and 3 supported" );
		}
	}

	auto ClampColor = []( Point< double , 3 > c )
	{
		for( unsigned int i=0 ; i<3 ; i++ ) c[i] = std::max< double >( 0. , std::min< double >( 255. , c[i] ) );
		return c;
	};


	switch( ValueType.value )
	{
	case VALUE_POSITION:
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].get<0>() = values[i];
		break;
	case VALUE_NORMAL:
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].get<1>() = values[i] = values[i] / Point< double , 3 >::Length( values[i] );
		break;
	case VALUE_COLOR:
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].get<2>() = ClampColor( values[i] );
		break;
	}

	if( Out.set )
	{
		if( hasColor && hasNormal ) PLY::WritePolygons( Out.value , factory , vertices , polygons , file_type );
		else if( hasColor )
		{
			typedef VertexFactory::Factory< double , VertexFactory::PositionFactory< double , 3 > , VertexFactory::RGBColorFactory< double > > _Factory;
			typedef typename _Factory::VertexType _Vertex;

			_Factory _factory;
			std::vector< _Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].get<0>() = vertices[i].get<0>() , _vertices[i].get<1>() = vertices[i].get<2>();
			PLY::WritePolygons( Out.value , _factory , _vertices , polygons , file_type );
		}
		else if( hasNormal )
		{
			typedef VertexFactory::Factory< double , VertexFactory::PositionFactory< double , 3 > , VertexFactory::NormalFactory< double , 3 > > _Factory;
			typedef typename _Factory::VertexType _Vertex;

			_Factory _factory;
			std::vector< _Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].get<0>() = vertices[i].get<0>() , _vertices[i].get<1>() = vertices[i].get<1>();
			PLY::WritePolygons( Out.value , _factory , _vertices , polygons , file_type );
		}
		else
		{
			typedef VertexFactory::PositionFactory< double , 3 > _Factory;
			typedef typename _Factory::VertexType _Vertex;

			_Factory _factory;
			std::vector< _Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i] = vertices[i].get<0>();
			PLY::WritePolygons( Out.value , _factory , _vertices , polygons , file_type );
		}
	}
	return EXIT_SUCCESS;
}