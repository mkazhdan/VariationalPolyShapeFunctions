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
#include "Misha/CmdLineParser.h"
#include "Misha/PlyVertexData.h"
#include "Misha/Miscellany.h"
#include "GeodesicsInHeat.inl"

const unsigned int ManifoldDimension = 2;

Misha::CmdLineParameter< std::string > In( "in" );
Misha::CmdLineParameter< int > Degree( "degree" , 2 ) , CoarseNodeDimension( "coarseDim" , 1 ) , Width( "width" , 512 ) , Height( "height" , 512 ) , SubdivisionIterations( "sub" , 0 );
Misha::CmdLineParameter< double > TimeStep( "time" , 1.e-3 ) , StiffnessRegularizer( "sRegularizer" , 1e-8 );
Misha::CmdLineReadable Multigrid( "mg" ) , Verbose( "verbose" ) , NoHelp( "noHelp" );
Misha::CmdLineReadable* params[] = { &In , &Degree , &TimeStep , &Verbose , &CoarseNodeDimension , &Multigrid , &Width , &Height , &SubdivisionIterations , &StiffnessRegularizer , &NoHelp , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name.c_str() );
	printf( "\t[--%s <element degree>=%d]\n" , Degree.name.c_str() , Degree.value );
	printf( "\t[--%s <coarse node dimensions>=%d]\n" , CoarseNodeDimension.name.c_str() , CoarseNodeDimension.value );
	printf( "\t[--%s <diffusion time step>=%e]\n" , TimeStep.name.c_str() , TimeStep.value );
	printf( "\t[--%s <stiffness regularizer>=%e]\n" , StiffnessRegularizer.name.c_str() , StiffnessRegularizer.value );
	printf( "\t[--%s <subdivision iterations>=%d]\n" , SubdivisionIterations.name.c_str() , SubdivisionIterations.value );
	printf( "\t[--%s <window width>=%d]\n" , Width.name.c_str() , Width.value );
	printf( "\t[--%s <window height>=%d]\n" , Height.name.c_str() , Height.value );
	printf( "\t[--%s]\n" , Multigrid.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
	printf( "\t[--%s]\n" , NoHelp.name.c_str() );
}

template< unsigned int Degree >
void Execute( int argc , char *argv[] , std::vector< Point3D< double > > &vertices , const std::vector< std::vector< unsigned int > > &polygons )
{
	if( Multigrid.set )
	{
		GeodesicsInHeatVisualization< Degree , true > v( vertices , polygons , CoarseNodeDimension.value , TimeStep.value , StiffnessRegularizer.value , SubdivisionIterations.value , Width.value , Height.value );

		std::stringstream sStream;
		sStream << "Geodesics in Heat: Degree=" << Degree << ", multigrid (v. " << VERSION << ")]";
		Misha::Viewable< GeodesicsInHeatVisualization< Degree , true > >::Viewer::Run( &v , argc , argv , sStream.str().c_str() );
	}
	else
	{
		GeodesicsInHeatVisualization< Degree , false > v( vertices , polygons , CoarseNodeDimension.value , TimeStep.value , StiffnessRegularizer.value , SubdivisionIterations.value , Width.value , Height.value );

		std::stringstream sStream;
		sStream << "Geodesics in Heat: Degree=" << Degree << ", direct (v. " << VERSION << ")]";
		Misha::Viewable< GeodesicsInHeatVisualization< Degree , false > >::Viewer::Run( &v , argc , argv , sStream.str().c_str() );
	}
}

int main( int argc , char* argv[] )
{
	typedef VertexFactory::PositionFactory< double , 3 > Factory;
	typedef typename Factory::VertexType Vertex;

	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	if( !NoHelp.set )
	{
		printf( "+----------------------------------------+\n" );
		printf( "| Interface Controls:                    |\n" );
		printf( "|    [Left Mouse]:            rotate     |\n" );
		printf( "|    [Left Mouse] + [CTRL]:   pan        |\n" );
		printf( "|    [Left Mouse] + [SHIFT]:  set seed   |\n" );
		printf( "+----------------------------------------+\n" );
	}

	std::vector< std::vector< unsigned int > > polygons;
	std::vector< Vertex > vertices;

	Factory factory;
	int file_type;
	PLY::ReadPolygons< Factory >( In.value , factory , vertices , polygons , NULL , file_type );
	if( Verbose.set ) std::cout << "Vertices / Polygons: " << vertices.size() << " / " << polygons.size() << std::endl;


	switch( Degree.value )
	{
		case 1: Execute< 1 >( argc , argv , vertices , polygons ) ; break;
		case 2: Execute< 2 >( argc , argv , vertices , polygons ) ; break;
		case 3: Execute< 3 >( argc , argv , vertices , polygons ) ; break;
		default: ERROR_OUT( "Only degrees 1, 2, and 3 supported" );
	}

	return EXIT_SUCCESS;
}