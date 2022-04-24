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
#include <omp.h>
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include <GL/glew.h>
#include <GL/glut.h>
#include <Misha/cmdLineParser.h>
#include "DeformablePolygonMesh.inl"

static const unsigned int Dim = 2;

Misha::CmdLineParameter< std::string > Input( "in" );
Misha::CmdLineParameter< int > Degree( "degree" , 2 ) , Width( "width" , 512 ) , Height( "height" , 512 ) , CoarseNodeDimension( "coarseDim" , 1 ) , RefinementResolution( "refine" , 8 ) , VCycles( "vCycles" , 1 ) , GSIterations( "gsIters" , 5 );
Misha::CmdLineParameterArray< float , Dim * Dim > AffineTransform( "xForm" );
Misha::CmdLineParameter< float > Gravity( "gravity" , -5e8f ) , TimeStep( "timeStep" );
Misha::CmdLineReadable Multigrid( "mg" ) , Lock( "lock" ) , NoHelp( "noHelp" );
Misha::CmdLineReadable* params[] = { &Input , &Width , &Height , &Degree , &AffineTransform , &Lock , &Gravity , &TimeStep , &RefinementResolution , &CoarseNodeDimension , &VCycles , &GSIterations , &Multigrid , &NoHelp , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s: <resolution>\n" , ex );
	printf( "\t --%s <input mesh>\n" , Input.name.c_str() );
	printf( "\t[--%s <screen width>=%d]\n" , Width.name.c_str() , Width.value );
	printf( "\t[--%s <screen height>=%d]\n" , Height.name.c_str() , Height.value );
	printf( "\t[--%s <fem degree>=%d]\n" , Degree.name.c_str() , Degree.value );
	printf( "\t[--%s <refinement resolution>=%d]\n" , RefinementResolution.name.c_str() , RefinementResolution.value );
	printf( "\t[--%s <xForm(0,0) , xForm(1,0) , xForm(0,1) , xForm(1,1) >=%f %f %f %f]\n" , AffineTransform.name.c_str() , 1. , 0. , 0. , 1. );
	printf( "\t[--%s <acceleration=%f>]\n" , Gravity.name.c_str() , Gravity.value );
	printf( "\t[--%s <time-step>]\n" , TimeStep.name.c_str() );
	printf( "\t[--%s <coarse node dimensions>=%d]\n" , CoarseNodeDimension.name.c_str() , CoarseNodeDimension.value );
	printf( "\t[--%s <v-cycles>=%d]\n" , VCycles.name.c_str() , VCycles.value );
	printf( "\t[--%s <Gauss-Seidel iterations>=%d]\n" , GSIterations.name.c_str() , GSIterations.value );
	printf( "\t[--%s]\n" , Multigrid.name.c_str() );
	printf( "\t[--%s]\n" , Lock.name.c_str() );
	printf( "\t[--%s]\n" , NoHelp.name.c_str() );
}

template< unsigned int Degree , unsigned int RefinementLevels >
void Execute( int argc , char *argv[] , const std::vector< Point< double , Dim > > &vertices , const std::vector< std::vector< unsigned int > > &polygons , const std::vector< bool > &lockedVertices , SquareMatrix< double , Dim > xForm , unsigned int width , unsigned int height , unsigned int refinementResolution )
{
	char windowName[1024];

	if( Multigrid.set )
	{
		sprintf( windowName , "Deformable Polygon Mesh Viewer (Degree=%d, multigrid)" , Degree );
		DeformablePolygonMeshVisualization< Degree , RefinementLevels , true > v;
		v.vCycles = VCycles.value;
		v.gsIters = GSIterations.value;
		v.init( vertices , polygons , lockedVertices , xForm , width , height , Gravity.value , CoarseNodeDimension.value );
		if( TimeStep.set ) v.setTimeStep( TimeStep.value );

		Misha::Viewable< DeformablePolygonMeshVisualization< Degree , RefinementLevels , true > >::Viewer::Run( &v , argc , argv , windowName );
	}
	else
	{
		sprintf( windowName , "Deformable Polygon Mesh Viewer (Degree=%d, direct)" , Degree );
		DeformablePolygonMeshVisualization< Degree , RefinementLevels , false > v;
		v.init( vertices , polygons , lockedVertices , xForm , width , height , Gravity.value , CoarseNodeDimension.value );
		if( TimeStep.set ) v.setTimeStep( TimeStep.value );

		Misha::Viewable< DeformablePolygonMeshVisualization< Degree , RefinementLevels , false > >::Viewer::Run( &v , argc , argv , windowName );
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

	if( !NoHelp.set )
	{
		printf( "+------------------------------------+\n" );
		printf( "| Interface Controls:                |\n" );
		printf( "|    [Left Mouse]: select vertex     |\n" );
		printf( "+------------------------------------+\n" );
	}

	SquareMatrix< double , Dim > xForm = SquareMatrix< double , Dim >::Identity();
	if( AffineTransform.set ) for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ )  xForm(i,j) = AffineTransform.values[j*Dim+i];
	std::vector< Point< double , Dim > > vertices;
	std::vector< SimplexIndex< Dim , unsigned int  > > triangles;
	std::vector< std::vector< unsigned int > > polygons;
	std::vector< bool > lockedVertices;

	typedef VertexFactory::PositionFactory< double , Dim > Factory;
	typedef typename Factory::VertexType Vertex;

	Factory factory;
	int file_type;
	PLY::ReadPolygons< Factory >( Input.value , factory , vertices , polygons , NULL , file_type , NULL );
	double xMin = std::numeric_limits< double >::infinity();
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) xMin = std::min< double >( xMin , vertices[i][0] );
	lockedVertices.resize( vertices.size() , false );
	if( Lock.set ) for( unsigned int i=0 ; i<vertices.size() ; i++ ) if( fabs( vertices[i][0]-xMin )<1e-5 ) lockedVertices[i] = true;

	static const unsigned int RefinementLevels = 2;

	switch( Degree.value )
	{
		case 1: Execute< 1 , RefinementLevels >( argc , argv , vertices , polygons , lockedVertices , xForm , Width.value , Height.value , RefinementResolution.value ) ; break;
		case 2: Execute< 2 , RefinementLevels >( argc , argv , vertices , polygons , lockedVertices , xForm , Width.value , Height.value , RefinementResolution.value ) ; break;
		case 3: Execute< 3 , RefinementLevels >( argc , argv , vertices , polygons , lockedVertices , xForm , Width.value , Height.value , RefinementResolution.value ) ; break;
		case 4: Execute< 4 , RefinementLevels >( argc , argv , vertices , polygons , lockedVertices , xForm , Width.value , Height.value , RefinementResolution.value ) ; break;
		default: ERROR_OUT( "Only degrees 1, 2, 3, or 4 supported" );
	}

	return EXIT_SUCCESS;
}