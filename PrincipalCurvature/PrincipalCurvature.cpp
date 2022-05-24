#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <set>
#include <omp.h>
#include "Misha/PreProcess.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Misha/Exceptions.h"
#include "Misha/Geometry.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Misha/Miscellany.h"

const unsigned int ManifoldDimension = 2;

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineReadable Verbose( "verbose" ) , MinCurvature( "kMin" );
Misha::CmdLineReadable* params[] = { &In , &Out , &MinCurvature , &Verbose , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name.c_str() );
	printf( "\t[--%s <output mesh>]\n" , Out.name.c_str() );
	printf( "\t[--%s]\n" , MinCurvature.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

typedef VertexFactory::Factory< double , VertexFactory::PositionFactory< double , 3 > , VertexFactory::NormalFactory< double , 3 > > Factory;
typedef typename Factory::VertexType Vertex;

void WriteVectorField( const std::vector< Point< double , 3 > > & vec , const char * fileName )
{
	FILE * file;
	file = fopen( fileName , "wb" );
	if( !file )
	{
		fprintf( stderr , "Unable to open file for writing: %s\n" , fileName );
		exit( 1 );
	}
	int vecSize = (int)vec.size();
	fwrite( &vecSize , sizeof(int) , 1 , file );
	fwrite( &vec[0] , sizeof( Point< double , 3 > ) , vecSize , file );
	fclose( file );
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
	sStream << "Running Time [Principal Curvature (v. " << VERSION << ")]";
	Miscellany::NestedTimer timer( sStream.str() , Verbose.set );

	std::vector< std::vector< unsigned int > > polygons;
	std::vector< Vertex > vertices;

	Factory factory;
	int file_type;

	bool *readFlags = new bool[ factory.plyReadNum() ];
	PLY::ReadPolygons< Factory >( In.value , factory , vertices , polygons , readFlags , file_type );
	if( !factory.plyValidReadProperties<1>( readFlags ) ) ERROR_OUT( "No normals provided" );

	if( Verbose.set ) std::cout << "Vertices / Polygons: " << vertices.size() << " / " << polygons.size() << std::endl;


	std::vector< Point< double , 3 > > vectorField( polygons.size() );
#pragma omp parallel for
	for( int i=0 ; i<(int)polygons.size() ; i++ )
	{
		const std::vector< unsigned int > &p = polygons[i];

		// Compute the frame of the polygon
		Point< double , 3 > f[3];

		for( unsigned int j=0 ; j<p.size() ; j++ ) f[0] += vertices[ p[j] ].get<1>();
		f[0] /= Point< double , 3 >::Length( f[0] );

		f[1] = Point< double , 3 >::CrossProduct( f[0] , Point< double , 3 >( 1 , 0 , 0 ) );
		if( Point< double , 3 >::SquareNorm( f[1] )<1e-8 ) f[1] = Point< double , 3 >::CrossProduct( f[0] , Point< double , 3 >( 0 , 1 , 0 ) );
		f[1] /= Point< double , 3 >::Length( f[1] );

		f[2] = Point< double , 3 >::CrossProduct( f[0] , f[1] );
		f[2] /= Point< double , 3 >::Length( f[2] );

		// Looking for a matrix D minimizing:
		//	E(D) = \sum_i || D * dv - dn ||^2
		//       = \sum_i Tr( D^t * D  * dv * dv^t ) - \sum_i Tr( D^t * dn * dv^t ) - \sum_i Tr( D * dv * dn^t ) + ...
		//       = \sum_i Tr( D^t * D  * dv * dv^t ) - \sum_i Tr( D^t * dn * dv^t ) - \sum_i Tr( dn * dv^t * D^t ) + ...
		//       = \sum_i Tr( D^t * D  * dv * dv^t ) - \sum_i Tr( D^t * dn * dv^t ) - \sum_i Tr( D^t * dn * dv^t ) + ...
		// <=> D * \sum_i dv * dv^t = \sum_i dn * dv^t
		SquareMatrix< double , 2 > A , b;

		for( unsigned int j=0 ; j<p.size() ; j++ )
		{
			unsigned int j0 = p[j] , j1 = p[ (j+1)%p.size() ];
			Vertex d = vertices[j1] - vertices[j0];
			Point< double , 2 > dv , dn;
			for( unsigned int k=0 ; k<2 ; k++ )
			{
				dv[k] = Point< double , 3 >::Dot( d.get<0>() , f[k+1] );
				dn[k] = Point< double , 3 >::Dot( d.get<1>() , f[k+1] );
			}
			for( unsigned int k=0 ; k<2 ; k++ ) for( unsigned int l=0 ; l<2 ; l++ ) A(k,l) += dv[k] * dv[l] , b(k,l) += dn[k] * dv[l];
		}
		A = A.inverse() * b;
		A = ( A + A.transpose() ) / 2.;

		Eigen::Matrix< double , 2 , 2 > D;
		for( unsigned int j=0 ; j<2 ; j++ ) for( unsigned int k=0 ; k<2 ; k++ ) D(j,k) = A(j,k);

		Eigen::SelfAdjointEigenSolver< Eigen::Matrix< double , 2 , 2 > > es( D );
		Eigen::Vector< double , 2 > eVector = es.eigenvectors().col( MinCurvature.set ? 0 : 1 );
		double eValue = es.eigenvalues()( MinCurvature.set ? 0 : 1 );
		Point< double , 2 > v( eVector[0] * eValue , eVector[1] * eValue );

		vectorField[i] = f[1] * v[0] + f[2] * v[1];
	}

	if( Out.set ) WriteVectorField( vectorField , Out.value.c_str() );

	return EXIT_SUCCESS;
}