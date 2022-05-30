#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <set>
#include <sstream>
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

Misha::CmdLineParameter< std::string > In( "in" ) , InVectorField( "inVF" ) , Out( "out" );
Misha::CmdLineParameter< int > Degree( "degree" , 2 ) , CoarseNodeDimension( "coarseDim" , 1 ) , VCycles( "vCycles" , 20 ) , GSIters( "gsIters" , 5 );
Misha::CmdLineParameter< double > StepSize( "stepSize" , 1.e-4 ) , SharpeningStepSize( "sStepSize" , 1e-5 ) , MetricScale( "mScale" , 1e4 ) , PlanarityEpsilon( "pEps" , 0 ) , Sharpen( "sharpen" , 16. );
Misha::CmdLineReadable Multigrid( "mg" ) , Verbose( "verbose" ) , PartitionOfUnity( "pou" ) , LinearPrecision( "linear" );
Misha::CmdLineReadable* params[] = { &In , &InVectorField , &Out , &MetricScale , &Degree , &StepSize , &Verbose , &CoarseNodeDimension , &Multigrid , &VCycles , &GSIters , &PartitionOfUnity , &LinearPrecision , &PlanarityEpsilon , &Sharpen , &SharpeningStepSize , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name.c_str() );
	printf( "\t --%s <input vector field>\n" , InVectorField.name.c_str() );
	printf( "\t[--%s <output mesh>]\n" , Out.name.c_str() );
	printf( "\t[--%s <element degree>=%d]\n" , Degree.name.c_str() , Degree.value );
	printf( "\t[--%s <coarse node dimensions>=%d]\n" , CoarseNodeDimension.name.c_str() , CoarseNodeDimension.value );
	printf( "\t[--%s <step size>=%e]\n" , StepSize.name.c_str() , StepSize.value );
	printf( "\t[--%s <sharpening step size>=%e]\n" , SharpeningStepSize.name.c_str() , SharpeningStepSize.value );
	printf( "\t[--%s <metric scale>=%e]\n" , MetricScale.name.c_str() , MetricScale.value );
	printf( "\t[--%s <sharpening gradient scale>=%e]\n" , Sharpen.name.c_str() , Sharpen.value );
	printf( "\t[--%s <v-cycles>=%d]\n" , VCycles.name.c_str() , VCycles.value );
	printf( "\t[--%s <Gauss-Seidel iterations>=%d]\n" , GSIters.name.c_str() , GSIters.value );
	printf( "\t[--%s <planarity epsilon>=%g]\n" , PlanarityEpsilon.name.c_str() , PlanarityEpsilon.value );
	printf( "\t[--%s]\n" , Multigrid.name.c_str() );
	printf( "\t[--%s]\n" , PartitionOfUnity.name.c_str() );
	printf( "\t[--%s]\n" , LinearPrecision.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

typedef VertexFactory::PositionFactory< double , 3 > Factory;
typedef typename Factory::VertexType Vertex;

typedef VertexFactory::Factory< double , VertexFactory::PositionFactory< double , 3 > , VertexFactory::RGBColorFactory< double > > ColorFactory;
typedef typename ColorFactory::VertexType ColorVertex;

void ReadVectorField( std::vector< Point< double , 3 > > & vec , const char * fileName )
{
	FILE * file;
	file = fopen( fileName , "rb" );
	if( !file )
	{
		fprintf( stderr , "Unable to open file for reading: %s\n" , fileName );
		exit( 1 );
	}
	int vecSize;
	fread( &vecSize , sizeof(int) , 1 , file );
	vec.resize( vecSize );
	fread( &vec[0] , sizeof( Point< double , 3 > ) , vecSize , file );
	fclose( file );
}

template< unsigned int Degree , bool Hierarchical >
using SimplexRefinableMesh = typename std::conditional< Hierarchical , HierarchicalSimplexRefinableCellMesh< ManifoldDimension , Degree > , SimplexRefinableCellMesh< ManifoldDimension , Degree > >::type;

template< unsigned int Degree >
std::vector< SquareMatrix< double , ManifoldDimension > > LICMetric( const SimplexMesh< ManifoldDimension , Degree > &sMesh , const Meshes::FullVertexPositionFunction< 3 > &fullVertexPositionFunction , const std::vector< Point< double , 3 > > &vectorField )
{
	std::vector< SquareMatrix< double , ManifoldDimension > > metric( sMesh.simplices() );
	for( unsigned int i=0 ; i<sMesh.simplices() ; i++ )
	{
		SimplexIndex< ManifoldDimension , unsigned int > s = sMesh.simplex( i );
		Point< double , 3 > frame[2];
		Point< double , 3 > v[] = { fullVertexPositionFunction( s[0] ) , fullVertexPositionFunction( s[1] ) , fullVertexPositionFunction( s[2] ) };
		Point< double , 3 > d[] = { v[1]-v[0] , v[2]-v[0] };
		Point< double , 3 > n = Point< double , 3 >::CrossProduct( d[0] , d[1] );
		double scale = Point< double , 3 >::Length( vectorField[i] );
		n /= Point< double , 3 >::Length( n );

		if( scale ) frame[0] = vectorField[i] / scale;
		else
		{
			frame[0] = Point< double , 3 >::CrossProduct( n , Point< double , 3 >(1,0,0) );
			if( Point< double , 3 >::SquareNorm( frame[0] ) < 1e-4 ) Point< double , 3 >::CrossProduct( n , Point< double , 3 >(0,1,0) );
			frame[0] /= Point< double , 3 >::Length( frame[0] );
		}
		frame[1] = Point< double , 3 >::CrossProduct( n , frame[0] );
		frame[1] /= Point< double , 3 >::Length( frame[1] );

		Point< double , 2 > _d[2];
		for( unsigned int j=0 ; j<2 ; j++ ) for( unsigned int k=0 ; k<2 ; k++ ) _d[j][k] = Point< double , 3 >::Dot( d[j] , frame[k] );
		SquareMatrix< double , 2 > A;
		A(0,0) = 1./( 1. + scale*MetricScale.value );
		A(1,1) = 1.;
		for( unsigned int j=0 ; j<2 ; j++ ) for( unsigned int k=0 ; k<2 ; k++ ) metric[i](j,k) = Point< double , 2 >::Dot( _d[j] , A*_d[k] );
	}
	return metric;
}

template< unsigned int Degree >
std::vector< SquareMatrix< double , ManifoldDimension > > EmbeddingMetric( const SimplexMesh< ManifoldDimension , Degree > &sMesh , const Meshes::FullVertexPositionFunction< 3 > &fullVertexPositionFunction )
{
	std::vector< SquareMatrix< double , ManifoldDimension > > metric( sMesh.simplices() );
	for( unsigned int i=0 ; i<sMesh.simplices() ; i++ )
	{
		SimplexIndex< ManifoldDimension , unsigned int > s = sMesh.simplex( i );
		Point< double , 3 > v[] = { fullVertexPositionFunction( s[0] ) , fullVertexPositionFunction( s[1] ) , fullVertexPositionFunction( s[2] ) };
		Point< double , 3 > d[] = { v[1]-v[0] , v[2]-v[0] };
		for( unsigned int j=0 ; j<2 ; j++ ) for( unsigned int k=0 ; k<2 ; k++ ) metric[i](j,k) = Point< double , 3 >::Dot( d[j] , d[k] );
	}
	return metric;
}

#if 0
template< unsigned int Dim > using PointVector = Eigen::Matrix< Point< double , Dim > , Eigen::Dynamic , 1 >;
template< unsigned int Dim >
PointVector< Dim > operator * ( const PointVector< Dim > &v , double s )
{
	PointVector< Dim > w( v.rows() );
	for( unsigned int i=0 ; i<v.rows() ; i++ ) w[i] = v[i] * s;
	return w;
}

template< unsigned int Dim > using PointVector = Eigen::Matrix< Point< double , Dim > , Eigen::Dynamic , 1 >;
template< unsigned int Dim >
PointVector< Dim > operator + ( const PointVector< Dim > &v1 , const PointVector< Dim > &v2 )
{
	assert( v1.rows()==v2.rows() );
	PointVector< Dim > w( v1.rows() );
	for( unsigned int i=0 ; i<v1.rows() ; i++ ) w[i] = v1[i] + v2[i];
	return w;
}

template< unsigned int Dim >
PointVector< Dim > operator * ( const Eigen::SparseMatrix< double > &M , const PointVector< Dim > &v )
{
	PointVector< Dim > w( M.rows() );
	for( unsigned int i=0 ; i<M.rows() ; i++ ) w[i] = Point< double , Dim >();
	for( unsigned int i=0 ; i<M.outerSize() ; i++ )
		for( Eigen::SparseMatrix< double >::InnerIterator iter(M,i) ; iter ; ++iter )
			w[i] += v[ iter.col() ] * iter.value();
	return w;
}
#endif

template< unsigned int Degree , bool Hierarchical >
std::vector< Point< double , 3 > > Execute( const std::vector< Vertex > &vertices , const std::vector< std::vector< unsigned int > > &polygons , const std::vector< Point< double , 3 > > &vectorField )
{
	using MGSolver::PointVector;
	using MGSolver::operator +;
	using MGSolver::operator -;
	using MGSolver::operator *;
	using MGSolver::operator /;

	typedef typename std::conditional< Hierarchical , HierarchicalSimplexRefinableCellMesh< ManifoldDimension , Degree > , SimplexRefinableCellMesh< ManifoldDimension , Degree > >::type SimplexRefinableMesh;

	std::vector< Point< double , 3 > > colors( vertices.size() );
	std::vector< Point< double , 3 > > simplexVectorField;

	typedef typename SimplexMesh< ManifoldDimension , Degree >::NodeMultiIndex NodeMultiIndex;

	Meshes::VertexPositionFunction< 3 , unsigned int > vertexPositionFunction = [&]( unsigned int v ){ return Point3D< double >( vertices[v] ); };

	SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
	eWeights.kWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1;

	Meshes::PolygonMesh< unsigned int > polygonMesh = Meshes::PolygonMesh< unsigned int >( polygons );
	Meshes::FullVertexPositionFunction< 3 > fullVertexPositionFunction = polygonMesh.fullVertexPositionFunction( vertexPositionFunction , true );

	for( unsigned int i=0 ; i<polygonMesh.polygons() ; i++ ) for( unsigned int j=0 ; j<polygonMesh.polygon(i).size() ; j++ )
		simplexVectorField.push_back( vectorField[i] );

	SimplexRefinableMesh pMesh;
	if constexpr( Hierarchical )
	{
		pMesh = polygonMesh.template hierarchicalSimplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , CoarseNodeDimension.value , PartitionOfUnity.set , LinearPrecision.set , PlanarityEpsilon.value , Verbose.set );
	}
	else
	{
		pMesh = polygonMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , CoarseNodeDimension.value , PartitionOfUnity.set , LinearPrecision.set , PlanarityEpsilon.value , Verbose.set );
	}

	Eigen::SparseMatrix< double > Pt , P;
	std::vector< Eigen::SparseMatrix< double > > Ps;

	{
		Miscellany::NestedTimer timer( "Gathered prolongations" , Verbose.set );
		if constexpr( Hierarchical )
		{
			P = pMesh.P( pMesh.maxLevel() , CoarseNodeDimension.value );
			Ps.resize( CoarseNodeDimension.value );
			for( unsigned int d=0 ; d<(unsigned int)CoarseNodeDimension.value ; d++ ) Ps[d] = pMesh.P( d+1 , d );
		}
		else P = pMesh.P();
		Pt = P.transpose();
	}

	const SimplexMesh< ManifoldDimension , Degree > &sMesh = pMesh.simplexMesh();

	Eigen::SparseMatrix< double > licS , licM , licL , S , M , L;

	PointVector< 3 > nodeColors;
	{
		Miscellany::NestedTimer timer( "Got system" , Verbose.set );
		std::vector< SquareMatrix< double , ManifoldDimension > > licMetric = LICMetric( sMesh , fullVertexPositionFunction , simplexVectorField );
		pMesh.setMetric( [&]( unsigned int idx ){ return licMetric[idx]; } );
		pMesh.makeUnitVolume();
		licS = Pt * pMesh.simplexMesh().stiffness() * P;
		licM = Pt * pMesh.simplexMesh().mass() * P;
		licL = licM + licS * StepSize.value;

		if( Sharpen.value>0 )
		{
			std::vector< SquareMatrix< double , ManifoldDimension > > metric = LICMetric( sMesh , fullVertexPositionFunction , simplexVectorField );
			pMesh.setMetric( [&]( unsigned int idx ){ return metric[idx]; } );
			pMesh.makeUnitVolume();
			S = Pt * pMesh.simplexMesh().stiffness() * P;
			M = Pt * pMesh.simplexMesh().mass() * P;
			L = M + S * SharpeningStepSize.value;
		}

		nodeColors = PointVector< 3 >( pMesh.nodes() );

		for( unsigned int i=0 ; i<nodeColors.size() ; i++ ) nodeColors[i] = Point< double , 3 >( Random< double >() , Random< double >() , Random< double >() );
	}

	if( Ps.size() )
	{
		MGSolver::Solver< MGSolver::ParallelGaussSeidelRelaxer< 20u > > *licMGSolver = nullptr;
		MGSolver::Solver< MGSolver::ParallelGaussSeidelRelaxer< 20u > > *mgSolver = nullptr;
		{
			Miscellany::NestedTimer timer( "Factorized system" , Verbose.set );
			licMGSolver = new MGSolver::Solver< MGSolver::ParallelGaussSeidelRelaxer< 20u > > ( licL , Ps , false );
			if( Sharpen.value>0 ) mgSolver = new MGSolver::Solver< MGSolver::ParallelGaussSeidelRelaxer< 20u > >( L , Ps , false );
		}
		{
			Miscellany::NestedTimer timer( "Solved mg system" , Verbose.set );
			PointVector< 3 > b;
			b = licM * nodeColors;
			nodeColors = licMGSolver->solve( b , nodeColors , VCycles.value , GSIters.value , GSIters.value , false );
			if( mgSolver )
			{
				b = M * nodeColors + S * nodeColors * SharpeningStepSize.value * Sharpen.value;
				nodeColors = mgSolver->solve( b , nodeColors , VCycles.value , GSIters.value , GSIters.value , false );
			}

			if( mgSolver ) delete mgSolver;
			delete licMGSolver;
		}
	}
	else
	{
		SparseSolver::LLT *licSolver = nullptr;
		SparseSolver::LLT *solver = nullptr;
		{
			Miscellany::NestedTimer timer( "Factorized system" , Verbose.set );
			licSolver = new SparseSolver::LLT( licL );
			if( Sharpen.value>0 ) solver = new SparseSolver::LLT( L );
		}
		{
			Miscellany::NestedTimer timer( "Solved system" , Verbose.set );
			Eigen::MatrixXd x( pMesh.nodes() , 3 ) , b;
			for( unsigned int d=0 ; d<3 ; d++ ) for( unsigned int i=0 ; i<pMesh.nodes() ; i++ ) x(i,d) = nodeColors[i][d];
			b = licM * x;
			x = licSolver->solve( b );
			if( solver )
			{
				b = M * x + S * x * SharpeningStepSize.value * Sharpen.value;
				x = solver->solve( b );
			}
			for( unsigned int d=0 ; d<3 ; d++ ) for( unsigned int i=0 ; i<pMesh.nodes() ; i++ ) nodeColors[i][d] = x(i,d);
		}
		if( solver ) delete solver;
		delete licSolver;
	}
	for( unsigned int i=0 ; i<colors.size() ; i++ )
	{
		unsigned int idx[Degree];
		for( unsigned int d=0 ; d<Degree ; d++ ) idx[d] = i;
		colors[i] = nodeColors[ pMesh.nodeIndex( NodeMultiIndex(idx) ) ];
	}

	return colors;
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set || !InVectorField.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	std::stringstream sStream;
	sStream << "Running Time [Line Integral Convolution (V" << VERSION << ")]";
	Miscellany::NestedTimer nestedTimer( "Running time" , Verbose.set );

	std::vector< std::vector< unsigned int > > polygons;
	std::vector< Vertex > vertices;
	std::vector< Point< double , 3 > > vectorField;

	Factory factory;
	int file_type;
	PLY::ReadPolygons< Factory >( In.value , factory , vertices , polygons , NULL , file_type );
	ReadVectorField( vectorField , InVectorField.value.c_str() );
	if( polygons.size()!=vectorField.size() ) ERROR_OUT( "Polygon/vector-field size mismatch: " , polygons.size() , " != " , vectorField.size() );

	if( Verbose.set ) std::cout << "Vertices / Polygons: " << vertices.size() << " / " << polygons.size() << std::endl;

	std::vector< std::string > comments( 2 );
	Miscellany::Timer timer;
	std::vector< Point< double , 3 > > colors;
	if( Multigrid.set )
	{
		Miscellany::NestedTimer timer( "Computed LIC" , true );
		switch( Degree.value )
		{
			case 1: colors = Execute< 1 , true >( vertices , polygons , vectorField ) ; break;
			case 2: colors = Execute< 2 , true >( vertices , polygons , vectorField ) ; break;
			case 3: colors = Execute< 3 , true >( vertices , polygons , vectorField ) ; break;
			default: ERROR_OUT( "Only degrees 1, 2, and 3 supported" );
		}
	}
	else
	{
		Miscellany::NestedTimer timer( "Computed LIC" , true );
			switch( Degree.value )
		{
			case 1: colors = Execute< 1 , false >( vertices , polygons , vectorField ) ; break;
			case 2: colors = Execute< 2 , false >( vertices , polygons , vectorField ) ; break;
			case 3: colors = Execute< 3 , false >( vertices , polygons , vectorField ) ; break;
			default: ERROR_OUT( "Only degrees 1, 2, and 3 supported" );
		}
	}
	std::stringstream ss;
	ss << "Running time: " << timer.elapsed() << " (s)";
	comments[0] = ss.str();

	ss.str( std::string() );
	ss << "Peak memory usage: " << Miscellany::MemoryInfo::PeakMemoryUsageMB() << " (MB)";
	comments[1] = ss.str();

	auto ClampColor = []( Point< double , 3 > c )
	{
		for( unsigned int i=0 ; i<3 ; i++ ) c[i] = std::max< double >( 0. , std::min< double >( 1. , c[i] ) );
		return c;
	};

	if( Out.set )
	{
		Miscellany::NestedTimer timer( "Wrote file" , true );
		ColorFactory colorFactory;
		std::vector< ColorVertex > colorVertices( vertices.size() );
		for( unsigned int i=0 ; i<vertices.size() ; i++ )
		{
			colorVertices[i].get<0>() = vertices[i];
			colorVertices[i].get<1>() = ClampColor( colors[i] ) * 255.;
		}
		PLY::WritePolygons( Out.value , colorFactory , colorVertices , polygons , file_type , &comments );
	}
	if( Verbose.set ) std::cout << comments[1] << std::endl;

	return EXIT_SUCCESS;
}