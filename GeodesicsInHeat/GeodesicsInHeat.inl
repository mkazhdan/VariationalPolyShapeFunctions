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

#include <map>
#include <set>
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include <Misha/Visualization3D.h>
#include <Misha/SimplexRefinableMesh.h>
#include <Misha/Meshes.h>
#include <Misha/MGSolver.h>
#include <Misha/RightTriangleQuadrature.h>


struct MeshRefiner
{
	MeshRefiner( unsigned int res ) : _res(res){}
	unsigned int subTriangles( void ) const { return _res*_res; }
	unsigned int subVertices( void ) const { return ((_res+2)*(_res+1))/2; }
	unsigned int subTriangles( unsigned int pSize ) const { return subTriangles()*pSize; }
	unsigned int subVertices( unsigned int pSize ) const { return subVertices()*pSize - _res*pSize - (pSize-1); }

	static unsigned int Index( unsigned int i , unsigned int j , unsigned int res )
	{
		// 0 <= i,j < res
		// i+j < res
		// I(i,j) = i + \sum_{k=0}^{j-1}(res-k)
		//        = i + j*R - \sum_{k=0}^{j-1}k
		//        = i + j*R - \sum_{k=1}^j(k-1)
		//        = i + j*R - j*(j-1)/2
		return i + j*res - ( j*(j-1) )/2;
	};

	void setRefinedTriangles( unsigned int sNum , std::vector< SimplexIndex< 2 , unsigned int > > &triangles , std::vector< typename SimplexMesh< 2 >::Sample > &vertexSamples ) const
	{
		triangles.resize( sNum * subTriangles() );
		vertexSamples.resize( sNum * subVertices() );

#pragma omp parallel for
		for( int s=0 ; s<(int)sNum ; s++ ) _addRefinedTriangle( s , subVertices()*s , &triangles[s*subTriangles()] , &vertexSamples[s*subVertices()] );
	}

	void setRefinedPolygons( unsigned int pNum , std::function< unsigned int (unsigned int) > PSize , std::vector< SimplexIndex< 2 , unsigned int > > &triangles , std::vector< typename SimplexMesh< 2 >::Sample > &vertexSamples ) const
	{
		// Cumulative counts of simplex, refined vertex, and refined triangle counts
		std::vector< unsigned int > pSizes(pNum) , vSizes(pNum) , tSizes(pNum);

		auto StartP = [&]( unsigned int p ){ return p==0 ? 0 : pSizes[p-1]; };
		auto StartV = [&]( unsigned int p ){ return p==0 ? 0 : vSizes[p-1]; };
		auto StartT = [&]( unsigned int p ){ return p==0 ? 0 : tSizes[p-1]; };

		for( unsigned int p=0 ; p<pNum ; p++ )
		{
			unsigned int pSize = PSize(p);
			pSizes[p] = pSize , vSizes[p] = subVertices( pSize ) , tSizes[p] = subTriangles( pSize );
			if( p ) pSizes[p] += pSizes[p-1] , vSizes[p] += vSizes[p-1] , tSizes[p] += tSizes[p-1];
		}

		triangles.resize( tSizes.back() );
		vertexSamples.resize( vSizes.back() );

#pragma omp parallel for
		for( int p=0 ; p<(int)pNum ; p++ )
		{
			unsigned int startV = StartV(p) , startT = StartT(p) , startP = StartP(p);
			_addRefinedPolygon( PSize(p) , startP , startV , &triangles[startT] , &vertexSamples[startV] );
		}

	}
protected:
	unsigned int _res;

	void _addRefinedTriangle( unsigned int s , unsigned int vStart , SimplexIndex< 2 , unsigned int > *triangles , typename SimplexMesh< 2 >::Sample *vertexSamples ) const
	{
		for( unsigned int i=0 , t=0 ; i<_res ; i++ ) for( unsigned int j=0 ; j<_res-i ; j++ )
		{
			triangles[t++] = SimplexIndex< 2 , unsigned int >( vStart + Index(i,j,_res+1) , vStart + Index(i+1,j,_res+1) , vStart + Index(i,j+1,_res+1) );
			if( i+j<_res-1 ) triangles[t++] = SimplexIndex< 2 , unsigned int >( vStart + Index(i+1,j+1,_res+1) , vStart + Index(i,j+1,_res+1) , vStart + Index(i+1,j,_res+1) );
		}

		typename SimplexMesh< 2 >::Sample sample;
		sample.sIdx = s;
		for( unsigned int i=0 ; i<_res+1 ; i++ ) for( unsigned int j=0 ; j<_res+1-i ; j++ )
		{
			double x = (i+0.)/_res , y = (j+0.)/_res;
			sample.bcCoordinates = Point< double , 3 >( 1.0-x-y , x , y );
			vertexSamples[ Index(i,j,_res+1) ] = sample;
		}
	}

	struct _WeightedVertex
	{
		int v;
		unsigned int w;

		_WeightedVertex( void ){}
		_WeightedVertex( unsigned int _v , unsigned int _w ) : v(_w==0?-1:_v) , w(_w){}
		bool operator < ( const _WeightedVertex &s ) const
		{
			if( v<s.v ) return true;
			else if( v>s.v ) return false;
			else return w<s.w;
		}
	};
	friend std::ostream &operator << ( std::ostream &os , const _WeightedVertex &v ){ return os << "( " << v.v << " , " << v.w << " )"; }

	struct _TriangleSample
	{
		unsigned int tIdx;
		std::tuple< _WeightedVertex , _WeightedVertex , _WeightedVertex > v;

		_TriangleSample( void ){}
		_TriangleSample( unsigned int i1 , unsigned int i2 , unsigned int i3 , unsigned int w1 , unsigned int w2 , unsigned int w3 )
		{
			_WeightedVertex _v[] = { _WeightedVertex( i1 , w1 ) , _WeightedVertex( i2 , w2 ) , _WeightedVertex( i3 , w3 ) };
			std::sort( _v , _v+3 , []( _WeightedVertex v1 , _WeightedVertex v2 ){ return v1.v<v2.v; } );
			std::get<0>( v ) = _v[0];
			std::get<1>( v ) = _v[1];
			std::get<2>( v ) = _v[2];
		}
		bool operator < ( const _TriangleSample &s ) const { return v<s.v; }
	};
	friend std::ostream &operator << ( std::ostream &os , const _TriangleSample &v ){ return os << "[ " << v.tIdx << " : " << std::get<0>(v.v) << " " << std::get<1>(v.v) << " " << std::get<2>(v.v) << " ]"; }

	void _addRefinedPolygon( unsigned int pSize , unsigned int s , unsigned int vStart , SimplexIndex< 2 , unsigned int > *triangles , typename SimplexMesh< 2 >::Sample *vertexSamples ) const
	{
		// Start by computing the refined geometry, with duplicate vertices
		std::vector< SimplexIndex< 2 , unsigned int > > _triangles( pSize * subTriangles() );
		std::vector< _TriangleSample > _vertexSamples( pSize * subVertices() );

		{
			for( unsigned int p=0 , t=0 ; p<pSize ; p++ )
			{
				unsigned int vStart = p * subVertices();

				// Create the (local) list of triangles
				for( unsigned int i=0 ; i<_res ; i++ ) for( unsigned int j=0 ; j<_res-i ; j++ )
				{
					_triangles[t++] = SimplexIndex< 2 , unsigned int >( vStart + Index(i,j,_res+1) , vStart + Index(i+1,j,_res+1) , vStart + Index(i,j+1,_res+1) );
					if( i+j<_res-1 ) _triangles[t++] = SimplexIndex< 2 , unsigned int >( vStart + Index(i+1,j+1,_res+1) , vStart + Index(i,j+1,_res+1) , vStart + Index(i+1,j,_res+1) );
				}

				// Create the (local) list of vertices
				for( unsigned int i=0 ; i<_res+1 ; i++ ) for( unsigned int j=0 ; j<_res+1-i ; j++ )
				{
#pragma message( "[WARNING] Assuming particular simplex ordering" )
					_TriangleSample sample( pSize , p , (p+1)%pSize , _res-i-j , i , j );
					sample.tIdx = p;
					_vertexSamples[ vStart + Index(i,j,_res+1) ] = sample;
				}
			}
		}

		// Identify the duplicate triangle samples
		std::map< _TriangleSample , unsigned int > sampleMap;
		for( unsigned int i=0 ; i<_vertexSamples.size() ; i++ ) sampleMap[ _vertexSamples[i] ] = i;

		// Assign the unique triangle samples an index, compute the associated simplex sample, and set the correct triangle indices
		{
			unsigned int idx = 0;
			for( auto &[key,value] : sampleMap )
			{
				typename SimplexMesh< 2 >::Sample sample;
				auto SetBC = [&]( unsigned int tIdx , _WeightedVertex wv )
				{
					if( wv.w )
					{
						double coord = (double)wv.w/_res;
						if( wv.v==pSize )     sample.bcCoordinates[2] = coord;
						else if( wv.v==tIdx ) sample.bcCoordinates[0] = coord;
						else                  sample.bcCoordinates[1] = coord;
					}
				};
				SetBC( key.tIdx , std::get<0>(key.v) );
				SetBC( key.tIdx , std::get<1>(key.v) );
				SetBC( key.tIdx , std::get<2>(key.v) );

				sample.sIdx = s + key.tIdx;
				vertexSamples[idx] = sample;
				value = idx++;
			}
			for( unsigned int t=0 ; t<_triangles.size() ; t++ ) for( unsigned int j=0 ; j<3 ; j++ ) triangles[t][j] = vStart + sampleMap[ _vertexSamples[ _triangles[t][j] ] ];
		}
	}
};

template< unsigned int Degree , bool Hierarchical >
struct GeodesicsInHeatVisualization : public Misha::Viewable3D< GeodesicsInHeatVisualization< Degree , Hierarchical > >
{
	typedef typename Misha::Viewable< GeodesicsInHeatVisualization >::KeyboardCallBack KeyboardCallBack;
	typedef typename KeyboardCallBack::Modifiers KeyboardCallBackModifiers;
	using Misha::Viewable< GeodesicsInHeatVisualization >::promptCallBack;

	static const unsigned int Dim = 2;
	typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	using Misha::Viewable< GeodesicsInHeatVisualization >::info;
	using Misha::Viewable< GeodesicsInHeatVisualization >::callBacks;
	using Misha::Viewable< GeodesicsInHeatVisualization >::screenWidth;
	using Misha::Viewable< GeodesicsInHeatVisualization >::screenHeight;
	using Misha::Viewable3D< GeodesicsInHeatVisualization >::camera;

	enum
	{
		EDGES_NONE ,
		EDGES_POLYGON ,
		EDGES_TRIANGLES ,
		EDGES_COUNT
	};
	unsigned int edgeMode;
	bool lightOn , showColor , showGeodesic , useStripeTexture;
	static const unsigned int TEXTURE_WIDTH = 4096;
	unsigned int stripeShift;

	GeodesicsInHeatVisualization( const std::vector< Point< double , 3 > > &vertices , const std::vector< std::vector< unsigned int > > &polygons , unsigned int coarseNodeDim , double diffusionTime , double stiffnessRegularizer , unsigned int subdivisionIterations , unsigned int width , unsigned int height );
	~GeodesicsInHeatVisualization( void );

	void draw3D( void );
	void outputMesh( std::string fileName ) const;
	void idle( void );
	void mouseFunc( int button , int state , int x , int y );

	static void DecreaseStripes( GeodesicsInHeatVisualization* v , const char * )
	{
		if( (1<<v->stripeShift)<TEXTURE_WIDTH )
		{
			v->stripeShift++;
			if( v->_stripeTextureHandle ) glDeleteTextures( 1 , &v->_stripeTextureHandle );
			v->_stripeTextureHandle = 0;
		}
	}
	static void IncreaseStripes( GeodesicsInHeatVisualization* v , const char * )
	{
		if( v->stripeShift>0 )
		{
			v->stripeShift--;
			if( v->_stripeTextureHandle ) glDeleteTextures( 1 , &v->_stripeTextureHandle );
			v->_stripeTextureHandle = 0;
		}
	}
	static void AdvanceIteration( GeodesicsInHeatVisualization* v , const char * )
	{
		v->_iterate = false;
		if( Hierarchical && v->_sourceVertex!=-1 )
		{
			v->_update();
			glutPostRedisplay();
		}
	}
	static void ToggleIteration( GeodesicsInHeatVisualization* v , const char * ){ v->_iterate=!v->_iterate; }
	static void ToggleStripeTexture( GeodesicsInHeatVisualization* v , const char * ){ v->useStripeTexture=!v->useStripeTexture; }
	static void ToggleColor( GeodesicsInHeatVisualization* v , const char * ){ v->showColor=!v->showColor; }
	static void CycleEdges( GeodesicsInHeatVisualization* v , const char * ){ v->edgeMode=(v->edgeMode+1)%EDGES_COUNT; }
	static void ToggleGeodesic( GeodesicsInHeatVisualization* v , const char * ){ v->showGeodesic=!v->showGeodesic; }
	static void ToggleLight( GeodesicsInHeatVisualization* v , const char * ){ v->lightOn=!v->lightOn; }
	static void SetDiffusionTime( GeodesicsInHeatVisualization* v , const char *prompt );
	static void SetSourceVertex( GeodesicsInHeatVisualization* v , const char *prompt );
	static void SetQuadratureSamples( GeodesicsInHeatVisualization* v , const char *prompt );
	static void SetRefinementResolution( GeodesicsInHeatVisualization* v , const char *prompt );
	static void ReadCamera( GeodesicsInHeatVisualization* v , const char *prompt );
	static void WriteCamera( GeodesicsInHeatVisualization* v , const char *prompt );
	static void OutputMesh( GeodesicsInHeatVisualization* v , const char *prompt ){ if( prompt ) v->outputMesh( std::string( prompt ) ); }

protected:
	typedef typename std::conditional< Hierarchical , HierarchicalSimplexRefinableCellMesh< Dim , Degree > , SimplexRefinableCellMesh< Dim , Degree > >::type SimplexRefinableMesh;
	typedef typename std::conditional< Hierarchical , MGSolver::Solver< MGSolver::ParallelGaussSeidelRelaxer< 20u > > , SparseSolver::LLT >::type SolverType;

	bool _verbose , _iterate;
	GLuint _valueTextureHandle , _stripeTextureHandle;
	unsigned int _sourceVertex;
	unsigned int _vCycles , _gsIters;
	unsigned int _vNum;
	unsigned int _quadratureSamples;
	double _scale;
	double _diffusionTime;
	unsigned int _diffusionTimeLine , _meshSizeLine , _sourceVertexLine , _quadratureSamplesLine , _vCyclesLine;
	unsigned int _systemTimeLine , _constraintTimeLine , _solutionTimeLine;

	enum
	{
		VBO_VERTEX ,
		VBO_TRIANGLE ,
		VBO_DIFFUSED_DELTA_TEXTURE ,
		VBO_GEODESIC_TEXTURE ,
		VBO_COUNT
	};
	GLuint *_vbo;
	std::vector< Point3D< double > > _verticesAndNormals;
	std::vector< double > _diffusedDeltaTextureCoordinates , _geodesicTextureCoordinates;
	std::vector< SimplexIndex< 2 , unsigned int > > _triangles;
	Eigen::SparseMatrix< double > _eMatrix;

	Meshes::PolygonMesh< unsigned int > _polygonMesh;
	SimplexRefinableMesh _simplexRefinableCellMesh;
	std::vector< std::vector< unsigned int > > _polygons;
	std::vector< Point3D< double > > _fullVertices , _polygonNormals;
	std::vector< unsigned int > _vertexToFineNode , _triangleToPolygon;
	Point< Polynomial::Polynomial< Dim , Degree-1 , double > , 2 > _dElements[ SimplexElements< Dim , Degree >::NodeNum ];

	SolverType *_diffusionSolver , *_geodesicSolver;
	Eigen::SparseMatrix< double >  _M , _S , _L , _Pt , _P;
	std::vector< Eigen::SparseMatrix< double > > __P;
	Eigen::VectorXd _coarseDiffusedDeltaB , _coarseDiffusedDeltaX , _coarseGeodesicB , _coarseGeodesicX;
	Eigen::VectorXd _fineDiffusedDeltaB , _fineDiffusedDeltaX , _fineGeodesicB , _fineGeodesicX;

	void _setSamplingMesh( unsigned int res );
	void _drawMesh( bool useColor );
	void _setVBO( void );

	void _setDiffusionTime( void );
	void _setSource( void );
	void _setDiffuseDeltaB( void );
	void _setDiffuseDeltaX( void );
	template< unsigned int Samples > void _setGeodesicB( void );
	void _setGeodesicX( void );
	void _update( void );
};

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::SetDiffusionTime( GeodesicsInHeatVisualization *v , const char *prompt )
{
	if( prompt )
	{
		v->_diffusionTime = atof(prompt);
		v->_setDiffusionTime();
		glutPostRedisplay();
	}
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::SetRefinementResolution( GeodesicsInHeatVisualization* v , const char *prompt )
{
	int res = atoi(prompt);
	if( res<=0 ) WARN( "Resolution must be positive: " , res );
	else
	{
		v->_setSamplingMesh( res );

		v->_setVBO();

		if( v->_sourceVertex!=-1 )
		{
			{
				Eigen::VectorXd eValues = v->_eMatrix * v->_fineDiffusedDeltaX;
				glBindBuffer( GL_ARRAY_BUFFER , v->_vbo[VBO_DIFFUSED_DELTA_TEXTURE] );
				glBufferSubData( GL_ARRAY_BUFFER , 0 , eValues.size() * sizeof( double ) , &eValues[0] );
				glBindBuffer( GL_ARRAY_BUFFER , 0 );
			}
			{
				Eigen::VectorXd eValues = v->_eMatrix * v->_fineGeodesicX;
				glBindBuffer( GL_ARRAY_BUFFER , v->_vbo[VBO_GEODESIC_TEXTURE] );
				glBufferSubData( GL_ARRAY_BUFFER , 0 , eValues.size() * sizeof( double ) , &eValues[0] );
				glBindBuffer( GL_ARRAY_BUFFER , 0 );
			}
		}
		glutPostRedisplay();
	}
}


template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::SetQuadratureSamples( GeodesicsInHeatVisualization *v , const char *prompt )
{
	if( prompt )
	{
		unsigned int s = atoi( prompt );
		switch( s )
		{
		case  1:
		case  3:
		case  4:
		case  6:
		case  7:
		case 12:
		case 13:
		case 24:
		case 27:
		case 32:
			v->_quadratureSamples = s;
			v->_setSource();
			glutPostRedisplay();
			break;
		default:
			std::cerr << "Unsupported number of integration samples: " << s << std::endl;
			std::cerr << "\tValid options are 1, 3, 4, 6, 7, 12, 13, 24, 27, or 32" << std::endl;
		}
	}
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::SetSourceVertex( GeodesicsInHeatVisualization *v , const char *prompt )
{
	if( prompt )
	{
		v->_sourceVertex = atoi(prompt);
		if( v->_sourceVertex>=v->_vNum )
		{
			WARN( "Source index excceeds number of vertices: " , v->_sourceVertex ,  " >= " , v->_vNum );
			v->_sourceVertex = -1;
		}
		v->_setSource();
		glutPostRedisplay();
	}
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::WriteCamera( GeodesicsInHeatVisualization *v , const char *prompt )
{
	FILE *fp = fopen( prompt , "w" );
	if( !fp ) WARN( "Failed to open file for writing: " , std::string( prompt ) );
	else
	{
		v->camera().write( fp );
		fclose( fp );
	}
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::ReadCamera( GeodesicsInHeatVisualization *v , const char *prompt )
{
	FILE *fp = fopen( prompt , "r" );
	if( !fp ) WARN( "Failed to open file for reading: " , std::string( prompt ) );
	else
	{
		v->camera().read( fp );
		fclose( fp );
		glutPostRedisplay();
	}
}

template< unsigned int Degree , bool Hierarchical >
GeodesicsInHeatVisualization< Degree , Hierarchical >::~GeodesicsInHeatVisualization( void )
{
	delete _diffusionSolver;
	delete _geodesicSolver;
	delete[] _vbo;
}

template< unsigned int Degree , bool Hierarchical >
GeodesicsInHeatVisualization< Degree , Hierarchical >::GeodesicsInHeatVisualization( const std::vector< Point3D< double > > &vertices , const std::vector< std::vector< unsigned int > > &polygons , unsigned int coarseNodeDim , double diffusionTime , double stiffnessRegularizer , unsigned int subdivisionIterations , unsigned int width , unsigned int height )
	: _diffusionTime(diffusionTime) , edgeMode(EDGES_NONE) , lightOn(true) , showColor(true) , _gsIters(2) , _diffusionSolver(NULL) , _geodesicSolver(NULL) , _verbose(false) , _sourceVertex(-1) , _valueTextureHandle(0) , _stripeTextureHandle(0) , showGeodesic(true) , useStripeTexture(true) , stripeShift(5) , _quadratureSamples(12) ,
	_vbo(NULL) , _iterate(false)
{
	// Set up the 3D visualization stuff
	{
		Point3D< double > bBox[2];
		bBox[0] = bBox[1] = vertices[0];
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) for( unsigned int d=0 ; d<3 ; d++ )
		{
			bBox[0][d] = std::min< double >( bBox[0][d] , vertices[i][d] );
			bBox[1][d] = std::max< double >( bBox[1][d] , vertices[i][d] );
		}
		_scale = Point3D< double >::Length( bBox[1]-bBox[0] );

		Misha::Viewable3D< GeodesicsInHeatVisualization< Degree , Hierarchical > >::init( bBox , width , height );
	}

	// Set up visualizaiton stuff
	{
		callBacks.push_back( KeyboardCallBack( this , '2' , KeyboardCallBackModifiers() , "read camera" , "Camera" , ReadCamera ) );
		callBacks.push_back( KeyboardCallBack( this , '1' , KeyboardCallBackModifiers() , "write camera" , "Camera" , WriteCamera ) );
		callBacks.push_back( KeyboardCallBack( this , '[' , KeyboardCallBackModifiers() , "decrease stripes" , DecreaseStripes ) );
		callBacks.push_back( KeyboardCallBack( this , ']' , KeyboardCallBackModifiers() , "increase stripes" , IncreaseStripes ) );
		callBacks.push_back( KeyboardCallBack( this , 'n' , KeyboardCallBackModifiers() , "set quadrature samples" , "Quadrature samples" , SetQuadratureSamples ) );
		callBacks.push_back( KeyboardCallBack( this , 'r' , KeyboardCallBackModifiers() , "set refinement resolution" , "Refinement resolution" , SetRefinementResolution ) );
		callBacks.push_back( KeyboardCallBack( this , 't' , KeyboardCallBackModifiers() , "set diffusion time" , "Diffusion time" , SetDiffusionTime ) );
		callBacks.push_back( KeyboardCallBack( this , 'v' , KeyboardCallBackModifiers() , "set source vertex" , "Source vertex" , SetSourceVertex ) );
		callBacks.push_back( KeyboardCallBack( this , '+' , KeyboardCallBackModifiers() , "advance iteration" , AdvanceIteration ) );
		callBacks.push_back( KeyboardCallBack( this , ' ' , KeyboardCallBackModifiers() , "toggle iteration" , ToggleIteration ) );
		callBacks.push_back( KeyboardCallBack( this , 'T' , KeyboardCallBackModifiers() , "toggle stripes" , ToggleStripeTexture ) );
		callBacks.push_back( KeyboardCallBack( this , 'e' , KeyboardCallBackModifiers() , "cycle edges" , CycleEdges ) );
		callBacks.push_back( KeyboardCallBack( this , 'g' , KeyboardCallBackModifiers() , "toggle Geodesic" , ToggleGeodesic ) );
		callBacks.push_back( KeyboardCallBack( this , 'c' , KeyboardCallBackModifiers() , "toggle color" , ToggleColor ) );
		callBacks.push_back( KeyboardCallBack( this , 'l' , KeyboardCallBackModifiers() , "toggle light" , ToggleLight ) );
		callBacks.push_back( KeyboardCallBack( this , 'o' , KeyboardCallBackModifiers() , "output mesh" , "File name" , OutputMesh ) );

		if( Hierarchical ) _vCyclesLine = (int)info.size() , info.push_back( new char[512] );
		_quadratureSamplesLine = (int)info.size() , info.push_back( new char[512] );
		_sourceVertexLine      = (int)info.size() , info.push_back( new char[512] );
		_diffusionTimeLine     = (int)info.size() , info.push_back( new char[512] );
		_solutionTimeLine      = (int)info.size() , info.push_back( new char[512] );
		_constraintTimeLine    = (int)info.size() , info.push_back( new char[512] );
		_systemTimeLine        = (int)info.size() , info.push_back( new char[512] );
		_meshSizeLine          = (int)info.size() , info.push_back( new char[512] );
	}

	// Set up the system stuff
	_vNum = (unsigned int)vertices.size();
	_polygons = polygons;
	_fullVertices.resize( vertices.size() + polygons.size() );

	for( unsigned int n=0 ; n<SimplexElements< Dim , Degree >::NodeNum ; n++ ) _dElements[n] = SimplexElements< Dim , Degree >::Differential( n );

	SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
	eWeights.kWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1;

	_polygonMesh = Meshes::PolygonMesh< unsigned int >( _polygons );
	std::function< Point< double , 3 > ( unsigned int ) > vertexPositionFunction = [&]( unsigned int v ){ return Point3D< double >( vertices[v] ); };
	std::function< Point< double , 3 > ( unsigned int ) > fullVertexPositionFunction = _polygonMesh.fullVertexPositionFunction( vertexPositionFunction , true );
	for( unsigned int i=0 ; i<_fullVertices.size() ; i++ ) _fullVertices[i] = fullVertexPositionFunction( i );

	_polygonNormals.resize( _polygons.size() );
	for( unsigned int i=0 ; i<_polygons.size() ; i++ )
	{
		const std::vector< unsigned int > &poly = _polygons[i];
		unsigned int sz = (unsigned int)poly.size();
		Point3D< double > n;
		for( unsigned int v=0 ; v<sz ; v++ ) n += Point3D< double >::CrossProduct( _fullVertices[ poly[v] ] , _fullVertices[ poly[(v+1)%sz] ] );
		_polygonNormals[i] = n / Point3D< double >::Length( n );
	}

	if constexpr( Hierarchical ) _simplexRefinableCellMesh = _polygonMesh.template hierarchicalSimplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarseNodeDim , false , false , 0 , true );
	else                         _simplexRefinableCellMesh = _polygonMesh.template             simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarseNodeDim , false , false , 0 , true );
	_simplexRefinableCellMesh.hashSimplexMeshLocalToGlobalNodeIndex();

	_simplexRefinableCellMesh.makeUnitVolume();

	const SimplexMesh< Dim , Degree > &simplexMesh = _simplexRefinableCellMesh.simplexMesh();

	_vertexToFineNode.resize( _fullVertices.size() );
	for( unsigned int v=0 ; v<_fullVertices.size() ; v++ )
	{
		unsigned int idx[Degree];
		for( unsigned int d=0 ; d<Degree ; d++ ) idx[d] = v;
		_vertexToFineNode[v] = simplexMesh.nodeIndex( NodeMultiIndex(idx) );
	}

	static const unsigned int NodeNum = SimplexElements< Dim , Degree >::NodeNum;

	_triangleToPolygon.resize( simplexMesh.simplices() );
	for( unsigned int p=0 , t=0 ; p<polygons.size() ; p++ ) for( unsigned int i=0 ; i<polygons[p].size() ; i++ ,t++ ) _triangleToPolygon[t] = p;

	_coarseDiffusedDeltaX = Eigen::VectorXd::Zero( _simplexRefinableCellMesh.nodes() );
	_coarseDiffusedDeltaB = Eigen::VectorXd::Zero( _simplexRefinableCellMesh.nodes() );
	_coarseGeodesicX = Eigen::VectorXd::Zero( _simplexRefinableCellMesh.nodes() );
	_coarseGeodesicB = Eigen::VectorXd::Zero( _simplexRefinableCellMesh.nodes() );
	_fineDiffusedDeltaX = Eigen::VectorXd::Zero( simplexMesh.nodes() );
	_fineDiffusedDeltaB = Eigen::VectorXd::Zero( simplexMesh.nodes() );
	_fineGeodesicX = Eigen::VectorXd::Zero( simplexMesh.nodes() );
	_fineGeodesicB = Eigen::VectorXd::Zero( simplexMesh.nodes() );

	if constexpr( Hierarchical )
	{
		_P = _simplexRefinableCellMesh.P( _simplexRefinableCellMesh.maxLevel() , coarseNodeDim );
		__P.resize( coarseNodeDim );
		for( unsigned int d=0 ; d<coarseNodeDim ; d++ ) __P[d] = _simplexRefinableCellMesh.P( d+1 , d );
	}
	else _P = _simplexRefinableCellMesh.P();
	_Pt = _P.transpose();
	_M = _Pt * simplexMesh.mass() * _P;
	_S = _Pt * simplexMesh.stiffness() * _P;

	Miscellany::Timer timer;
	if constexpr( Hierarchical )
	{
		_geodesicSolver = new SolverType( _S + _M * stiffnessRegularizer , __P , false );
		_geodesicSolver->state.vCycles = 1;
		_geodesicSolver->state.rIters  = 0;
		_geodesicSolver->state.pIters  = _gsIters;
		_geodesicSolver->state.verbose = _verbose;
	}
	else
	{
		_geodesicSolver = new SolverType( _S + _M * stiffnessRegularizer );
		switch( _geodesicSolver->info() )
		{
			case Eigen::NumericalIssue: ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- numerical issue" );
			case Eigen::NoConvergence:  ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- no convergence" );
			case Eigen::InvalidInput:   ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- invalid input" );
			case Eigen::Success: ;
		}
	}

	_setDiffusionTime();

	_setSamplingMesh( 1<<subdivisionIterations );

	if( Hierarchical ) sprintf( info[ _vCyclesLine ] , "V-Cycles:" );
	sprintf( info[ _solutionTimeLine ] , "Solve time:" );
	sprintf( info[ _constraintTimeLine ] , "Constraint time:" );
	sprintf( info[ _systemTimeLine ] , "System time: %.2f(s)" , timer.elapsed() );
	sprintf( info[_meshSizeLine] , "Nodes: %d" , (int)_simplexRefinableCellMesh.nodes() );
	sprintf( info[_sourceVertexLine] , "Source vertex: %d\n" , (int)_sourceVertex );
	sprintf( info[_quadratureSamplesLine] , "Quadrature samples: %d\n" , (int)_quadratureSamples );
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_setSamplingMesh( unsigned int res )
{
	const SimplexMesh< Dim , Degree > &simplexMesh = _simplexRefinableCellMesh.simplexMesh();

	Miscellany::Timer timer;
	MeshRefiner refiner(res);
	std::vector< typename SimplexMesh< 2 >::Sample > vertexSamples;

	refiner.setRefinedPolygons( (unsigned int)_polygons.size() , [&]( unsigned int p ){ return (unsigned int)_polygons[p].size(); } , _triangles , vertexSamples );
	_verticesAndNormals.resize( vertexSamples.size()*2 );

#pragma omp parallel for
	for( int i=0 ; i<(int)vertexSamples.size() ; i++ )
	{
		SimplexIndex< 2 , unsigned int > tri = simplexMesh.simplex( vertexSamples[i].sIdx );
		_verticesAndNormals[i] = Point3D< double >();
		for( unsigned int j=0 ; j<3 ; j++ ) _verticesAndNormals[i] += _fullVertices[ tri[j] ] * vertexSamples[i].bcCoordinates[j];
		_verticesAndNormals[ i+vertexSamples.size() ] = _polygonNormals[ _triangleToPolygon[ vertexSamples[i].sIdx ] ];
	}
	_eMatrix = simplexMesh.evaluationMatrix( vertexSamples );

	std::cout << "Set sampling mesh (" << res << "): " << timer.elapsed() << std::endl;
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_setDiffusionTime( void )
{
	_L = _M + _S * _diffusionTime;

	if constexpr( Hierarchical )
	{
		if( _diffusionSolver ) delete _diffusionSolver;
		_diffusionSolver = new SolverType( _L , __P , false );
		_diffusionSolver->state.vCycles = 1;
		_diffusionSolver->state.rIters  = 0;
		_diffusionSolver->state.pIters  = _gsIters;
		_diffusionSolver->state.verbose = _verbose;
	}
	else
	{
		if( !_diffusionSolver )
		{
			_diffusionSolver = new SolverType();
			_diffusionSolver->analyzePattern( _L );
		}
		_diffusionSolver->factorize( _L );

		switch( _diffusionSolver->info() )
		{
			case Eigen::NumericalIssue: ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- numerical issue" );
			case Eigen::NoConvergence:  ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- no convergence" );
			case Eigen::InvalidInput:   ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- invalid input" );
			case Eigen::Success: ;
		}
	}
	_setSource();

	sprintf( info[ _diffusionTimeLine ] , "Diffusion time: %.2e" , _diffusionTime );
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_setDiffuseDeltaB( void )
{
	for( unsigned int i=0 ; i<_coarseDiffusedDeltaX.size() ; i++ ) _coarseDiffusedDeltaX[i] = 0;
	for( unsigned int i=0 ; i<_fineDiffusedDeltaB.size() ; i++ ) _fineDiffusedDeltaB[i] = 0;
	_fineDiffusedDeltaB[ _vertexToFineNode[ _sourceVertex ] ] = 1;
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_setDiffuseDeltaX( void )
{
	if constexpr( Hierarchical ) _coarseDiffusedDeltaX = _diffusionSolver->solve( (Eigen::VectorXd)( _Pt * _fineDiffusedDeltaB ) , _coarseDiffusedDeltaX );
	else                         _coarseDiffusedDeltaX = _diffusionSolver->solve( _Pt * _fineDiffusedDeltaB );
	_fineDiffusedDeltaX = _P * _coarseDiffusedDeltaX;

	// Normalize the diffused delta (which will only affect the gradients by scaling, which will be normalized out in any case)
	double min = std::numeric_limits< double >::infinity() , max = -std::numeric_limits< double >::infinity();
	for( unsigned int i=0 ; i<_fineDiffusedDeltaX.size() ; i++ ) min = std::min< double >( min , _fineDiffusedDeltaX[i] ) , max = std::max< double >( max , _fineDiffusedDeltaX[i] );
	for( unsigned int i=0 ; i<_fineDiffusedDeltaX.size() ; i++ ) _fineDiffusedDeltaX[i] = ( _fineDiffusedDeltaX[i] - min ) / ( max - min );

	if( _vbo )
	{
		Eigen::VectorXd eValues = _eMatrix * _fineDiffusedDeltaX;
		glBindBuffer( GL_ARRAY_BUFFER , _vbo[VBO_DIFFUSED_DELTA_TEXTURE] );
		glBufferSubData( GL_ARRAY_BUFFER , 0 , eValues.size() * sizeof( double ) , &eValues[0] );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );
	}
	else WARN( "vbo has not been set" );
}

template< unsigned int Degree , bool Hierarchical >
template< unsigned int Samples >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_setGeodesicB( void )
{
	static const unsigned int NodeNum = SimplexElements< Dim , Degree >::NodeNum;
	typedef TriangleIntegrator< Samples > Integrator;

	Point2D< double > d[ NodeNum ][ Samples ];
	for( unsigned int n=0 ; n<NodeNum ; n++ ) for( unsigned int s=0 ; s<Samples ; s++ )
		d[n][s] = Point2D< double >( _dElements[n][0]( Integrator::Positions[s] ) , _dElements[n][1]( Integrator::Positions[s] ) );

	for( unsigned int i=0 ; i<_coarseGeodesicX.size() ; i++ ) _coarseGeodesicX[i] = 0;
	for( unsigned int i=0 ; i<_fineGeodesicB.size() ; i++ ) _fineGeodesicB[i] = 0;

	const SimplexMesh< Dim , Degree > &simplexMesh = _simplexRefinableCellMesh.simplexMesh();

	static std::vector< std::vector< double > > __fineGeodesicB( omp_get_max_threads() );
	for( unsigned int t=0 ; t<__fineGeodesicB.size() ; t++ )
	{
		__fineGeodesicB[t].resize( _fineGeodesicB.size() );
#pragma omp parallel for
		for( int i=0 ; i<(int)__fineGeodesicB[t].size() ; i++ ) __fineGeodesicB[t][i] = 0;
	}

#pragma omp parallel for
	for( int s=0 ; s<(int)simplexMesh.simplices() ; s++ )
	{
		Point2D< double > _d[ Samples ];
		SquareMatrix< double , 2 > g = simplexMesh.metric( s );
		SquareMatrix< double , 2 > gInv = g.inverse();
		double area = sqrt( g.determinant() )/2.;

		unsigned nodeIndices[NodeNum];
		for( unsigned int n=0 ; n<NodeNum ; n++ ) nodeIndices[n] = simplexMesh.nodeIndex(s,n);

		// Get the differential
		for( unsigned int i=0 ; i<Samples ; i++ )
		{
			_d[i] = Point2D< double >();
			for( unsigned int n=0 ; n<NodeNum ; n++ ) _d[i] += d[n][i] * _fineDiffusedDeltaX[ nodeIndices[n] ];
		}

		// Normalize, dualize, and quadrature-weight
		for( unsigned int i=0 ; i<Samples ; i++ )
		{
			_d[i] = gInv * _d[i];
			double l = Point2D< double >::Dot( g * _d[i] , _d[i] );
			if( l>0 ) _d[i] = _d[i] * Integrator::Weights[i] * area / sqrt(l);
			else _d[i] = Point2D< double >();
		}

		std::vector< double > &_fineGeodesicB = __fineGeodesicB[ omp_get_thread_num() ];
		for( unsigned int n=0 ; n<NodeNum ; n++ ) for( unsigned int i=0 ; i<Samples ; i++ )
			_fineGeodesicB[ nodeIndices[n] ] -= Point2D< double >::Dot( d[n][i] , _d[i] );
	}
	for( unsigned int t=0 ; t<__fineGeodesicB.size() ; t++ ) for( unsigned int i=0 ; i<_fineGeodesicB.size() ; i++ ) _fineGeodesicB[i] += __fineGeodesicB[t][i];

	sprintf( info[_quadratureSamplesLine] , "Quadrature samples: %d\n" , (int)_quadratureSamples );
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_setGeodesicX( void )
{
	if constexpr( Hierarchical ) _coarseGeodesicX = _geodesicSolver->solve( (Eigen::VectorXd)( _Pt * _fineGeodesicB ) , _coarseGeodesicX );
	else                         _coarseGeodesicX = _geodesicSolver->solve( _Pt * _fineGeodesicB );
	_fineGeodesicX = _P * _coarseGeodesicX;

	double sourceEDT = _fineGeodesicX[ _vertexToFineNode[ _sourceVertex ] ];
	for( unsigned int i=0 ; i<_fineGeodesicX.size() ; i++ ) _fineGeodesicX[i] -= sourceEDT;
	double min = std::numeric_limits< double >::infinity() , max = -std::numeric_limits< double >::infinity();
	for( unsigned int i=0 ; i<_fineGeodesicX.size() ; i++ ) min = std::min< double >( min , _fineGeodesicX[i] ) , max = std::max< double >( max , _fineGeodesicX[i] );
	for( unsigned int i=0 ; i<_fineGeodesicX.size() ; i++ ) _fineGeodesicX[i] /= max;

	if( _vbo )
	{
		Eigen::VectorXd eValues = _eMatrix * _fineGeodesicX;
		glBindBuffer( GL_ARRAY_BUFFER , _vbo[VBO_GEODESIC_TEXTURE] );
		glBufferSubData( GL_ARRAY_BUFFER , 0 , eValues.size() * sizeof( double ) , &eValues[0] );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );
	}
	else WARN( "vbo has not been set" );
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_update( void )
{
	double constraintTime=0 , solveTime=0;
	Miscellany::Timer timer;

	// Solve for the diffused delta
	timer.reset();
	_setDiffuseDeltaX();
	solveTime += timer.elapsed();

	// Set the geodesic constraints
	timer.reset();
	switch( _quadratureSamples )
	{
		case  1: _setGeodesicB<  1 >() ; break;
		case  3: _setGeodesicB<  3 >() ; break;
		case  4: _setGeodesicB<  4 >() ; break;
		case  6: _setGeodesicB<  6 >() ; break;
		case  7: _setGeodesicB<  7 >() ; break;
		case 12: _setGeodesicB< 12 >() ; break;
		case 13: _setGeodesicB< 13 >() ; break;
		case 24: _setGeodesicB< 24 >() ; break;
		case 27: _setGeodesicB< 27 >() ; break;
		case 32: _setGeodesicB< 32 >() ; break;
		default: WARN( "Number of quadrature samples not supported: " , _quadratureSamples );
	}
	constraintTime += timer.elapsed();

	// Solve for the geodesic
	timer.reset();
	_setGeodesicX();
	solveTime += timer.elapsed();

	if( Hierarchical ) sprintf( info[ _vCyclesLine ] , "V-Cycles: %d" , ++_vCycles );
	sprintf( info[ _solutionTimeLine ] , "Solve time: %.2f(s)" , solveTime );
	sprintf( info[ _constraintTimeLine ] , "Constraint time: %.2f(s)" , constraintTime );
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_setSource( void )
{

	if( _sourceVertex!=-1 )
	{
		_vCycles = 0;
		_setDiffuseDeltaB();
		_update();
	}
	sprintf( info[_sourceVertexLine] , "Source vertex: %d\n" , (int)_sourceVertex );
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::idle( void )
{
	if( !promptCallBack )
	{ 
		if( Hierarchical && _iterate && _sourceVertex!=-1 )
		{
			_update();
			glutPostRedisplay();
		}
	}
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::mouseFunc( int button , int state , int x , int y )
{
	if( button==GLUT_LEFT_BUTTON && state==GLUT_DOWN && ( glutGetModifiers() & GLUT_ACTIVE_SHIFT ) )
	{
		_sourceVertex = -1;
		Point3D< double > p;
		if( Misha::Viewable3D< GeodesicsInHeatVisualization< Degree , Hierarchical > >::select( x , screenHeight-1-y , p ) )
		{
			double d = std::numeric_limits< double >::infinity();
			for( unsigned int i=0 ; i<_vNum ; i++ ) 
			{
				double _d = Point3D< double >::SquareNorm( p - _fullVertices[i] );
				if( _d<d ) d = _d , _sourceVertex = i;
			}
		}
		_setSource();
		glutPostRedisplay();
	}
	else Misha::Viewable3D< GeodesicsInHeatVisualization< Degree , Hierarchical > >::mouseFunc( button , state , x , y );
}


template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_setVBO( void )
{
	if( _vbo ) glDeleteBuffers( VBO_COUNT , _vbo );
	else _vbo = new GLuint[ VBO_COUNT ];

	glGenBuffers( VBO_COUNT , _vbo );

	glBindBuffer( GL_ARRAY_BUFFER , _vbo[VBO_VERTEX] );
	glBufferData( GL_ARRAY_BUFFER , _verticesAndNormals.size() * sizeof( Point3D< double > ) , &_verticesAndNormals[0] , GL_STATIC_DRAW );

	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _vbo[VBO_TRIANGLE] );
	glBufferData( GL_ELEMENT_ARRAY_BUFFER , _triangles.size() * sizeof( int ) * 3 , &_triangles[0] , GL_STATIC_DRAW );

	std::vector< double > textureCoordinates( _verticesAndNormals.size()/2 , 0.75 );

	glBindBuffer( GL_ARRAY_BUFFER , _vbo[VBO_DIFFUSED_DELTA_TEXTURE] );
	glBufferData( GL_ARRAY_BUFFER , textureCoordinates.size() * sizeof( double ) , &textureCoordinates[0] , GL_DYNAMIC_DRAW );

	glBindBuffer( GL_ARRAY_BUFFER , _vbo[VBO_GEODESIC_TEXTURE] );
	glBufferData( GL_ARRAY_BUFFER , textureCoordinates.size() * sizeof( double ) , &textureCoordinates[0] , GL_DYNAMIC_DRAW );
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::outputMesh( std::string fileName ) const
{
	std::vector< SimplexIndex< 2 , int > > triangles( _triangles.size() );
	for( unsigned int i=0 ; i<_triangles.size() ; i++ ) for( unsigned int j=0 ; j<3 ; j++ ) triangles[i][j] = (int)_triangles[i][j];
	if( _sourceVertex==-1 )
	{
		typedef VertexFactory::Factory< double , VertexFactory::PositionFactory< double , 3 > , VertexFactory::NormalFactory< double , 3 > > Factory;
		typedef typename Factory::VertexType Vertex;

		std::vector< Vertex > vertices( _verticesAndNormals.size()/2 );
		for( unsigned int i=0 ; i<vertices.size() ; i++ )
		{
			vertices[i].get<0>() = _verticesAndNormals[i];
			vertices[i].get<1>() = _verticesAndNormals[ vertices.size() + i ];
		}
		Factory factory;
		PLY::WriteTriangles< Factory >( fileName , factory , vertices , triangles , PLY_BINARY_NATIVE );
	}
	else
	{
		typedef VertexFactory::Factory< double , VertexFactory::PositionFactory< double , 3 > , VertexFactory::NormalFactory< double , 3 > , VertexFactory::TextureFactory< double , 1 > > Factory;
		typedef typename Factory::VertexType Vertex;

		Eigen::VectorXd eValues = _eMatrix * _fineGeodesicX;
		std::vector< Vertex > vertices( _verticesAndNormals.size()/2 );
		for( unsigned int i=0 ; i<vertices.size() ; i++ )
		{
			vertices[i].get<0>() = _verticesAndNormals[i];
			vertices[i].get<1>() = _verticesAndNormals[ vertices.size() + i ];
			vertices[i].get<2>() = eValues[i];
		}

		Factory factory;
		PLY::WriteTriangles< Factory >( fileName , factory , vertices , triangles , PLY_BINARY_NATIVE );
//		PLY::WriteTriangles< Factory >( fileName , factory , vertices , triangles , PLY_ASCII );
	}
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::_drawMesh( bool useColor )
{
	if( !_vbo ) _setVBO();

	if( !_valueTextureHandle )
	{
		float texture[TEXTURE_WIDTH*3];
		for( unsigned int i=0 ; i<TEXTURE_WIDTH ; i++ ) texture[3*i+0] = texture[3*i+1] = texture[3*i+2] = (float)i/(TEXTURE_WIDTH-1);
		glGenTextures( 1, &_valueTextureHandle );
		glBindTexture( GL_TEXTURE_1D , _valueTextureHandle );
		glTexEnvf( GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE ,GL_MODULATE );
		glTexParameterf( GL_TEXTURE_1D , GL_TEXTURE_MIN_FILTER , GL_LINEAR );
		glTexParameterf( GL_TEXTURE_1D , GL_TEXTURE_WRAP_S , GL_CLAMP_TO_EDGE );
		glTexImage1D( GL_TEXTURE_1D , 0 , GL_RGB , TEXTURE_WIDTH , 0 , GL_RGB , GL_FLOAT , texture );
	}
	if( !_stripeTextureHandle )
	{
		float texture[TEXTURE_WIDTH*3];
		for( unsigned int i=0 ; i<TEXTURE_WIDTH ; i++ )
		{
			if( (i>>stripeShift)&1 ) texture[3*i+0] = 1.f , texture[3*i+1] = texture[3*i+2] = 0.f;
			else                     texture[3*i+2] = 1.f , texture[3*i+1] = texture[3*i+0] = 0.f;
		}
		glGenTextures( 1, &_stripeTextureHandle );
		glBindTexture( GL_TEXTURE_1D , _stripeTextureHandle );
		glTexEnvf( GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE ,GL_MODULATE );
		glTexParameterf( GL_TEXTURE_1D , GL_TEXTURE_MIN_FILTER , GL_LINEAR );
		glTexParameterf( GL_TEXTURE_1D , GL_TEXTURE_WRAP_S , GL_CLAMP_TO_EDGE );
		glTexImage1D( GL_TEXTURE_1D , 0 , GL_RGB , TEXTURE_WIDTH , 0 , GL_RGB , GL_FLOAT , texture );
	}

	if( _sourceVertex!=-1 ) glEnable( GL_TEXTURE_1D ) , glBindTexture( GL_TEXTURE_1D , useStripeTexture ? _stripeTextureHandle : _valueTextureHandle );
	else glDisable( GL_TEXTURE_1D );

	glBindBuffer( GL_ARRAY_BUFFER , _vbo[VBO_VERTEX] );
	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_NORMAL_ARRAY );
	glVertexPointer  ( 3 , GL_DOUBLE , 0 , (GLubyte*)NULL );
	glNormalPointer  (     GL_DOUBLE , 0 , (GLubyte*)NULL + sizeof( Point3D< double > ) * _verticesAndNormals.size()/2 );

	glBindBuffer( GL_ARRAY_BUFFER , showGeodesic ? _vbo[VBO_GEODESIC_TEXTURE] : _vbo[VBO_DIFFUSED_DELTA_TEXTURE] );
	glEnableClientState( GL_TEXTURE_COORD_ARRAY );
	glTexCoordPointer( 1 , GL_DOUBLE , 0 , (GLubyte*)NULL );

	glColor3d( 0.75 , 0.75 , 0.75 );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _vbo[VBO_TRIANGLE] );
	glDrawElements( GL_TRIANGLES , (GLsizei)(_triangles.size()*3) , GL_UNSIGNED_INT , NULL );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );

	glDisable( GL_TEXTURE_1D );
}

template< unsigned int Degree , bool Hierarchical >
void GeodesicsInHeatVisualization< Degree , Hierarchical >::draw3D( void )
{
	if( lightOn ) glEnable( GL_LIGHTING );
	else          glDisable( GL_LIGHTING );
	glColor3f( 0.75f , 0.75f , 0.75f );
	_drawMesh( showColor );

	if( edgeMode==EDGES_POLYGON || edgeMode==EDGES_TRIANGLES )
	{
		Point3D< double > f = camera().forward;
		f /= Point3D< double >::Length( f );
		f *= _scale / 1000.;
		glPushMatrix();
		glTranslated( -f[0] , -f[1] , -f[2] );
		glDisable( GL_LIGHTING );
		if( edgeMode==EDGES_POLYGON )
		{
			glColor3d( 0. , 0. , 0. );
			glLineWidth( 1.5f );
		}
		else if( edgeMode==EDGES_TRIANGLES )
		{
			glColor3d( 0.25 , 0.25 , 0.25 );
			glLineWidth( 0.75f );
		}
		glDepthMask( GL_FALSE );
		glEnable( GL_BLEND );
		glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
		glEnable( GL_LINE_SMOOTH );

		if( edgeMode==EDGES_POLYGON )
		{
			for( unsigned int p=0 ; p<_polygons.size() ; p++ )
			{
				glBegin( GL_LINE_LOOP );
				for( unsigned int i=0 ; i<_polygons[p].size() ; i++ ) glVertex3dv( &_fullVertices[ _polygons[p][i] ][0] );
				glEnd();
			}
		}
		else if( edgeMode==EDGES_TRIANGLES )
		{
			const SimplexMesh< Dim , Degree > &simplexMesh = _simplexRefinableCellMesh.simplexMesh();
			for( unsigned int s=0 ; s<simplexMesh.simplices() ; s++ )
			{
				const SimplexIndex< 2 , unsigned int > &tri = simplexMesh.simplex(s);
				glBegin( GL_LINE_LOOP );
				for( unsigned int i=0 ; i<3 ; i++ ) glVertex3dv( &_fullVertices[ tri[i] ][0] );
				glEnd();
			}
		}

		glDepthMask( GL_TRUE );
		glDisable( GL_LINE_SMOOTH );
		glDisable( GL_BLEND );
		glPopMatrix();
	}
	else if( edgeMode==EDGES_TRIANGLES )
	{
	}
}
