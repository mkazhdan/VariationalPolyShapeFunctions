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
#include "Misha/Visualization.h"
#include "Misha/Camera.h"
#include "Misha/Meshes.h"
#include "Misha/DeformableSolid.h"
#include "Misha/Visualization3D.h"


template< unsigned int Degree , bool Hierarchical >
struct DeformablePolyhedralMeshVisualization : public Misha::Viewable3D< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >
{
	typedef typename Misha::Viewable< DeformablePolyhedralMeshVisualization >::KeyboardCallBack KeyboardCallBack;
	typedef typename KeyboardCallBack::Modifiers KeyboardCallBackModifiers;
	using Misha::Viewable< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::callBacks;
	using Misha::Viewable< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::info;
	using Misha::Viewable< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::screenWidth;
	using Misha::Viewable< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::screenHeight;
	using Misha::Viewable< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::promptCallBack;
	using Misha::Viewable3D< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::cameraType;
	using Misha::Viewable3D< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::PERSPECTIVE_CAMERA;
	using Misha::Viewable3D< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::perspectiveCameraInfo;
	using Misha::Viewable3D< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::ORTHOGRAPHIC_CAMERA;
	using Misha::Viewable3D< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::orthographicCameraInfo;

	static const unsigned int Dim = 3;

	enum
	{
		DRAW_FACES ,
		DRAW_EDGES ,
	};

	bool showInitial , showEdges , showFaces;
	bool animate;

	unsigned int vCycles , gsIters;

	DeformablePolyhedralMeshVisualization( void );
	void init( const std::vector< Point< double , Dim > > &restVertices , const std::vector< std::vector< unsigned int > > &polygons , const std::vector< std::vector< std::pair< unsigned int , bool > > > &polyhedra , const std::vector< bool > &lockedVertices , SquareMatrix< double , Dim > xForm , unsigned int width , unsigned int height , double gravity , unsigned int subResolution , unsigned int coarseNodeDim , bool linearPrecision );

	void draw3D( void );
	void idle( void );

	static void Reset( DeformablePolyhedralMeshVisualization* v , const char * ){ v->reset() ; glutPostRedisplay(); }
	static void Advance( DeformablePolyhedralMeshVisualization* v , const char * ){ v->advance(); }
	static void ToggleAnimate( DeformablePolyhedralMeshVisualization* v , const char * ){ v->animate=!v->animate; }
	static void ToggleInitial( DeformablePolyhedralMeshVisualization* v , const char * ){ v->showInitial=!v->showInitial; }
	static void ToggleEdges( DeformablePolyhedralMeshVisualization* v , const char * ){ v->showEdges=!v->showEdges; }
	static void ToggleFaces( DeformablePolyhedralMeshVisualization* v , const char * ){ v->showFaces=!v->showFaces; }
	static void DecreaseGravity( DeformablePolyhedralMeshVisualization* v , const char * ){ v->_deformableSolid->gravity /= 2.; }
	static void IncreaseGravity( DeformablePolyhedralMeshVisualization* v , const char * ){ v->_deformableSolid->gravity *= 2.; }
	static void SetTimeStep( DeformablePolyhedralMeshVisualization* v , const char *prompt );
	static void SetYoungsModulus( DeformablePolyhedralMeshVisualization* v , const char *prompt );
	static void SetPoissonRatio( DeformablePolyhedralMeshVisualization* v , const char *prompt );

	void advance( void );
	void reset( void );
	void setTimeStep( double timeStep ){ _deformableSolid->setTimeStep( timeStep ); }

protected:
	typedef typename std::conditional< Hierarchical , HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree > , SolidSimplexRefinableCellMesh< Dim , Degree > >::type SolidSimplexRefinableCellMeshType;
	typedef DeformableSolid< Dim , Degree , SolidSimplexRefinableCellMeshType , Hierarchical > DeformableSolidType;

	SquareMatrix< double , Dim > _xForm;
	typedef SimplexIndex< 2 , unsigned int > Triangle;

	GLuint _vbo , _ebo;
	std::vector< Triangle > _tris;
	Meshes::PolyhedronMesh< unsigned int > _polyMesh;
	std::vector< SimplexIndex< 1 , unsigned int > > _edges;
	std::vector< Point< double , Dim > > _restVertices;
	std::function< Point< double , Dim > ( unsigned int ) > _vertexPositionFunction , _fullVertexPositionFunction;

	typename DeformableSolidType::AdvanceStats _aStats;

	int _youngsModulusLine , _poissonRatioLine , _timeStepLine , _timeLine , _meshSizeLine , _selectedVertex , _energyLine , _gravityLine , _advanceLine;
	double _time;
	SolidSimplexRefinableCellMeshType _solidSimplexRefinableCellMesh;
	DeformableSolidType *_deformableSolid;
	std::vector< bool > _lockedVertices;
	Eigen::VectorXd _x;
	Eigen::SparseMatrix< double , Eigen::RowMajor > _triE , _edgeE;

	void _drawMesh( bool rest , int drawType , Point< double , 3 > outsideColor , Point< double , 3 > insideColor , double lineWidth ) const;
};

template< unsigned int Degree , bool Hierarchical >
void DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::SetTimeStep( DeformablePolyhedralMeshVisualization *v , const char *prompt )
{
	if( prompt ) v->_deformableSolid->setTimeStep( atof(prompt) );
}

template< unsigned int Degree , bool Hierarchical >
void DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::SetYoungsModulus( DeformablePolyhedralMeshVisualization* v , const char *prompt )
{
	if( prompt )
	{
		typename DeformableSolidType::MaterialProperties mp = v->_deformableSolid->materialProperties();
		mp.youngsModulus = atof( prompt );
		v->_deformableSolid->setMaterialProperties( mp );
	}
}

template< unsigned int Degree , bool Hierarchical >
void DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::SetPoissonRatio( DeformablePolyhedralMeshVisualization* v , const char *prompt )
{
	if( prompt )
	{
		typename DeformableSolidType::MaterialProperties mp = v->_deformableSolid->materialProperties();
		mp.poissonRatio = atof( prompt );
		v->_deformableSolid->setMaterialProperties( mp );
	}
}

template< unsigned int Degree , bool Hierarchical >
DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::DeformablePolyhedralMeshVisualization( void ) : animate(false) , _time(0) , _selectedVertex(-1) , showInitial(true) , showEdges(true) , showFaces(true) , _vbo(0) , _ebo(0)
{
	callBacks.push_back( KeyboardCallBack( this , '[' , KeyboardCallBackModifiers() , "decrease gravity" , DecreaseGravity ) );
	callBacks.push_back( KeyboardCallBack( this , ']' , KeyboardCallBackModifiers() , "increase gravity" , IncreaseGravity ) );
	callBacks.push_back( KeyboardCallBack( this , ' ' , KeyboardCallBackModifiers() , "animate" , ToggleAnimate ) );
	callBacks.push_back( KeyboardCallBack( this , 'r' , KeyboardCallBackModifiers() , "reset" , Reset ) );
	callBacks.push_back( KeyboardCallBack( this , '+' , KeyboardCallBackModifiers() , "advance" , Advance ) );
	callBacks.push_back( KeyboardCallBack( this , 'y' , KeyboardCallBackModifiers() , "set young\'s modulus" , "Young\'s modulus" , SetYoungsModulus ) );
	callBacks.push_back( KeyboardCallBack( this , 'p' , KeyboardCallBackModifiers() , "set poisson ratio" , "Poisson ratio" , SetPoissonRatio ) );
	callBacks.push_back( KeyboardCallBack( this , 't' , KeyboardCallBackModifiers() , "set time step" , "Time-step" , SetTimeStep ) );
	callBacks.push_back( KeyboardCallBack( this , 'b' , KeyboardCallBackModifiers() , "toggle initial" , ToggleInitial ) );
	callBacks.push_back( KeyboardCallBack( this , 'e' , KeyboardCallBackModifiers() , "toggle edges" , ToggleEdges ) );
	callBacks.push_back( KeyboardCallBack( this , 'f' , KeyboardCallBackModifiers() , "toggle faces" , ToggleFaces ) );

	_youngsModulusLine = (int)info.size() , info.push_back( new char[512] );
	_poissonRatioLine  = (int)info.size() , info.push_back( new char[512] );
	_timeStepLine      = (int)info.size() , info.push_back( new char[512] );
	_gravityLine       = (int)info.size() , info.push_back( new char[512] );
	_energyLine        = (int)info.size() , info.push_back( new char[512] );
	_timeLine          = (int)info.size() , info.push_back( new char[512] );
	_advanceLine       = (int)info.size() , info.push_back( new char[512] );
	_meshSizeLine      = (int)info.size() , info.push_back( new char[512] );

	sprintf( info[_advanceLine] , "Update / Solve time:" );
}

template< unsigned int Degree , bool Hierarchical >
void DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::reset( void )
{
	_deformableSolid->setGeometricState( [&]( unsigned int idx ){ return _xForm * _fullVertexPositionFunction( idx ); } , []( unsigned int ){ return Point< double , Dim >(); } );
	_time = 0;
	animate = false;
}

template< unsigned int Degree , bool Hierarchical >
void DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::init( const std::vector< Point< double , Dim > > &restVertices , const std::vector< std::vector< unsigned int > > &polygons , const std::vector< std::vector< std::pair< unsigned int , bool > > > &polyhedra , const std::vector< bool > &lockedVertices , SquareMatrix< double , Dim > xForm , unsigned int width , unsigned int height , double gravity , unsigned int subResolution , unsigned int coarseNodeDim , bool linearPrecision )
{
	typedef MultiIndex< Degree , unsigned int , false > NodeMultiIndex;
	Point3D< double > bBox[2];
	bBox[0] = bBox[1] = restVertices[0];
	for( unsigned int i=0 ; i<restVertices.size() ; i++ )
	{
		Point< double , Dim > v0 = restVertices[i];
		Point< double , Dim > v1 = xForm * restVertices[i];
		for( unsigned int d=0 ; d<Dim ; d++ )
		{
			bBox[0][d] = std::min< double >( bBox[0][d] , std::min< double >( v0[d] , v1[d] ) );
			bBox[1][d] = std::max< double >( bBox[1][d] , std::max< double >( v0[d] , v1[d] ) );
		}
	}

	Misha::Viewable3D< DeformablePolyhedralMeshVisualization< Degree , Hierarchical > >::init( bBox , width , height );

	_restVertices = restVertices;
	_xForm = xForm;
	_lockedVertices = lockedVertices;

	Point< double , Dim > g;
	g[1] = gravity;

	SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
	eWeights.kWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1;

	_polyMesh = Meshes::PolyhedronMesh< unsigned int >( polyhedra , polygons );
	_vertexPositionFunction = [&]( unsigned int idx )
	{
		if( idx<_restVertices.size() ) return _restVertices[idx];
		ERROR_OUT( "Bad vertex index: " , idx , " / " ,  _restVertices.size() );
		return Point< double , Dim >();
	};
	_fullVertexPositionFunction = _polyMesh.fullVertexPositionFunction( _vertexPositionFunction , true );

	if constexpr( Hierarchical )
	{
		_solidSimplexRefinableCellMesh = _polyMesh.template hierarchicalSolidSimplexRefinableCellMesh< Degree >( _fullVertexPositionFunction , eWeights , coarseNodeDim , false , linearPrecision , 0 , true );
	}
	else
	{
		_solidSimplexRefinableCellMesh = _polyMesh.template solidSimplexRefinableCellMesh< Degree >( _fullVertexPositionFunction , eWeights , coarseNodeDim , false , linearPrecision , 0 , true );
	}
	_deformableSolid = new DeformableSolidType( _solidSimplexRefinableCellMesh , (unsigned int)restVertices.size() , _fullVertexPositionFunction , g , true );
	_deformableSolid->gsIters = gsIters;
	_deformableSolid->vCycles = vCycles;
	_deformableSolid->verbose = false;

	_deformableSolid->setTimeStep( 1e-6 );
	_x = _deformableSolid->x();

	reset();
	_deformableSolid->lock( lockedVertices );

	{
		std::vector< typename SimplexMesh< Dim >::Sample > edgeSamples , triangleSamples;

		// assuming that x <= subResolution
		auto EdgeVertexIndex = [&]( unsigned int x ){ return x + (unsigned int)edgeSamples.size(); };

		// assuming that x+y <= subResolution
		auto TriangleVertexIndex = [&]( unsigned int x , unsigned int y )
		{
			// Idx = x + Offset(y)
			// Offset(y) = Sum_{0<=i<y}( subResolution + 1 - i )
			//           = Sum_{0<=i<y}( subResolution + 1 ) - Sum_{0<=i<y}( i )
			//           = ( subResolution + 1 ) * y - y*(y-1)/2
			return x + ( subResolution + 1 ) * y - (y*(y-1))/2 + (unsigned int)triangleSamples.size();
		};

		std::map< unsigned int , unsigned int > polygonCount;
		for( unsigned int i=0 ; i<polyhedra.size() ; i++ ) for( unsigned int j=0 ; j<polyhedra[i].size() ; j++ ) polygonCount[ polyhedra[i][j].first ]++;

		unsigned int sOffset = 0;
		// Iterate over the polyhedra
		for( unsigned int i=0 ; i<_polyMesh.polyhedra() ; i++ )
		{
			Meshes::SimplexRefinablePolyhedron srp = _polyMesh.template simplexRefinable< Dim >( i , []( unsigned in ){ return Point< double , Dim >(); } );

			unsigned int _sOffset = 0;
			// For each face of the polyhedra
			for( unsigned int f=0 ; f<srp.faces() ; f++ )
			{
				const SimplexRefinableCell< Dim-1 > &face = srp.face( f );

				// If it's a boundary face
				if( polygonCount[ polyhedra[i][f].first ]==1 )
				{
					// For each triangle comprising the boundary face
					for( unsigned int s=0 ; s<face.size() ; s++ )
					{
						{
							typename SimplexMesh< 3 >::Sample sample;
							sample.sIdx = sOffset + _sOffset + s;

							SimplexIndex< Dim , unsigned int > tetIndex = srp[s+_sOffset];
							{
								unsigned int cutOff = (unsigned int)restVertices.size();
								unsigned int i1=-1 , i2=-1;
								for( unsigned int _i1=0 ; _i1<=Dim ; _i1++ ) for( unsigned int _i2=_i1+1 ; _i2<=Dim ; _i2++ )
									if( tetIndex[_i1]<cutOff && tetIndex[_i2]<cutOff ) i1 = _i1 , i2 = _i2;
								if( i1==-1 || i2==-1 ) ERROR_OUT( "Could not find edge" );
								for( unsigned int j=0 ; j<subResolution ; j++ )
									_edges.push_back( SimplexIndex< 1 , unsigned int >( EdgeVertexIndex(j+0) , EdgeVertexIndex(j+1) ) );
								for( unsigned int j=0 ; j<=subResolution ; j++ )
								{
									double x = (double)j/subResolution;
									sample.bcCoordinates = Point< double , Dim+1 >();
									sample.bcCoordinates[i1] = 1.-x;
									sample.bcCoordinates[i2] = x;
									edgeSamples.push_back( sample );
								}
							}
							{
								unsigned int cutOff = (unsigned int)( restVertices.size() + polygons.size() );
								unsigned int i1=-1 , i2=-1 , i3=-1;
								for( unsigned int _i1=0 ; _i1<=Dim ; _i1++ ) for( unsigned int _i2=_i1+1 ; _i2<=Dim ; _i2++ ) for( unsigned int _i3=_i2+1 ; _i3<=Dim ; _i3++ )
									if( tetIndex[_i1]<cutOff && tetIndex[_i2]<cutOff && tetIndex[_i3]<cutOff ) i1 = _i1 , i2 = _i2 , i3 = _i3;
								if( i1==-1 || i2==-1 || i3==-1 ) ERROR_OUT( "Could not find triangle" );
								for( unsigned int k=0 ; k<subResolution ; k++ ) for( unsigned int j=0 ; j+k<subResolution ; j++ )
								{
									_tris.push_back( SimplexIndex< 2 , unsigned int >( TriangleVertexIndex(j+0,k+1) , TriangleVertexIndex(j+1,k+0) , TriangleVertexIndex(j+0,k+0) ) );
									if( j+k<subResolution-1 )
										_tris.push_back( SimplexIndex< 2 , unsigned int >( TriangleVertexIndex(j+1,k+0) , TriangleVertexIndex(j+0,k+1) , TriangleVertexIndex(j+1,k+1) ) );
								}
								for( unsigned int k=0 ; k<=subResolution ; k++ ) for( unsigned int j=0 ; j+k<=subResolution ; j++ )
								{
									double x = (double)j/subResolution;
									double y = (double)k/subResolution;
									sample.bcCoordinates = Point< double , Dim+1 >();
									sample.bcCoordinates[i1] = 1.-x-y;
									sample.bcCoordinates[i2] = x;
									sample.bcCoordinates[i3] = y;
									triangleSamples.push_back( sample );
								}
							}
						}
					}
				}
				_sOffset += face.size();
			}
			sOffset += _sOffset;
		}

		if constexpr( Hierarchical )
		{
			Eigen::SparseMatrix< double > P = _solidSimplexRefinableCellMesh.P( _solidSimplexRefinableCellMesh.maxLevel() , coarseNodeDim );
			_edgeE = _solidSimplexRefinableCellMesh.solidSimplexMesh().evaluationMatrix( edgeSamples ) * P;
			_triE = _solidSimplexRefinableCellMesh.solidSimplexMesh().evaluationMatrix( triangleSamples ) * P;
		}
		else
		{
			_edgeE = _solidSimplexRefinableCellMesh.solidSimplexMesh().evaluationMatrix( edgeSamples ) * _solidSimplexRefinableCellMesh.P();
			_triE = _solidSimplexRefinableCellMesh.solidSimplexMesh().evaluationMatrix( triangleSamples ) * _solidSimplexRefinableCellMesh.P();
		}
	}
}

template< unsigned int Degree , bool Hierarchical >
void DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::_drawMesh( bool rest , int drawType , Point< double , 3 > outsideColor , Point< double , 3 > insideColor , double lineWidth ) const
{
	Camera camera;
	if     ( cameraType== PERSPECTIVE_CAMERA ) camera =  perspectiveCameraInfo.camera;
	else if( cameraType==ORTHOGRAPHIC_CAMERA ) camera = orthographicCameraInfo.camera;

	Eigen::VectorXd x;

	if( drawType==DRAW_FACES )
	{
		x = _triE * ( rest ? _x : _deformableSolid->x() );

		static std::vector< double > _vertexData( 2 * x.size() );
		static std::vector< Point< double  , 3 > > normals( _tris.size() );
		memset( &normals[0] , 0 , sizeof(double)*3*_tris.size() );
		Point3D< double > *_vertices  = ( Point< double , 3 >* )( &_vertexData[0] );
		Point3D< double > *_normals   = ( Point< double , 3 >* )( &_vertexData[0] + x.size() );

#pragma omp parallel for
		for( int i=0 ; i<(int)x.size()/3 ; i++ ) _vertices[i] = Point< double , 3 >( x[i*3+0] , x[i*3+1] , x[i*3+2] );
#pragma omp parallel for
		for( int t=0 ; t<(int)_tris.size() ; t++ ) normals[t] = Point< double , 3 >::CrossProduct( _vertices[ _tris[t][1] ] - _vertices[ _tris[t][0] ] , _vertices[ _tris[t][2] ] - _vertices[ _tris[t][0] ] );

		for( unsigned int t=0 ; t<_tris.size() ; t++ ) for( unsigned int i=0 ; i<3 ; i++ ) _normals[ _tris[t][i] ] += normals[t];

		glBindBuffer( GL_ARRAY_BUFFER , _vbo );
		glBufferData( GL_ARRAY_BUFFER , 2 * x.size() * sizeof( double ) , &_vertexData[0] , GL_DYNAMIC_DRAW );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );

		glEnable( GL_NORMALIZE  );
		glDisable( GL_COLOR_MATERIAL );
		GLfloat frontAmbientAndDiffuse[] = { (float)outsideColor[0] , (float)outsideColor[1] , (float)outsideColor[2] , 1. };
		GLfloat backAmbientAndDiffuse[] = { (float)insideColor[0] , (float)insideColor[1] , (float)insideColor[2] , 1. };
		glMaterialfv( GL_FRONT , GL_AMBIENT_AND_DIFFUSE , frontAmbientAndDiffuse );
		glMaterialfv( GL_BACK , GL_AMBIENT_AND_DIFFUSE , backAmbientAndDiffuse );
		
		glBindBuffer( GL_ARRAY_BUFFER , _vbo );
		glEnableClientState( GL_VERTEX_ARRAY );
		glEnableClientState( GL_NORMAL_ARRAY );
		glVertexPointer  ( 3 , GL_DOUBLE , 0 , (GLubyte*)NULL );
		glNormalPointer  (     GL_DOUBLE , 0 , (GLubyte*)NULL + sizeof( double ) * x.size() );

		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _ebo );
		glDrawElements( GL_TRIANGLES , (GLsizei)(_tris.size()*3) , GL_UNSIGNED_INT , NULL );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );

		glDisableClientState( GL_NORMAL_ARRAY );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );
	}
	else if( drawType==DRAW_EDGES )
	{
		auto DrawVertex = []( Point< double , Dim > p ){ glVertex3d( p[0] , p[1] , p[2] ); };
		auto Vertex = [&]( unsigned int v ){ return Point< double , Dim >( x[3*v+0] , x[3*v+1] , x[3*v+2] ); };

		x = _edgeE * ( rest ? _x : _deformableSolid->x() );
		bool lighting = glIsEnabled( GL_LIGHTING );
		if( lighting ) glDisable( GL_LIGHTING );

		glEnable( GL_COLOR_MATERIAL );
		glColor3d( outsideColor[0] , outsideColor[1] , outsideColor[2] );
		glLineWidth( (float)lineWidth );

		Point< double , Dim > offset = -camera.forward / Point< double , Dim >::Length( camera.forward ) / 100;
		glBegin( GL_LINES );
		for( unsigned int e=0 ; e<_edges.size() ; e++ )
		{
			DrawVertex( Vertex( _edges[e][0] ) + offset );
			DrawVertex( Vertex( _edges[e][1] ) + offset );
		}
		glEnd();
		if( lighting ) glEnable( GL_LIGHTING );
	}
}

template< unsigned int Degree , bool Hierarchical >
void DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::draw3D( void )
{
	if( !_ebo )
	{
		glGenBuffers( 1 , &_ebo );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _ebo );
		glBufferData( GL_ELEMENT_ARRAY_BUFFER , _tris.size() * sizeof( unsigned int ) * 3 , &_tris[0] , GL_STATIC_DRAW );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );
	}
	if( !_vbo )
	{
		glGenBuffers( 1 , &_vbo );
		glBindBuffer( GL_ARRAY_BUFFER , _vbo );
		glBufferData( GL_ARRAY_BUFFER , 2 * _triE.rows() * sizeof( double ) , NULL , GL_DYNAMIC_DRAW );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );
	}

	Point< double , 3 > black(0,0,0) , gray(0.5,0.5,0.5) , red(1.0,0.75,0.75) , blue(0.75,0.75,1.0);
	// Draw the rest state
	if( showInitial ) _drawMesh( true , DRAW_EDGES , gray , gray , 2. );

	// Draw the faces of deformed state
	if( showFaces ) _drawMesh( false , DRAW_FACES , red , blue , 0. );

	// Draw the wire-frame of deformed state
	if( showEdges ) _drawMesh( false , DRAW_EDGES , black , black , 2. );

	// Draw the locked vertices
	glDisable( GL_LIGHTING );
	glColor3f( 0. , 0. , 0. );
	glPointSize( 8. );
	glBegin( GL_POINTS );
	for( unsigned int v=0 ; v<_deformableSolid->vertices() ; v++ ) if( _deformableSolid->lockedVertex( v ) )
	{
		Point< double , Dim > vertex = _deformableSolid->vertex( v );
		glVertex3d( vertex[0] , vertex[1] , vertex[2] );
	}
	glEnd();
	glEnable( GL_LIGHTING );

	// Print out the properties
	typename DeformableSolidType::MaterialProperties mp = _deformableSolid->materialProperties();
	sprintf( info[_youngsModulusLine] , "Young\'s modulus: %g" , mp.youngsModulus );
	sprintf( info[_poissonRatioLine] , "Poisson ratio: %g" , mp.poissonRatio );
	sprintf( info[_timeStepLine] , "Time-step: %g" , _deformableSolid->timeStep() );
	sprintf( info[_timeLine] , "Time: %g" , _time );
	sprintf( info[_meshSizeLine] , "Nodes: %d" , (int)_solidSimplexRefinableCellMesh.nodes() );
	sprintf( info[_energyLine] , "Energy: %g" , _deformableSolid->energy() );
	sprintf( info[_gravityLine] , "Gravity: %g" , _deformableSolid->gravity );
}

template< unsigned int Degree , bool Hierarchical >
void DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::idle( void )
{
	if( !promptCallBack )
	{ 
		if( animate ) advance();
	}
}

template< unsigned int Degree , bool Hierarchical >
void DeformablePolyhedralMeshVisualization< Degree , Hierarchical >::advance( void )
{
	_deformableSolid->advance( _aStats );
	sprintf( info[_advanceLine] , "Update / Solve time: %.2f / %.2f" , _aStats.updateTime , _aStats.solveTime );
	_time += _deformableSolid->timeStep();
	glutPostRedisplay();
}
