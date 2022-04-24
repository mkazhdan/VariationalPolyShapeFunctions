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
#include "Misha/DeformableSolid.h"

template< unsigned int Level >
using RefinedVertexIndex = MultiIndex< (1<<Level) , unsigned int , false >;

struct PolygonEdgePosition
{
	unsigned int p;
	unsigned int e;
	double x;
};

template< unsigned int Level >
void RefinePolygonMesh( const std::vector< std::vector< unsigned int > > &coarseMesh , std::vector< PolygonEdgePosition > &vertices , std::vector< std::vector< unsigned int > > &fineMesh )
{
	static const unsigned int Res = 1<<Level;

	std::vector< std::vector< std::pair< RefinedVertexIndex< Level > , PolygonEdgePosition > > > pMesh( coarseMesh.size() );
	fineMesh.resize( coarseMesh.size() );

	std::map< RefinedVertexIndex< Level > , unsigned int > vMap;
	PolygonEdgePosition pep;
	unsigned int idx[ Res ];
	// Iterate over the coarse polygons
	for( unsigned int p=0 ; p<coarseMesh.size() ; p++ )
	{
		pep.p = p;
		pMesh[p].resize( coarseMesh[p].size()*Res );
		fineMesh[p].resize( coarseMesh[p].size()*Res );

		// Iterate over each edge
		for( unsigned int e=0 ; e<coarseMesh[p].size() ; e++ )
		{
			pep.e = e;

			// The end-points of the edge
			unsigned int v1 = coarseMesh[p][e] , v2 = coarseMesh[p][(e+1)%coarseMesh[p].size()];

			// Iterate over each point on the edge
			for( unsigned int r=0 ; r<Res ; r++ )
			{
				pep.x = (double)r/Res;
				for( unsigned int k=0 ; k<Res-r ; k++ ) idx[k] = v1;
				for( unsigned int k=Res-r ; k<Res ; k++ ) idx[k] = v2;
				RefinedVertexIndex< Level > vIdx( idx );
				pMesh[p][e*Res+r] = std::make_pair( vIdx , pep );
				vMap[vIdx] = 0;
			}
		}
	}
	{
		unsigned int vCount = 0;
		for( auto &iter : vMap ) iter.second = vCount++;
	}

	vertices.resize( vMap.size() );
	for( unsigned int p=0 ; p<pMesh.size() ; p++ ) for( unsigned int v=0 ; v<pMesh[p].size() ; v++ )
	{
		fineMesh[p][v] = vMap[ pMesh[p][v].first ];
		vertices[ fineMesh[p][v] ] = pMesh[p][v].second;
	}
}


template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
struct DeformablePolygonMeshVisualization : public Misha::Viewable< DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical > >
{
	typedef typename Misha::Viewable< DeformablePolygonMeshVisualization >::KeyboardCallBack KeyboardCallBack;
	typedef typename KeyboardCallBack::Modifiers KeyboardCallBackModifiers;
	using Misha::Viewable< DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical > >::callBacks;
	using Misha::Viewable< DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical > >::info;
	using Misha::Viewable< DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical > >::screenWidth;
	using Misha::Viewable< DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical > >::screenHeight;
	using Misha::Viewable< DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical > >::promptCallBack;

	static const unsigned int Dim = 2;
	static const unsigned int EdgeSamples = 1 + (1<<RefinementLevels);

	enum
	{
		DRAW_FACES ,
		DRAW_EDGES ,
		DRAW_BOUNDARY
	};

	int oldX , oldY , newX , newY;
	bool animate;
	bool showBoundary , showEdges , showFaces;

	unsigned int vCycles , gsIters;

	DeformablePolygonMeshVisualization( void );
	void init( const std::vector< Point< double , Dim > > &restVertices , const std::vector< std::vector< unsigned int > > &polygons , const std::vector< bool > &lockedVertices , SquareMatrix< double , Dim > xForm , unsigned int width , unsigned int height , double gravity , unsigned int coarseNodeDim );

	void idle( void );
	void keyboardFunc( unsigned char key , int x , int y );
	void specialFunc( int key, int x, int y );
	void display( void );
	void mouseFunc( int button , int state , int x , int y );
	void motionFunc( int x , int y );
	void reshape( int w , int h );

	static void Reset( DeformablePolygonMeshVisualization* v , const char * ){ v->reset() ; glutPostRedisplay(); }
	static void Advance( DeformablePolygonMeshVisualization* v , const char * ){ v->advance(); }
	static void ToggleAnimate( DeformablePolygonMeshVisualization* v , const char * ){ v->animate=!v->animate; }
	static void ToggleBoundary( DeformablePolygonMeshVisualization* v , const char * ){ v->showBoundary=!v->showBoundary; }
	static void ToggleEdges( DeformablePolygonMeshVisualization* v , const char * ){ v->showEdges=!v->showEdges; }
	static void ToggleFaces( DeformablePolygonMeshVisualization* v , const char * ){ v->showFaces=!v->showFaces; }
	static void DecreaseGravity( DeformablePolygonMeshVisualization* v , const char * ){ v->_deformableSolid->gravity /= 2.; }
	static void IncreaseGravity( DeformablePolygonMeshVisualization* v , const char * ){ v->_deformableSolid->gravity *= 2.; }
	static void SetTimeStep( DeformablePolygonMeshVisualization* v , const char *prompt );
	static void SetYoungsModulus( DeformablePolygonMeshVisualization* v , const char *prompt );
	static void SetPoissonRatio( DeformablePolygonMeshVisualization* v , const char *prompt );

	void advance( void );
	void reset( void );
	void setTimeStep( double timeStep ){ _deformableSolid->setTimeStep( timeStep ); }

protected:
	typedef typename std::conditional< Hierarchical , HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree > , SolidSimplexRefinableCellMesh< Dim , Degree > >::type SolidSimplexRefinableCellMeshType;
	typedef DeformableSolid< Dim , Degree , SolidSimplexRefinableCellMeshType , Hierarchical > DeformableSolidType;

	// The world-to-screen transformaiton maps the mesh into the square [-1,1]x[-1,1]
	XForm< double , Dim+1 > _worldToCamera , _cameraToWorld;
	SquareMatrix< double , Dim > _xForm;

	// camera-to-screen transformation maps the square [-1,1]x[-1,1] to the square [0,screenWidth]x[0,screenHeight]
	XForm< double , 3 > _cameraToScreen , _screenToCamera;


	Meshes::PolygonMesh< unsigned int > _polygonMesh;
	std::vector< std::vector< unsigned int > > _polygons , _refinedPolygons;
	std::vector< std::pair< unsigned int , unsigned int > > _boundaryEdges;
	std::vector< Point< double , Dim > > _restVertices;
	std::function< Point< double , Dim > ( unsigned int ) > _vertexPositionFunction , _fullVertexPositionFunction;
	typedef MultiIndex< EdgeSamples , unsigned int , false > VertexIndex;
	std::vector< VertexIndex > _vertices;

	int _youngsModulusLine , _poissonRatioLine , _timeStepLine , _timeLine , _meshSizeLine , _selectedVertex , _energyLine , _gravityLine , _advanceLine;
	double _time;
	SolidSimplexRefinableCellMeshType _solidSimplexRefinableCellMesh;
	DeformableSolidType *_deformableSolid;
	std::vector< bool > _lockedVertices;
	Eigen::VectorXd _x;
	Eigen::SparseMatrix< double > _E;

	void _setCameraToScreenXForm();
	void _drawMesh( bool rest , int drawType ) const;
	void _drawDisk( Point< double , Dim > center , double radius , unsigned int res ) const;
};

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::SetTimeStep( DeformablePolygonMeshVisualization *v , const char *prompt )
{
	if( prompt ) v->_deformableSolid->setTimeStep( atof(prompt) );
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::SetYoungsModulus( DeformablePolygonMeshVisualization* v , const char *prompt )
{
	if( prompt )
	{
		typename DeformableSolidType::MaterialProperties mp = v->_deformableSolid->materialProperties();
		mp.youngsModulus = atof( prompt );
		v->_deformableSolid->setMaterialProperties( mp );
	}
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::SetPoissonRatio( DeformablePolygonMeshVisualization* v , const char *prompt )
{
	if( prompt )
	{
		typename DeformableSolidType::MaterialProperties mp = v->_deformableSolid->materialProperties();
		mp.poissonRatio = atof( prompt );
		v->_deformableSolid->setMaterialProperties( mp );
	}
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::DeformablePolygonMeshVisualization( void ) : animate(false) , _time(0) , _selectedVertex(-1) , showBoundary(true) , showEdges(true) , showFaces(true)
{
	callBacks.push_back( KeyboardCallBack( this , '[' , KeyboardCallBackModifiers() , "decrease gravity" , DecreaseGravity ) );
	callBacks.push_back( KeyboardCallBack( this , ']' , KeyboardCallBackModifiers() , "increase gravity" , IncreaseGravity ) );
	callBacks.push_back( KeyboardCallBack( this , ' ' , KeyboardCallBackModifiers() , "animate" , ToggleAnimate ) );
	callBacks.push_back( KeyboardCallBack( this , 'r' , KeyboardCallBackModifiers() , "reset" , Reset ) );
	callBacks.push_back( KeyboardCallBack( this , '+' , KeyboardCallBackModifiers() , "advance" , Advance ) );
	callBacks.push_back( KeyboardCallBack( this , 'y' , KeyboardCallBackModifiers() , "set young\'s modulus" , "Young\'s modulus" , SetYoungsModulus ) );
	callBacks.push_back( KeyboardCallBack( this , 'p' , KeyboardCallBackModifiers() , "set poisson ratio" , "Poisson ratio" , SetPoissonRatio ) );
	callBacks.push_back( KeyboardCallBack( this , 't' , KeyboardCallBackModifiers() , "set time step" , "Time-step" , SetTimeStep ) );
	callBacks.push_back( KeyboardCallBack( this , 'b' , KeyboardCallBackModifiers() , "toggle boundary" , ToggleBoundary ) );
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

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::reset( void )
{
	std::function< Point< double , Dim > ( unsigned int  ) > vFunction = [&]( unsigned int idx ){ return _xForm * _fullVertexPositionFunction( idx ); };
	_deformableSolid->setGeometricState( vFunction , []( unsigned int ){ return Point< double , Dim >(); } );
	_time = 0;
	animate = false;
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::init( const std::vector< Point< double , Dim > > &restVertices , const std::vector< std::vector< unsigned int > > &polygons , const std::vector< bool > &lockedVertices , SquareMatrix< double , Dim > xForm , unsigned int width , unsigned int height , double gravity , unsigned int coarseNodeDim )
{
	typedef MultiIndex< Degree , unsigned int , false > NodeMultiIndex;

	_restVertices = restVertices;
	_xForm = xForm;
	_lockedVertices = lockedVertices;
	Point< double , Dim > g;
	g[Dim-1] = gravity;
	_polygons = _refinedPolygons = polygons;

	_polygonMesh = Meshes::PolygonMesh< unsigned int >( _polygons );

	SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
	eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

	_vertexPositionFunction = [&]( unsigned int idx )
	{
		if( idx<_restVertices.size() ) return _restVertices[idx];
		ERROR_OUT( "Bad vertex index: " , idx , " / " ,  _restVertices.size() );
		return Point< double , Dim >();
	};
	_fullVertexPositionFunction = _polygonMesh.fullVertexPositionFunction( _vertexPositionFunction , true );

	if constexpr( Hierarchical )
	{
		_solidSimplexRefinableCellMesh = _polygonMesh.template hierarchicalSolidSimplexRefinableCellMesh< Degree >( _fullVertexPositionFunction , eWeights , coarseNodeDim , false , false , 0 , true );
	}
	else
	{
		_solidSimplexRefinableCellMesh = _polygonMesh.template solidSimplexRefinableCellMesh< Degree >( _fullVertexPositionFunction , eWeights , coarseNodeDim , false , false , 0 , true );
	}
	_deformableSolid = new DeformableSolidType( _solidSimplexRefinableCellMesh , (unsigned int)restVertices.size() , _fullVertexPositionFunction , g );
	_deformableSolid->gsIters = gsIters;
	_deformableSolid->vCycles = vCycles;
	_deformableSolid->verbose = false;

	_deformableSolid->setTimeStep( 1e-6 );
	_x = _deformableSolid->x();

	screenWidth = width , screenHeight = height;

	Point< double , Dim > center;
	double radius;
	{
		Point< double , Dim > min , max;
		for( unsigned int d=0 ; d<Dim ; d++ ) min[d] = std::numeric_limits< double >::infinity() , max[d] = -std::numeric_limits< double >::infinity();
		for( unsigned int i=0 ; i<restVertices.size() ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ )
			min[d] = std::min< double >( min[d] , restVertices[i][d] ) , max[d] = std::max< double >( max[d] , restVertices[i][d] );
		center = (min+max)/2;
		radius = 0;
		for( unsigned int i=0 ; i<restVertices.size() ; i++ ) radius = std::max< double >( radius , ( restVertices[i]-center ).squareNorm() );
		radius = sqrt( radius );
	}
	_cameraToWorld = XForm< double , Dim+1 >::Identity();
	for( unsigned int d=0 ; d<Dim ; d++ ) _cameraToWorld(d,d) = radius*1.25 , _cameraToWorld(Dim,d) = center[d];
	_worldToCamera = _cameraToWorld.inverse();
	_setCameraToScreenXForm();

	reset();
	_deformableSolid->lock( lockedVertices );

	{
		std::vector< unsigned int > simplexStartIndex( polygons.size() );
		simplexStartIndex[0] = 0;
		for( unsigned int p=1 ; p<polygons.size() ; p++ ) simplexStartIndex[p] = simplexStartIndex[p-1] + (unsigned int)polygons[p-1].size();

		std::vector< PolygonEdgePosition > vertices;
		RefinePolygonMesh< RefinementLevels >( polygons , vertices , _refinedPolygons );

		std::vector< typename SimplexMesh< Dim >::Sample > samples( vertices.size() );
		for( unsigned int i=0 ; i<vertices.size() ; i++ )
		{
			const Meshes::Polygon< unsigned int > &polygon = polygons[ vertices[i].p ];
			std::function< Point< double , 2 > ( unsigned int ) > vFunction = []( unsigned int ){ return Point< double , 2 >(); };
			Meshes::SimplexRefinablePolygon srp( polygon , (unsigned int)restVertices.size() + vertices[i].p , vFunction );
#pragma message( "[WARNING] Assuming that the simplex index is the index of the edge" )
			SimplexIndex< 2 , unsigned int > simplexIndex = srp[ vertices[i].e ];
			samples[i].sIdx = simplexStartIndex[ vertices[i].p ] + vertices[i].e;
			bool found = false;
			for( unsigned int j=0 ; j<3 ; j++ )
				if( simplexIndex[j]==polygon[ vertices[i].e ] && simplexIndex[(j+1)%3]==polygon[ (vertices[i].e+1)%polygon.size() ] )
				{
					samples[i].bcCoordinates[j] = 1.-vertices[i].x;
					samples[i].bcCoordinates[(j+1)%3] = vertices[i].x;
					found = true;
				}
			if( !found ) ERROR_OUT( "Could not find edge in simplex" );
		}
		{
			typedef MultiIndex< 2 , unsigned int , true > EdgeIndex;
			std::map< EdgeIndex , unsigned int > eMap;

			for( unsigned int p=0 ; p<_refinedPolygons.size() ; p++ ) for( unsigned int v=0 ; v<_refinedPolygons[p].size() ; v++ )
				eMap[ EdgeIndex( _refinedPolygons[p][v] , _refinedPolygons[p][(v+1)%_refinedPolygons[p].size()] ) ]++;

			for( auto &i : eMap ) if( i.second==1 ) _boundaryEdges.push_back( std::make_pair( i.first[0] , i.first[1] ) );
		}

		if constexpr( Hierarchical )
		{
			_E = _solidSimplexRefinableCellMesh.solidSimplexMesh().evaluationMatrix( samples ) * _solidSimplexRefinableCellMesh.P( _solidSimplexRefinableCellMesh.maxLevel() , coarseNodeDim );
		}
		else
		{
			_E = _solidSimplexRefinableCellMesh.solidSimplexMesh().evaluationMatrix( samples ) * _solidSimplexRefinableCellMesh.P();
		}
	}
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::reshape( int w , int h )
{
	Misha::Viewable< DeformablePolygonMeshVisualization >::reshape( w , h );
	_setCameraToScreenXForm();
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::_setCameraToScreenXForm( void )
{
	_cameraToScreen = XForm< double , 3 >::Identity();
	_cameraToScreen(0,0) = _cameraToScreen(1,1) = std::min< double >( screenWidth , screenHeight )/2;
	_cameraToScreen(2,0) = (double)screenWidth/2;
	_cameraToScreen(2,1) = (double)screenHeight/2;
	_screenToCamera = _cameraToScreen.inverse();
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::_drawMesh( bool rest , int drawType ) const
{
	Eigen::VectorXd x = _E * ( rest ? _x : _deformableSolid->x() );

	auto drawVertexIndex = [&]( unsigned int v )
	{
		Point< double , Dim > p = _cameraToScreen * ( _worldToCamera * Point< double , Dim >( x[2*v+0] , x[2*v+1] ) );
		glVertex2d( p[0] , p[1] );
	};
	auto drawVertex = [&]( Point< double , Dim > p )
	{
		p = _cameraToScreen * ( _worldToCamera * p );
		glVertex2d( p[0] , p[1] );
	};

	if( drawType==DRAW_EDGES )
	{
		for( unsigned int p=0 ; p<_refinedPolygons.size() ; p++ )
		{
			glBegin( GL_LINE_LOOP );
			for( unsigned int v=0 ; v<_refinedPolygons[p].size() ; v++ ) drawVertexIndex( _refinedPolygons[p][v] );
			glEnd();
		}
	}
	else if( drawType==DRAW_BOUNDARY )
	{
		glBegin( GL_LINES );
		for( unsigned int b=0 ; b<_boundaryEdges.size() ; b++ )
		{
			drawVertexIndex( _boundaryEdges[b].first );
			drawVertexIndex( _boundaryEdges[b].second );
		}
		glEnd();
	}
	else if( drawType==DRAW_FACES )
	{
		for( unsigned int p=0 ; p<_refinedPolygons.size() ; p++ )
		{
			Point2D< double > center;
			for( unsigned int v=0 ; v<_refinedPolygons[p].size() ; v++ ) center += Point2D< double >( x[2*_refinedPolygons[p][v]+0] , x[2*_refinedPolygons[p][v]+1] );
			center /= (double)_refinedPolygons[p].size();
			glBegin( GL_TRIANGLE_FAN );
			drawVertex( center );
			for( unsigned int v=0 ; v<_refinedPolygons[p].size() ; v++ ) drawVertexIndex( _refinedPolygons[p][v] );
			drawVertexIndex( _refinedPolygons[p][0] );
			glEnd();
		}
	}
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::_drawDisk( Point< double , Dim > center , double radius , unsigned int res ) const
{
	glBegin( GL_POLYGON );
	for( unsigned int i=0 ; i<res ; i++ )
	{
		double theta = 2.*M_PI*i/res;
		Point< double , Dim > p = center + Point< double , Dim >( cos(theta) , sin(theta) ) * radius;
		glVertex2f( (float)p[0] , (float)p[1] );
	}
	glEnd();
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::display( void )
{
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glOrtho( 0 , screenWidth , 0 , screenHeight , -1.5 , 1.5 );

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();

	// Draw the rest state
	glColor3f( 0.75f , 0.75f , 0.75f );
	if( showFaces ) _drawMesh( true , DRAW_FACES );
	glLineWidth( 1.f );
	glColor3f( 0.5f , 0.5f , 0.5f );
	if( showEdges ) _drawMesh( true , DRAW_EDGES );
	glLineWidth( 2.f );
	if( showBoundary ) _drawMesh( true , DRAW_BOUNDARY );

	// Draw the deformed state
	glColor3f( 0.5f , 0.5f , 0.5f );
	if( showFaces ) _drawMesh( false , DRAW_FACES );
	glLineWidth( 1.f );
	glColor3f( 0.f , 0.f , 0.f );
	if( showEdges ) _drawMesh( false , DRAW_EDGES );
	glLineWidth( 2.f );
	if( showBoundary ) _drawMesh( false , DRAW_BOUNDARY );

	// Draw the locked vertices
	glColor3f( 0.f , 0.f , 0.f );
	for( unsigned int v=0 ; v<_deformableSolid->vertices() ; v++ ) if( _deformableSolid->lockedVertex( v ) )
		_drawDisk( _cameraToScreen * ( _worldToCamera * _deformableSolid->vertex( v ) ) , 5 , 16 );

	// Draw the selected vertex
	glColor3f( 1.f , 0.f , 0.f );
	if( _selectedVertex!=-1 ) _drawDisk( _cameraToScreen * ( _worldToCamera * _deformableSolid->vertex( _selectedVertex ) ) , 5 , 16 );

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

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::mouseFunc( int button , int state , int x , int y )
{
	newX = x ; newY = y;

	static const unsigned int SelectionRadius = 50;

	if( button==GLUT_LEFT_BUTTON )
	{
		if( state==GLUT_DOWN )
		{
			if( _selectedVertex!=-1 ) _lockedVertices[ _selectedVertex ] = false;
			//			_deformableSolid->lock( lockedVertices );

			double d = std::numeric_limits< double >::infinity();
			unsigned int idx = -1;

			Point2D< double > q( x , screenHeight-y );
			for( unsigned int v=0 ; v<_deformableSolid->vertices() ; v++ )
			{
				Point2D< double > p = _cameraToScreen * ( _worldToCamera * _deformableSolid->vertex( v ) );
				double d2 = ( p - q ).squareNorm() ;
				if( d2<d ) d = d2 , idx = v;
			}
			_selectedVertex = d<(SelectionRadius*SelectionRadius) ? idx : -1;
			if( _selectedVertex!=-1 )
			{
				if( _lockedVertices[ _selectedVertex ] ) _selectedVertex = -1;
				else _lockedVertices[ _selectedVertex ] = true;
			}
			_deformableSolid->lock( _lockedVertices );
		}
		else
		{
			if( _selectedVertex!=-1 ) _lockedVertices[ _selectedVertex ] = false;
			_selectedVertex = -1;
			_deformableSolid->lock( _lockedVertices );
		}
	}
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::motionFunc( int x , int y )
{
	if( _selectedVertex!=-1 )
	{
		Point< double , Dim > q = _cameraToWorld * ( _screenToCamera * Point2D< double >( x , screenHeight-y ) );
		_deformableSolid->setVertex( _selectedVertex , q );
		glutPostRedisplay();
	}
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::advance( void )
{
	typename DeformableSolidType::AdvanceStats aStats = _deformableSolid->advance();
	sprintf( info[_advanceLine] , "Update / Solve time: %.2f / %.2f" , aStats.updateTime , aStats.solveTime );
	_time += _deformableSolid->timeStep();
	glutPostRedisplay();
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::idle( void )
{
	if( !promptCallBack )
	{ 
		if( animate ) advance();
	}
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::keyboardFunc( unsigned char key , int x , int y )
{
}

template< unsigned int Degree , unsigned int RefinementLevels , bool Hierarchical >
void DeformablePolygonMeshVisualization< Degree , RefinementLevels , Hierarchical >::specialFunc( int key, int x, int y )
{
}