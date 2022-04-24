/* -*- C++ -*-
Copyright (c) 2019, Michael Kazhdan
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

namespace PLY
{
#ifdef NEW_PLY
	template<> inline int Type< int           >( void ){ return PLY_INT   ; }
	template<> inline int Type<          char >( void ){ return PLY_CHAR  ; }
	template<> inline int Type< unsigned char >( void ){ return PLY_UCHAR ; }
	template<> inline int Type<        float  >( void ){ return PLY_FLOAT ; }
	template<> inline int Type<        double >( void ){ return PLY_DOUBLE; }
	template< class Real > inline int Type( void )
	{
		ERROR_OUT( "Unrecognized type" );
		return -1;
	}

	template<> const std::string Traits<          int >::name="int";
	template<> const std::string Traits< unsigned int >::name="unsigned int";
	template<> const std::string Traits<          long >::name="long";
	template<> const std::string Traits< unsigned long >::name="unsigned long";
	template<> const std::string Traits<          long long >::name="long long";
	template<> const std::string Traits< unsigned long long >::name="unsigned long long";

	template<>
	PlyProperty Face<          int       >::Properties[] = { PlyProperty( "vertex_indices" , PLY_INT       , PLY_INT       , offsetof( Face , vertices ) , 1 , PLY_INT , PLY_INT , offsetof( Face , nr_vertices ) ) };
	template<>
	PlyProperty Face< unsigned int       >::Properties[] = { PlyProperty( "vertex_indices" , PLY_UINT      , PLY_UINT      , offsetof( Face , vertices ) , 1 , PLY_INT , PLY_INT , offsetof( Face , nr_vertices ) ) };
	template<>
	PlyProperty Face<          long long >::Properties[] = { PlyProperty( "vertex_indices" , PLY_LONGLONG  , PLY_LONGLONG  , offsetof( Face , vertices ) , 1 , PLY_INT , PLY_INT , offsetof( Face , nr_vertices ) ) };
	template<>
	PlyProperty Face< unsigned long long >::Properties[] = { PlyProperty( "vertex_indices" , PLY_ULONGLONG , PLY_ULONGLONG , offsetof( Face , vertices ) , 1 , PLY_INT , PLY_INT , offsetof( Face , nr_vertices ) ) };

#else // !NEW_PLY
	struct Face
	{
		int nr_vertices;
		int *vertices;
	};
	const PlyProperty FaceProps[] =
	{
		PlyProperty( "vertex_indices" , PLY_INT , PLY_INT , offsetof( Face , vertices ) , 1 , PLY_INT, PLY_INT , offsetof( Face , nr_vertices ) ) ,
	};
#endif // NEW_PLY

	struct Edge{ int v1 , v2; };
	const PlyProperty EdgeProps[] =
	{ 
		{ "v1" , PLY_INT , PLY_INT , (int)offsetof( Edge , v1 ) , 0 , 0 , 0 , 0 },
		{ "v2" , PLY_INT , PLY_INT , (int)offsetof( Edge , v2 ) , 0 , 0 , 0 , 0 }
	};

	// Read
	inline void ReadHeader( std::string fileName , const PlyProperty *properties , int propertyNum , bool *readFlags , int &file_type )
	{
		std::vector< std::string > elist;
		float version;

		PlyFile *ply = PlyFile::Read( fileName , elist , file_type , version );
		if( !ply ) THROW( "could not create read ply file: " , fileName );

		for( int i=0 ; i<(int)elist.size() ; i++ ) if( elist[i]=="vertex" ) for( int j=0 ; j<propertyNum ; j++ ) if( readFlags ) readFlags[j] = ply->get_property( elist[i] , &properties[j] )!=0;

		delete ply;
	}

	inline void ReadHeader( std::string fileName , const PlyProperty *properties , int propertyNum , bool *readFlags )
	{
		int file_type;
		ReadHeader( fileName , properties , propertyNum , readFlags , file_type );
	}

#ifdef NEW_PLY
#ifdef USE_PLY_FACTORY
	template< class VertexFactory , typename Index >
	void Read
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		std::vector< typename VertexFactory::VertexType >& vertices , 
		std::vector< std::pair< Index , Index > >* edges ,
		std::vector< std::vector< Index > >* polygons ,
		bool *vertexPropertiesFlag ,
		int &file_type ,
		std::vector< std::string > *comments
	)
#else // !USE_PLY_FACTORY
	template< class Vertex , typename Index >
	void Read
	(
		std::string fileName ,
		std::vector< Vertex >& vertices , 
		std::vector< std::pair< Index , Index > >* edges ,
		std::vector< std::vector< Index > >* polygons ,
		const PlyProperty* vertexProperties ,
		bool* vertexPropertiesFlag ,
		int vertexPropertyNum ,
		int& file_type ,
		std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
#else // !NEW_PLY
	template< class Vertex >
	void Read
	(
		std::string fileName ,
		std::vector< Vertex >& vertices , 
		std::vector< std::pair< int , int > >* edges ,
		std::vector< std::vector< int > >* polygons ,
		const PlyProperty* vertexProperties ,
		bool* vertexPropertiesFlag ,
		int vertexPropertyNum ,
		int& file_type ,
		std::vector< std::string > *comments
	)
#endif // NEW_PLY
	{

		float version;
		std::vector< std::string > elist;

		PlyFile *ply = PlyFile::Read( fileName , elist , file_type , version );
		if( !ply ) THROW( "could not create read ply file: " , fileName );

		if( comments )
		{
			comments->reserve( comments->size() + ply->comments.size() );
			for( int i=0 ; i<(int)ply->comments.size() ; i++ ) comments->push_back( ply->comments[i] );
		}

		for( int i=0 ; i<elist.size() ; i++ )
		{
			std::string &elem_name = elist[i];
			size_t num_elems;
			std::vector< PlyProperty * > plist = ply->get_element_description( elem_name , num_elems );
			if( !num_elems ) continue;
			else if( !plist.size() )
			{
				delete ply;
				THROW( "could not get element description for: " , elem_name );
			}

			if( elem_name=="vertex" )
			{
#ifdef USE_PLY_FACTORY
				for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++ )
				{
					PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(i) : vFactory.plyReadProperty(i);
					int hasProperty = ply->get_property( elem_name , &prop );
					if( vertexPropertiesFlag ) vertexPropertiesFlag[i] = (hasProperty!=0);
				}
				vertices.resize( num_elems , vFactory() );
#else // !USE_PLY_FACTORY
				for( int i=0 ; i<vertexPropertyNum ; i++)
				{
					int hasProperty = ply->get_property( elem_name , &vertexProperties[i] );
					if( vertexPropertiesFlag ) vertexPropertiesFlag[i] = (hasProperty!=0);
				}
				vertices.resize( num_elems );
#endif // USE_PLY_FACTORY
				for( int j=0 ; j<num_elems ; j++ ) ply->get_element( (void *)&vertices[j] );
			}
			else if( elem_name=="face" && polygons )
			{
#ifdef NEW_PLY
				ply->get_property( elem_name , &Face< Index >::Properties[0] );
#else // !NEW_PLY
				ply->get_property( elem_name , &FaceProps[0] );
#endif // NEW_PLY
				polygons->resize( num_elems );
				for( int j=0 ; j<num_elems ; j++ )
				{
					Face< Index > ply_face;
					ply->get_element( (void *)&ply_face );
					(*polygons)[j].resize( ply_face.nr_vertices );
					for( unsigned int k=0 ; k<ply_face.nr_vertices ; k++ ) (*polygons)[j][k] = ply_face.vertices[k];
					free( ply_face.vertices );
				}  // for, read faces
			}  // if face
			else if( elem_name=="edge" && edges )
			{
				ply->get_property( elem_name , &EdgeProps[0] );
				ply->get_property( elem_name , &EdgeProps[1] );
				edges->resize( num_elems );
				for( int j=0 ; j<num_elems ; j++ )
				{
					Edge ply_edge;
					ply->get_element( (void*)&ply_edge );
					(*edges)[j].first = ply_edge.v1 , (*edges)[j].second = ply_edge.v2;
				}
			}
			else ply->get_other_element( elem_name , num_elems );

			for( int j=0 ; j<plist.size() ; j++ ) delete plist[j];
		}  // for each type of element
		delete ply;
	}

#ifdef USE_PLY_FACTORY
	template< class VertexFactory >
	void ReadVertices
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		std::vector< typename VertexFactory::VertexType > &vertices ,
		bool* vertexPropertiesFlag ,
		int& file_type ,
		std::vector< std::string > *comments
	)
	{
		return Read< VertexFactory , unsigned int >( fileName , vFactory , vertices , NULL , NULL , vertexPropertiesFlag , file_type , comments );
	}
#else // !USE_PLY_FACTORY
	template< class Vertex >
	void ReadVertices
	(
		std::string fileName ,
		std::vector< Vertex > &vertices ,
		const PlyProperty *vertexProperties ,
		bool* vertexPropertiesFlag ,
		int vertexPropertyNum ,
		int& file_type ,
		std::vector< std::string > *comments
	)
	{
		return Read( fileName , vertices , NULL , NULL , vertexProperties , vertexPropertiesFlag , vertexPropertyNum , file_type , comments );
	}
#endif // USE_PLY_FACTORY

#ifdef NEW_PLY
#ifdef USE_PLY_FACTORY
	template< typename VertexFactory , typename Real , unsigned int Dim , typename Index >
	void ReadTriangles
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		std::vector< typename VertexFactory::VertexType > &vertices ,
		std::vector< SimplexIndex< 2 , Index > > &triangles ,
		bool* vertexPropertiesFlag ,
		int& file_type ,
		std::function< Point< Real , Dim > ( typename VertexFactory::VertexType ) > VertexToPointFunctor ,
		std::vector< std::string > *comments
	)
#else // !USE_PLY_FACTORY
	template< class Vertex , typename Real , unsigned int Dim , typename Index >
	void ReadTriangles
	(
		std::string fileName ,
		std::vector< Vertex > &vertices ,
		std::vector< SimplexIndex< 2 , Index > > &triangles ,
		const PlyProperty *vertexProperties ,
		bool* vertexPropertiesFlag ,
		int vertexPropertyNum ,
		int& file_type ,
		std::function< Point< Real , Dim > (Vertex) > VertexToPointFunctor ,
		std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
#else // !NEW_PLY
	template< class Vertex >
	void ReadTriangles
	(
		std::string fileName ,
		std::vector< Vertex > &vertices ,
		std::vector< TriangleIndex > &triangles ,
		const PlyProperty *vertexProperties ,
		bool* vertexPropertiesFlag ,
		int vertexPropertyNum ,
		int& file_type ,
		std::function< Point3D< double > (Vertex) > VertexToPointFunctor ,
		std::vector< std::string > *comments
	)
#endif // NEW_PLY
	{
#ifdef NEW_PLY
		MinimalAreaTriangulation< Real , Dim > MAT;
		std::vector< std::vector< Index > > polygons;
#else // !NEW_PLY
//		MinimalAreaTriangulation< double > MAT;
		MinimalAreaTriangulation< double , 3 > MAT;
		std::vector< std::vector< int > > polygons;
#endif // NEW_PLY
#ifdef USE_PLY_FACTORY
		ReadPolygons( fileName , vFactory , vertices , polygons , vertexPropertiesFlag , file_type , comments );
#else // !USE_PLY_FACTORY
		ReadPolygons( fileName , vertices , polygons , vertexProperties , vertexPropertiesFlag , vertexPropertyNum , file_type , comments );
#endif // USE_PLY_FACTORY
#ifdef NEW_PLY
		std::vector< Point3D< Real > > poly;
		std::vector< SimplexIndex< 2 , Index > > tris;
#else // !NEW_PLY
		std::vector< Point3D< double > > poly;
		std::vector< TriangleIndex > tris;
#endif // NEW_PLY

		triangles.clear();
		for( unsigned int i=0 ; i<polygons.size() ; i++ )
		{
			poly.resize( polygons[i].size( ) );
			for( unsigned int j=0 ; j<polygons[i].size() ; j++ ) poly[j] = VertexToPointFunctor( vertices[ polygons[i][j] ] );
			MAT.GetTriangulation( poly , tris );
			for( unsigned int j=0 ; j<tris.size() ; j++ )
			{
#ifdef NEW_PLY
				SimplexIndex< 2 , Index > tri;
#else // !NEW_PLY
				TriangleIndex tri;
#endif // NEW_PLY
				tri[0] = polygons[i][ tris[j][0] ];
				tri[1] = polygons[i][ tris[j][1] ];
				tri[2] = polygons[i][ tris[j][2] ];
				triangles.push_back( tri );
			}
		}
	}


#ifdef USE_PLY_FACTORY
	template< typename VertexFactory , typename Index >
	void ReadTriangles
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		std::vector< typename VertexFactory::VertexType > &vertices ,
		std::vector< SimplexIndex< 2 , Index > > &triangles ,
		bool* vertexPropertiesFlag ,
		int& file_type ,
		std::vector< std::string > *comments
	)
	{
		std::vector< std::vector< Index > > polygons;
		ReadPolygons( fileName , vFactory , vertices , polygons , vertexPropertiesFlag , file_type , comments );
		triangles.resize( polygons.size() );
		for( unsigned int i=0 ; i<polygons.size() ; i++ )
		{
			if( polygons[i].size()!=3 ) ERROR_OUT( "Polygon is not a triangle: " , polygons[i].size() , " != " , 3 );
			for( int j=0 ; j<3 ; j++ ) triangles[i][j] = polygons[i][j];
		}
	}

#endif // USE_PLY_FACTORY
#ifdef NEW_PLY
#ifdef USE_PLY_FACTORY
	template< class VertexFactory , typename Index >
	void ReadPolygons
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		std::vector< typename VertexFactory::VertexType > &vertices ,
		std::vector< std::vector< Index > > &polygons ,
		bool *readFlags ,
		int &file_type ,
		std::vector< std::string > *comments
	)
#else // !USE_PLY_FACTORY
	template< class Vertex , typename Index >
	void ReadPolygons
	(
		std::string fileName ,
		std::vector< Vertex > &vertices ,
		std::vector< std::vector< Index > > &polygons ,
		const PlyProperty *properties ,
		bool *readFlags ,
		int propertyNum ,
		int &file_type ,
		std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
#else // !NEW_PLY
	template< class Vertex >
	void ReadPolygons
	(
		std::string fileName ,
		std::vector< Vertex > &vertices ,
		std::vector< std::vector< int > > &polygons ,
		const PlyProperty *properties ,
		bool *readFlags ,
		int propertyNum ,
		int &file_type ,
		std::vector< std::string > *comments
	)
#endif // NEW_PLY
	{
		std::vector< std::string > elist;
		float version;
		PlyFile *ply = PlyFile::Read( fileName , elist , file_type , version );
		if( !ply ) THROW( "could not create read ply file: " , fileName );

		if( comments )
		{
			comments->reserve( comments->size() + ply->comments.size() );
			for( int i=0 ; i<(int)ply->comments.size() ; i++ ) comments->push_back( ply->comments[i] );
		}

		for( int i=0 ; i<(int)elist.size() ; i++ )
		{
			std::string &elem_name = elist[i];
			size_t num_elems;
			std::vector< PlyProperty * > plist = ply->get_element_description( elem_name , num_elems );
			if( !num_elems ) continue;
			else if( !plist.size() )
			{
				delete ply;
				THROW( "could not get element description for: " , elem_name );
			}
			if( elem_name=="vertex" )
			{
#ifdef USE_PLY_FACTORY
				for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++)
				{
					PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(i) : vFactory.plyReadProperty(i);
					int hasProperty = ply->get_property( elem_name , &prop );
					if( readFlags ) readFlags[i] = (hasProperty!=0);
				}
				vertices.resize( num_elems , vFactory() );

				char *buffer = new char[ vFactory.bufferSize() ];
				for( size_t j=0 ; j<num_elems ; j++ )
				{
					if( vFactory.isStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
					else
					{
						ply->get_element( (void *)buffer );
						vFactory.fromBuffer( buffer , vertices[j] );
					}
				}
				delete[] buffer;
#else // !USE_PLY_FACTORY
				for( int i=0 ; i<propertyNum ; i++)
				{
					int hasProperty = ply->get_property( elem_name , &properties[i] );
					if( readFlags ) readFlags[i] = (hasProperty!=0);
				}
				vertices.resize( num_elems );
				for( int j=0 ; j<num_elems ; j++ ) ply->get_element( (void *)&vertices[j] );
#endif // USE_PLY_FACTORY
			}
			else if( elem_name=="face" )
			{
#ifdef NEW_PLY
				ply->get_property( elem_name , &Face< Index >::Properties[0] );
#else // !NEW_PLY
				ply->get_property( elem_name , &FaceProps[0] );
#endif // NEW_PLY
				polygons.resize( num_elems );
				for( unsigned int j=0 ; j<num_elems ; j++ )
				{
#ifdef NEW_PLY
					Face< Index > ply_face;
#else // !NEW_PLY
					Face ply_face;
#endif // NEW_PLY
					ply->get_element( (void *)&ply_face );
					polygons[j].resize( ply_face.nr_vertices );
					for( unsigned int k=0 ; k<ply_face.nr_vertices ; k++ ) polygons[j][k] = ply_face.vertices[k];
					free( ply_face.vertices );
				}  // for, read faces
			}  // if face
			else ply->get_other_element( elem_name , num_elems );

			for( int j=0 ; j<(int)plist.size() ; j++ ) delete plist[j];
		}  // for each type of element

		delete ply;
	}


#ifdef USE_PLY_FACTORY
	template< class VertexFactory , class Polygon >
	int ReadPolygons
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		std::vector< typename VertexFactory::VertexType >& vertices ,
		std::vector< Polygon >& polygons ,
		bool *vertexPropertiesFlag ,
		PlyProperty* polygonProperties , bool* polygonPropertiesFlag , int polygonPropertyNum ,
		int& file_type,
		std::vector< std::string > *comments
	)
#else // !USE_PLY_FACTORY
	template< class Vertex , class Polygon >
	int ReadPolygons
	(
		std::string fileName ,
		std::vector< Vertex >& vertices ,
		std::vector< Polygon >& polygons ,
		PlyProperty*  vertexProperties , bool*  vertexPropertiesFlag , int  vertexPropertyNum ,
		PlyProperty* polygonProperties , bool* polygonPropertiesFlag , int polygonPropertyNum ,
		int& file_type,
		std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
	{
		std::vector< std::string > elist = { std::string( "vertex" ) , std::string( "face" ) };
		float version;

		PlyFile *ply = PlyFile::Read( fileName , elist , file_type , version );
		if(!ply) return 0;

		if( comments )
		{
			comments->reserve( comments->size() + ply->comments.size() );
			for( int i=0 ; i<ply->comments.size() ; i++ ) comments->push_back( ply->comments[i] );
		}

		for( int i=0 ; i<elist.size() ; i++ )
		{
			std::string &elem_name = elist[i];
			size_t num_elems;
			std::vector< PlyProperty * > plist = ply->get_element_description( elem_name , num_elems );
			if( !plist.size() )
			{
				delete ply;
				return 0;
			}		
			if( elem_name=="vertex" )
			{
#ifdef USE_PLY_FACTORY
#if 1
				for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++)
				{
					PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(i) : vFactory.plyReadProperty(i);
					int hasProperty = ply->get_property( elem_name , &prop );
					if( vertexPropertiesFlag ) vertexPropertiesFlag[i] = (hasProperty!=0);
				}
				vertices.resize( num_elems , vFactory() );

				char *buffer = new char[ vFactory.bufferSize() ];
				for( size_t j=0 ; j<num_elems ; j++ )
				{
					if( vFactory.isStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
					else
					{
						ply->get_element( (void *)buffer );
						vFactory.fromBuffer( buffer , vertices[j] );
					}
				}
				delete[] buffer;
#else
				for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++)
				{
					PlyProperty property = vFactory.plyReadProperty( i );
					int hasProperty = ply->get_property( elem_name , &property );
					if( vertexPropertiesFlag ) vertexPropertiesFlag[i] = (hasProperty!=0);
				}
				vertices.resize( num_elems , vFactory() );
#endif
#else // !USE_PLY_FACTORY
				for( int i=0 ; i<vertexPropertyNum ; i++)
				{
					int hasProperty = ply->get_property( elem_name , &vertexProperties[i] );
					if( vertexPropertiesFlag ) vertexPropertiesFlag[i] = (hasProperty!=0);
				}
				vertices.resize( num_elems );
#endif // USE_PLY_FACTORY
				for( size_t j=0 ; j<num_elems ; j++ ) ply->get_element( (void *)&vertices[j] );
			}
			else if( elem_name=="face" )
			{
				for( int i=0 ; i<polygonPropertyNum ; i++ )
				{
					int hasProperty = ply->get_property( elem_name , &polygonProperties[i] );
					if( polygonPropertiesFlag ) polygonPropertiesFlag[i] = (hasProperty!=0);
				}
				polygons.resize( num_elems );
				for( size_t j=0 ; j<num_elems ; j++ ) ply->get_element( (void *)&polygons[j] );
			}
			else ply->get_other_element( elem_name , num_elems );

			for( int j=0 ; j<plist.size() ; j++ ) delete plist[j];
		}
		delete ply;
		return 1;
	}

#ifdef NEW_PLY
#ifdef USE_PLY_FACTORY
	template< class VertexFactory , typename Index >
	void ReadTetrahedra
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		std::vector< typename VertexFactory::VertexType > &vertices ,
		std::vector< SimplexIndex< 3 , Index > > &tetrahedra ,
		bool* vertexPropertiesFlag ,
		int& file_type ,
		std::vector< std::string > *comments
	)
#else // !USE_PLY_FACTORY
	template< class Vertex , typename Index >
	void ReadTetrahedra
	(
		std::string fileName ,
		std::vector< Vertex > &vertices ,
		std::vector< SimplexIndex< 3 , Index > > &tetrahedra ,
		const PlyProperty *vertexProperties ,
		bool* vertexPropertiesFlag ,
		int vertexPropertyNum ,
		int& file_type ,
		std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
#else // !NEW_PLY
	template< class Vertex >
	void ReadTetrahedra
	(
		std::string fileName ,
		std::vector< Vertex > &vertices ,
		std::vector< TetrahedronIndex > &tetrahedra ,
		const PlyProperty *vertexProperties ,
		bool* vertexPropertiesFlag ,
		int vertexPropertyNum ,
		int& file_type ,
		std::vector< std::string > *comments
	)
#endif // NEW_PLY
	{
#ifdef NEW_PLY
		std::vector< std::vector< Index > > polygons;
#else // !NEW_PLY
		std::vector< std::vector< int > > polygons;
#endif // NEW_PLY
#ifdef USE_PLY_FACTORY
		ReadPolygons( fileName , vFactory , vertices , polygons , vertexPropertiesFlag , file_type , comments );
#else // !USE_PLY_FACTORY
		ReadPolygons( fileName , vertices , polygons , vertexProperties , vertexPropertiesFlag , vertexPropertyNum , file_type , comments );
#endif // USE_PLY_FACTORY

		for( int i=0 ; i<polygons.size() ; i++ ) if( polygons[i].size()!=4 ) ERROR_OUT( "Expected polygon with four vertices" );
		tetrahedra.resize( polygons.size() );
		for( unsigned int i=0 ; i<polygons.size() ; i++ ) for( int j=0 ; j<4 ; j++ ) tetrahedra[i][j] = polygons[i][j];
	}

#ifdef NEW_PLY
#ifdef USE_PLY_FACTORY
	template< class VertexFactory , unsigned int Dim , typename Index >
	void ReadSimplices
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		std::vector< typename VertexFactory::VertexType > &vertices ,
		std::vector< SimplexIndex< Dim , Index > > &simplices ,
		bool* vertexPropertiesFlag ,
		int& file_type ,
		std::vector< std::string > *comments
	)

	{
		std::vector< std::vector< Index > > polygons;
		ReadPolygons( fileName , vFactory , vertices , polygons , vertexPropertiesFlag , file_type , comments );

		for( int i=0 ; i<polygons.size() ; i++ ) if( polygons[i].size()!=Dim+1 ) ERROR_OUT( "Expected polygon with " , Dim+1 , " vertices" );
		simplices.resize( polygons.size() );
		for( unsigned int i=0 ; i<polygons.size() ; i++ ) for( int j=0 ; j<=Dim ; j++ ) simplices[i][j] = polygons[i][j];
	}
#endif // USE_PLY_FACTORY
#endif // NEW_PLY

	// Write
#ifdef NEW_PLY
#ifdef USE_PLY_FACTORY
	template< typename VertexFactory , typename Index >
	void Write
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		const std::vector< typename VertexFactory::VertexType > &vertices , 
		const std::vector< std::pair< Index , Index > > *edges , 
		const std::vector< std::vector< Index > > *polygons,
		int file_type ,
		const std::vector< std::string > *comments
	)
#else // !USE_PLY_FACTORY
	template< class Vertex , typename Index >
	void Write
	(
		std::string fileName ,
		const std::vector< Vertex > &vertices , 
		const std::vector< std::pair< Index , Index > > *edges , 
		const std::vector< std::vector< Index > > *polygons,
		const PlyProperty *vertexProperties ,
		int vertexPropertyNum ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
#else // !NEW_PLY
	template< class Vertex >
	void Write
	(
		std::string fileName ,
		const std::vector< Vertex > &vertices , 
		const std::vector< std::pair< int , int > > *edges , 
		const std::vector< std::vector< int > > *polygons,
		const PlyProperty *vertexProperties ,
		int vertexPropertyNum ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#endif // NEW_PLY
	{
		int nr_vertices =            (int) vertices.size()    ;
		int nr_edges    = edges    ? (int)   edges->size() : 0;
		int nr_faces    = polygons ? (int)polygons->size() : 0;
		float version;
		std::vector< std::string > elist = { std::string( "vertex" ) , std::string( "edge" ) , std::string( "face" ) };

		PlyFile *ply = PlyFile::Write( fileName , elist , file_type , version );
		if( !ply ) THROW( "could not create write ply file: " , fileName );

		//
		// describe vertex, edge, and face properties
		//
		{
			ply->element_count( "vertex", nr_vertices );
#ifdef USE_PLY_FACTORY
			for( int i=0 ; i<vFactory.writeNum() ; i++ )
			{
				PlyProperty property = vFactory.writeProperty( i );
				ply->describe_property( "vertex" , &property );
			}
#else // !USE_PLY_FACTORY
			for( int i=0 ; i<vertexPropertyNum ; i++ ) ply->describe_property( "vertex" , &vertexProperties[i] );
#endif // USE_PLY_FACTORY
		}
		{
			ply->element_count( "edge" , nr_edges );
			ply->describe_property( "edge" , &EdgeProps[0] );
			ply->describe_property( "edge" , &EdgeProps[1] );
		}
		{
			ply->element_count( "face" , nr_faces );
#ifdef NEW_PLY
			ply->describe_property( "face" , &Face< Index >::Properties[0] );
#else // !NEW_PLY
			ply->describe_property( "face" , &FaceProps[0] );
#endif // NEW_PLY
		}

		// Write in the comments
		if( comments ) for( int i=0 ; i<comments->size() ; i++ ) ply->put_comment( (*comments)[i] );

		ply->header_complete();

		// write vertices
		ply->put_element_setup( "vertex" );
		for( int i=0 ; i<nr_vertices ; i++ ) ply->put_element( (void*)&vertices[i] );

		// write edges
		if( nr_edges )
		{
			Edge ply_edge;
			ply->put_element_setup( "edge" );
			for( int i=0 ; i<nr_edges ; i++ )
			{
				ply_edge.v1 = (*edges)[i].first , ply_edge.v2 = (*edges)[i].second;
				ply->put_element( (void*)&ply_edge );
			}
		}

		// write faces
		if( nr_faces )
		{
#ifdef NEW_PLY
			Face< Index > ply_face;
#else // !NEW_PLY
			Face ply_face;
#endif // NEW_PLY
			int maxFaceVerts=3;
			ply_face.nr_vertices = 3;
#ifdef NEW_PLY
			ply_face.vertices = new Index[3];
#else // !NEW_PLY
			ply_face.vertices = new int[3];
#endif // NEW_PLY

			ply->put_element_setup( "face" );
			for( int i=0 ; i<nr_faces ; i++ )
			{
				int face_size = (int)(*polygons)[i].size();
				if( face_size>maxFaceVerts )
				{
					delete[] ply_face.vertices;
					maxFaceVerts = face_size;
#ifdef NEW_PLY
					ply_face.vertices = new Index[face_size];
#else // !NEW_PLY
					ply_face.vertices = new int[face_size];
#endif // NEW_PLY
				}
				ply_face.nr_vertices = face_size;
				for( int j=0 ; j<ply_face.nr_vertices ; j++ ) ply_face.vertices[j] = (*polygons)[i][j];
				ply->put_element( (void*)&ply_face );
			}
			delete[] ply_face.vertices;
		}
		delete ply;
	}

#ifdef USE_PLY_FACTORY
		template< class VertexFactory >
		void WriteVertices
		(
			std::string fileName ,
			const VertexFactory &vFactory ,
			const std::vector< typename VertexFactory::VertexType > &vertices ,
			int file_type ,
			const std::vector< std::string > *comments
		)
#else // !USE_PLY_FACTORY
	template< class Vertex >
	void WriteVertices
	(
		std::string fileName ,
		const std::vector< Vertex > &vertices ,
		const PlyProperty *properties ,
		int propertyNum ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
	{
		int nr_vertices = (int)vertices.size();
		float version;
		std::vector< std::string > elem_names = { std::string( "vertex" ) };
		PlyFile *ply = PlyFile::Write( fileName , elem_names , file_type , version );
		if( !ply ) THROW( "could not create write ply file: " , fileName );

		//
		// describe vertex properties
		//
		ply->element_count( "vertex", nr_vertices );
#ifdef USE_PLY_FACTORY
		for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++ )
		{
			PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticWriteProperty(i) : vFactory.plyWriteProperty(i);
			ply->describe_property( "vertex" , &prop );
		}
#else // !USE_PLY_FACTORY
		for( int i=0 ; i<propertyNum ; i++ ) ply->describe_property( "vertex" , &properties[i] );
#endif // USE_PLY_FACTORY

		// Write in the comments
		if( comments ) for( int i=0 ; i<comments->size() ; i++ ) ply->put_comment( (*comments)[i] );
		ply->header_complete();

		// write vertices
		ply->put_element_setup( elem_names[0] );
#ifdef USE_PLY_FACTORY
		char *buffer = new char[ vFactory.bufferSize() ];
		for( size_t j=0 ; j<(int)vertices.size() ; j++ )
		{
			if( vFactory.isStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
			else
			{
				vFactory.toBuffer( vertices[j] , buffer );
				ply->put_element( (void *)buffer );
			}
		}
		delete[] buffer;
#else // !USE_PLY_FACTORY
		for( int i=0 ; i<(int)vertices.size() ; i++ ) ply->put_element( (void *)&vertices[i] );
#endif // USE_PLY_FACTORY

		delete ply;
	}

#ifdef NEW_PLY
#ifdef USE_PLY_FACTORY
	template< class VertexFactory , typename Index >
	void WriteTriangles
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		const std::vector< typename VertexFactory::VertexType > &vertices ,
		const std::vector< SimplexIndex< 2 , Index > > &triangles ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#else // !USE_PLY_FACTORY
	template< class Vertex , typename Index >
	void WriteTriangles
	(
		std::string fileName ,
		const std::vector< Vertex > &vertices ,
		const std::vector< SimplexIndex< 2 , Index > > &triangles ,
		const PlyProperty *properties ,
		int propertyNum ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
#else // !NEW_PLY
	template< class Vertex >
	void WriteTriangles
	(
		std::string fileName ,
		const std::vector< Vertex > &vertices ,
		const std::vector< TriangleIndex > &triangles ,
		const PlyProperty *properties ,
		int propertyNum ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#endif // NEW_PLY
	{
#ifdef NEW_PLY
		std::vector< std::vector< Index > > polygons( triangles.size() );
#else // !NEW_PLY
		std::vector< std::vector< int > > polygons( triangles.size() );
#endif // NEW_PLY
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			polygons[i].resize( 3 );
			for( int j=0 ; j<3 ; j++ ) polygons[i][j] = triangles[i][j];
		}
#ifdef USE_PLY_FACTORY
		WritePolygons( fileName , vFactory , vertices , polygons , file_type , comments );
#else // !USE_PLY_FACTORY
		WritePolygons( fileName , vertices , polygons , properties , propertyNum , file_type , comments );
#endif // USE_PLY_FACTORY
	}

#ifdef NEW_PLY
#ifdef USE_PLY_FACTORY
	template< class VertexFactory , typename Index >
	void WritePolygons
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		const std::vector< typename VertexFactory::VertexType > &vertices ,
		const std::vector< std::vector< Index > > &polygons ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#else // !USE_PLY_FACTORY
	template< class Vertex , typename Index >
	void WritePolygons
	(
		std::string fileName ,
		const std::vector< Vertex > &vertices ,
		const std::vector< std::vector< Index > > &polygons ,
		const PlyProperty *properties ,
		int propertyNum ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
#else // !NEW_PLY
	template< class Vertex >
	void WritePolygons
	(
		std::string fileName ,
		const std::vector< Vertex > &vertices ,
		const std::vector< std::vector< int > > &polygons ,
		const PlyProperty *properties ,
		int propertyNum ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#endif // NEW_PLY
	{
		int nr_vertices = int(vertices.size());
		int nr_faces = int(polygons.size());
		float version;
		std::vector< std::string > elem_names = { std::string( "vertex" ) , std::string( "face" ) };
		PlyFile *ply = PlyFile::Write( fileName , elem_names , file_type , version );
		if( !ply ) THROW( "could not create write ply file: " , fileName );

		//
		// describe vertex and face properties
		//
		ply->element_count( "vertex", nr_vertices );
#ifdef USE_PLY_FACTORY
#if 1
		for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++)
		{
			PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticWriteProperty(i) : vFactory.plyWriteProperty(i);
			ply->describe_property( "vertex" , &prop );
		}
#else
		for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++ )
		{
			PlyProperty property = vFactory.plyWriteProperty( i );
			ply->describe_property( "vertex" , &property );
		}
#endif
#else // !USE_PLY_FACTORY
		for( int i=0 ; i<propertyNum ; i++ ) ply->describe_property( "vertex" , &properties[i] );
#endif // USE_PLY_FACTORY
		ply->element_count( "face" , nr_faces );
#ifdef NEW_PLY
		ply->describe_property( "face" , &Face< Index >::Properties[0] );
#else // !NEW_PLY
		ply->describe_property( "face" , &FaceProps[0] );
#endif // NEW_PLY

		// Write in the comments
		if( comments ) for( size_t i=0 ; i<comments->size() ; i++ ) ply->put_comment( (*comments)[i] );
		ply->header_complete();

		// write vertices
		ply->put_element_setup( elem_names[0] );
#ifdef USE_PLY_FACTORY
		char *buffer = new char[ vFactory.bufferSize() ];
		for( size_t j=0 ; j<(int)vertices.size() ; j++ )
		{
			if( vFactory.isStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
			else
			{
				vFactory.toBuffer( vertices[j] , buffer );
				ply->put_element( (void *)buffer );
			}
		}
		delete[] buffer;
#else // !USE_PLY_FACTORY
		for( int i=0 ; i<(int)vertices.size() ; i++ ) ply->put_element( (void *)&vertices[i] );
#endif // USE_PLY_FACTORY

		// write faces
#ifdef NEW_PLY
		Face< Index > ply_face;
#else // !NEW_PLY
		Face ply_face;
#endif // NEW_PLY
		int maxFaceVerts = 3;
		ply_face.nr_vertices = maxFaceVerts;
#ifdef NEW_PLY
		ply_face.vertices = new Index[ maxFaceVerts ];
#else // !NEW_PLY
		ply_face.vertices = new int[ maxFaceVerts ];
#endif // NEW_PLY

		ply->put_element_setup( elem_names[1] );
		for( int i=0 ; i<nr_faces ; i++ )
		{
			if( (int)polygons[i].size()>maxFaceVerts )
			{
				delete[] ply_face.vertices;
				maxFaceVerts = (int)polygons[i].size();
#ifdef NEW_PLY
				ply_face.vertices=new Index[ maxFaceVerts ];
#else // !NEW_PLY
				ply_face.vertices=new int[ maxFaceVerts ];
#endif // NEW_PLY
			}
			ply_face.nr_vertices = (int)polygons[i].size();
			for( unsigned int j=0 ; j<ply_face.nr_vertices ; j++ ) ply_face.vertices[j] = polygons[i][j];
			ply->put_element( (void *)&ply_face );
		}
		delete[] ply_face.vertices;
		delete ply;
	}

#ifdef NEW_PLY
#ifdef USE_PLY_FACTORY
	template< class VertexFactory , typename Index >
	void WriteTetrahedra
	(
		std::string fileName ,
		const VertexFactory &vFactory ,
		const std::vector< typename VertexFactory::VertexType > &vertices ,
		const std::vector< SimplexIndex< 3 , Index > > &tetrahedra ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#else // !USE_PLY_FACTORY
	template< class Vertex , typename Index >
	void WriteTetrahedra
	(
		std::string fileName ,
		const std::vector< Vertex > &vertices ,
		const std::vector< SimplexIndex< 3 , Index > > &tetrahedra ,
		const PlyProperty *properties ,
		int propertyNum ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#endif // USE_PLY_FACTORY
#else // !NEW_PLY
	template< class Vertex >
	void WriteTetrahedra
	(
		std::string fileName ,
		const std::vector< Vertex > &vertices ,
		const std::vector< TetrahedronIndex > &tetrahedra ,
		const PlyProperty *properties ,
		int propertyNum ,
		int file_type ,
		const std::vector< std::string > *comments
	)
#endif // NEW_PLY
	{
#ifdef NEW_PLY
		std::vector< std::vector< Index > > polygons( tetrahedra.size() );
#else // !NEW_PLY
		std::vector< std::vector< int > > polygons( tetrahedra.size() );
#endif // NEW_PLY
		for( int i=0 ; i<tetrahedra.size() ; i++ )
		{
			polygons[i].resize( 4 );
			for( int j=0 ; j<4 ; j++ ) polygons[i][j] = tetrahedra[i][j];
		}
#ifdef USE_PLY_FACTORY
		WritePolygons( fileName , vFactory , vertices , polygons , file_type , comments );
#else // !USE_PLY_FACTORY
		WritePolygons( fileName , vertices , polygons , properties , propertyNum , file_type , comments );
#endif // USE_PLY_FACTORY
	}

	inline int DefaultFileType( void ){ return PLY_ASCII; }
}
