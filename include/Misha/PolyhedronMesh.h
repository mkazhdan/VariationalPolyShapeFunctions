#ifndef POLYHEDRON_MESH_INCLUDED
#define POLYHEDRON_MESH_INCLUDED

#include "Polygon.h"
#include "Polyhedron.h"
#include "SimplexMesh.h"
#include "SimplexRefinableMesh.h"
#include "SolidSimplexRefinableMesh.h"

namespace Meshes
{
	////////////////////
	// PolyhedronMesh //
	////////////////////
	template< typename VIndex >
	struct PolyhedronMesh
	{
		static const unsigned int Dim = 3;

		template< unsigned int EmbeddingDim >
		static void ReadOVM( const std::string &fileName , std::vector< std::vector< std::pair< unsigned int , bool > > > &polyhedra , std::vector< Meshes::Polygon< unsigned int > > &polygons , std::vector< Point< double , EmbeddingDim > > &vertices )
		{
			typedef std::pair< unsigned int , unsigned int > Edge;
			std::vector< Edge > edges;
			std::ifstream file( fileName );
			if( !file.is_open() ) ERROR_OUT( "Could not open file for reading: " , fileName );

			auto CheckString = [&]( const std::string& s )
			{
				std::string buf;
				file >> buf;
				if( s.compare(buf) )
				{
					std::cout << "reading ovm: string " << s << " not found. Returning." << std::endl;
					return false;
				}
				return true;
			};

			if( !CheckString( "OVM" ) ) ERROR_OUT( "Failed to parse 'OVM'" );
			if( !CheckString( "ASCII" ) ) ERROR_OUT( "Failed to parse 'ASCII'" );

			{
				if( !CheckString( "Vertices" ) ) ERROR_OUT( "Failed to parse 'Vertices'" );

				int nv;
				file >> nv;
				vertices.reserve( nv );

				for( int i=0 ; i<nv ; i++ )
				{
					Point< double , EmbeddingDim > p;
					for( unsigned int d=0 ; d<EmbeddingDim ; d++ ) file >> p[d];
					vertices.push_back( p );
				}
			}

			{
				if( !CheckString( "Edges" ) ) ERROR_OUT( "Could not find 'Edges'" );
				int ne;
				file >> ne;
				edges.reserve( ne );

				Edge e;
				for( int i=0 ; i<ne ; i++ )
				{
					file >> e.first >> e.second;
					edges.push_back(e);
				}
			}

			{
				if( !CheckString( "Faces" ) ) ERROR_OUT( "Could not find 'Faces'" );
				int nf;
				file >> nf;
				polygons.resize( nf );

				for( int i=0 ; i<nf ; i++ )
				{
					int nei;
					file >> nei;
					polygons[i].reserve(nei);
					unsigned int id;

					for( int j=0 ; j<nei ; j++ )
					{
						file >> id;
						if( id&1 ) polygons[i].push_back( edges[ id>>1 ].second );
						else       polygons[i].push_back( edges[ id>>1 ].first );
					}
				}
			}

			{
				if( !CheckString( "Polyhedra" ) ) ERROR_OUT( "Could not find 'Polyhedra'" );
				int nc;
				file >> nc;
				polyhedra.resize( nc );

				for( int i=0 ; i<nc ; i++ )
				{
					int nfi;
					file >> nfi;
					polyhedra[i].reserve( nfi );

					unsigned int id;
					for( int j=0 ; j<nfi ; j++ )
					{
						file >> id;
						polyhedra[i].push_back( std::pair< unsigned int , bool >( id>>1 , !(id&1) ) );
					}
				}
			}

			file.close();
		}

		// Given a function giving the position of the vertices, generate the function giving the position of the center of all faces
		template< unsigned int EmbeddingDimension >
		FullVertexPositionFunction< EmbeddingDimension > fullVertexPositionFunction( const VertexPositionFunction< EmbeddingDimension , VIndex > &vertexPositionFunction , bool squaredAreaMinimizer ) const
		{
			return [&,squaredAreaMinimizer]( unsigned int idx )
			{
				if( idx<_vertices.size() ) return vertexPositionFunction( _vertices[idx] );
				idx -= (unsigned int)_vertices.size();
				if( idx<_polygons.size() )
				{
					if( squaredAreaMinimizer ) return SquaredAreaMinimizingCenter( polygon(idx) , vertexPositionFunction );
					else                       return                CenterOfMass( polygon(idx) , vertexPositionFunction );
				}
				idx -= (unsigned int)_polygons.size();
				if( idx<_polyhedra.size() )
				{
					VertexPositionFunction< EmbeddingDimension , unsigned int > facePositionFunction = [&]( unsigned int f ){ return CenterOfMass( polygon( _polyhedra[idx][f].first ) , vertexPositionFunction ); };

					if( squaredAreaMinimizer ) return SquaredVolumeMinimizingCenter( polyhedron(idx) , vertexPositionFunction , facePositionFunction );
					else                       return                  CenterOfMass( polyhedron(idx) , vertexPositionFunction , facePositionFunction );
				}
				ERROR_OUT( "Bad index: " , idx );
				return Point< double , EmbeddingDimension >();
			};
		}

		unsigned int polyhedra( void ) const { return (unsigned int)_polyhedra.size(); }
		Polyhedron< VIndex > polyhedron( unsigned int p ) const
		{
			const std::vector< std::pair< unsigned int , bool > > &_polyhedron = _polyhedra[p];
			Polyhedron< VIndex > polyhedron( _polyhedron.size() );
			for( unsigned int i=0 ; i<_polyhedron.size() ; i++ )
			{
				const Polygon< unsigned int > &_face = _polygons[ _polyhedron[i].first ];
				Polygon< VIndex > &face = polyhedron[i];
				face.resize( _face.size() );
				if( _polyhedron[i].second ) for( unsigned int j=0 ; j<face.size() ; j++ ) face[j] = _vertices[ _face[j] ];
				else                        for( unsigned int j=0 ; j<face.size() ; j++ ) face[j] = _vertices[ _face[ _face.size()-1-j ] ];
			}
			return polyhedron;
		};

		unsigned int polyhedronFaces( unsigned int p ) const { return (unsigned int)_polyhedra[p].size(); }
		std::pair< unsigned int , bool > polyhedronFace( unsigned int p , unsigned int f ) const { return _polyhedra[p][f]; }

		unsigned int polygons( void ) const { return (unsigned int)_polygons.size(); }
		Polygon< VIndex > polygon( unsigned int p ) const
		{
			Polygon< VIndex > polygon( _polygons[p].size() );
			for( unsigned int v=0 ; v<_polygons[p].size() ; v++ ) polygon[v] = _vertices[ _polygons[p][v] ];
			return polygon;
		}

		unsigned int vertices( void ) const { return (unsigned int)_vertices.size(); }
		VIndex vertex( unsigned int v ) const { return _vertices[v]; }

		PolyhedronMesh( void ){}

		// Construct the polyhedron mesh from a list of polyhedra and the associated faces
		PolyhedronMesh( const std::vector< std::vector< std::pair< unsigned int , bool > > > &polyhedra , const std::vector< Polygon< VIndex > > &polygons )
		{
			std::map< VIndex , unsigned int > vMap;
			for( unsigned int i=0 ; i<polygons.size() ; i++ ) for( unsigned int j=0 ; j<polygons[i].size() ; j++ ) vMap[ polygons[i][j] ] = 0;

			_vertices.resize( vMap.size() );
			unsigned int count = 0;
			for( auto & [ vIdx , idx ] : vMap )
			{
				_vertices[count] = vIdx;
				idx = count++;
			}
			_polygons.resize( polygons.size() );
			for( unsigned int i=0 ; i<polygons.size() ; i++ )
			{
				_polygons[i].resize( polygons[i].size() );
				for( unsigned int j=0 ; j<polygons[i].size() ; j++ ) _polygons[i][j] = vMap[ polygons[i][j] ];
			}
			_polyhedra = polyhedra;
		}

		template< unsigned int EmbeddingDimension >
		SimplexRefinablePolyhedron simplexRefinable( unsigned int p , const FullVertexPositionFunction< EmbeddingDimension > &vFunction ) const
		{
			std::vector< unsigned int > faceCenterIndices( _polyhedra[p].size() );
			for( unsigned int i=0 ; i<_polyhedra[p].size() ; i++ ) faceCenterIndices[i] = (unsigned int)_vertices.size() + _polyhedra[p][i].first;
			return SimplexRefinablePolyhedron( polyhedron(p) , (unsigned int)( _vertices.size() + _polygons.size() ) + p , faceCenterIndices , vFunction );
		}

		// Create a simplex-refinable cell mesh, with:
		//		metric defined by the embedding,
		//		prolongation defining basis function at the vertices, defined by:
		//			iterating over the K-dimensional meshes (with K going from lower to higher)
		//				using the energy weights to define a prolongation coefficients for nodes on the K-dimensional meshes that fixes the coefficients from the (K-1)-dimensional solution
		template< unsigned int Degree , unsigned int EmbeddingDimension >
		HierarchicalSimplexRefinableCellMesh< Dim , Degree > hierarchicalSimplexRefinableCellMesh
		(
			FullVertexPositionFunction< EmbeddingDimension > fullVertexPositionFunction ,
			typename SimplexRefinableElements<>::EnergyWeights eWeights ,
			unsigned int finestNodeDim ,
			bool pou , 
			bool forceLinearPrecision ,
			double planarityEpsilon ,
			bool verbose
		) const
		{
			std::function< SimplexRefinablePolyhedron (unsigned int) > cellFunctor = [&]( unsigned int c ){ return simplexRefinable( c , fullVertexPositionFunction ); };
			typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinablePolyhedron > cellList( (unsigned int)_polyhedra.size() , cellFunctor );
			if( forceLinearPrecision ) return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::Init( cellList , eWeights , pou , fullVertexPositionFunction , planarityEpsilon , finestNodeDim , verbose );
			else                       return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::Init( cellList , eWeights , pou , finestNodeDim , verbose );
		}

		template< unsigned int Degree , unsigned int EmbeddingDimension >
		SimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh
		(
			FullVertexPositionFunction< EmbeddingDimension > fullVertexPositionFunction ,
			typename SimplexRefinableElements<>::EnergyWeights eWeights ,
			unsigned int finestNodeDim ,
			bool pou , 
			bool forceLinearPrecision ,
			double planarityEpsilon ,
			bool verbose
		) const
		{
			std::function< SimplexRefinablePolyhedron (unsigned int) > cellFunctor = [&]( unsigned int c ){ return simplexRefinable( c , fullVertexPositionFunction ); };
			typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinablePolyhedron > cellList( (unsigned int)_polyhedra.size() , cellFunctor );
			if( forceLinearPrecision ) return SimplexRefinableCellMesh< Dim , Degree >::Init( cellList , eWeights , pou , fullVertexPositionFunction , planarityEpsilon , finestNodeDim , verbose );
			else                       return SimplexRefinableCellMesh< Dim , Degree >::Init( cellList , eWeights , pou , finestNodeDim , verbose );
		}

		// Create a solid simplex-refinable cell mesh, with:
		//		metric defined by the embedding,
		//		prolongation defining basis function at the vertices, defined by:
		//			iterating over the K-dimensional meshes (with K going from lower to higher)
		//				using the energy weights to define a prolongation coefficients for nodes on the K-dimensional meshes that fixes the coefficients from the (K-1)-dimensional solution
		template< unsigned int Degree >
		HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree > hierarchicalSolidSimplexRefinableCellMesh
		(
			FullVertexPositionFunction< Dim > fullVertexPositionFunction ,
			typename SimplexRefinableElements<>::EnergyWeights eWeights ,
			unsigned int finestNodeDim ,
			bool pou ,
			bool forceLinearPrecision ,
			double planarityEpsilon ,
			bool verbose
		) const
		{
			std::function< SimplexRefinablePolyhedron (unsigned int) > cellFunctor = [&]( unsigned int c ){ return simplexRefinable( c , fullVertexPositionFunction ); };
			typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinablePolyhedron > cellList( (unsigned int)_polyhedra.size() , cellFunctor );
			return HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::Init( cellList , fullVertexPositionFunction , eWeights , finestNodeDim , pou , forceLinearPrecision , planarityEpsilon , verbose );
		}

		template< unsigned int Degree >
		SolidSimplexRefinableCellMesh< Dim , Degree > solidSimplexRefinableCellMesh
		(
			FullVertexPositionFunction< Dim > fullVertexPositionFunction ,
			typename SimplexRefinableElements<>::EnergyWeights eWeights ,
			unsigned int finestNodeDim ,
			bool pou ,
			bool forceLinearPrecision ,
			double planarityEpsilon ,
			bool verbose
		) const
		{
			std::function< SimplexRefinablePolyhedron (unsigned int) > cellFunctor = [&]( unsigned int c ){ return simplexRefinable( c , fullVertexPositionFunction ); };
			typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinablePolyhedron > cellList( (unsigned int)_polyhedra.size() , cellFunctor );
			return SolidSimplexRefinableCellMesh< Dim , Degree >::Init( cellList , fullVertexPositionFunction , eWeights , finestNodeDim , pou , forceLinearPrecision , planarityEpsilon , verbose );
		}

	protected:
		std::vector< VIndex > _vertices;
		std::vector< Polygon< unsigned int > > _polygons;
		std::vector< std::vector< std::pair< unsigned int , bool > > > _polyhedra;
	};

}

#endif // POLYHEDRON_MESH_INCLUDED