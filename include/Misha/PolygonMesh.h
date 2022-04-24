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

#ifndef POLYGON_MESH_INCLUDED
#define POLYGON_MESH_INCLUDED

#include "Polygon.h"
#include "SimplexMesh.h"
#include "SimplexRefinableMesh.h"
#include "SolidSimplexRefinableMesh.h"

namespace Meshes
{
	/////////////////
	// PolygonMesh //
	/////////////////
	template< typename VIndex >
	struct PolygonMesh
	{
		static const unsigned int Dim = 2;

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
				ERROR_OUT( "Bad index: " , idx );
				return Point< double , EmbeddingDimension >();
			};
		}

		unsigned int polygons( void ) const { return (unsigned int)_polygons.size(); }
		Polygon< VIndex > polygon( unsigned int p ) const
		{
			Polygon< VIndex > polygon( _polygons[p].size() );
			for( unsigned int v=0 ; v<_polygons[p].size() ; v++ ) polygon[v] = _vertices[ _polygons[p][v] ];
			return polygon;
		}

		unsigned int vertices( void ) const { return (unsigned int)_vertices.size(); }
		VIndex vertex( unsigned int v ) const { return _vertices[v]; }

		PolygonMesh( void ){}

		// Construct the polygon mesh from a list of polygons
		PolygonMesh( const std::vector< Polygon< VIndex > > &polygons )
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
#ifdef INTERPOLATION_CONSTRAINTS
			bool forceLinearPrecision ,
			double planarityEpsilon ,
			bool verbose
#else // !INTERPOLATION_CONSTRAINTS
			bool verbose=false
#endif // INTERPOLATION_CONSTRAINTS
		) const
		{
			std::function< SimplexRefinablePolygon (unsigned int) > cellFunctor = [&]( unsigned int c ){ return _simplexRefinable( c , fullVertexPositionFunction ); };
			typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinablePolygon > cellList( (unsigned int)_polygons.size() , cellFunctor );
#ifdef INTERPOLATION_CONSTRAINTS
			if( forceLinearPrecision )
				if( pou ) return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template Init< true  >( cellList , eWeights , fullVertexPositionFunction , planarityEpsilon , finestNodeDim , verbose );
				else      return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template Init< false >( cellList , eWeights , fullVertexPositionFunction , planarityEpsilon , finestNodeDim , verbose );
			else
				if( pou ) return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template Init< true  >( cellList , eWeights , finestNodeDim , verbose );
				else      return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template Init< false > ( cellList , eWeights , finestNodeDim , verbose );
#else // !INTERPOLATION_CONSTRAINTS
			return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::Init( cellList , eWeights , finestNodeDim , verbose );
#endif // INTERPOLATION_CONSTRAINTS
		}

		template< unsigned int Degree , unsigned int EmbeddingDimension >
		SimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh
		(
			FullVertexPositionFunction< EmbeddingDimension > fullVertexPositionFunction ,
			typename SimplexRefinableElements<>::EnergyWeights eWeights ,
			unsigned int finestNodeDim ,
			bool pou ,
#ifdef INTERPOLATION_CONSTRAINTS
			bool forceLinearPrecision ,
			double planarityEpsilon ,
			bool verbose
#else // !INTERPOLATION_CONSTRAINTS
			bool verbose=false
#endif // INTERPOLATION_CONSTRAINTS
		) const
		{
			std::function< SimplexRefinablePolygon (unsigned int) > cellFunctor = [&]( unsigned int c ){ return _simplexRefinable( c , fullVertexPositionFunction ); };
			typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinablePolygon > cellList( (unsigned int)_polygons.size() , cellFunctor );
#ifdef INTERPOLATION_CONSTRAINTS
			if( forceLinearPrecision )
				if( pou ) return SimplexRefinableCellMesh< Dim , Degree >::template Init< true  >( cellList , eWeights , fullVertexPositionFunction , planarityEpsilon , finestNodeDim , verbose );
				else      return SimplexRefinableCellMesh< Dim , Degree >::template Init< false >( cellList , eWeights , fullVertexPositionFunction , planarityEpsilon , finestNodeDim , verbose );
			else
				if( pou ) return SimplexRefinableCellMesh< Dim , Degree >::template Init< true  >( cellList , eWeights , finestNodeDim , verbose );
				else      return SimplexRefinableCellMesh< Dim , Degree >::template Init< false >( cellList , eWeights , finestNodeDim , verbose );
#else // !INTERPOLATION_CONSTRAINTS
			return SimplexRefinableCellMesh< Dim , Degree >::Init( cellList , eWeights , finestNodeDim , verbose );
#endif // INTERPOLATION_CONSTRAINTS
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
#ifdef INTERPOLATION_CONSTRAINTS
			bool forceLinearPrecision ,
			double planarityEpsilon ,
			bool verbose
#else // !INTERPOLATION_CONSTRAINTS
			bool verbose=false
#endif // INTERPOLATION_CONSTRAINTS
		) const
		{
			std::function< SimplexRefinablePolygon (unsigned int) > cellFunctor = [&]( unsigned int c ){ return _simplexRefinable( c , fullVertexPositionFunction ); };
			typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinablePolygon > cellList( (unsigned int)_polygons.size() , cellFunctor );
#ifdef INTERPOLATION_CONSTRAINTS
			if( pou ) return HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::template Init< true  >( cellList , fullVertexPositionFunction , eWeights , finestNodeDim , forceLinearPrecision , planarityEpsilon , verbose );
			else      return HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::template Init< false >( cellList , fullVertexPositionFunction , eWeights , finestNodeDim , forceLinearPrecision , planarityEpsilon , verbose );
#else // !INTERPOLATION_CONSTRAINTS
			return HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::Init( cellList , fullVertexPositionFunction , eWeights , finestNodeDim , verbose );
#endif // INTERPOLATION_CONSTRAINTS
		}

		template< unsigned int Degree >
		SolidSimplexRefinableCellMesh< Dim , Degree > solidSimplexRefinableCellMesh
		(
			FullVertexPositionFunction< Dim > fullVertexPositionFunction ,
			typename SimplexRefinableElements<>::EnergyWeights eWeights ,
			unsigned int finestNodeDim ,
			bool pou ,
#ifdef INTERPOLATION_CONSTRAINTS
			bool forceLinearPrecision ,
			double planarityEpsilon ,
			bool verbose
#else // !INTERPOLATION_CONSTRAINTS
			bool verbose=false
#endif // INTERPOLATION_CONSTRAINTS
		) const
		{
			std::function< SimplexRefinablePolygon (unsigned int) > cellFunctor = [&]( unsigned int c ){ return _simplexRefinable( c , fullVertexPositionFunction ); };
			typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinablePolygon > cellList( (unsigned int)_polygons.size() , cellFunctor );
#ifdef INTERPOLATION_CONSTRAINTS
			if( pou ) return SolidSimplexRefinableCellMesh< Dim , Degree >::template Init< true  >( cellList , fullVertexPositionFunction , eWeights , finestNodeDim , forceLinearPrecision , planarityEpsilon , verbose );
			else      return SolidSimplexRefinableCellMesh< Dim , Degree >::template Init< false >( cellList , fullVertexPositionFunction , eWeights , finestNodeDim , forceLinearPrecision , planarityEpsilon , verbose );
#else // !INTERPOLATION_CONSTRAINTS
			return SolidSimplexRefinableCellMesh< Dim , Degree >::Init( cellList , fullVertexPositionFunction , eWeights , finestNodeDim , verbose );
#endif // INTERPOLATION_CONSTRAINTS
		}

	protected:
		std::vector< VIndex > _vertices;
		std::vector< Polygon< unsigned int > > _polygons;

		template< unsigned int EmbeddingDimension >
		SimplexRefinablePolygon _simplexRefinable( unsigned int p , const FullVertexPositionFunction< EmbeddingDimension > &vFunction ) const
		{
			return SimplexRefinablePolygon( _polygons[p] , (unsigned int)_vertices.size() + p , vFunction );
		}
	};

}

#endif // POLYGON_MESH_INCLUDED