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

#ifndef POLYHEDRON_INCLUDED
#define POLYHEDRON_INCLUDED

#include <functional>
#include "Polygon.h"
#include "SimplexRefinableBasis.h"

namespace Meshes
{
	// Aliasing a polyhedron as a vector of polygons
	// (Not checking for consistency)
	template< typename VIndex > using Polyhedron = std::vector< std::vector< VIndex > >;

	// Function for returning the point minimizing the squared are of the triangulation
	template< unsigned int EmbeddingDimension , typename VIndex >
	Point< double , EmbeddingDimension > SquaredVolumeMinimizingCenter( const Polyhedron< VIndex > &polyhedron , const VertexPositionFunction< EmbeddingDimension , VIndex > &vertexPositionFunction , const VertexPositionFunction< EmbeddingDimension , unsigned int > &facePositionFunction )
	{
		typedef Point< double , EmbeddingDimension > Pt;
		auto OuterProduct = []( Pt p , Pt q )
		{
			SquareMatrix< double , EmbeddingDimension > M;
			for( unsigned int i=0 ; i<EmbeddingDimension ; i++ ) for( unsigned int j=0 ; j<EmbeddingDimension ; j++ ) M(i,j) = p[i] * q[j];
			return M;
		};
		SquareMatrix< double , EmbeddingDimension > M;
		Pt b;
		for( unsigned int f=0 ; f<polyhedron.size() ; f++ )
		{
			const Polygon< VIndex > &polygon = polyhedron[f];
			Pt c = facePositionFunction( f );
			Pt p[2];
			p[0] = vertexPositionFunction( polygon[0] ) - c;
			for( unsigned int e=0 ; e<polygon.size() ; e++ )
			{
				SquareMatrix< double , EmbeddingDimension > _M;
				p[1] = vertexPositionFunction( polygon[ (e+1)%polygon.size() ] ) - c;
				double inner[2][2];
				SquareMatrix< double , EmbeddingDimension > outer[2][2];
				for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int j=0 ; j<2 ; j++ ) inner[i][j] = Pt::Dot( p[i] , p[j] ) , outer[i][j]  = OuterProduct( p[i] , p[j] );

				_M = SquareMatrix< double , EmbeddingDimension >::Identity() * ( inner[0][0] * inner[1][1] - inner[0][1] * inner[1][0] )
					- outer[0][0] * inner[1][1] - outer[1][1] * inner[0][0]
					+ outer[0][1] * inner[1][0] + outer[1][0] * inner[0][1];
				b += _M * c;
				M += _M;

				p[0] = p[1];
			}
		}
		return M.inverse() * b;
	}
	template< unsigned int EmbeddingDimension , typename VIndex >
	Point< double , EmbeddingDimension > CenterOfMass( const Polyhedron< VIndex > &polyhedron , const VertexPositionFunction< EmbeddingDimension , VIndex > &vertexPositionFunction , const VertexPositionFunction< EmbeddingDimension , unsigned int > &facePositionFunction )
	{
		typedef Point< double , EmbeddingDimension > Pt;
		Pt c;
		double area = 0;
		for( unsigned int f=0 ; f<polyhedron.size() ; f++ )
		{
			const Polygon< VIndex > &polygon = polyhedron[f];
			Pt _c = facePositionFunction( f );
			Pt p1 = vertexPositionFunction( polygon[0] ) - _c;
			for( unsigned int e=0 ; e<polygon.size() ; e++ )
			{
				Pt p2 = vertexPositionFunction( polygon[ (e+1)%polygon.size() ] ) - _c;
				double a = sqrt( Pt::SquareNorm(p1) * Pt::SquareNorm(p2) - Pt::Dot(p1,p2) * Pt::Dot(p2,p1) );
				c += ( _c*3. + p1 + p2 ) * a / 3.;
				area += a;
				p1 = p2;
			}
		}
		return c / area;
	}


	////////////////////////////////
	// SimplexRefinablePolyhedron //
	////////////////////////////////

	// This structure represents a single polygon that can be refined into simplices (triangles)
	struct SimplexRefinablePolyhedron : public SimplexRefinableCell< 3 >
	{
		static const unsigned int Dim = 3;
		typedef SimplexRefinablePolygon SimplexRefinableBoundaryCellType;

		// Constructs a simplex-refinable polyhedron given:
		//		the polyhedron
		//		the indices of the (virtual) polyhedron and face centers
		//		the function embedding the vertices, face, and polyhedron centers
		template< unsigned int EmbeddingDimension >
		SimplexRefinablePolyhedron( const Polyhedron< unsigned int > &p , unsigned int cellCenterIndex , const std::vector< unsigned int > &faceCenterIndices , const FullVertexPositionFunction< EmbeddingDimension > &fullVertexPositionFunction );

		// Returns the number of (Dim-1)-dimensional cells on the boundary
		unsigned int faces( void ) const { return (unsigned int)_faces.size(); }

		// Returns the prescribed face from the boundary
		const SimplexRefinableCell< Dim-1 > &face( unsigned int faceIndex ) const{ return _faces[faceIndex]; }

		// Returns the index of the (virtual) polygon center
		unsigned int centerIndex( void ) const { return _centerIndex; }

		// Returns the metric associated with s-th simplex
		SquareMatrix< double , Dim > metric( unsigned int s ) const { return _metrics[s]; }
	protected:
		unsigned int _centerIndex;
		std::vector< SimplexRefinablePolygon > _faces;
		std::vector< SquareMatrix< double , Dim > > _metrics;
	};

	template< unsigned int EmbeddingDimension >
	SimplexRefinablePolyhedron::SimplexRefinablePolyhedron( const Polyhedron< unsigned int > &p , unsigned int cellCenterIndex , const std::vector< unsigned int > &faceCenterIndices , const FullVertexPositionFunction< EmbeddingDimension > &fullVertexPositionFunction )
		: _centerIndex(cellCenterIndex)
	{
		// First set SimplexRefinablePolygon::_faces
		_faces.resize( p.size() );
		for( unsigned int i=0 ; i<p.size() ; i++ ) _faces[i] = SimplexRefinablePolygon( p[i] , faceCenterIndices[i] , fullVertexPositionFunction );

		// Now that SimplexRefinablePolygon::_faces is set, we can invoke
		//		SimplexRefinableCell< 3 >::operator[], which invokes
		//		SimplexRefinablePolyhedron::face method, which uses
		//		SimplexRefinablePolyhedron::_faces,
		// to get the simplex indices, to get the simplices, to set the metrics
		_metrics.resize( size() );
		for( unsigned int i=0 ; i<size() ; i++ )
		{
			SimplexIndex< Dim , unsigned int > simplexIndex = (*this)[i];
			Simplex< double , EmbeddingDimension , Dim > simplex;
			for( unsigned int j=0 ; j<=Dim ; j++ ) simplex[j] = fullVertexPositionFunction( simplexIndex[j] );
			_metrics[i] = RightSimplex< Dim >::Metric( simplex );
		}
	}
}

#endif // POLYHEDRON_INCLUDED