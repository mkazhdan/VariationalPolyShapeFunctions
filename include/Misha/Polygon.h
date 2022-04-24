#ifndef POLYGON_INCLUDED
#define POLYGON_INCLUDED

#include <functional>
#include "SimplexRefinableBasis.h"
#include "Misha/Geometry.h"


namespace Meshes
{
	// Aliasing a polygon as a vector of indices
	template< typename VIndex > using Polygon = std::vector< VIndex >;

	// Function for returning the point minimizing the squared are of the triangulation
	template< unsigned int EmbeddingDimension , typename VIndex >
	Point< double , EmbeddingDimension > SquaredAreaMinimizingCenter( const Polygon< VIndex > &polygon , const VertexPositionFunction< EmbeddingDimension , VIndex > &vertexPositionFunction )
	{
		auto OuterProduct = []( Point< double , EmbeddingDimension > p )
		{
			SquareMatrix< double , EmbeddingDimension > M;
			for( unsigned int i=0 ; i<EmbeddingDimension ; i++ ) for( unsigned int j=0 ; j<EmbeddingDimension ; j++ ) M(i,j) = p[i] * p[j];
			return M;
		};
		SquareMatrix< double , EmbeddingDimension > M;
		Point< double , EmbeddingDimension > b;
		Point< double , EmbeddingDimension > p = vertexPositionFunction( polygon[0] );
		for( unsigned int e=0 ; e<polygon.size() ; e++ )
		{
			Point< double , EmbeddingDimension > _p = vertexPositionFunction( polygon[ (e+1)%polygon.size() ] );
			Point< double , EmbeddingDimension > d = _p - p;
			double innerProduct = Point< double , EmbeddingDimension >::SquareNorm( d );
			SquareMatrix< double , EmbeddingDimension > outerProduct = OuterProduct( d );
			M += SquareMatrix< double , EmbeddingDimension >::Identity() * innerProduct - outerProduct;
			b += p * innerProduct - outerProduct * p;
			p = _p;
		}
		return M.inverse() * b;
	}

	template< unsigned int EmbeddingDimension , typename VIndex >
	Point< double , EmbeddingDimension > CenterOfMass( const Polygon< VIndex > &polygon , const VertexPositionFunction< EmbeddingDimension , VIndex > &vertexPositionFunction )
	{
		Point< double , EmbeddingDimension > c;
		double l = 0;
		Point< double , EmbeddingDimension > p = vertexPositionFunction( polygon[0] );
		for( unsigned int e=0 ; e<polygon.size() ; e++ )
		{
			Point< double , EmbeddingDimension > _p = vertexPositionFunction( polygon[ (e+1)%polygon.size() ] );
			double _l = Point< double , EmbeddingDimension >::Length( p - _p );
			c += ( p + _p ) * _l / 2.;
			l += _l;
			p = _p;
		}
		return c / l;
	}

	/////////////////////////////
	// SimplexRefinablePolygon //
	/////////////////////////////

	// This structure represents a single polygon that can be refined into simplices (triangles)
	struct SimplexRefinablePolygon : public SimplexRefinableCell< 2 >
	{
		static const unsigned int Dim = 2;
		typedef SimplexRefinableCell< Dim-1 > SimplexRefinableBoundaryCellType;

		SimplexRefinablePolygon( void ){}

		// Constructs a simplex-refinable polygon given:
		//		the polygon
		//		the index of the (virtual) polygon center
		//		the function embedding the vertices and polygon centers
		template< unsigned int EmbeddingDimension >
		SimplexRefinablePolygon( const Polygon< unsigned int > &p , unsigned int centerIndex , const FullVertexPositionFunction< EmbeddingDimension > &fullVertexPositionFunction );

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
		std::vector< SimplexRefinableCell< Dim-1 > > _faces;
		std::vector< SquareMatrix< double , Dim > > _metrics;
	};

	template< unsigned int EmbeddingDimension >
	SimplexRefinablePolygon::SimplexRefinablePolygon( const Polygon< unsigned int > &p , unsigned int centerIndex , const FullVertexPositionFunction< EmbeddingDimension > &fullVertexPositionFunction )
		: _centerIndex(centerIndex)
	{
		// First set SimplexRefinablePolygon::_faces
		_faces.resize( p.size() );
		for( unsigned int i=0 ; i<p.size() ; i++ )
		{
			unsigned int v1 = p[i] , v2 = p[ (i+1)%p.size() ];
			double squareLength = Point< double , EmbeddingDimension >::SquareNorm( fullVertexPositionFunction(v1) - fullVertexPositionFunction(v2) );
			_faces[i] = SimplexRefinableCell< 1 >( SimplexIndex< Dim-1 , unsigned int >( v1 , v2 ) , squareLength );
		}

		// Now that SimplexRefinablePolygon::_faces is set, we can invoke
		//		SimplexRefinableCell< 2 >::operator[], which invokes
		//		SimplexRefinablePolygon::face method, which uses
		//		SimplexRefinablePolygon::_faces,
		// to get the simplex indices, to get the simplices, to set the metrics
		_metrics.resize( size() );
		for( unsigned int i=0 ; i<p.size() ; i++ )
		{
			SimplexIndex< Dim , unsigned int > simplexIndex = (*this)[i];
			Simplex< double , EmbeddingDimension , Dim > simplex;
			for( unsigned int j=0 ; j<=Dim ; j++ ) simplex[j] = fullVertexPositionFunction( simplexIndex[j] );
			_metrics[i] = RightSimplex< Dim >::Metric( simplex );
		}
	}
}

#endif // POLYGON_INCLUDED