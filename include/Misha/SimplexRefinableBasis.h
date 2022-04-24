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

#ifndef SIMPLEX_REFINABLE_BASIS_INCLUDED
#define SIMPLEX_REFINABLE_BASIS_INCLUDED
#include <map>
#include <set>
#include <vector>
#include "Eigen/Dense"
#include "Misha/Exceptions.h"
#include "Misha/Geometry.h"
#include "SimplexBasis.h"
#include "LCQO.h"

#define INTERPOLATION_CONSTRAINTS

// An abstract class representing topology that can be decomposed into Dim-dimensional simplices
template< unsigned int Dim >
struct SimplexIndexRefinable
{
	// The number of sub-simplices
	virtual unsigned int size( void ) const = 0;

	// The idx-th sub-simplex
	virtual SimplexIndex< Dim , unsigned int > operator[]( unsigned int idx ) const = 0;
};

// An abstract class representing geometry that can be decomposed into Dim-dimensional simplices
template< unsigned int Dim >
struct SimplexRefinable : public SimplexIndexRefinable< Dim >
{
	// The metric on the idx-th sub-simplex.
	virtual SquareMatrix< double , Dim > metric( unsigned int idx ) const = 0;
};

// A (abstract) class that represents cells which:
//		Can be decomposed into Dim-dimensional simplices
//		Can enumerate the (Dim-1)-dimensional cells on their boundary
template< unsigned int Dim > struct SimplexRefinableCell;

// For the case Dim=1 the cell is simply an edge
template<>
struct SimplexRefinableCell<1> : public SimplexRefinable<1>
{
	static const unsigned int Dim = 1;
	SimplexRefinableCell( void ) : _si( SimplexIndex< Dim , unsigned int >(-1,-1) ) { _metric(0,0) = -1; }
	SimplexRefinableCell( SimplexIndex< Dim , unsigned int > si , double squareLength ) : _si(si) { _metric(0,0) = squareLength; }
	unsigned int size( void ) const { return 1u; }
	SimplexIndex< Dim , unsigned int > operator[]( unsigned int ) const { return _si; }
	SquareMatrix< double , Dim > metric( unsigned int ) const { return _metric; }
	void setNodeIndexSet( unsigned int cellDim , std::set< unsigned int > &nodeIndexSet ) const { nodeIndexSet.insert( _si[0] ) , nodeIndexSet.insert( _si[1] ); }
protected:
	SimplexIndex< Dim , unsigned int > _si;
	SquareMatrix< double , Dim > _metric;
};

// For the case Dim>1 the cells are more complex and this class abstract
template< unsigned int Dim >
struct SimplexRefinableCell : public SimplexRefinable< Dim >
{
	// Inherits the pure virtual method SimplexRefinable::metric( unsigned int ) const
	// Additional pure virtual methods of SimplexRefinableCell
	virtual unsigned int faces( void ) const = 0;
	virtual const SimplexRefinableCell< Dim-1 > &face( unsigned int faceIndex ) const = 0;
	virtual unsigned int centerIndex( void ) const = 0;

	// Implementation of SimplexRefinable::size
	unsigned int size( void ) const;
	// Implementation of SimplexRefinable::operator[]
	SimplexIndex< Dim , unsigned int > operator[]( unsigned int idx ) const;
	// Get the indices of the vertices by recursing through the faces
	void setNodeIndexSet( unsigned int cellDim , std::set< unsigned int > &nodeIndexSet ) const
	{
		if( Dim>cellDim ) for( unsigned int f=0 ; f<faces() ; f++ ) face(f).setNodeIndexSet( cellDim , nodeIndexSet );
		else
			for( unsigned int i=0 ; i<size() ; i++ )
			{
				SimplexIndex< Dim , unsigned int > s = (*this)[i];
				for( unsigned int d=0 ; d<=Dim ; d++ ) nodeIndexSet.insert( s[d] );
			}
	}
};

// Support for computing finite elements for simplex-refinable geometry
template< unsigned int ... Args > struct SimplexRefinableElements;

template<>
struct SimplexRefinableElements<>
{
	struct EnergyWeights
	{
		enum
		{
			MASS ,
			GRADIENT_SQUARE_NORM ,
			HESSIAN_SQUARE_NORM ,
			LAPLACIAN_SQUARE_NORM ,
			CROSS_FACE_GRADIENT_DIFFERENCE ,
			COUNT
		};
		static const std::string Names[];

		double weights[ COUNT ];
		EnergyWeights( void ){ memset( weights , 0 , sizeof(weights) ); }
		EnergyWeights( unsigned int eType ){ memset( weights , 0 , sizeof(weights) ) ; weights[eType] = 1; }
		double &operator[]( unsigned int idx ){ return weights[idx]; }
		const double &operator[]( unsigned int idx ) const { return weights[idx]; }
	};
};


template< unsigned int Dim >
struct SimplexRefinableElements< Dim >
{
	// The representation of a face (over which to integrate cross-face gradient differences)
	typedef MultiIndex< Dim , unsigned int , false > FaceMultiIndex;

	struct EnergyType
	{
		typename SimplexRefinableElements<>::EnergyWeights weights;
		std::function< bool ( FaceMultiIndex ) > isIntegrationFaceFunctor;
		EnergyType( void ) : isIntegrationFaceFunctor( [&]( FaceMultiIndex ){ return true; } ) {}
		EnergyType( typename SimplexRefinableElements<>::EnergyWeights eWeights ) : weights(eWeights) , isIntegrationFaceFunctor( [&]( FaceMultiIndex ){ return true; } ) {}
		EnergyType( unsigned int eType ) : weights(eType) , isIntegrationFaceFunctor( [&]( FaceMultiIndex ){ return true; } ) {}
	};
};

template< unsigned int Degree >
struct SimplexRefinableElements< 0 , Degree >
{
	// The representation of nodes at which basis functions are defined
	// [NOTE] The definition of NodeMultiIndex is independent of the value of Dim
	typedef MultiIndex< Degree , unsigned int , false > NodeMultiIndex;
};

template< unsigned int Dim , unsigned int Degree >
struct SimplexRefinableElements< Dim , Degree >
{
	// The representation of nodes at which basis functions are defined
	// [NOTE] The definition of NodeMultiIndex is independent of the value of Dim
	typedef typename SimplexRefinableElements< 0 , Degree >::NodeMultiIndex NodeMultiIndex;

	// The representation of a face (over which to integrate cross-face gradient differences)
	typedef typename SimplexRefinableElements< Dim >::FaceMultiIndex FaceMultiIndex;

	SimplexRefinableElements( const SimplexRefinable< Dim > &simplexRefinable );
	SimplexRefinableElements( SimplexRefinableElements &&sre ){ std::swap( _simplexRefinable , sre._simplexRefinable ) , std::swap( _nodeIndex , sre.nodeIndex ); }
	SimplexRefinableElements &operator = ( SimplexRefinableElements &&sre ){ std::swap( _simplexRefinable , sre._simplexRefinable ) , std::swap( _nodeIndex , sre.nodeIndex ) ; return *this; }

	// The number of nodes/basis-functions
	unsigned int size( void ) const { return _nodeIndex.size(); }

	// The multi-index describing the idx-th node
	const NodeMultiIndex &operator[]( unsigned int idx ) const { return _nodeIndex.toMultiIndex( idx ); }

	unsigned int nodeIndex( NodeMultiIndex nodeIndex ) const;
	unsigned int nodeIndex( SimplexIndex< Dim , unsigned int > subSimplexIndex , unsigned int nIdx ) const;

	Eigen::MatrixXd massMatrix( void ) const;
	Eigen::MatrixXd gradientSquareNormMatrix( void ) const;
	Eigen::MatrixXd hessianSquareNormMatrix( void ) const;
	Eigen::MatrixXd laplacianSquareNormMatrix( void ) const;
	virtual Eigen::MatrixXd crossFaceGradientDifferenceMatrix( void ) const { return crossFaceGradientDifferenceMatrix( [&]( FaceMultiIndex ){ return true; } ); }
	Eigen::MatrixXd systemMatrix( typename SimplexRefinableElements< Dim >::EnergyType eType ) const;
	template< unsigned int EmbeddingDimension >
	Eigen::MatrixXd systemMatrix( typename SimplexRefinableElements< Dim >::EnergyType eType , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor ) const;

	Eigen::MatrixXd crossFaceGradientDifferenceMatrix( std::function< bool ( FaceMultiIndex ) > faceSelectionFunctor ) const;

	struct NodeMultiIndex_Index
	{
		template< typename SubSimplexIndexFunctor >
		NodeMultiIndex_Index( SubSimplexIndexFunctor subSimplexIndexFunctor , unsigned int sz );
		unsigned int fromMultiIndex( SimplexRefinableElements::NodeMultiIndex nmi ) const;
		const SimplexRefinableElements::NodeMultiIndex &toMultiIndex( unsigned int idx ) const { return _nodeList[idx]; }
		unsigned int size( void ) const { return (unsigned int)_nodeMap.size(); }

		static SimplexRefinableElements::NodeMultiIndex NodeMultiIndex( SimplexIndex< Dim , unsigned int > subSimplexIndex , unsigned int n );

	protected:
		std::map< SimplexRefinableElements::NodeMultiIndex , unsigned int > _nodeMap;
		std::vector< SimplexRefinableElements::NodeMultiIndex > _nodeList;

	};

protected:
	const SimplexRefinable< Dim > &_simplexRefinable;
	NodeMultiIndex_Index _nodeIndex;

	template< typename SystemMatrixFunctor > Eigen::MatrixXd _systemMatrix( SystemMatrixFunctor F ) const;
};


// A struct for computing the prolongation matrix giving the expression of coarse basis function as a linear combination of finer basis function
// that minimizes the prescribed energy.
template< bool PoU >
struct InterpolatingProlongationSystem
{
	struct Constraint
	{
		Constraint( unsigned int cIndex=0 , unsigned int fIndex=0 , double v=0 ) : coarseIndex(cIndex) , fineIndex(fIndex) , value(v){}
		unsigned int coarseIndex , fineIndex;
		double value;
	};

	// Given a SPD matrix E and a subset of indices, construct the prolongation matrix taking coarse basis functions
	// to finer basis functions, that:
	// -- Is interpolating (if index idx is in the subset, it will be mapped to itself)
	// -- Forms a partition of unity

	// CoarseSelectionFunctor:
	//		A functor that takes a node index and returns a boolean indicating if that node is a degree of freedom in the coarse system
	template< typename CoarseSelectionFunctor >
	InterpolatingProlongationSystem( const Eigen::MatrixXd &E , CoarseSelectionFunctor coarseSelectionFunctor );
	InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices );
	InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , const std::vector< Constraint > &constraints );
#ifdef INTERPOLATION_CONSTRAINTS
	template< unsigned int InterpolationDim >
	InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , const Point< double , InterpolationDim > *interpolationConstraints  );
	template< unsigned int InterpolationDim >
	InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , const std::vector< Constraint > &constraints , const Point< double , InterpolationDim > *interpolationConstraints  );
#endif // INTERPOLATION_CONSTRAINTS

	template< unsigned int Degree >
	struct ProlongationInfo
	{
		Eigen::MatrixXd P;
		std::vector< typename SimplexRefinableElements< 0 , Degree >::NodeMultiIndex > coarseMultiIndices;
	};

	template< unsigned int Dim , unsigned int Degree >
	static void HierarchicalProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , typename SimplexRefinableElements<>::EnergyWeights eWeights , ProlongationInfo< Degree > pInfo[Dim] , unsigned int finestDim=Dim );

#ifdef INTERPOLATION_CONSTRAINTS
	template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
	static void HierarchicalProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , ProlongationInfo< Degree > pInfo[Dim] , unsigned int finestDim );
#endif // INTERPOLATION_CONSTRAINTS

	template< bool StableSolve >
	Eigen::MatrixXd prolongation( void ) const;
	double energy( const Eigen::MatrixXd &P ) const;

	const LCQO &lcqo( void ) const { return _lcqo; }
protected:
	unsigned int _coarseDim , _fineDim;
	Eigen::VectorXd _c;
	Eigen::SparseMatrix< double > _Q , _C;
	Eigen::VectorXd _q;
	std::vector< unsigned int > _coarseIndices , _fineIndices;
	LCQO _lcqo;

#ifdef INTERPOLATION_CONSTRAINTS
	template< unsigned int InterpolationDim >
	void _init( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , const std::vector< Constraint > &constraints , const Point< double , InterpolationDim > *interpolationConstraints );
#else // !INTERPOLATION_CONSTRAINTS
	void _init( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , const std::vector< Constraint > &constraints );
#endif // INTERPOLATION_CONSTRAINTS

	unsigned int _index( unsigned fine , unsigned int coarse ) const;
	Eigen::MatrixXd _toMatrix( const Eigen::VectorXd &v ) const;
	Eigen::VectorXd _toVector( const Eigen::MatrixXd &P ) const;

	template< unsigned int Dim , unsigned int Degree >
	static void _SetNodeMaps( const SimplexRefinableCell< Dim > &simplexRefinableCell , std::map< typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex , unsigned int > nodeMaps[Dim+1] );
	template< unsigned int Dim , unsigned int Degree >
	static void __SetNodeMaps( const SimplexRefinableCell< Dim > &simplexRefinableCell , std::map< typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex , unsigned int > nodeMaps[Dim+1] );

	template< unsigned int Dim , unsigned int Degree >
	static void _HierarchicalProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		ProlongationInfo< Degree > *pInfo , 
		unsigned int finestDim
	);

	template< unsigned int Dim , unsigned int Degree >
	static Eigen::MatrixXd _BoundaryProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		const std::vector< unsigned int > &coarseIndices
	);

	template< unsigned int Dim , unsigned int Degree >
	static Eigen::MatrixXd _BoundaryProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		const std::vector< unsigned int > &coarseIndices ,
		const Eigen::MatrixXd &coarseP
	);

#ifdef INTERPOLATION_CONSTRAINTS
	template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
	static void _HierarchicalProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor ,
		double planarityEpsilon ,
		ProlongationInfo< Degree > *pInfo , 
		unsigned int finestDim
	);

	template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
	static Eigen::MatrixXd _BoundaryProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor ,
		double planarityEpsilon ,
		const std::vector< unsigned int > &coarseIndices
	);

	template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
	static Eigen::MatrixXd _BoundaryProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor ,
		double planarityEpsilon ,
		const std::vector< unsigned int > &coarseIndices ,
		const Eigen::MatrixXd &coarseP
	);
#endif // INTERPOLATION_CONSTRAINTS
};

#include "SimplexRefinableBasis.inl"
#endif // SIMPLEX_REFINABLE_BASIS_INCLUDED