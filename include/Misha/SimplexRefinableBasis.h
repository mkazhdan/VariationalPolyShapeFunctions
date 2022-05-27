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
#ifdef USE_UNORDERED_SET_MAP
	void setNodeIndexSet( unsigned int cellDim , std::unordered_set< unsigned int > &nodeIndexSet ) const { nodeIndexSet.insert( _si[0] ) , nodeIndexSet.insert( _si[1] ); }
#else // !USE_UNORDERED_SET_MAP
	void setNodeIndexSet( unsigned int cellDim , std::set< unsigned int > &nodeIndexSet ) const { nodeIndexSet.insert( _si[0] ) , nodeIndexSet.insert( _si[1] ); }
#endif // USE_UNORDERED_SET_MAP
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
#ifdef USE_UNORDERED_SET_MAP
	void setNodeIndexSet( unsigned int cellDim , std::unordered_set< unsigned int > &nodeIndexSet ) const
#else // !USE_UNORDERED_SET_MAP
	void setNodeIndexSet( unsigned int cellDim , std::set< unsigned int > &nodeIndexSet ) const
#endif // USE_UNORDERED_SET_MAP
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

		double weights[ COUNT ] , kWeights[ COUNT ];
		EnergyWeights( void ){ memset( weights , 0 , sizeof(weights) ) , memset( kWeights , 0 , sizeof(kWeights) ); }
		EnergyWeights( unsigned int eType ) : EnergyWeights() { weights[eType] = 1; }
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
	Eigen::MatrixXd systemMatrix( const	double *weights , std::function< bool ( FaceMultiIndex ) > isIntegrationFaceFunctor=[]( FaceMultiIndex ){ return true; } ) const;

	Eigen::MatrixXd crossFaceGradientDifferenceMatrix( std::function< bool ( FaceMultiIndex ) > faceSelectionFunctor ) const;

	struct NodeMultiIndex_Index
	{
		template< typename SubSimplexIndexFunctor >
		NodeMultiIndex_Index( SubSimplexIndexFunctor subSimplexIndexFunctor , unsigned int sz );
		unsigned int fromMultiIndex( NodeMultiIndex nmi ) const;
		const NodeMultiIndex &toMultiIndex( unsigned int idx ) const { return _nodeList[idx]; }
		unsigned int size( void ) const { return (unsigned int)_nodeMap.size(); }

		static NodeMultiIndex GetNodeMultiIndex( SimplexIndex< Dim , unsigned int > subSimplexIndex , unsigned int n );

	protected:
		typename NodeMultiIndex::map _nodeMap;
		std::vector< NodeMultiIndex > _nodeList;
	};

protected:
	const SimplexRefinable< Dim > &_simplexRefinable;
	NodeMultiIndex_Index _nodeIndex;

	template< typename SystemMatrixFunctor > Eigen::MatrixXd _systemMatrix( SystemMatrixFunctor F ) const;
};


// A struct for computing the prolongation matrix giving the expression of coarse basis function as a linear combination of finer basis function
// that minimizes the prescribed energy.
struct InterpolatingProlongationSystem
{
	template< unsigned int Dim , unsigned int EmbeddingDimension >
	static bool IsPlanar( const SimplexRefinableCell< Dim > &simplexRefinableCell , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon );

	// Given a SPD matrix E and a subset of indices, construct the prolongation matrix taking coarse basis functions
	// to finer basis functions, that:
	// -- Is interpolating (if index idx is in the subset, it will be mapped to itself)
	// -- Forms a partition of unity

	InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const Eigen::MatrixXd &kE , const std::vector< unsigned int > &coarseIndices , bool forcePoU );
	template< unsigned int InterpolationDim >
	InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const Eigen::MatrixXd &kE , const std::vector< unsigned int > &coarseIndices , bool forcePoU , const Point< double , InterpolationDim > *interpolationConstraints );

	// Trivial kernel
	InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , bool forcePoU ) : InterpolatingProlongationSystem( E , Eigen::MatrixXd::Zero( E.rows() , E.cols() ) , coarseIndices , forcePoU ) {}
	template< unsigned int InterpolationDim >
	InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , bool forcePoU , const Point< double , InterpolationDim > *interpolationConstraints ) : InterpolatingProlongationSystem( E , Eigen::MatrixXd::Zero( E.rows() , E.cols() ) , coarseIndices , forcePoU , interpolationConstraints ) {}

	template< unsigned int Degree >
	struct ProlongationInfo
	{
		Eigen::MatrixXd P;
		std::vector< typename SimplexRefinableElements< 0 , Degree >::NodeMultiIndex > coarseMultiIndices;
	};

	template< unsigned int Dim , unsigned int Degree >
	static void HierarchicalProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool forcePoU , ProlongationInfo< Degree > pInfo[Dim] , unsigned int finestDim=Dim );

	template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
	static void HierarchicalProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool forcePoU , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , ProlongationInfo< Degree > pInfo[Dim] , unsigned int finestDim );

	Eigen::MatrixXd prolongation( void ) const;
	double energy( const Eigen::MatrixXd &P ) const;
	double kernelEnergy( const Eigen::MatrixXd &P ) const;

protected:
	unsigned int _coarseDim , _fineDim;
	std::vector< unsigned int > _coarseIndices , _fineIndices;
	// The energy is defined as E(x) = Q(x,x) + < L , x > (not E(x) = Q(x,x) - < L , x >)
	Eigen::MatrixXd _Q , _kQ;
	Eigen::VectorXd _L , _kL;
	bool _hasKernelRegularizer , _hasKernel , _hasImage;
	Eigen::MatrixXd _kernel , _image;

	struct _DependentSystem
	{
		Eigen::VectorXd c;
		Eigen::SparseMatrix< double > C;
	};
	bool _isIndependent;
	_DependentSystem _dependentSystem;

	template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
	static std::vector< Point< double , EmbeddingDimension > > _InterpolationConstraints( const SimplexRefinableElements< Dim , Degree > &sre , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor );

	template< unsigned int InterpolationDim >
	void _init( const Eigen::MatrixXd &E , const Eigen::MatrixXd &kE , const std::vector< unsigned int > &coarseIndices , bool pou , const Point< double , InterpolationDim > *interpolationConstraints );

	unsigned int _index( unsigned fine , unsigned int coarse ) const;
	Eigen::MatrixXd _toMatrix( const Eigen::VectorXd &v ) const;
	Eigen::VectorXd _toVector( const Eigen::MatrixXd &P ) const;

	template< unsigned int Dim , unsigned int Degree >
	static void _SetNodeMaps( const SimplexRefinableCell< Dim > &simplexRefinableCell , typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex::map nodeMaps[Dim+1] );
	template< unsigned int Dim , unsigned int Degree >
	static void __SetNodeMaps( const SimplexRefinableCell< Dim > &simplexRefinableCell , typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex::map nodeMaps[Dim+1] );

	template< unsigned int Dim , unsigned int Degree >
	static void _HierarchicalProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		bool forcePoU ,
		ProlongationInfo< Degree > *pInfo , 
		unsigned int finestDim
	);

	template< unsigned int Dim , unsigned int Degree >
	static Eigen::MatrixXd _BoundaryProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		bool forcePoU ,
		const std::vector< unsigned int > &coarseIndices
	);

	template< unsigned int Dim , unsigned int Degree >
	static Eigen::MatrixXd _BoundaryProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		bool forcePoU ,
		const std::vector< unsigned int > &coarseIndices ,
		const Eigen::MatrixXd &coarseP
	);

	template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
	static void _HierarchicalProlongation
	(
		const SimplexRefinableCell< Dim > &simplexRefinableCell ,
		const SimplexRefinableElements< Dim , Degree > &sre ,
		typename SimplexRefinableElements<>::EnergyWeights eWeights ,
		bool forcePoU ,
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
		bool forcePoU ,
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
		bool forcePoU ,
		std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor ,
		double planarityEpsilon ,
		const std::vector< unsigned int > &coarseIndices ,
		const Eigen::MatrixXd &coarseP
	);
};

#include "SimplexRefinableBasis.inl"
#endif // SIMPLEX_REFINABLE_BASIS_INCLUDED