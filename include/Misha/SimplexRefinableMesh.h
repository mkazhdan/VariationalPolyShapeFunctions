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

#ifndef SIMPLEX_REFINABLE_MESH_INCLUDED
#define SIMPLEX_REFINABLE_MESH_INCLUDED
#include <iostream>
#include <map>
#include <vector>
#include <functional>
#include <type_traits>
#include <atomic>
#include <omp.h>
#include "Misha/Geometry.h"
#include "Misha/Timer.h"
#include "SimplexMesh.h"
#include "SimplexRefinableBasis.h"

template< unsigned int Dim , unsigned int Degree >
struct HierarchicalSimplexRefinableCellMesh
{
	typedef MultiIndex< Degree , unsigned int , false > NodeMultiIndex;
	typedef MultiIndex< Dim , unsigned int , false > FaceMultiIndex;

	HierarchicalSimplexRefinableCellMesh( void ){}

	// This struct represents a list of cells. 
	template< typename SimplexRefinableCellType >
	struct CellList
	{
		CellList( unsigned int cellNum , const std::function< SimplexRefinableCellType (unsigned int) > &cellFunctor ) : _cellNum(cellNum) , _cellFunctor(cellFunctor){}
		unsigned int size( void ) const { return _cellNum; }
		SimplexRefinableCellType operator[]( unsigned int idx ) const { return _cellFunctor(idx); }

		// This struct represents an object that can return the list of boundary cells
		struct BoundaryCellList
		{
			typedef typename SimplexRefinableCellType::SimplexRefinableBoundaryCellType BoundaryCellType;

			BoundaryCellList( const CellList &cellList ) : _cellList( cellList )
			{
				std::unordered_map< unsigned int , std::pair< unsigned int , unsigned int > > boundaryMap;
				for( unsigned int c=0 ; c<_cellList.size() ; c++ )
				{
					SimplexRefinableCellType src = _cellList[c];
					for( unsigned int f=0 ; f<src.faces() ; f++ ) boundaryMap[ src.face[f].centerIndex() ] = std::make_pair( c , f );
				}
				_boundaryFaces.reserve( boundaryMap.size() );
				for( const auto & [ idx , f ] : boundaryMap ) _boundaryFaces.push_back(f);

				_boundaryCellFunctor = [&]( unsigned int idx ){ return _cellList[ _boundaryFaces[idx].first ].face[ _boundaryFaces[idx].second ]; };
				_boundaryCellList = new typename HierarchicalSimplexRefinableCellMesh< Dim-1 , Degree >::template CellList< BoundaryCellType >( _boundaryFaces.size() , _boundaryCellFunctor );
			}

			~BoundaryCellList( void ) { delete _boundaryCellList; }

			const typename HierarchicalSimplexRefinableCellMesh< Dim-1 , Degree >::template CellList< BoundaryCellType > &operator()( void ) const
			{
				return *_boundaryCellList;
			}
		protected:
			const CellList &_cellList;
			std::function< BoundaryCellType ( unsigned int ) > _boundaryCellFunctor;
			std::vector< std::pair< unsigned int , unsigned int > > _boundaryFaces;
			typename HierarchicalSimplexRefinableCellMesh< Dim-1 , Degree >::template CellList< BoundaryCellType > *_boundaryCellList;
		};

	protected:
		unsigned int _cellNum;
		const std::function< SimplexRefinableCellType (unsigned int) > &_cellFunctor;
	};

	template< typename SimplexRefinableCellType >
	static HierarchicalSimplexRefinableCellMesh Init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou ,                          bool verbose=false );

	template< typename SimplexRefinableCellType >
	static HierarchicalSimplexRefinableCellMesh Init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , unsigned int finestDim , bool verbose=false );

	template< typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
	static HierarchicalSimplexRefinableCellMesh Init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon ,                          bool verbose=false );

	template< typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
	static HierarchicalSimplexRefinableCellMesh Init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestDim , bool verbose=false );

	HierarchicalSimplexRefinableCellMesh( HierarchicalSimplexRefinableCellMesh &&srcm ){ std::swap( _simplexMesh , srcm._simplexMesh ) , std::swap( _prolongationAndNodeMap , srcm._prolongationAndNodeMap ); }
	HierarchicalSimplexRefinableCellMesh &operator = ( HierarchicalSimplexRefinableCellMesh &&srcm ){ std::swap( _simplexMesh , srcm._simplexMesh ) , std::swap( _prolongationAndNodeMap , srcm._prolongationAndNodeMap ) ; return *this; }

	unsigned int maxLevel( void ) const { return (unsigned int)_prolongationAndNodeMap.size(); };
	size_t nodes( unsigned int l ) const { return l<_prolongationAndNodeMap.size() ? _prolongationAndNodeMap[l].second.size() : _simplexMesh.nodes(); }
	const std::map< NodeMultiIndex , unsigned int > &nodeMap( unsigned int l ) const { return l<_prolongationAndNodeMap.size() ? _prolongationAndNodeMap[l].second : _simplexMesh.nodeMap(); }
	const Eigen::SparseMatrix< double > &P( unsigned int l ) const { return _prolongationAndNodeMap[l].first; }
	unsigned int nodeIndex( unsigned int l , NodeMultiIndex ni ) const;

	Eigen::SparseMatrix< double > P( unsigned int lOut , unsigned int lIn ) const;

	size_t nodes( void ) const { return nodes( maxLevel()-1 ); }
	const std::map< NodeMultiIndex , unsigned int > &nodeMap( void ) const { return nodeMap( maxLevel()-1 ); }
	unsigned int nodeIndex( NodeMultiIndex ni ) const { return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::nodeIndex( maxLevel()-1 , ni ); }

	const SimplexMesh< Dim , Degree > &simplexMesh( void ) const { return _simplexMesh; }
	void hashSimplexMeshLocalToGlobalNodeIndex( void ){ _simplexMesh.hashLocalToGlobalNodeIndex(); }


	double volume( void ) const { return _simplexMesh.volume(); }
	void makeUnitVolume( void ) { _simplexMesh.makeUnitVolume(); }
	void setMetric( std::function< SquareMatrix< double , Dim > (unsigned int) > metricFunction ){ _simplexMesh.setMetric( metricFunction ); }
protected:
	SimplexMesh< Dim , Degree > _simplexMesh;
	std::vector< std::pair< Eigen::SparseMatrix< double > , std::map< NodeMultiIndex , unsigned int > > > _prolongationAndNodeMap;

	template< typename SimplexRefinableCellType >
	void _init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , unsigned int finestDim , bool verbose );

	template< typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
	void _init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestDim , bool verbose );

	template< typename SimplexRefinableCellType >
	void _setSimplexMesh( const CellList< SimplexRefinableCellType > &cellList , bool verbose );

	template< typename SimplexRefinableCellType >
	void _setProlongationAndNodeMap( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , unsigned int finestDim , bool verbose );

	template< typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
	void _setProlongationAndNodeMap( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestDim , bool verbose );
};

template< unsigned int Dim , unsigned int Degree >
struct SimplexRefinableCellMesh : protected HierarchicalSimplexRefinableCellMesh< Dim , Degree >
{
	typedef MultiIndex< Degree , unsigned int , false > NodeMultiIndex;
	template< typename SimplexRefinableCellType >
	using CellList = typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinableCellType >;

	SimplexRefinableCellMesh( void ) : HierarchicalSimplexRefinableCellMesh< Dim , Degree >(){}

	SimplexRefinableCellMesh( SimplexRefinableCellMesh &&srcm ) : HierarchicalSimplexRefinableCellMesh< Dim , Degree >( std::move(srcm) ) { _level = srcm._level , std::swap( _P , srcm._P ) , std::swap( _Pt , srcm._Pt ); }
	SimplexRefinableCellMesh &operator = ( SimplexRefinableCellMesh &&srcm ){ HierarchicalSimplexRefinableCellMesh< Dim , Degree >::operator=( std::move(srcm) ) , _level = srcm._level , std::swap( _P , srcm._P ) , std::swap( _Pt , srcm._Pt ) ; return *this; }

	template< typename SimplexRefinableCellType >
	static SimplexRefinableCellMesh Init( const typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , unsigned int finestNodeDim , bool verbose=false )
	{
		SimplexRefinableCellMesh srcm;
		srcm.template _init< SimplexRefinableCellType >( cellList , eWeights , pou , finestNodeDim , verbose );
		return srcm;
	}

	template< typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
	static SimplexRefinableCellMesh Init( const typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , std::function< Point< double , EmbeddingDimension > (unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestNodeDim , bool verbose )
	{
		SimplexRefinableCellMesh srcm;
		srcm.template _init< SimplexRefinableCellType , EmbeddingDimension >( cellList , eWeights , pou , positionFunctor , planarityEpsilon , finestNodeDim , verbose );
		return srcm;
	}

	size_t nodes( void ) const { return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::nodes( _level ); }
	const std::map< NodeMultiIndex , unsigned int > &nodeMap( void ) const { return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::nodeMap( _level ); }
	const Eigen::SparseMatrix< double > &P ( void ) const { return _P ; }
	const Eigen::SparseMatrix< double > &Pt( void ) const { return _Pt; }
	unsigned int nodeIndex( NodeMultiIndex ni ) const { return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::nodeIndex( _level , ni ); }

	const SimplexMesh< Dim , Degree > &simplexMesh( void ) const { return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::simplexMesh(); }
	void hashSimplexMeshLocalToGlobalNodeIndex( void ){ HierarchicalSimplexRefinableCellMesh< Dim , Degree >::hashSimplexMeshLocalToGlobalNodeIndex(); }

	double volume( void ) const { return HierarchicalSimplexRefinableCellMesh< Dim , Degree >::volume(); }
	void makeUnitVolume( void ) { HierarchicalSimplexRefinableCellMesh< Dim , Degree >::makeUnitVolume(); }
	void setMetric( std::function< SquareMatrix< double , Dim > (unsigned int) > metricFunction ){ HierarchicalSimplexRefinableCellMesh< Dim , Degree >::setMetric( metricFunction ); }
protected:
	unsigned int _level;
	Eigen::SparseMatrix< double > _P , _Pt;

	template< typename SimplexRefinableCellType >
	void _init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , unsigned int finestNodeDim , bool verbose )
	{
		HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_init( cellList , eWeights , pou , finestNodeDim , verbose );
		_level = finestNodeDim;
		_P  = HierarchicalSimplexRefinableCellMesh< Dim , Degree >::P( _level+1 , _level );
		_Pt = HierarchicalSimplexRefinableCellMesh< Dim , Degree >::P( _level , _level+1 );
	}

	template< typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
	void _init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool pou , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestNodeDim , bool verbose )
	{
		HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_init( cellList , eWeights , pou , positionFunctor , planarityEpsilon , finestNodeDim , verbose );
		_level = finestNodeDim;
		_P  = HierarchicalSimplexRefinableCellMesh< Dim , Degree >::P( _level+1 , _level );
		_Pt = HierarchicalSimplexRefinableCellMesh< Dim , Degree >::P( _level , _level+1 );
	}
};

#include "SimplexRefinableMesh.inl"
#endif // SIMPLEX_REFINABLE_MESH_INCLUDED