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

#ifndef SOLID_SIMPLEX_REFINABLE_MESH_INCLUDED
#define SOLID_SIMPLEX_REFINABLE_MESH_INCLUDED

#include "SolidSimplexMesh.h"
#include "SimplexRefinableMesh.h"

template< unsigned int Dim , unsigned int Degree >
struct HierarchicalSolidSimplexRefinableCellMesh
{
	typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	template< typename SimplexRefinableCellType >
	using CellList = typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinableCellType >;

	HierarchicalSolidSimplexRefinableCellMesh( void ){}

	HierarchicalSolidSimplexRefinableCellMesh( HierarchicalSolidSimplexRefinableCellMesh &&ssrcm ){ std::swap( _solidSimplexMesh , ssrcm._solidSimplexMesh ) , std::swap( _prolongationAndNodeMap , ssrcm._prolongationAndNodeMap ); }
	HierarchicalSolidSimplexRefinableCellMesh &operator = ( HierarchicalSolidSimplexRefinableCellMesh &&ssrcm ){ std::swap( _solidSimplexMesh , ssrcm._solidSimplexMesh ) , std::swap( _prolongationAndNodeMap , ssrcm._prolongationAndNodeMap ) ; return *this; }

	template< typename SimplexRefinableCellType >
	static HierarchicalSolidSimplexRefinableCellMesh Init( const CellList< SimplexRefinableCellType > &cellList , std::function< Point< double , Dim > ( unsigned int ) > vFunction , typename SimplexRefinableElements<>::EnergyWeights eWeights ,                          bool forcePoU , bool linearPrecision , double planarityEpsilon , bool verbose );

	template< typename SimplexRefinableCellType >
	static HierarchicalSolidSimplexRefinableCellMesh Init( const CellList< SimplexRefinableCellType > &cellList , std::function< Point< double , Dim > ( unsigned int ) > vFunction , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestDim , bool forcePoU , bool linearPrecision , double planarityEpsilon , bool verbose );

	unsigned int maxLevel( void ) const { return (unsigned int)_prolongationAndNodeMap.size(); };
	size_t nodes( unsigned int l ) const { return l<_prolongationAndNodeMap.size() ? _prolongationAndNodeMap[l].second.size() : _solidSimplexMesh.nodes(); }
	const std::map< NodeMultiIndex , unsigned int > &nodeMap( unsigned int l ) const{ return l<_prolongationAndNodeMap.size() ? _prolongationAndNodeMap[l].second : _solidSimplexMesh.nodeMap(); }
	const Eigen::SparseMatrix< double > &P( unsigned int l ) const { return _prolongationAndNodeMap[l].first; }
	unsigned int nodeIndex( unsigned int l , NodeMultiIndex ni ) const;

	Eigen::SparseMatrix< double > P( unsigned int lOut , unsigned int lIn ) const;

	size_t nodes( void ) const { return nodes( maxLevel()-1 ); }
	const std::map< NodeMultiIndex , unsigned int > &nodeMap( void ) const { return nodeMap( maxLevel()-1 ); }
	unsigned int nodeIndex( NodeMultiIndex ni ) const { return HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::nodeIndex( maxLevel()-1 , ni ); }

	const SolidSimplexMesh< Dim , Degree > &solidSimplexMesh( void ) const { return _solidSimplexMesh; }
protected:
	SolidSimplexMesh< Dim , Degree > _solidSimplexMesh;
	std::vector< std::pair< Eigen::SparseMatrix< double > , std::map< NodeMultiIndex , unsigned int > > > _prolongationAndNodeMap;

	template< typename SimplexRefinableCellType >
	void _init( const CellList< SimplexRefinableCellType > &cellList , std::function< Point< double , Dim > ( unsigned int ) > vFunction , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestNodeDim , bool forcePoU , bool linearPrecision , double planarityEpsilon , bool verbose );
};

template< unsigned int Dim , unsigned int Degree >
struct SolidSimplexRefinableCellMesh : protected HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >
{
	typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	template< typename SimplexRefinableCellType >
	using CellList = typename HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinableCellType >;

	SolidSimplexRefinableCellMesh( void ) : HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >() {}

	SolidSimplexRefinableCellMesh( SolidSimplexRefinableCellMesh &&ssrcm ) : HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >( std::move(ssrcm) ) { _level = ssrcm._level , std::swap( _P , ssrcm._P ) , std::swap( _Pt , ssrcm._Pt ); }
	SolidSimplexRefinableCellMesh &operator = ( SolidSimplexRefinableCellMesh &&ssrcm ){ HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::operator=( std::move(ssrcm) ) , _level = ssrcm._level , std::swap( _P , ssrcm._P ) , std::swap( _Pt , ssrcm._Pt ) ; return *this; }

	template< typename SimplexRefinableCellType >
	static SolidSimplexRefinableCellMesh Init( typename HierarchicalSimplexRefinableCellMesh< Dim , Degree >::template CellList< SimplexRefinableCellType > &cellList , std::function< Point< double , Dim > ( unsigned int ) > vFunction , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestNodeDim , bool forcePoU , bool linearPrecision , double planarityEpsilon , bool verbose )
	{
		SolidSimplexRefinableCellMesh ssrcm;
		ssrcm.template _init< SimplexRefinableCellType >( cellList , vFunction , eWeights , finestNodeDim , forcePoU , linearPrecision , planarityEpsilon , verbose ); 
		return ssrcm;
	}

	size_t nodes( void ) const { return HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::nodes( _level ); }
	const std::map< NodeMultiIndex , unsigned int > &nodeMap( void ) const { return HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::nodeMap( _level ); }
	unsigned int nodeIndex( NodeMultiIndex ni ) const { return HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::nodeIndex( _level , ni ); }

	const Eigen::SparseMatrix< double > &P ( void ) const { return _P ; }
	const Eigen::SparseMatrix< double > &Pt( void ) const { return _Pt; }

	Eigen::SparseMatrix< double > massMatrix( void ) const { return _Pt * solidSimplexMesh().massMatrix() * _P; }
	Eigen::SparseMatrix< double > frobeniusStiffnessMatrix( void ) const { return _Pt * solidSimplexMesh().frobeniusStiffnessMatrix() * _P; }
	Eigen::SparseMatrix< double > traceStiffnessMatrix( void ) const { return _Pt * solidSimplexMesh().traceStiffnessMatrix() * _P; }
	void setMassFrobeniusStiffnessAndTraceStiffnessMatrices( Eigen::SparseMatrix< double > &M , Eigen::SparseMatrix< double > &F , Eigen::SparseMatrix< double > &T ) const
	{
		Eigen::SparseMatrix< double > _M , _F , _T;
		solidSimplexMesh().setMassFrobeniusStiffnessAndTraceStiffnessMatrices( _M , _F , _T );
		M = _Pt * _M * _P;
		F = _Pt * _F * _P;
		T = _Pt * _T * _P;
	}
	Eigen::VectorXd dcVector( Point< double , Dim > v ) const { return _Pt * solidSimplexMesh().dcVector( v ); }
	Eigen::VectorXd stiffnessVector( void ) const { return _Pt * solidSimplexMesh().stiffnessVector(); }

	const SolidSimplexMesh< Dim , Degree > &solidSimplexMesh( void ) const { return HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::solidSimplexMesh(); }
protected:
	unsigned int _level;
	Eigen::SparseMatrix< double > _P , _Pt;

	template< typename SimplexRefinableCellType >
	void _init( const CellList< SimplexRefinableCellType > &cellList , std::function< Point< double , Dim > ( unsigned int ) > vFunction , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestNodeDim , bool forcePoU , bool linearPrecision , double planarityEpsilon , bool verbose )
	{
		HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::_init( cellList , vFunction , eWeights , finestNodeDim , forcePoU , linearPrecision , planarityEpsilon , verbose );
		_level = finestNodeDim;
		_P  = HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::P( _level+1 , _level );
		_Pt = HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::P( _level , _level+1 );
	}
};


#include "SolidSimplexRefinableMesh.inl"
#endif // SOLID_SIMPLEX_REFINABLE_MESH_INCLUDED