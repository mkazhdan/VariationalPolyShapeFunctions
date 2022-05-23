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

///////////////////////////////////////////////
// HierarchicalSolidSimplexRefinableCellMesh //
///////////////////////////////////////////////
template< unsigned int Dim , unsigned int Degree >
template< typename SimplexRefinableCellType >
HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree > HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::Init( const CellList< SimplexRefinableCellType > &cellList , std::function< Point< double , Dim > ( unsigned int ) > vFunction , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool forcePoU , bool linearPrecision , double planarityEpsilon , bool verbose )
{
	return Init( cellList , vFunction , eWeights , Dim , forcePoU , linearPrecision , planarityEpsilon , verbose );
}

template< unsigned int Dim , unsigned int Degree >
template< typename SimplexRefinableCellType >
HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree > HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::Init( const CellList< SimplexRefinableCellType > &cellList , std::function< Point< double , Dim > ( unsigned int ) > vFunction , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestDim , bool forcePoU , bool linearPrecision , double planarityEpsilon , bool verbose )
{
	static_assert( std::is_base_of< SimplexRefinable< Dim > , SimplexRefinableCellType >::value , "[ERROR] SimplexRefinableCellType must derive from SimplexRefinable" );
	HierarchicalSolidSimplexRefinableCellMesh ssrcm;
	ssrcm._init( cellList , vFunction , eWeights , finestDim , forcePoU , linearPrecision , planarityEpsilon , verbose );
	return ssrcm;
}

template< unsigned int Dim , unsigned int Degree >
template< typename SimplexRefinableCellType >
void HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::_init( const CellList< SimplexRefinableCellType > &cellList , std::function< Point< double , Dim > ( unsigned int ) > vFunction , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestDim , bool forcePoU , bool linearPrecision , double planarityEpsilon , bool verbose )
{
	// Start by constructing the underlying SolidSimplexMesh
	{
		Miscellany::NestedTimer timer( "Got solid simplex mesh" , verbose );

		std::vector< SimplexIndex< Dim , unsigned int > > simplices;

		unsigned int sCount = 0;
		for( unsigned int c=0 ; c<cellList.size() ; c++ ) sCount += cellList[c].size();

		// Get the simplices
		simplices.reserve( sCount );
		for( unsigned int c=0 ; c<cellList.size() ; c++ )
		{
			SimplexRefinableCellType simplexRefinable = cellList[c];
			for( unsigned int s=0 ; s<simplexRefinable.size() ; s++ ) simplices.push_back( simplexRefinable[s] );
		}

		_solidSimplexMesh = SolidSimplexMesh< Dim , Degree >::Init( simplices , vFunction );
	}

	// Next, compute the prolongation information
	{
		HierarchicalSimplexRefinableCellMesh< Dim , Degree > srm;
		{
			Miscellany::NestedTimer( "Got simplex refinable cell mesh" , verbose );
			if( linearPrecision ) srm = HierarchicalSimplexRefinableCellMesh< Dim , Degree >::Init( cellList , eWeights , forcePoU , vFunction , planarityEpsilon , finestDim , verbose );
			else srm = HierarchicalSimplexRefinableCellMesh< Dim , Degree >::Init( cellList , eWeights , forcePoU , finestDim , verbose );
		}
		{
			Miscellany::NestedTimer( "Expanded solid prolongation and node map" , verbose );
			_prolongationAndNodeMap.resize( srm.maxLevel() );
			for( unsigned int l=0 ; l<srm.maxLevel() ; l++ )
			{
				_prolongationAndNodeMap[l].first = BlockExpand< Dim >( srm.P(l+1,l) );
				_prolongationAndNodeMap[l].second = srm.nodeMap(l);
			}
		}
	}
}


template< unsigned int Dim , unsigned int Degree >
unsigned int HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::nodeIndex( unsigned int l , NodeMultiIndex ni ) const
{
	if( l<_prolongationAndNodeMap.size() )
	{
		auto iter = _prolongationAndNodeMap[l].second.find( ni );
		if( iter==_prolongationAndNodeMap[l].second.end() ) ERROR_OUT( "Could not find node index: " , ni );
		return iter->second;
	}
	else return _solidSimplexMesh.nodeIndex( ni );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > HierarchicalSolidSimplexRefinableCellMesh< Dim , Degree >::P( unsigned int lOut , unsigned int lIn ) const
{
	if( lOut>_prolongationAndNodeMap.size() ) ERROR_OUT( "Output prolongation level exceeded: " , lOut , " >= " , _prolongationAndNodeMap.size() );
	if( lIn >_prolongationAndNodeMap.size() ) ERROR_OUT(  "Input prolongation level exceeded: " , lIn  , " >= " , _prolongationAndNodeMap.size() );
	if( lOut<lIn ) return P( lIn , lOut ).transpose();
	else if( lOut==lIn )
	{
		Eigen::SparseMatrix< double > I( nodes(lOut) , nodes(lIn) );
		I.setIdentity();
		return I;
	}
	else
	{
		Eigen::SparseMatrix< double > P = _prolongationAndNodeMap[lIn].first;
		for( unsigned int l=lIn+1 ; l<lOut ; l++ ) P = _prolongationAndNodeMap[l].first * P;
		return P;
	}
}
