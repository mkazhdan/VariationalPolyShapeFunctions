//////////////////////////////////////////
// HierarchicalSimplexRefinableCellMesh //
//////////////////////////////////////////
template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType >
HierarchicalSimplexRefinableCellMesh< Dim , Degree > HierarchicalSimplexRefinableCellMesh< Dim , Degree >::Init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , bool verbose )
{
	return Init< PoU >( cellList , eWeights , Dim , verbose );
}

template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType >
HierarchicalSimplexRefinableCellMesh< Dim , Degree > HierarchicalSimplexRefinableCellMesh< Dim , Degree >::Init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestDim , bool verbose )
{
	static_assert( std::is_base_of< SimplexRefinableCell< Dim > , SimplexRefinableCellType >::value , "[ERROR] SimplexRefinableCellType must derive from SimplexRefinableCell" );
	HierarchicalSimplexRefinableCellMesh srcm;
	srcm.template _init< PoU >( cellList , eWeights , finestDim , verbose );
	return srcm;
}

template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
HierarchicalSimplexRefinableCellMesh< Dim , Degree > HierarchicalSimplexRefinableCellMesh< Dim , Degree >::Init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , bool verbose )
{
	return Init< PoU >( cellList , eWeights , positionFunctor , planarityEpsilon , Dim , verbose );
}

template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
HierarchicalSimplexRefinableCellMesh< Dim , Degree > HierarchicalSimplexRefinableCellMesh< Dim , Degree >::Init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestDim , bool verbose )
{
	static_assert( std::is_base_of< SimplexRefinableCell< Dim > , SimplexRefinableCellType >::value , "[ERROR] SimplexRefinableCellType must derive from SimplexRefinableCell" );
	HierarchicalSimplexRefinableCellMesh srcm;
	srcm.template _init< PoU >( cellList , eWeights , positionFunctor , planarityEpsilon , finestDim , verbose );
	return srcm;
}

template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType >
void HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestDim , bool verbose )
{
	_setSimplexMesh           <       SimplexRefinableCellType >( cellList ,                        verbose );
	_setProlongationAndNodeMap< PoU , SimplexRefinableCellType >( cellList , eWeights , finestDim , verbose );
}

template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
void HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestDim , bool verbose )
{
	_setSimplexMesh           <       SimplexRefinableCellType >( cellList ,                                                             verbose );
	_setProlongationAndNodeMap< PoU , SimplexRefinableCellType >( cellList , eWeights , positionFunctor , planarityEpsilon , finestDim , verbose );
}

template< unsigned int Dim , unsigned int Degree >
template< typename SimplexRefinableCellType >
void HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_setSimplexMesh( const CellList< SimplexRefinableCellType > &cellList , bool verbose )
{
	// The total list of simplices
	std::vector< SimplexIndex< Dim , unsigned int > > simplices;
	// The total list of metrics
	std::vector< SquareMatrix< double , Dim > > metrics;

	Timer timer;
	unsigned int sCount = 0;
	for( unsigned int c=0 ; c<cellList.size() ; c++ ) sCount += cellList[c].size();
	simplices.resize( sCount );
	metrics.resize( sCount );

	{
		std::vector< std::vector< std::pair< SimplexIndex< Dim , unsigned int > , SquareMatrix< double , Dim > > > > _simplicesAndMetrics( omp_get_max_threads() );
		for( unsigned int t=0 ; t<_simplicesAndMetrics.size() ; t++ ) _simplicesAndMetrics[t].reserve( sCount );

#pragma omp parallel for
		for( int c=0 ; c<(int)cellList.size() ; c++ )
		{
			SimplexRefinableCellType simplexRefinable = cellList[c];
			std::vector< std::pair< SimplexIndex< Dim , unsigned int > , SquareMatrix< double , Dim > > > &simplicesAndMetrics = _simplicesAndMetrics[ omp_get_thread_num() ];
			for( unsigned int s=0 ; s<simplexRefinable.size() ; s++ ) simplicesAndMetrics.push_back( std::make_pair( simplexRefinable[s] , simplexRefinable.metric(s) ) );
		}
		unsigned int idx = 0;
		for( unsigned int t=0 ; t<_simplicesAndMetrics.size() ; t++ ) for( unsigned int i=0 ; i<_simplicesAndMetrics[t].size() ; i++ )
		{
			simplices[idx] = _simplicesAndMetrics[t][i].first , metrics[idx] = _simplicesAndMetrics[t][i].second;
			idx++;
		}
	}
	std::function< SquareMatrix< double , Dim > ( unsigned int ) > gFunction = [&]( unsigned int idx ){ return metrics[idx]; };
	_simplexMesh = SimplexMesh< Dim , Degree >::Init( simplices , gFunction );

	if( verbose ) std::cout << "Got simplicial mesh: " << timer.elapsed() << std::endl;
}

template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType >
void HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_setProlongationAndNodeMap( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestDim , bool verbose )
{
	_prolongationAndNodeMap.resize( finestDim+1 );

	struct PEntry
	{
		NodeMultiIndex coarse , fine;
		double value;
		PEntry( void ) : value(0) {}
		PEntry( NodeMultiIndex c , NodeMultiIndex f , double v ) : coarse(c) , fine(f) , value(v) {}
	};

	// The prolongation entries (represented using MultiNodeIndex)
	std::vector< PEntry > pEntries[Dim];

	Timer timer;
	{
		std::vector< std::vector< PEntry > > _pEntries[Dim];
		std::vector< std::map< NodeMultiIndex , unsigned int > > _nodeMaps[Dim];
		for( unsigned int d=0 ; d<Dim ; d++ ) _pEntries[d].resize( omp_get_max_threads() ) , _nodeMaps[d].resize( omp_get_max_threads() );
#pragma omp parallel for
		for( int c=0 ; c<(int)cellList.size() ; c++ )
		{
			unsigned int tIdx = omp_get_thread_num();
			SimplexRefinableCellType simplexRefinable = cellList[c];
			SimplexRefinableElements< Dim , Degree > sre( simplexRefinable );
			typename InterpolatingProlongationSystem< PoU >::template ProlongationInfo< Degree > pInfo[Dim];
			InterpolatingProlongationSystem< PoU >::template HierarchicalProlongation< Dim , Degree >( simplexRefinable , eWeights , pInfo , finestDim );
			for( unsigned int d=0 ; d<Dim && d<=finestDim ; d++ )
			{
				const std::vector< NodeMultiIndex > &coarseMultiIndices = pInfo[d].coarseMultiIndices;
				const Eigen::MatrixXd &P = pInfo[d].P;
				for( unsigned int i=0 ; i<coarseMultiIndices.size() ; i++ ) _nodeMaps[d][tIdx][ coarseMultiIndices[i] ] = 0;
				for( unsigned int r=0 ; r<P.rows() ; r++ ) for( unsigned int c=0 ; c<P.cols() ; c++ ) if( P(r,c) ) _pEntries[d][tIdx].push_back( PEntry( coarseMultiIndices[c] , sre[r] , P(r,c) ) );
			}
		}

		for( unsigned int d=0 ; d<Dim && d<=finestDim ; d++ )
		{
			unsigned int count=0;
			for( unsigned int t=0 ; t<_pEntries[d].size() ; t++ ) count += (unsigned int)_pEntries[d][t].size();
			pEntries[d].reserve( count );
			for( unsigned int t=0 ; t<_pEntries[d].size() ; t++ )
			{
				for( unsigned int i=0 ; i<_pEntries[d][t].size() ; i++ ) pEntries[d].push_back( _pEntries[d][t][i] );
				for( const auto &[nmi,idx] : _nodeMaps[d][t] ) _prolongationAndNodeMap[d].second[ nmi ] = 0;
			}
		}
	}

	for( unsigned int d=0 ; d<Dim && d<=finestDim ; d++ )
	{
		unsigned int nodeCount = 0;
		for( auto &iter : _prolongationAndNodeMap[d].second ) iter.second = nodeCount++;
	}

	if( finestDim==Dim )
	{
		_prolongationAndNodeMap[Dim].second = _simplexMesh.nodeMap();
		_prolongationAndNodeMap[Dim].first.resize( nodes(Dim) , nodes(Dim) );
		_prolongationAndNodeMap[Dim].first.setIdentity();
	}

	for( unsigned int d=0 ; d<Dim && d<=finestDim ; d++ )
	{
		std::vector< double > rowSums( _simplexMesh.nodes() , 0. );
		for( unsigned int i=0 ; i<pEntries[d].size() ; i++ ) rowSums[ _simplexMesh.nodeIndex( pEntries[d][i].fine ) ] += pEntries[d][i].value;
		for( unsigned int i=0 ; i<pEntries[d].size() ; i++ ) pEntries[d][i].value /= rowSums[ _simplexMesh.nodeIndex( pEntries[d][i].fine ) ];

		std::vector< Eigen::Triplet< double > > entries( pEntries[d].size() );
		for( unsigned int i=0 ; i<pEntries[d].size() ; i++ ) entries[i] = Eigen::Triplet< double >( nodeIndex( d+1 , pEntries[d][i].fine ) , nodeIndex( d , pEntries[d][i].coarse ) , pEntries[d][i].value );

		_prolongationAndNodeMap[d].first.resize( nodes(d+1) , nodes(d) );
		_prolongationAndNodeMap[d].first.setFromTriplets( entries.begin() , entries.end() );		
	}

	if( verbose ) std::cout << "Got prolongation/restriction: " << timer.elapsed() << std::endl;
}

template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
void HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_setProlongationAndNodeMap( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestDim , bool verbose )
{
	_prolongationAndNodeMap.resize( finestDim+1 );

	struct PEntry
	{
		NodeMultiIndex coarse , fine;
		double value;
		PEntry( void ) : value(0) {}
		PEntry( NodeMultiIndex c , NodeMultiIndex f , double v ) : coarse(c) , fine(f) , value(v) {}
	};

	// The prolongation entries (represented using MultiNodeIndex)
	std::vector< PEntry > pEntries[Dim];

	Timer timer;
	{
		std::vector< std::vector< PEntry > > _pEntries[Dim];
		std::vector< std::map< NodeMultiIndex , unsigned int > > _nodeMaps[Dim];
		for( unsigned int d=0 ; d<Dim ; d++ ) _pEntries[d].resize( omp_get_max_threads() ) , _nodeMaps[d].resize( omp_get_max_threads() );
#pragma omp parallel for
		for( int c=0 ; c<(int)cellList.size() ; c++ )
		{
			unsigned int tIdx = omp_get_thread_num();
			SimplexRefinableCellType simplexRefinable = cellList[c];
			SimplexRefinableElements< Dim , Degree > sre( simplexRefinable );
			typename InterpolatingProlongationSystem< PoU >::template ProlongationInfo< Degree > pInfo[Dim];
			InterpolatingProlongationSystem< PoU >::template HierarchicalProlongation< Dim , Degree >( simplexRefinable , eWeights , positionFunctor, planarityEpsilon , pInfo , finestDim );
			for( unsigned int d=0 ; d<Dim && d<=finestDim ; d++ )
			{
				const std::vector< NodeMultiIndex > &coarseMultiIndices = pInfo[d].coarseMultiIndices;
				const Eigen::MatrixXd &P = pInfo[d].P;
				for( unsigned int i=0 ; i<coarseMultiIndices.size() ; i++ ) _nodeMaps[d][tIdx][ coarseMultiIndices[i] ] = 0;
				for( unsigned int r=0 ; r<P.rows() ; r++ ) for( unsigned int c=0 ; c<P.cols() ; c++ ) if( P(r,c) ) _pEntries[d][tIdx].push_back( PEntry( coarseMultiIndices[c] , sre[r] , P(r,c) ) );
			}
		}

		for( unsigned int d=0 ; d<Dim && d<=finestDim ; d++ )
		{
			unsigned int count=0;
			for( unsigned int t=0 ; t<_pEntries[d].size() ; t++ ) count += (unsigned int)_pEntries[d][t].size();
			pEntries[d].reserve( count );
			for( unsigned int t=0 ; t<_pEntries[d].size() ; t++ )
			{
				for( unsigned int i=0 ; i<_pEntries[d][t].size() ; i++ ) pEntries[d].push_back( _pEntries[d][t][i] );
				for( const auto &[nmi,idx] : _nodeMaps[d][t] ) _prolongationAndNodeMap[d].second[ nmi ] = 0;
			}
		}
	}

	for( unsigned int d=0 ; d<Dim && d<=finestDim ; d++ )
	{
		unsigned int nodeCount = 0;
		for( auto &iter : _prolongationAndNodeMap[d].second ) iter.second = nodeCount++;
	}

	if( finestDim==Dim )
	{
		_prolongationAndNodeMap[Dim].second = _simplexMesh.nodeMap();
		_prolongationAndNodeMap[Dim].first.resize( nodes(Dim) , nodes(Dim) );
		_prolongationAndNodeMap[Dim].first.setIdentity();
	}

	for( unsigned int d=0 ; d<Dim && d<=finestDim ; d++ )
	{
		std::vector< double > rowSums( _simplexMesh.nodes() , 0. );
		for( unsigned int i=0 ; i<pEntries[d].size() ; i++ ) rowSums[ _simplexMesh.nodeIndex( pEntries[d][i].fine ) ] += pEntries[d][i].value;
		for( unsigned int i=0 ; i<pEntries[d].size() ; i++ ) pEntries[d][i].value /= rowSums[ _simplexMesh.nodeIndex( pEntries[d][i].fine ) ];

		std::vector< Eigen::Triplet< double > > entries( pEntries[d].size() );
		for( unsigned int i=0 ; i<pEntries[d].size() ; i++ ) entries[i] = Eigen::Triplet< double >( nodeIndex( d+1 , pEntries[d][i].fine ) , nodeIndex( d , pEntries[d][i].coarse ) , pEntries[d][i].value );

		_prolongationAndNodeMap[d].first.resize( nodes(d+1) , nodes(d) );
		_prolongationAndNodeMap[d].first.setFromTriplets( entries.begin() , entries.end() );		
	}

	if( verbose ) std::cout << "Got prolongation/restriction: " << timer.elapsed() << std::endl;
}

template< unsigned int Dim , unsigned int Degree >
unsigned int HierarchicalSimplexRefinableCellMesh< Dim , Degree >::nodeIndex( unsigned int l , NodeMultiIndex ni ) const
{
	if( l<_prolongationAndNodeMap.size() )
	{
		auto iter = _prolongationAndNodeMap[l].second.find( ni );
		if( iter==_prolongationAndNodeMap[l].second.end() ) ERROR_OUT( "Could not find node index: " , ni );
		return iter->second;
	}
	else return _simplexMesh.nodeIndex(ni);
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > HierarchicalSimplexRefinableCellMesh< Dim , Degree >::P( unsigned int lOut , unsigned int lIn ) const
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
