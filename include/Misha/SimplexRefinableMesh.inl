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

#ifdef INTERPOLATION_CONSTRAINTS
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
#endif // INTERPOLATION_CONSTRAINTS

template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType >
void HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , unsigned int finestDim , bool verbose )
{
	_setSimplexMesh           <       SimplexRefinableCellType >( cellList ,                        verbose );
	_setProlongationAndNodeMap< PoU , SimplexRefinableCellType >( cellList , eWeights , finestDim , verbose );
}

#ifdef INTERPOLATION_CONSTRAINTS
template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
void HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_init( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestDim , bool verbose )
{
	_setSimplexMesh           <       SimplexRefinableCellType >( cellList ,                                                             verbose );
	_setProlongationAndNodeMap< PoU , SimplexRefinableCellType >( cellList , eWeights , positionFunctor , planarityEpsilon , finestDim , verbose );
}
#endif // INTERPOLATION_CONSTRAINTS

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

#ifdef NEW_SIMPLEX_MESH
	_cellSimplices.resize( cellNum );
#endif // NEW_SIMPLEX_MESH

	{
#ifdef NEW_SIMPLEX_MESH
		struct SimplexData
		{
			unsigned int cellIndex;
			SimplexIndex< Dim , unsigned int > simplexIndex;
			SquareMatrix< double , Dim > simplexMetric;

			SimplexData( void ){}
			SimplexData( unsigned int ci , SimplexIndex< Dim , unsigned int > si , SquareMatrix< double , Dim > sm )
				: cellIndex(ci) , simplexIndex(si) , simplexMetric(sm) {}
		};
		std::vector< std::vector< SimplexData > > _simplexData( omp_get_max_threads() );
		for( unsigned int t=0 ; t<_simplexData.size() ; t++ ) _simplexData.reserve( sCount );
#else // !NEW_SIMPLEX_MESH
		std::vector< std::vector< std::pair< SimplexIndex< Dim , unsigned int > , SquareMatrix< double , Dim > > > > _simplicesAndMetrics( omp_get_max_threads() );
		for( unsigned int t=0 ; t<_simplicesAndMetrics.size() ; t++ ) _simplicesAndMetrics[t].reserve( sCount );
#endif // NEW_SIMPLEX_MESH

#pragma omp parallel for
		for( int c=0 ; c<(int)cellList.size() ; c++ )
		{
			SimplexRefinableCellType simplexRefinable = cellList[c];
#ifdef NEW_SIMPLEX_MESH
			_cellSimplices[c].reserve( simplexRefinable.size() );
			std::vector< SimplexData > &simplexData = _simplexData[ omp_get_thread_num() ];
			for( unsigned int s=0 ; s<simplexRefinable.size() ; s++ ) simplexData.push_back( SimplexData( c , simplexRefinable[s] , simplexRefinable.metric(s) ) );
#else // !NEW_SIMPLEX_MESH
			std::vector< std::pair< SimplexIndex< Dim , unsigned int > , SquareMatrix< double , Dim > > > &simplicesAndMetrics = _simplicesAndMetrics[ omp_get_thread_num() ];
			for( unsigned int s=0 ; s<simplexRefinable.size() ; s++ ) simplicesAndMetrics.push_back( std::make_pair( simplexRefinable[s] , simplexRefinable.metric(s) ) );
#endif // NEW_SIMPLEX_MESH
		}
#ifdef NEW_SIMPLEX_MESH
		unsigned int idx = 0;
		for( unsigned int t=0 ; t<_simplexData.size() ; t++ ) for( unsigned int i=0 ; i<_simplexData[t].size() ; i++ )
		{
			simplices[idx] = _simplexData[t][i].simplexIndex;
			metrics[idx] = _simplexData[t][i].simplexMetric;
			_cellSimplices[ _simplexData[t][i].cellIndex ].push_back( idx );
			idx++;
		}
#else // !NEW_SIMPLEX_MESH
		unsigned int idx = 0;
		for( unsigned int t=0 ; t<_simplicesAndMetrics.size() ; t++ ) for( unsigned int i=0 ; i<_simplicesAndMetrics[t].size() ; i++ )
		{
			simplices[idx] = _simplicesAndMetrics[t][i].first , metrics[idx] = _simplicesAndMetrics[t][i].second;
			idx++;
		}
#endif // NEW_SIMPLEX_MESH
	}
	std::function< SquareMatrix< double , Dim > ( unsigned int ) > gFunction = [&]( unsigned int idx ){ return metrics[idx]; };
	_simplexMesh = SimplexMesh< Dim , Degree >::Init( simplices , gFunction );

#ifdef NEW_SIMPLEX_MESH
	_vertexDim.resize( _simplexMesh.vertices() , 0 );
	if constexpr( Dim>1 )
		for( int c=0 ; c<(int)cellNum ; c++ )
		{
			SimplexRefinableCellType simplexRefinable = cellFunctor( c );
			_setVertexDim( simplexRefinable );
		}
#endif // NEW_SIMPLEX_MESH

	if( verbose ) std::cout << "Got simplicial mesh: " << timer.elapsed() << std::endl;
}

#ifdef NEW_SIMPLEX_MESH
template< unsigned int Dim , unsigned int Degree >
template< unsigned int _Dim >
void HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_setVertexDim( const SimplexRefinableCell< _Dim > &src )
{
	_vertexDim[ src.centerIndex() ] = _Dim;
	if constexpr( _Dim-1>1 ) for( unsigned int f=0 ; f<src.faces() ; f++ ) _setVertexDim( src.face(f) );
}

template< unsigned int Dim , unsigned int Degree >
unsigned int HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_nodeDim( const NodeMultiIndex &multiIndex ) const
{
	unsigned int dim = SimplexMesh< Dim , Degree >::NodeDim( multiIndex );
	for( unsigned int d=0 ; d<Degree ; d++ ) dim = std::max< unsigned int >( dim , _vertexDim[ multiIndex[d] ] );
	return dim;
}
#endif // NEW_SIMPLEX_MESH

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

#ifdef INTERPOLATION_CONSTRAINTS
template< unsigned int Dim , unsigned int Degree >
template< bool PoU , typename SimplexRefinableCellType , unsigned int EmbeddingDimension >
void HierarchicalSimplexRefinableCellMesh< Dim , Degree >::_setProlongationAndNodeMap( const CellList< SimplexRefinableCellType > &cellList , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , unsigned int finestDim , bool verbose )
{
	WARN_ONCE( "Enforcing linear precision" );

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
#endif // INTERPOLATION_CONSTRAINTS

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
