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

//////////////////////////
// SimplexRefinableCell //
//////////////////////////
template< unsigned int Dim >
unsigned int SimplexRefinableCell< Dim >::size( void ) const
{
	unsigned int sz = 0;
	for( unsigned int f=0 ; f<faces() ; f++ ) sz += face( f ).size();
	return sz;
}

template< unsigned int Dim >
SimplexIndex< Dim , unsigned int > SimplexRefinableCell< Dim >::operator[]( unsigned int idx ) const
{
	SimplexIndex< Dim-1 , unsigned int > faceSimplexIndex;
	SimplexIndex< Dim , unsigned int > simplexIndex;
	for( unsigned int f=0 ; f<faces() ; f++ )
	{
		const SimplexRefinableCell< Dim-1 > &faceCell = face(f);
		unsigned int faceSize = faceCell.size();
		if( idx<faceSize )
		{
			faceSimplexIndex = faceCell[idx];
			break;
		}
		else idx -= faceSize;
	}
	for( unsigned int d=0 ; d<Dim ; d++ ) simplexIndex[d] = faceSimplexIndex[d];
	simplexIndex[Dim] = centerIndex();
	return simplexIndex;
};

////////////////////////////////////////////////////
// SimplexRefinableElements::NodeMultiIndex_Index //
////////////////////////////////////////////////////
template< unsigned int Dim , unsigned int Degree >
template< typename SubSimplexIndexFunctor >
SimplexRefinableElements< Dim , Degree >::NodeMultiIndex_Index::NodeMultiIndex_Index( SubSimplexIndexFunctor subSimplexIndexFunctor , unsigned int sz )
{
	for( unsigned int i=0 ; i<sz ; i++ )
	{
		SimplexIndex< Dim , unsigned int > subSimplexIndex = subSimplexIndexFunctor(i);
		for( unsigned int n=0 ; n<SimplexElements< Dim , Degree >::NodeNum ; n++ ) _nodeMap[ NodeMultiIndex( subSimplexIndex , n ) ] = 0;
	}
	unsigned int nodeNum = 0;
	for( auto & [ nodeMultiIndex , idx ] : _nodeMap ) idx = nodeNum++;

	_nodeList.resize( nodeNum );
	for( auto & [ nodeMultiIndex , idx ] : _nodeMap ) _nodeList[idx] = nodeMultiIndex;
}

template< unsigned int Dim , unsigned int Degree >
typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex SimplexRefinableElements< Dim , Degree >::NodeMultiIndex_Index::NodeMultiIndex( SimplexIndex< Dim , unsigned int > subSimplexIndex , unsigned int n )
{
	unsigned int v[Degree];
	SimplexElements< Dim , Degree >::FactorNodeIndex( n , v );
	for( unsigned int d=0 ; d<Degree ; d++ ) v[d] = subSimplexIndex[ v[d] ];
	return SimplexRefinableElements::NodeMultiIndex(v);
}

template< unsigned int Dim , unsigned int Degree >
unsigned int SimplexRefinableElements< Dim , Degree >::NodeMultiIndex_Index::fromMultiIndex( SimplexRefinableElements::NodeMultiIndex nodeMultiIndex ) const
{
	auto iter = _nodeMap.find( nodeMultiIndex );
	if( iter==_nodeMap.cend() ) ERROR_OUT( "node multi-index not in map: " , nodeMultiIndex );
	return iter->second;
}

/////////////////////////////////////////////
// SimplexRefinableElements::EnergyWeights //
/////////////////////////////////////////////
const std::string SimplexRefinableElements<>::EnergyWeights::Names[] =
{
	"mass" ,
	"gradient square norm" ,
	"hessian square norm" ,
	"laplacian square norm" ,
	"cross-face gradient difference" ,
};

//////////////////////////////
// SimplexRefinableElements //
//////////////////////////////
template< unsigned int Dim , unsigned int Degree >
SimplexRefinableElements< Dim , Degree >::SimplexRefinableElements( const SimplexRefinable< Dim > &simplexRefinable )
	: _simplexRefinable(simplexRefinable) , _nodeIndex( [&]( unsigned int idx ){ return simplexRefinable[idx]; } , (unsigned int)_simplexRefinable.size() ) {}

template< unsigned int Dim , unsigned int Degree >
unsigned int SimplexRefinableElements< Dim , Degree >::nodeIndex( NodeMultiIndex nodeIndex ) const { return _nodeIndex.fromMultiIndex( nodeIndex ); }

template< unsigned int Dim , unsigned int Degree >
unsigned int SimplexRefinableElements< Dim , Degree >::nodeIndex( SimplexIndex< Dim , unsigned int > subSimplexIndex , unsigned int n ) const
{
	return nodeIndex( NodeMultiIndex_Index::NodeMultiIndex( subSimplexIndex , n ) );
}

template< unsigned int Dim , unsigned int Degree >
template< typename SystemMatrixFunctor >
Eigen::MatrixXd SimplexRefinableElements< Dim , Degree >::_systemMatrix( SystemMatrixFunctor F ) const
{
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero( size() , size() );
	for( unsigned int s=0 ; s<_simplexRefinable.size() ; s++ )
	{
		SimplexIndex< Dim , unsigned int > subSimplexIndex = _simplexRefinable[s];
		SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum > _A = F( _simplexRefinable.metric(s) );
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum ; i++ ) for( unsigned int j=0 ; j<SimplexElements< Dim , Degree >::NodeNum ; j++ )
			A( nodeIndex( NodeMultiIndex_Index::NodeMultiIndex( subSimplexIndex , i ) ) , nodeIndex( NodeMultiIndex_Index::NodeMultiIndex( subSimplexIndex , j ) ) ) += _A(i,j);
	}
	return A;
}

template< unsigned int Dim , unsigned int Degree >
Eigen::MatrixXd SimplexRefinableElements< Dim , Degree >::systemMatrix( typename SimplexRefinableElements< Dim >::EnergyType eType ) const
{
	typedef SimplexRefinableElements<>::EnergyWeights EnergyWeights;
	Eigen::MatrixXd E = Eigen::MatrixXd::Zero( size() , size() );
	if( eType.weights[ EnergyWeights::MASS                           ]>0 ) E += massMatrix()                                                        * eType.weights[ EnergyWeights::MASS                           ];
	if( eType.weights[ EnergyWeights::GRADIENT_SQUARE_NORM           ]>0 ) E += gradientSquareNormMatrix()                                          * eType.weights[ EnergyWeights::GRADIENT_SQUARE_NORM           ];
	if( eType.weights[ EnergyWeights::HESSIAN_SQUARE_NORM            ]>0 ) E += hessianSquareNormMatrix()                                           * eType.weights[ EnergyWeights::HESSIAN_SQUARE_NORM            ];
	if( eType.weights[ EnergyWeights::LAPLACIAN_SQUARE_NORM          ]>0 ) E += laplacianSquareNormMatrix()                                         * eType.weights[ EnergyWeights::LAPLACIAN_SQUARE_NORM          ];
	if( eType.weights[ EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE ]>0 ) E += crossFaceGradientDifferenceMatrix( eType.isIntegrationFaceFunctor ) * eType.weights[ EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE ];
	return E;
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int EmbeddingDimension >
Eigen::MatrixXd SimplexRefinableElements< Dim , Degree >::systemMatrix( typename SimplexRefinableElements< Dim >::EnergyType eType , std::function< Point< double , EmbeddingDimension > ( unsigned int idx ) > positionFunctor ) const
{
	typedef SimplexRefinableElements<>::EnergyWeights EnergyWeights;
	Eigen::MatrixXd E = Eigen::MatrixXd::Zero( size() , size() );
	if( eType.weights[ EnergyWeights::MASS                                     ]>0 ) E += massMatrix()                                                                                   * eType.weights[ EnergyWeights::MASS                                     ];
	if( eType.weights[ EnergyWeights::GRADIENT_SQUARE_NORM                     ]>0 ) E += gradientSquareNormMatrix()                                                                     * eType.weights[ EnergyWeights::GRADIENT_SQUARE_NORM                     ];
	if( eType.weights[ EnergyWeights::HESSIAN_SQUARE_NORM                      ]>0 ) E += hessianSquareNormMatrix()                                                                      * eType.weights[ EnergyWeights::HESSIAN_SQUARE_NORM                      ];
	if( eType.weights[ EnergyWeights::LAPLACIAN_SQUARE_NORM                    ]>0 ) E += laplacianSquareNormMatrix()                                                                    * eType.weights[ EnergyWeights::LAPLACIAN_SQUARE_NORM                    ];
	if( eType.weights[ EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE           ]>0 ) E += crossFaceGradientDifferenceMatrix( eType.isIntegrationFaceFunctor )                            * eType.weights[ EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE           ];
	return E;
}

template< unsigned int Dim , unsigned int Degree >
Eigen::MatrixXd SimplexRefinableElements< Dim , Degree >::massMatrix( void ) const
{
	return _systemMatrix( SimplexElements< Dim , Degree >::MassMatrix );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::MatrixXd SimplexRefinableElements< Dim , Degree >::gradientSquareNormMatrix( void ) const
{
	return _systemMatrix( SimplexElements< Dim , Degree >::GradientSquareNormMatrix );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::MatrixXd SimplexRefinableElements< Dim , Degree >::hessianSquareNormMatrix( void ) const
{
	return _systemMatrix( SimplexElements< Dim , Degree >::HessianSquareNormMatrix );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::MatrixXd SimplexRefinableElements< Dim , Degree >::laplacianSquareNormMatrix( void ) const
{
	return _systemMatrix( SimplexElements< Dim , Degree >::LaplacianSquareNormMatrix );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::MatrixXd SimplexRefinableElements< Dim , Degree >::crossFaceGradientDifferenceMatrix( std::function< bool ( FaceMultiIndex ) >  faceSelectionFunctor ) const
{
	typedef Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > FaceD;
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero( size() , size() );
	if constexpr( Dim==1 )
	{
//		WARN_ONCE( "Shouldn't be called for Dim==1" );
		return C;
	}
	else
	{
		// Functionality for returning the multi-index (and parity) associated with a face of a simplex
		auto GetFaceMultiIndex = [&]( SimplexIndex< Dim , unsigned int > simplexIndex , unsigned int f , Permutation< Dim > &p )
		{
			SimplexIndex< Dim-1 , unsigned int > faceIndex = RightSimplex< Dim >::Face( f );
			p = Permutation< Dim >( [&]( unsigned int i , unsigned int j ){ return simplexIndex[ faceIndex[i] ]<simplexIndex[ faceIndex[j] ]; } );
			for( unsigned int d=0 ; d<Dim ; d++ ) faceIndex[d] = simplexIndex[ faceIndex[d] ];
			return FaceMultiIndex( &faceIndex[0] );
		};

		// A mapping from all simplex faces (including those we will not integrate over) to indices
		std::map< FaceMultiIndex , unsigned int > faceMap;
		{
			// Iterate over the simplices
			for( unsigned int i=0 ; i<_simplexRefinable.size() ; i++ )
			{
				Permutation< Dim > p;
				SimplexIndex< Dim , unsigned int > simplexIndex = _simplexRefinable[i];
				// Iterate over the faces of each simplex
				for( unsigned int f=0 ; f<=Dim ; f++ ) faceMap[ GetFaceMultiIndex( simplexIndex , f , p ) ] = 0;
			}
			unsigned int faceCount = 0;
			for( auto & [ mIdx , idx ] : faceMap ) idx = faceCount++;
		}

		// Per face, a map giving the gradient components of a node
		std::vector< std::map< unsigned int , FaceD > > faceGradientOrthogonalComponents( faceMap.size() );
		// Per face, the metric on the face
		std::vector< SquareMatrix< double , Dim-1 > > faceMetrics( faceMap.size() );

		// Iterate over the simplices
		for( unsigned int i=0 ; i<_simplexRefinable.size() ; i++ )
		{
			SimplexIndex< Dim , unsigned int > simplexIndex = _simplexRefinable[i];
			SquareMatrix< double , Dim > g = _simplexRefinable.metric(i);

			// Compute the gradient components and metrics for the faces of the simplex
			Matrix< FaceD , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > _faceGradientOrthogonalComponents;
			Point< SquareMatrix< double , Dim-1 > , Dim+1 > _faceMetrics;

			_faceGradientOrthogonalComponents = SimplexElements< Dim , Degree >::FaceGradientOrthogonalComponents( g , &simplexIndex[0] );
			_faceMetrics = SimplexElements< Dim , Degree >::FaceMetrics( g , &simplexIndex[0] );

			for( unsigned int f=0 ; f<=Dim ; f++ )
			{
				Permutation< Dim > p;
				FaceMultiIndex fmi = GetFaceMultiIndex( simplexIndex , f , p );
				if( faceSelectionFunctor( fmi ) )
				{
					bool evenParity = ( p.parity()&1 )==0;
					unsigned int fi = faceMap[fmi];
					Matrix< double , Dim , Dim-1 > A = RightSimplex< Dim-1 >::AffineTransform( p );

					// Accumulate the metrics and (signed) gradient components
					faceMetrics[fi] += _faceMetrics[f]/2;
					for( unsigned int n=0 ; n<SimplexElements< Dim , Degree >::NodeNum ; n++ )
					{
						Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > &gOrthoComponents = faceGradientOrthogonalComponents[fi][ nodeIndex( simplexIndex , n ) ];
						for( unsigned int d=0 ; d<Dim ; d++ )
							if( evenParity ) gOrthoComponents[d] += _faceGradientOrthogonalComponents(n,f)[d].template operator()< Dim >(A);
							else             gOrthoComponents[d] -= _faceGradientOrthogonalComponents(n,f)[d].template operator()< Dim >(A);
					}
				}
			}
		}

		// Iterate over each face
		for( unsigned int f=0 ; f<faceGradientOrthogonalComponents.size() ; f++ )
			// Iterate over each (node indices , gradient components) pair associated to the face
			for( const auto & [ i1 , g1 ] : faceGradientOrthogonalComponents[f] ) for( const auto & [ i2 , g2 ] : faceGradientOrthogonalComponents[f] )
				for( unsigned int d=0 ; d<Dim ; d++ ) C(i1,i2) += RightSimplex< Dim-1 >::Integral( g1[d] * g2[d] , faceMetrics[f] );

		return C;
	}
}

/////////////////////////////////////
// InterpolatingProlongationSystem //
/////////////////////////////////////
template< bool PoU >
template< typename CoarseSelectionFunctor >
InterpolatingProlongationSystem< PoU >::InterpolatingProlongationSystem( const Eigen::MatrixXd &E , CoarseSelectionFunctor coarseSelectionFunctor )
{
	std::vector< unsigned int > coarseIndices;
	for( unsigned int i=0 ; i<E.rows() ; i++ ) if( coarseSelectionFunctor(i) ) coarseIndices.push_back( i );
	std::vector< Constraint > constraints;
	_init< 0 >( E , coarseIndices , constraints , NULL );
}

template< bool PoU >
InterpolatingProlongationSystem< PoU >::InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices )
{
	// Do some sanity checking
	for( unsigned int i=0 ; i<coarseIndices.size() ; i++ )
	{
		if( coarseIndices[i]>=E.rows() ) ERROR_OUT( "Coarse index cannot exceed fine index: " , coarseIndices[i] , " >= " , E.rows() );
		for( unsigned int j=0 ; j<i ; j++ )
			if( coarseIndices[i]==coarseIndices[j] ) ERROR_OUT( "Duplicated coarse index: " , i , " , " , j , " ->"  , coarseIndices[i] );
	}
	std::vector< Constraint > constraints;
	_init< 0 >( E , coarseIndices , constraints , NULL );
}

template< bool PoU >
InterpolatingProlongationSystem< PoU >::InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , const std::vector< Constraint > &constraints )
{
	// Do some sanity checking
	for( unsigned int i=0 ; i<coarseIndices.size() ; i++ )
	{
		if( coarseIndices[i]>=E.rows() ) ERROR_OUT( "Coarse index cannot exceed fine index: " , coarseIndices[i] , " >= " , E.rows() );
		for( unsigned int j=0 ; j<i ; j++ )
			if( coarseIndices[i]==coarseIndices[j] ) ERROR_OUT( "Duplicated coarse index: " , i , " , " , j , " ->"  , coarseIndices[i] );
	}
	_init< 0 >( E , coarseIndices , constraints , NULL );
}

template< bool PoU >
template< unsigned int InterpolationDim >
InterpolatingProlongationSystem< PoU >::InterpolatingProlongationSystem( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , const Point< double , InterpolationDim > *interpolationConstraints )
{
	// Do some sanity checking
	for( unsigned int i=0 ; i<coarseIndices.size() ; i++ )
	{
		if( coarseIndices[i]>=E.rows() ) ERROR_OUT( "Coarse index cannot exceed fine index: " , coarseIndices[i] , " >= " , E.rows() );
		for( unsigned int j=0 ; j<i ; j++ )
			if( coarseIndices[i]==coarseIndices[j] ) ERROR_OUT( "Duplicated coarse index: " , i , " , " , j , " ->"  , coarseIndices[i] );
	}
	std::vector< Constraint > constraints;
	_init< InterpolationDim >( E , coarseIndices , constraints , interpolationConstraints );
}

template< bool PoU >
template< unsigned int InterpolationDim >
void InterpolatingProlongationSystem< PoU >::_init( const Eigen::MatrixXd &E , const std::vector< unsigned int > &coarseIndices , const std::vector< Constraint > &constraints , const Point< double , InterpolationDim > *interpolationConstraints )
{
	_coarseIndices = coarseIndices;
	_fineIndices.reserve( E.rows() - _coarseIndices.size() );

	std::vector< bool > isCoarse( E.rows() , false );
	std::vector< unsigned int > fineIndexMap( E.rows() , -1 );
	for( unsigned int i=0 ; i<_coarseIndices.size() ; i++ ) isCoarse[ _coarseIndices[i] ] = true;
	for( unsigned int i=0 ; i<E.rows() ; i++ ) if( !isCoarse[i] )
	{
		fineIndexMap[i] = (unsigned int)_fineIndices.size();
		_fineIndices.push_back( i );
	}
	_coarseDim = (unsigned int)_coarseIndices.size();
	_fineDim   = (unsigned int)  _fineIndices.size();

	// We would like to solve for the prolongation coeffiicents P(i,j) (with 0<=i<_fineDim and 0<=j<_coarseDim)
	// that minimizes the energy:
	//		E(P) =   \sum_{i=0}^{_coarseDim-1} \sum_{j=0}^{  _fineDim-1} \sum_{k=0}^{  _fineDim-1}  E(   _fineIndices[j] ,   _fineIndices[k] ) * P(j,i) * P(k,i)
	// 	         + 2 \sum_{i=0}^{_coarseDim-1} \sum_{j=0}^{  _fineDim-1} \sum_{k=0}^{_coarseDim-1}  E(   _fineIndices[j] , _coarseIndices[k] ) * P(j,i) * d(k,i)
	// 	         +   \sum_{i=0}^{_coarseDim-1} \sum_{j=0}^{_coarseDim-1} \sum_{k=0}^{_coarseDim-1}  E( _coarseIndices[j] , _coarseIndices[k] ) * d(j,i) * d(k,i)
	// Subject to the (_fineDim) PoU constraints:
	//		\sum_{i=0}^{coarseDim-1} P(j,i) = 1 for all j
	// And the additional (constraints.size) constraints provided as input

	std::vector< Eigen::Triplet< double > > qEntries , cEntries;

	unsigned int dim = _coarseDim * _fineDim;

	// The quadratic energy term
	_Q.resize( dim , dim );
	_q.resize( dim );

	{
		unsigned int count = 0;
		for( unsigned int j=0 ; j<_fineDim ; j++ ) for( unsigned int k=0 ; k<_fineDim ; k++ ) if( E( _fineIndices[j] , _fineIndices[k] ) ) count += _coarseDim;
		qEntries.reserve( _coarseDim );
	}

	for( unsigned int j=0 ; j<_fineDim ; j++ ) for( unsigned int k=0 ; k<_fineDim ; k++ ) if( E( _fineIndices[j] , _fineIndices[k] ) )
		for( unsigned int i=0 ; i<_coarseDim ; i++ ) qEntries.push_back( Eigen::Triplet< double >( _index(j,i) , _index(k,i) , E( _fineIndices[j] , _fineIndices[k] ) ) );
	for( unsigned int j=0 ; j<_fineDim ; j++ ) for( unsigned int k=0 ; k<_coarseDim ; k++ ) _q[ _index(j,k) ] = E( _fineIndices[j] , _coarseIndices[k] );

	// The linear interpolation constraint
	if( interpolationConstraints )
	{
		cEntries.reserve( _fineDim * _coarseDim + constraints.size() + _fineDim*_coarseDim*InterpolationDim );
		if constexpr( PoU )
		{
			_C.resize( _fineDim + constraints.size() + _fineDim*InterpolationDim , dim );
			_c = Eigen::VectorXd::Zero( _fineDim + constraints.size() + _fineDim*InterpolationDim );
		}
		else
		{
			_C.resize( constraints.size() + _fineDim*InterpolationDim , dim );
			_c = Eigen::VectorXd::Zero( constraints.size() + _fineDim*InterpolationDim );
		}
	}
	else
	{
		cEntries.reserve( _fineDim * _coarseDim + constraints.size() );
		if constexpr( PoU )
		{
			_C.resize( _fineDim + constraints.size() , dim );
			_c = Eigen::VectorXd::Zero( _fineDim + constraints.size() );
		}
		else
		{
			_C.resize( constraints.size() , dim );
			_c = Eigen::VectorXd::Zero( constraints.size() );
		}
	}

	if constexpr( PoU )
	{
		// The PoU constraints
		for( unsigned int i=0 ; i<_fineDim ; i++ )
		{
			for( unsigned int j=0 ; j<_coarseDim ; j++ ) cEntries.push_back( Eigen::Triplet< double >( i , _index(i,j) , 1 ) );
			_c[i] = 1;
		}
	}

	// Locked weight constraints
	for( unsigned int i=0 ; i<constraints.size() ; i++ )
	{
		unsigned int row = PoU ? _fineDim + i : i;
		if( isCoarse[ constraints[i].fineIndex ] )
		{
			if     ( constraints[i].fineIndex==constraints[i].coarseIndex && constraints[i].value!=1. ) ERROR_OUT( "Bad constraint" );
			else if( constraints[i].fineIndex!=constraints[i].coarseIndex && constraints[i].value!=0. ) ERROR_OUT( "Bad constraint" );
		}
		else
		{
			cEntries.push_back( Eigen::Triplet< double >( row , _index( fineIndexMap[ constraints[i].fineIndex ] , constraints[i].coarseIndex ) , 1 ) );
			_c[ row ] = constraints[i].value;
		}
	}

	// Interpolation constraints
	// For the i-th node, we would like to satisfy:
	//		\sum_j P( _coarseIndices[j] , _fineIndices[i] ) * internalConstraints[ _coarseIndices[j] ] = internalConstraints[ _fineIndices[i] ]
	if constexpr( InterpolationDim!=0 ) if( interpolationConstraints )
	{
		unsigned int start = PoU ? _fineDim + (unsigned int)constraints.size() : (unsigned int)constraints.size();
		for( unsigned int f=0 ; f<_fineDim ; f++ ) for( unsigned int d=0 ; d<InterpolationDim ; d++ )
		{
			unsigned int row = start + f*InterpolationDim + d;
			for( unsigned int c=0 ; c<_coarseDim ; c++ )
				cEntries.push_back( Eigen::Triplet< double >( row , _index(f,c) , interpolationConstraints[ _coarseIndices[c] ][d] ) );
			_c[row] = interpolationConstraints[ _fineIndices[f] ][d];
		}
	}

	_Q.setFromTriplets( qEntries.begin() , qEntries.end() );
	_C.setFromTriplets( cEntries.begin() , cEntries.end() );

	if( _c.size() ) _lcqo = LCQO( _Q , _q , _C , _c , true );
	else            _lcqo = LCQO( _Q , _q );
}

template< bool PoU >
unsigned int InterpolatingProlongationSystem< PoU >::_index( unsigned int fine , unsigned int coarse ) const
{
	return fine + coarse*_fineDim;
}

template< bool PoU >
Eigen::MatrixXd InterpolatingProlongationSystem< PoU >::_toMatrix( const Eigen::VectorXd &v ) const
{
	Eigen::MatrixXd P( _fineDim+_coarseDim , _coarseDim );
	for( unsigned int i=0 ; i<_fineDim ; i++ ) for( unsigned int j=0 ; j<_coarseDim ; j++ ) P( _fineIndices[i] , j ) = v[ _index(i,j) ];
	for( unsigned int i=0 ; i<_coarseDim ; i++ ) for( unsigned int j=0 ; j<_coarseDim ; j++ ) P( _coarseIndices[i] , j ) = i==j ? 1 : 0;
	return P;
}

template< bool PoU >
Eigen::VectorXd InterpolatingProlongationSystem< PoU >::_toVector( const Eigen::MatrixXd &P ) const
{
	Eigen::VectorXd v( _fineDim * _coarseDim );
	for( unsigned int i=0 ; i<_fineDim ; i++ ) for( unsigned int j=0 ; j<_coarseDim ; j++ ) v[ _index(i,j) ] = P( _fineIndices[i] , j );
	return v;
}

template< bool PoU >
template< bool StableSolve >
Eigen::MatrixXd InterpolatingProlongationSystem< PoU >::prolongation( void ) const
{
	return _toMatrix( _lcqo.solve< StableSolve >() );
}

template< bool PoU >
double InterpolatingProlongationSystem< PoU >::energy( const Eigen::MatrixXd &P ) const
{
	if( P.rows()!=(_fineDim+_coarseDim) || P.cols()!=_coarseDim ) ERROR_OUT( "Bad dimensions: " , P.rows() , " x " , P.cols() , " != " , (_fineDim+_coarseDim) , " x " , _coarseDim );
	Eigen::VectorXd v = _toVector( P );
	Eigen::VectorXd Qv = _Q * v;
	double e = 0;
	for( unsigned int i=0 ; i<v.size() ; i++ ) e += v[i] * Qv[i] + 2 * v[i] * _q[i];
	return e;
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree >
void InterpolatingProlongationSystem< PoU >::HierarchicalProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , typename SimplexRefinableElements<>::EnergyWeights eWeights , ProlongationInfo< Degree > pInfo[Dim] , unsigned int finestDim )
{
	SimplexRefinableElements< Dim , Degree > sre( simplexRefinableCell );
	_HierarchicalProlongation< Dim , Degree >( simplexRefinableCell , sre , eWeights , &pInfo[0] , finestDim );
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
void InterpolatingProlongationSystem< PoU >::HierarchicalProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , ProlongationInfo< Degree > pInfo[Dim] , unsigned int finestDim )
{
	SimplexRefinableElements< Dim , Degree > sre( simplexRefinableCell );
	_HierarchicalProlongation< Dim , Degree , EmbeddingDimension >( simplexRefinableCell , sre , eWeights , positionFunctor , planarityEpsilon , &pInfo[0] , finestDim );
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree >
void InterpolatingProlongationSystem< PoU >::_HierarchicalProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , const SimplexRefinableElements< Dim , Degree > &sre , typename SimplexRefinableElements<>::EnergyWeights eWeights , ProlongationInfo< Degree > *pInfo , unsigned int finestDim )
{
	typedef typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	static const unsigned int CoarseDim = Dim-1;
	// Our goals are:
	// 1. To compute the prolognation from CoarseDim -> Dim
	// 2. Merge the prolongations from coarser dimensions

	// The map mapping multi-indices to indices for nodes residing on the faces of the cell
	std::map< NodeMultiIndex , unsigned int > boundaryNodeMap;

	std::vector< NodeMultiIndex > _coarseMultiIndices;
	std::vector< NodeMultiIndex > &coarseMultiIndices = CoarseDim<=finestDim ? pInfo[CoarseDim].coarseMultiIndices : _coarseMultiIndices;

	// (Recursively) compute the prolongation at the coarser resolution
	if constexpr( Dim!=1 )
	{
		// The hierarchical prolongation information for each face
		std::vector< std::vector< ProlongationInfo< Degree > > > _pInfo( simplexRefinableCell.faces() );

		// The maps mapping multi-indices to indices for nodes residing on the sub-faces of the cell
		std::vector< std::map< NodeMultiIndex , unsigned int > > nodeMaps( finestDim+1 );

		// For each face, the list of the nodes defined on that face
		std::vector< std::vector< NodeMultiIndex > > boundaryMultiIndices( simplexRefinableCell.faces() );

		for( unsigned int f=0 ; f<simplexRefinableCell.faces() ; f++ )
		{
			const SimplexRefinableCell< CoarseDim > &_simplexRefinableCell = simplexRefinableCell.face(f);
			SimplexRefinableElements< CoarseDim , Degree > _sre( _simplexRefinableCell );

			// Compute the prolongation information for the coarser resolutions
			_pInfo[f].resize( finestDim+1 );
			_HierarchicalProlongation( _simplexRefinableCell , _sre , eWeights , &_pInfo[f][0] , finestDim );

			// Merge the list of node multi-indices
			for( unsigned int d=0 ; d<CoarseDim && d<=finestDim ; d++ )	for( auto nmi : _pInfo[f][d].coarseMultiIndices ) nodeMaps[d][nmi] = 0;

			boundaryMultiIndices[f].resize( _sre.size() );
			for( unsigned int i=0 ; i<_sre.size() ; i++ ) boundaryNodeMap[ _sre[i] ] = 0 , boundaryMultiIndices[f][i] = _sre[i];
		}

		// Assign indices to the nodes at different levels
		{
			unsigned int count = 0;
			for( auto &[nmi,idx] : boundaryNodeMap ) idx = count++;
			for( unsigned int d=0 ; d<CoarseDim && d<=finestDim ; d++ ) 
			{
				unsigned int count = 0;
				for( auto &[nmi,idx] : nodeMaps[d] ) idx = count++ ;
			}
		}

		// Set the accumulated list of coarse node indices
		for( unsigned int d=0 ; d<CoarseDim && d<=finestDim ; d++ )
		{
			pInfo[d].coarseMultiIndices.resize( nodeMaps[d].size() );
			for( const auto &[nmi,idx] : nodeMaps[d] ) pInfo[d].coarseMultiIndices[idx] = nmi;
		}

		// Accumulate the prolongation matrix information
		for( unsigned int d=0 ; d<CoarseDim && d<=finestDim ; d++ )
		{
			auto Index = []( const std::map< NodeMultiIndex , unsigned int > &nodeMap , NodeMultiIndex nmi )
			{
				auto iter = nodeMap.find(nmi);
				if( iter==nodeMap.end() ) ERROR_OUT( "Could not find node multi index: " , nmi );
				return iter->second;
			};

			const std::map< NodeMultiIndex , unsigned int > &coarseNodeMap = nodeMaps[d];
			const std::map< NodeMultiIndex , unsigned int > &  fineNodeMap = ( ( d+1<CoarseDim && d+1<=finestDim ) ? nodeMaps[d+1] : boundaryNodeMap );
			unsigned int  lowDim = (unsigned int)coarseNodeMap.size();
			unsigned int highDim = (unsigned int)  fineNodeMap.size();
			Eigen::MatrixXd &P = pInfo[d].P;
			Eigen::MatrixXi  C = Eigen::MatrixXi::Zero( highDim , lowDim );
			P = Eigen::MatrixXd::Zero( highDim , lowDim );

			for( unsigned f=0 ; f<_pInfo.size() ; f++ )
			{
				const std::vector< ProlongationInfo< Degree > > &facePInfo = _pInfo[f];
				const ProlongationInfo< Degree > &__pInfo = facePInfo[d];
				const std::vector< NodeMultiIndex > &coareNodeMultiIndices = __pInfo.coarseMultiIndices;
				const std::vector< NodeMultiIndex > & fineNodeMultiIndices = ( ( d+1<CoarseDim && d+1<=finestDim ) ? facePInfo[d+1].coarseMultiIndices : boundaryMultiIndices[f] );
				for( unsigned int c=0 ; c<__pInfo.P.cols() ; c++ ) for( unsigned int r=0 ; r<__pInfo.P.rows() ; r++ )
				{
					unsigned int _c = Index( coarseNodeMap , coareNodeMultiIndices[c] ) , _r = Index( fineNodeMap , fineNodeMultiIndices[r] );
					P(_r,_c) += __pInfo.P(r,c);
					C(_r,_c)++;
				}
			}

			for( unsigned int c=0 ; c<C.cols() ; c++ ) for( unsigned int r=0 ; r<C.rows() ; r++ ) if( C(r,c) ) P(r,c) /= C(r,c);
		}
	}
	else
	{
		auto VertexNodeIndex = []( unsigned int v )
		{
			unsigned int idx[Degree];
			for( unsigned int d=0 ; d<Degree ; d++ ) idx[d] = v;
			return NodeMultiIndex( idx );
		};
		for( unsigned int i=0 ; i<2 ; i++ ) boundaryNodeMap[ VertexNodeIndex( simplexRefinableCell[0][i] ) ] = i;
	}
	{
		coarseMultiIndices.resize( boundaryNodeMap.size() );
		for( const auto &[nmi,idx] : boundaryNodeMap ) coarseMultiIndices[idx] = nmi;	
	}

	// Then compute the prolongation at the finest level
	std::vector< unsigned int > coarseIndices( coarseMultiIndices.size() );
	for( unsigned int i=0 ; i<coarseMultiIndices.size() ; i++ ) coarseIndices[i] = sre.nodeIndex( coarseMultiIndices[i] );
	if( CoarseDim<=finestDim ) pInfo[CoarseDim].P = _BoundaryProlongation( simplexRefinableCell , sre , eWeights , coarseIndices );
	else                       pInfo[finestDim].P = _BoundaryProlongation( simplexRefinableCell , sre , eWeights , coarseIndices , pInfo[finestDim].P );
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree >
Eigen::MatrixXd InterpolatingProlongationSystem< PoU >::_BoundaryProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , const SimplexRefinableElements< Dim , Degree > &sre , typename SimplexRefinableElements<>::EnergyWeights eWeights , const std::vector< unsigned int > &coarseIndices )
{
	typedef typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex NodeMultiIndex;

	// Construct the energy description from the weights
	typename SimplexRefinableElements< Dim >::EnergyType eType(eWeights);
	unsigned int centerIndex = -1;
	if constexpr( Dim>1 ) centerIndex = simplexRefinableCell.centerIndex();
	eType.isIntegrationFaceFunctor = [&]( const typename SimplexRefinableElements< Dim >::FaceMultiIndex &fmi )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) if( fmi[d]==centerIndex ) return true;
		return false;
	};

	// Get the energy matrix
	Eigen::MatrixXd E = sre.systemMatrix( eType );

	// Get the prolongation
	return InterpolatingProlongationSystem( E , coarseIndices ).prolongation< false >();
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree >
Eigen::MatrixXd InterpolatingProlongationSystem< PoU >::_BoundaryProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , const SimplexRefinableElements< Dim , Degree > &sre , typename SimplexRefinableElements<>::EnergyWeights eWeights , const std::vector< unsigned int > &coarseIndices , const Eigen::MatrixXd &coarseP )
{
	// Conceptually, we can break up the set of nodes into three parts:
	// 1. base nodes (# = coarseP.cols)
	// 2. fine nodes (# = sre.size - coarseP.rows)
	// 3. middle nodes (# = coarseP.rows - coarseP.cols)
	// We denote by "coarse" the union of "base" and "middle" nodes (# = coarseIndices.size = coarseP.rows)

	typedef typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	if( coarseP.rows()!=coarseIndices.size() ) ERROR_OUT( "Dimensional mismatch: " , coarseP.rows() , " != "  , coarseIndices.size() );

	// Construct the energy description from the weights
	typename SimplexRefinableElements< Dim >::EnergyType eType(eWeights);
	unsigned int centerIndex = -1;
	if constexpr( Dim>1 ) centerIndex = simplexRefinableCell.centerIndex();
	eType.isIntegrationFaceFunctor = [&]( const typename SimplexRefinableElements< Dim >::FaceMultiIndex &fmi )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) if( fmi[d]==centerIndex ) return true;
		return false;
	};

	Eigen::MatrixXd fineE = sre.systemMatrix( eType );

	// Compute the subset of element indices that are not in the coarse index set.
	std::vector< unsigned int > fineIndices;
	{
		fineIndices.reserve( sre.size() - coarseIndices.size() );
		std::vector< bool > isCoarse( sre.size() , false );
		for( unsigned int i=0 ; i<coarseIndices.size() ; i++ ) isCoarse[ coarseIndices[i] ] = true;
		for( unsigned int i=0 ; i<isCoarse.size() ; i++ ) if( !isCoarse[i] ) fineIndices.push_back( i );
	}

	// Construct the prolongation matrix: base + fine -> base + middle + fine
	Eigen::MatrixXd _P = Eigen::MatrixXd::Zero( fineE.cols() , coarseP.cols() + fineIndices.size() );
	for( unsigned int c=0 ; c<coarseP.cols() ; c++ ) for( unsigned int r=0 ; r<coarseP.rows() ; r++ ) _P( coarseIndices[r] , c ) = coarseP( r , c );
	for( unsigned int i=0 ; i<fineIndices.size() ; i++ ) _P( fineIndices[i] , coarseP.cols()+i ) = 1.;

	// Get the energy matrix for: base + fine
	Eigen::MatrixXd E = _P.transpose() * fineE * _P;

	// Get the relative prolongation: coarse -> coase + fine
	std::vector< unsigned int > _coarseIndices( coarseP.cols() );
	for( unsigned int i=0 ; i<_coarseIndices.size() ; i++ ) _coarseIndices[i] = i;
	_P = InterpolatingProlongationSystem( E , _coarseIndices ).prolongation< false >();

	// Extend to the complete prolongation: coarse -> coarse + middle + fine
	Eigen::MatrixXd P = Eigen::MatrixXd::Zero( sre.size() , _P.cols() );
	for( unsigned int c=0 ; c<_P.cols() ; c++ )
	{
		for( unsigned int r=0 ; r<coarseP.rows() ; r++ ) P( coarseIndices[r] , c ) = coarseP( r , c );
		for( unsigned int r=0 ; r<fineIndices.size() ; r++ ) P( fineIndices[r] , c ) = _P( coarseP.cols()+r , c );
	}

	return P;
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
void InterpolatingProlongationSystem< PoU >::_HierarchicalProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , const SimplexRefinableElements< Dim , Degree > &sre , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , ProlongationInfo< Degree > *pInfo , unsigned int finestDim )
{
	typedef typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	static const unsigned int CoarseDim = Dim-1;
	// Our goals are:
	// 1. To compute the prolognation from CoarseDim -> Dim
	// 2. Merge the prolongations from coarser dimensions

	// The map mapping multi-indices to indices for nodes residing on the faces of the cell
	std::map< NodeMultiIndex , unsigned int > boundaryNodeMap;

	std::vector< NodeMultiIndex > _coarseMultiIndices;
	std::vector< NodeMultiIndex > &coarseMultiIndices = CoarseDim<=finestDim ? pInfo[CoarseDim].coarseMultiIndices : _coarseMultiIndices;

	// (Recursively) compute the prolongation at the coarser resolution
	if constexpr( Dim!=1 )
	{
		// The hierarchical prolongation information for each face
		std::vector< std::vector< ProlongationInfo< Degree > > > _pInfo( simplexRefinableCell.faces() );

		// The maps mapping multi-indices to indices for nodes residing on the sub-faces of the cell
		std::vector< std::map< NodeMultiIndex , unsigned int > > nodeMaps( finestDim+1 );

		// For each face, the list of the nodes defined on that face
		std::vector< std::vector< NodeMultiIndex > > boundaryMultiIndices( simplexRefinableCell.faces() );

		for( unsigned int f=0 ; f<simplexRefinableCell.faces() ; f++ )
		{
			const SimplexRefinableCell< CoarseDim > &_simplexRefinableCell = simplexRefinableCell.face(f);
			SimplexRefinableElements< CoarseDim , Degree > _sre( _simplexRefinableCell );

			// Compute the prolongation information for the coarser resolutions
			_pInfo[f].resize( finestDim+1 );
			_HierarchicalProlongation( _simplexRefinableCell , _sre , eWeights , positionFunctor , planarityEpsilon , &_pInfo[f][0] , finestDim );

			// Merge the list of node multi-indices
			for( unsigned int d=0 ; d<CoarseDim && d<=finestDim ; d++ )	for( auto nmi : _pInfo[f][d].coarseMultiIndices ) nodeMaps[d][nmi] = 0;

			boundaryMultiIndices[f].resize( _sre.size() );
			for( unsigned int i=0 ; i<_sre.size() ; i++ ) boundaryNodeMap[ _sre[i] ] = 0 , boundaryMultiIndices[f][i] = _sre[i];
		}

		// Assign indices to the nodes at different levels
		{
			unsigned int count = 0;
			for( auto &[nmi,idx] : boundaryNodeMap ) idx = count++;
			for( unsigned int d=0 ; d<CoarseDim && d<=finestDim ; d++ ) 
			{
				unsigned int count = 0;
				for( auto &[nmi,idx] : nodeMaps[d] ) idx = count++ ;
			}
		}

		// Set the accumulated list of coarse node indices
		for( unsigned int d=0 ; d<CoarseDim && d<=finestDim ; d++ )
		{
			pInfo[d].coarseMultiIndices.resize( nodeMaps[d].size() );
			for( const auto &[nmi,idx] : nodeMaps[d] ) pInfo[d].coarseMultiIndices[idx] = nmi;
		}

		// Accumulate the prolongation matrix information
		for( unsigned int d=0 ; d<CoarseDim && d<=finestDim ; d++ )
		{
			auto Index = []( const std::map< NodeMultiIndex , unsigned int > &nodeMap , NodeMultiIndex nmi )
			{
				auto iter = nodeMap.find(nmi);
				if( iter==nodeMap.end() ) ERROR_OUT( "Could not find node multi index: " , nmi );
				return iter->second;
			};

			const std::map< NodeMultiIndex , unsigned int > &coarseNodeMap = nodeMaps[d];
			const std::map< NodeMultiIndex , unsigned int > &  fineNodeMap = ( ( d+1<CoarseDim && d+1<=finestDim ) ? nodeMaps[d+1] : boundaryNodeMap );
			unsigned int  lowDim = (unsigned int)coarseNodeMap.size();
			unsigned int highDim = (unsigned int)  fineNodeMap.size();
			Eigen::MatrixXd &P = pInfo[d].P;
			Eigen::MatrixXi  C = Eigen::MatrixXi::Zero( highDim , lowDim );
			P = Eigen::MatrixXd::Zero( highDim , lowDim );

			for( unsigned f=0 ; f<_pInfo.size() ; f++ )
			{
				const std::vector< ProlongationInfo< Degree > > &facePInfo = _pInfo[f];
				const ProlongationInfo< Degree > &__pInfo = facePInfo[d];
				const std::vector< NodeMultiIndex > &coareNodeMultiIndices = __pInfo.coarseMultiIndices;
				const std::vector< NodeMultiIndex > & fineNodeMultiIndices = ( ( d+1<CoarseDim && d+1<=finestDim ) ? facePInfo[d+1].coarseMultiIndices : boundaryMultiIndices[f] );
				for( unsigned int c=0 ; c<__pInfo.P.cols() ; c++ ) for( unsigned int r=0 ; r<__pInfo.P.rows() ; r++ )
				{
					unsigned int _c = Index( coarseNodeMap , coareNodeMultiIndices[c] ) , _r = Index( fineNodeMap , fineNodeMultiIndices[r] );
					P(_r,_c) += __pInfo.P(r,c);
					C(_r,_c)++;
				}
			}

			for( unsigned int c=0 ; c<C.cols() ; c++ ) for( unsigned int r=0 ; r<C.rows() ; r++ ) if( C(r,c) ) P(r,c) /= C(r,c);
		}
	}
	else
	{
		auto VertexNodeIndex = []( unsigned int v )
		{
			unsigned int idx[Degree];
			for( unsigned int d=0 ; d<Degree ; d++ ) idx[d] = v;
			return NodeMultiIndex( idx );
		};
		for( unsigned int i=0 ; i<2 ; i++ ) boundaryNodeMap[ VertexNodeIndex( simplexRefinableCell[0][i] ) ] = i;
	}
	{
		coarseMultiIndices.resize( boundaryNodeMap.size() );
		for( const auto &[nmi,idx] : boundaryNodeMap ) coarseMultiIndices[idx] = nmi;	
	}

	// Then compute the prolongation at the finest level
	std::vector< unsigned int > coarseIndices( coarseMultiIndices.size() );
	for( unsigned int i=0 ; i<coarseMultiIndices.size() ; i++ ) coarseIndices[i] = sre.nodeIndex( coarseMultiIndices[i] );
	if( CoarseDim<=finestDim ) pInfo[CoarseDim].P = _BoundaryProlongation( simplexRefinableCell , sre , eWeights , positionFunctor , planarityEpsilon , coarseIndices );
	else                       pInfo[finestDim].P = _BoundaryProlongation( simplexRefinableCell , sre , eWeights , positionFunctor , planarityEpsilon , coarseIndices , pInfo[finestDim].P );
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
Eigen::MatrixXd InterpolatingProlongationSystem< PoU >::_BoundaryProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , const SimplexRefinableElements< Dim , Degree > &sre , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , const std::vector< unsigned int > &coarseIndices )

{
	typedef typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex NodeMultiIndex;

	// Construct the energy description from the weights
	typename SimplexRefinableElements< Dim >::EnergyType eType(eWeights);
	unsigned int centerIndex = -1;
	if constexpr( Dim>1 ) centerIndex = simplexRefinableCell.centerIndex();
	eType.isIntegrationFaceFunctor = [&]( const typename SimplexRefinableElements< Dim >::FaceMultiIndex &fmi )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) if( fmi[d]==centerIndex ) return true;
		return false;
	};

	// Get the energy matrix
	Eigen::MatrixXd E = sre.template systemMatrix< EmbeddingDimension >( eType , positionFunctor );

	// Get the prolongation
	if constexpr( Dim<EmbeddingDimension )
	{
		if( planarityEpsilon>0 )
		{
			unsigned int count = 0;
			Point< double , EmbeddingDimension > c;
			SquareMatrix< double , EmbeddingDimension > cov;
			for( unsigned int i=0 ; i<simplexRefinableCell.size() ; i++ )
			{
				SimplexIndex< Dim , unsigned int > si = simplexRefinableCell[i];
				for( unsigned int d=0 ; d<=Dim ; d++ ) c += positionFunctor( si[d] ) , count++;
			}
			c /= (double) count;
			for( unsigned int i=0 ; i<simplexRefinableCell.size() ; i++ )
			{
				SimplexIndex< Dim , unsigned int > si = simplexRefinableCell[i];
				for( unsigned int d=0 ; d<=Dim ; d++ )
				{
					Point< double , EmbeddingDimension > p = positionFunctor( si[d] ) - c;
					for( unsigned int j=0 ; j<EmbeddingDimension ; j++ ) for( unsigned int k=0 ; k<EmbeddingDimension ; k++ ) cov(j,k) += p[j] * p[k];
				}
			}
			Polynomial::Polynomial< 1 , Dim , double > cPoly = cov.characteristicPolynomial();
			bool planar = true;
			for( unsigned int i=0 ; i<(EmbeddingDimension-Dim) ; i++ ) if( fabs( cPoly.coefficient(i) )>planarityEpsilon ) planar = false;
			if( planar ) return InterpolatingProlongationSystem( E , coarseIndices ).prolongation< false >();
			else
			{
				std::vector< Point< double , EmbeddingDimension > > interpolationConstraints( sre.size() );
				for( unsigned int i=0 ; i<sre.size() ; i++ )
				{
					const NodeMultiIndex &nmi = sre[i];
					Point< double , EmbeddingDimension > p;
					for( unsigned int d=0 ; d<Degree ; d++ ) p += positionFunctor( nmi[d] );
					interpolationConstraints[i] = p / Degree;
				}
				return InterpolatingProlongationSystem( E , coarseIndices , &interpolationConstraints[0] ).prolongation< true >();
			}
		}
		else
		{
			std::vector< Point< double , EmbeddingDimension > > interpolationConstraints( sre.size() );
			for( unsigned int i=0 ; i<sre.size() ; i++ )
			{
				const NodeMultiIndex &nmi = sre[i];
				Point< double , EmbeddingDimension > p;
				for( unsigned int d=0 ; d<Degree ; d++ ) p += positionFunctor( nmi[d] );
				interpolationConstraints[i] = p / Degree;
			}
			return InterpolatingProlongationSystem( E , coarseIndices , &interpolationConstraints[0] ).prolongation< true >();
		}
	}
	else return InterpolatingProlongationSystem( E , coarseIndices ).prolongation< false >();
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree , unsigned int EmbeddingDimension >
Eigen::MatrixXd InterpolatingProlongationSystem< PoU >::_BoundaryProlongation( const SimplexRefinableCell< Dim > &simplexRefinableCell , const SimplexRefinableElements< Dim , Degree > &sre , typename SimplexRefinableElements<>::EnergyWeights eWeights , std::function< Point< double , EmbeddingDimension > ( unsigned int ) > positionFunctor , double planarityEpsilon , const std::vector< unsigned int > &coarseIndices , const Eigen::MatrixXd &coarseP )
{
	typedef typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex NodeMultiIndex;
	if( coarseP.rows()!=coarseIndices.size() ) ERROR_OUT( "Dimensional mismatch: " , coarseP.rows() , " != "  , coarseIndices.size() );

	// Construct the energy description from the weights
	typename SimplexRefinableElements< Dim >::EnergyType eType(eWeights);
	unsigned int centerIndex = -1;
	if constexpr( Dim>1 ) centerIndex = simplexRefinableCell.centerIndex();
	eType.isIntegrationFaceFunctor = [&]( const typename SimplexRefinableElements< Dim >::FaceMultiIndex &fmi )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) if( fmi[d]==centerIndex ) return true;
		return false;
	};

	Eigen::MatrixXd fineE = sre.systemMatrix( eType , positionFunctor );

	// Compute the subset of element indices that are not in the coarse index set.
	// [NOTE] The set of indices in the coarse index set can be larger than the number of columns in coarseP
	// In particular, the coarse index set should contain _all_ Dim-1 nodes while coarseP will only index a subset.
	std::vector< unsigned int > fineIndices;
	{
		fineIndices.reserve( sre.size() - coarseIndices.size() );
		std::vector< bool > isCoarse( sre.size() , false );
		for( unsigned int i=0 ; i<coarseIndices.size() ; i++ ) isCoarse[ coarseIndices[i] ] = true;
		for( unsigned int i=0 ; i<isCoarse.size() ; i++ ) if( !isCoarse[i] ) fineIndices.push_back( i );
	}

	Eigen::MatrixXd _P = Eigen::MatrixXd::Zero( fineE.cols() , coarseP.cols() + fineIndices.size() );
	for( unsigned int c=0 ; c<coarseP.cols() ; c++ ) for( unsigned int r=0 ; r<coarseP.rows() ; r++ ) _P( coarseIndices[r] , c ) = coarseP( r , c );
	for( unsigned int i=0 ; i<fineIndices.size() ; i++ ) _P( fineIndices[i] , coarseP.cols()+i ) = 1.;

	// Get the energy matrix
	Eigen::MatrixXd E = _P.transpose() * fineE * _P;

	// Get the relative prolongation
	std::vector< unsigned int > _coarseIndices( coarseP.cols() );
	for( unsigned int i=0 ; i<_coarseIndices.size() ; i++ ) _coarseIndices[i] = i;
	_P = InterpolatingProlongationSystem( E , _coarseIndices ).prolongation< false >();

	// Extend the relative prolongation to the complete prolongation
	Eigen::MatrixXd P = Eigen::MatrixXd::Zero( sre.size() , _P.cols() );
	for( unsigned int c=0 ; c<_P.cols() ; c++ )
	{
		for( unsigned int r=0 ; r<coarseP.rows() ; r++ ) P( coarseIndices[r] , c ) = coarseP( r , c );
		for( unsigned int r=0 ; r<fineIndices.size() ; r++ ) P( fineIndices[r] , c ) = _P( coarseP.cols()+r , c );
	}

	return P;
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree >
void InterpolatingProlongationSystem< PoU >::_SetNodeMaps( const SimplexRefinableCell< Dim > &simplexRefinableCell , std::map< typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex , unsigned int > nodeMaps[Dim+1] )
{
	__SetNodeMaps< Dim , Degree >( simplexRefinableCell , nodeMaps );

	// Set the indices for (merged) lower resolutions nodes.
	// Note that for d=Dim, the indices are already set in correspondences with the node indices.
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		unsigned int count = 0;
		for( auto &[ nmi , idx ] : nodeMaps[d] ) idx = count++;
	}
}

template< bool PoU >
template< unsigned int Dim , unsigned int Degree >
void InterpolatingProlongationSystem< PoU >::__SetNodeMaps( const SimplexRefinableCell< Dim > &simplexRefinableCell , std::map< typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex , unsigned int > nodeMaps[Dim+1] )
{
	typedef typename SimplexRefinableElements< Dim , Degree >::NodeMultiIndex NodeMultiIndex;

	SimplexRefinableElements< Dim , Degree > sre( simplexRefinableCell );
	for( unsigned int i=0 ; i<sre.size() ; i++ ) nodeMaps[Dim][ sre[i] ] = i;
	if constexpr( Dim>1 ) for( unsigned int f=0 ; f<simplexRefinableCell.faces() ; f++ ) __SetNodeMaps< Dim-1 , Degree >( simplexRefinableCell.face(f) , nodeMaps );
	else if constexpr( Dim==1 )
	{
		auto VertexNodeIndex = []( unsigned int v )
		{
			unsigned int idx[Degree];
			for( unsigned int d=0 ; d<Degree ; d++ ) idx[d] = v;
			return NodeMultiIndex( idx );
		};

		// Start by getting the nodes on the boundary
		for( unsigned int i=0 ; i<2 ; i++ ) nodeMaps[0][ VertexNodeIndex( simplexRefinableCell[0][i] ) ] = i;
	}
}