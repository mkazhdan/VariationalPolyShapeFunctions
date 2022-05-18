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

template< unsigned int BlockSize >
Eigen::SparseMatrix< double > BlockExpand( const Eigen::SparseMatrix< double > &A )
{
	Eigen::SparseMatrix< double > _A( A.rows()*BlockSize , A.cols()*BlockSize );

	std::vector< Eigen::Triplet< double > > entries;
	entries.reserve( A.nonZeros()*BlockSize );
	for( unsigned int i=0 ; i<A.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(A,i) ; iter ; ++iter )
		for( unsigned int b=0 ; b<BlockSize ; b++ )
			entries.push_back( Eigen::Triplet< double >( (unsigned int)(iter.row()*BlockSize+b) , (unsigned int)(iter.col()*BlockSize+b) , iter.value() ) );

	_A.setFromTriplets( entries.begin() , entries.end() );
	return _A;
}


//////////////////////
// SolidSimplexMesh //
//////////////////////
template< unsigned int Dim , unsigned int Degree >
template< typename Index >
SolidSimplexMesh< Dim , Degree > SolidSimplexMesh< Dim , Degree >::Init( const std::vector< SimplexIndex< Dim , Index > > &simplices , std::function< Point< double , Dim > ( Index ) > vFunction )
{
	SolidSimplexMesh sm;
	sm._init( simplices , vFunction );

	sm._A.resize( sm._simplices.size() );
	for( unsigned int s=0 ; s<sm._simplices.size() ; s++ ) for( unsigned int d=1 ; d<=Dim ; d++ )
	{
		Point< double , Dim > v = Point< double , Dim >( vFunction( sm._simplices[s][d] ) - vFunction( sm._simplices[s][0] ) );
		for( unsigned int _d=0 ; _d<Dim ; _d++ ) sm._A[s](d-1,_d) = v[_d];
	}

	return sm;
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SolidSimplexMesh< Dim , Degree >::massMatrix( void ) const
{

	std::vector< Eigen::Triplet< double > > entries( _simplices.size() * NodesPerSimplex * NodesPerSimplex * Dim * Dim );
#pragma omp parallel for
	for( int s=0 ; s<(int)_simplices.size() ; s++ )
	{
		unsigned int indices[ SimplexElements< Dim , Degree >::NodeNum*Dim ];
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum ; i++ )
		{
			unsigned int idx = nodeIndex( s , i );
			for( unsigned int d=0 ; d<Dim ; d++ ) indices[i*Dim+d] = idx*Dim+d;
		}

		SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > M = SolidSimplexElements< Dim , Degree >::MassMatrix( _A[s] );
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum*Dim ; i++ )
			for( unsigned int j=0 ; j<SimplexElements< Dim , Degree >::NodeNum*Dim ; j++ )
					entries[ s*NodesPerSimplex*NodesPerSimplex*Dim*Dim + i * NodesPerSimplex*Dim + j ] = Eigen::Triplet< double >( indices[i] , indices[j] , M( i , j ) );
	}
	Eigen::SparseMatrix< double > M( nodes()*Dim , nodes()*Dim );
	M.setFromTriplets( entries.begin() , entries.end() );
	return M;
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SolidSimplexMesh< Dim , Degree >::frobeniusStiffnessMatrix( void ) const
{
	std::vector< Eigen::Triplet< double > > entries( _simplices.size() * NodesPerSimplex * NodesPerSimplex * Dim * Dim );
#pragma omp parallel for
	for( int s=0 ; s<(int)_simplices.size() ; s++ )
	{
		unsigned int indices[ SimplexElements< Dim , Degree >::NodeNum*Dim ];
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum ; i++ )
		{
			unsigned int idx = nodeIndex( s , i );
			for( unsigned int d=0 ; d<Dim ; d++ ) indices[i*Dim+d] = idx*Dim+d;
		}

		SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > S = SolidSimplexElements< Dim , Degree >::FrobeniusStiffnessMatrix( _A[s] );
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum*Dim ; i++ )
			for( unsigned int j=0 ; j<SimplexElements< Dim , Degree >::NodeNum*Dim ; j++ )
					entries[ s*NodesPerSimplex*NodesPerSimplex*Dim*Dim + i * NodesPerSimplex*Dim + j ] = Eigen::Triplet< double >( indices[i] , indices[j] , S( i , j ) );
	}
	Eigen::SparseMatrix< double > S( nodes()*Dim , nodes()*Dim );
	S.setFromTriplets( entries.begin() , entries.end() );
	return S;
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SolidSimplexMesh< Dim , Degree >::traceStiffnessMatrix( void ) const
{
	std::vector< Eigen::Triplet< double > > entries( _simplices.size() * NodesPerSimplex * NodesPerSimplex * Dim * Dim );
#pragma omp parallel for
	for( int s=0 ; s<(int)_simplices.size() ; s++ )
	{
		unsigned int indices[ SimplexElements< Dim , Degree >::NodeNum*Dim ];
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum ; i++ )
		{
			unsigned int idx = nodeIndex( s , i );
			for( unsigned int d=0 ; d<Dim ; d++ ) indices[i*Dim+d] = idx*Dim+d;
		}

		SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > S = SolidSimplexElements< Dim , Degree >::TraceStiffnessMatrix( _A[s] );
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum*Dim ; i++ )
			for( unsigned int j=0 ; j<SimplexElements< Dim , Degree >::NodeNum*Dim ; j++ )
				entries[ s*NodesPerSimplex*NodesPerSimplex*Dim*Dim + i * NodesPerSimplex*Dim + j ] = Eigen::Triplet< double >( indices[i] , indices[j] , S( i , j ) );
	}
	Eigen::SparseMatrix< double > S( nodes()*Dim , nodes()*Dim );
	S.setFromTriplets( entries.begin() , entries.end() );
	return S;
}

template< unsigned int Dim , unsigned int Degree >
void SolidSimplexMesh< Dim , Degree >::setMassFrobeniusStiffnessAndTraceStiffnessMatrices( Eigen::SparseMatrix< double > &M , Eigen::SparseMatrix< double > &F , Eigen::SparseMatrix< double > & T ) const
{
	std::vector< Eigen::Triplet< double > > entriesM( _simplices.size() * NodesPerSimplex * NodesPerSimplex * Dim * Dim );
	std::vector< Eigen::Triplet< double > > entriesF( _simplices.size() * NodesPerSimplex * NodesPerSimplex * Dim * Dim );
	std::vector< Eigen::Triplet< double > > entriesT( _simplices.size() * NodesPerSimplex * NodesPerSimplex * Dim * Dim );
#pragma omp parallel for
	for( int s=0 ; s<(int)_simplices.size() ; s++ )
	{
		unsigned int indices[ NodesPerSimplex*Dim ];
		for( unsigned int i=0 ; i<NodesPerSimplex ; i++ )
		{
			unsigned int idx = nodeIndex( s , i );
			for( unsigned int d=0 ; d<Dim ; d++ ) indices[i*Dim+d] = idx*Dim+d;
		}

		SquareMatrix< double , NodesPerSimplex*Dim > M , F , T;
		SolidSimplexElements< Dim , Degree >::SetMassFrobeniusStiffnessAndTraceStiffnessMatrices( _A[s] , M , F , T );

		unsigned int idx = s*NodesPerSimplex*NodesPerSimplex*Dim*Dim;
		for( unsigned int i=0 ; i<NodesPerSimplex*Dim ; i++ ) for( unsigned int j=0 ; j<NodesPerSimplex*Dim ; j++ , idx++ )
		{
			entriesM[ idx ] = Eigen::Triplet< double >( indices[i] , indices[j] , M(i,j) );
			entriesF[ idx ] = Eigen::Triplet< double >( indices[i] , indices[j] , F(i,j) );
			entriesT[ idx ] = Eigen::Triplet< double >( indices[i] , indices[j] , T(i,j) );
		}
	}
	M.resize( nodes()*Dim , nodes()*Dim );
	F.resize( nodes()*Dim , nodes()*Dim );
	T.resize( nodes()*Dim , nodes()*Dim );
	M.setFromTriplets( entriesM.begin() , entriesM.end() );
	F.setFromTriplets( entriesF.begin() , entriesF.end() );
	T.setFromTriplets( entriesT.begin() , entriesT.end() );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::VectorXd SolidSimplexMesh< Dim , Degree >::stiffnessVector( void ) const
{
	Eigen::VectorXd S = Eigen::VectorXd::Zero ( nodes()*Dim );

	for( unsigned int s=0 ; s<_simplices.size() ; s++ )
	{
		Point< double , SimplexElements< Dim , Degree >::NodeNum * Dim > _stiffness = SolidSimplexElements< Dim , Degree >::StiffnessVector( _A[s] );
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ )
			S[ nodeIndex( s , i )*Dim + d ] += _stiffness[ i*Dim + d ];
	}
	return S;
}

template< unsigned int Dim , unsigned int Degree >
Eigen::VectorXd SolidSimplexMesh< Dim , Degree >::dcVector( Point< double , Dim > v ) const
{
	Eigen::VectorXd DC = Eigen::VectorXd::Zero ( nodes()*Dim );

	for( unsigned int s=0 ; s<_simplices.size() ; s++ )
	{
		Point< double , SimplexElements< Dim , Degree >::NodeNum * Dim > _dc = SolidSimplexElements< Dim , Degree >::DC( _A[s] ) * v;
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ )
			DC[ nodeIndex( s , i )*Dim + d ] += _dc[ i*Dim + d ];
	}
	return DC;
}


template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SolidSimplexMesh< Dim , Degree >::evaluationMatrix( const std::vector< typename SimplexMesh< Dim >::Sample > &samples ) const
{
	return BlockExpand< Dim >( SimplexMesh< Dim , Degree >::evaluationMatrix( samples ) );
}