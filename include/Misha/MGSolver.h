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

#ifndef MG_SOLVER_INCLUDED
#define MG_SOLVER_INCLUDED

#include <omp.h>
#include <list>
#include "Misha/PreProcess.h"
#include "Miscellany.h"

//#define USE_PARALLEL_GS_SORT



namespace MGSolver
{
	template< unsigned int Dim > using PointVector = Eigen::Matrix< Point< double , Dim > , Eigen::Dynamic , 1 >;

	template< unsigned int Dim >
	PointVector< Dim > operator + ( const PointVector< Dim > &v1 , const PointVector< Dim > &v2 )
	{
		assert( v1.rows()==v2.rows() );
		PointVector< Dim > w( v1.rows() );
		for( unsigned int i=0 ; i<v1.rows() ; i++ ) w[i] = v1[i] + v2[i];
		return w;
	}
	template< unsigned int Dim >
	PointVector< Dim > operator - ( const PointVector< Dim > &v1 , const PointVector< Dim > &v2 )
	{
		assert( v1.rows()==v2.rows() );
		PointVector< Dim > w( v1.rows() );
		for( unsigned int i=0 ; i<v1.rows() ; i++ ) w[i] = v1[i] - v2[i];
		return w;
	}

	template< unsigned int Dim >
	PointVector< Dim > operator * ( const PointVector< Dim > &v , double s )
	{
		PointVector< Dim > w( v.rows() );
		for( unsigned int i=0 ; i<v.rows() ; i++ ) w[i] = v[i] * s;
		return w;
	}

	template< unsigned int Dim >
	PointVector< Dim > operator * ( double s , const PointVector< Dim > &v ){ return v * s; }

	template< unsigned int Dim >
	PointVector< Dim > operator / ( const PointVector< Dim > &v , double s ){ return v * 1./s; }

	template< unsigned int Dim >
	PointVector< Dim > operator * ( const Eigen::SparseMatrix< double > &M , const PointVector< Dim > &v )
	{
		PointVector< Dim > w( M.rows() );
		for( unsigned int i=0 ; i<M.rows() ; i++ ) w[i] = Point< double , Dim >();
		for( unsigned int i=0 ; i<M.outerSize() ; i++ )
			for( Eigen::SparseMatrix< double >::InnerIterator iter(M,i) ; iter ; ++iter )
				w[ iter.row() ] += v[ iter.col() ] * iter.value();
		return w;
	}

	template< typename RelaxerType > struct Solver;

	enum
	{
		RELAXER_JACOBI ,
		RELAXER_GAUSS_SEIDEL ,
		RELAXER_PARALLEL_GAUSS_SEIDEL ,
		RELAXER_COUNT
	};
	static const std::string RelaxerTypeNames[] = { "Jacobi" , "Gauss-Seidel (serial)" , "GaussSeidel (parallel)" };

	template< typename V > struct VectorWrapper;
	template< typename V > struct ConstVectorWrapper;

	template< unsigned int D >
	struct VectorWrapper< std::vector< Point< double , D > > >
	{
		typedef std::vector< Point< double , D > > type;
		typedef Point< double , D > element;
		typedef Eigen::MatrixXd base_eigen_type;

		element &operator[]( unsigned int idx ) { return _values[idx]; }
		const element &operator[]( unsigned int idx ) const { return _values[idx]; }
		double squaredNorm( void ) const { double n2=0 ; for( unsigned int i=0 ; i<_values.size() ; i++ ) n2 += Point< double , D >::SquareNorm( _values[i] ) ; return n2; }
		size_t size( void ) const { return _values.size(); }
		unsigned int dim( void ) const { return D; }
		static unsigned int Dim( const type &t ){ return D; }
		void resize( size_t sz , unsigned int d=D ){ assert(d==D) ; _values.resize( sz ); }
		void setZero( void ){ for( unsigned int i=0 ; i<_values.size() ; i++ ) _values[i] = Point< double , D >(); }
		VectorWrapper &operator = ( const base_eigen_type &b )
		{
			assert( b.rows()==_values.size() && b.cols()==D );
			for( unsigned int i=0 ; i<b.rows() ; i++ ) for( unsigned int j=0 ; b.cols() ; j++ ) _values[i][j] = b(i,j);
			return *this;
		}
		operator base_eigen_type() const
		{
			base_eigen_type b(_values.size() , D );
			for( unsigned int i=0 ; i<_values.size() ; i++ ) for( unsigned int j=0 ; j<D ; j++ ) b(i,j) = _values[i][j];
			return b;
		}

		VectorWrapper( type &values ) : _values( values ){}
	protected:
		type &_values;
	};

	template<>
	struct VectorWrapper< Eigen::VectorXd >
	{
		typedef Eigen::VectorXd type;
		typedef double element;
		typedef Eigen::VectorXd base_eigen_type;

		element &operator[]( unsigned int idx ) { return _values[idx]; }
		const element &operator[]( unsigned int idx ) const { return _values[idx]; }
		double squaredNorm( void ) const { return _values.squaredNorm(); }
		size_t size( void ) const { return _values.size(); }
		unsigned int dim( void ) const { return 1; }
		static unsigned int Dim( const type &t ){ return 1; }
		void resize( size_t sz , unsigned int d=1 ){ assert(d==1) ; _values.resize( sz ); }
		void setZero( void ){ _values.setZero(); }
		VectorWrapper &operator = ( const base_eigen_type &b ){ _values = b ; return *this; }
		operator base_eigen_type() const { return _values; }

		VectorWrapper( type &values ) : _values( values ){}
	protected:
		type &_values;
	};

	template< unsigned int D >
	struct VectorWrapper< Eigen::Matrix< Point< double , D > , Eigen::Dynamic , 1 > >
	{
		typedef Eigen::Matrix< Point< double , D > , Eigen::Dynamic , 1 > type;
		typedef Point< double , D > element;
		typedef Eigen::MatrixXd base_eigen_type;

		element &operator[]( unsigned int idx ) { return _values[idx]; }
		const element &operator[]( unsigned int idx ) const { return _values[idx]; }
		double squaredNorm( void ) const { double n2=0 ; for( unsigned int i=0 ; i<_values.size() ; i++ ) n2 += Point< double , D >::SquareNorm( _values[i] ) ; return n2; }
		size_t size( void ) const { return _values.size(); }
		unsigned int dim( void ) const { return D; }
		static unsigned int Dim( const type &t ){ return D; }
		void resize( size_t sz , unsigned int d=D ){ assert(d==D) ; _values.resize( sz ); }
		void setZero( void ){ for( unsigned int i=0 ; i<_values.size() ; i++ ) _values[i] = Point< double , D >(); }
		VectorWrapper &operator = ( const base_eigen_type &b )
		{
			assert( b.rows()==_values.size() && b.cols()==D );
			for( unsigned int i=0 ; i<b.rows() ; i++ ) for( unsigned int j=0 ; j<b.cols() ; j++ ) _values[i][j] = b(i,j);
			return *this;
		}
		operator base_eigen_type() const
		{
			base_eigen_type b(_values.size() , D );
			for( unsigned int i=0 ; i<_values.size() ; i++ ) for( unsigned int j=0 ; j<D ; j++ ) b(i,j) = _values[i][j];
			return b;
		}

		VectorWrapper( type &values ) : _values( values ){}
	protected:
		type &_values;
	};

	template<>
	struct VectorWrapper< Eigen::MatrixXd >
	{
		typedef Eigen::MatrixXd type;
		typedef Eigen::Block< Eigen::MatrixXd , 1 , Eigen::Dynamic , false > element;
		typedef Eigen::MatrixXd base_eigen_type;

		element operator[]( unsigned int idx ) { return _values.row( (Eigen::Index)idx ); }
		const element operator[]( unsigned int idx ) const { return _values.row( idx ); }
		double squaredNorm( void ) const { return _values.squaredNorm(); }
		size_t size( void ) const { return _values.rows(); }
		unsigned int dim( void ) const { return (unsigned int)_values.cols(); }
		static unsigned int Dim( const type &t ){ return (unsigned int)t.cols(); }
		void resize( size_t sz , unsigned int d ){ _values.resize( sz , d ); }
		void setZero( void ){ _values.setZero(); }
		VectorWrapper &operator = ( const base_eigen_type &b ){ _values = b ; return *this; }
		operator base_eigen_type() const { return _values; }

		VectorWrapper( type &values ) : _values( values ){}
	protected:
		type &_values;
	};

	template< unsigned int D >
	struct ConstVectorWrapper< std::vector< Point< double , D > > >
	{
		typedef std::vector< Point< double , D > > type;
		typedef Point< double , D > element;
		typedef Eigen::MatrixXd base_eigen_type;

		const element &operator[]( unsigned int idx ) const { return _values[idx]; }
		double squaredNorm( void ) const { double n2=0 ; for( unsigned int i=0 ; i<_values.size() ; i++ ) n2 += Point< double , D >::SquareNorm( _values[i] ) ; return n2; }
		size_t size( void ) const { return _values.size(); }
		unsigned int dim( void ) const { return D; }
		static unsigned int Dim( const type &t ){ return D; }
		operator base_eigen_type() const
		{
			base_eigen_type b(_values.size() , D );
			for( unsigned int i=0 ; i<_values.size() ; i++ ) for( unsigned int j=0 ; j<D ; j++ ) b(i,j) = _values[i][j];
			return b;
		}

		ConstVectorWrapper( const type &values ) : _values( values ){}
	protected:
		const type &_values;
	};

	template<>
	struct ConstVectorWrapper< Eigen::VectorXd >
	{
		typedef Eigen::VectorXd type;
		typedef double element;
		typedef Eigen::VectorXd base_eigen_type;

		const element &operator[]( unsigned int idx ) const { return _values[idx]; }
		double squaredNorm( void ) const { return _values.squaredNorm(); }
		size_t size( void ) const { return _values.size(); }
		unsigned int dim( void ) const { return 1; }
		static unsigned int Dim( const type &t ){ return 1; }
		operator base_eigen_type() const { return _values; }

		ConstVectorWrapper( const type &values ) : _values( values ){}
	protected:
		const type &_values;
	};

	template< unsigned int D >
	struct ConstVectorWrapper< Eigen::Matrix< Point< double , D > , Eigen::Dynamic , 1 > >
	{
		typedef Eigen::Matrix< Point< double , D > , Eigen::Dynamic , 1 > type;
		typedef Point< double , D > element;
		typedef Eigen::MatrixXd base_eigen_type;

		const element &operator[]( unsigned int idx ) const { return _values[idx]; }
		double squaredNorm( void ) const { double n2=0 ; for( unsigned int i=0 ; i<_values.size() ; i++ ) n2 += Point< double , D >::SquareNorm( _values[i] ) ; return n2; }
		size_t size( void ) const { return _values.size(); }
		unsigned int dim( void ) const { return D; }
		static unsigned int Dim( const type &t ){ return D; }
		operator base_eigen_type() const
		{
			base_eigen_type b(_values.size() , D );
			for( unsigned int i=0 ; i<_values.size() ; i++ ) for( unsigned int j=0 ; j<D ; j++ ) b(i,j) = _values[i][j];
			return b;
		}

		ConstVectorWrapper( const type &values ) : _values( values ){}
	protected:
		const type &_values;
	};

	template<>
	struct ConstVectorWrapper< Eigen::MatrixXd >
	{
		typedef Eigen::MatrixXd type;
		typedef Eigen::Block< const Eigen::MatrixXd , 1 , Eigen::Dynamic , false > element;
		typedef Eigen::MatrixXd base_eigen_type;

		const element operator[]( unsigned int idx ) const { return _values.row( idx ); }
		double squaredNorm( void ) const { return _values.squaredNorm(); }
		size_t size( void ) const { return _values.rows(); }
		unsigned int dim( void ) const { return (unsigned int)_values.cols(); }
		static unsigned int Dim( const type &t ){ return (unsigned int)t.cols(); }
		operator base_eigen_type() const { return _values; }

		ConstVectorWrapper( const type &values ) : _values( values ){}
	protected:
		const type &_values;
	};

	struct RowRelaxer
	{
		// Jacobi / Gauss-Seidel relaxation:
		// For each i, we want:
		//		\sum_j S(i,j) * x[j] = b[i]
		// => S(i,i) * x[i] + \sum_j S(i,j) * x[j] - S(i,i) * x[i] - b[i] = 0
		// => x[i]  = ( b[i] - \sum_j S(i,j) * x[j] + S(i,i) * x[i] ) / S(i,i)
		// => x[i]  = ( b[i] - \sum_j S(i,j) * x[j] ) / S(i,i) + x[i]
		// => x[i] += ( b[i] - \sum_j S(i,j) * x[j] ) / S(i,i)

		RowRelaxer( void ){}
		RowRelaxer( RowRelaxer &&relaxer ){ std::swap( _diagonalR , relaxer._diagonalR ); }
		RowRelaxer &operator = ( RowRelaxer &&relaxer ){ std::swap( _diagonalR , relaxer._diagonalR ) ; return *this; }

		void setDiagonalR( const Eigen::SparseMatrix< double > &M );
	protected:
		std::vector< double > _diagonalR;
	};

	struct JacobiRelaxer : public RowRelaxer
	{
		JacobiRelaxer( const Eigen::SparseMatrix< double > &M , bool verbose=false ){ setDiagonalR( M ); }
		template< typename V >
		void relax( const Eigen::SparseMatrix< double > &M , const V &b , V &x , unsigned int iters , bool forward );
	protected:
		JacobiRelaxer( void ){}
		friend struct Solver< JacobiRelaxer >;
	};

	struct GaussSeidelRelaxer : public RowRelaxer
	{
		GaussSeidelRelaxer( const Eigen::SparseMatrix< double > &M , bool verbose=false ){ setDiagonalR( M ); }
		template< typename V >
		void relax( const Eigen::SparseMatrix< double > &M , const V &b , V &x , unsigned int iters , bool forward );
	protected:
		GaussSeidelRelaxer( void ){}
		friend struct Solver< GaussSeidelRelaxer >;
	};

	template< unsigned int ParallelizationThreshold >
	struct ParallelGaussSeidelRelaxer : public RowRelaxer
	{
#ifdef USE_PARALLEL_GS_SORT
		enum
		{
			SORT_NONE ,
			SORT_FLOOD_FILL ,
			SORT_MIN_DEGREE ,
			SORT_MAX_DEGREE ,
			SORT_COUNT
		};
		const static std::string SortTypeNames[];
		ParallelGaussSeidelRelaxer( const Eigen::SparseMatrix< double > &M , unsigned int sortType=SORT_MAX_DEGREE ){ setDiagonalR( M ) ; setMultiColorIndices( M , sortType ); }
#else // !USE_PARALLEL_GS_SORT
		ParallelGaussSeidelRelaxer( const Eigen::SparseMatrix< double > &M , bool verbose=false ){ setDiagonalR( M ) ; setMultiColorIndices( M , verbose ); }
#endif // USE_PARALLEL_GS_SORT
		template< typename V >
		void relax( const Eigen::SparseMatrix< double > &M , const V &b , V &x , unsigned int iters , bool forward );
#ifdef USE_PARALLEL_GS_SORT
		void setMultiColorIndices( const Eigen::SparseMatrix< double > &M , unsigned int sortType );
#else // !USE_PARALLEL_GS_SORT
		void setMultiColorIndices( const Eigen::SparseMatrix< double > &M , bool verbose );
#endif // USE_PARALLEL_GS_SORT
	protected:
		ParallelGaussSeidelRelaxer( void ){}

		std::vector< std::vector< unsigned int > > _multiColorIndices;
		std::vector< unsigned int > _serialIndices;

		friend struct Solver< ParallelGaussSeidelRelaxer< ParallelizationThreshold > >;
	};

	template< typename RelaxerType >
	struct Solver
	{
		struct State
		{
			unsigned int vCycles;
			unsigned int rIters , pIters;
			bool verbose;
			State( void ) : vCycles(1) , rIters(0) , pIters(0) , verbose(false) {}
		};

		State state;

		~Solver( void ){ delete[] _relaxers; }
		Solver( const Eigen::SparseMatrix< double > &S , const std::vector< Eigen::SparseMatrix< double > > &P , bool verbose );
		template< typename V > V solve( const V &b ){ return solve( b , state.vCycles , state.rIters , state.pIters , state.verbose ); }
		template< typename V > V solve( const V &b , unsigned int vCycles , unsigned int rIters , unsigned int pIters , bool verbose );
		template< typename V > V solve( const V &b , const V &x ){ return solve( b , x , state.vCycles , state.rIters , state.pIters , state.verbose ); }
		template< typename V > V solve( const V &b , const V &x , unsigned int vCycles , unsigned int rIters , unsigned int pIters , bool verbose );
	protected:
		std::vector< Eigen::SparseMatrix< double > > _S , _P , _Pt;
		SparseSolver::LLT _solver;
		RelaxerType *_relaxers;

		template< typename V >
		void _solve( std::vector< V > &b , std::vector< V > &x , unsigned int vCycles , unsigned int rIters , unsigned int pIters , bool verbose );
	};

	typedef Solver< JacobiRelaxer > Jacobi;
	typedef Solver< GaussSeidelRelaxer > GaussSeidel;
	template< unsigned int ParallelizationThreshold >
	using ParallelGaussSeidel = Solver< ParallelGaussSeidelRelaxer< ParallelizationThreshold > >;

	////////////////
	// RowRelaxer //
	////////////////
	void RowRelaxer::setDiagonalR( const Eigen::SparseMatrix< double > &M )
	{
		_diagonalR.resize( M.rows() , 0 );
		for( unsigned int j=0 ; j<M.outerSize() ; j++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(M,j) ; iter ; ++iter )
			if( iter.row()==iter.col() ) _diagonalR[ iter.row() ] += iter.value();
		for( unsigned int j=0 ; j<_diagonalR.size() ; j++ ) _diagonalR[j] = 1./_diagonalR[j];
	}

	///////////////////
	// JacobiRelaxer //
	///////////////////
	template< typename V >
	void JacobiRelaxer::relax( const Eigen::SparseMatrix< double > &M , const V &b , V &x , unsigned int iters , bool forward )
	{
		assert( b.rows()==x.rows() && b.cols()==x.cols() && M.cols()==b.rows() );
		Eigen::MatrixXd _x( b.rows() , b.cols() );

		auto ProcessIndex = [&]( unsigned int i )
		{
			Eigen::VectorXd delta = b.row(i);
			for( Eigen::SparseMatrix< double >::InnerIterator iter(M,i) ; iter ; ++iter ) delta -= iter.value() * x.row( iter.row() );
			_x.row(i) = x.row(i) + delta * _diagonalR[i];
		};

		for( unsigned int iter=0 ; iter<iters ; iter++ )
		{
#pragma omp parallel for
			for( int i=0 ; i<(int)b.size() ; i++ ) ProcessIndex(i);
#pragma omp parallel for
			for( int i=0 ; i<(int)x.size() ; i++ ) x.row(i) = _x.row(i);
		}
	}

	////////////////////////
	// GaussSeidelRelaxer //
	////////////////////////
	template< typename V >
	void GaussSeidelRelaxer::relax( const Eigen::SparseMatrix< double > &M , const V &b , V &x , unsigned int iters , bool forward )
	{
		auto ProcessIndex = [&]( unsigned int i )
		{
			Eigen::VectorXd delta = b.row(i);
			for( Eigen::SparseMatrix< double >::InnerIterator iter(M,i) ; iter ; ++iter ) delta -= iter.value() * x.row( iter.row() );
			x.row(i) += delta * _diagonalR[i];
		};

		for( unsigned int iter=0 ; iter<iters ; iter++ )
			if( forward ) for( unsigned int i=0 ; i<b.size() ; i++ ) ProcessIndex(i);
			else          for( int i=(int)b.size()-1 ; i>=0  ; i-- ) ProcessIndex(i);
	}

	////////////////////////////////
	// ParallelGaussSeidelRelaxer //
	////////////////////////////////
#ifdef USE_PARALLEL_GS_SORT
	template< unsigned int ParallelizationThreshold >
	const std::string ParallelGaussSeidelRelaxer< ParallelizationThreshold >::SortTypeNames[] = { "none" , "flood fill" , "min degree" , "max degree" };
#endif // USE_PARALLEL_GS_SORT
	template< unsigned int ParallelizationThreshold >
	template< typename V >
	void ParallelGaussSeidelRelaxer< ParallelizationThreshold >::relax( const Eigen::SparseMatrix< double > &M , const V &b , V &x , unsigned int iters , bool forward )
	{
		VectorWrapper< V > _x(x);
		ConstVectorWrapper< V > _b(b);
		auto ProcessIndex = [&]( unsigned int i )
		{
			typename VectorWrapper< V >::element delta = _b[i];
			for( Eigen::SparseMatrix< double >::InnerIterator iter(M,i) ; iter ; ++iter ) delta -= iter.value() * _x[ (unsigned int)iter.row() ];
			_x[i] += delta * _diagonalR[i];
		};

		auto ProcessColor = [&]( const std::vector< unsigned int > &color )
		{
#pragma omp parallel for
			for( int i=0 ; i<(int)color.size() ; i++ ) ProcessIndex( color[i] );
		};

		for( unsigned int iter=0 ; iter<iters ; iter++ )
		{
			if( !forward ) std::for_each( _serialIndices.rbegin() , _serialIndices.rend() , [&]( unsigned int i ){ ProcessIndex(i); } );
			if( forward ) std::for_each( _multiColorIndices. begin() , _multiColorIndices. end() , [&]( const std::vector< unsigned int > &clr ){ ProcessColor( clr ); } );
			else          std::for_each( _multiColorIndices.rbegin() , _multiColorIndices.rend() , [&]( const std::vector< unsigned int > &clr ){ ProcessColor( clr ); } );
			if(  forward ) std::for_each( _serialIndices. begin() , _serialIndices. end() , [&]( unsigned int i ){ ProcessIndex(i); } );
		}
	}

	template< unsigned int ParallelizationThreshold >
#ifdef USE_PARALLEL_GS_SORT
	void ParallelGaussSeidelRelaxer< ParallelizationThreshold >::setMultiColorIndices( const Eigen::SparseMatrix< double > &M , unsigned int sortType )
#else // !USE_PARALLEL_GS_SORT
	void ParallelGaussSeidelRelaxer< ParallelizationThreshold >::setMultiColorIndices( const Eigen::SparseMatrix< double > &M , bool verbose )
#endif // USE_PARALLEL_GS_SORT
	{
		Miscellany::Timer timer;
#ifdef USE_PARALLEL_GS_SORT
		std::vector< unsigned int > indices( M.rows() );

		if( sortType==SORT_NONE )
		{
			for( unsigned int i=0 ; i<M.rows() ; i++ ) indices[i] = i;
		}
		else if( sortType==SORT_MIN_DEGREE || sortType==SORT_MAX_DEGREE ) 
		{
			std::vector< unsigned int > rowSizes( M.rows() , 0 );
			for( unsigned int i=0 ; i<M.rows() ; i++ )
			{
				indices[i] = i;
				for( Eigen::SparseMatrix< double >::InnerIterator iter(M,i) ; iter ; ++iter ) rowSizes[i]++;
			}
			if     ( sortType==SORT_MIN_DEGREE ) std::sort( indices.begin() , indices.end() , [&]( unsigned int i , unsigned int j ){ return rowSizes[i]<rowSizes[j]; } );
			else if( sortType==SORT_MAX_DEGREE ) std::sort( indices.begin() , indices.end() , [&]( unsigned int i , unsigned int j ){ return rowSizes[i]>rowSizes[j]; } );
		}
		else if( sortType==SORT_FLOOD_FILL )
		{
			unsigned int idx = 0;
			std::vector< bool > added( M.rows() , false );
			added[idx] = true;
			indices[idx++] = 0;
			for( unsigned int i=0 ; i<M.rows() ; i++ )
			{
				for( Eigen::SparseMatrix< double >::InnerIterator iter(M,indices[i]) ; iter ; ++iter ) if( !added[ iter.row() ] )
				{
					added[ iter.row() ] = true;
					indices[idx++] = (unsigned int)iter.row();
				}
			}
		}
#endif // USE_PARALLEL_GS_SORT

		// Repeatedly:
		//		Iterate through unmarked and unprocessed nodes
		//		Add each such node to the next color groups
		//		Mark neighbors

		// Linked-list tracking the set of indices that have not been assigned a color
		std::list< unsigned int > indexList;

		// Flags tracking whether an index is "independent"/"dependent"
		std::vector< bool > dependentFlags( M.rows() , false );

		// Construct the initial (total) linked list
#ifdef USE_PARALLEL_GS_SORT
		for( unsigned int i=0 ; i<M.rows() ; i++ ) indexList.push_back( indices[i] );
#else // !USE_PARALLEL_GS_SORT
		for( unsigned int i=0 ; i<M.rows() ; i++ ) indexList.insert( indexList.begin() , i );
#endif //  USE_PARALLEL_GS_SORT

		// The next color group
		std::vector< unsigned int > colorIndices;
		colorIndices.reserve( M.rows() );

		// Iterate until each index has been assigned a color
		while( !indexList.empty() )
		{
			// Set the states of all uncolored indices to "independent"
			for( auto iter=indexList.begin() ; iter!=indexList.end() ; iter++ ) dependentFlags[*iter] = false;

			// Iterate over all "independent" indices
			for( auto iter=indexList.begin() ; iter!=indexList.end() ; ) if( !dependentFlags[*iter] )
			{
				unsigned int i = *iter;

				// Add the index to the current color group
				colorIndices.push_back( i );

				// Mark neighbors as "dependent"
				for( Eigen::SparseMatrix< double >::InnerIterator iter(M,i) ; iter ; ++iter ) dependentFlags[ iter.row() ] = true;

				// Remove the index from the linked list
				iter = indexList.erase( iter );
			}
			else iter++;

			// Add the color group to the list of color groups
			_multiColorIndices.push_back( colorIndices );
			colorIndices.clear();
		}

		// Merge the small color groups into a single color group that will be processed serially
		{
			// The count indices to be processed serially
			unsigned int serialCount = 0;
			for( unsigned int i=0 ; i<_multiColorIndices.size() ; i++ ) if( _multiColorIndices[i].size()<ParallelizationThreshold ) serialCount += (unsigned int)_multiColorIndices[i].size();

			// Generate the serial group
			_serialIndices.reserve( serialCount );
			for( int i=(int)_multiColorIndices.size()-1 ; i>=0 ; i-- ) if( _multiColorIndices[i].size()<ParallelizationThreshold )
			{
				for( unsigned int j=0 ; j<_multiColorIndices[i].size() ; j++ ) _serialIndices.push_back( _multiColorIndices[i][j] );
				std::swap( _multiColorIndices[i] , _multiColorIndices.back() );
				_multiColorIndices.pop_back();
			}
		}
		if( verbose ) std::cout << "Got " << _multiColorIndices.size() << " colors in: " << timer.elapsed() << "(s)" << std::endl;
	}

	////////////
	// Solver //
	////////////
	template< typename RelaxerType >
	Solver< RelaxerType >::Solver( const Eigen::SparseMatrix< double > &S , const std::vector< Eigen::SparseMatrix< double > > &P , bool verbose ) : _P(P)
	{
		_Pt.resize( _P.size() );
		for( unsigned int i=0 ; i<_P.size() ; i++ ) _Pt[i] = _P[i].transpose();
		_S.resize( _P.size()+1 );
		_S.back() = S;
		for( int i=(int)_P.size()-1 ; i>=0 ; i-- ) _S[i] = _Pt[i] * _S[i+1] * _P[i];
		_solver.compute( _S[0] );
		switch( _solver.info() )
		{
			case Eigen::NumericalIssue: ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- numerical issue" );
			case Eigen::NoConvergence:  ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- no convergence" );
			case Eigen::InvalidInput:   ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- invalid input" );
			case Eigen::Success: ;
		}
		_relaxers = new RelaxerType[ _S.size() ];

		for( unsigned int i=0 ; i<_S.size() ; i++ ) _relaxers[i] = RelaxerType( _S[i] , verbose );

		if( verbose )
		{
			for( unsigned int l=0 ; l<_S.size() ; l++ )
			{
				std::cout << "Level: " << l << std::endl;
				std::cout << "\tDoFs / Non-zero matrix entries / Entries per row: " << _S[l].rows() << " / " << _S[l].nonZeros() << " / " << _S[l].nonZeros()/_S[l].rows() << std::endl;
			}
		}
	}

	template< typename RelaxerType >
	template< typename V >
	V  Solver< RelaxerType >::solve( const V &b , unsigned int vCycles , unsigned int rIters , unsigned int pIters , bool verbose )
	{
		std::vector< V > _x , _b;
		_x.resize( _S.size() );
		_b.resize( _S.size() );

		for( unsigned int i=0 ; i<_S.size() ; i++ )
		{
			VectorWrapper< V > __x(_x[i]) , __b(_b[i]);
			__x.resize( _S[i].rows() , VectorWrapper< V >::Dim( b ) );
			__b.resize( _S[i].rows() , VectorWrapper< V >::Dim( b ) );
		}

		_b.back() = b;
		VectorWrapper< V >( _x.back() ).setZero();
		_solve( _b , _x , vCycles , rIters , pIters , verbose );
		return _x.back();
	}

	template< typename RelaxerType >
	template< typename V >
	V Solver< RelaxerType >::solve( const V &b , const V &x , unsigned int vCycles , unsigned int rIters , unsigned int pIters , bool verbose )
	{
		assert( b.rows()==x.rows() && b.cols()==x.cols() );
		std::vector< V > _x , _b;
		_x.resize( _S.size() );
		_b.resize( _S.size() );

		for( unsigned int i=0 ; i<_S.size() ; i++ )
		{
			VectorWrapper< V > __x(_x[i]) , __b(_b[i]);
			__x.resize( _S[i].rows() , VectorWrapper< V >::Dim( b ) );
			__b.resize( _S[i].rows() , VectorWrapper< V >::Dim( b ) );
		}

		_b.back() = b;
		_x.back() = x;
		_solve< V >( _b , _x , vCycles , rIters , pIters , verbose );
		return _x.back();
	}

	template< typename RelaxerType >
	template< typename V >
	void Solver< RelaxerType >::_solve( std::vector< V > &_b , std::vector< V > &_x , unsigned int vCycles , unsigned int rIters , unsigned int pIters , bool verbose )
	{
		auto PrintError = [&]( unsigned int d )
		{
			double b1 = 0 , b2 = 0;
			V r = _b[d] - _S[d] * _x[d];
			b1 = VectorWrapper< V >( _b[d] ).squaredNorm();
			b2 = VectorWrapper< V >( r ).squaredNorm();
			for( unsigned int i=0 ; i<=d ; i++ ) std::cout << "  ";
			std::cout << b1 << " -> " << b2 << std::endl;
		};

		for( unsigned int v=0 ; v<vCycles ; v++ )
		{
			// The restriction phase
			for( unsigned int d=(unsigned int)_S.size()-1 ; d!=0 ; d-- )
			{
				// Do the relaxation
				_relaxers[d].template relax< typename VectorWrapper< V >::type >( _S[d] , _b[d] , _x[d] , rIters , true );
				if( verbose ) PrintError(d);
				// Compute the residual and restrict
				_b[d-1] = _Pt[d-1] * ( _b[d] - _S[d] * _x[d] );
				VectorWrapper< V >( _x[d-1] ).setZero();
			}

			// The base solve
			typename VectorWrapper< V >::base_eigen_type __b = VectorWrapper< V >( _b[0] );
			typename VectorWrapper< V >::base_eigen_type __x = _solver.solve( __b );
			VectorWrapper< V > _wx( _x[0] );
			_wx = __x;
			if( verbose ) PrintError(0);

			// The prolongation phase
			for( unsigned int d=1 ; d<_S.size() ; d++ )
			{
				// Prolong the coarse solution
				_x[d] += _P[d-1] * _x[d-1];

				// Do the relaxation
				_relaxers[d].template relax< typename VectorWrapper< V >::type >( _S[d] , _b[d] , _x[d] , pIters , false );
				if( verbose ) PrintError(d);
			}
		}
	}
}
#endif // MG_SOLVER_INCLUDED