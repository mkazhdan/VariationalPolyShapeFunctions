/*
Copyright (c) 2020, Michael Kazhdan
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

#ifndef TENSORS_INCLUDED
#define TENSORS_INCLUDED

#include <iostream>
#include "Misha/Algebra.h"
#include "Windows.h"

namespace AutoDiff
{
	template< typename Pack > struct _Tensor;

	// A zero-tensor is the same as a double value
	template<>
	struct _Tensor< UIntPack<> > : public InnerProductSpace< double , _Tensor< UIntPack<> > >
	{
		typedef UIntPack<> Pack;
		static const unsigned int Size = 1;

		double data;

		_Tensor( double d=0 ) : data(d) {}
		explicit operator double &( void ){ return data; }
		explicit operator const double &( void ) const { return data; }

		void Add( const _Tensor &t ){ data += t.data; }
		void Scale( double s ){ data *= s; }
		double InnerProduct( const _Tensor &t ) const { return data * t.data; }
		template< unsigned int ... _Dims >
		_Tensor< UIntPack< _Dims ... > > operator * ( const _Tensor< UIntPack< _Dims ... > > &t ) const { return t * data; }
		template< unsigned int I , unsigned int ... _Dims >
		_Tensor< UIntPack< _Dims ... > > contractedOuterProduct( const _Tensor< UIntPack< _Dims ... > > &t ) const 
		{
			static_assert( I==0 , "[ERROR] Contraction suffix/prefix don't match" );
			return *this * t;
		}

		// Permute indices
		template< unsigned int ... PermutationValues >
		_Tensor< Permutation< Pack , UIntPack< PermutationValues ... > > > permute( UIntPack< PermutationValues ... > ) const
		{
			static_assert( sizeof ... ( PermutationValues ) == Size , "[ERROR] Permutation size doesn't match dimension" );
			return *this;
		}

		friend std::ostream &operator << ( std::ostream &os , const _Tensor &t ){ return os << t.data; }
	};

	// A general tensor
	template< unsigned int ... Dims >
	struct _Tensor< UIntPack< Dims ... > > : public Window< double , UIntPack< Dims ... > > , public InnerProductSpace< double , _Tensor< UIntPack< Dims ... > > >
	{
		typedef UIntPack< Dims ... > Pack;
		static const unsigned int Size = Pack::Size;

		_Tensor( void ){ memset( Window< double , Pack >::data , 0 , sizeof( double ) * WindowSize< Pack >::Size ); }

		template< typename ... UInts >
		double &operator()( unsigned int index , UInts ... indices )
		{
			static_assert( sizeof...(indices)==Pack::Size-1 , "[ERROR] Wrong number of indices" );
			unsigned int idx[] = { index , indices ... };
			return Window< double , UIntPack< Dims ... > >::operator()( idx );
		}

		template< typename ... UInts >
		const double &operator()( unsigned int index , UInts ... indices ) const
		{
			static_assert( sizeof...(indices)==Pack::Size-1 , "[ERROR] Wrong number of indices" );
			unsigned int idx[] = { index , indices ... };
			return Window< double , UIntPack< Dims ... > >::operator()( idx );
		}

		double &operator()( const unsigned int indices[] ){ return Window< double , UIntPack< Dims ... > >::operator()( indices ); }

		const double &operator()( const unsigned int indices[] ) const { return Window< double , UIntPack< Dims ... > >::operator()( indices ); }

		// Inner-product space methods
		void Add( const _Tensor &t )
		{
			WindowLoop< Size >::Run
			(
				ZeroUIntPack< Size >() , Pack() ,
				[]( int d , int i ){} ,
				[]( double &v1 , const double &v2 ){ v1 += v2; } ,
				*this , t
			);
		}
		void Scale( double s )
		{
			WindowLoop< Size >::Run
			(
				ZeroUIntPack< Size >() , Pack() ,
				[]( int d , int i ){} ,
				[&]( double &v ){ v *= s; } ,
				*this
			);
		}
		double InnerProduct( const _Tensor &t ) const
		{
			double innerProduct = 0;
			WindowLoop< Size >::Run
			(
				ZeroUIntPack< Size >() , Pack() ,
				[]( int d , int i ){} ,
				[&]( double v1 , double v2 ){ innerProduct += v1*v2; } ,
				*this , t
			);
			return innerProduct;
		}

		// Permute indices
		template< unsigned int ... PermutationValues >
		_Tensor< Permutation< Pack , UIntPack< PermutationValues ... > > > permute( UIntPack< PermutationValues ... > ) const
		{
			static_assert( sizeof ... ( PermutationValues ) == Size , "[ERROR] Permutation size doesn't match dimension" );
			typedef UIntPack< PermutationValues ... > PPack;
			const unsigned int PValues[] = { PermutationValues ... };
			unsigned int IPValues[ Size ];
			for( unsigned int i=0 ; i<Size ; i++ ) IPValues[ PValues[i] ] = i;

			_Tensor< Permutation< Pack , PPack > > t;
			unsigned int idx[ Size ];
			WindowLoop< Size >::Run
			(
				ZeroUIntPack< Size >() , Pack() ,
				[&]( int d , int i ){ idx[ IPValues[d] ] = i; } ,
				[&]( const double &v ){ t( idx ) = v; } ,
				*this
			);
			return t;
		}

		// Extract slice
#if 0
		template< unsigned int ... Indices >
		_Tensor< typename Split< sizeof ... ( Indices ) , Pack >::Second > extract( UIntPack< Indices ... > ) const
		{
			typedef typename Split< sizeof ... ( Indices ) , Pack >::Second Remainder;
			_Tensor< Remainder> t;

			WindowLoop< Remainder::Size >::Run
			(
				ZeroUIntPack< Remainder::Size >() , Remainder() ,
				[&]( int d , int i ){; } ,
				[&]( const double &v , double &v ){ t( idx ) = v; } ,
				*this , t
			);
			return t;
		}
#endif

		static auto TransposeTensor( void )
		{
			_Tensor< Concatenation< typename Pack::Transpose , Pack > > t;
			unsigned int idx[ 2*Size ];
			WindowLoop< Size >::Run
			(
				ZeroUIntPack< Size >() , Pack::Transpose() ,
				[&]( int d , int i ){ idx[d] = idx[ 2*Size - 1 - d ] = i; } ,
				[&]( void ){ t( idx ) = 1; }
			);
			return t;
		}

		// Transpose operator
		_Tensor< typename Pack::Transpose > transpose( void ) const
		{
			_Tensor< typename Pack::Transpose > t;
			unsigned int idx[ Size ];
			WindowLoop< Size >::Run
			(
				ZeroUIntPack< Size >() , Pack() ,
				[&]( int d , int i ){ idx[ Size - 1 - d ] = i; } ,
				[&]( const double &v ){ t( idx ) = v; } ,
				*this
			);
			return t;
		}

		// Outer product
		template< unsigned int ... _Dims >
		_Tensor< Concatenation< Pack , UIntPack< _Dims ... > > > operator * ( const _Tensor< UIntPack< _Dims ... > > &t ) const 
		{
			typedef UIntPack< _Dims ... > _Pack;
			_Tensor< Concatenation< Pack , _Pack > > _t;

			WindowLoop< Pack::Size >::Run
			(
				ZeroUIntPack< Pack::Size >() , Pack() ,
				[]( int d , int i ){} ,
				[&]( WindowWrapper< double , _Pack > __t , const double &v )
				{
					WindowLoop< _Pack::Size >::Run
					(
						ZeroUIntPack< _Pack::Size >() , _Pack() ,
						[]( int d , int i ){} ,
						[&]( double &v1 , const double &v2 ){ v1 += v*v2; } ,
						__t , t
					);
				} ,
				_t , *this
			);
			return _t;
		}

		_Tensor< Pack > operator * ( const _Tensor< UIntPack<> > &t ) const { return *this * t.data; }

	protected:
		template< unsigned int D1 , unsigned int D2 >
		static auto _ContractionTensor( void )
		{
			static_assert( D1<D2 , "[ERROR] Contraction indices are the same" );
			static_assert( D1<Pack::Size , "[ERROR] First contraction index too large" );
			static_assert( D2<Pack::Size , "[ERROR] Second contraction index too large" );
			static_assert( Select< D1 , Pack >::Value==Select< D2 , Pack >::Value , "[ERROR] Contraction dimensions differ" );
			typedef typename Select< D1 , typename Select< D2 , Pack >::Complement >::Complement OutPack;

			_Tensor< Concatenation< OutPack , Pack > > t;

			unsigned int index[ Pack::Size+OutPack::Size ];
			if constexpr( OutPack::Size==0 )
				for( unsigned int i=0 ; i<Pack::template Get<D1>() ; i++ )
				{
					index[D1] = index[D2] = i;
					t( index ) = 1;
				}
			else
			{
				unsigned int out2in[ OutPack::Size ];
				{
					unsigned int count = 0;
					for( unsigned int i=0 ; i<Pack::Size ; i++ ) if( i!=D1 && i!=D2 ) out2in[ count++ ] = i;
				}

				WindowLoop< OutPack::Size >::Run
				(
					ZeroUIntPack< OutPack::Size >() , OutPack() ,
					[&]( int d , int i ){ index[d] = index[ out2in[d] ] = i; } ,
					[&]( void )
					{
						for( unsigned int i=0 ; i<Pack::template Get<D1>() ; i++ )
						{
							index[ OutPack::Size+D1 ] = i;
							index[ OutPack::Size+D2 ] = i;
							t( index ) = 1;
						}
					}
				);
			}

			return t;
		}

	public:
		template< unsigned int D1 , unsigned int D2 >
		static auto ContractionTensor( void )
		{
			if constexpr( D1<D2 ) return _ContractionTensor< D1 , D2 >();
			else                  return _ContractionTensor< D2 , D1 >();
		}

		// Tensor contraction
		template< unsigned int D1 , unsigned int D2 >
		_Tensor< typename Select< (D1<D2?D2:D1) , typename Select< (D1<D2?D1:D2) , Pack >::Complement >::Complement > contract( void ) const
		{
			return ContractionTensor< D1 , D2 >().contractedOuterProduct< Pack::Size-2 >( *this );
		}

		// In1 := [ N{1} , ... , N{I} , N{I+1} , ... , N{K} ]
		// In2 :=                     [ N{I+1} , ... , N{K} , N{K+1} , ... N{M} ]
		// Out := [ N{1} , ... , N{I}             ,           N{K+1} , ... N{M} ]
		template< unsigned int I , unsigned int ... _Dims >
		_Tensor< Concatenation< typename Split< Size-I , Pack >::First , typename Split< I , UIntPack< _Dims ... > >::Second > > contractedOuterProduct( const _Tensor< UIntPack< _Dims ... > > &t ) const 
		{
			static_assert( Compare< typename Split< Size-I , Pack >::Second , typename Split< I , UIntPack< _Dims ... > >::First >::Equal , "[ERROR] Contraction suffix/prefix don't match" );
			typedef UIntPack< _Dims ... > _Pack;
			static const unsigned int _Size = _Pack::Size;
			typedef typename Split< Size-I ,  Pack >:: First P1;
			typedef typename Split< Size-I ,  Pack >::Second P2;
			typedef typename Split<      I , _Pack >::Second P3;

			typedef typename std::conditional< P2::Size!=0 , ConstWindowWrapper< double , P2 > , const double & >::type In1SliceType;
			typedef typename std::conditional< P3::Size!=0 , ConstWindowWrapper< double , P3 > , const double & >::type In2SliceType;
			typedef typename std::conditional< Concatenation< P1 , P3 >::Size!=0 , double , _Tensor< UIntPack<> > >::type OutBaseType;
			typedef typename std::conditional< P3::Size!=0 ,      WindowWrapper< double , P3 > , OutBaseType & >::type OutSliceType;

			const _Tensor<  Pack > &in1 = *this;
			const _Tensor< _Pack > &in2 = t;
			_Tensor< Concatenation< P1 , P3 > > out;

			// Iterate over {1,...,I} of in1 and out
			WindowLoop< P1::Size >::Run
			(
				ZeroUIntPack< P1::Size >() , P1() ,
				[]( int d , int i ){} ,
				[&]( In1SliceType _in1 , OutSliceType _out )
				{
					// Iterate over {I,...,K} of in1 and in2
					WindowLoop< P2::Size >::Run
					(
						ZeroUIntPack< P2::Size >() , P2() ,
						[]( int d , int i ){} ,
						[&]( double __in1 , In2SliceType _in2 )
						{
							// Iterate over {K+1,...,M} of in2 and out
							WindowLoop< P3::Size >::Run
							(
								ZeroUIntPack< P3::Size >() , P3() ,
								[]( int d , int i ){} ,
								[&]( double __in2 , OutBaseType &_out_ ){ _out_ += __in1 * __in2; } ,
								_in2 , _out
							);
						} ,
						_in1 , in2
					);
				} ,
				in1 , out
			);
			return out;
		}

		template< unsigned int I >
		_Tensor< Pack > contractedOuterProduct( const _Tensor< UIntPack<> > &t ) const { return *this * t; }
	};

	template< unsigned int ... Dims > using Tensor = _Tensor< UIntPack< Dims ... > >;
}
#endif // TENSORS_INCLUDED
