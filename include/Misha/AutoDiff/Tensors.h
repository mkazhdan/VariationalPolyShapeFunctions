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

namespace Misha
{
	template< typename Pack > struct Tensor;

	// A zero-tensor is the same as a double value
	template<>
	struct Tensor< UIntPack<> > : public StaticWindow< double , UIntPack<> > , public InnerProductSpace< double , Tensor< UIntPack<> > >
	{
		typedef  UIntPack<> Pack;
		static const unsigned int Size = Pack::Size;

		Tensor( double d=0 ) : StaticWindow< double , UIntPack<> >(d) {}

		void Add( const Tensor &t ){ data += t.data; }
		void Scale( double s ){ data *= s; }
		double InnerProduct( const Tensor &t ) const { return data * t.data; }
		template< unsigned int ... _Dims >
		Tensor< UIntPack< _Dims ... > > operator * ( const Tensor< UIntPack< _Dims ... > > &t ) const { return t * data; }
		template< unsigned int I , unsigned int ... _Dims >
		Tensor< UIntPack< _Dims ... > > contractedOuterProduct( const Tensor< UIntPack< _Dims ... > > &t ) const 
		{
			static_assert( I==0 , "[ERROR] Contraction suffix/prefix don't match" );
			return *this * t;
		}

		// Permute indices
		template< unsigned int ... PermutationValues >
		Tensor< Permutation< Pack , UIntPack< PermutationValues ... > > > permute( UIntPack< PermutationValues ... > ) const
		{
			static_assert( sizeof ... ( PermutationValues ) == Size , "[ERROR] Permutation size doesn't match dimension" );
			return *this;
		}

		friend std::ostream &operator << ( std::ostream &os , const Tensor &t ){ return os << t.data; }
	};


	// A general tensor
	template< unsigned int ... Dims >
	struct Tensor< UIntPack< Dims ... > > : public StaticWindow< double , UIntPack< Dims ... > > , public InnerProductSpace< double , Tensor< UIntPack< Dims ... > > >
	{
		typedef UIntPack< Dims ... > Pack;
		static const unsigned int Size = Pack::Size;

		Tensor( void ){ memset( StaticWindow< double , Pack >::data , 0 , sizeof( double ) * WindowSize< Pack >::Size ); }

		// Inner-product space methods
		void Add( const Tensor &t )
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
		double InnerProduct( const Tensor &t ) const
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
		Tensor< Permutation< Pack , UIntPack< PermutationValues ... > > > permute( UIntPack< PermutationValues ... > ) const
		{
			static_assert( sizeof ... ( PermutationValues ) == Size , "[ERROR] Permutation size doesn't match dimension" );
			typedef UIntPack< PermutationValues ... > PPack;
			const unsigned int PValues[] = { PermutationValues ... };
			unsigned int IPValues[ Size ];
			for( unsigned int i=0 ; i<Size ; i++ ) IPValues[ PValues[i] ] = i;

			Tensor< Permutation< Pack , PPack > > t;
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
		Tensor< typename Split< sizeof ... ( Indices ) , Pack >::Second > extract( UIntPack< Indices ... > ) const
		{
			typedef typename Split< sizeof ... ( Indices ) , Pack >::Second Remainder;
			Tensor< Remainder> t;

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

		// Transpose operator
		Tensor< typename Pack::Transpose > transpose( void ) const
		{
			Tensor< typename Pack::Transpose > t;
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
		Tensor< Concatenation< Pack , UIntPack< _Dims ... > > > operator * ( const Tensor< UIntPack< _Dims ... > > &t ) const 
		{
			typedef UIntPack< _Dims ... > _Pack;
			Tensor< Concatenation< Pack , _Pack > > _t;

			WindowLoop< Pack::Size >::Run
			(
				ZeroUIntPack< Pack::Size >() , Pack() ,
				[]( int d , int i ){} ,
				[&]( WindowSlice< double , _Pack > __t , const double &v )
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

		Tensor< Pack > operator * ( const Tensor< UIntPack<> > &t ) const { return *this * t.data; }

		// Tensor contraction
		template< unsigned int D1 , unsigned int D2 >
		Tensor< typename Select< (D1<D2?D2:D1) , typename Select< (D1<D2?D1:D2) , Pack >::Complement >::Complement > contract( void ) const
		{
			static_assert( Select< D1 , Pack >::Value == Select< D2 , Pack >::Value , "[ERROR] Values differ" );
			if constexpr( D2<D1 ){ return contract< D2 , D1 >(); }
			else
			{
				Tensor< typename Select< D2 , typename Select< D1 , Pack >::Complement >::Complement > t;

				unsigned int idx[ Size ] , _idx[ Size==2 ?  1 : Size - 2 ];
				WindowLoop< Size >::Run
				(
					ZeroUIntPack< Size >() , Pack() ,
					[&]( int d , int i ){ idx[d] = i; } ,
					[&]( const double &v )
					{
						if( idx[D1]==idx[D2] )
						{
							int ii=0;
							for( int d=0 ; d<sizeof ... (Dims) ; d++ ) if( d!=D1 && d!=D2 ) _idx[ii++] = idx[d];
							t( _idx ) += v;
						}
					} ,
					*this
				);
				return t;
			}
		};

		template< unsigned int I , unsigned int ... _Dims >
		Tensor< Concatenation< typename Split< Size-I , Pack >::First , typename Split< I , UIntPack< _Dims ... > >::Second > > contractedOuterProduct( const Tensor< UIntPack< _Dims ... > > &t ) const 
		{
			static_assert( Compare< typename Split< Size-I , Pack >::Second , typename Split< I , UIntPack< _Dims ... > >::First >::Equal , "[ERROR] Contraction suffix/prefix don't match" );
			typedef UIntPack< _Dims ... > _Pack;
			static const unsigned int _Size = _Pack::Size;
			typedef typename Split< Size-I ,  Pack >:: First P1;
			typedef typename Split< Size-I ,  Pack >::Second P2;
			typedef typename Split<      I , _Pack >::Second P3;

			const Tensor<  Pack > &in1 = *this;
			const Tensor< _Pack > &in2 = t;
			Tensor< Concatenation< P1 , P3 > > out;
			WindowLoop< P1::Size >::Run
			(
				ZeroUIntPack< P1::Size >() , P1() ,
				[]( int d , int i ){} ,
				[&]( ConstWindowSlice< double , P2 > _in1 , WindowSlice< double , P3 > _out )
				{
					WindowLoop< P2::Size >::Run
					(
						ZeroUIntPack< P2::Size >() , P2() ,
						[]( int d , int i ){} ,
						[&]( double __in1 , ConstWindowSlice< double , P3 > _in2 )
						{
							WindowLoop< P3::Size >::Run
							(
								ZeroUIntPack< P3::Size >() , P3() ,
								[]( int d , int i ){} ,
								[&]( double __in2 , double &_out_ ){ _out_ += __in1 * __in2; } ,
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
		Tensor< Pack > contractedOuterProduct( const Tensor< UIntPack<> > &t ) const { return *this * t; }
	};

	template< unsigned int ... Dims > Tensor< UIntPack< Dims ... > >& operator *= ( Tensor< UIntPack< Dims ... > > &t , int s ){ return t *= (double)s; };
	template< unsigned int ... Dims > Tensor< UIntPack< Dims ... > > operator * ( const Tensor< UIntPack< Dims ... > > &t , int s ){ return t * (double)s; };
	template< unsigned int ... Dims > Tensor< UIntPack< Dims ... > > operator * ( int s , const Tensor< UIntPack< Dims ... > > &t ){ return t * (double)s; };
}
#endif // TENSORS_INCLUDED
