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
#include <random>
#include "Misha/Algebra.h"
#include "MultiDimensionalArray.h"
#include "UIntPack.h"

namespace AutoDiff
{
	template< typename Pack > struct Tensor;

	// A zero-tensor is the same as a double value
	template<>
	struct Tensor< UIntPack::Pack<> > : public InnerProductSpace< double , Tensor< UIntPack::Pack<> > >
	{
		typedef UIntPack::Pack<> Pack;
		static const unsigned int Size = 1;

		double data;

		Tensor( double d=0 ) : data(d) {}
		explicit operator double &( void ){ return data; }
		explicit operator const double &( void ) const { return data; }

		void Add( const Tensor &t ){ data += t.data; }
		void Scale( double s ){ data *= s; }
		double InnerProduct( const Tensor &t ) const { return data * t.data; }
		template< unsigned int ... _Dims >
		Tensor< UIntPack::Pack< _Dims ... > > operator * ( const Tensor< UIntPack::Pack< _Dims ... > > &t ) const { return t * data; }
		template< unsigned int I , unsigned int ... _Dims >
		Tensor< UIntPack::Pack< _Dims ... > > contractedOuterProduct( const Tensor< UIntPack::Pack< _Dims ... > > &t ) const 
		{
			static_assert( I==0 , "[ERROR] Contraction suffix/prefix don't match" );
			return *this * t;
		}

		// Permute indices
		template< unsigned int ... PermutationValues >
		Tensor< UIntPack::Permutation< Pack , UIntPack::Pack< PermutationValues ... > > > permute( UIntPack::Pack< PermutationValues ... > ) const
		{
			static_assert( sizeof ... ( PermutationValues ) == Size , "[ERROR] Permutation size doesn't match dimension" );
			return *this;
		}

		static Tensor Random( std::default_random_engine &generator )
		{
			// From https://www.cplusplus.com/reference/random/uniform_real_distribution/
			std::uniform_real_distribution< double > distribution( 0.0 , 1.0 );

			return Tensor( distribution( generator ) );
		}

		friend std::ostream &operator << ( std::ostream &os , const Tensor &t ){ return os << t.data; }
	};

	// A general tensor
	template< unsigned int ... Dims >
	struct Tensor< UIntPack::Pack< Dims ... > > : public MultiDimensionalArray::Array< double , Dims ... > , public InnerProductSpace< double , Tensor< UIntPack::Pack< Dims ... > > >
	{
		typedef UIntPack::Pack< Dims ... > Pack;
		static const unsigned int Size = Pack::Size;

		Tensor( void ){ memset( MultiDimensionalArray::Array< double , Dims ... >::data , 0 , sizeof( double ) * MultiDimensionalArray::ArraySize< Dims ... >() ); }

		template< typename ... UInts >
		double &operator()( unsigned int index , UInts ... indices )
		{
			static_assert( sizeof...(indices)==Pack::Size-1 , "[ERROR] Wrong number of indices" );
			unsigned int idx[] = { index , indices ... };
			return MultiDimensionalArray::Array< double , Dims ... >::operator()( idx );
		}

		template< typename ... UInts >
		const double &operator()( unsigned int index , UInts ... indices ) const
		{
			static_assert( sizeof...(indices)==Pack::Size-1 , "[ERROR] Wrong number of indices" );
			unsigned int idx[] = { index , indices ... };
			return MultiDimensionalArray::Array< double , Dims ... >::operator()( idx );
		}

		double &operator()( const unsigned int indices[] ){ return MultiDimensionalArray::Array< double , Dims ... >::operator()( indices ); }
		const double &operator()( const unsigned int indices[] ) const { return MultiDimensionalArray::Array< double , Dims ... >::operator()( indices ); }

		// Inner-product space methods
		void Add( const Tensor &t )
		{
			MultiDimensionalArray::Loop< Size >::Run
			(
				UIntPack::IsotropicPack< Size >::Values , Pack::Values ,
				[]( int d , int i ){} ,
				[]( double &v1 , const double &v2 ){ v1 += v2; } ,
				*this , t
			);
		}
		void Scale( double s )
		{
			MultiDimensionalArray::Loop< Size >::Run
			(
				UIntPack::IsotropicPack< Size >::Values , Pack::Values ,
				[]( int d , int i ){} ,
				[&]( double &v ){ v *= s; } ,
				*this
			);
		}
		double InnerProduct( const Tensor &t ) const
		{
			double innerProduct = 0;
			MultiDimensionalArray::Loop< Size >::Run
			(
				UIntPack::IsotropicPack< Size >::Values , Pack::Values ,
				[]( int d , int i ){} ,
				[&]( double v1 , double v2 ){ innerProduct += v1*v2; } ,
				*this , t
			);
			return innerProduct;
		}

		static Tensor Random( std::default_random_engine &generator )
		{
			// From https://www.cplusplus.com/reference/random/uniform_real_distribution/
			std::uniform_real_distribution< double > distribution( 0.0 , 1.0 );

			Tensor t;

			MultiDimensionalArray::Loop< Size >::Run
			(
				UIntPack::IsotropicPack< Size >::Values , Pack::Values ,
				[&]( int , int ){} ,
				[&]( double &v ){ v = distribution( generator ); } ,
				t
			);
			return t;
		}


		template< unsigned int ... PermutationValues >
		static auto PermutationTensor( UIntPack::Pack< PermutationValues ... > )
		{
#pragma message( "[WARNING] Should avoid using PermutationTensor" )
			WARN_ONCE( "Invoking PermutationTensor" );
			Tensor< UIntPack::Concatenation< UIntPack::Permutation< Pack , UIntPack::Pack< PermutationValues ... > > , Pack > > t;
			const unsigned int permutation[] = { PermutationValues ... };
			unsigned int idx[ 2*Size ];
			MultiDimensionalArray::Loop< Size >::Run
			(
				UIntPack::IsotropicPack< Size >::Values , Pack::Values ,
				[&]( int d , int i ){ idx[ permutation[d] ] = idx[ Size+d ] = i; } ,
				[&]( void ){ t( idx ) = 1; }
			);
			return t;
		}

		// Permute indices
		template< unsigned int ... PermutationValues >
		Tensor< UIntPack::Permutation< Pack , UIntPack::Pack< PermutationValues ... > > > permute( UIntPack::Pack< PermutationValues ... > ) const
		{
			static_assert( sizeof ... ( PermutationValues ) == Size , "[ERROR] Permutation size doesn't match dimension" );
			typedef UIntPack::Pack< PermutationValues ... > PPack;
			const unsigned int PValues[] = { PermutationValues ... };
			unsigned int IPValues[ Size ];
			for( unsigned int i=0 ; i<Size ; i++ ) IPValues[ PValues[i] ] = i;

			Tensor< UIntPack::Permutation< Pack , PPack > > t;
			unsigned int idx[ Size ];
			MultiDimensionalArray::Loop< Size >::Run
			(
				UIntPack::IsotropicPack< Size >::Values , Pack::Values ,
				[&]( int d , int i ){ idx[ IPValues[d] ] = i; } ,
				[&]( const double &v ){ t( idx ) = v; } ,
				*this
			);
			return t;
		}

		// Extract slice
		template< unsigned int I >
		auto extract( const unsigned int indices[/*I*/] ) const
		{
			typedef typename UIntPack::Partition< I , Pack >::Second Remainder;
			Tensor< Remainder> t;

			if constexpr( Remainder::Size!=0 )
			{
				unsigned int _indices[ Pack::Size ];
				for( unsigned int i=0 ; i<I ; i++ ) _indices[i] = indices[i];

				MultiDimensionalArray::Loop< Remainder::Size >::Run
				(
					UIntPack::IsotropicPack< Remainder::Size >::Values , Remainder::Values ,
					[&]( int d , int i ){ _indices[d+I] = i; } ,
					[&]( double &_t ){ _t = operator()( _indices ); } ,
					t
				);
			}
			else static_cast< double & >( t ) = operator()( indices );
			return t;
		}

		static auto TransposeTensor( void )
		{
#pragma message( "[WARNING] Should avoid using TransposeTensor" )
			WARN_ONCE( "Invoking TransposeTensor" );
			Tensor< UIntPack::Concatenation< typename Pack::Transpose , Pack > > t;
			unsigned int idx[ 2*Size ];
			MultiDimensionalArray::Loop< Size >::Run
			(
				UIntPack::IsotropicPack< Size >::Values , Pack::Transpose::Values ,
				[&]( int d , int i ){ idx[d] = idx[ 2*Size - 1 - d ] = i; } ,
				[&]( void ){ t( idx ) = 1; }
			);
			return t;
		}

		// Transpose operator
		Tensor< typename Pack::Transpose > transpose( void ) const
		{
			return permute( UIntPack::SequentialPack< Pack::Size >::Transpose() );
		}

		// Outer product
		template< unsigned int ... _Dims >
		Tensor< UIntPack::Concatenation< Pack , UIntPack::Pack< _Dims ... > > > operator * ( const Tensor< UIntPack::Pack< _Dims ... > > &t ) const 
		{
			typedef UIntPack::Pack< _Dims ... > _Pack;
			Tensor< UIntPack::Concatenation< Pack , _Pack > > _t;

			MultiDimensionalArray::Loop< Pack::Size >::Run
			(
				UIntPack::IsotropicPack< Pack::Size >::Values , Pack::Values ,
				[]( int d , int i ){} ,
				[&]( MultiDimensionalArray::ArrayWrapper< double , _Dims ... > __t , const double &v )
				{
					MultiDimensionalArray::Loop< _Pack::Size >::Run
					(
						UIntPack::IsotropicPack< _Pack::Size >::Values , _Pack::Values ,
						[]( int d , int i ){} ,
						[&]( double &v1 , const double &v2 ){ v1 += v*v2; } ,
						__t , t
					);
				} ,
				_t , *this
					);
			return _t;
		}

		Tensor< Pack > operator * ( const Tensor< UIntPack::Pack<> > &t ) const { return *this * t.data; }

	protected:
		template< unsigned int D1 , unsigned int D2 >
		static auto _ContractionTensor( void )
		{
			static_assert( D1<D2 , "[ERROR] Contraction indices are the same" );
			static_assert( D1<Pack::Size , "[ERROR] First contraction index too large" );
			static_assert( D2<Pack::Size , "[ERROR] Second contraction index too large" );
			static_assert( UIntPack::Selection< D1 , Pack >::Value==UIntPack::Selection< D2 , Pack >::Value , "[ERROR] Contraction dimensions differ" );
			typedef typename UIntPack::Selection< D1 , typename UIntPack::Selection< D2 , Pack >::Complement >::Complement OutPack;

			Tensor< UIntPack::Concatenation< OutPack , Pack > > t;

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

				MultiDimensionalArray::Loop< OutPack::Size >::Run
				(
					UIntPack::IsotropicPack< OutPack::Size >::Values , OutPack::Values ,
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
#pragma message( "[WARNING] Should avoid using ContractionTensor" )
			WARN_ONCE( "Invoking ContractionTensor" );
			if constexpr( D1<D2 ) return _ContractionTensor< D1 , D2 >();
			else                  return _ContractionTensor< D2 , D1 >();
		}

		// Tensor contraction
		template< unsigned int I1 , unsigned int I2 >
		auto contract( void ) const
		{
			static_assert( I1!=I2 , "[ERROR] Contraction indices must differ" );
			static_assert( Pack::template Get< I1 >()==Pack::template Get< I2 >() , "[ERROR] Contraction dimensions don't match" );
			static_assert( I1<Pack::Size && I2<Pack::Size , "[ERROR] Contraction indices out of bounds" );
			if constexpr( I2<I1 ) return this->template contract< I2 , I1 >();
			typedef typename UIntPack::Selection< I1 , typename UIntPack::Selection< I2 , Pack >::Complement >::Complement OutPack;
			Tensor< OutPack > out;
			if constexpr( Pack::Size>2 )
			{
				unsigned int indices[ OutPack::Size ];
				MultiDimensionalArray::Loop< OutPack::Size >::Run
				(
					UIntPack::IsotropicPack< OutPack::Size >::Values , OutPack::Values ,
					[&]( int d , int i ){ indices[d] = i; } ,
					[&]( double &_out )
					{
						unsigned int _indices[ Pack::Size ];
						unsigned int idx=0;
						for( unsigned int i=0 ; i<Pack::Size ; i++ ) if( i!=I1 && i!=I2 ) _indices[i] = indices[idx++];
						_out = 0;
						for( unsigned int i=0 ; i<Pack::template Get< I1 >() ; i++ )
						{
							_indices[I1] = _indices[I2] = i;
							_out += operator()( _indices );
						}
					} ,
					out
						);
			}
			else
			{
				double &_out = static_cast< double & >( out );
				unsigned int _indices[2];
				for( unsigned int i=0 ; i<Pack::template Get<I1>() ; i++ )
				{
					_indices[I1] = _indices[I2] = i;
					_out += operator()( _indices );
				}
			}
			return out;
		}

		// In1 := [ N{1} , ... , N{I} , N{I+1} , ... , N{K} ]
		// In2 :=                     [ N{I+1} , ... , N{K} , N{K+1} , ... N{M} ]
		// Out := [ N{1} , ... , N{I}             ,           N{K+1} , ... N{M} ]
		template< unsigned int I , unsigned int ... _Dims >
		Tensor< UIntPack::Concatenation< typename UIntPack::Partition< Size-I , Pack >::First , typename UIntPack::Partition< I , UIntPack::Pack< _Dims ... > >::Second > > contractedOuterProduct( const Tensor< UIntPack::Pack< _Dims ... > > &t ) const 
		{
			static_assert( UIntPack::Comparison< typename UIntPack::Partition< Size-I , Pack >::Second , typename UIntPack::Partition< I , UIntPack::Pack< _Dims ... > >::First >::Equal , "[ERROR] Contraction suffix/prefix don't match" );
			typedef UIntPack::Pack< _Dims ... > _Pack;
			static const unsigned int _Size = _Pack::Size;
			typedef typename UIntPack::Partition< Size-I ,  Pack >:: First P1;
			typedef typename UIntPack::Partition< Size-I ,  Pack >::Second P2;
			typedef typename UIntPack::Partition<      I , _Pack >::Second P3;

			typedef typename MultiDimensionalArray::SliceType< P1::Size , double ,  Dims ... >::const_type In1SliceType;
			typedef typename MultiDimensionalArray::SliceType< P2::Size , double , _Dims ... >::const_type In2SliceType;
			// In the case that we are collapsing completely, out is of type Tensor< UIntPack::Pack<> >
			// -- Then the first and last loops are trivial and we never access the contents of out using operator[]
			typedef typename std::conditional< UIntPack::Concatenation< P1 , P3 >::Size!=0 , double , Tensor< UIntPack::Pack<> > >::type OutBaseType;
			typedef typename std::conditional< P3::Size!=0 , typename MultiDimensionalArray::SliceType< P2::Size , double , _Dims ... >::type , OutBaseType & >::type OutSliceType;

			const Tensor<  Pack > &in1 = *this;
			const Tensor< _Pack > &in2 = t;
			Tensor< UIntPack::Concatenation< P1 , P3 > > out;

			// Iterate over {1,...,I} of in1 and out
			MultiDimensionalArray::Loop< P1::Size >::Run
			(
				UIntPack::IsotropicPack< P1::Size >::Values , P1::Values ,
				[]( int d , int i ){} ,
				[&]( In1SliceType _in1 , OutSliceType _out )
				{
					// Iterate over {I,...,K} of in1 and in2
					MultiDimensionalArray::Loop< P2::Size >::Run
					(
						UIntPack::IsotropicPack< P2::Size >::Values , P2::Values ,
						[]( int d , int i ){} ,
						[&]( double __in1 , In2SliceType _in2 )
						{
							// Iterate over {K+1,...,M} of in2 and out
							MultiDimensionalArray::Loop< P3::Size >::Run
							(
								UIntPack::IsotropicPack< P3::Size >::Values , P3::Values ,
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
		Tensor< Pack > contractedOuterProduct( const Tensor< UIntPack::Pack<> > &t ) const { return *this * t; }
	};
}
#endif // TENSORS_INCLUDED
