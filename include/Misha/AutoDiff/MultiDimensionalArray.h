/*
Copyright (c) 2016, Michael Kazhdan
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

#ifndef MULTI_DIMENSIONAL_ARRAY_INCLUDED
#define MULTI_DIMENSIONAL_ARRAY_INCLUDED

#include <functional>
#include <iostream>
#include "Misha/Exceptions.h"
#include "Misha/Array.h"

//#define USE_DYNAMIC_MULTI_DIMENSIONAL_ARRAY

namespace MultiDimensionalArray
{
	///////////////////////////
	// MultiDimensionalArray //
	///////////////////////////

	template< unsigned int Res , unsigned int ... Ress >
	constexpr unsigned int ArraySize( void )
	{
		if constexpr( sizeof...(Ress)!=0 ) return Res * ArraySize< Ress ... >();
		else                               return Res;
	}

	template< unsigned int Res , unsigned int ... Ress >
	constexpr unsigned int ArrayIndex( const unsigned int idx[] )
	{
		if constexpr( sizeof...(Ress)!=0 ) return idx[0] * ArraySize< Ress ... >() + ArrayIndex< Ress ... >( idx+1 );
		else                               return idx[0];
	};

	template< unsigned int Res , unsigned int ... Ress > 
	constexpr unsigned int ArrayIndex( const int idx[] )
	{
		if constexpr( sizeof...(Ress)!=0 ) return idx[0] * ArraySize< Ress ... >() + ArrayIndex< Ress ... >( idx+1 );
		else                               return idx[0];
	};

	// Wrappers for accessing a pointer as a multi-dimensional array
	template< class Data , unsigned int ... Ress > struct ConstArrayWrapper{};
	template< class Data , unsigned int ... Ress > struct      ArrayWrapper{};


	// The types representing an I-th dimensional slice of a multi-dimensional array
	template< unsigned int I , class Data , unsigned int ... Ress > struct SliceType;

	// Base slice: Return the wrapper itself
	template< class Data , unsigned int Res , unsigned int ... Ress >
	struct SliceType< 0 , Data , Res , Ress ... >
	{
		typedef ConstArrayWrapper< Data , Res , Ress ... > const_type;
		typedef      ArrayWrapper< Data , Res , Ress ... >       type;
	};

	// Specialized base slice: Return a reference to the data
	template< class Data >
	struct SliceType< 0 , Data >
	{
		typedef const Data & const_type;
		typedef       Data &       type;
	};

	// Partial slice: Recurse
	template< unsigned int I , class Data , unsigned int Res , unsigned int ... Ress >
	struct SliceType< I , Data , Res , Ress ... >
	{
		typedef typename SliceType< I-1 , Data , Ress ... >::const_type const_type;
		typedef typename SliceType< I-1 , Data , Ress ... >::      type       type;
	};

	template< class Data , unsigned int Res , unsigned int ... Ress >
	struct ConstArrayWrapper< Data , Res , Ress ... >
	{
		static const unsigned int Size = ArraySize< Res , Ress ... >();
		typedef       Data                  data_type;
		typedef const Data&       data_reference_type;
		typedef const Data& const_data_reference_type;

		// The constructors wrap the pointer
		ConstArrayWrapper(      Pointer( Data ) d ) : _data(d) {}
		ConstArrayWrapper( ConstPointer( Data ) d ) : _data(d) {}

		// Operator [] for slicing off a sub-array
		// -- When wrapper is a one-dimensional array, return the indexed data.
		// -- Otherwise return the sub-array.
		typename std::conditional< sizeof ... ( Ress )==0 , const Data & , ConstArrayWrapper< Data , Ress ... > >::type operator[]( unsigned int idx ) const
		{
			if constexpr( sizeof ... ( Ress )==0 ) return _data[idx];
			else return ConstArrayWrapper< Data , Ress ... >( _data + ArraySize< Ress ... >() * idx );
		}

		// Operator () for accessing elements of the multi-dimensional array using a multi-dimensional index
		const_data_reference_type operator()( const          int idx[] ) const { return _data[ ArrayIndex< Res , Ress ... >( idx ) ]; }
		const_data_reference_type operator()( const unsigned int idx[] ) const { return _data[ ArrayIndex< Res , Ress ... >( idx ) ]; }
		data_reference_type       operator()( const          int idx[] )       { return _data[ ArrayIndex< Res , Ress ... >( idx ) ]; }
		data_reference_type       operator()( const unsigned int idx[] )       { return _data[ ArrayIndex< Res , Ress ... >( idx ) ]; }

		// Operator << for writing the contents of the multi-dimensional array to a stream
		friend std::ostream &operator << ( std::ostream &os , const ConstArrayWrapper &slice ) 
		{
			os << "{ ";
			for( unsigned int i=0 ; i<Res ; i++ )
			{
				if( i ) os << " , ";
				os << slice[i];
			}
			os << " }";
			return os;
		}
	protected:
		// A pointer to the multi-dimensional array data
		ConstPointer( Data ) _data;
	};

	template< class Data , unsigned int Res , unsigned int ... Ress >
	struct ArrayWrapper< Data , Res , Ress ... >
	{
		static const unsigned int Size = ArraySize< Res , Ress ... >();
		typedef       Data                  data_type;
		typedef       Data&       data_reference_type;
		typedef const Data& const_data_reference_type;

		ArrayWrapper( Pointer( Data ) d ) : _data(d) {}

		// Operator [] for slicing off a sub-array
		typename std::conditional< sizeof ... ( Ress )==0 , Data & , ArrayWrapper< Data , Ress ... > >::type operator[]( unsigned int idx )
		{
			if constexpr( sizeof ... ( Ress )==0 ) return _data[idx];
			else return ArrayWrapper< Data , Ress ... >( _data + ArraySize< Ress ... >() * idx );
		}
		typename std::conditional< sizeof ... ( Ress )==0 , const Data & , ConstArrayWrapper< Data , Ress ... > >::type operator[]( unsigned int idx ) const
		{
			if constexpr( sizeof ... ( Ress )==0 ) return _data[idx];
			else return ConstArrayWrapper< Data , Ress ... >( _data + ArraySize< Ress ... >() * idx );
		}

		data_reference_type       operator()( const          int idx[] )       { return _data[ ArrayIndex< Res , Ress ... >( idx ) ]; }
		data_reference_type       operator()( const unsigned int idx[] )       { return _data[ ArrayIndex< Res , Ress ... >( idx ) ]; }
		const_data_reference_type operator()( const          int idx[] ) const { return _data[ ArrayIndex< Res , Ress ... >( idx ) ]; }
		const_data_reference_type operator()( const unsigned int idx[] ) const { return _data[ ArrayIndex< Res , Ress ... >( idx ) ]; }

		operator ConstArrayWrapper< Data , Res , Ress ... >() const { return ConstArrayWrapper< Data , Res , Ress ... >( ( ConstPointer( Data ) )_data ); }

		friend std::ostream &operator << ( std::ostream &os , const ArrayWrapper &slice ){ return os << ( ConstArrayWrapper< Data , Res , Ress ... > )slice; }
	protected:
		// A pointer to the zero-dimensional array data
		Pointer( Data ) _data;
#ifdef USE_DYNAMIC_MULTI_DIMENSIONAL_ARRAY
		ArrayWrapper( void ) : _data( NullPointer< Data >() ) {}
		void _init( Pointer( Data ) data ){ _data = data; }
		template< class _Data , unsigned int _Res , unsigned int ... _Ress > friend struct Array;
#endif // USE_DYNAMIC_MULTI_DIMENSIONAL_ARRAY
	};

	template< typename Data , unsigned int Res , unsigned int ... Ress >
	struct Array : public ArrayWrapper< Data , Res , Ress ... >
	{
		using ArrayWrapper< Data , Res , Ress ... >::operator();

#ifdef USE_DYNAMIC_MULTI_DIMENSIONAL_ARRAY
		Pointer( Data ) data;

		Array( void )
		{
			data = NewPointer< Data >( ArraySize< Res , Ress ... >() );
			ArrayWrapper< Data , Res , Ress ... >::_init( data );
		}
		Array( const Array &mda )
		{
			data = NewPointer< Data >( ArraySize< Res , Ress ... >() );
			ArrayWrapper< Data , Res , Ress ... >::_init( data );
#ifdef _WIN64
#pragma message( "[WARNING] Why does Windows make me pull the size into a separate variable?" )
			static const unsigned int Size = ArraySize< Res , Ress ... >();
			memcpy( data , mda.data , sizeof(Data) * Size );
#else // !_WIN64
			memcpy( data , mda.data , sizeof(Data) * ArraySize< Res , Ress ... >();
#endif // _WIN64
		}
		Array &operator=( const Array &mda ){ memcpy( data , mda.data , sizeof(Data) * ArraySize< Res , Ress ... >() ) ; return *this; }
		~Array( void ){ DeletePointer( data ); }
#else // !USE_DYNAMIC_MULTI_DIMENSIONAL_ARRAY
		Data data[ ArraySize< Res , Ress ... >() ];

		Array( void ) : ArrayWrapper< Data , Res , Ress ... >( data ){}
		Array( const Array &mda ) : ArrayWrapper< Data , Res , Ress ... >( data )
		{
#ifdef _WIN64
#pragma message( "[WARNING] Why does Windows make me pull the size into a separate variable?" )
			static const unsigned int Size = ArraySize< Res , Ress ... >();
			memcpy( data , mda.data , sizeof(Data) * Size );
#else // !_WIN64
			memcpy( data , mda.data , sizeof(Data) * ArraySize< Res , Ress ... >() );
#endif // _WIN64
		}
		Array &operator=( const Array &mda ){ memcpy( data , mda.data , sizeof(Data) * ArraySize< Res , Ress ... >() ) ; return *this; }
#endif // USE_DYNAMIC_MULTI_DIMENSIONAL_ARRAY

		ArrayWrapper< Data , Res , Ress ... > operator()( void ){ return ArrayWrapper< Data , Res , Ress ... >( data ); }
		ConstArrayWrapper< Data , Res , Ress ... > operator()( void ) const { return ConstArrayWrapper< Data , Res , Ress ... >( data ); }
	};

	template< unsigned int N >
	struct Loop
	{
		template< typename UpdateFunction , typename ProcessFunction , class ... Arrays >
		static void Run( const unsigned int *begin , const unsigned int *end , UpdateFunction updateState , ProcessFunction function , Arrays&& ... mda )
		{
			_Run< 0 >( begin , end , updateState , function , std::forward< Arrays >( mda ) ... ); 
		}
	protected:
		template< unsigned int I , typename UpdateFunction , typename ProcessFunction , class ... Values >
		static void _Run( const unsigned int *begin , const unsigned int *end , UpdateFunction& updateState , ProcessFunction& function , Values&& ... v )
		{
			if constexpr( I<N )
				for( unsigned int i=begin[0] ; i<end[0] ; i++ ){ updateState( I , i ) ; _Run< I+1 >( begin+1 , end+1 , updateState , function , std::forward< decltype(v[i]) >( v[i] ) ... ); }
			else
				function( std::forward< Values >( v ) ... );
		}
	};
}
#endif // MULTI_DIMENSIONAL_ARRAY_INCLUDED
