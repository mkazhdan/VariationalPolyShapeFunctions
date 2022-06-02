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

#ifndef WINDOWS_INCLUDED
#define WINDOWS_INCLUDED

#include <functional>
#include <iostream>
#include "Misha/Exceptions.h"
#include "Misha/Array.h"
#include "UIntPack.h"

//#define USE_DYNAMIC_WINDOW

namespace AutoDiff
{
	////////////
	// Window //
	////////////

	template< typename T > struct WindowSize{};
	template< typename T1 , typename T2 > struct WindowIndex{};

	template< unsigned int Res , unsigned int ... Ress > struct WindowSize< UIntPack< Res , Ress ... > >{ static const unsigned int Size = WindowSize< UIntPack< Ress ... > >::Size * Res; };
	template< unsigned int Res                         > struct WindowSize< UIntPack< Res            > >{ static const unsigned int Size = Res; };

	template< unsigned int Res , unsigned int ... Ress , unsigned int Idx , unsigned int ... Idxs > struct WindowIndex< UIntPack< Res , Ress ... > , UIntPack< Idx , Idxs ... > >{ static const unsigned int Index = Idx * WindowSize< UIntPack< Ress ... > >::Size + WindowIndex< UIntPack< Ress ... > , UIntPack< Idxs ... > >::Index; };
	template< unsigned int Res                         , unsigned int Idx                         > struct WindowIndex< UIntPack< Res            > , UIntPack< Idx            > >{ static const unsigned int Index = Idx; };

	template< unsigned int Res , unsigned int ... Ress > typename std::enable_if< (sizeof...(Ress)!=0) , unsigned int >::type GetWindowIndex( UIntPack< Res , Ress ... > , const unsigned int idx[] ){ return idx[0] * WindowSize< UIntPack< Ress ... > >::Size + GetWindowIndex( UIntPack< Ress ... >() , idx+1 ); };
	template< unsigned int Res                         >                                                 unsigned int         GetWindowIndex( UIntPack< Res >            , const unsigned int idx[] ){ return idx[0]; }

	template< unsigned int Res , unsigned int ... Ress > typename std::enable_if< (sizeof...(Ress)!=0) , unsigned int >::type GetWindowIndex( UIntPack< Res , Ress ... > , const int idx[] ){ return idx[0] * WindowSize< UIntPack< Ress ... > >::Size + GetWindowIndex( UIntPack< Ress ... >() , idx+1 ); };
	template< unsigned int Res                         >                                                 unsigned int         GetWindowIndex( UIntPack< Res >            , const int idx[] ){ return idx[0]; }

	template< typename Data , typename Pack > struct   ConstWindowWrapper{};
	template< typename Data , typename Pack > struct        WindowWrapper{};
	template< typename Data , typename Pack > struct        Window       {};

	// A wrapper for accessing a pointer as a multi-dimensional array
	template< class Data , unsigned int ... Ress >
	struct ConstWindowWrapper< Data , UIntPack< Ress ... > >
	{
		typedef UIntPack< Ress ... > Pack;
		static const unsigned int Size = WindowSize< Pack >::Size;
		typedef       Data                  data_type;
		typedef const Data&       data_reference_type;
		typedef const Data& const_data_reference_type;

		// The constructors wrap the pointer
		ConstWindowWrapper(      Pointer( Data ) d ) : _data(d) {}
		ConstWindowWrapper( ConstPointer( Data ) d ) : _data(d) {}

		// Operator [] for slicing off a sub-array
		// -- When wrapper is a one-dimensional array, return the indexed data.
		// -- Otherwise return the sub-array.
		typename std::conditional< sizeof ... ( Ress )==1 , const Data & , ConstWindowWrapper< Data , typename Pack::Rest > >::type operator[]( int idx ) const
		{
			if constexpr( sizeof ... ( Ress )==1 ) return _data[idx];
			else return ConstWindowWrapper< Data , typename Pack::Rest >( _data + WindowSize< typename Pack::Rest >::Size * idx );
		}

		// Operator () for accessing elements of the multi-dimensional array using a multi-dimensional index
		const_data_reference_type operator()( const          int idx[/*Pack::Size*/] ) const { return _data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		const_data_reference_type operator()( const unsigned int idx[/*Pack::Size*/] ) const { return _data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		data_reference_type operator()( const          int idx[/*Pack::Size*/] ){ return _data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		data_reference_type operator()( const unsigned int idx[/*Pack::Size*/] ){ return _data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }

		// Operator << for writing the contents of the multi-dimensional array to a stream
		friend std::ostream &operator << ( std::ostream &os , const ConstWindowWrapper &slice ) 
		{
			os << "{ ";
			for( unsigned int i=0 ; i<Pack::First ; i++ )
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


	template< class Data , unsigned int ... Ress >
	struct WindowWrapper< Data , UIntPack< Ress ... > >
	{
		typedef UIntPack< Ress ... > Pack;
		static const unsigned int Size = WindowSize< Pack >::Size;
		typedef       Data                  data_type;
		typedef       Data&       data_reference_type;
		typedef const Data& const_data_reference_type;

		WindowWrapper( Pointer( Data ) d ) : _data(d) {}

		// Operator [] for slicing off a sub-array
		typename std::conditional< sizeof ... ( Ress )==1 , Data & , WindowWrapper< Data , typename Pack::Rest > >::type operator[]( int idx )
		{
			if constexpr( sizeof ... ( Ress )==1 ) return _data[idx];
			else return WindowWrapper< Data , typename Pack::Rest >( _data + WindowSize< typename Pack::Rest >::Size * idx );
		}
		typename std::conditional< sizeof ... ( Ress )==1 , const Data & , ConstWindowWrapper< Data , typename Pack::Rest > >::type operator[]( int idx ) const
		{
			if constexpr( sizeof ... ( Ress )==1 ) return _data[idx];
			else return ConstWindowWrapper< Data , typename Pack::Rest >( _data + WindowSize< typename Pack::Rest >::Size * idx );
		}

		data_reference_type       operator()( const          int idx[/*Pack::Size*/] )       { return _data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		data_reference_type       operator()( const unsigned int idx[/*Pack::Size*/] )       { return _data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		const_data_reference_type operator()( const          int idx[/*Pack::Size*/] ) const { return _data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		const_data_reference_type operator()( const unsigned int idx[/*Pack::Size*/] ) const { return _data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }

		operator ConstWindowWrapper< Data , Pack >() const { return ConstWindowWrapper< Data , Pack >( ( ConstPointer( Data ) )_data ); }

		friend std::ostream &operator << ( std::ostream &os , const WindowWrapper &slice ){ return os << ( ConstWindowWrapper< Data , Pack > )slice; }
	protected:
		// A pointer to the zero-dimensional array data
		Pointer( Data ) _data;
#ifdef USE_DYNAMIC_WINDOW
		WindowWrapper( void ) : _data( NullPointer< Data >() ) {}
		void _init( Pointer( Data ) data ){ _data = data; }
		template< class _Data , typename Pack > friend struct Window;
#endif // USE_DYNAMIC_WINDOW
	};

	template< class Data , unsigned int ... Ress >
	struct Window< Data , UIntPack< Ress ... > > : public WindowWrapper< Data , UIntPack< Ress ... > >
	{
		using WindowWrapper< Data , UIntPack< Ress ... > >::operator();

#ifdef USE_DYNAMIC_WINDOW
		Pointer( Data ) data;

		Window( void )
		{
			data = NewPointer< Data >( WindowSize< UIntPack< Ress ... > >::Size );
			WindowWrapper< Data , UIntPack< Ress ... > >::_init( data );
		}
		Window( const Window &w )
		{
			data = NewPointer< Data >( WindowSize< UIntPack< Ress ... > >::Size );
			WindowWrapper< Data , UIntPack< Ress ... > >::_init( data );
			memcpy( data , w.data , sizeof(Data) * WindowSize< UIntPack< Ress ... > >::Size );
		}
		Window &operator=( const Window &w ){ memcpy( data , w.data , sizeof(Data) * WindowSize< UIntPack< Ress ... > >::Size ) ; return *this; }
		~Window( void ){ DeletePointer( data ); }
#else // !USE_DYNAMIC_WINDOW
		Data data[ WindowSize< UIntPack< Ress ... > >::Size ];

		Window( void ) : WindowWrapper< Data , UIntPack< Ress ... > >( data ){}
		Window( const Window &w ) : WindowWrapper< Data , UIntPack< Ress ... > >( data ){ memcpy( data , w.data , sizeof(Data) * WindowSize< UIntPack< Ress ... > >::Size ); }
		Window &operator=( const Window &w ){ memcpy( data , w.data , sizeof(Data) * WindowSize< UIntPack< Ress ... > >::Size ) ; return *this; }
#endif // USE_DYNAMIC_WINDOW

		WindowWrapper< Data , UIntPack< Ress ... > > operator()( void ){ return WindowWrapper< Data , UIntPack< Ress ... > >( data ); }
		ConstWindowWrapper< Data , UIntPack< Ress ... > > operator()( void ) const { return ConstWindowWrapper< Data , UIntPack< Ress ... > >( data ); }
	};

	// Recursive loop iterations for processing window slices
	//		N: the number of dimensions to process
	//		I: the current processing iteration of the window ( <=N )

	template< unsigned int I , unsigned int N > struct _WindowLoop;

	template< unsigned int N >
	struct WindowLoop
	{
		template< typename UpdateFunction , typename ProcessFunction , class ... Windows >
		static void Run( int begin , int end , UpdateFunction updateState , ProcessFunction function , Windows&& ... w )
		{
			_WindowLoop< 0 , N >::Run( begin , end , updateState , function , std::forward< Windows >( w ) ... ); 
		}
		template< typename UpdateFunction , typename ProcessFunction , class ... Windows >
		static void Run( const int* begin , const int* end , UpdateFunction updateState , ProcessFunction function , Windows&& ... w )
		{
			_WindowLoop< 0 , N >::Run( begin , end , updateState , function , std::forward< Windows >( w ) ... ); 
		}
		template< unsigned int ... Begin , unsigned int ... End , typename UpdateFunction , typename ProcessFunction , class ... Windows >
		static void Run( UIntPack< Begin ... > begin , UIntPack< End ... > end , UpdateFunction updateState , ProcessFunction function , Windows&& ... w )
		{
			_WindowLoop< 0 , N >::Run( begin , end , updateState , function , std::forward< Windows >( w ) ... ); 
		}
	};

	template< unsigned int I , unsigned int N >
	struct _WindowLoop
	{
		static_assert( I<=N , "[ERROR] Index exceeds bound" );
	protected:
		friend struct WindowLoop< N >;
		template< unsigned int _I , unsigned int _N > friend struct _WindowLoop;

		template< typename UpdateFunction , typename ProcessFunction , class ... Values >
		static void Run( int begin , int end , UpdateFunction& updateState , ProcessFunction& function , Values&& ... v )
		{
			if constexpr( I<N )
				for( int i=begin ; i<end ; i++ ){ updateState( I , i ) ; _WindowLoop< I+1 , N >::Run( begin , end , updateState , function , std::forward< decltype(v[i]) >( v[i] ) ... ); }
			else
				function( std::forward< Values >( v ) ... );
		}
		template< typename UpdateFunction , typename ProcessFunction , class ... Values >
		static void Run( const int* begin , const int* end , UpdateFunction& updateState , ProcessFunction& function , Values&& ... v )
		{
			if constexpr( I<N )
				for( int i=begin[0] ; i<end[0] ; i++ ){ updateState( I , i ) ; _WindowLoop< I+1 , N >::Run( begin+1 , end+1 , updateState , function , std::forward< decltype(v[i]) >( v[i] ) ... ); }
			else
				function( std::forward< Values >( v ) ... );
		}
		template< unsigned int ... Begin , unsigned int ... End , typename UpdateFunction , typename ProcessFunction , class ... Values >
		static void Run( UIntPack< Begin ... > begin , UIntPack< End ... > end , UpdateFunction& updateState , ProcessFunction& function , Values&& ... v )
		{
			if constexpr( I<N )
				for( int i=UIntPack< Begin ... >::First ; i<UIntPack< End ... >::First ; i++ ){ updateState( I , i ) ; _WindowLoop< I+1 , N >::Run( typename UIntPack< Begin ... >::Rest() , typename UIntPack< End ... >::Rest() , updateState , function , std::forward< decltype(v[i]) >( v[i] ) ... ); }
			else
				function( std::forward< Values >( v ) ... );
		}
	};
}
#endif // WINDOWS_INCLUDED
