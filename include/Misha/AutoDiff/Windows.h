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

namespace Misha
{
	////////////
	// Window //
	////////////


	template< typename T > struct WindowSize{};
	template< typename T1 , typename T2 > struct WindowIndex{};

	template< unsigned int Res , unsigned int ... Ress > struct WindowSize< UIntPack< Res , Ress ... > >{ static const unsigned int Size = WindowSize< UIntPack< Ress ... > >::Size * Res; };
	template< unsigned int Res                         > struct WindowSize< UIntPack< Res            > >{ static const unsigned int Size = Res; };
	template<                                          > struct WindowSize< UIntPack<                > >{ static const unsigned int Size = 1; };

	template< unsigned int Res , unsigned int ... Ress , unsigned int Idx , unsigned int ... Idxs > struct WindowIndex< UIntPack< Res , Ress ... > , UIntPack< Idx , Idxs ... > >{ static const unsigned int Index = Idx * WindowSize< UIntPack< Ress ... > >::Size + WindowIndex< UIntPack< Ress ... > , UIntPack< Idxs ... > >::Index; };
	template< unsigned int Res                         , unsigned int Idx                         > struct WindowIndex< UIntPack< Res            > , UIntPack< Idx            > >{ static const unsigned int Index = Idx; };
	template<                                                                                     > struct WindowIndex< UIntPack<                > , UIntPack<                > >{ static const unsigned int Index = 0; };

	template< unsigned int Res , unsigned int ... Ress > typename std::enable_if< (sizeof...(Ress)!=0) , unsigned int >::type GetWindowIndex( UIntPack< Res , Ress ... > , const unsigned int idx[] ){ return idx[0] * WindowSize< UIntPack< Ress ... > >::Size + GetWindowIndex( UIntPack< Ress ... >() , idx+1 ); };
	template< unsigned int Res                         >                                                 unsigned int         GetWindowIndex( UIntPack< Res >            , const unsigned int idx[] ){ return idx[0]; }

	template< unsigned int Res , unsigned int ... Ress > typename std::enable_if< (sizeof...(Ress)!=0) , unsigned int >::type GetWindowIndex( UIntPack< Res , Ress ... > , const int idx[] ){ return idx[0] * WindowSize< UIntPack< Ress ... > >::Size + GetWindowIndex( UIntPack< Ress ... >() , idx+1 ); };
	template< unsigned int Res                         >                                                 unsigned int         GetWindowIndex( UIntPack< Res >            , const int idx[] ){ return idx[0]; }

	template< typename Data , typename Pack > struct   ConstWindowSlice{};
	template< typename Data , typename Pack > struct        WindowSlice{};
	template< typename Data , typename Pack > struct  StaticWindow     {};
	template< typename Data , typename Pack > struct DynamicWindow     {};

	// A wrapper for storing multi-dimensional array data
	template< class Data , unsigned int ... Ress >
	struct ConstWindowSlice< Data , UIntPack< Ress ... > >
	{
		typedef UIntPack< Ress ... > Pack;
		static const unsigned int Size = WindowSize< Pack >::Size;
		typedef       Data                  data_type;
		typedef const Data&       data_reference_type;
		typedef const Data& const_data_reference_type;

		// A pointer to the multi-dimensional array data
		ConstPointer( Data ) data;

		ConstWindowSlice(      Pointer( Data ) d ) : data(d) {}
		ConstWindowSlice( ConstPointer( Data ) d ) : data(d) {}

		// Operator [] for slicing off a sub-array
		typename std::conditional< sizeof ... ( Ress )==1 , const Data & , ConstWindowSlice< Data , typename Pack::Rest > >::type operator[]( int idx ) const
		{
			if constexpr( sizeof ... ( Ress )==1 ) return data[idx];
			else return ConstWindowSlice< Data , typename Pack::Rest >( data + WindowSize< typename Pack::Rest >::Size * idx );
		}

		// Operator () for accessing elements of the multi-dimensional array using a multi-dimensional index
		const_data_reference_type operator()( const          int idx[/*Pack::Size*/] ) const { return data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		const_data_reference_type operator()( const unsigned int idx[/*Pack::Size*/] ) const { return data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		data_reference_type operator()( const          int idx[/*Pack::Size*/] ){ return data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		data_reference_type operator()( const unsigned int idx[/*Pack::Size*/] ){ return data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }

		// Operator << for writing the contents of the multi-dimensional array to a stream
		friend std::ostream &operator << ( std::ostream &os , const ConstWindowSlice &slice ) 
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
	};

	template< class Data >
	struct ConstWindowSlice< Data , UIntPack<> >
	{
		typedef UIntPack<> Pack;
		static const unsigned int Size = WindowSize< Pack >::Size;
		typedef       Data                  data_type;
		typedef const Data&       data_reference_type;
		typedef const Data& const_data_reference_type;

		ConstWindowSlice(       Data &d ) : data( GetPointer(d) ){}
		ConstWindowSlice( const Data &d ) : data( GetPointer(d) ){}

		// A pointer to the zero-dimensional array data
		ConstPointer( Data ) data;

		// operator for transforming the ConstWindowSlice to the data
		operator const Data &( void ) const { return *data; }

		// In the case that the datatype is a pointer
		auto operator ->( void ) const
		{
			static_assert( std::is_pointer< Data >::value , "[ERROR] Can only dereference pointer types" );
			return *data;
		}

		ConstWindowSlice(      Pointer( Data ) d ) : data(d) {}
		ConstWindowSlice( ConstPointer( Data ) d ) : data(d) {}

		// Operator () for accessing elements of the multi-dimensional array using a multi-dimensional index
		data_reference_type       operator()( const          int idx[/*Pack::Size*/] )       { return data[0]; }
		data_reference_type       operator()( const unsigned int idx[/*Pack::Size*/] )       { return data[0]; }
		const_data_reference_type operator()( const          int idx[/*Pack::Size*/] ) const { return data[0]; }
		const_data_reference_type operator()( const unsigned int idx[/*Pack::Size*/] ) const { return data[0]; }


		// Operator << for writing the contents of the multi-dimensional array to a stream
		friend std::ostream &operator << ( std::ostream &os , const ConstWindowSlice &slice ) { return os << *slice.data; }
	};


	template< class Data , unsigned int ... Ress >
	struct WindowSlice< Data , UIntPack< Ress ... > >
	{
		typedef UIntPack< Ress ... > Pack;
		static const unsigned int Size = WindowSize< Pack >::Size;
		typedef       Data                  data_type;
		typedef       Data&       data_reference_type;
		typedef const Data& const_data_reference_type;

		// A pointer to the multi-dimensional array data
		Pointer( Data ) data;

		WindowSlice( Pointer( Data ) d ) : data(d) {}

		// Operator [] for slicing off a sub-array
		typename std::conditional< sizeof ... ( Ress )==1 , Data & , WindowSlice< Data , typename Pack::Rest > >::type operator[]( int idx )
		{
			if constexpr( sizeof ... ( Ress )==1 ) return data[idx];
			else return WindowSlice< Data , typename Pack::Rest >( data + WindowSize< typename Pack::Rest >::Size * idx );
		}
		typename std::conditional< sizeof ... ( Ress )==1 , const Data & , ConstWindowSlice< Data , typename Pack::Rest > >::type operator[]( int idx ) const
		{
			if constexpr( sizeof ... ( Ress )==1 ) return data[idx];
			else return ConstWindowSlice< Data , typename Pack::Rest >( data + WindowSize< typename Pack::Rest >::Size * idx );
		}

		data_reference_type       operator()( const          int idx[/*Pack::Size*/] )       { return data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		data_reference_type       operator()( const unsigned int idx[/*Pack::Size*/] )       { return data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		const_data_reference_type operator()( const          int idx[/*Pack::Size*/] ) const { return data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }
		const_data_reference_type operator()( const unsigned int idx[/*Pack::Size*/] ) const { return data[ GetWindowIndex( UIntPack< Ress ... >() , idx ) ]; }

		operator ConstWindowSlice< Data , Pack >() const { return ConstWindowSlice< Data , Pack >( ( ConstPointer( Data ) )data ); }

		friend std::ostream &operator << ( std::ostream &os , const WindowSlice &slice ){ return os << ( ConstWindowSlice< Data , Pack > )slice; }
	};

	template< class Data >
	struct WindowSlice< Data , UIntPack<> >
	{
		typedef UIntPack<> Pack;
		static const unsigned int Size = WindowSize< Pack >::Size;
		typedef       Data                  data_type;
		typedef       Data&       data_reference_type;
		typedef const Data& const_data_reference_type;

		// A pointer to the zero-dimensional array data
		Pointer( Data ) data;

		WindowSlice( Data& d ) : data( GetPointer(d) ){}
		WindowSlice( Pointer( Data ) d ) : data(d) {}

		// Assignment directly from the data
		Data &operator = ( const Data &d ){ *data = d ; return *data; }

		// operator for transforming the WindowSlice to the data
		operator Data &( void ){ return *data; }
		operator const Data &( void ) const { return *data; }

		// In the case that the datatype is a pointer
		auto operator ->( void )
		{
			static_assert( std::is_pointer< Data >::value , "[ERROR] Can only dereference pointer types" );
			return *data;
		}
		auto operator ->( void ) const
		{
			static_assert( std::is_pointer< Data >::value , "[ERROR] Can only dereference pointer types" );
			return *data;
		}

		data_reference_type       operator()( const          int idx[/*Pack::Size*/] )       { return data[0]; }
		data_reference_type       operator()( const unsigned int idx[/*Pack::Size*/] )       { return data[0]; }
		const_data_reference_type operator()( const          int idx[/*Pack::Size*/] ) const { return data[0]; }
		const_data_reference_type operator()( const unsigned int idx[/*Pack::Size*/] ) const { return data[0]; }

		operator ConstWindowSlice< Data , Pack >() const { return ConstWindowSlice< Data , Pack >( ( ConstPointer( Data ) )data ); }

		friend std::ostream &operator << ( std::ostream &os , const WindowSlice &slice ){ return os << *slice.data; }
	};

	template< class Data , unsigned int ... Ress >
	struct StaticWindow< Data , UIntPack< Ress ... > > : public WindowSlice< Data , UIntPack< Ress ... > >
	{
		using WindowSlice< Data , UIntPack< Ress ... > >::operator();

		Data data[ WindowSize< UIntPack< Ress ... > >::Size ];

		StaticWindow( void ) : WindowSlice< Data , UIntPack< Ress ... > >( data ){}
		StaticWindow( const StaticWindow &w ) : WindowSlice< Data , UIntPack< Ress ... > >( data ){ memcpy( data , w.data , sizeof(Data) * WindowSize< UIntPack< Ress ... > >::Size ); }
		StaticWindow &operator=( const StaticWindow &w ){ memcpy( data , w.data , sizeof(Data) * WindowSize< UIntPack< Ress ... > >::Size ) ; return *this; }

		WindowSlice< Data , UIntPack< Ress ... > > operator()( void ){ return WindowSlice< Data , UIntPack< Ress ... > >( data ); }
		ConstWindowSlice< Data , UIntPack< Ress ... > > operator()( void ) const { return ConstWindowSlice< Data , UIntPack< Ress ... > >( data ); }
	};

	template< class Data >
	struct StaticWindow< Data , UIntPack<> > : public WindowSlice< Data , UIntPack<> >
	{
		using WindowSlice< Data , UIntPack<> >::operator();

		Data data;

		StaticWindow( void ) : WindowSlice< Data , UIntPack<> >( &data ){}
		StaticWindow( const Data &d ) : WindowSlice< Data , UIntPack<> >( &data ) , data(d) {}
		StaticWindow( const StaticWindow &w ) : WindowSlice< Data , UIntPack<> >( &data ){ data = w.data; }
		StaticWindow &operator=( const StaticWindow &w ){ data = w.data ; return *this; }

		WindowSlice< Data , UIntPack<> > operator()( void ){ return WindowSlice< Data , UIntPack<> >( &data ); }
		ConstWindowSlice< Data , UIntPack<> > operator()( void ) const { return ConstWindowSlice< Data , UIntPack<> >( &data ); }
	};

	template< typename Data >
	struct DynamicData
	{
		Pointer( Data ) data;

		DynamicData( size_t sz ){ data = data = NewPointer< Data >( sz ); }
		~DynamicData( void ){ DeletePointer( data ); }
	};

	template< class Data , unsigned int ... Ress >
	struct DynamicWindow< Data , UIntPack< Ress ... > > : public DynamicData< Data > , WindowSlice< Data , UIntPack< Ress ... > >
	{
		using DynamicData< Data >::data;
		using WindowSlice< Data , UIntPack< Ress ... > >::operator();

		DynamicWindow( void ) : DynamicData< Data >( WindowSize< UIntPack< Ress ... > >::Size ) , WindowSlice< Data , UIntPack< Ress ... > >( DynamicData< Data >::data ){}
		DynamicWindow( const DynamicWindow &w ) : DynamicData< Data >( WindowSize< UIntPack< Ress ... > >::Size ) , WindowSlice< Data , UIntPack< Ress ... > >( data ){ memcpy( data , w.data , sizeof(Data) * WindowSize< UIntPack< Ress ... > >::Size ); }
		DynamicWindow &operator=( const DynamicWindow &w ){ memcpy( data , w.data , sizeof(Data) * WindowSize< UIntPack< Ress ... > >::Size ) ; return *this; }

		WindowSlice< Data , UIntPack< Ress ... > > operator()( void ){ return WindowSlice< Data , UIntPack< Ress ... > >( DynamicData< Data >::data ); }
		ConstWindowSlice< Data , UIntPack< Ress ... > > operator()( void ) const { return ConstWindowSlice< Data , UIntPack< Ress ... > >( DynamicData< Data >::data ); }
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
#include "Windows.inl"
}
#endif // WINDOWS_INCLUDED
