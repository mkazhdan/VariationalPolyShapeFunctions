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

#ifndef UINT_PACK_INCLUDED
#define UINT_PACK_INCLUDED

//////////////////////////////////////
// Unsigned integer parameter packs //
//////////////////////////////////////

namespace UIntPack
{
	////////////////////////////
	////////////////////////////
	//// Short declarations ////
	////////////////////////////
	////////////////////////////

	// A wrapper class for passing unsigned integer parameter packs
	template< unsigned int  ... Values > struct Pack;

	// A class that identifies a single entry and its complement
	template< unsigned int I , typename Pack > struct Selection;

	// A class that splits a Pack into two sub-packs
	template< unsigned int I , typename Pack > struct Partition;

	// A class for comparing two Packs
	template< typename Pack1 , typename Pack2 > struct Comparison;

	// A helper class for defining a concatenation of multiple Packs
	template< typename ... Packs > struct _Concatenation;

	// A helper class for defining the permutation of a Pack
	template< typename Pack , typename PermutationPack > struct _Permutation;

	// A helper class for defining a Pack with the same value repeated Dim times
	template< unsigned int Dim , unsigned int Value > struct _IsotropicPack;

	// A helper class for defining a Pack with sequentially increasing values
	template< unsigned int Dim , unsigned int StartValue > struct _SequentialPack;

	// A Pack that is the concatenation of multiple Packs
	template< typename ... Packs > using Concatenation = typename _Concatenation< Packs ... >::type;

	// A Pack that is the permtuation of a Pack
	// [NOTE] The entry in the i-th position  of PermutationPack indicates where the i-th position comes from (not goes to)
	template< typename Pack , typename PermutationPack > using Permutation = typename _Permutation< Pack , PermutationPack >::type;

	// A Pack that has the same value repeated Dim times
	template< unsigned int Dim , unsigned int Value=0 > using IsotropicPack = typename _IsotropicPack< Dim , Value >::type;

	// A Pack with sequentially increasing values
	template< unsigned int Dim , unsigned int StartValue=0 > using SequentialPack = typename _SequentialPack< Dim , StartValue >::type;


	/////////////////////
	/////////////////////
	//// Definitions ////
	/////////////////////
	/////////////////////


	//////////
	// Pack //
	//////////

	// The general case
	template< unsigned int _Value , unsigned int ... _Values > struct Pack< _Value , _Values ... >
	{
		static const unsigned int First = _Value;
		typedef Pack< _Values ... > Rest;
		typedef typename Rest::Transpose::template Append< First > Transpose;

		static const unsigned int Size = 1 + sizeof ... ( _Values );

		template< unsigned int ... __Values > using  Append = Pack< _Value , _Values ... , __Values ... >;
		template< unsigned int ... __Values > using Prepend = Pack< __Values ... , _Value , _Values ... >;

		static const unsigned int Values[];

		template< unsigned int I > constexpr static unsigned int Get( void )
		{
			if constexpr( I==0 ) return _Value;
			else return Rest::template Get< I-1 >();
		}

		friend std::ostream &operator << ( std::ostream &os , Pack )
		{
			os << "< ";
			for( unsigned int i=0 ; i<Size ; i++ )
			{
				if( i ) os << " , ";
				os << Values[i];
			}
			return os << " >";
		}
	};

	// The specialized case with one entry
	template< unsigned int _Value > struct Pack< _Value >
	{
		static const unsigned int First = _Value;
		typedef Pack<> Rest;
		typedef Pack< _Value > Transpose;

		static const unsigned int Size = 1;

		template< unsigned int ... __Values > using  Append = Pack< _Value , __Values ... >;
		template< unsigned int ... __Values > using Prepend = Pack< __Values ... , _Value >;

		static const unsigned int Values[];

		template< unsigned int I > constexpr static unsigned int Get( void )
		{
			static_assert( I==0 , "[ERROR] Pack< Value >::Get called with non-zero index" );
			return _Value;
		}

		friend std::ostream &operator << ( std::ostream &os , Pack )
		{
			return os << "< " << First << " >";
		}
	};

	// The specialized case with no entries
	template<> struct Pack<>
	{
		typedef Pack<> Rest;
		static const unsigned int Size = 0;
		static constexpr unsigned int Values[] = { 0 };
		typedef Pack<> Transpose;
		template< unsigned int ... __Values > using  Append = Pack< __Values ... >;
		template< unsigned int ... __Values > using Prepend = Pack< __Values ... >;
		friend std::ostream &operator << ( std::ostream &os , Pack ){ return os << "< >"; }
	};

	template< unsigned int _Value , unsigned int ... _Values > const unsigned int Pack< _Value , _Values ... >::Values[] = { _Value , _Values ... };
	template< unsigned int _Value > const unsigned int Pack< _Value >::Values[] = { _Value };

	///////////////
	// Selection //
	///////////////
	template< unsigned int I , unsigned _Value , unsigned int ... _Values >
	struct Selection< I , Pack< _Value , _Values ... > >
	{
		static const unsigned int Value = Selection< I-1 , Pack< _Values ... > >::Value;
		typedef typename Selection< I-1 , Pack< _Values ... > >::Complement::template Prepend< _Value > Complement;
	};

	template< unsigned _Value , unsigned int ... _Values >
	struct Selection< 0 , Pack< _Value , _Values ... > >
	{
		static const unsigned int Value = _Value;
		typedef Pack< _Values ... > Complement;
	};

	///////////////
	// Partition //
	///////////////
	template< unsigned int ... Values >
	struct Partition< 0 , Pack< Values ... > >
	{
		typedef Pack<> First;
		typedef Pack< Values ... > Second;
	};

	template< unsigned int I , unsigned int ... Values >
	struct Partition< I , Pack< Values ... > >
	{
		typedef Concatenation< Pack< Pack< Values ... >::First > , typename Partition< I-1 , typename Pack< Values ... >::Rest >::First > First;
		typedef typename Partition< I-1 , typename Pack< Values ... >::Rest >::Second Second;
	};

	////////////////
	// Comparison //
	////////////////
	template< unsigned int ... Values1 , unsigned int ... Values2 >
	struct Comparison< Pack< Values1 ... > , Pack< Values2 ... > >
	{
		typedef Pack< Values1 ... > Pack1;
		typedef Pack< Values2 ... > Pack2;
		static const bool              Equal = Pack1::First==Pack2::First && Comparison< typename Pack1::Rest , typename Pack2::Rest >::Equal;
		static const bool           NotEqual = Pack1::First!=Pack2::First || Comparison< typename Pack1::Rest , typename Pack2::Rest >::NotEqual;
		static const bool    LessThan        = Pack1::First< Pack2::First && Comparison< typename Pack1::Rest , typename Pack2::Rest >::LessThan;
		static const bool    LessThanOrEqual = Pack1::First<=Pack2::First && Comparison< typename Pack1::Rest , typename Pack2::Rest >::LessThanOrEqual;
		static const bool GreaterThan        = Pack1::First> Pack2::First && Comparison< typename Pack1::Rest , typename Pack2::Rest >::GreaterThan;
		static const bool GreaterThanOrEqual = Pack1::First>=Pack2::First && Comparison< typename Pack1::Rest , typename Pack2::Rest >::GreaterThanOrEqual;
	};

	template< unsigned int Value1 , unsigned int Value2 >
	struct Comparison< Pack< Value1 > , Pack< Value2 > >
	{
		static const bool Equal = Value1==Value2;
		static const bool NotEqual = Value1!=Value2;
		static const bool LessThan = Value1<Value2;
		static const bool LessThanOrEqual = Value1<=Value2;
		static const bool GreaterThan = Value1>Value2;
		static const bool GreaterThanOrEqual = Value1>=Value2;
	};

	template<>
	struct Comparison< Pack<> , Pack<> >
	{
		static const bool Equal = true;
		static const bool NotEqual = false;
		static const bool LessThan = false;
		static const bool LessThanOrEqual = true;
		static const bool GreaterThan = false;
		static const bool GreaterThanOrEqual = true;
	};

	////////////////////
	// _Concatenation //
	////////////////////
	template< unsigned int ... Values1 , unsigned int ... Values2 , typename ... Packs >
	struct _Concatenation< Pack< Values1 ... > , Pack< Values2 ... > , Packs ... >
	{
		typedef typename _Concatenation< typename Pack< Values1 ... >::template Append< Values2 ... > , Packs ... >::type type;
	};
	template< unsigned int ... Values >
	struct _Concatenation< Pack< Values ... > >
	{
		typedef Pack< Values ... > type;
	};

	//////////////////
	// _Permutation //
	//////////////////
	template< unsigned int ... Values , unsigned int ... PermutationValues >
	struct _Permutation< Pack< Values ... > , Pack< PermutationValues ... > >
	{
		typedef Pack< PermutationValues ... > PPack;
		typedef Concatenation< Pack< Selection< PPack::First , Pack< Values ... > >::Value > , typename _Permutation< Pack< Values ... > , typename PPack::Rest >::type > type;
	};
	template< unsigned int ... Values >
	struct _Permutation< Pack< Values ... > , Pack<> >
	{
		typedef Pack<> type;
	};

	////////////////////
	// _IsotropicPack //
	////////////////////
	template< unsigned int Dim , unsigned int Value > struct _IsotropicPack             { typedef typename _IsotropicPack< Dim-1 , Value >::type::template Append< Value > type; };
	template<                    unsigned int Value > struct _IsotropicPack< 1 , Value >{ typedef Pack< Value > type; };
	template<                    unsigned int Value > struct _IsotropicPack< 0 , Value >{ typedef Pack< > type; };

	/////////////////////
	// _SequentialPack //
	/////////////////////
	template< unsigned int Dim , unsigned int Value >   struct _SequentialPack             { typedef Concatenation< Pack< Value > , typename _SequentialPack< Dim-1 , Value+1 >::type > type; };
	template<                    unsigned int Value >   struct _SequentialPack< 0 , Value >{ typedef Pack<> type; };
}
#endif // UINT_PACK_INCLUDED
