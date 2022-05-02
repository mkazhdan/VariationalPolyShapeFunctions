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

namespace AutoDiff
{

	// A wrapper class for passing unsigned integer parameter packs
	template< unsigned int  ... Values > struct UIntPack{};

	template< unsigned int _Value , unsigned int ... _Values > struct UIntPack< _Value , _Values ... >
	{
		static const unsigned int First = _Value;
		typedef UIntPack< _Values ... > Rest;
		typedef typename Rest::Transpose::template Append< First > Transpose;

		static const unsigned int Size = 1 + sizeof ... ( _Values );
		template< unsigned int ... __Values > using  Append = UIntPack< _Value , _Values ... , __Values ... >;
		template< unsigned int ... __Values > using Prepend = UIntPack< __Values ... , _Value , _Values ... >;

		static const unsigned int Values[];
		static constexpr unsigned int Min( void ){ return _Value < Rest::Min() ? _Value : Rest::Min(); }
		static constexpr unsigned int Max( void ){ return _Value > Rest::Max() ? _Value : Rest::Max(); }

		template< typename T > struct Plus{};
		template< typename T > struct Minus{};
		template< typename T > struct Compare{};
		template< unsigned int __Value , unsigned int ... __Values > struct Plus < UIntPack< __Value , __Values ... > >{ typedef typename Rest::template Plus < UIntPack< __Values ... > >::type::template Prepend< _Value + __Value > type; };
		template< unsigned int __Value , unsigned int ... __Values > struct Minus< UIntPack< __Value , __Values ... > >{ typedef typename Rest::template Minus< UIntPack< __Values ... > >::type::template Prepend< _Value - __Value > type; };
		template< unsigned int __Value , unsigned int ... __Values > struct Compare< UIntPack< __Value , __Values ... > >
		{
			static const bool              Equal = _Value==__Value && Rest::template Compare< UIntPack< __Values ... > >::             Equal;
			static const bool           NotEqual = _Value!=__Value || Rest::template Compare< UIntPack< __Values ... > >::          NotEqual;
			static const bool    LessThan        = _Value< __Value && Rest::template Compare< UIntPack< __Values ... > >::   LessThan       ;
			static const bool    LessThanOrEqual = _Value<=__Value && Rest::template Compare< UIntPack< __Values ... > >::   LessThanOrEqual;
			static const bool GreaterThan        = _Value> __Value && Rest::template Compare< UIntPack< __Values ... > >::GreaterThan       ;
			static const bool GreaterThanOrEqual = _Value>=__Value && Rest::template Compare< UIntPack< __Values ... > >::GreaterThanOrEqual;
		};

		template< unsigned int I > constexpr static typename std::enable_if< I==0 , unsigned int >::type Get( void ){ return _Value; }
		template< unsigned int I > constexpr static typename std::enable_if< I!=0 , unsigned int >::type Get( void ){ return Rest::template Get< I-1 >(); }

		template< unsigned int __Value , unsigned int ... __Values > constexpr bool operator <  ( UIntPack< __Value , __Values ... > ) const { return _Value< __Value && Rest()< UIntPack< __Values ... >(); }
		template< unsigned int __Value , unsigned int ... __Values > constexpr bool operator <= ( UIntPack< __Value , __Values ... > ) const { return _Value<=__Value && Rest()<=UIntPack< __Values ... >(); }
		template< unsigned int __Value , unsigned int ... __Values > constexpr bool operator >  ( UIntPack< __Value , __Values ... > ) const { return _Value> __Value && Rest()> UIntPack< __Values ... >(); }
		template< unsigned int __Value , unsigned int ... __Values > constexpr bool operator >= ( UIntPack< __Value , __Values ... > ) const { return _Value>=__Value && Rest()>=UIntPack< __Values ... >(); }
		template< unsigned int __Value , unsigned int ... __Values > constexpr bool operator == ( UIntPack< __Value , __Values ... > ) const { return _Value==__Value && Rest()==UIntPack< __Values ... >(); }
		template< unsigned int __Value , unsigned int ... __Values > constexpr bool operator != ( UIntPack< __Value , __Values ... > ) const { return _Value!=__Value && Rest()!=UIntPack< __Values ... >(); }

		friend std::ostream &operator << ( std::ostream &os , UIntPack )
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
	template< unsigned int _Value > struct UIntPack< _Value >
	{
		static const unsigned int First = _Value;
		typedef UIntPack<> Rest;
		typedef UIntPack< _Value > Transpose;

		static const unsigned int Size = 1;
		template< unsigned int ... __Values > using  Append = UIntPack< _Value , __Values ... >;
		template< unsigned int ... __Values > using Prepend = UIntPack< __Values ... , _Value >;
		static const unsigned int Values[];
		static constexpr unsigned int Min( void ){ return _Value; }
		static constexpr unsigned int Max( void ){ return _Value; }

		template< typename T > struct Plus{};
		template< typename T > struct Minus{};
		template< typename T > struct Compare{};
		template< unsigned int __Value > struct Plus < UIntPack< __Value > >{ typedef UIntPack< _Value + __Value > type; };
		template< unsigned int __Value > struct Minus< UIntPack< __Value > >{ typedef UIntPack< _Value - __Value > type; };
		template< unsigned int __Value > struct Compare< UIntPack< __Value > >
		{
			static const bool              Equal = _Value==__Value;
			static const bool           NotEqual = _Value!=__Value;
			static const bool    LessThan        = _Value< __Value;
			static const bool    LessThanOrEqual = _Value<=__Value;
			static const bool GreaterThan        = _Value> __Value;
			static const bool GreaterThanOrEqual = _Value>=__Value;
		};

		template< unsigned int I > constexpr static unsigned int Get( void ){ static_assert( I==0 , "[ERROR] UIntPack< Value >::Get called with non-zero index" ) ; return _Value; }

		template< unsigned int __Value > constexpr bool operator <  ( UIntPack< __Value > ) const { return _Value< __Value; }
		template< unsigned int __Value > constexpr bool operator <= ( UIntPack< __Value > ) const { return _Value<=__Value; }
		template< unsigned int __Value > constexpr bool operator >  ( UIntPack< __Value > ) const { return _Value> __Value; }
		template< unsigned int __Value > constexpr bool operator >= ( UIntPack< __Value > ) const { return _Value>=__Value; }
		template< unsigned int __Value > constexpr bool operator == ( UIntPack< __Value > ) const { return _Value==__Value; }
		template< unsigned int __Value > constexpr bool operator != ( UIntPack< __Value > ) const { return _Value!=__Value; }

		friend std::ostream &operator << ( std::ostream &os , UIntPack )
		{
			return os << "< " << First << " >";
		}
	};


	template<> struct UIntPack<>
	{
		typedef UIntPack<> Rest;
		static const unsigned int Size = 0;
		static constexpr unsigned int Values[] = { 0 };
		typedef UIntPack<> Transpose;
		template< unsigned int ... __Values > using  Append = UIntPack< __Values ... >;
		template< unsigned int ... __Values > using Prepend = UIntPack< __Values ... >;
		friend std::ostream &operator << ( std::ostream &os , UIntPack ){ return os << "< >"; }
	};

	template< unsigned int _Value , unsigned int ... _Values > const unsigned int UIntPack< _Value , _Values ... >::Values[] = { _Value , _Values ... };
	template< unsigned int _Value > const unsigned int UIntPack< _Value >::Values[] = { _Value };
	template< unsigned int ... V1 , unsigned int ... V2 > typename UIntPack< V1 ... >::template Plus < UIntPack< V2 ... > >::type operator + ( UIntPack< V1 ... > , UIntPack< V2 ... > ){ return typename UIntPack< V1 ... >::template Plus < UIntPack< V2 ... > >::type(); }
	template< unsigned int ... V1 , unsigned int ... V2 > typename UIntPack< V1 ... >::template Minus< UIntPack< V2 ... > >::type operator - ( UIntPack< V1 ... > , UIntPack< V2 ... > ){ return typename UIntPack< V1 ... >::template Minus< UIntPack< V2 ... > >::type(); }

	/////////////////////////////////////////////////////////
	// Selection of an individual index and its complement //
	/////////////////////////////////////////////////////////
	template< unsigned int I , typename Pack > struct Select;

	template< unsigned int I , unsigned _Value , unsigned int ... _Values >
	struct Select< I , UIntPack< _Value , _Values ... > >
	{
		static const unsigned int Value = Select< I-1 , UIntPack< _Values ... > >::Value;
		typedef typename Select< I-1 , UIntPack< _Values ... > >::Complement::template Prepend< _Value > Complement;
	};

	template< unsigned _Value , unsigned int ... _Values >
	struct Select< 0 , UIntPack< _Value , _Values ... > >
	{
		static const unsigned int Value = _Value;
		typedef UIntPack< _Values ... > Complement;
	};

	///////////////////
	// Concatenation //
	///////////////////
	template< typename ... Packs > struct _Concatenation;
	template< unsigned int ... Values1 , unsigned int ... Values2 , typename ... Packs >
	struct _Concatenation< UIntPack< Values1 ... > , UIntPack< Values2 ... > , Packs ... >
	{
		typedef typename _Concatenation< typename UIntPack< Values1 ... >::template Append< Values2 ... > , Packs ... >::type type;
	};
	template< unsigned int ... Values >
	struct _Concatenation< UIntPack< Values ... > >
	{
		typedef UIntPack< Values ... > type;
	};

	template< typename ... Packs > using Concatenation = typename _Concatenation< Packs ... >::type;

	///////////
	// Power //
	///////////
	template< typename Pack , unsigned int P > struct _Power;
	template< unsigned int ... Values , unsigned int P >
	struct _Power< UIntPack< Values ... > , P >
	{
		typedef Concatenation< UIntPack< Values ... > , typename _Power< UIntPack< Values ... > , P-1 >::type > type;
	};
	template< unsigned int ... Values >
	struct _Power< UIntPack< Values ... > , 0 >
	{
		typedef UIntPack<> type;
	};

	template< typename Pack , unsigned int P > using Power = typename _Power< Pack , P >::type;

	///////////
	// Split //
	///////////
	template< unsigned int I , typename Pack > struct Split;

	template< unsigned int ... Values >
	struct Split< 0 , UIntPack< Values ... > >
	{
		typedef UIntPack<> First;
		typedef UIntPack< Values ... > Second;
	};

	template< unsigned int I , unsigned int ... Values >
	struct Split< I , UIntPack< Values ... > >
	{
		typedef UIntPack< Values ... > Pack;
		typedef Concatenation< UIntPack< Pack::First > , typename Split< I-1 , typename Pack::Rest >::First > First;
		typedef typename Split< I-1 , typename Pack::Rest >::Second Second;
	};

	///////////////
	// Insertion //
	///////////////
	template< unsigned int I , typename Pack1 , typename Pack2 > struct _Insertion;
	template< unsigned int I , unsigned int ... Values1 , unsigned int ... Values2 >
	struct _Insertion< I , UIntPack< Values1 ... > , UIntPack< Values2 ... > >
	{
		typedef UIntPack< Values1 ... > Pack1;
		typedef UIntPack< Values2 ... > Pack2;
		typedef Concatenation< Concatenation< typename Split< I , Pack1 >::First , Pack2 > , typename Split< I , Pack2 >::Second > type;
	};
	template< unsigned int I , typename Pack1 , typename Pack2 >
	using Insertion = typename _Insertion< I , Pack1 , Pack2 >::type;

	/////////////////
	// Permutation //
	/////////////////
	template< typename Pack , typename PermutationPack > struct _Permutation;
	template< unsigned int ... Values , unsigned int ... PermutationValues >
	struct _Permutation< UIntPack< Values ... > , UIntPack< PermutationValues ... > >
	{
		typedef UIntPack< Values ... > Pack;
		typedef UIntPack< PermutationValues ... > PPack;
		typedef Concatenation< UIntPack< Select< PPack::First , Pack >::Value > , typename _Permutation< Pack , typename PPack::Rest >::type > type;
	};
	template< unsigned int ... Values >
	struct _Permutation< UIntPack< Values ... > , UIntPack<> >
	{
		typedef UIntPack<> type;
	};

	template< typename Pack , typename PermutationPack > using Permutation = typename _Permutation< Pack , PermutationPack >::type;


	template< typename Pack1 , typename Pack2 > struct Compare;

	template< unsigned int ... Values1 , unsigned int ... Values2 >
	struct Compare< UIntPack< Values1 ... > , UIntPack< Values2 ... > >
	{
		typedef UIntPack< Values1 ... > Pack1;
		typedef UIntPack< Values2 ... > Pack2;
		static const bool              Equal = Pack1::First==Pack2::First && Compare< typename Pack1::Rest , typename Pack2::Rest >::Equal;
		static const bool           NotEqual = Pack1::First!=Pack2::First || Compare< typename Pack1::Rest , typename Pack2::Rest >::NotEqual;
		static const bool    LessThan        = Pack1::First< Pack2::First && Compare< typename Pack1::Rest , typename Pack2::Rest >::LessThan;
		static const bool    LessThanOrEqual = Pack1::First<=Pack2::First && Compare< typename Pack1::Rest , typename Pack2::Rest >::LessThanOrEqual;
		static const bool GreaterThan        = Pack1::First> Pack2::First && Compare< typename Pack1::Rest , typename Pack2::Rest >::GreaterThan;
		static const bool GreaterThanOrEqual = Pack1::First>=Pack2::First && Compare< typename Pack1::Rest , typename Pack2::Rest >::GreaterThanOrEqual;
	};
	template< unsigned int Value1 , unsigned int Value2 >
	struct Compare< UIntPack< Value1 > , UIntPack< Value2 > >
	{
		static const bool Equal = Value1==Value2;
		static const bool NotEqual = Value1!=Value2;
		static const bool LessThan = Value1<Value2;
		static const bool LessThanOrEqual = Value1<=Value2;
		static const bool GreaterThan = Value1>Value2;
		static const bool GreaterThanOrEqual = Value1>=Value2;
	};
	template<>
	struct Compare< UIntPack<> , UIntPack<> >
	{
		static const bool Equal = true;
		static const bool NotEqual = false;
		static const bool LessThan = false;
		static const bool LessThanOrEqual = true;
		static const bool GreaterThan = false;
		static const bool GreaterThanOrEqual = true;
	};


	///////////////////////////
	// The isotropic variant //
	///////////////////////////
	template< unsigned int Dim , unsigned int Value > struct _IsotropicUIntPack             { typedef typename _IsotropicUIntPack< Dim-1 , Value >::type::template Append< Value > type; };
	template<                    unsigned int Value > struct _IsotropicUIntPack< 1 , Value >{ typedef UIntPack< Value > type; };
	template<                    unsigned int Value > struct _IsotropicUIntPack< 0 , Value >{ typedef UIntPack< > type; };
	template< unsigned int Dim , unsigned int Value > using IsotropicUIntPack = typename _IsotropicUIntPack< Dim , Value >::type;
	template< unsigned int Dim > using ZeroUIntPack = IsotropicUIntPack< Dim , 0 >;

	////////////////////////
	// Sequential variant //
	////////////////////////
	template< unsigned int Dim , unsigned int Value >   struct _SequentialUIntPack             { typedef Concatenation< UIntPack< Value > , typename _SequentialUIntPack< Dim-1 , Value+1 >::type > type; };
	template<                    unsigned int Value >   struct _SequentialUIntPack< 0 , Value >{ typedef UIntPack<> type; };
	template< unsigned int Dim , unsigned int Value=0 > using SequentialUIntPack = typename _SequentialUIntPack< Dim , Value >::type;
}
#endif // UINT_PACK_INCLUDED
