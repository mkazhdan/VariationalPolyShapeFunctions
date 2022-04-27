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
#pragma message( "TODO:" )
#pragma message( "  1. Extract a slice of a function" )
#pragma message( "  2. Determinant, inverse, and adjugate" )

#ifndef AUTO_DIFF_INCLUDED
#define AUTO_DIFF_INCLUDED

#include <iostream>
#include "Tensors.h"
#include "Factorial.h"

namespace Misha
{
	template< unsigned int D , typename OutPack , typename InPack > struct _OutDPack;

	template< unsigned int D , unsigned int ... OutDims , unsigned int ... InDims >
	struct _OutDPack< D , UIntPack< OutDims ... > , UIntPack< InDims ... > >
	{
		typedef UIntPack< InDims ... > InPack;
		typedef UIntPack< OutDims ... > OutPack;
		typedef Concatenation< typename _OutDPack< D-1 , OutPack , InPack >::type , InPack > type;
	};
	template< unsigned int ... OutDims , unsigned int ... InDims >
	struct _OutDPack< 0 , UIntPack< OutDims ... > , UIntPack< InDims ... > >
	{
		typedef UIntPack< OutDims ... > type;
	};
	template< unsigned int D , typename OutPack , typename InPack > using OutDPack = typename _OutDPack< D , OutPack , InPack >::type;

	template< unsigned int D >
	struct ContinuousDifferentiator
	{
		template< typename F >
		static Tensor< OutDPack< D , typename F::OutPack , typename F::InPack > > Derivative( const F &f , const Tensor< typename F::InPack > &t )
		{
			return ContinuousDifferentiator< D-1 >::Derivative( f.d() , t );
		}
		template< typename F > static auto Derivative( const F &f ){ return ContinuousDifferentiator< D-1 >::Derivative( f.d() ); }
	};
	template<>
	struct ContinuousDifferentiator< 0 >
	{
		template< typename F >
		static Tensor< OutDPack< 0 , typename F::OutPack , typename F::InPack > > Derivative( const F &f , const Tensor< typename F::InPack > &t ){ return f(t); }
		template< typename F > static auto Derivative( const F &f ){ return f; }
	};

	template< unsigned int D >
	struct DiscreteDifferentiator
	{
		template< typename F >
		static Tensor< OutDPack< D , typename F::OutPack , typename F::InPack > > Derivative( const F &f , const Tensor< typename F::InPack > &t , double eps )
		{
			typedef OutDPack< D-1 , typename F::OutPack , typename F::InPack > _OutPack;

			Tensor< OutDPack< D , typename F::OutPack , typename F::InPack > > out;

			unsigned int index[ ( _OutPack::Size + F::InPack::Size )>0 ? ( _OutPack::Size + F::InPack::Size ) : 1 ] , idx[ F::InPack::Size>0 ? F::InPack::Size : 1 ];
			WindowLoop< F::InPack::Size >::Run
			(
				ZeroUIntPack< F::InPack::Size >() , F::InPack() ,
				[&]( int d , int i ){ index[ _OutPack::Size+d ] = idx[d] = i; } ,
				[&]( void )
				{
					Tensor< typename F::InPack > delta;
					delta( idx ) = eps;

					Tensor< _OutPack > d = ( DiscreteDifferentiator< D-1 >::template Derivative< F >( f , t+delta , eps ) - DiscreteDifferentiator< D-1 >::template Derivative< F >( f , t-delta , eps ) ) / ( 2*eps );
					WindowLoop< _OutPack::Size >::Run
					(
						ZeroUIntPack< _OutPack::Size >() , _OutPack() ,
						[&]( int d , int i ){ index[d] = i; } ,
						[&]( double v ){ out( index ) = v; } ,
						d
					);
				}
			);
			return out;
		}
	};

	template<>
	struct DiscreteDifferentiator< 0 >
	{
		template< typename F >
		static Tensor< OutDPack< 0 , typename F::OutPack , typename F::InPack > > Derivative( const F &f , const Tensor< typename F::InPack > &t , double eps ){ return f( t ); }
	};

	template< typename OutPack , typename InPack , typename F > struct Function;
	template< typename F1 , typename F2 > struct Composition;

	template< unsigned int ... OutDims , unsigned int ... InDims , typename F >
	struct Function< UIntPack < OutDims ... > , UIntPack< InDims ... > , F >
	{
		// Type of the tensor taken as input
		typedef UIntPack< InDims ... > InPack;
		// Type of the tensor returned as output
		typedef UIntPack< OutDims ... > OutPack;

		// These two methods are expected to be supported by derived classes

		// Gives the D-th derivative
		// [NOTE] We can't actually declare this one as virtual because it is templated.
		// template< unsigned int D >
		// virtual Tensor< OutDPack< D , OutPack , InPack > > d( const Tensor< InPack > &t ) const = 0;

		template< typename _F >
		typename std::enable_if< !std::is_arithmetic_v< _F > , struct Composition< F , _F > >::type operator()( const _F &_f ) const { return Composition< F , _F >( static_cast< const F & >( *this ) , _f ); }
		template< typename V >
		typename std::enable_if< std::is_arithmetic_v< V > && InPack::Size==0 , Tensor< OutPack > >::type operator()( const V &v ) const { return static_cast< const F & >( *this ).template d<0>( Tensor< UIntPack<> >( v ) ); }
		Tensor< OutPack > operator()( const Tensor< InPack > &p ) const { return static_cast< const F & >( *this ).template d<0>( p ); }
	};

	template< typename F >
	struct Scale : public Function< typename F::OutPack , typename F::InPack , Scale< F > >
	{
		typedef Function< typename F::OutPack , typename F::InPack , Scale > _Function;
		template< typename _F > friend struct Scale;

		Scale( const F &f , double s ) : _f(f) , _s(s) {}

		template< unsigned int D >
		Tensor< OutDPack< D , typename _Function::OutPack , typename _Function::InPack > > d( const Tensor< typename _Function::InPack > &t ) const	{ return _f.template d<D>(t) * _s; }

		auto d( void ) const { return Scale< decltype( std::declval< F >().d() ) >( _f.d() , _s ); }

		friend std::ostream &operator <<( std::ostream &os , const Scale &scale )
		{
			if( scale._s==-1 ) return os << "( -" << scale._f << " )";
			else               return os << "( " << scale._s << " * " << scale._f << " )";
		}

	protected:
		const F _f;
		const double _s;
	};

	template< typename F > struct Scale< Scale< F > > : public Scale< F >{ Scale( const Scale< F > &f , double s ) : Scale< F >( f._f , s*f._s ){} };

	template< typename F1 , typename F2 >
	struct Add : public Function< typename F1::OutPack , typename F1::InPack , Add< F1 , F2 > >
	{
		static_assert( Compare< typename F1::OutPack , typename F2::OutPack >::Equal , "[ERROR] Output types differ" );
		static_assert( Compare< typename F2::InPack , typename F2::InPack >::Equal , "[ERROR] Input types differ" );

		typedef Function< typename F1::OutPack , typename F1::InPack , Add > _Function;

		Add( const F1 &f1 , const F2 &f2 ) : _f1(f1) , _f2(f2) {}

		template< unsigned int D >
		Tensor< OutDPack< D , typename _Function::OutPack , typename _Function::InPack > > d( const Tensor< typename _Function::InPack > &t ) const { return _f1.template d<D>(t) + _f2.template d<D>(t); }

		auto d( void ) const { return Add< decltype( std::declval< F1 >().d() ) , decltype( std::declval< F2 >().d() ) >( _f1.d() , _f2.d() ); }

		friend std::ostream &operator <<( std::ostream &os , const Add &add ){ return os << "( " << add._f1 << " + " << add._f2 << " )"; }
	protected:
		const F1 _f1;
		const F2 _f2;
	};

	template< unsigned int I , typename F1 , typename F2 >
	struct ContractedOuterProduct : public Function< Concatenation< typename Split< F1::OutPack::Size-I , typename F1::OutPack >::First , typename Split< I , typename F2::OutPack >::Second > , typename F1::InPack , ContractedOuterProduct< I , F1 , F2 > >
	{
		static_assert( Compare< typename F2::InPack , typename F2::InPack >::Equal , "[ERROR] Input types differ" );
		typedef typename F1::OutPack OutPack1;
		typedef typename F2::OutPack OutPack2;

		typedef Function< Concatenation< typename Split< F1::OutPack::Size-I , typename F1::OutPack >::First , typename Split< I , typename F2::OutPack >::Second > , typename F1::InPack , ContractedOuterProduct > _Function;

		ContractedOuterProduct( const F1 &f1 , const F2 &f2 ) : _f1(f1) , _f2(f2) {}

		template< unsigned int D >
		Tensor< OutDPack< D , typename _Function::OutPack , typename _Function::InPack > > d( const Tensor< typename _Function::InPack > &t ) const
		{
			return _d< D >( t );
		}
		friend std::ostream &operator <<( std::ostream &os , const ContractedOuterProduct &op )
		{
			if( I==0 ) return os << "( " << op._f1 << " * " << op._f2 << " )";
			else       return os << "( " << op._f1 << " *_" << I << " " << op._f2 << " )";
		}

		auto d( void ) const
		{
			return Add
				<
				ContractedOuterProduct< I , decltype( std::declval< F1 >().d() ) , F2 > ,
				ContractedOuterProduct< I , F1 , decltype( std::declval< F2 >().d() ) >
				>
				(
					ContractedOuterProduct< I , decltype( std::declval< F1 >().d() ) , F2 >( _f1.d() , _f2 ) ,
					ContractedOuterProduct< I , F1 , decltype( std::declval< F2 >().d() ) >( _f1 , _f2.d() )
				);
		}

	protected:
		const F1 _f1;
		const F2 _f2;

		template< unsigned int D , unsigned int J=D >
		Tensor< OutDPack< D , typename _Function::OutPack , typename _Function::InPack > > _d( const Tensor< typename _Function::InPack > &t ) const
		{
			if constexpr( J!=0 )
			{
				typedef SequentialUIntPack< OutPack1::Size , 0 > P1;
				typedef SequentialUIntPack< J*_Function::InPack::Size , P1::Size > P2;
				typedef SequentialUIntPack< OutPack2::Size + (D-J)*_Function::InPack::Size , P1::Size + P2::Size > P3;

				return ( _f1.template d< J >(t).template contractedOuterProduct<I>( _f2.template d< D-J >(t) ) ).permute( Concatenation< P1 , P3 , P2 >() ) * (double)Choose< D , J >::Value + _d< D , J-1 >( t );
			}
			else return _f1.template d< J >(t).template contractedOuterProduct< I >( _f2.template d< D-J >(t) );
		}
	};

	template< typename F1 , typename F2 >
	struct Composition : public Function< typename F1::OutPack , typename F2::InPack , Composition< F1 , F2 > >
	{
		static_assert( Compare< typename F1::InPack , typename F2::OutPack >::Equal , "[ERROR] Input/Output types differ" );

		typedef Function< typename F1::OutPack , typename F2::InPack , Composition > _Function;

		Composition( const F1 &f1 , const F2 &f2 ) : _f1(f1) , _f2(f2) {}

		template< unsigned int D >
		Tensor< OutDPack< D , typename _Function::OutPack , typename _Function::InPack > > d( const Tensor< typename _Function::InPack > &t ) const
		{
			if constexpr( D==0 ) return _f1( _f2( t ) );
			else return _d< D , UIntPack< 1 > >( t );
		}

		auto d( void ) const
		{
			typedef decltype( std::declval< F1 >().d() ) DF1;
			typedef decltype( std::declval< F2 >().d() ) DF2;
			return ContractedOuterProduct< F1::InPack::Size , Composition< DF1 , F2 > , DF2 >
				(
					ContractedOuterProduct< F1::InPack::Size , Composition< DF1 , F2 > , DF2 >( Composition< DF1 , F2 >( _f1.d() , _f2 ) , _f2.d() )
				);
		}

		friend std::ostream &operator <<( std::ostream &os , const Composition &composition ){ return os << composition._f1 << "( " << composition._f2 << " )"; }
	protected:

		template< typename JPack >
		static constexpr unsigned int _OutPower( void )
		{
			if constexpr( JPack::Size==1 ) return JPack::First;
			else return JPack::First + _OutPower< JPack::Rest >();
		}
		template< typename JPack , unsigned int DStart=1 >
		static constexpr unsigned int _InPower( void )
		{
			if constexpr( JPack::Size==1 ) return JPack::First * DStart;
			else return JPack::First * DStart + _InPower< JPack::Rest , DStart+1 >();
		}

		template< unsigned int D , typename JPack , unsigned int Index=(unsigned int)-1 >
		Tensor< OutDPack< D , typename _Function::OutPack , typename _Function::InPack > > _d( const Tensor< typename _Function::InPack > &t ) const
		{
			if constexpr( D==_InPower< JPack >() )
			{
				static const unsigned int I = _OutPower< JPack >();
				Tensor< Concatenation< typename F1::OutPack , Power< typename F1::InPack , I > > > t1 = _f1.template d< I >( _f2( t ) );
				Tensor< Concatenation< Power< typename F2::OutPack , _OutPower< JPack >() > , Power< typename F2::InPack , _InPower< JPack >() > > > t2 = _multiplyPower2< JPack >( t );
				return t1.template contractedOuterProduct< F2::OutPack::Size * I >( t2 );
			}
			else
			{
				if      constexpr( Index==-1 ) return _d< D , Concatenation< UIntPack< JPack::First+1 > , typename JPack::Rest > >( t ) + _d< D , JPack , Index+1 >( t );
				else if constexpr( Index<JPack::Size-1 )
				{
					if constexpr( Select< Index , JPack >::Value )
					{
						typedef Concatenation< typename Split< Index , JPack >::First , UIntPack< Select< Index , JPack >::Value-1 > , UIntPack< Select< Index+1 , JPack >::Value+1 > , typename Split< Index+2 , JPack >::Second > _JPack;
						return _d< D , _JPack >( t ) * (double)Select< Index , JPack >::Value + _d< D , JPack , Index+1 >( t );
					}
					else return _d< D , JPack , Index+1 >( t );
				}
				else if constexpr( Index==JPack::Size-1 )
				{
					if constexpr( Select< Index , JPack >::Value!=0 )
					{
						typedef Concatenation< typename Split< Index , JPack >::First , UIntPack< Select< Index , JPack >::Value-1 > , UIntPack< 1 > > _JPack;
						return _d< D , _JPack >( t ) * (double)Select< Index , JPack >::Value;
					}
					else return Tensor< OutDPack< D , typename _Function::OutPack , typename _Function::InPack > >();
				}
			}
		}

		// Comptes the D-th derivate, raised to the P-th power (and reorders)
		template< unsigned int D , unsigned int P >
		Tensor< Concatenation< Power< typename F2::OutPack , P > , Power< typename F2::InPack , D*P > > > _power2( const Tensor< typename F2::InPack > &t ) const
		{
			if      constexpr( P==0 ) return Tensor< UIntPack<> >(1);
			else if constexpr( P==1 ) return _f2.template d< D >( t );
			else
			{
				typedef SequentialUIntPack< F2::OutPack::Size*(P-1)   , 0 > P1;
				typedef SequentialUIntPack< F2:: InPack::Size*(P-1)*D , P1::Size > P2;
				typedef SequentialUIntPack< F2::OutPack::Size         , P1::Size + P2::Size > P3;
				typedef SequentialUIntPack< F2:: InPack::Size*D       , P1::Size + P2::Size + P3::Size > P4;

				Tensor< Concatenation< Power< typename F2::OutPack , P-1 > , Power< typename F2::InPack , D*(P-1) > > > t1 = _power2< D , P-1 >( t );
				Tensor< Concatenation< typename F2::OutPack , Power< typename F2::InPack , D > > > t2 = _f2.template d< D >( t );
				return ( t1 * t2 ).permute( Concatenation< P1 , P3 , P2 , P4 >() );
			}
		}

		// Multiplies a product of derivatives by the D-th derivative, raised to the P-th power (and reorders)
		template< unsigned int D , unsigned int P , unsigned int POut , unsigned int PIn >
		Tensor< Concatenation< Power< typename F2::OutPack , POut + P > , Power< typename F2::InPack , PIn + D*P > > > _leftMultiply2( const Tensor< Concatenation< Power< typename F2::OutPack , POut > , Power< typename F2::InPack , PIn > > > &t2 , const Tensor< typename F2::InPack > &t ) const
		{
			typedef SequentialUIntPack< F2::OutPack::Size*P    , 0 > P1;
			typedef SequentialUIntPack< F2:: InPack::Size*D*P  , P1::Size > P2;
			typedef SequentialUIntPack< F2::OutPack::Size*POut , P1::Size+P2::Size > P3;
			typedef SequentialUIntPack< F2:: InPack::Size*PIn  , P1::Size + P2::Size + P3::Size > P4;

			Tensor< Concatenation< Power< typename F2::OutPack , P > , Power< typename F2::InPack , D*P > > > t1 = _power2< D , P >( t );
			return ( t1 * t2 ).permute( Concatenation< P1 , P3 , P2 , P4 >() );
		}

		// Computes the product of the powers of the derivatives
		template< typename JPack , unsigned int DStart=1 >
		Tensor< Concatenation< Power< typename F2::OutPack , _OutPower< JPack >() > , Power< typename F2::InPack , _InPower< JPack , DStart >() > > > _multiplyPower2( const Tensor< typename F2::InPack > &t ) const
		{
			if constexpr( JPack::Size==1 ) return _power2< DStart , JPack::First >( t );
			// Assuming JPack = < 0 , 1 > and DStart = 1
			// F2::OutPack = <3>
			// F2::InPack = <3>
			// F1::OutPack = <>
			// F1::InPack = <3>
			// The output of:
			//		_leftMultiply2< D=DStart , P=JPack::First , POut=_OutPower< JPack::Rest > , PIn=_InPower< JPack::Rest > >
			//		_leftMultiply2< D=1 , P=0 , POut=_OutPower< < 1 > > , PIn=_InPower< < 1 > > >
			//		_leftMultiply2< D=1 , P=0 , POut=1 , PIn=1 >
			// wil be of type: 
			//		Tensor< Concatenation< Power< typename F2::OutPack , POut + P > , Power< typename F2::InPack , PIn + D*P > > >
			//		Tensor< Concatenation< Power< <3> , 1 > , Power< <3> , 1 > > >
			//		Tensor< Concatenation< <>  , <3> > >
			//		Tensor< < 3 > >
			else return _leftMultiply2< DStart , JPack::First , _OutPower< JPack::Rest >() , _InPower< JPack::Rest , DStart+1 >() >( _multiplyPower2< JPack::Rest , DStart+1 >( t ) , t );
		}

		const F1 _f1;
		const F2 _f2;
	};

	template< typename SquareMatrixF >
	struct Adjugate
	{
		// Need to define this functionality
	};
	template< typename SquareMatrixF >
	struct Determinant
	{
		// Need to define this functionality
	};

	template< typename OutPack , typename InPack > struct Constant;
	template< typename OutPack , typename InPack > struct Linear;

	template< unsigned int ... OutDims , unsigned int ... InDims >
	struct Constant< UIntPack< OutDims ... > , UIntPack< InDims ... > > : public Function< UIntPack< OutDims ... > , UIntPack< InDims ... > , Constant< UIntPack< OutDims ... > , UIntPack< InDims ... > > >
	{
		Constant( void ){}
		Constant( const Tensor< UIntPack< OutDims ... > > &c ) : _c(c){}
		template< unsigned int D >
		Tensor< OutDPack< D , UIntPack< OutDims ... > , UIntPack< InDims ... > > > d( const Tensor< UIntPack< InDims ... > > &t ) const { return _d<D>(t); }
		Constant< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > , UIntPack< InDims ... > > d( void ) const { return Constant< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > , UIntPack< InDims ... > >(); }
		friend std::ostream &operator <<( std::ostream &os , const Constant &c ){ return os << c._c; }
	protected:
		const Tensor< UIntPack< OutDims ... > > _c;

		template< unsigned int D >
		Tensor< OutDPack< D , UIntPack< OutDims ... > , UIntPack< InDims ... > > >_d( const Tensor< UIntPack< InDims ... > > & ) const
		{
			if constexpr( D==0 ) return _c;
			else return Tensor< OutDPack< D , UIntPack< OutDims ... > , UIntPack< InDims ... > > >();
		}
	};

	template< unsigned int ... OutDims , unsigned int ... InDims >
	struct Linear< UIntPack< OutDims ... > , UIntPack< InDims ... > > : public Function< UIntPack< OutDims ... > , UIntPack< InDims ... > , Linear< UIntPack< OutDims ... > , UIntPack< InDims ... > > >
	{
		Linear( void ){}
		Linear( const Tensor< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > > &l ) : _l(l){}
		template< unsigned int D >
		Tensor< OutDPack< D , UIntPack< OutDims ... > , UIntPack< InDims ... > > > d( const Tensor< UIntPack< InDims ... > > &t ) const { return _d<D>(t); }
		Constant< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > , UIntPack< InDims ... > > d( void ) const { return Constant< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > , UIntPack< InDims ... > >(_l); }
		friend std::ostream &operator <<( std::ostream &os , const Linear &l ){ return os << l._l; }
	protected:
		Tensor< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > > _l;

		template< unsigned int D >
		Tensor< OutDPack< D , UIntPack< OutDims ... > , UIntPack< InDims ... > > > _d( const Tensor< UIntPack< InDims ... > > &t ) const
		{
			if      constexpr( D==0 ) return _l.template contractedOuterProduct< sizeof ... ( InDims ) >( t );
			else if constexpr( D==1 ) return _l;
			else return Tensor< OutDPack< D , UIntPack< OutDims ... > , UIntPack< InDims ... > > >();
		}
	};

	template< typename F1 , typename F2 >
	Add< F1 , F2 > operator + ( const F1 &f1 , const F2 &f2 ){ return Add< F1 , F2 >(f1,f2); }

#if 1
	// The things I need to do to keep clang happy
	template< typename F >
	typename std::enable_if< std::is_class< F >::value , Scale< F > >::type operator * ( const F &f , double s ){ return Scale< F >(f,s); }

	template< typename F >
	typename std::enable_if< std::is_class< F >::value , Scale< F > >::type operator * ( const F &f , int s ){ return Scale< F >(f,s); }

	template< typename F >
	typename std::enable_if< std::is_class< F >::value , Scale< F > >::type operator * ( double s , const F &f ){ return Scale< F >(f,s); }

	template< typename F >
	typename std::enable_if< std::is_class< F >::value , Scale< F > >::type operator * ( int s , const F &f ){ return Scale< F >(f,s); }

	template< typename F >
	typename std::enable_if< std::is_class< F >::value , Scale< F > >::type operator - ( const F &f ){ return f * -1.; }

	template< typename F >
	typename std::enable_if< std::is_class< F >::value , Scale< F > >::type operator / ( const F &f , double s ){ return Scale< F >(f,1./s); }
#else
	template< typename F >
	Scale< F > operator * ( const F &f , double s ){ return Scale< F >(f,s); }

	template< typename F >
	Scale< F > operator * ( const F &f , int s ){ return Scale< F >(f,s); }

	template< typename F >
	Scale< F > operator * ( double s , const F &f ){ return Scale< F >(f,s); }

	template< typename F >
	Scale< F > operator * ( int s , const F &f ){ return Scale< F >(f,s); }

	template< typename F >
	Scale< F > operator - ( const F &f ){ return f * -1.; }

	template< typename F >
	Scale< F > operator / ( const F &f , double s ){ return Scale< F >(f,1./s); }
#endif 
	template< typename F1 , typename F2 >
	Add< F1 , Scale< F2 > > operator - ( const F1 &f1 , const F2 &f2 ){ return f1 + ( -f2 ); }


	template< typename F1 , typename F2 >
	ContractedOuterProduct< 0 , F1 , F2 > operator * ( const F1 &f1 , const F2 &f2 ){ return ContractedOuterProduct< 0 , F1 , F2 >(f1,f2); }

	template< typename F >
	auto operator + ( const F &f , const Tensor< typename F::OutPack > &t )
	{
		return Add< F , Constant< typename F::OutPack , typename F::InPack > >( f , Constant< typename F::OutPack , typename F::InPack >( t ) );
	}

	template< typename F >
	auto operator + ( const Tensor< typename F::OutPack > &t , const F &f )
	{
		return Add< Constant< typename F::OutPack , typename F::InPack > , F >( Constant< typename F::OutPack , typename F::InPack >( t ) , f );
	}

	template< typename F >
	auto operator - ( const F &f , const Tensor< typename F::OutPack > &t )
	{
		return Add< F , Constant< typename F::OutPack , typename F::InPack > >( f , Constant< typename F::OutPack , typename F::InPack >( -t ) );
	}

	template< typename F >
	auto operator - ( const Tensor< typename F::OutPack > &t , const F &f )
	{
		return Add< Constant< typename F::OutPack , typename F::InPack > , Scale< F > >( Constant< typename F::OutPack , typename F::InPack >( t ) , -f );
	}

	template< typename F >
	auto operator + ( const F &f , double s )
	{
		static_assert( Compare< typename F::OutPack , UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
		return Add< F , Constant< typename F::OutPack , typename F::InPack > >( f , Constant< typename F::OutPack , typename F::InPack >( s ) );
	}

	template< typename F >
	auto operator + ( double s , const F &f )
	{
		static_assert( Compare< typename F::OutPack , UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
		return Add< Constant< typename F::OutPack , typename F::InPack > , F >( Constant< typename F::OutPack , typename F::InPack >( s ) , f );
	}

	template< typename F >
	auto operator - ( const F &f , double s )
	{
		static_assert( Compare< typename F::OutPack , UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
		return Add< F , Constant< typename F::OutPack , typename F::InPack > >( f , Constant< typename F::OutPack , typename F::InPack >( -s ) );
	}

	template< typename F >
	auto operator - ( double s , const F &f )
	{
		static_assert( Compare< typename F::OutPack , UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
		return Add< Constant< typename F::OutPack , typename F::InPack > , Scale< F > >( Constant< typename F::OutPack , typename F::InPack >( s ) , -f );
	}

	struct _Pow;
	struct _Exp;
	struct _Log;
	struct _Sin;
	struct _Cos;
	template< typename Pack > struct Identity;
	template< typename Pack > struct SquareNorm;

	struct _Pow : public Function< UIntPack<> , UIntPack<> , _Pow >
	{
		_Pow( double e ) : _e(e) {}
		template< unsigned int D >
		Tensor< UIntPack<> > d( const Tensor< UIntPack<> > &t ) const { return Tensor< UIntPack<> >( pow( t() , _e-D ) * _Scalar<D>( _e ) ); }
		Scale< _Pow > d( void ) const { return _Pow( _e-1 ) * _e; }
		friend std::ostream &operator<<( std::ostream &os , const _Pow &p ){ return os << "pow( " << p._e << " )"; }
	protected:
		double _e;

		template< unsigned int D >
		static double _Scalar( double e )
		{
			if constexpr( D==0 ) return e;
			else return e * _Scalar< D-1 >( e-1 );
		}
	};

	template< typename F >
	auto Pow( const F &f , double e ){ return Composition< _Pow , F >( _Pow(e) , f ); }

	template< typename F >
	auto operator / ( double n , const F &f ){ return n * Pow(f,-1.); }

	struct _Exp : public Function< UIntPack<> , UIntPack<> , _Exp >
	{
		template< unsigned int D >
		Tensor< UIntPack<> > d( const Tensor< UIntPack<> > &t ) const { return Tensor< UIntPack<> >( exp( t() ) ); }
		_Exp d( void ) const { return _Exp(); }
		friend std::ostream &operator<<( std::ostream &os , const _Exp & ){ return os << "exp"; }
	};

	template< typename F >
	auto Exp( const F &f ){ return Composition< _Exp , F >( _Exp() , f ); }

	struct _Log : public Function< UIntPack<> , UIntPack<> , _Log >
	{
		template< unsigned int D >
		Tensor< UIntPack<> > d( const Tensor< UIntPack<> > &t ) const
		{
			if constexpr( D==0 ) return Tensor< UIntPack<> >( log( t() ) );
			else return _Pow( -1. ).template d< D-1  >( t );
		}
		_Pow d( void ) const { return _Pow(-1.); }
		friend std::ostream &operator<<( std::ostream &os , const _Log & ){ return os << "log"; }
	};
	template< typename F >
	auto Log( const F &f ){ return Composition< _Log , F >( _Log() , f ); }


	struct _Sin : public Function< UIntPack<> , UIntPack<> , _Sin >
	{
		template< unsigned int D >
		Tensor< UIntPack<> > d( const Tensor< UIntPack<> > &t ) const
		{
			if( D&1 )
				if( (D>>1)&1 ) return Tensor< UIntPack<> >( -cos( t() ) );
				else           return Tensor< UIntPack<> >(  cos( t() ) );
			else
				if( (D>>1)&1 ) return Tensor< UIntPack<> >( -sin( t() ) );
				else           return Tensor< UIntPack<> >(  sin( t() ) );
		}
		struct _Cos d( void ) const;
		friend std::ostream &operator<<( std::ostream &os , const _Sin & ){ return os << "sin"; }
	};
	template< typename F >
	auto Sin( const F &f ){ return Composition< _Sin , F >( _Sin() , f ); }

	struct _Cos : public Function< UIntPack<> , UIntPack<> , _Cos >
	{
		template< unsigned int D >
		Tensor< UIntPack<> > d( const Tensor< UIntPack<> > &t ) const
		{
			if( D&1 )
				if( (D>>1)&1 ) return Tensor< UIntPack<> >(  sin( t() ) );
				else           return Tensor< UIntPack<> >( -sin( t() ) );
			else
				if( (D>>1)&1 ) return Tensor< UIntPack<> >( -cos( t() ) );
				else           return Tensor< UIntPack<> >(  cos( t() ) );
		}
		Scale< _Sin > d( void ) const { return -_Sin(); }
		friend std::ostream &operator <<( std::ostream &os , const _Cos & ){ return os << "cos"; }
	};
	template< typename F >
	auto Cos( const F &f ){ return Composition< _Cos , F >( _Cos() , f ); }

	_Cos _Sin::d( void ) const { return _Cos(); }


	template< unsigned int ... Dims >
	struct Identity< UIntPack< Dims ... > > : public Linear< UIntPack< Dims ... > , UIntPack< Dims ... > >
	{
		typedef UIntPack< Dims ... > Pack;
		using Linear< Pack , Pack >::_l;
		Identity( void )
		{
			unsigned int index[ 2*Pack::Size ];
			WindowLoop< Pack::Size >::Run
			(
				ZeroUIntPack< Pack::Size >() , Pack() ,
				[&]( int d , int i ){ index[d] = index[d+Pack::Size] = i; } ,
				[&]( void ){ _l( index ) = 1; }
			);
		}
	};

	template< unsigned int ... Dims >
	struct SquareNorm< UIntPack< Dims ... > > : public Function< UIntPack<> , UIntPack< Dims ... > , SquareNorm< UIntPack< Dims ... > > >
	{
		typedef UIntPack< Dims ... > Pack;
		template< unsigned int D >
		Tensor< OutDPack< D , UIntPack<> , Pack > > d( const Tensor< Pack > &t ) const { return _d< D >(t); }
		auto d( void ) const { return Identity< UIntPack< Dims ... > >() * 2.; }
		friend std::ostream &operator <<( std::ostream &os , const SquareNorm & ){ return os << "SquareNorm[" << Pack() << "]"; }
	protected:
		template< unsigned int D >
		Tensor< OutDPack< D , UIntPack<> , Pack > > _d( const Tensor< Pack > &t ) const
		{
			Tensor< OutDPack< D , UIntPack<> , Pack > > d;
			if constexpr( D==0 )
			{
				WindowLoop< Pack::Size >::Run
				(
					ZeroUIntPack< Pack::Size >() , Pack() ,
					[]( int d , int i ){} ,
					[&]( double v ){ d() += v*v; } ,
					t
				);
			}
			else if constexpr( D==1 )
			{
				WindowLoop< Pack::Size >::Run
				(
					ZeroUIntPack< Pack::Size >() , Pack() ,
					[]( int d , int i ){} ,
					[]( double &d , double v ){ d = 2*v; } ,
					d , t
				);
			}
			else if constexpr( D==2 )
			{
				unsigned int index[ 2*Pack::Size ];
				WindowLoop< Pack::Size >::Run
				(
					ZeroUIntPack< Pack::Size >() , Pack() ,
					[&]( int d , int i ){ index[d] = index[ d+Pack::Size ] = i; } ,
					[&]( void ){ d( index ) = 2.; }
				);
			}
			return d;
		}
	};

}
#endif // AUTO_DIFF_INCLUDED
