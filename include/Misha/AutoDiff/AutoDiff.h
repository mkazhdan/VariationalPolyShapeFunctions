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
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMIRTED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef AUTO_DIFF_INCLUDED
#define AUTO_DIFF_INCLUDED

#include <iostream>
#include "Tensors.h"
#include "Misha/Exceptions.h"

namespace AutoDiff
{
	template< typename T > bool constexpr IsScalar( void );

	////////////////////////////
	////////////////////////////
	//// Short declarations ////
	////////////////////////////
	////////////////////////////

	//////////////
	// Function //
	//////////////
	// A general class describing a function taking in a tensor and outputting a tensor
	// (Parameters OutPack and InPack are assumed to be of variadic UIntPack<...> type giving the
	// dimensions along the different ranks of the input/output tensors.)
	template< typename OutPack , typename InPack , typename F > struct Function;

	/////////////////////
	// Basic Functions //
	/////////////////////

	// A Function that is constant in its input
	// InPack: a UIntPacK< ... > describing the dimensions of the input tensor
	// OutPack: a UIntPacK< ... > describing the dimensions of the output tensor
	template< typename OutPack , typename InPack > struct Constant;

	// A Function that is linear in its input
	// InPack: a UIntPacK< ... > describing the dimensions of the input tensor
	// OutPack: a UIntPacK< ... > describing the dimensions of the output tensor
	template< typename OutPack , typename InPack > struct Linear;

	// The identity Function
	// InOutPack: a UIntPacK< ... > describing the dimensions of the input/output tensor
	template< typename InOutPack > struct Identity;


	//////////////////////////
	// Function composition //
	//////////////////////////

	// A function returning the negation of a Function
	// Output type derives from Function< F::OutPack , F::InPack >
	template< typename F > auto operator - ( const F &f );

	// A function returning the sum of two Function's (and specializations)
	// Assumes:
	//		F1::InPack==F2::InPack
	//		F1::OutPack==F2::OutPack
	// Output type derives from Function< F1::OutPack , F1::InPack >
	template< typename F1 , typename F2 > auto operator + ( const F1 &f1 , const F2 &f2 );

	// A function returning the difference of two Function's (and specializations)
	// Assumes:
	//		F1::InPack==F2::InPack
	//		F1::OutPack==F2::OutPack
	// Output type derives from Function< F1::OutPack , F1::InPack >
	template< typename F1 , typename F2 > auto operator - ( const F1 &f1 , const F2 &f2 );

	// A function returning the outer product of a Function with a Function (and specialization)
	// Assumes:
	//		F1::InPack==F2::InPack
	// Output type derives from Function< Concatenation< F1::OutPack , F2::OutPack > , F1::InPack >
	template< typename F1 , typename F2 > auto operator * ( const F1 &f , const F2 &f2 );

	// A function returning the quotient of a Function by a scalar-valued Function
	// Assumes:
	//		F1::InPack==F2::InPack
	//		F2::OutPack==UIntPack<>
	// Output type derives from Function< F1::OutPack , F1::InPack >
	template< typename F1 , typename F2 > auto operator / ( const F1 &f1 , const F2 &f2 );


	// A function returning the permutation a Function's output
	// Assumes:
	//		PermutationPack is a permutation
	//		PermutationPack::Size==F::OutPack::Size
	// Output type derives from Function< Permutation< F::OutPack , PermutationPack > >
	template< typename PermutationPack , typename F > auto Permutation( const F &f );

	// A function returning the contracted outer product of two Function's
	// Assumes:
	//		F1::InPack==F2::InPack
	//		F1::OutPack==Concatenation< Pack1 , Pack2 >
	//		F2::OutPack==Concatenation< Pack2 , Pack3 >
	//		Pack::Size = I
	// Output type derives from Function< Concatenation< Pack1 , Pack3 > , F1::InPack >
	template< unsigned int I , typename F1 , typename F2 > auto ContractedOuterProduct( const F1 &f1 , const F2 &f2 );

	// A function returning the composition two Function's
	// Assumes:
	//		F2::InPack==F1::OutPack
	// Output type derives from Function< F2::OutPack , F1::InPack >
	template< typename F1 , typename F2 > auto Composition( const F1 &f1 , const F2 &f2 );

	// A function extracting a single coefficient of the output
	// Assumes:
	//		sizeof...(indices)==F::OutPack::Size
	//		indices[i]<F::OutPack[i]
	// Output type derives from Function< UIntPack<> , F1::InPack >
	template< typename F > auto Coefficient( const unsigned int indices[] , const F &f );
	template< unsigned int ... CoefficientIndices , typename F > auto Coefficient( UIntPack< CoefficientIndices ... > , const F &f );

	// A function returning the determinant of the output of a Function (assumed to return a square 2-tensor)
	// Assumes:
	//		F::OutPack==UIntPack< Dim , Dim >
	// Output type derives from Function< UIntPack<> , F::InPack >
	template< typename F > auto Determinant( const F &f );

	// A function returning the cofactor of the output of a Function (assumed to return a square 2-tensor)
	// Assumes:
	//		F::OutPack==UIntPack< Dim , Dim >
	// Output type derives from Function< UIntPack<> , F::InPack >
	template< typename F > auto Cofactor( const F &f );

	// A function returning the transpose of the output of a Function
	// Output type derives from Function< F::OutPack::Transpose , F::InPack >
	template< typename F > auto Transpose( const F &f );

	// A function returning the adjugate of the output of a function (assumed to return a square 2-tensor)
	// Assumes:
	//		F::OutPack==UIntPack< Dim , Dim >
	// Output type derives from Function< UIntPack<> , F::InPack >
	template< typename F > auto Adjugate( const F &f );

	// A function returning the inverse of the output of a function (assumed to return a square 2-tensor)
	// Assumes:
	//		F::OutPack==UIntPack< Dim , Dim >
	// Output type derives from Function< UIntPack<> , F::InPack >
	template< typename F > auto Inverse( const F &f );

	// A function returning the cross-product of the columns of the output of a Function (assumed to return a 2-tensor with one more row than columns)
	// Assumes:
	//		F::OutPack==UIntPack< Dim , Dim-1 >
	// Output type derives from Function< UIntPack<Dim> , F::InPack >
	template< typename F > auto CrossProduct( const F &f );

	// A function returning the contraction of a Function
	// Assumes:
	//		F::OutPack::Size>=2
	// Output type derives from Function< F::OutPack::Remove(I1,I2) , F::InPack >
	template< unsigned int I1 , unsigned int I2 , typename F > auto Contract( const F &f );

	// A function returning the square norm of the output of a Function
	// Output type derives from Function< UIntPack<> , F::InPack >
	template< typename F > auto SquareNorm( const F &f );

	// Some common transcendental Function's
	// Assumes:
	//		F::OutPack==UIntPack<>
	// Output type derives from Function< UIntPack<> , F::InPack >
	template< typename F > auto Pow( const F &f , double e );
	template< typename F > auto Exp( const F &f );
	template< typename F > auto Log( const F &f );
	template< typename F > auto Sin( const F &f );
	template< typename F > auto Cos( const F &f );
	template< typename F > auto Sinh( const F &f );
	template< typename F > auto Cosh( const F &f );

	// A function returning the discrete derivative of a Function at a particular input
	template< typename F > auto DiscreteDerivative( const F &f , const Tensor< typename F::InPack > &in , double eps );

	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	//// Full declarations (including helper functions) ////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////

	template< typename T >
	bool constexpr IsScalar( void )
	{
		if constexpr( std::is_arithmetic_v< T > || std::is_base_of< Tensor< UIntPack<> > , T >::value ) return true;
		else return false;
	}

	// Some combinatoric functions that will be useful
	template< unsigned int D > struct Factorial      { static const unsigned int Value = Factorial< D-1 >::Value * D; };
	template<>                 struct Factorial< 0 > { static const unsigned int Value = 1; };

	template< unsigned int D , unsigned int K > struct Choose;
	template< unsigned int D , unsigned int K > struct Choose         { static const unsigned int Value = ( Choose< D-1 , K-1 >::Value * D ) / K; };
	template< unsigned int D >                  struct Choose< D , 0 >{ static const unsigned int Value = 1; };

	// A class for describing a function
	template< unsigned int ... OutDims , unsigned int ... InDims , typename F >
	struct Function< UIntPack < OutDims ... > , UIntPack< InDims ... > , F >
	{
		// Type of the tensor taken as input
		typedef UIntPack< InDims ... > InPack;
		// Type of the tensor returned as output
		typedef UIntPack< OutDims ... > OutPack;

		// [NOTE] There are three ways the function call operator can be invoked with a single argument:
		//	1. [Composition] Where the argument is a function, implying composition
		//	2. [Evaluation] Where the argument is a scalar and the input is a zero-order tensor
		//	3. [Evaluation] Where the argument is a tensor of the same type as the input
		template< typename V > auto operator()( const V &v ) const;
	};


	// A class for describing the product of a function with a scalar
	template< typename F >
	struct _Scale : public Function< typename F::OutPack , typename F::InPack , _Scale< F > >
	{
		typedef Function< typename F::OutPack , typename F::InPack , _Scale > _Function;
		template< typename _F > friend struct _Scale;

		_Scale( const F &f , double s ) : _f(f) , _s(s) {}

		auto value( const Tensor< typename _Function::InPack > &t ) const;
		auto d( void ) const;
		template< typename _F > friend std::ostream &operator << ( std::ostream &os , const _Scale< _F > &scale );
		const F &f( void ) const { return f; }
	protected:
		template< typename _F > friend auto operator * ( const _Scale< _F > & , double );
		template< typename _F > friend auto operator * ( double , const _Scale< _F > & );
		const F _f;
		const double _s;
	};

	template< bool B , bool ... Bs >
	constexpr bool AND( void )
	{
		if constexpr( sizeof...(Bs)==0 ) return B;
		else return B && AND< Bs... >();
	}
	template< bool B , bool ... Bs >
	constexpr bool OR( void )
	{
		if constexpr( sizeof...(Bs)==0 ) return B;
		else return B || OR< Bs... >();
	}

	// A class for describing the sum of two or more functions (with the same order input and the same order output)
	template< typename ... Fs > struct _Add;

	template< typename F , typename ... Fs >
	struct _Add< F , Fs ... > : public Function< typename F::OutPack , typename F::InPack , _Add< F , Fs ... > >
	{
		static_assert( AND< Compare< typename F::OutPack , typename Fs::OutPack >::Equal ... >() , "[ERROR] Output types differ" );
		static_assert( AND< Compare< typename F:: InPack , typename Fs:: InPack >::Equal ... >() , "[ERROR] Input types differ" );

		typedef Function< typename F::OutPack , typename F::InPack , _Add > _Function;

		_Add( const std::tuple< F , Fs ... > f ) : _f(f) {}

		auto value( const Tensor< typename _Function::InPack > &t ) const;
		auto d( void ) const;
		template< typename _F , typename ... _Fs >
		friend std::ostream &operator << ( std::ostream &os , const _Add< _F , _Fs... > &_Add );
		const std::tuple< F , Fs ... > &f_tuple( void ) const { return _f; }
	protected:
		template< unsigned int I > void _toStream( std::ostream &os ) const;
		template< unsigned int I > auto _d( void ) const;
		template< unsigned int I > auto _value( const Tensor< typename _Function::InPack > &t ) const;

		const std::tuple< F , Fs... > _f;
	};

	template< typename PermutationPack , typename F > struct _Permutation;

	template< unsigned int ... PermutationIndices , typename F >
	struct _Permutation< UIntPack< PermutationIndices ... > , F > : public Function< PermutedUIntPack< typename F::OutPack , UIntPack< PermutationIndices ... > > , typename F::InPack , _Permutation< UIntPack< PermutationIndices ... > , F > >
	{
		typedef UIntPack< PermutationIndices ... > PermutationPack;
		static_assert( PermutationPack::Size==F::OutPack::Size , "[ERROR] Sizes don't match" );
		typedef Function< PermutedUIntPack< typename F::OutPack , PermutationPack > , typename F::InPack , _Permutation< PermutationPack , F > > _Function;

		_Permutation( const F &f ) : _f(f) {}
		auto value( const Tensor< typename _Function::InPack > &t ) const;
		auto d( void ) const;
		template< typename _PermutationPack , typename _F > friend std::ostream &operator << ( std::ostream & , const _Permutation< _PermutationPack , _F > & );

	protected:
		const F _f;
	};

	// A class for describing the product of two functions (with the same order input)
	template< unsigned int I , typename F1 , typename F2 >
	struct _ContractedOuterProduct : public Function< Concatenation< typename Split< F1::OutPack::Size-I , typename F1::OutPack >::First , typename Split< I , typename F2::OutPack >::Second > , typename F1::InPack , _ContractedOuterProduct< I , F1 , F2 > >
	{
		static_assert( Compare< typename F2::InPack , typename F2::InPack >::Equal , "[ERROR] Input types differ" );
		typedef typename F1::OutPack OutPack1;
		typedef typename F2::OutPack OutPack2;

		typedef Function< Concatenation< typename Split< F1::OutPack::Size-I , typename F1::OutPack >::First , typename Split< I , typename F2::OutPack >::Second > , typename F1::InPack , _ContractedOuterProduct > _Function;

		_ContractedOuterProduct( const F1 &f1 , const F2 &f2 ) : _f1(f1) , _f2(f2) {}

		auto value( const Tensor< typename _Function::InPack > &t ) const;
		auto d( void ) const;
		template< unsigned int _I , typename _F1 , typename _F2 > friend std::ostream &operator << ( std::ostream &os , const _ContractedOuterProduct< _I , _F1 , _F2 > &op );

	protected:
		const F1 _f1;
		const F2 _f2;
	};


	// A class for describing the composition of two functions
	template< typename F1 , typename F2 >
	struct _Composition : public Function< typename F1::OutPack , typename F2::InPack , _Composition< F1 , F2 > >
	{
		static_assert( Compare< typename F1::InPack , typename F2::OutPack >::Equal , "[ERROR] Input/Output types differ" );

		typedef Function< typename F1::OutPack , typename F2::InPack , _Composition > _Function;

		_Composition( const F1 &f1 , const F2 &f2 ) : _f1(f1) , _f2(f2) {}

		auto value( const Tensor< typename _Function::InPack > &t ) const;
		auto d( void ) const;
		template< typename _F1 , typename _F2 >
		friend std::ostream &operator << ( std::ostream &os , const _Composition< _F1 , _F2 > &composition );
	protected:
		const F1 _f1;
		const F2 _f2;
	};

	// A class for describing a constant function
	template< unsigned int ... OutDims , unsigned int ... InDims >
	struct Constant< UIntPack< OutDims ... > , UIntPack< InDims ... > > : public Function< UIntPack< OutDims ... > , UIntPack< InDims ... > , Constant< UIntPack< OutDims ... > , UIntPack< InDims ... > > >
	{
		Constant( void ){}
		Constant( const Tensor< UIntPack< OutDims ... > > &c ) : _c(c){}

		auto value( const Tensor< UIntPack< InDims ... > > &t ) const;
		auto d( void ) const;
		template< unsigned int ... _OutDims , unsigned int ... _InDims >
		friend std::ostream &operator << ( std::ostream &os , const Constant< UIntPack< _OutDims ... > , UIntPack< _InDims ... > > &c );
	protected:
		const Tensor< UIntPack< OutDims ... > > _c;
	};

	// A class for describing a linear function
	template< unsigned int ... OutDims , unsigned int ... InDims >
	struct Linear< UIntPack< OutDims ... > , UIntPack< InDims ... > > : public Function< UIntPack< OutDims ... > , UIntPack< InDims ... > , Linear< UIntPack< OutDims ... > , UIntPack< InDims ... > > >
	{
		// A constructor generating the zero linear function
		Linear( void ){}

		// A constructor generating a linear function with the prescribed tensor taking input tensors to output tensors
		Linear( const Tensor< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > > &l ) : _l(l){}

		// A constructor generating a simple Linear function with one in the entry {out,in} and zero for all other entries.
		Linear( std::initializer_list< unsigned int > out , std::initializer_list< unsigned int > in );

		// A constructor generating a simple Linear function with one in the entry {_OutDims,_InDims} and zero for all other entries.
		template< unsigned int ... _OutDims , unsigned int ... _InDims > Linear( UIntPack< _OutDims ... > , UIntPack< _InDims ... > );

		auto value( const Tensor< UIntPack< InDims ... > > &t ) const;
		auto d( void ) const;
		template< unsigned int ... _OutDims , unsigned int ... _InDims >
		friend std::ostream &operator << ( std::ostream &os , const Linear< UIntPack< _OutDims ... > , UIntPack< _InDims ... > > &l );

	protected:
		Tensor< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > > _l;
	};

	template< unsigned int ... Dims >
	struct Identity< UIntPack< Dims ... > > : public Linear< UIntPack< Dims ... > , UIntPack< Dims ... > >
	{
		typedef UIntPack< Dims ... > Pack;
		using Linear< Pack , Pack >::_l;
		Identity( void );
	};

	template< typename InPack > struct _Coefficient;

	template< unsigned int ... InDims >
	struct _Coefficient< UIntPack< InDims ... > > : public Linear< UIntPack<> , UIntPack< InDims ... > >
	{
		typedef UIntPack< InDims ... > InPack;
		using Linear< UIntPack<> , InPack >::_l;
		_Coefficient( const unsigned int indices[] );
		template< unsigned int ... CoefficientIndices >
		_Coefficient( UIntPack< CoefficientIndices ... > );
	};

	template< unsigned int I1 , unsigned int I2 , typename InPack > struct _Contraction;

	template< unsigned int I1 , unsigned int I2 , unsigned int ... Dims >
	struct _Contraction< I1 , I2 , UIntPack< Dims ... > > : public Linear< typename Select< I1 , typename Select< I2 , UIntPack< Dims ... > >::Complement >::Complement , UIntPack< Dims ... > >
	{
		typedef UIntPack< Dims ... > InPack;
		typedef typename Select< I1 , typename Select< I2 , InPack >::Complement >::Complement OutPack;
		using Linear< OutPack , InPack >::_l;
		_Contraction( void );
	};

	////////////////////////
	////////////////////////
	//// Implementation ////
	////////////////////////
	////////////////////////

	//////////////
	// Function //
	//////////////
	template< unsigned int ... OutDims , unsigned int ... InDims , typename F >
	template< typename V >
	auto Function< UIntPack < OutDims ... > , UIntPack< InDims ... > , F >::operator()( const V &v ) const
	{
		if constexpr( std::is_base_of< Tensor< InPack > , V >::value ) return static_cast< const F & >( *this ).value( v );
		else if constexpr( std::is_arithmetic_v< V > && InPack::Size==0 ) return static_cast< const F & >( *this ).value( Tensor< UIntPack<> >( v ) );
		else return Composition< F , V >( static_cast< const F & >( *this ) , v );
	}

	////////////
	// _Scale //
	////////////
	template< typename F >
	auto _Scale< F >::value( const Tensor< typename _Function::InPack > &t ) const { return _f(t) * _s; }

	template< typename F >
	auto _Scale< F >::d( void ) const { return _Scale< decltype( std::declval< F >().d() ) >( _f.d() , _s ); }

	template< typename F >
	std::ostream &operator << ( std::ostream &os , const _Scale< F > &scale )
	{
		if( scale._s==-1 ) return os << "( -" << scale._f << " )";
		else               return os << "( " << scale._s << " * " << scale._f << " )";
	}

	//////////
	// _Add //
	//////////
	template< typename F , typename ... Fs >
	auto _Add< F , Fs ... >::value( const Tensor< typename _Function::InPack > &t ) const { return this->template _value<0>(t); }

	template< typename F , typename ... Fs >
	auto _Add< F , Fs ... >::d( void ) const { return this->template _d<0>(); }

	template< typename F , typename ... Fs >
	std::ostream &operator << ( std::ostream &os , const _Add< F , Fs ... > &add )
	{
		os << "( ";
		add.template _toStream<0>( os );
		os << " )";
		return os;
	}

	template< typename F , typename ... Fs >
	template< unsigned int I >
	void _Add< F , Fs ... >::_toStream( std::ostream &os ) const
	{
		if constexpr( I==0 ) os << std::get<I>( _f );
		else os << " + " << std::get<I>( _f );
		if constexpr( I<sizeof...(Fs) ) this->template _toStream< I+1 >( os );
	}

	template< typename F , typename ... Fs >
	template< unsigned int I >
	auto _Add< F , Fs ... >::_d( void ) const
	{
		if constexpr( I==sizeof...(Fs) ) return std::get<I>(_f).d();
		else return std::get<I>(_f).d() + this->template _d<I+1>();
	}

	template< typename F , typename ... Fs >
	template< unsigned int I >
	auto _Add< F , Fs ... >::_value( const Tensor< typename _Function::InPack > &t ) const
	{
		if constexpr( I==sizeof...(Fs) ) return std::get<I>(_f)(t);
		else return std::get<I>(_f)(t) + this->template _value<I+1>(t);
	}

	//////////////
	// Negation //
	//////////////
	template< typename F >
	auto operator - ( const F &f ){ return f * -1.; }

	//////////////
	// Addition //
	//////////////
	template< typename F1 , typename F2 >
	auto operator + ( const F1 &f1 , const F2 &f2 ){ return _Add< F1 , F2 >( std::make_tuple(f1,f2) ); }

	template< typename F1 , typename F2 , typename ... Fs >
	auto operator + ( const F1 &f , const _Add< F2 , Fs ... > &add ){ return _Add< F1 , F2 , Fs ... >( std::tuple_cat( std::make_tuple(f) , add.f_tuple() ) ); }

	template< typename F1 , typename ... Fs , typename F2 >
	auto operator + ( const _Add< F1 , Fs ... > &add , const F2 &f ){ return _Add< F1 , Fs ... , F2 >( std::tuple_cat( add.f_tuple() , std::make_tuple(f) ) ); }

	template< typename F1 , typename ... F1s , typename F2 , typename ... F2s >
	auto operator + ( const _Add< F1 , F1s ... > &add1 , const _Add< F2 , F2s ... > &add2 ){ return _Add< F1 , F1s ... , F2 , F2s ... >( std::tuple_cat( add1.f_tuple() , add2.f_tuple() ) ); }

	template< typename F >
	auto operator + ( const F &f , const Tensor< typename F::OutPack > &t ){ return f + Constant< typename F::OutPack , typename F::InPack >( t ); }

	template< typename F >
	auto operator + ( const Tensor< typename F::OutPack > &t , const F &f ){ return Constant< typename F::OutPack , typename F::InPack >( t ) + f; }

	template< typename F >
	auto operator + ( const F &f , double s )
	{
		static_assert( Compare< typename F::OutPack , UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
		return f + Constant< typename F::OutPack , typename F::InPack >( s );
	}

	template< typename F >
	auto operator + ( double s , const F &f )
	{
		static_assert( Compare< typename F::OutPack , UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
		return Constant< typename F::OutPack , typename F::InPack >( s ) + f;
	}

	/////////////////
	// Subtraction //
	/////////////////
	template< typename F1 , typename F2 >
	auto operator - ( const F1 &f1 , const F2 &f2 ){ return f1 + ( -f2 ); }

	template< typename F >
	auto operator - ( const F &f , const Tensor< typename F::OutPack > &t ){ return f + Constant< typename F::OutPack , typename F::InPack >( -t ); }

	template< typename F >
	auto operator - ( const Tensor< typename F::OutPack > &t , const F &f ){ Constant< typename F::OutPack , typename F::InPack >( t ) - f; }

	template< typename F >
	auto operator - ( const F &f , double s )
	{
		static_assert( Compare< typename F::OutPack , UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
		return f + Constant< typename F::OutPack , typename F::InPack >( -s );
	}

	template< typename F >
	auto operator - ( double s , const F &f )
	{
		static_assert( Compare< typename F::OutPack , UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
		return Constant< typename F::OutPack , typename F::InPack >( s ) - f;
	}

	////////////////////
	// Scalar product //
	////////////////////
	template< typename F >
	auto operator * ( const F &f , const double &s ){ return _Scale< F >(f,s); }

	template< typename F >
	auto operator * ( const _Scale< F > &f , double s ){ return _Scale< F >( f._f , f._s*s ); }

	template< typename F >
	auto operator * ( const double &s , const F &f ){ return f*s; }

	/////////////////////
	// Scalar quotient //
	/////////////////////
	template< typename F >
	auto operator / ( const F &f , const double &s ){ return f * (1./s); }

	////////////////////
	// Tensor product //
	////////////////////
	template< typename F1 , typename F2 >
	auto operator * ( const F1 &f1 , const F2 &f2 ){ return _ContractedOuterProduct< 0 , F1 , F2 >(f1,f2); }

	/////////////////////
	// Tensor quotient //
	/////////////////////
	template< typename F1 , typename F2 >
	auto operator / ( const F1 &f1 , const F2 &f2 )
	{
		static_assert( Compare< typename F2::OutPack , UIntPack<> >::Equal , "[ERROR] denominator is not scalar-valued" );

		if constexpr( std::is_same< F1 , double >::value ) return Constant< UIntPack<> , typename F2::InPack >( f1 );
		else if constexpr( std::is_same< F1 , Tensor< UIntPack<> > >::value ) return Constant< UIntPack<> , typename F2::InPack >( f1 );
		else
		{
			static_assert( Compare< typename F1::InPack , typename F2::InPack >::Equal , "[ERROR] Input types differ" );
			return f1 * Pow( f2 , -1. );
		}
	}

	/////////////////
	// Permutation //
	/////////////////
	template< unsigned int ... PermutationIndices , typename F >
	auto _Permutation< UIntPack< PermutationIndices ... > , F >::value( const Tensor< typename _Function::InPack > &t ) const
	{
		return _f( t ).permute( PermutationPack() );
	}

	template< unsigned int ... PermutationIndices , typename F >
	auto _Permutation< UIntPack< PermutationIndices ... > , F >::d( void ) const
	{
		return Permutation< Concatenation< PermutationPack , SequentialUIntPack< _Function::InPack::Size , PermutationPack::Size > > >( _f.d() ); 
	}

	template< unsigned int ... PermutationIndices , typename F >
	std::ostream &operator << ( std::ostream &os , const _Permutation< UIntPack< PermutationIndices ... > , F > &p )
	{
		os << "Permutation" << _Permutation< UIntPack< PermutationIndices ... > , F >::PermutationPack() << "( " << p._f << " )";
	}

	template< typename PermutationPack , typename F >
	auto Permutation( const F &f ){ return _Permutation< PermutationPack , F >( f ); }

	////////////////////////////
	// ContractedOuterProduct //
	////////////////////////////
	template< unsigned int I , typename F1 , typename F2 >
	auto _ContractedOuterProduct< I , F1 , F2 >::value( const Tensor< typename _Function::InPack > &t ) const
	{
		return _f1(t).template contractedOuterProduct< I >( _f2(t) );
	}

	template< unsigned int I , typename F1 , typename F2 >
	auto _ContractedOuterProduct< I , F1 , F2 >::d( void ) const
	{
		// [Note] As a function is differentiated, the input type is appended to the output type
		// 
		// Get the PermutationPack permuting _f1.d() so that the InPack is on the left
		typedef SequentialUIntPack< F1::InPack::Size , F1::OutPack::Size > PreP1Pack;
		typedef SequentialUIntPack< F1::OutPack::Size , 0 > PreP2Pack;
		typedef Concatenation< PreP1Pack , PreP2Pack > PrePermutationPack;

		// Get the PermutationPack permuting _f1.d() * f2 so that the InPack is on the right
		typedef Concatenation< typename Split< F1::OutPack::Size-I , typename F1::OutPack >::First , typename Split< I , typename F2::OutPack >::Second > ContractionPack;
		typedef SequentialUIntPack< ContractionPack::Size , F1::InPack::Size > PostP1Pack;
		typedef SequentialUIntPack< F1::InPack::Size , 0 > PostP2Pack;
		typedef Concatenation< PostP1Pack , PostP2Pack > PostPermutationPack;

		return Permutation< PostPermutationPack >( ContractedOuterProduct< I >( Permutation< PrePermutationPack >( _f1.d() ) , _f2 ) ) + ContractedOuterProduct< I >( _f1 , _f2.d() );
	}

	template< unsigned int I , typename F1 , typename F2 >
	std::ostream &operator << ( std::ostream &os , const _ContractedOuterProduct< I , F1 , F2 > &op )
	{
		if( I==0 ) return os << "( " << op._f1 << " * " << op._f2 << " )";
		else       return os << "( " << op._f1 << " *_" << I << " " << op._f2 << " )";
	}

	template< unsigned int I , typename F1 , typename F2 >
	auto ContractedOuterProduct( const F1 &f1 , const F2 &f2 ){ return _ContractedOuterProduct< I , F1 , F2 >( f1 , f2 ); }


	/////////////////
	// Composition //
	/////////////////
	template< typename F1 , typename F2 >
	auto _Composition< F1 , F2 >::value( const Tensor< typename _Function::InPack > &t ) const { return _f1( _f2(t) ); }

	template< typename F1 , typename F2 >
	auto _Composition< F1 , F2 >::d( void ) const { return ContractedOuterProduct< F1::InPack::Size >( Composition( _f1.d() , _f2 ) , _f2.d() ); }

	template< typename F1 , typename F2 >
	std::ostream &operator << ( std::ostream &os , const _Composition< F1 , F2 > &composition ){ return os << composition._f1 << "( " << composition._f2 << " )"; }

	template< typename F1 , typename F2 > auto Composition( const F1 &f1 , const F2 &f2 ){ return _Composition< F1 , F2 >( f1 , f2 ); }

	//////////////
	// Constant //
	//////////////
	// A class for reprenting a constant function
	template< unsigned int ... OutDims , unsigned int ... InDims >
	auto Constant< UIntPack< OutDims ... > , UIntPack< InDims ... > >::value( const Tensor< UIntPack< InDims ... > > &t ) const { return _c; }

	template< unsigned int ... OutDims , unsigned int ... InDims >
	auto Constant< UIntPack< OutDims ... > , UIntPack< InDims ... > >::d( void ) const { return Constant< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > , UIntPack< InDims ... > >(); }

	template< unsigned int ... OutDims , unsigned int ... InDims >
	std::ostream &operator << ( std::ostream &os , const Constant< UIntPack< OutDims ... > , UIntPack< InDims ... > > &c ){ return os << c._c; }

	////////////
	// Linear //
	////////////
	// A constructor generating a simple Linear function with one in the entry {out,in} and zero for all other entries.
	template< unsigned int ... OutDims , unsigned int ... InDims >
	Linear< UIntPack< OutDims ... > , UIntPack< InDims ... > >::Linear( std::initializer_list< unsigned int > out , std::initializer_list< unsigned int > in )
	{
		if( out.size()!=sizeof...(OutDims) || in.size()!=sizeof...(InDims) ) ERROR_OUT( "Output dimensions don't match" );
		if constexpr( sizeof...(OutDims)+sizeof...(InDims)==0 ) _l = 1.;
		else
		{
			unsigned int idx[ sizeof...(OutDims)+sizeof...(InDims) ];
			unsigned int c = 0;
			for( auto it=out.begin() ; it!=out.end() ; it++ ) idx[c++] = *it;
			for( auto it= in.begin() ; it!= in.end() ; it++ ) idx[c++] = *it;
			_l(idx) = 1;
		}
	}

	// A constructor generating a simple Linear function with one in the entry {_OutDims,_InDims} and zero for all other entries.
	template< unsigned int ... OutDims , unsigned int ... InDims >
	template< unsigned int ... _OutDims , unsigned int ... _InDims >
	Linear< UIntPack< OutDims ... > , UIntPack< InDims ... > >::Linear( UIntPack< _OutDims ... > , UIntPack< _InDims ... > )
	{
		static_assert( sizeof...(_OutDims)==sizeof...(OutDims) && sizeof...(_InDims)==sizeof...(InDims) , "[ERROR] Size mismatch" );
		unsigned int idx[] = { _OutDims ... , _InDims ... };
		_l(idx) = 1;
	}

	template< unsigned int ... OutDims , unsigned int ... InDims >
	auto Linear< UIntPack< OutDims ... > , UIntPack< InDims ... > >::value( const Tensor< UIntPack< InDims ... > > &t ) const
	{
		return _l.template contractedOuterProduct< sizeof ... ( InDims ) >( t );
	}

	template< unsigned int ... OutDims , unsigned int ... InDims >
	auto Linear< UIntPack< OutDims ... > , UIntPack< InDims ... > >::d( void ) const { return Constant< Concatenation< UIntPack< OutDims ... > , UIntPack< InDims ... > > , UIntPack< InDims ... > >(_l); }

	template< unsigned int ... OutDims , unsigned int ... InDims >
	std::ostream &operator << ( std::ostream &os , const Linear< UIntPack< OutDims ... > , UIntPack< InDims ... > > &l ){ return os << l._l; }

	//////////////
	// Identity //
	//////////////
	template< unsigned int ... Dims >
	Identity< UIntPack< Dims ... > >::Identity( void )
	{
		unsigned int index[ 2*Pack::Size ];
		WindowLoop< Pack::Size >::Run
		(
			ZeroUIntPack< Pack::Size >() , Pack() ,
			[&]( int d , int i ){ index[d] = index[d+Pack::Size] = i; } ,
			[&]( void ){ _l( index ) = 1; }
		);
	}

	//////////////////
	// _Coefficient //
	//////////////////
	template< unsigned int ... InDims >
	_Coefficient< UIntPack< InDims ... > >::_Coefficient( const unsigned int indices[] ){ _l( indices ) = 1; }
	template< unsigned int ... InDims >
	template< unsigned int ... CoefficientIndices >
	_Coefficient< UIntPack< InDims ... > >::_Coefficient( UIntPack< CoefficientIndices ... > )
	{
		typedef UIntPack< CoefficientIndices ... > CoefficientIndexPack;
		typedef UIntPack< InDims ... > InPack;
		static_assert( CoefficientIndexPack::Size==InPack::Size , "[ERROR] Coefficient count must match dimensions" );
		static_assert( Compare< CoefficientIndexPack , InPack >::LessThan , "[ERROR] Coefficient index cannot exceed dimension" );
		const unsigned int indices[] = { CoefficientIndices ... };
		_l( indices ) = 1;
	}

	/////////////////
	// Coefficient //
	/////////////////
	template< typename F >
	auto Coefficient( const unsigned int indices[] , const F &f ){ return _Coefficient< typename F::OutPack >( indices )( f ); }
	template< unsigned int ... CoefficientIndices , typename F >
	auto Coefficient( UIntPack< CoefficientIndices ... > , const F &f ){ return _Coefficient< typename F::OutPack >( UIntPack< CoefficientIndices ... >() )( f ); }

	/////////////////
	// Determinant //
	/////////////////
	template< unsigned int R , unsigned int C ,typename F >
	auto _SubMatrix( const F &f )
	{
		typedef typename F::OutPack OutPack;
		static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
		static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
		static const unsigned int Dim = OutPack::template Get<0>();

		// Create the (Dim-1)x(Dim-1) matrix
		Tensor< UIntPack< Dim-1 , Dim-1 , Dim , Dim > > l;
		unsigned int index[4];
		for( unsigned int i=0 ; i<Dim-1 ; i++ ) for( unsigned int j=0 ; j<Dim-1 ; j++ )
		{
			index[0] = i;
			index[1] = j;
			index[2] = i<R ? i : i+1;
			index[3] = j<C ? j : j+1;
			l( index ) = 1;
		}
		Linear< UIntPack< Dim-1 , Dim-1 > , UIntPack< Dim , Dim > > L( l );
		return L(f);
	}

	template< unsigned int C ,typename F >
	auto _Determinant0( const F &f )
	{
		typedef typename F::OutPack OutPack;
		static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
		static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
		static const unsigned int Dim = OutPack::template Get<0>();

		auto _det = Determinant( _SubMatrix<0,C>(f) ) * Coefficient( UIntPack<0,C>() , f );

		if constexpr( C==0 ) return _det;
		else
		{
			if constexpr( C&1 ) return _Determinant0< C-1 >( f ) - _det;
			else                return _Determinant0< C-1 >( f ) + _det;
		}
	}

	template< typename F >
	auto Determinant( const F &f )
	{
		typedef typename F::OutPack OutPack;
		static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
		static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
		static const unsigned int Dim = OutPack::template Get<0>();
		if      constexpr( Dim==1 ) return Coefficient( UIntPack<0,0>() , f );
		else if constexpr( Dim==2 ) return Coefficient( UIntPack<0,0>() , f ) * Coefficient( UIntPack<1,1>() , f ) - Coefficient( UIntPack<1,0>() , f ) * Coefficient( UIntPack<0,1>() , f );
		else return _Determinant0< Dim-1 >( f );
	}

	//////////////
	// Cofactor //
	//////////////

	template< unsigned int R , unsigned int C , typename F >
	auto _Cofactor( const F &f )
	{
		typedef typename F::InPack InPack;
		typedef typename F::OutPack OutPack;
		static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
		static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
		static const unsigned int Dim = OutPack::template Get<0>();

		// Create the Dim x Dim matrix
		Tensor< UIntPack< Dim , Dim > > e;
		{
			unsigned int index[2];
			index[0] = R;
			index[1] = C;
			e( index ) = 1;
		}

		auto _det = Constant< UIntPack< Dim , Dim > , InPack >( e ) * Determinant( _SubMatrix< R , C >( f ) );

		if constexpr( R==0 && C==0 ) return _det;
		else if constexpr( C==0 )
		{
			if constexpr( (R+C)&1 ) return _Cofactor< R-1 , Dim-1 >(f) - _det;
			else                    return _Cofactor< R-1 , Dim-1 >(f) + _det;
		}
		else
		{
			if constexpr( (R+C)&1 ) return _Cofactor< R , C-1 >( f ) - _det;
			else                    return _Cofactor< R , C-1 >( f ) + _det;
		}
	}

	template< typename F >
	auto Cofactor( const F& f )
	{
		typedef typename F::OutPack OutPack;
		static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
		static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
		static const unsigned int Dim = OutPack::template Get<0>();
		return _Cofactor< Dim-1 , Dim-1 >( f );
	}

	//////////////////
	// CrossProduct //
	//////////////////

	template< unsigned int R , typename F >
	auto _CrossProduct( const F &f )
	{
		typedef typename F::InPack InPack;
		typedef typename F::OutPack OutPack;
		static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
		static_assert( OutPack::template Get<0>()==OutPack::template Get<1>()+1 , "[ERROR] Output 2-tensor must have one more row than column" );
		static const unsigned int Dim = OutPack::template Get<0>();

		// Create the Dim-dimensional vector
		Tensor< UIntPack< Dim > > c;
		{
			unsigned int index[1];
			index[0] = R;
			c( index ) = 1;
		}

		// Create the (Dim-1)x(Dim-1) matrix
		Tensor< UIntPack< Dim-1 , Dim-1 , Dim , Dim-1 > > l;
		{
			unsigned int index[4];
			for( unsigned int i=0 ; i<Dim-1 ; i++ ) for( unsigned int j=0 ; j<Dim-1 ; j++ )
			{
				index[0] = i;
				index[1] = j;
				index[2] = i<R ? i : i+1;
				index[3] = j;
				l( index ) = 1;
			}
		}
		Constant< UIntPack< Dim > , InPack > C( c );
		Linear< UIntPack< Dim-1 , Dim-1 > , UIntPack< Dim , Dim-1 > > L( l );

		auto _det = C * Determinant( L( f ) );

		if constexpr( R==0 ) return _det;
		else
		{
			if constexpr( R&1 ) return _CrossProduct< R-1 >( f ) - _det;
			else                return _CrossProduct< R-1 >( f ) + _det;
		}
	}

	template< typename F >
	auto CrossProduct( const F& f )
	{
		typedef typename F::OutPack OutPack;
		static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
		static_assert( OutPack::template Get<0>()==OutPack::template Get<1>()+1 , "[ERROR] Output 2-tensor must have one more row than column" );
		static const unsigned int Dim = OutPack::template Get<0>();
		return _CrossProduct< Dim-1 >( f );
	}

	//////////////////
	// _Contraction //
	//////////////////
	template< unsigned int I1 , unsigned int I2 , unsigned int ... Dims >
	_Contraction< I1 , I2 , UIntPack< Dims ... > >::_Contraction( void ){ _l = Tensor< UIntPack< Dims ... > >::template ContractionTensor< I1 , I2 >(); }

	/////////////////
	// Contraction //
	/////////////////
	template< unsigned I1 , unsigned I2 , typename F >
	auto Contraction( const F &f )
	{
		if constexpr( I1<I2 ) return _Contraction< I1 , I2 , typename F::OutPack >()( f );
		else                  return _Contraction< I2 , I1 , typename F::OutPack >()( f );
	}

	///////////////
	// Transpose //
	///////////////
	template< typename F >
	auto Transpose( const F &f ){ return Permutation< ReverseUIntPack< SequentialUIntPack< F::OutPack::Size > > >( f ); }

	//////////////
	// Adjugate //
	//////////////
	template< typename F >
	auto Adjugate( const F &f ){ return Transpose( Cofactor( f ) ); }

	/////////////
	// Inverse //
	/////////////
	template< typename F >
	auto Inverse( const F &f ){ return Adjugate( f ) / Determinant( f ); }

	////////////////
	// SquareNorm //
	////////////////
	template< typename F >
	auto SquareNorm( const F &f ){ return ContractedOuterProduct< F::OutPack::Size >( f , f ); }

	////////////////////////
	// DiscreteDerivative //
	////////////////////////
	template< typename F >
	auto DiscreteDerivative( const F &f , const Tensor< typename F::InPack > &in , double eps )
	{
		typedef Tensor< typename F::InPack > InTensor;
		typedef Tensor< typename F::OutPack > OutTensor;

		if constexpr( F::InPack::Size==0 ) return ( F( in + InTensor( eps ) ) - F( in - InTensor( eps ) ) ) / ( 2. * eps );
		else
		{
			Tensor< Concatenation< typename F::OutPack , typename F::InPack > > out;

			unsigned int index[ F::InPack::Size ];
			WindowLoop< F::InPack::Size >::Run
			(
				ZeroUIntPack< F::InPack::Size >() , F::InPack() ,
				[&]( int d , int i ){ index[d] = i; } ,
				[&]( void )
				{
					InTensor in1 = in , in2 = in;
					in1( index ) += eps;
					in2( index ) -= eps;
					OutTensor _out = ( f( in1 ) - f( in2 ) ) / (2.*eps);
					if constexpr( F::OutPack::Size==0 ) out( index ) = (double)_out;
					else
					{
						unsigned int _index[ F::OutPack::Size + F::InPack::Size ];
						for( unsigned int i=0 ; i<F::InPack::Size ; i++ ) _index[i+F::OutPack::Size] = index[i];
						WindowLoop< F::OutPack::Size >::Run
						(
							ZeroUIntPack< F::OutPack::Size >() , F::OutPack() ,
							[&]( int d , int i ){ _index[d] = i; } ,
							[&]( const double &v ){ out( _index ) = v; } ,
							_out
						);
					}
				}
			);
			return out;
		}
	}

#include "AutoDiff.Transcendental.inc"
}
#endif // AUTO_DIFF_INCLUDED
