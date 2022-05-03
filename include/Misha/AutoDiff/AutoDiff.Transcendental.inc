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


struct _Pow : public Function< UIntPack<> , UIntPack<> , _Pow >
{
	_Pow( double e ) : _e(e) {}
	template< unsigned int D >
	Tensor<> d( const Tensor<> &t ) const { return Tensor<>( pow( static_cast< double >( t ) , _e-D ) * _Scalar<D>( _e ) ); }
	auto d( void ) const { return _Pow( _e-1. ) * _e; }
	friend std::ostream &operator<< ( std::ostream &os , const _Pow &p ){ return os << "pow( " << p._e << " )"; }
protected:
	double _e;

	template< unsigned int D >
	static double _Scalar( double e )
	{
		if constexpr( D==0 ) return 1;
		else return e * _Scalar< D-1 >( e-1 );
	}
};

template< typename F > auto Pow( const F &f , double e )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Pow , F >( _Pow(e) , f );
}

template< typename F >
auto operator / ( double n , const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return n * Pow(f,-1.);
}


struct _Exp : public Function< UIntPack<> , UIntPack<> , _Exp >
{
	template< unsigned int D >
	Tensor<> d( const Tensor<> &t ) const { return Tensor<>( exp( static_cast< double >( t ) ) ); }
	_Exp d( void ) const { return _Exp(); }
	friend std::ostream &operator<< ( std::ostream &os , const _Exp & ){ return os << "exp"; }
};

template< typename F > auto Exp( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Exp , F >( _Exp() , f );
}

struct _Log : public Function< UIntPack<> , UIntPack<> , _Log >
{
	template< unsigned int D >
	Tensor<> d( const Tensor<> &t ) const
	{
		if constexpr( D==0 ) return Tensor<>( log( static_cast< double >( t ) ) );
		else return _Pow( -1. ).template d< D-1  >( t );
	}
	_Pow d( void ) const { return _Pow(-1.); }
	friend std::ostream &operator<< ( std::ostream &os , const _Log & ){ return os << "log"; }
};

template< typename F > auto Log( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Log , F >( _Log() , f );
}


struct _Sin : public Function< UIntPack<> , UIntPack<> , _Sin >
{
	template< unsigned int D >
	Tensor<> d( const Tensor<> &t ) const
	{
		if( D&1 )
			if( (D>>1)&1 ) return Tensor<>( -cos( static_cast< double >( t ) ) );
			else           return Tensor<>(  cos( static_cast< double >( t ) ) );
		else
			if( (D>>1)&1 ) return Tensor<>( -sin( static_cast< double >( t ) ) );
			else           return Tensor<>(  sin( static_cast< double >( t ) ) );
	}
	struct _Cos d( void ) const;
	friend std::ostream &operator<< ( std::ostream &os , const _Sin & ){ return os << "sin"; }
};

template< typename F > auto Sin( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Sin , F >( _Sin() , f );
}


struct _Cos : public Function< UIntPack<> , UIntPack<> , _Cos >
{
	template< unsigned int D >
	Tensor<> d( const Tensor<> &t ) const
	{
		if( D&1 )
			if( (D>>1)&1 ) return Tensor<>(  sin( static_cast< double >( t ) ) );
			else           return Tensor<>( -sin( static_cast< double >( t ) ) );
		else
			if( (D>>1)&1 ) return Tensor<>( -cos( static_cast< double >( t ) ) );
			else           return Tensor<>(  cos( static_cast< double >( t ) ) );
	}
	auto d( void ) const { return -_Sin(); }
	friend std::ostream &operator << ( std::ostream &os , const _Cos & ){ return os << "cos"; }
};

template< typename F > auto Cos( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Cos , F >( _Cos() , f );
}

_Cos _Sin::d( void ) const { return _Cos(); }


struct _Sinh : public Function< UIntPack<> , UIntPack<> , _Sinh >
{
	template< unsigned int D >
	Tensor<> d( const Tensor<> &t ) const
	{
		if( D&1 ) return Tensor<>( cosh( static_cast< double >( t ) ) );
		else      return Tensor<>( sinh( static_cast< double >( t ) ) );
	}
	struct _Cosh d( void ) const;
	friend std::ostream &operator<< ( std::ostream &os , const _Sinh & ){ return os << "sinh"; }
};

template< typename F > auto Sinh( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Sinh , F >( _Sinh() , f );
}


struct _Cosh : public Function< UIntPack<> , UIntPack<> , _Cosh >
{
	template< unsigned int D >
	Tensor<> d( const Tensor<> &t ) const
	{
		if( D&1 ) return Tensor<>( sinh( static_cast< double >( t ) ) );
		else      return Tensor<>( cosh( static_cast< double >( t ) ) );
	}
	auto d( void ) const { return _Sinh(); }
	friend std::ostream &operator << ( std::ostream &os , const _Cosh & ){ return os << "cosh"; }
};
template< typename F > auto Cosh( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Cosh , F >( _Cosh() , f );
}

_Cosh _Sinh::d( void ) const { return _Cosh(); }