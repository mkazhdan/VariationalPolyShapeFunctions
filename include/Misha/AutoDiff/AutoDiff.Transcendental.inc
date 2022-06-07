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


struct _Pow : public Function< UIntPack::Pack<> , UIntPack::Pack<> , _Pow >
{
	_Pow( double e ) : _e(e) {}
	Tensor< UIntPack::Pack<> > value( const Tensor< UIntPack::Pack<> > &t ) const { return Tensor< UIntPack::Pack<> >( pow( static_cast< double >( t ) , _e ) ); }
	auto d( void ) const { return _Pow( _e-1. ) * _e; }
	friend std::ostream &operator<< ( std::ostream &os , const _Pow &p ){ return os << "pow( " << p._e << " )"; }
protected:
	double _e;
};

template< typename F > auto Pow( const F &f , double e )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Pow , F >( _Pow(e) , f );
}

template< typename F > auto Sqrt( const F &f ){ return Pow( f , 0.5 ); }

template< typename F >
auto operator / ( double n , const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return n * Pow(f,-1.);
}


struct _Exp : public Function< UIntPack::Pack<> , UIntPack::Pack<> , _Exp >
{
	Tensor< UIntPack::Pack<> > value( const Tensor< UIntPack::Pack<> > &t ) const { return Tensor< UIntPack::Pack<> >( exp( static_cast< double >( t ) ) ); }
	auto d( void ) const { return _Exp(); }
	friend std::ostream &operator<< ( std::ostream &os , const _Exp & ){ return os << "exp"; }
};

template< typename F > auto Exp( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Exp , F >( _Exp() , f );
}

struct _Log : public Function< UIntPack::Pack<> , UIntPack::Pack<> , _Log >
{
	Tensor< UIntPack::Pack<> > value( const Tensor< UIntPack::Pack<> > &t ) const { return Tensor< UIntPack::Pack<> >( log( static_cast< double >( t ) ) ); }
	auto d( void ) const { return _Pow(-1.); }
	friend std::ostream &operator<< ( std::ostream &os , const _Log & ){ return os << "log"; }
};

template< typename F > auto Log( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Log , F >( _Log() , f );
}


struct _Sin : public Function< UIntPack::Pack<> , UIntPack::Pack<> , _Sin >
{
	Tensor< UIntPack::Pack<> > value( const Tensor< UIntPack::Pack<> > &t ) const { return Tensor< UIntPack::Pack<> >( sin( static_cast< double >( t ) ) ); }
	auto d( void ) const;
	friend std::ostream &operator<< ( std::ostream &os , const _Sin & ){ return os << "sin"; }
};

template< typename F > auto Sin( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Sin , F >( _Sin() , f );
}


struct _Cos : public Function< UIntPack::Pack<> , UIntPack::Pack<> , _Cos >
{
	Tensor< UIntPack::Pack<> > value( const Tensor< UIntPack::Pack<> > &t ) const { return Tensor< UIntPack::Pack<> >( cos( static_cast< double >( t ) ) ); }
	auto d( void ) const { return -_Sin(); }
	friend std::ostream &operator << ( std::ostream &os , const _Cos & ){ return os << "cos"; }
};

template< typename F > auto Cos( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Cos , F >( _Cos() , f );
}

auto _Sin::d( void ) const { return _Cos(); }


struct _Sinh : public Function< UIntPack::Pack<> , UIntPack::Pack<> , _Sinh >
{
	Tensor< UIntPack::Pack<> > value( const Tensor< UIntPack::Pack<> > &t ) const { return Tensor< UIntPack::Pack<> >( sinh( static_cast< double >( t ) ) ); }
	auto d( void ) const;
	friend std::ostream &operator<< ( std::ostream &os , const _Sinh & ){ return os << "sinh"; }
};

template< typename F > auto Sinh( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Sinh , F >( _Sinh() , f );
}


struct _Cosh : public Function< UIntPack::Pack<> , UIntPack::Pack<> , _Cosh >
{
	Tensor< UIntPack::Pack<> > value( const Tensor< UIntPack::Pack<> > &t ) const { return Tensor< UIntPack::Pack<> >( cosh( static_cast< double >( t ) ) ); };
	auto d( void ) const { return _Sinh(); }
	friend std::ostream &operator << ( std::ostream &os , const _Cosh & ){ return os << "cosh"; }
};
template< typename F > auto Cosh( const F &f )
{
	static_assert( F::OutPack::Size==0 , "[ERROR] Function should be scalar-valued" );
	return Composition< _Cosh , F >( _Cosh() , f );
}

auto _Sinh::d( void ) const { return _Cosh(); }
