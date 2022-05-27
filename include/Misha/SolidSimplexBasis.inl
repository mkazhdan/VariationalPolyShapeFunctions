/*
Copyright (c) 2022, Michael Kazhdan
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

////////////////////////////////////////////////
// SolidSimplexElements::_StiffnessIntegrator //
////////////////////////////////////////////////
template< unsigned int Dim , unsigned int Degree >
bool SolidSimplexElements< Dim , Degree >::Integrator::_Initialized = false;

template< unsigned int Dim , unsigned int Degree >
Point< double , SimplexElements< Dim , Degree >::NodeNum  > SolidSimplexElements< Dim , Degree >::Integrator::_DC;

template< unsigned int Dim , unsigned int Degree >
SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum  > SolidSimplexElements< Dim , Degree >::Integrator::_Mass;

template< unsigned int Dim , unsigned int Degree >
SquareMatrix< double , Dim > SolidSimplexElements< Dim , Degree >::Integrator::_Stiffness[ SimplexElements< Dim , Degree >::NodeNum ][ SimplexElements< Dim , Degree >::NodeNum ];

template< unsigned int Dim , unsigned int Degree >
Point< double , Dim > SolidSimplexElements< Dim , Degree >::Integrator::_B[ SimplexElements< Dim , Degree >::NodeNum ];

template< unsigned int Dim , unsigned int Degree >
template< typename Real >
void SolidSimplexElements< Dim , Degree >::Integrator::_Init( void )
{
	if( !_Initialized )
	{
		Polynomial::Polynomial< Dim , Degree , Real > elements[ SimplexElements< Dim , Degree >::NodeNum ];
		Polynomial::Polynomial< Dim , Degree-1 , Real > dElements[ SimplexElements< Dim , Degree >::NodeNum ][ Dim ];
		SimplexElements< Dim , Degree >::SetElements( elements );
		for( unsigned int n=0 ; n<SimplexElements< Dim , Degree >::NodeNum ; n++ ) for( int d=0 ; d<Dim ; d++ ) dElements[n][d] = elements[n].d( d );

		for( unsigned int n1=0 ; n1<SimplexElements< Dim , Degree >::NodeNum ; n1++ )
		{
			_DC[n1] = RightSimplex< Dim >::Integral( elements[n1] );
			for( unsigned int n2=0 ; n2<SimplexElements< Dim , Degree >::NodeNum ; n2++ )
				_Mass(n1,n2) = RightSimplex< Dim >::Integral( elements[n1] * elements[n2] );
		}

		for( unsigned int n1=0 ; n1<SimplexElements< Dim , Degree >::NodeNum ; n1++ ) for( unsigned int d1=0 ; d1<Dim ; d1++ )
		{
			_B[n1][d1] = RightSimplex< Dim >::Integral( dElements[n1][d1] );
			for( unsigned int n2=0 ; n2<SimplexElements< Dim , Degree >::NodeNum ; n2++ ) for( unsigned int d2=0 ; d2<Dim ; d2++ )
				_Stiffness[n1][n2](d1,d2) = RightSimplex< Dim >::Integral( dElements[n1][d1] * dElements[n2][d2] );
		}

		_Initialized = true;
	}
}

template< unsigned int Dim , unsigned int Degree >
SolidSimplexElements< Dim , Degree >::Integrator::Integrator( SquareMatrix< double , Dim > A )
{
	_Init< double >();

	_det = sqrt( fabs( A.determinant() ) );

	SquareMatrix< double , Dim > A_inv = A.inverse();

	for( unsigned int d1=0 ; d1<Dim ; d1++ ) for( unsigned int d2=0 ; d2<Dim ; d2++ ) 
		for( unsigned int _d1=0 ; _d1<Dim ; _d1++ ) for( unsigned int _d2=0 ; _d2<Dim ; _d2++ )
			_A[d1][d2](_d1,_d2) = A_inv( d1 , _d1 ) * A_inv( d2 , _d2 );

	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		Point< double , Dim > e;
		e[d] = 1.;
		_e[d] = A_inv * e;
	}
}

template< unsigned int Dim , unsigned int Degree >
double SolidSimplexElements< Dim , Degree >::Integrator::dc( unsigned int n1 , unsigned int d1 , unsigned int d2 ) const
{
	return d1==d2 ? _DC[n1] * _det : 0;
}

template< unsigned int Dim , unsigned int Degree >
double SolidSimplexElements< Dim , Degree >::Integrator::mass( unsigned int n1 , unsigned int d1 , unsigned int n2 , unsigned int d2 ) const
{
	return d1==d2 ? _Mass(n1,n2) * _det : 0;
}

template< unsigned int Dim , unsigned int Degree >
double SolidSimplexElements< Dim , Degree >::Integrator::_stiffness( unsigned int n , unsigned int d ) const
{
	return Point< double , Dim >::Dot( _e[d] , _B[n] ) * _det;
}

template< unsigned int Dim , unsigned int Degree >
double SolidSimplexElements< Dim , Degree >::Integrator::_stiffness( unsigned int n1 , unsigned int d1 , unsigned int n2 , unsigned int d2 ) const
{
	return SquareMatrix< double , Dim >::Dot( _A[d1][d2] , _Stiffness[n1][n2] ) * _det;
}

template< unsigned int Dim , unsigned int Degree >
double SolidSimplexElements< Dim , Degree >::Integrator::frobeniusStiffness( unsigned int n1 , unsigned int d1 , unsigned int n2 , unsigned int d2 ) const
{
	double dot = 0;
	if( d1==d2 ) for( unsigned int d=0 ; d<Dim ; d++ ) dot += _stiffness( n1 , d , n2 , d );
	dot += _stiffness( n1 , d2 , n2 , d1 );
	return dot/2;
}

template< unsigned int Dim , unsigned int Degree >
double SolidSimplexElements< Dim , Degree >::Integrator::traceStiffness( unsigned int n1 , unsigned int d1 , unsigned int n2 , unsigned int d2 ) const
{
	return _stiffness( n1 , d1 , n2 , d2 );
}

template< unsigned int Dim , unsigned int Degree >
double SolidSimplexElements< Dim , Degree >::Integrator::stiffness( unsigned int n , unsigned int d ) const
{
	return _stiffness( n , d );
}

//////////////////////////
// SolidSimplexElements //
//////////////////////////
template< unsigned int Dim , unsigned int Degree >
const std::string SolidSimplexElements< Dim , Degree >::EnergyNames[] = { "mass" , "frobenius-stiffness" , "trace-stiffness" };

template< unsigned int Dim , unsigned int Degree >
Matrix< double , Dim , SimplexElements< Dim , Degree >::NodeNum * Dim > SolidSimplexElements< Dim , Degree >::DC( SquareMatrix< double , Dim > A )
{
	Integrator integrator( A );
	Matrix< double , Dim , SimplexElements< Dim , Degree >::NodeNum*Dim > m;
	for( unsigned int n1=0 ; n1<SimplexElements< Dim , Degree >::NodeNum ; n1++ ) for( unsigned int d1=0 ; d1<Dim ; d1++ ) for( unsigned int d2=0 ; d2<Dim ; d2++ )
		m( d2 , n1*Dim+d1 ) = integrator.dc(n1,d1,d2);
	return m;
}

template< unsigned int Dim , unsigned int Degree >
SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > SolidSimplexElements< Dim , Degree >::MassMatrix( SquareMatrix< double , Dim > A )
{
	Integrator integrator( A );
	SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum*Dim > m;
	for( unsigned int n1=0 ; n1<SimplexElements< Dim , Degree >::NodeNum ; n1++ ) for( int unsigned n2=0 ; n2<SimplexElements< Dim , Degree >::NodeNum ; n2++ )
		for( unsigned int d1=0 ; d1<Dim ; d1++ ) for( unsigned int d2=0 ; d2<Dim ; d2++ ) m( n1*Dim+d1 , n2*Dim+d2 ) = integrator.mass(n1,d1,n2,d2);
	return m;
}

template< unsigned int Dim , unsigned int Degree >
SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > SolidSimplexElements< Dim , Degree >::FrobeniusStiffnessMatrix( SquareMatrix< double , Dim > A )
{
	Integrator integrator( A );

	SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum*Dim > s;
	for( unsigned int n1=0 ; n1<SimplexElements< Dim , Degree >::NodeNum ; n1++ ) for( int unsigned n2=0 ; n2<SimplexElements< Dim , Degree >::NodeNum ; n2++ )
		for( unsigned int d1=0 ; d1<Dim ; d1++ ) for( unsigned int d2=0 ; d2<Dim ; d2++ ) s( n1*Dim+d1 , n2*Dim+d2 ) = integrator.frobeniusStiffness(n1,d1,n2,d2);
	return s;
}

template< unsigned int Dim , unsigned int Degree >
SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > SolidSimplexElements< Dim , Degree >::TraceStiffnessMatrix( SquareMatrix< double , Dim > A )
{
	Integrator integrator( A );

	SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum*Dim > s;
	for( unsigned int n1=0 ; n1<SimplexElements< Dim , Degree >::NodeNum ; n1++ ) for( int unsigned n2=0 ; n2<SimplexElements< Dim , Degree >::NodeNum ; n2++ )
		for( unsigned int d1=0 ; d1<Dim ; d1++ ) for( unsigned int d2=0 ; d2<Dim ; d2++ ) s( n1*Dim+d1 , n2*Dim+d2 ) = integrator.traceStiffness(n1,d1,n2,d2);
	return s;
}

template< unsigned int Dim , unsigned int Degree >
void SolidSimplexElements< Dim , Degree >::SetMassFrobeniusStiffnessAndTraceStiffnessMatrices( SquareMatrix< double , Dim > A , SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > &M , SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > &F , SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > &T )
{
	Integrator integrator( A );

	for( unsigned int n1=0 ; n1<SimplexElements< Dim , Degree >::NodeNum ; n1++ ) for( int unsigned n2=0 ; n2<SimplexElements< Dim , Degree >::NodeNum ; n2++ )
		for( unsigned int d1=0 ; d1<Dim ; d1++ ) for( unsigned int d2=0 ; d2<Dim ; d2++ )
		{
			M( n1*Dim+d1 , n2*Dim+d2 ) = integrator.mass(n1,d1,n2,d2);
			F( n1*Dim+d1 , n2*Dim+d2 ) = integrator.frobeniusStiffness(n1,d1,n2,d2);
			T( n1*Dim+d1 , n2*Dim+d2 ) = integrator.traceStiffness(n1,d1,n2,d2);
		}
}

template< unsigned int Dim , unsigned int Degree >
Point< double , SimplexElements< Dim , Degree >::NodeNum * Dim > SolidSimplexElements< Dim , Degree >::StiffnessVector( SquareMatrix< double , Dim > A )
{
	Integrator integrator( A );
	Point< double , SimplexElements< Dim , Degree >::NodeNum*Dim > s;
	for( unsigned int n=0 ; n<SimplexElements< Dim , Degree >::NodeNum ; n++ ) for( unsigned int d=0 ; d<Dim ; d++ )
		s[ n*Dim+d ] = integrator.stiffness(n,d);
	return s;
}


template< unsigned int Dim , unsigned int Degree >
SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > SolidSimplexElements< Dim , Degree >::SystemMatrix( SquareMatrix< double , Dim > A , unsigned int type )
{
	switch( type )
	{
		case SYSTEM_MASS:                return               MassMatrix( A );
		case SYSTEM_FROBENIUS_STIFFNESS: return FrobeniusStiffnessMatrix( A );
		case SYSTEM_TRACE_STIFFNESS:     return     TraceStiffnessMatrix( A );
		default: ERROR_OUT( "Unrecognized energy type: " , type );
	}
}
