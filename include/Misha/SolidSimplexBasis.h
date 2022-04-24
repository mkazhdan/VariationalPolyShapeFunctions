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

#ifndef SOLID_SIMPLEX_BASIS_INCLUDED
#define SOLID_SIMPLEX_BASIS_INCLUDED
#include "SimplexBasis.h"

template< unsigned int Dim , unsigned int Degree >
struct SolidSimplexElements
{
	enum
	{
		SYSTEM_MASS ,
		SYSTEM_FROBENIUS_STIFFNESS ,
		SYSTEM_TRACE_STIFFNESS ,
		SYSTEM_COUNT
	};
	static const std::string EnergyNames[];

	static SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > SystemMatrix( SquareMatrix< double , Dim > A , unsigned int type );

	static SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim >               MassMatrix( SquareMatrix< double , Dim > A );
	static SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > FrobeniusStiffnessMatrix( SquareMatrix< double , Dim > A );
	static SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim >     TraceStiffnessMatrix( SquareMatrix< double , Dim > A );
	static void SetMassFrobeniusStiffnessAndTraceStiffnessMatrices( SquareMatrix< double , Dim > A , SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > &M , SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > &F , SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum * Dim > & T );
	static Matrix< double , Dim , SimplexElements< Dim , Degree >::NodeNum * Dim >                       DC( SquareMatrix< double , Dim > A );
	static Point< double , SimplexElements< Dim , Degree >::NodeNum * Dim >                 StiffnessVector( SquareMatrix< double , Dim > A );

	struct Integrator
	{
		Integrator( SquareMatrix< double , Dim > A );
		double dc( unsigned int n1 , unsigned int d1 , unsigned int d2 ) const;
		double mass( unsigned int n1 , unsigned int d1 , unsigned int n2 , unsigned int d2 ) const;
		double frobeniusStiffness( unsigned int n1 , unsigned int d1 , unsigned int n2 , unsigned int d2 ) const;
		double traceStiffness( unsigned int n1 , unsigned int d1 , unsigned int n2 , unsigned int d2 ) const;
		double stiffness( unsigned int n , unsigned int d ) const;
	protected:
		double _det;
		SquareMatrix< double , Dim  > _A[Dim][Dim];
		Point< double , Dim > _e[Dim];
		double _stiffness( unsigned int n , unsigned int d ) const;
		double _stiffness( unsigned int n1 , unsigned int d1 , unsigned int n2 , unsigned int d2 ) const;

		static bool _Initialized;

		template< typename Real >
		static void _Init( void );
		static Point< double , Dim > _B[ SimplexElements< Dim , Degree >::NodeNum ];
		static Point< double , SimplexElements< Dim , Degree >::NodeNum > _DC;
		static SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum  > _Mass;
		static SquareMatrix< double , Dim > _Stiffness[ SimplexElements< Dim , Degree >::NodeNum ][ SimplexElements< Dim , Degree >::NodeNum ];
	};
};

#include "SolidSimplexBasis.inl"
#endif // SOLID_SIMPLEX_BASIS_INCLUDED