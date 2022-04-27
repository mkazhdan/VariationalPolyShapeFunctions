/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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

/////////////////////
// I != N: Recurse //
/////////////////////
template< unsigned int I , unsigned int N >
struct _WindowLoop
{
	static_assert( I<N , "[ERROR] Index exceeds bound" );
protected:
	friend struct WindowLoop< N >;
	template< unsigned int _I , unsigned int _N > friend struct _WindowLoop;

	template< typename UpdateFunction , typename ProcessFunction , class ... Windows >
	static void Run( int begin , int end , UpdateFunction& updateState , ProcessFunction& function , Windows&& ... w )
	{
		for( int i=begin ; i<end ; i++ ){ updateState( I , i ) ; _WindowLoop< I+1 , N >::Run( begin , end , updateState , function , std::forward< decltype(w[i]) >( w[i] ) ... ); }
	}
	template< typename UpdateFunction , typename ProcessFunction , class ... Windows >
	static void Run( const int* begin , const int* end , UpdateFunction& updateState , ProcessFunction& function , Windows&& ... w )
	{
		for( int i=begin[0] ; i<end[0] ; i++ ){ updateState( I , i ) ; _WindowLoop< I+1 , N >::Run( begin+1 , end+1 , updateState , function , std::forward< decltype(w[i]) >( w[i] ) ... ); }
	}
	template< unsigned int ... Begin , unsigned int ... End , typename UpdateFunction , typename ProcessFunction , class ... Windows >
	static void Run( UIntPack< Begin ... > begin , UIntPack< End ... > end , UpdateFunction& updateState , ProcessFunction& function , Windows&& ... w )
	{
		for( int i=UIntPack< Begin ... >::First ; i<UIntPack< End ... >::First ; i++ ){ updateState( I , i ) ; _WindowLoop< I+1 , N >::Run( typename UIntPack< Begin ... >::Rest() , typename UIntPack< End ... >::Rest() , updateState , function , std::forward< decltype(w[i]) >( w[i] ) ... ); }
	}
};

///////////////////////
// I == N: Base case //
///////////////////////
template< unsigned int N >
struct _WindowLoop< N , N >
{
protected:
	static const unsigned int I = N;
	friend struct WindowLoop< N >;
	template< unsigned int _I , unsigned int _N > friend struct _WindowLoop;

	template< typename UpdateFunction , typename ProcessFunction , class ... Values >
	static void Run( int , int , UpdateFunction & , ProcessFunction &function , Values&& ... v )
	{
		function( std::forward< Values >( v ) ... );
	}
	template< typename UpdateFunction , typename ProcessFunction , class ... Values >
	static void Run( const int * , const int * , UpdateFunction & , ProcessFunction &function , Values&& ... v )
	{
		function( std::forward< Values >( v ) ... );
	}
	template< unsigned int ... Begin , unsigned int ... End , typename UpdateFunction , typename ProcessFunction , class ... Values >
	static void Run( UIntPack< Begin ... > , UIntPack< End ... > , UpdateFunction & , ProcessFunction &function , Values&& ... v )
	{
		function( std::forward< Values >( v ) ... );
	}
};

