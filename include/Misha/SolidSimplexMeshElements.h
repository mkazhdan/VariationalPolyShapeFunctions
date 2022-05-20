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

#ifndef SOLID_SIMPLEX_MESH_ELEMENTS_INCLUDED
#define SOLID_SIMPLEX_MESH_ELEMENTS_INCLUDED
#include <iostream>
#include <map>
#include <vector>
#include <functional>
#include "Misha/Geometry.h"
#include "SolidSimplexBasis.h"
#include "SimplexMeshElements.h"

template< unsigned int BlockSize >
Eigen::SparseMatrix< double > BlockExpand( const Eigen::SparseMatrix< double > &A );

template< unsigned int Dim , unsigned int Degree >
struct SolidSimplexMesh : protected SimplexMesh< Dim , Degree >
{
	typedef MultiIndex< Degree , unsigned int , false > NodeIndex;
	static const unsigned int NodesPerSimplex = SimplexElements< Dim , Degree >::NodeNum;

	SolidSimplexMesh( void ){}

	template< typename Index >
	static SolidSimplexMesh Init( const std::vector< SimplexIndex< Dim , Index > > &simplices , std::function< Point< double , Dim > ( Index ) > vFunction );

	size_t simplices( void ) const { return _simplices.size(); }
	size_t nodes( void ) const { return _nodeMap.size(); }

	Eigen::SparseMatrix< double > massMatrix( void ) const;
	Eigen::SparseMatrix< double > frobeniusStiffnessMatrix( void ) const;
	Eigen::SparseMatrix< double > traceStiffnessMatrix( void ) const;
	Eigen::VectorXd dcVector( Point< double , Dim > v ) const;
	Eigen::VectorXd stiffnessVector( void ) const;

	Eigen::SparseMatrix< double > evaluationMatrix( const std::vector< typename SimplexMesh< Dim , Degree >::Sample > &samples ) const;

	SolidSimplexMesh( SolidSimplexMesh &&sm ){ std::swap( _simplices , sm._simplices ) , std::swap( _A , sm._A ) , std::swap( _nodeMap , sm._nodeMap ); }
	SolidSimplexMesh & operator = ( SolidSimplexMesh &&sm ){ std::swap( _simplices , sm._simplices ) , std::swap( _A , sm._A ) , std::swap( _nodeMap , sm._nodeMap ) ; return *this; }

	NodeIndex nodeMultiIndex( unsigned int s , unsigned int n ) const;
	unsigned int nodeIndex( const NodeIndex &multiIndex ) const;
	unsigned int nodeIndex( unsigned int s , unsigned int n ) const { return nodeIndex( nodeMultiIndex( s , n ) ); }

	static Eigen::SparseMatrix< double > ScalarToSolid( const Eigen::SparseMatrix< double > &A );

protected:
	std::vector< SquareMatrix< double , Dim > > _A;
};

#include "SolidSimplexMesh.inl"
#endif // SOLID_SIMPLEX_MESH_ELEMENTS_INCLUDED