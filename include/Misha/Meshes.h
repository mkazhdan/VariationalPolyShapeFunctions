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

#ifndef MESHES_INCLUDED
#define MESHES_INCLUDED

#include "SimplexMesh.h"
#include "SimplexRefinableMesh.h"
#include "SolidSimplexRefinableMesh.h"

namespace Meshes
{
	// A function that takes an index and returns the spatial position of the point
	template< unsigned int EmbeddingDimension , typename VIndex > using     VertexPositionFunction = std::function< Point< double , EmbeddingDimension > ( VIndex ) >;
	template< unsigned int EmbeddingDimension                   > using FullVertexPositionFunction = std::function< Point< double , EmbeddingDimension > ( unsigned int ) >;

	// A function that takes a representation of a node and returns a boolean value indicating if it should exist within the coarsened representation
	template< unsigned int Dim , unsigned int Degree > using IsCoarseNodeFunction = std::function< bool ( typename SimplexMesh< Dim , Degree >::NodeMultiIndex ) >;
}

#include "PolygonMesh.h"
#include "PolyhedronMesh.h"

#endif // MESHES_INCLUDED