<center><h2>Variational Polygonal/Polyhedral Shape Functions (Version 2.5)</h2></center>
<center>
<a href="#LINKS">links</a>
<a href="#EXECUTABLES">executables</a>
<a href="#USAGE">usage</a>
<a href="#COMPILATION">compilation</a>
<a href="#CHANGES">changes</a>
<!--
<a href="#SUPPORT">support</a>
-->
</center>
<hr>
This software supports finite-elements-type calculations over polygonal and polyhedral meshes. Supported applications include:
<UL>
<LI>Simulation of deformable solids (with linear elasticity) in 2D and 3D,</LI>
<LI>Solution of the Franke test in 2D and 3D,</LI>
<LI>Calculation of Geodesics in Heat on polygonal meshes, and</LI>
<LI>Gradient domain processing of signals on polygonal meshes.</LI>
</UL>
<hr>
<a name="LINKS"><b>LINKS</b></a><br>
<ul>
<b>Papers:</b>
<a href="http://www.cs.jhu.edu/~misha/MyPapers/SIG22.pdf">[Bunge, Herholz, Sorkine-Hornung, Botsch, and Kazhdan, 2022]</a>,
<a href="https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/">[Crane, Weischedel, and Wardetzky, 2013]</a>
<br>
<b>Executables: </b>
<a href="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/VPSF.x64.zip">Win64</a><br>
<b>Source Code:</b>
<a href="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/VPSF.Source.zip">ZIP</a> <a href="https://github.com/mkazhdan/VariationalPolyShapeFunctions">GitHub</a><br>
<B>Data:</B>
<A HREF="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/VPSF.Data.zip">ZIP</A><br>
<b>Older Versions:</b>
<a href="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/Version2.00/">V2</a>
<a href="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/Version1.00/">V1</a>
</ul>
<hr>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>

<ul>
<dl>
<details>
<summary>
<font size="+1"><b>DeformableSolids2D/DeformableSolids3D</b></font>:
Supports the simulation of deformable solids in 2D and 3D using linear elasticity. The executable launches an interactive viewer that provides a visualization of the deforming solid. The solid can deform either through the action of gravity or by applying a prescribed linear transformation and having the solid evolve towards its rest state. (In 2D, the applications supports selecting and dragging of individual vertices.)<BR>
Hit [SPACE] to start the simulation or "+" to advance one time-step.
</summary>
<dt><b>--in</b> &lt;<i>input polygonal/polyhedral mesh</i>&gt;</dt>
<dd> This string specifies the the name of the mesh.<br>
For 2D simulations, the input polygonal mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i> and the set of polygons encoded by a list of vertex indices.<br>
For 3D simulations, the input polyhedral mesh is assumed to be in <a HREF="https://www.graphics.rwth-aachen.de/software/openvolumemesh/">OVM</a> format.
</dd>

<dt>[<b>--xForm</b> &lt;<i>linear transform</i>&gt;]</dt>
<dd> This 2x2 (resp. 3x3) set of floating point values describes the entries of the linear transformation initially applied to the solid.<BR>
The default values spcify the identity transformation.
</dd>

</dd><dt>[<b>--lock</B>]</dt>
<dd> If enabled, this flag specifies that the values on the <i>y</i>-axis (resp. <i>yz</i>-plane) should be locked during the course of the animation.
</dd>

<dt>[<b>--gravity</b> &lt;<i>gravitational force</i>&gt;]</dt>
<dd> This floating point value describes the force of gravity acting on the solid. (Note that without the <b>--lock</b> parameter, using a non-zero value for gravity will have the solid fall off the screen.)<BR>
The default value for this parameter is -500,000,000.
</dd>

</dd><dt>[<b>--mg</B>]</dt>
<dd> If enabled, this flag specifies that a multigrid solver should be used (instead of the default sparse Cholesky solver).
</dd>

<dt>[<b>--vCycles</b> &lt;<i>number of v-cycles per animation step</i>&gt;]</dt>
<dd> If a multigrid solver is used, ths integer value specifies the number of v-cycles to be performed at each step of the animation.<BR>
The default value for this parameter is 1.
</dd>

<dt>[<b>--gsIters</b> &lt;<i>number of Gauss-Seidel iterations per level</i>&gt;]</dt>
<dd> If a multigrid solver is used, ths integer value specifies the number of Gauss-Seidel iterations to be done within each level of the v-cycle.<BR>
The default value for this parameter is 5.
</dd>

</details>
</dl>
</ul>



<ul>
<dl>
<details>
<summary>
<font size="+1"><b>FrankeTest2D/FrankeTest3D</b></font>:
Supports the evaluation of function space quality by solving a Poisson equation over the unit square/cube, with boundary values fixed to the analytic values of the Franke test function. The executable takes in geometry, and outputs the RMS of the solution (compared to the analytic solution).
</summary>
<dt><b>--in</b> &lt;<i>input polygonal/polyhedral mesh</i>&gt;</dt>
<dd> This string specifies the the name of the mesh.<br>
For 2D simulations, the input polygonal mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i> and the set of polygons encoded by a list of vertex indices.<br>
For 3D simulations, the input polyhedral mesh is assumed to be in <a HREF="https://www.graphics.rwth-aachen.de/software/openvolumemesh/">OVM</a> format.
</dd>

</dd><dt>[<b>--mg</B>]</dt>
<dd> If enabled, this flag specifies that a multigrid solver should be used (instead of the default sparse Cholesky solver).
</dd>

<dt>[<b>--vCycles</b> &lt;<i>number of v-cycles per animation step</i>&gt;]</dt>
<dd> If a multigrid solver is used, ths integer value specifies the number of v-cycles to be performed at each step of the animation.<BR>
The default value for this parameter is 3.
</dd>

<dt>[<b>--gsIters</b> &lt;<i>number of Gauss-Seidel iterations per level</i>&gt;]</dt>
<dd> If a multigrid solver is used, ths integer value specifies the number of Gauss-Seidel iterations to be done within each level of the v-cycle.<BR>
The default value for this parameter is 5.
</dd>

</details>
</dl>
</ul>



<ul>
<dl>
<details>
<summary>
<font size="+1"><b>GeodesicsInHeat</b></font>:
Supports the interactive visualization of single-source geodesics on the surface of polygonal mesh using the <A HREF="https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/">heat method</A>.<BR>
Left-clicking while holding down the [SHIFT] key selects the source.
</summary>
<dt><b>--in</b> &lt;<i>input polygonal mesh</i>&gt;</dt>
<dd> This string specifies the the name of the mesh.<br>
The input polygonal mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i> and the set of polygons encoded by a list of vertex indices.
</dd>

<dt>[<b>--time</b> &lt;<i>diffusion time</i>&gt;]</dt>
<dd> This floating point values specifies the time for diffusing the source delta function .<BR>
The default value for this parameter is 1e-3.
</dd>

</dd><dt>[<b>--mg</B>]</dt>
<dd> If enabled, this flag specifies that a multigrid solver should be used (instead of the default sparse Cholesky solver).
</dd>

<dt>[<b>--vCycles</b> &lt;<i>number of v-cycles per animation step</i>&gt;]</dt>
<dd> If a multigrid solver is used, ths integer value specifies the number of v-cycles to be performed at each step of the animation.<BR>
The default value for this parameter is 1.
</dd>

<dt>[<b>--gsIters</b> &lt;<i>number of Gauss-Seidel iterations per level</i>&gt;]</dt>
<dd> If a multigrid solver is used, ths integer value specifies the number of Gauss-Seidel iterations to be done within each level of the v-cycle.<BR>
The default value for this parameter is 5.
</dd>

</details>
</dl>
</ul>



<ul>
<dl>
<details>
<summary>
<font size="+1"><b>GradientDomainProcessing</b></font>:
Supports the gradient domain smoothing and sharpening of surface geometry by solving a screened Poisson equation where the target values are given by the input geometry and the the target gradients are given by the dampened/amplified gradients of the input.
</summary>

<dt><b>--in</b> &lt;<i>input polygonal mesh</i>&gt;</dt>
<dd> This string specifies the the name of the input polygonal mesh.<br>
The polygonal mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i> and the set of polygons encoded by a list of vertex indices.
</dd>

<dt>[<b>--out</b> &lt;<i>output polygonal mesh</i>&gt;]</dt>
<dd> This string specifies the the name of the output (processed) polygonal mesh.<br>
The polygonal mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i> and the set of polygons encoded by a list of vertex indices.
</dd>

</dd><dt>[<b>--gWeight</B> &lt;<i>gradient interpolation weight</i>&gt;]</dt>
<dd> This floating point value specifies the weight that should be given to gradient interpolation.<BR>
The default value for this parameter is 1e-5.
</dd>

</dd><dt>[<b>--gScale</B> &lt;<i>gradient dampening/amplification factor</i>&gt;]</dt>
<dd> This floating point value specifies the scale that is to be appled to the gradients.<BR>
The default value for this parameter is 1, specifying unmodified output.
</dd>

</dd><dt>[<b>--mg</B>]</dt>
<dd> If enabled, this flag specifies that a multigrid solver should be used (instead of the default sparse Cholesky solver).
</dd>

<dt>[<b>--vCycles</b> &lt;<i>number of v-cycles per animation step</i>&gt;]</dt>
<dd> If a multigrid solver is used, ths integer value specifies the number of v-cycles to be performed at each step of the animation.<BR>
The default value for this parameter is 3.
</dd>

<dt>[<b>--gsIters</b> &lt;<i>number of Gauss-Seidel iterations per level</i>&gt;]</dt>
<dd> If a multigrid solver is used, this integer value specifies the number of Gauss-Seidel iterations to be done within each level of the v-cycle.<BR>
The default value for this parameter is 5.
</dd>

</dd><dt>[<b>--value</b> &lt;<i>data type</i>&gt;]</dt>
<dd> This integer value specifies which type of per-vertex-data is to be processed. If the value is <B>0</B>, the vertices' positions will be processed. If the value is <B>1</B>, the vertices' normals will be processed. If the value is <B>2</B>, the vertices' colors will be processed. (Note that if the value is <B>1</B> and the vertices do not have normals, the code will synthesize them by averaging from the indicent faces.)<BR>
The default value for this parameter is 0.
</dd>

</details>
</dl>
</ul>




<ul>
<dl>
<details>
<summary>
<font size="+1"><b>PrincipalCurvature</b></font>:
Supports the computation of per-face, extrinsic principal curvature directions by approximating the derivative of the normals across edges.
</summary>

<dt><b>--in</b> &lt;<i>input polygonal mesh</i>&gt;</dt>
<dd> This string specifies the the name of the input polygonal mesh.<br>
The polygonal mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i>, the x-, y-, and z-coordinates of the normals encoded by the properties <i>nx</i>, <i>ny</i>, and <i>nz</i>, and the set of polygons encoded by a list of vertex indices.
</dd>

<dt>[<b>--out</b> &lt;<i>output principal curvature directions</i>&gt;]</dt>
<dd> This string specifies the the name of the file to which the (extrinsic) principal curvature directions will be written, scaled by the associated curvature values.
</dd>

<dt>[<b>--kMin</b>]</dt>
<dd> If enabled, this flag specifies that the output tangent vectors will correspond to the directions of minimal curvature. Otherwise they will correspond to the directions of maximal curvature.
</dd>

</details>
</dl>
</ul>




<ul>
<dl>
<details>
<summary>
<font size="+1"><b>LineIntegralConvolution</b></font>:
Supports the generation of a line-integral-convolution visualization of a tangent vector field by diffusing a random scalar field with respect to an anisotropic metric that compresses distances along the the prescribed tangent directions.
</summary>

<dt><b>--in</b> &lt;<i>input polygonal mesh</i>&gt;</dt>
<dd> This string specifies the the name of the input polygonal mesh.<br>
The polygonal mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i> and the set of polygons encoded by a list of vertex indices.
</dd>

<dt><b>--inVF</b> &lt;<i>input tangent vector field</i>&gt;</dt>
<dd> This string specifies the the name of the file prescribing the per-face tangent vectors.
</dd>

<dt>[<b>--out</b> &lt;<i>output polygonal mesh</i>&gt;]</dt>
<dd> This string specifies the the name of the output polygonal mesh.<br>
The polygonal mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i>, the colors encoded by the properties <i>red</i>, <i>green</i>, and <i>blue</i>, and the set of polygons encoded by a list of vertex indices.
</dd>

</dd><dt>[<b>--mg</B>]</dt>
<dd> If enabled, this flag specifies that a multigrid solver should be used (instead of the default sparse Cholesky solver).
</dd>

<dt>[<b>--vCycles</b> &lt;<i>number of v-cycles per animation step</i>&gt;]</dt>
<dd> If a multigrid solver is used, ths integer value specifies the number of v-cycles to be performed at each step of the animation.<BR>
The default value for this parameter is 20.
</dd>

<dt>[<b>--gsIters</b> &lt;<i>number of Gauss-Seidel iterations per level</i>&gt;]</dt>
<dd> If a multigrid solver is used, ths integer value specifies the number of Gauss-Seidel iterations to be done within each level of the v-cycle.<BR>
The default value for this parameter is 5.
</dd>

</details>
</dl>
</ul>




<hr>
<a name="USAGE"><b>USAGE EXAMPLES (WITH SAMPLE DATA)</b></a><br>
For testing purposes, a number of <A HREF="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/VPSF.Data.zip">polgonal/polyhedral models</A> are provided (using the <U>.ply</U> and <U>.ovm</U> extensions respectively).

<ul>

<dl>
<details>
<summary>
<font size="+1"><b>Linear Elasticity in 2D</b></font>
</summary>
To run this executable you must specify the input polygonal mesh. For example, to see the deformation of the unit square, tessellated by a Voronoi diagram, deforming under the action of gravity, with the vertices on the left side locked, and using a direct solver to advance time-steps, execute:
<blockquote><code>% Bin/*/DeformableSolids2D --in ../VPSF.Data/square.voronoi.3.ply --lock</code></blockquote>
To see the deformation of the unit square, tessellated using concave polygons, evolving to its rest state after an initial anisotropic scaling is applied, using a hierarchical solver to advance time-steps, execute:
<blockquote><code>% Bin/*/DeformableSolids2D --in ../VPSF.Data/square.concave.3.ply --gravity 0 --xForm 2 0  0 0.5 --mg</code></blockquote>
You can toggle the animtation by hitting [SPACE] and you can step through the animation by hitting "+".<BR>
You can also interact with the animation by left-clicking to drag a vertex.
</details>
</dl>

<dl>
<details>
<summary>
<font size="+1"><b>Linear Elasticity in 3D</b></font>
</summary>
To run this executable you must specify the input polyhedral mesh. For example, to see the deformation of a unit cube, tessellated by a Voronoi diagram, deforming under the action of gravity, with the vertices on the left side locked, and using a direct solver to advance time-steps, execute:
<blockquote><code>% Bin/*/DeformableSolids3D --in ../VPSF.Data/cube.voronoi.3.ovm --lock</code></blockquote>
To see the deformation of the unit cube, tessellated using truncated cells, evolving to its rest state after an initial anisotropic scaling is applied, using a hierarchical solver to advance time-steps, execute:
<blockquote><code>% Bin/*/DeformableSolids3D --in ../VPSF.Data/cube.truncated.3.ovm --gravity 0 --xForm 2 0 0  0 1 0  0 0 0.5 --mg</code></blockquote>
You can toggle the animtation by hitting [SPACE] and you can step through the animation by hitting "+".<BR>
You can pan by by dragging with the left mouse button while holding down the [CTRL] key.<BR>
You can rotate by dragging with the left mouse button.<BR>
You can also rotate by using the "q", "w" , "a", "z", "s", and "x" keys.
</details>
</dl>

<dl>
<details>
<summary>
<font size="+1"><b>Franke Test in 2D</b></font>
</summary>
To run this executable you must specify the input polygonal mesh. For example, to run  the test on the unit square tessellated by a Voronoi diagram and using a direct solver, execute:
<blockquote><code>% Bin/*/FrankeTest2D --in ../VPSF.Data/square.voronoi.3.ply </code></blockquote>
To run the test on the unit square tessellated using concave polygons and using a hierarchical solver, execute:
<blockquote><code>% Bin/*/FrankeTest2D --in ../VPSF.Data/square.concave.3.ply --mg</code></blockquote>
</details>
</dl>

<dl>
<details>
<summary>
<font size="+1"><b>Franke Test in 3D</b></font>
</summary>
To run this executable you must specify the input polyhedral mesh. For example, to run  the test on the unit cube, tessellated by a Voronoi diagram, and using a direct solver, execute:
<blockquote><code>% Bin/*/FrankeTest3D --in ../VPSF.Data/cube.voronoi.3.ovm </code></blockquote>
To run the test on the unit cube, tessellated using truncated cells, and using a hierarchical solver, execute:
<blockquote><code>% Bin/*/FrankeTest3D --in ../VPSF.Data/cube.truncated.3.ovm --mg</code></blockquote>
</details>
</dl>


<dl>
<details>
<summary>
<font size="+1"><b>Geodesics in Heat</b></font>
</summary>
To run this executable you must specify the input polygonal mesh. For example, to visualize single-source geodesics on the model of the Armadillo Man, using a direct solver, execute:
<blockquote><code>% Bin/*/GeodesicsInHeat --in ../VPSF.Data/armadillo.ply</code></blockquote>
To visualize single-source geodesics on the model of the Fanblade, using a hierarchical solver, execute:
<blockquote><code>% Bin/*/GeodesicsInHeat --in ../VPSF.Data/fanblade.ply --mg</code></blockquote>
You can specify the geodesic source by left-clicking while holding down the [SHIFT] key.<BR>
You can pan by by dragging with the left mouse button while holding down the [CTRL] key.<BR>
You can rotate by dragging with the left mouse button.<BR>
You can also rotate by using the "q", "w" , "a", "z", "s", and "x" keys.
</details>
</dl>

<dl>
<details>
<summary>
<font size="+1"><b>Gradient Domain Processing</b></font>
</summary>
To run this executable you must specify the input and output polygonal meshes as well as the gradient interpolation weight and the gradient dampening/amplification scale. For example, to smooth the Bunny model using a direct solver, execute:
<blockquote><code>% Bin/*/GradientDomainProcessing --in ../VPSF.Data/bunny.ply --gScale 0 --out bunny.smooth.ply</code></blockquote>
For more aggressive smoothing, you can increase the gradient interpolation weight:
<blockquote><code>% Bin/*/GradientDomainProcessing --in ../VPSF.Data/bunny.ply --gScale 0 --out bunny.smooth.ply --gWeight 1e-3</code></blockquote>
To sharpen the Armadillo Man model using a hierarchical solver, execute:
<blockquote><code>% Bin/*/GradientDomainProcessing --in ../VPSF.Data/armadillo.ply --gScale 2 --out armadillo.sharp.ply --mg</code></blockquote>
</details>
</dl>


<dl>
<details>
<summary>
<font size="+1"><b>Line Integral Convolution</b></font>
</summary>
Generating a line-integral-convolution visualizing the maximal curvature directions requires three steps -- defining a smooth normal field, computing the maximal curvature direction, and diffusing a random noise signal using a metric defined by the maximal curvature directions. We provide the steps for synthesizing a visualization of the maximal curvature directions on the Fertility model, using a multigrid solver.
<OL>
<LI> To compute a smoothed normal field, execute:
<blockquote><code>% Bin/*/GradientDomainProcessing --in ../VPSF.Data/fertility.ply --gScale 0 --gWeight 1e-3 --out fertility.normal.ply --value 1 --mg</code></blockquote>
(For this example we smooth more aggressively, setting <i>--gWeight=1e-3</i>.)
<LI> To compute the maximal curvature directions, execute:
<blockquote><code>% Bin/*/PrincipalCurvature --in fertility.normal.ply --out fertility.kmax</code></blockquote>
<LI> To generate the line-integral-convolution visualization, execute:
<blockquote><code>% Bin/*/LineIntegralConvolution --in fertility.normal.ply --inVF fertility.kmax --out fertility.lic.ply --mg</code></blockquote>
</OL>
We make two observations about this executable. First, to compute the smooth normal field we diffuse more aggressively, using a gradient weight of <i>--gWeight=1e-3</i>, ten times larger than the default value. (This is neccessary for capturing a more global notion of curvature.) Second, the default number of v-cycles for the multigrid solver is 20. (This is necessary because the severe anisotropy makes the multigrid solver converge less efficiently.)

</details>
</dl>


</ul>

<hr>
<details>
<summary>
<a name="COMPILATION"><b>COMPILATION AND EXECUTION</b></a><br>
</summary>
<UL>
<LI>The Windows executables require both the <B>glew</B> and <B>glut</B> dynamically linked libraries to run. These can be found <A HREF="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/VPSF.DLLs.zip">here</A> and should be included either in the directory with the executables, or in the directory from which the executables are run.</LI>
<LI>Compiling under Windows requires both the <B>glew</B> and <B>glut</B> libraries. These can be found <A HREF="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/VPSF.LIBs.zip">here</A> and should be placed in the output directory for linkage.</LI></LI>
<LI> The implementation uses the standard <A HREF="https://eigen.tuxfamily.org/">Eigen</A> source code. If your system supports MKL it may be possible to develop faster executables by <CODE>#define</CODE>'ing the <B>EIGEN_USE_MKL_ALL</B> flag in <I>Misha/PreProcess.h</I> and ensuring that the associated include/library directories are set correctly.
</UL>
</details>

<hr>
<details>
<summary>
<a name="CHANGES"><b>HISTORY OF CHANGES</b></a><br>
</summary>
<a href="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/Version2.00/">Version 2</a>:
<ul><li> Added support for regularizing within the kernel.</li></ul>
<a href="http://www.cs.jhu.edu/~misha/Code/VariationalPolyShapeFunctions/Version2.50/">Version 2.5</a>:
<ul><li> Added line-integral-convolution code.</li></ul>
</details>


<!--
<hr>
<a name="SUPPORT"><b>SUPPORT</b></a><br>
-->

<hr>
<a href="http://www.cs.jhu.edu/~misha">HOME</a>
