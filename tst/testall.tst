#############################################################################
##
#W  testall.tst                                              Oeyvind Solberg
##
##
#Y  Copyright (C)  2015,  The QPA-team
##
##  
##
gap> START_TEST("qpa.tst");

# Entering a quiver.
gap> Q:= Quiver(3,[[1,2,"a"],[1,2,"b"],[2,2,"c"],[2,3,"d"],[3,1,"e"]]);
<quiver with 3 vertices and 5 arrows>

# Listing the vertices in the quiver.
gap> VerticesOfQuiver(Q);
[ v1, v2, v3 ]

# Listing the arrows in the quiver.
gap> ArrowsOfQuiver(Q);
[ a, b, c, d, e ]

#
# Constructing the opposite quiver.
gap> Qop := OppositeQuiver(Q);
<quiver with 3 vertices and 5 arrows>

#
# Listing the vertices.
gap> VerticesOfQuiver(Qop);
[ v1_op, v2_op, v3_op ]

#
# Listing the arrows.
gap> ArrowsOfQuiver(Qop);
[ a_op, b_op, c_op, d_op, e_op ]

#
# Constructing the path algebra over this quiver over the field of rationals.
gap> KQ:= PathAlgebra(Rationals, Q);
<Rationals[<quiver with 3 vertices and 5 arrows>]>

#
# The arrows are called a, b, ... , in the definition of the quiver, but they do not
# exist as variables in GAP yet.  The following defines these and variables for the
# vertices. 
gap> AssignGeneratorVariables(KQ);
#I  Assigned the global variables [ v1, v2, v3, a, b, c, d, e ]

#
# Defining a set of relations.
gap> rels:= [d*e,c^2,a*c*d-b*d,e*a];
[ (1)*d*e, (1)*c^2, (-1)*b*d+(1)*a*c*d, (1)*e*a ]

#
# Constructing the quotient of the path algebra.
gap> A:= KQ/rels;
<Rationals[<quiver with 3 vertices and 5 arrows>]/
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (5 generators)>>
gap> Dimension(A);
17

#
# Constructing the tenor product of  A  with itself.
gap> TensorProductOfAlgebras(A,A);
<Rationals[<quiver with 9 vertices and 30 arrows>]/
<two-sided ideal in <Rationals[<quiver with 9 vertices and 30 arrows>]>, 
  (55 generators)>>

#
# Defining a module over the original algebra  A.
gap> mat := [["a",[[1,2],[0,3],[1,5]]],["b",[[2,0],[3,0],[5,0]]],
> ["c",[[0,0],[1,0]]],["d",[[1,2],[0,1]]],["e",[[0,0,0],[0,0,0]]]];
[ [ "a", [ [ 1, 2 ], [ 0, 3 ], [ 1, 5 ] ] ], 
  [ "b", [ [ 2, 0 ], [ 3, 0 ], [ 5, 0 ] ] ], [ "c", [ [ 0, 0 ], [ 1, 0 ] ] ], 
  [ "d", [ [ 1, 2 ], [ 0, 1 ] ] ], [ "e", [ [ 0, 0, 0 ], [ 0, 0, 0 ] ] ] ]
gap> N := RightModuleOverPathAlgebra(A,mat);  
<[ 3, 2, 2 ]>

#
# Listing the defining matrices for the module  N.
gap> MatricesOfPathAlgebraModule(N);
[ [ [ 1, 2 ], [ 0, 3 ], [ 1, 5 ] ], [ [ 2, 0 ], [ 3, 0 ], [ 5, 0 ] ], 
  [ [ 0, 0 ], [ 1, 0 ] ], [ [ 1, 2 ], [ 0, 1 ] ], 
  [ [ 0, 0, 0 ], [ 0, 0, 0 ] ] ]

#
# Recalling the dimension of  N.
gap> DimensionVector(N);
[ 3, 2, 2 ]

#
# Finding the radical, top, socle of the module, projective cover.
gap> RadicalOfModuleInclusion(N);
<<[ 0, 2, 2 ]> ---> <[ 3, 2, 2 ]>>

gap> TopOfModuleProjection(N);   
<<[ 3, 2, 2 ]> ---> <[ 3, 0, 0 ]>>

gap> SocleOfModuleInclusion(N);  
<<[ 1, 0, 2 ]> ---> <[ 3, 2, 2 ]>>

gap> p := ProjectiveCover(N);
<<[ 3, 12, 9 ]> ---> <[ 3, 2, 2 ]>>


# Constructing the direct sum of  N  with itself and finding naturally
# associated inclusions and projections.  
gap> M := DirectSumOfQPAModules([N,N]);
<[ 6, 4, 4 ]>
gap> DirectSumInclusions(M);
[ <<[ 3, 2, 2 ]> ---> <[ 6, 4, 4 ]>>
    , <<[ 3, 2, 2 ]> ---> <[ 6, 4, 4 ]>>
     ]
gap> DirectSumProjections(M);
[ <<[ 6, 4, 4 ]> ---> <[ 3, 2, 2 ]>>
    , <<[ 6, 4, 4 ]> ---> <[ 3, 2, 2 ]>>
     ]
gap> IsDirectSummand(N,M);
true
gap> IsomorphicModules(N,M);
false

#
# Construting the simple module  S1  in two different ways.
gap> S1 := SimpleModules(A)[1];
<[ 1, 0, 0 ]>
gap> U1 := RightModuleOverPathAlgebra(A,[1,0,0],[]);
<[ 1, 0, 0 ]>
gap> S1 = U1;
true
gap> IsIdenticalObj(S1,U1);
false

#
# Testing various operations on modules.
gap> IsExceptionalModule(S1);
true
gap> AnnihilatorOfModule(N);
[ [(1)*e], [(-1)*b+(1)*a*c], [(1)*b*c], [(1)*e*b], [(1)*b*c*d], [(1)*e*b*c], 
  [(1)*e*b*c*d] ]
gap> BlockDecompositionOfModule(N);
[ <[ 2, 2, 2 ]>, <[ 1, 0, 0 ]> ]
gap> IsProjectiveModule(N);
false
gap> ProjDimensionOfModule(N,10);
4
gap> IsRigidModule(N);
true
gap> IsSemisimpleModule(N);
false
gap> IsSimpleQPAModule(N);
false
gap> IsSimpleQPAModule(S1);
true
gap> gens := MinimalGeneratingSetOfModule(N);
[ [ [ 1, 0, 0 ], [ 0, 0 ], [ 0, 0 ] ], [ [ 0, 1, 0 ], [ 0, 0 ], [ 0, 0 ] ], 
  [ [ 0, 0, 1 ], [ 0, 0 ], [ 0, 0 ] ] ]
gap> f1 := SubRepresentationInclusion(N,[gens[1]]);
<<[ 1, 2, 2 ]> ---> <[ 3, 2, 2 ]>>

gap> f2 := SubRepresentationInclusion(N,[gens[2]]);
<<[ 1, 2, 2 ]> ---> <[ 3, 2, 2 ]>>

gap> SumOfSubmodules(f1,f2);
[ <<[ 2, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
    , <<[ 1, 2, 2 ]> ---> <[ 2, 2, 2 ]>>
    , <<[ 1, 2, 2 ]> ---> <[ 2, 2, 2 ]>>
     ]
gap> RadicalSeries(N);
[ [ 3, 0, 0 ], [ 0, 1, 0 ], [ 0, 1, 1 ], [ 0, 0, 1 ] ]
gap> SocleSeries(N);
[ [ 1, 0, 0 ], [ 1, 1, 0 ], [ 0, 1, 0 ], [ 1, 0, 2 ] ]
gap> CommonDirectSummand(N,N);
[ <[ 1, 0, 0 ]>, <[ 2, 2, 2 ]>, <[ 1, 0, 0 ]>, <[ 2, 2, 2 ]> ]
gap> IndecProjectiveModules(A);
[ <[ 1, 4, 3 ]>, <[ 0, 2, 2 ]>, <[ 1, 2, 2 ]> ]
gap> IndecInjectiveModules(A);
[ <[ 1, 0, 1 ]>, <[ 4, 2, 2 ]>, <[ 3, 2, 2 ]> ]
gap> DualOfModule(N);
<[ 3, 2, 2 ]>
gap> DTr(N);
<[ 17, 10, 7 ]>
gap> NakayamaFunctorOfModule(N);
<[ 0, 0, 0 ]>
gap> StarOfModule(N);
<[ 0, 0, 0 ]>
gap> TransposeOfModule(N);
<[ 17, 10, 7 ]>

#
# Testing various operations for homomorphisms of modules
gap> hom := HomOverAlgebra(N,N);
[ <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
    , <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
    , <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
    , <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
    , <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
     ]
gap> f := hom[1];
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>

gap> IsInjective(f);
false
gap> IsIsomorphism(f);
false
gap> IsLeftMinimal(f);
false
gap> IsRightMinimal(f);
false
gap> IsSplitEpimorphism(f);
false
gap> IsSplitMonomorphism(f);
false
gap> IsSurjective(f);
false
gap> IsZero(f);
false
gap> LeftInverseOfHomomorphism(f);
false
gap> CoKernel(f);
<[ 2, 2, 2 ]>
gap> CoKernelProjection(f);
<<[ 3, 2, 2 ]> ---> <[ 2, 2, 2 ]>>

gap> EndOverAlgebra(N);
<algebra-with-one of dimension 5 over Rationals>
gap> Image(f);
<[ 1, 0, 0 ]>
gap> ImageInclusion(f);
<<[ 1, 0, 0 ]> ---> <[ 3, 2, 2 ]>>

gap> ImageProjection(f);
<<[ 3, 2, 2 ]> ---> <[ 1, 0, 0 ]>>

gap> Kernel(f);
<[ 2, 2, 2 ]>
gap> LeftMinimalVersion(f);
[ <<[ 3, 2, 2 ]> ---> <[ 1, 0, 0 ]>>
    , [ <[ 2, 2, 2 ]> ] ]
gap> RightMinimalVersion(f);
[ <<[ 1, 0, 0 ]> ---> <[ 3, 2, 2 ]>>
    , [ <[ 2, 2, 2 ]> ] ]

#
# Testing homological algebra operations.
gap> CotiltingModule(N,2);
false
gap> ExtOverAlgebra(N,N);
[ <<[ 0, 10, 7 ]> ---> <[ 3, 12, 9 ]>>
    , [  ], [  ] ]
gap> GlobalDimensionOfAlgebra(A,4);
infinity
gap> GorensteinDimensionOfAlgebra(A,5);
4
gap> InjDimensionOfModule(S1,4);  
3
gap> P:=DirectSumOfQPAModules(IndecProjectiveModules(A));
<[ 2, 8, 7 ]>
gap> TiltingModule(P,1);                                              
[ 0, [ 0 -> 0:(1,4,3) -> -1:(1,4,3) -> 0, 0 -> 0:(0,2,2) -> -1:(0,2,2) -> 0, 
      0 -> 0:(1,2,2) -> -1:(1,2,2) -> 0 ] ]
gap> CotiltingModule(P,3);
false
gap> MinimalRightAddMApproximation(N,S1);
<<[ 1, 0, 0 ]> ---> <[ 1, 0, 0 ]>>

gap> MinimalLeftAddMApproximation(S1,N);   
<<[ 1, 0, 0 ]> ---> <[ 1, 0, 0 ]>>

gap> NthSyzygy(S1,4);
Computing syzygy number: 1 ...
Dimension vector for syzygy: [ 0, 4, 3 ]
Top of the 1th syzygy: [ 0, 2, 0 ]
Computing syzygy number: 2 ...
Dimension vector for syzygy: [ 0, 0, 1 ]
Top of the 2th syzygy: [ 0, 0, 1 ]
Computing syzygy number: 3 ...
Dimension vector for syzygy: [ 1, 2, 1 ]
Top of the 3th syzygy: [ 1, 0, 0 ]
Computing syzygy number: 4 ...
Dimension vector for syzygy: [ 0, 2, 2 ]
Top of the 4th syzygy: [ 0, 1, 0 ]
The module has projective dimension 4.
<[ 0, 2, 2 ]>
gap> ProjDimensionOfModule(S1,4);
4
gap> ProjectiveResolutionOfPathAlgebraModule(S1,4);
Computing f2's ...
Computing f3's ...
Computing f4's ...
ProjectiveResolutionOfPathAlgebraModule
gap> PullBack(hom[1],hom[2]);
[ <<[ 5, 4, 4 ]> ---> <[ 3, 2, 2 ]>>
    , <<[ 5, 4, 4 ]> ---> <[ 3, 2, 2 ]>>
     ]
gap> PushOut(hom[1],hom[2]);
[ <<[ 3, 2, 2 ]> ---> <[ 4, 4, 4 ]>>
    , <<[ 3, 2, 2 ]> ---> <[ 4, 4, 4 ]>>
     ]

#
# Testing some AR-theory.
gap> AlmostSplitSequence(S1);
[ <<[ 7, 4, 3 ]> ---> <[ 8, 4, 3 ]>>
    , <<[ 8, 4, 3 ]> ---> <[ 1, 0, 0 ]>>
     ]

#
# Testing chain complexes
gap> alg := A;
<Rationals[<quiver with 3 vertices and 5 arrows>]/
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (5 generators)>>
gap> f := hom[1];
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>

gap> h := CoKernelProjection(f);
<<[ 3, 2, 2 ]> ---> <[ 2, 2, 2 ]>>

gap> g := h*HomOverAlgebra(Range(h),N)[2];
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>

gap> L := ShallowCopy(N);
<[ 3, 2, 2 ]>
gap> M := ShallowCopy(N);
<[ 3, 2, 2 ]>
gap> cat := CatOfRightAlgebraModules(alg);
<cat: right modules over algebra>
gap> cat.zeroObj;
<[ 0, 0, 0 ]>
gap> cat.isZeroObj(M);
false
gap> cat.zeroMap(M,N);
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>

gap> cat.composeMaps(g,f);
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>

gap> cat.ker(g);
<[ 2, 2, 2 ]>
gap> cat.isExact(g,f);
false
gap> C := FiniteComplex(cat, 1, [g,f]);
0 -> 2:(3,2,2) -> 1:(3,2,2) -> 0:(3,2,2) -> 0
gap> Ms := StalkComplex(cat, M, 3);
0 -> 3:(3,2,2) -> 0
gap> IsFiniteComplex(C);
true
gap> UpperBound(C);
2
gap> LowerBound(C);
0
gap> LengthOfComplex(C);
3
gap> HighestKnownDegree(C);
+inf
gap> LowestKnownDegree(C);
-inf
gap> IsExactSequence(C);
false
gap> IsExactInDegree(C,1);
false
gap> IsExactInDegree(C,2);
false
gap> Shift(C,1);
0 -> 1:(3,2,2) -> 0:(3,2,2) -> -1:(3,2,2) -> 0
gap> D := Shift(C,-1);
0 -> 3:(3,2,2) -> 2:(3,2,2) -> 1:(3,2,2) -> 0
gap> dc := DifferentialOfComplex(C,3)!.maps;
[ [ [ 0, 0, 0 ] ], [ [ 0, 0 ] ], [ [ 0, 0 ] ] ]
gap> dd := DifferentialOfComplex(D,4)!.maps;
[ [ [ 0, 0, 0 ] ], [ [ 0, 0 ] ], [ [ 0, 0 ] ] ]

#
# Testing unit forms.
gap> u := TitsUnitFormOfAlgebra(NakayamaAlgebra(Rationals,[3,2,1]));
[ [ 2, -1, 0 ], [ -1, 2, -1 ], [ 0, -1, 2 ] ]
gap> IsWeaklyNonnegativeUnitForm(u);
true
gap> IsWeaklyPositiveUnitForm(u);
true
gap> PositiveRootsOfUnitForm(u);
[ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 1, 0 ], [ 0, 1, 1 ], 
  [ 1, 1, 1 ] ]
gap> BilinearFormOfUnitForm(u);
function( x, y ) ... end
gap> QuadraticFormOfUnitForm(u);
function( x ) ... end
gap> SymmetricMatrixOfUnitForm(u);
[ [ 2, -1, 0 ], [ -1, 2, -1 ], [ 0, -1, 2 ] ]
