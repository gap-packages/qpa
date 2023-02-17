############################################################################
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
gap> IsQuiver(Q);
true
gap> dq := DynkinQuiver( "A", 4, ["r","l","r"]);
<quiver with 4 vertices and 3 arrows>
gap> IsAcyclicQuiver(dq);
true
gap> IsUAcyclicQuiver(dq);
true
gap> IsConnectedQuiver(dq);
true
gap> IsTreeQuiver(dq);
true
gap> IsDynkinQuiver(dq);
true

# Listing the vertices in the quiver.
gap> VerticesOfQuiver(Q);
[ v1, v2, v3 ]

# Listing the arrows in the quiver.
gap> ArrowsOfQuiver(Q);
[ a, b, c, d, e ]
gap> Q.a; Q.b; Q.v1;
a
b
v1
gap> AdjacencyMatrixOfQuiver(Q);
[ [ 0, 2, 0 ], [ 0, 1, 1 ], [ 1, 0, 0 ] ]
gap> GeneratorsOfQuiver(Q);
[ v1, v2, v3, a, b, c, d, e ]
gap> NumberOfVertices(Q);
3
gap> NumberOfArrows(Q);
5
gap> OrderingOfQuiver(Q);
<length left lexicographic ordering>

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
gap> Qsub := FullSubquiver(Q, [Q.v1, Q.v2]);
<quiver with 2 vertices and 3 arrows>
gap> Q2 := Quiver(4, [[1,2,"a"], [3,4,"b"]]);
<quiver with 4 vertices and 2 arrows>
gap> ConnectedComponentsOfQuiver(Q2);
[ <quiver with 2 vertices and 1 arrows>, 
  <quiver with 2 vertices and 1 arrows> ]
gap> Qsep := SeparatedQuiver(Q);
<quiver with 6 vertices and 5 arrows>
gap> p := Q.a*Q.c;
a*c
gap> IsPath(p);
true
gap> IsQuiverVertex(Q.v1);
true
gap> IsArrow(Q.a);
true
gap> IsZeroPath(Q.v1);
false
gap> SourceOfPath(Q.a);
v1
gap> TargetOfPath(Q.b);
v2
gap> LengthOfPath(p);
2
gap> WalkOfPath(p);
[ a, c ]
gap> p = Q.a*Q.c;
true
gap> Q.a < p;
true
gap> IncomingArrowsOfVertex(Q.v1);
[ e ]
gap> OutgoingArrowsOfVertex(Q.v1);
[ a, b ]
gap> InDegreeOfVertex(Q.v1);
1
gap> OutDegreeOfVertex(Q.v1);
2
gap> NeighborsOfVertex(Q.v1);
[ v2 ]

#
# Constructing the path algebra over this quiver over the field of rationals.
gap> KQ:= PathAlgebra(Rationals, Q);
<Rationals[<quiver with 3 vertices and 5 arrows>]>
gap> IsPathAlgebra(KQ);             
true

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
# Constructing the path algebra but now using HighLevelGroebnerBasis.
gap> KQ:= PathAlgebra(Rationals, Q, HighLevelGroebnerBasis);
<Rationals[<quiver with 3 vertices and 5 arrows>]>
gap> AssignGeneratorVariables(KQ);
#I  Global variable `v1' is already defined and will be overwritten
#I  Global variable `v2' is already defined and will be overwritten
#I  Global variable `v3' is already defined and will be overwritten
#I  Global variable `a' is already defined and will be overwritten
#I  Global variable `b' is already defined and will be overwritten
#I  Global variable `c' is already defined and will be overwritten
#I  Global variable `d' is already defined and will be overwritten
#I  Global variable `e' is already defined and will be overwritten
#I  Assigned the global variables [ v1, v2, v3, a, b, c, d, e ]
gap> rels:= [d*e,c^2,a*c*d-b*d,e*a];
[ (1)*d*e, (1)*c^2, (-1)*b*d+(1)*a*c*d, (1)*e*a ]
gap> A:= KQ/rels;
<Rationals[<quiver with 3 vertices and 5 arrows>]/
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (5 generators)>>
gap> Dimension(A);
17

#
# Testing HighLevelGroebnerBasis on the first example from issue #20.
gap> Q_ := Quiver(1, [[1,1,"x"], [1,1,"y"]]);
<quiver with 1 vertices and 2 arrows>
gap> KQ_ := PathAlgebra(Rationals, Q_, HighLevelGroebnerBasis);
<Rationals[<quiver with 1 vertices and 2 arrows>]>
gap> x := KQ_.x;; y := KQ_.y;;
gap> rels_ := [ x^3, x^2*y, x*y^2+x^2, y*x^2, x^3+y^2+y*x, x*y^2+x*y*x, y*x*y+y*x^2, y^2*x+y*x^2, y^3+y*x^2, x^2*y*x, x*y*x^2, y*x^2*y ];
[ (1)*x^3, (1)*x^2*y, (1)*x^2+(1)*x*y^2, (1)*y*x^2, (1)*y*x+(1)*y^2+(1)*x^3, 
  (1)*x*y*x+(1)*x*y^2, (1)*y*x^2+(1)*y*x*y, (1)*y*x^2+(1)*y^2*x, 
  (1)*y*x^2+(1)*y^3, (1)*x^2*y*x, (1)*x*y*x^2, (1)*y*x^2*y ]
gap> A_ := KQ_/rels_;
<Rationals[<quiver with 1 vertices and 2 arrows>]/
<two-sided ideal in <Rationals[<quiver with 1 vertices and 2 arrows>]>, 
  (6 generators)>>
gap> RelationsOfAlgebra(A_);
[ (1)*y*x+(1)*y^2, (1)*x^3, (1)*x^2*y, (-1)*x^2+(1)*x*y*x, (1)*y*x^2, 
  (1)*y*x*y ]

#
#
gap> B := AssociatedMonomialAlgebra(A);
<Rationals[<quiver with 3 vertices and 5 arrows>]/
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (5 generators)>>
gap> qb := QuiverOfPathAlgebra(B);
<quiver with 3 vertices and 5 arrows>
gap> OrderingOfAlgebra(B);
<length left lexicographic ordering>
gap> KQ.a; A.a;
(1)*a
[(1)*a]
gap> aa := ElementOfPathAlgebra(KQ, Q.a);
(1)*a
gap> KQ.a = aa;
true
gap> KQ.b < KQ.a*KQ.c;
true
gap> IsLeftUniform(KQ.b);
true
gap> IsRightUniform(KQ.b);
true
gap> IsUniform(KQ.b);
true
gap> LeadingTerm(2*KQ.b + 3*KQ.a*KQ.c);
(3)*a*c
gap> LeadingCoefficient(2*KQ.b + 3*KQ.a*KQ.c);
3
gap> LeadingMonomial(2*KQ.b + 3*KQ.a*KQ.c);   
a*c
gap> MakeUniformOnRight([2*KQ.b + 3*KQ.a*KQ.c, KQ.c + KQ.e]);
[ (2)*b+(3)*a*c, (1)*c, (1)*e ]
gap> VertexPosition(KQ.v2);
2
gap> relations := RelationsOfAlgebra(A);
[ (1)*c^2, (1)*d*e, (1)*e*a, (-1)*b*d+(1)*a*c*d, (1)*e*b*d ]
gap> I := Ideal(KQ, relations);
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (5 generators)>
gap> II := IdealOfQuotient(A);
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (5 generators)>
gap> PathsOfLengthTwo(Q);
[ a*c, a*d, b*c, b*d, c^2, c*d, d*e, e*a, e*b ]
gap> NthPowerOfArrowIdeal(KQ,2);
[ (1)*a*c, (1)*a*d, (1)*b*c, (1)*b*d, (1)*c^2, (1)*c*d, (1)*d*e, (1)*e*a, 
  (1)*e*b ]
gap> newrels := ShallowCopy(relations);
[ (1)*c^2, (1)*d*e, (1)*e*a, (-1)*b*d+(1)*a*c*d, (1)*e*b*d ]
gap> AddNthPowerToRelations(KQ, newrels, 2);  
[ (1)*c^2, (1)*d*e, (1)*e*a, (-1)*b*d+(1)*a*c*d, (1)*e*b*d, (1)*a*c, (1)*a*d, 
  (1)*b*c, (1)*b*d, (1)*c^2, (1)*c*d, (1)*d*e, (1)*e*a, (1)*e*b ]
gap> newrels;
[ (1)*c^2, (1)*d*e, (1)*e*a, (-1)*b*d+(1)*a*c*d, (1)*e*b*d, (1)*a*c, (1)*a*d, 
  (1)*b*c, (1)*b*d, (1)*c^2, (1)*c*d, (1)*d*e, (1)*e*a, (1)*e*b ]
gap> Q.a in KQ;
false
gap> KQ.a in KQ;
true
gap> IsAdmissibleIdeal(I);
true
gap> IsIdealInPathAlgebra(I);    
true
gap> IsMonomialIdeal(I);
false
gap> IsQuadraticIdeal(newrels);
false
gap> quad := NthPowerOfArrowIdeal(KQ,2);
[ (1)*a*c, (1)*a*d, (1)*b*c, (1)*b*d, (1)*c^2, (1)*c*d, (1)*d*e, (1)*e*a, 
  (1)*e*b ]
gap> QuadraticPerpOfPathAlgebraIdeal(quad);
[ <Rationals[<quiver with 3 vertices and 5 arrows>]>, [  ] ]
gap> GroebnerBasisOfIdeal(I);
<complete two-sided Groebner basis containing 5 elements>
gap> IsAdmissibleQuotientOfPathAlgebra(A);
true
gap> IsQuotientOfPathAlgebra(A);
true
gap> IsFiniteDimensional(A);
true
gap> IsDistributiveAlgebra(A);
false
gap> IsGentleAlgebra(A);
false
gap> IsHereditaryAlgebra(A);
false
gap> IsMonomialAlgebra(A);
false
gap> IsQuiverAlgebra(A);
true
gap> IsRadicalSquareZeroAlgebra(A);
false
gap> IsSchurianAlgebra(A);
false
gap> IsSelfinjectiveAlgebra(A);
false
gap> IsSemicommutativeAlgebra(A);
false
gap> IsSemisimpleAlgebra(A);
false
gap> IsSpecialBiserialAlgebra(A);
false
gap> IsStringAlgebra(A);
false
gap> IsSymmetricAlgebra(A);
false
gap> IsTriangularReduced(A);
true
gap> IsWeaklySymmetricAlgebra(A);
false
gap> IsFiniteTypeAlgebra(A);
Quiver contains an (un)oriented cycle.
Quiver contains an (un)oriented cycle.
Have checked representation type of the algebra modulo the square of the radic\
al.
false
gap> Q1 := Quiver(4, [[1,2,"a"], [1,2,"b"], [2,3,"c"], [3,4,"d"], [3,4,"e"]]);
<quiver with 4 vertices and 5 arrows>
gap> k := Rationals;
Rationals
gap> kQ1 := PathAlgebra(k,Q1);
<Rationals[<quiver with 4 vertices and 5 arrows>]>
gap> rho := [kQ1.b*kQ1.c, kQ1.c*kQ1.d];
[ (1)*b*c, (1)*c*d ]
gap> A1 := kQ1/rho;
<Rationals[<quiver with 4 vertices and 5 arrows>]/
<two-sided ideal in <Rationals[<quiver with 4 vertices and 5 arrows>]>, (2 generators)>>
gap> IsValidString(A1,"eca");
true
gap> StringsLessThan(A1, 2);
[ "(1,1)", "b", "Ab", "(1,-1)", "a", "ca", "Ba", "(2,1)", "c", "B", "ec", "aB",
"(2,-1)", "A", "bA", "(3,1)", "e", "De", "(3,-1)", "d", "C", "Ed", "AC",
"(4,1)", "E", "dE", "CE", "(4,-1)", "D", "eD" ]
gap> IsABand(A1,"eca");
false
gap> IsABand(A1,"Ab");
true
gap> BandsLessThan(A1,3);
[ "Ab", "Ba", "aB", "bA", "De", "Ed", "dE", "eD" ]
gap> BandRepresentativesLessThan(A1,3);
[ "aB", "bA", "dE", "eD" ]
gap> IsDomesticStringAlgebra(A1);
true
gap> Q2 := BridgeQuiver(A1);
<quiver with 4 vertices and 2 arrows>
gap> Display(Q2);
Quiver( ["v1","v2","v3","v4"], [["v3","v2","CE"],["v1","v4","ec"]] )
gap> LocalARQuiver(A1,"eca");
[ "The canonical embedding StringModule(eca) to StringModule(eDeca)",
  "The canonical projection StringModule(ecaBa) to StringModule(eca)" ]
gap> CartanMatrix(A);
[ [ 1, 4, 3 ], [ 0, 2, 2 ], [ 1, 2, 2 ] ]
gap> Center(A);
<algebra-with-one of dimension 2 over Rationals>
gap> ComplexityOfAlgebra(A,10);
2
gap> CoxeterMatrix(A);
[ [ 1, 0, 0 ], [ 4, 3, 2 ], [ -6, -4, -3 ] ]
gap> CoxeterPolynomial(A);
x_1^3-x_1^2-x_1+1
gap> Dimension(A);
17
gap> LoewyLength(A);
5
gap> Nak := NakayamaAlgebra(Rationals, [3,3,3]);
<Rationals[<quiver with 3 vertices and 3 arrows>]/
<two-sided ideal in <Rationals[<quiver with 3 vertices and 3 arrows>]>, 
  (3 generators)>>
gap> IsNakayamaAlgebra(Nak);
true
gap> NakayamaAutomorphism(Nak);                 
function( x ) ... end
gap> NakayamaPermutation(Nak);                  
[ function( x ) ... end, function( x ) ... end ]
gap> OrderOfNakayamaAutomorphism(Nak);
3
gap> ComplexityOfModule(SimpleModules(Nak)[1],10);
1
gap> RadicalSeriesOfAlgebra(A);
[ <Rationals[<quiver with 3 vertices and 5 arrows>]/
    <two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
      (5 generators)>>, <algebra of dimension 14 over Rationals>, 
  <algebra of dimension 9 over Rationals>, 
  <algebra of dimension 4 over Rationals>, 
  <algebra of dimension 1 over Rationals>, 
  <algebra of dimension 0 over Rationals> ]
gap> IsElementOfQuotientOfPathAlgebra(A.a);
true
gap> OriginalPathAlgebra(A);
<Rationals[<quiver with 3 vertices and 5 arrows>]>
gap> CanonicalAlgebra(Rationals, [2,3,4], [-1]);
<Rationals[<quiver with 8 vertices and 9 arrows>]/
<two-sided ideal in <Rationals[<quiver with 8 vertices and 9 arrows>]>, 
  (1 generators)>>
gap> kron := KroneckerAlgebra(Rationals, 3);
<Rationals[<quiver with 2 vertices and 3 arrows>]>
gap> IsKroneckerAlgebra(kron);
true
gap> NakayamaAlgebra(Rationals, [3,3,2,1]);
<Rationals[<quiver with 4 vertices and 3 arrows>]/
<two-sided ideal in <Rationals[<quiver with 4 vertices and 3 arrows>]>, 
  (1 generators)>>
gap> TruncatedPathAlgebra(Rationals,Q,2); 
<Rationals[<quiver with 3 vertices and 5 arrows>]/
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (9 generators)>>
gap> IsSpecialBiserialQuiver(Q);
false
gap> OppositePath(Q.a*Q.c);
c_op*a_op
gap> OppositePathAlgebra(KQ);
<Rationals[<quiver with 3 vertices and 5 arrows>]>
gap> OppositePathAlgebraElement(KQ.a);
(1)*a_op
gap> qNak := QuiverOfPathAlgebra(Nak);      
<quiver with 3 vertices and 3 arrows>
gap> prodqNak := QuiverProduct(qNak,qNak);  
<quiver with 9 vertices and 18 arrows>
gap> TensorProductOfAlgebras(Nak,Nak);
<Rationals[<quiver with 9 vertices and 18 arrows>]/
<two-sided ideal in <Rationals[<quiver with 9 vertices and 18 arrows>]>, 
  (27 generators)>>
gap> T := TensorProductOfAlgebras(Nak,Nak);
<Rationals[<quiver with 9 vertices and 18 arrows>]/
<two-sided ideal in <Rationals[<quiver with 9 vertices and 18 arrows>]>, 
  (27 generators)>>
gap> TensorAlgebraInclusion(T,1);          
[ [(1)*v1+(1)*v2+(1)*v3], [(1)*v1], [(1)*v2], [(1)*v3], [(1)*a1], [(1)*a2], 
  [(1)*a3] ] -> 
[ [(1)*v1_v1+(1)*v1_v2+(1)*v1_v3+(1)*v2_v1+(1)*v2_v2+(1)*v2_v3+(1)*v3_v1+(
    1)*v3_v2+(1)*v3_v3], [(1)*v1_v1+(1)*v1_v2+(1)*v1_v3], 
  [(1)*v2_v1+(1)*v2_v2+(1)*v2_v3], [(1)*v3_v1+(1)*v3_v2+(1)*v3_v3], 
  [(1)*a1_v1+(1)*a1_v2+(1)*a1_v3], [(1)*a2_v1+(1)*a2_v2+(1)*a2_v3], 
  [(1)*a3_v1+(1)*a3_v2+(1)*a3_v3] ]
gap> SimpleTensor([Nak.v1, Nak.a1],T);
[(1)*v1_a1]
gap> TensorProductDecomposition(T);
[ <Rationals[<quiver with 3 vertices and 3 arrows>]/
    <two-sided ideal in <Rationals[<quiver with 3 vertices and 3 arrows>]>, 
      (3 generators)>>, <Rationals[<quiver with 3 vertices and 3 arrows>]/
    <two-sided ideal in <Rationals[<quiver with 3 vertices and 3 arrows>]>, 
      (3 generators)>> ]
gap> Nakenv := EnvelopingAlgebra(Nak);
<Rationals[<quiver with 9 vertices and 18 arrows>]/
<two-sided ideal in <Rationals[<quiver with 9 vertices and 18 arrows>]>, 
  (27 generators)>>
gap> IsEnvelopingAlgebra(Nakenv);
true
gap> AlgebraAsModuleOverEnvelopingAlgebra(Nak);
<[ 1, 1, 1, 1, 1, 1, 1, 1, 1 ]>
gap> TrivialExtensionOfQuiverAlgebra(Nak);
<Rationals[<quiver with 3 vertices and 6 arrows>]/
<two-sided ideal in <Rationals[<quiver with 3 vertices and 6 arrows>]>, 
  (9 generators)>>

#
# Defining a module over the original algebra  A.
gap> mat := [["a",[[1,2],[0,3],[1,5]]],["b",[[2,0],[3,0],[5,0]]],
> ["c",[[0,0],[1,0]]],["d",[[1,2],[0,1]]],["e",[[0,0,0],[0,0,0]]]];
[ [ "a", [ [ 1, 2 ], [ 0, 3 ], [ 1, 5 ] ] ], 
  [ "b", [ [ 2, 0 ], [ 3, 0 ], [ 5, 0 ] ] ], [ "c", [ [ 0, 0 ], [ 1, 0 ] ] ], 
  [ "d", [ [ 1, 2 ], [ 0, 1 ] ] ], [ "e", [ [ 0, 0, 0 ], [ 0, 0, 0 ] ] ] ]
gap> N := RightModuleOverPathAlgebra(A,mat);  
<[ 3, 2, 2 ]>
gap> IsPathAlgebraMatModule(N);
true
gap> Dimension(N);
7

#
# Listing the defining matrices for the module  N.
gap> MatricesOfPathAlgebraModule(N);
[ [ [ 1, 2 ], [ 0, 3 ], [ 1, 5 ] ], [ [ 2, 0 ], [ 3, 0 ], [ 5, 0 ] ], 
  [ [ 0, 0 ], [ 1, 0 ] ], [ [ 1, 2 ], [ 0, 1 ] ], 
  [ [ 0, 0, 0 ], [ 0, 0, 0 ] ] ]

#
# Recalling the dimension of  N.
gap> Dimension(N);
7
gap> DimensionVector(N);
[ 3, 2, 2 ]
gap> LoewyLength(N);
4
gap> IsZero(N);
false

#
# Finding the radical, top, socle of the module, projective cover.
gap> RadicalOfModule(N);
<[ 0, 2, 2 ]>
gap> RadicalOfModuleInclusion(N);
<<[ 0, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> TopOfModule(N);
<[ 3, 0, 0 ]>
gap> TopOfModuleProjection(N);   
<<[ 3, 2, 2 ]> ---> <[ 3, 0, 0 ]>>
gap> SocleOfModule(N);
<[ 1, 0, 2 ]>
gap> SocleOfModuleInclusion(N);  
<<[ 1, 0, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> p := ProjectiveCover(N);
<<[ 3, 12, 9 ]> ---> <[ 3, 2, 2 ]>>
gap> q := InjectiveEnvelope(N);             
<<[ 3, 2, 2 ]> ---> <[ 7, 4, 5 ]>>

# Constructing the direct sum of  N  with itself and finding naturally
# associated inclusions and projections.  
gap> M := DirectSumOfQPAModules([N,N]);
<[ 6, 4, 4 ]>
gap> IsDirectSumOfModules(M);
true
gap> DirectSumInclusions(M);
[ <<[ 3, 2, 2 ]> ---> <[ 6, 4, 4 ]>>, <<[ 3, 2, 2 ]> ---> <[ 6, 4, 4 ]>> ]
gap> DirectSumProjections(M);
[ <<[ 6, 4, 4 ]> ---> <[ 3, 2, 2 ]>>, <<[ 6, 4, 4 ]> ---> <[ 3, 2, 2 ]>> ]
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
gap> ZeroModule(A);
<[ 0, 0, 0 ]>
gap> Amod := RightAlgebraModule(A,\*,A);
<17-dimensional right-module over <Rationals[<quiver with 3 vertices and 
5 arrows>]/
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (5 generators)>>>
gap> RightAlgebraModuleToPathAlgebraMatModule(Amod);
<[ 2, 8, 7 ]>

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
gap> ProjDimension(N);
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
gap> gens[1]^A.a;   
[ [ 0, 0, 0 ], [ 1, 2 ], [ 0, 0 ] ]
gap> SubRepresentation(N,[gens[1]]);
<[ 1, 2, 2 ]>
gap> SupportModuleElement(gens[1]);
[ [(1)*v1] ]
gap> f1 := SubRepresentationInclusion(N,[gens[1]]);
<<[ 1, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> f2 := SubRepresentationInclusion(N,[gens[2]]);
<<[ 1, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> SumOfSubmodules(f1,f2);
[ <<[ 2, 2, 2 ]> ---> <[ 3, 2, 2 ]>>, <<[ 1, 2, 2 ]> ---> <[ 2, 2, 2 ]>>, 
  <<[ 1, 2, 2 ]> ---> <[ 2, 2, 2 ]>> ]
gap> IntersectionOfSubmodules([f1,f2]);
<<[ 0, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> RadicalSeries(N);
[ [ 3, 0, 0 ], [ 0, 1, 0 ], [ 0, 1, 1 ], [ 0, 0, 1 ] ]
gap> SocleSeries(N);
[ [ 1, 0, 0 ], [ 1, 1, 0 ], [ 0, 1, 0 ], [ 1, 0, 2 ] ]
gap> CommonDirectSummand(N,N);
[ <[ 1, 0, 0 ]>, <[ 2, 2, 2 ]>, <[ 1, 0, 0 ]>, <[ 2, 2, 2 ]> ]
gap> MaximalCommonDirectSummand(N,S1);
[ [ <[ 1, 0, 0 ]> ], <[ 2, 2, 2 ]>, <[ 0, 0, 0 ]> ]
gap> NumberOfNonIsoDirSummands(N);
[ 2, [ 1, 1 ] ]
gap> IndecProjectiveModules(A);
[ <[ 1, 4, 3 ]>, <[ 0, 2, 2 ]>, <[ 1, 2, 2 ]> ]
gap> IndecInjectiveModules(A);
[ <[ 1, 0, 1 ]>, <[ 4, 2, 2 ]>, <[ 3, 2, 2 ]> ]
gap> DualOfModule(N);
<[ 3, 2, 2 ]>
gap> DTr(N);
<[ 17, 10, 7 ]>
gap> TrD(N);
<[ 0, 2, 1 ]>
gap> NakayamaFunctorOfModule(N);
<[ 0, 0, 0 ]>
gap> StarOfModule(N);
<[ 0, 0, 0 ]>
gap> TransposeOfModule(N);
<[ 17, 10, 7 ]>
gap> TensorProductOfModules(N,DualOfModule(N));                       
[ ( Rationals^5 ), function( m, n ) ... end ]

#
# Testing various operations for homomorphisms of modules
gap> hom := HomOverAlgebra(N,N);
[ <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>, <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>, 
  <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>, <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>, 
  <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>> ]
gap> f := hom[1];
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> IsPathAlgebraMatModuleHomomorphism(f);
true
gap> MatricesOfPathAlgebraMatModuleHomomorphism(f);
[ [ [ -1, -1, 1 ], [ 0, 0, 0 ], [ 0, 0, 0 ] ], [ [ 0, 0 ], [ 0, 0 ] ], 
  [ [ 0, 0 ], [ 0, 0 ] ] ]
gap> PathAlgebraOfMatModuleMap(f);
<Rationals[<quiver with 3 vertices and 5 arrows>]/
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (5 generators)>>
gap> Range(f);
<[ 3, 2, 2 ]>
gap> Source(f);
<[ 3, 2, 2 ]>
gap> Zero(f);
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> Zero(f) = ZeroMapping(N,N);                                      
true
gap> g := hom[2];
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> f = g;
false
gap> f + g;
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> f * g;
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> IdentityMapping(N);
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> m := ImageElm(f, BasisVectors(Basis(N))[1]);
[ [ -1, -1, 1 ], [ 0, 0 ], [ 0, 0 ] ]
gap> ImagesSet(hom[4], BasisVectors(Basis(N)){[1..3]});
[ [ [ 0, -4/3, 0 ], [ 0, 0 ], [ 0, 0 ] ], [ [ 3, -4, 0 ], [ 0, 0 ], [ 0, 0 ] ]
    , [ [ 3, -16/3, 0 ], [ 0, 0 ], [ 0, 0 ] ] ]
gap> PreImagesRepresentative(f, m);                        
[ [ 1, 0, 0 ], [ 0, 0 ], [ 0, 0 ] ]
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
gap> RightInverseOfHomomorphism(f);
false
gap> CoKernel(f);
<[ 2, 2, 2 ]>
gap> CoKernelProjection(f);
<<[ 3, 2, 2 ]> ---> <[ 2, 2, 2 ]>>
gap> EndN := EndOverAlgebra(N);
<algebra-with-one of dimension 5 over Rationals>
gap> BEndN := BasisVectors(Basis(EndN));
[ [ [ 1, 1, -1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 0, 0, 0 ] ], 
  [ [ 0, 0, 0, 0, 0, 0, 0 ], [ 1, 1, -1, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 0, 0, 0 ] ], 
  [ [ 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], 
      [ 1, 1, -1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], 
      [ 0, 0, 0, 0, 0, 0, 0 ] ], 
  [ [ 0, 1, 0, 0, 0, 0, 0 ], [ 0, 21/4, -9/4, 0, 0, 0, 0 ], 
      [ 0, 25/4, -9/4, 0, 0, 0, 0 ], [ 0, 0, 0, 3/2, 0, 0, 0 ], 
      [ 0, 0, 0, -3/4, 3/2, 0, 0 ], [ 0, 0, 0, 0, 0, 3, 3 ], 
      [ 0, 0, 0, 0, 0, -3/4, 0 ] ], 
  [ [ 0, 0, 1, 0, 0, 0, 0 ], [ 0, 25/4, -9/4, 0, 0, 0, 0 ], 
      [ 0, 25/4, -5/4, 0, 0, 0, 0 ], [ 0, 0, 0, 5/2, 0, 0, 0 ], 
      [ 0, 0, 0, -3/4, 5/2, 0, 0 ], [ 0, 0, 0, 0, 0, 4, 3 ], 
      [ 0, 0, 0, 0, 0, -3/4, 1 ] ] ]
gap> BEndN[1] in EndOverAlgebra(N);
true
gap> FromEndMToHomMM(N,BEndN[1]);
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> FromHomMMToEndM(f);        
[ [ -1, -1, 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ],
  [ 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0 ] ]
gap> HomFactoringThroughProjOverAlgebra(N,N);                         
[  ]
gap> HomFromProjective(BasisVectors(Basis(N))[1],N);                  
<<[ 1, 4, 3 ]> ---> <[ 3, 2, 2 ]>>
gap> EndModuloProjOverAlgebra(N);
IdentityMapping( <algebra-with-one of dimension 5 over Rationals> )
gap> Image(f);
<[ 1, 0, 0 ]>
gap> ImageInclusion(f);
<<[ 1, 0, 0 ]> ---> <[ 3, 2, 2 ]>>
gap> ImageProjection(f);
<<[ 3, 2, 2 ]> ---> <[ 1, 0, 0 ]>>
gap> ImageProjectionInclusion(f);
[ <<[ 3, 2, 2 ]> ---> <[ 1, 0, 0 ]>>, <<[ 1, 0, 0 ]> ---> <[ 3, 2, 2 ]>> ]
gap> Kernel(f);
<[ 2, 2, 2 ]>
gap> KernelInclusion(f);
<<[ 2, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> LeftMinimalVersion(f);
[ <<[ 3, 2, 2 ]> ---> <[ 1, 0, 0 ]>>, [ <[ 2, 2, 2 ]> ] ]
gap> RightMinimalVersion(f);
[ <<[ 1, 0, 0 ]> ---> <[ 3, 2, 2 ]>>, [ <[ 2, 2, 2 ]> ] ]
gap> DualOfModuleHomomorphism(f);
<<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> NakayamaFunctorOfModuleHomomorphism(f);
<<[ 0, 0, 0 ]> ---> <[ 0, 0, 0 ]>>
gap> StarOfModuleHomomorphism(f);
<<[ 0, 0, 0 ]> ---> <[ 0, 0, 0 ]>>

#
# Testing homological algebra operations.
gap> CotiltingModule(N,2);
false
gap> ExtOverAlgebra(N,N);
[ <<[ 0, 10, 7 ]> ---> <[ 3, 12, 9 ]>>, [  ], [  ] ]
gap> ExtAlgebraGenerators(N,3);
[ [ 5, 0, 0, 3 ], [ 4, 0, 0, 1 ], 
  [ [ <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>, <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>, 
          <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>>, 
          <<[ 3, 2, 2 ]> ---> <[ 3, 2, 2 ]>> ], [  ], [  ], 
      [ <<[ 3, 6, 3 ]> ---> <[ 3, 2, 2 ]>> ] ] ]
gap> GlobalDimensionOfAlgebra(A,4);
infinity
gap> GlobalDimension(A);
infinity
gap> GorensteinDimensionOfAlgebra(A,5);
4
gap> GorensteinDimension(A);
4
gap> InjDimensionOfModule(S1,4);  
3
gap> InjDimension(S1);
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
gap> RejectOfModule(N,S1);
<<[ 0, 2, 2 ]> ---> <[ 3, 2, 2 ]>>
gap> TraceOfModule(S1,N);
<<[ 1, 0, 0 ]> ---> <[ 3, 2, 2 ]>>
gap> O4 := NthSyzygy(S1,4);
<[ 0, 2, 2 ]>
gap> IsNthSyzygy(O4,4);
true
gap> ProjDimensionOfModule(S1,4);
4
gap> ProjectiveResolutionOfPathAlgebraModule(S1,4);
Computing f2's ...
Computing f3's ...
Computing f4's ...
ProjectiveResolutionOfPathAlgebraModule
gap> PullBack(hom[1],hom[2]);
[ <<[ 5, 4, 4 ]> ---> <[ 3, 2, 2 ]>>, <<[ 5, 4, 4 ]> ---> <[ 3, 2, 2 ]>> ]
gap> PushOut(hom[1],hom[2]);
[ <<[ 3, 2, 2 ]> ---> <[ 4, 4, 4 ]>>, <<[ 3, 2, 2 ]> ---> <[ 4, 4, 4 ]>> ]
gap> P := DirectSumOfQPAModules(IndecProjectiveModules(A){[1,3]});             
<[ 2, 6, 5 ]>
gap> M := IndecProjectiveModules(A)[2];
<[ 0, 2, 2 ]>
gap> AllComplementsOfAlmostCompleteTiltingModule(P,M);                         
[ [  ], 0 -> 0:(0,2,2) -> -1:(4,16,12) -> -2:(4,14,10) -> 0 ]
gap> DominantDimensionOfAlgebra(A,3);
0
gap> DominantDimensionOfModule(N,3);
0
gap> FaithfulDimension(N);                                            
0
gap> HaveFiniteCoresolutionInAddM(SimpleModules(A)[1],N,2);
0 -> 0:(1,0,0) -> -1:(1,0,0) -> 0
gap> HaveFiniteResolutionInAddM(SimpleModules(A)[1],N,2);  
0 -> 1:(1,0,0) -> 0:(1,0,0) -> 0
gap> IsOmegaPeriodic(SimpleModules(Nak)[1],4);
Computing syzygy number: 1
Computing syzygy number: 2
2

#
# Testing some AR-theory.
gap> AlmostSplitSequence(S1);
[ <<[ 7, 4, 3 ]> ---> <[ 8, 4, 3 ]>>, <<[ 8, 4, 3 ]> ---> <[ 1, 0, 0 ]>> ]
gap> IsTauPeriodic(SimpleModules(Nak)[1],4);
3

# Algebras and modules over finite fields
gap> fNak := NakayamaAlgebra(GF(2),[3,2,1]);
<GF(2)[<quiver with 3 vertices and 2 arrows>]>
gap> fS1 := SimpleModules(fNak)[1];
<[ 1, 0, 0 ]>
gap> PredecessorsOfModule(fS1,4);
[ [ [ <[ 1, 0, 0 ]> ], [ <[ 1, 1, 0 ]> ], [ <[ 0, 1, 0 ]>, <[ 1, 1, 1 ]> ], 
      [ <[ 0, 1, 1 ]> ], [ <[ 0, 0, 1 ]> ] ], 
  [ [ [ 1, 1, [ 1, false ] ] ], [ [ 1, 1, [ 1, 1 ] ], [ 2, 1, [ 1, false ] ] ]
        , [ [ 1, 1, [ 1, 1 ] ], [ 1, 2, [ 1, 1 ] ] ], 
      [ [ 1, 1, [ false, 1 ] ] ] ] ]
gap> m1 := [[1,1],[0,1]]*One(GF(2));
[ [ Z(2)^0, Z(2)^0 ], [ 0*Z(2), Z(2)^0 ] ]
gap> m2 := [[0,1],[1,1]]*One(GF(2));
[ [ 0*Z(2), Z(2)^0 ], [ Z(2)^0, Z(2)^0 ] ]
gap> fmat := [m1,m2];
[ [ [ Z(2)^0, Z(2)^0 ], [ 0*Z(2), Z(2)^0 ] ], 
  [ [ 0*Z(2), Z(2)^0 ], [ Z(2)^0, Z(2)^0 ] ] ]
gap> M := RightModuleOverPathAlgebra(fNak,fmat);
<[ 2, 2, 2 ]>
gap> U := BasicVersionOfModule(M);
<[ 1, 1, 1 ]>
gap> IsIndecomposableModule(U);
true
gap> IsInAdditiveClosure(M,U);
true
gap> IsInjectiveModule(U);
true
gap> IsTauRigidModule(U);
true
gap> BlockSplittingIdempotents(M);
[ <<[ 2, 2, 2 ]> ---> <[ 2, 2, 2 ]>> ]
gap> DecomposeModule(M);
[ <[ 1, 1, 1 ]>, <[ 1, 1, 1 ]> ]
gap> DecomposeModuleWithMultiplicities(M);
[ [ <[ 1, 1, 1 ]> ], [ 2 ] ]

#
# Groebner basis commands
#
gap> I := IdealOfQuotient(A);
<two-sided ideal in <Rationals[<quiver with 3 vertices and 5 arrows>]>, 
  (5 generators)>
gap> grb := GroebnerBasisOfIdeal(I);
<complete two-sided Groebner basis containing 5 elements>
gap> IsGroebnerBasis(grb);
true
gap> IsTipReducedGroebnerBasis(grb);
true
gap> TipReduceGroebnerBasis(grb);
<complete two-sided Groebner basis containing 5 elements>
gap> IsCompletelyReducedGroebnerBasis(grb);
true
gap> AdmitsFinitelyManyNontips(grb);
true
gap> CompletelyReduce(grb, A.c^4);
[<zero> of ...]
gap> CompletelyReduceGroebnerBasis(grb);
<complete two-sided Groebner basis containing 5 elements>
gap> Nontips(grb);
[ v1, v2, v3, a, b, c, d, e, a*c, a*d, b*c, b*d, c*d, e*b, b*c*d, e*b*c, 
  e*b*c*d ]
gap> NontipSize(grb);
17
gap> TipReduce(grb,A.c^2 + A.a);
[(1)*a]
gap> rgrb := RightGroebnerBasis(I);
<complete right Groebner basis containing 13 elements>
gap> IsRightGroebnerBasis(rgrb);
true

#
# Vertex projective modules
#
gap> F := GF(11);
GF(11)
gap> Q := Quiver(["v","w", "x"],[["v","w","a"],["v","w","b"],
> ["w","x","c"]]);
<quiver with 3 vertices and 3 arrows>
gap> AA := PathAlgebra(F,Q);
<GF(11)[<quiver with 3 vertices and 3 arrows>]>
gap> P := RightProjectiveModule(AA,[AA.v,AA.v,AA.w]);
<right-module over <GF(11)[<quiver with 3 vertices and 3 arrows>]>>
gap> IsPathAlgebraModule(P);
true
gap> Dimension(P);
12
gap> rgrbP := RightGroebnerBasisOfModule(P);
<Groebner basis of <12-dimensional right-module over <GF(11)[<quiver with 
3 vertices and 3 arrows>]>>, [ [ (Z(11)^0)*v, <zero> of ..., <zero> of ... ], 
  [ <zero> of ..., (Z(11)^0)*v, <zero> of ... ], 
  [ <zero> of ..., <zero> of ..., (Z(11)^0)*w ] ] >
gap> CompletelyReduceGroebnerBasisForModule(rgrbP);
[ [ (Z(11)^0)*v, <zero> of ..., <zero> of ... ], 
  [ <zero> of ..., (Z(11)^0)*v, <zero> of ... ], 
  [ <zero> of ..., <zero> of ..., (Z(11)^0)*w ] ]
gap> p1 := Vectorize(P,[AA.b*AA.c,AA.a*AA.c,AA.c]);
[ (Z(11)^0)*b*c, (Z(11)^0)*a*c, (Z(11)^0)*c ]
gap> p2 := Vectorize(P,[AA.a,AA.b,AA.w]);
[ (Z(11)^0)*a, (Z(11)^0)*b, (Z(11)^0)*w ]
gap> IsPathAlgebraVector(p1);               
false
gap> IsPathAlgebraVector(ExtRepOfObj(p1));
true
gap> LeadingCoefficient(ExtRepOfObj(p1));
Z(11)^0
gap> LeadingComponent(ExtRepOfObj(p1));  
(Z(11)^0)*b*c
gap> LeadingPosition(ExtRepOfObj(p1)); 
1
gap> LeadingTerm(ExtRepOfObj(p1));    
[ (Z(11)^0)*b*c, <zero> of ..., <zero> of ... ]
gap> TargetVertex(ExtRepOfObj(p1));   
x
gap> 2*p1 + p2;
[ (Z(11)^0)*a+(Z(11))*b*c, (Z(11)^0)*b+(Z(11))*a*c, (Z(11)^0)*w+(Z(11))*c ]
gap> UniformGeneratorsOfModule(P);
[ [ <zero> of ..., (Z(11)^0)*v, <zero> of ... ], 
  [ (Z(11)^0)*v, <zero> of ..., <zero> of ... ], 
  [ <zero> of ..., <zero> of ..., (Z(11)^0)*w ] ]
gap> S := SubAlgebraModule(P,[p1,p2]);
<right-module over <GF(11)[<quiver with 3 vertices and 3 arrows>]>>
gap> Dimension(S);
3
gap> p2^(AA.c - AA.w);
[ (Z(11)^5)*a+(Z(11)^0)*a*c, (Z(11)^5)*b+(Z(11)^0)*b*c, 
  (Z(11)^5)*w+(Z(11)^0)*c ]
gap> p1 < p2;
false
gap> PS := P/S;
<9-dimensional right-module over <GF(11)[<quiver with 3 vertices and 
3 arrows>]>>
gap> ProjectivePathAlgebraPresentation(N);
[ <right-module over <Rationals[<quiver with 3 vertices and 5 arrows>]>>, 
  <right-module over <Rationals[<quiver with 3 vertices and 5 arrows>]>>, 
  [ [ (1)*v1, <zero> of ..., <zero> of ... ], 
      [ <zero> of ..., (1)*v1, <zero> of ... ], 
      [ <zero> of ..., <zero> of ..., (1)*v1 ] ], 
  [ [ (-1)*a, (-1)*a, (1)*a ], [ (3/4)*b, (5/2)*a, (-3/2)*a ], 
      [ <zero> of ..., (5)*a+(1)*b, (-3)*a ], 
      [ <zero> of ..., (-50/3)*a+(-10/3)*a*c, (10)*a ], 
      [ <zero> of ..., (-10/3)*a*d*e, <zero> of ... ], 
      [ <zero> of ..., (25/3)*a, (-5)*a+(1)*b ], 
      [ <zero> of ..., (10/3)*a, (-2)*a+(2/5)*a*c ], 
      [ <zero> of ..., <zero> of ..., (2/5)*a*d*e ] ], 
  [ [ (-1)*a, (3/4)*b, <zero> of ..., <zero> of ..., <zero> of ..., 
          <zero> of ..., <zero> of ..., <zero> of ... ], 
      [ (-1)*a, (5/2)*a, (5)*a+(1)*b, (-50/3)*a+(-10/3)*a*c, (-10/3)*a*d*e, 
          (25/3)*a, (10/3)*a, <zero> of ... ], 
      [ (1)*a, (-3/2)*a, (-3)*a, (10)*a, <zero> of ..., (-5)*a+(1)*b, 
          (-2)*a+(2/5)*a*c, (2/5)*a*d*e ] ],
  [ [ (-1)*a, (-1)*a, (1)*a ], [ (3/4)*b, (5/2)*a, (-3/2)*a ], 
      [ <zero> of ..., (5)*a+(1)*b, (-3)*a ], 
      [ <zero> of ..., (-50/3)*a+(-10/3)*a*c, (10)*a ], 
      [ <zero> of ..., (-10/3)*a*d*e, <zero> of ... ], 
      [ <zero> of ..., (25/3)*a, (-5)*a+(1)*b ], 
      [ <zero> of ..., (10/3)*a, (-2)*a+(2/5)*a*c ], 
      [ <zero> of ..., <zero> of ..., (2/5)*a*d*e ] ] ]

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
infinity
gap> LowestKnownDegree(C);
-infinity
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
gap> IsUnitForm(u);
true
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
gap> EulerBilinearFormOfAlgebra(NakayamaAlgebra(Rationals,[3,2,1]));
function( x, y ) ... end
gap> UnitForm([[2,1],[1,2]]);
[ [ 2, 1 ], [ 1, 2 ] ]


# 
# Testing ReadAlgebra/SaveAlgebra
gap> A := PathAlgebra(Rationals, DynkinQuiver("A", 11, ["r","r","r","r","r","r","r","r","r","r"]));
<Rationals[<quiver with 11 vertices and 10 arrows>]>
gap> SaveAlgebra(A, Filename(DirectoryCurrent(), "writealgebratest~"), "delete");
true
gap> ReadAlgebra("writealgebratest~");
<Rationals[<quiver with 11 vertices and 10 arrows>]>
gap> RemoveFile("writealgebratest~");
true

#
# Testing iterators/enumerators of path algebras
gap> e := Enumerator(Quiver(5, [[1,3],[2,3],[3, 4],[3,5]]));
<object>
gap> [e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], e[10], e[11], e[12], e[13], e[14]];
[ v1, v2, v3, v4, v5, a1, a2, a3, a4, a1*a3, a1*a4, a2*a3, a2*a4, fail ]
gap> e := Enumerator(Quiver(1, [[1,1],[1,1]]));
<object>
gap> [e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], e[10], e[11], e[12], e[13], e[14], e[15], e[16], e[17], e[18], e[19], e[20] ];
[ v1, a1, a2, a1^2, a1*a2, a2*a1, a2^2, a1^3, a1^2*a2, a1*a2*a1, a1*a2^2, a2*a1^2, a2*a1*a2, a2^2*a1, a2^3, a1^4, a1^3*a2, a1^2*a2*a1, a1^2*a2^2, a1*a2*a1^2 ]





#
# Testing isolated vertices
gap> Q := Quiver(3, [ [1, 2, "a"] ]);
<quiver with 3 vertices and 1 arrows>
gap> K := Rationals;
Rationals
gap> KQ := PathAlgebra(K, Q);
<Rationals[<quiver with 3 vertices and 1 arrows>]>

# Testing kernels and cokernels for isolated vertices:
gap> module := RightModuleOverPathAlgebra(KQ, [1, 1, 2], [["a", [[1]]]]);
<[ 1, 1, 2 ]>
gap> homomorphism := RightModuleHomOverAlgebra(module, module, [ [[1]], [[1]], [ [1, 0], [0, 0] ]  ]);
<<[ 1, 1, 2 ]> ---> <[ 1, 1, 2 ]>>
gap> CoKernel(homomorphism);
<[0, 0, 1]>
gap> Kernel(homomorphism);
<[0, 0, 1]>

# Testing submodules for isolated vertices:
gap> S := RightModuleOverPathAlgebra(KQ, [0, 0, 1], []);
<[ 0, 0, 1 ]>
gap> x := BasisVectors(Basis(S))[1];
[ [ 0 ], [ 0 ], [ 1 ] ]
gap> Source(SubRepresentationInclusion(S, [x]));
<[0, 0, 1]>

# Testing restrictions via algebra homomorphisms for isolated vertices:
gap> RestrictionViaAlgebraHomomorphism(IdentityMapping(KQ), module);
<[1, 1, 2]>

# Testing RightAlgebraModuleToPathAlgebraMatModule for isolated vertices:
gap> RightAlgebraModuleToPathAlgebraMatModule(module);
<[ 1, 1, 2 ]>


gap> STOP_TEST("qpa.tst");
