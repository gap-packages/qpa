Q := Quiver( ["u","v"] , [ ["u","u","a"], ["u","v","b"], 
["v","u","c"], ["v","v","d"] ] );
F := Rationals;
FQ := PathAlgebra(F,Q);
##########
IsPathAlgebra(FQ);
##########
QuiverOfPathAlgebra(FQ);
##########
FQ.a;
FQ.v;
elem := 2*FQ.a - 3*FQ.v;
##########
IsLeftUniform(elem);
IsRightUniform(elem);
IsUniform(elem);
another := FQ.a*FQ.b + FQ.b*FQ.d*FQ.c*FQ.b*FQ.d;
IsLeftUniform(another);
IsRightUniform(another);
IsUniform(another);
##########
elem := FQ.a*FQ.b*FQ.c + FQ.b*FQ.d*FQ.c+FQ.d*FQ.d;
LeadingTerm(elem);
LeadingCoefficient(elem);
mon := LeadingMonomial(elem);
mon in FQ;
mon in Q;
##########
Q := Quiver( 1, [ [1,1,"a"], [1,1,"b"] ] );
kQ := PathAlgebra(Rationals, Q);
gens := GeneratorsOfAlgebra(kQ);
a := gens[2];
b := gens[3];
relations := [a^2,a*b-b*a, b*b];
A := kQ/relations;
IndecProjectiveModules(A);  
##########
gb := GBNPGroebnerBasis(relations,kQ);
I := Ideal(kQ,gb);
GroebnerBasis(I,gb);
IndecProjectiveModules(A);
A := kQ/I;
IndecProjectiveModules(A);
##########
I := Ideal(FQ, [FQ.a - FQ.b*FQ.c, FQ.d*FQ.d]); 
GeneratorsOfIdeal(I);
IsIdealInPathAlgebra(I);
##########
rels := [FQ.a - FQ.b*FQ.c, FQ.d*FQ.d];
gb := GBNPGroebnerBasis(rels, FQ); 
I := Ideal(FQ, gb);
GroebnerBasis(I, gb);
quot := FQ/I;
##########
quot := FQ/I;
IsQuotientOfPathAlgebra(quot);
IsQuotientOfPathAlgebra(FQ);
##########
Q := Quiver(5, [ [1,2,"a"], [2,4,"b"], [3,2,"c"], [2,5,"d"] ]);
A := PathAlgebra(Rationals, Q);
IsFiniteTypeAlgebra(A);
quo := A/[A.a*A.b, A.c*A.d];;
IsFiniteTypeAlgebra(quo);
##########
elem := quot.a*quot.b;
IsElementOfQuotientOfPathAlgebra(elem);
IsElementOfQuotientOfPathAlgebra(FQ.a*FQ.b); 
##########
alg := NakayamaAlgebra([2,1], Rationals);
QuiverOfPathAlgebra(alg);
##########
Q := Quiver(1, [ [1,1,"a"], [1,1,"b"] ]);;  
A := PathAlgebra(Rationals, Q);;
IsSpecialBiserialAlgebra(A); IsStringAlgebra(A);
rel1 := [A.a*A.b, A.a^2, A.b^2];  
quo1 := A/rel1;;
IsSpecialBiserialAlgebra(quo1); IsStringAlgebra(quo1);
rel2 := [A.a*A.b-A.b*A.a, A.a^2, A.b^2];  
quo2 := A/rel2;;
IsSpecialBiserialAlgebra(quo2); IsStringAlgebra(quo2);
rel3 := [A.a*A.b+A.b*A.a, A.a^2, A.b^2, A.b*A.a];   
quo3 := A/rel3;;
IsSpecialBiserialAlgebra(quo3); IsStringAlgebra(quo3);
rel4 := [A.a*A.b, A.a^2, A.b^3]; 
quo4 := A/rel4;;
IsSpecialBiserialAlgebra(quo4); IsStringAlgebra(quo4);
##########
Q1 := Quiver(4, [[1,2,"a"], [1,2,"b"], [2,3,"c"], [3,4,"d"], [3,4,"e"]]);
k := Rationals;
kQ1 := PathAlgebra(k,Q1);
rho := [kQ1.b*kQ1.c, kQ1.c*kQ1.d];
A1 := kQ1/rho;
IsValidString(A1,"eca");
StringsLessThan(A1,2);
IsABand(A1,"eca");
IsABand(A1,"Ab");
BandsLessThan(A1,3);
BandRepresentativesLessThan(A1,3);
IsDomesticStringAlgebra(A1);
Q2 := BridgeQuiver(A1);
Display(Q2);
LocalARQuiver(A1,"eca");
##########
Q := Quiver( [ "u", "v" ], [ [ "u", "u", "a" ], 
             [ "u", "v", "b" ] ] );
Qop := OppositeQuiver(Q);
VerticesOfQuiver( Qop );
ArrowsOfQuiver( Qop );
OppositePath( Q.a * Q.b );
IsIdenticalObj( Q, OppositeQuiver( Qop ) );
OppositePath( Qop.b_op * Qop.a_op );
##########
Q := Quiver( [ "u", "v" ], [ [ "u", "u", "a" ], 
             [ "u", "v", "b" ] ] );
A := PathAlgebra( Rationals, Q );
OppositePathAlgebra( A );
OppositePathAlgebraElement( A.u + 2*A.a + 5*A.a*A.b );
IsIdenticalObj( A, 
        OppositePathAlgebra( OppositePathAlgebra( A ) ) );
##########
q1 := Quiver( [ "u1", "u2" ], [ [ "u1", "u2", "a" ] ] );
q2 := Quiver( [ "v1", "v2", "v3" ],
                      [ [ "v1", "v2", "b" ],
                        [ "v2", "v3", "c" ] ] );
q1_q2 := QuiverProduct( q1, q2 );
q1_q2.u1_b * q1_q2.a_v2;
IncludeInProductQuiver( [ q1.a, q2.b * q2.c ], q1_q2 );
ProjectFromProductQuiver( 2, q1_q2.a_v1 * q1_q2.u2_b * q1_q2.u2_c );
q1_q2_dec := QuiverProductDecomposition( q1_q2 );
q1_q2_dec[ 1 ];
q1_q2_dec[ 1 ] = q1;
##########
q1 := Quiver( [ "u1", "u2" ], [ [ "u1", "u2", "a" ] ] );
q2 := Quiver( [ "v1", "v2", "v3", "v4" ],
                      [ [ "v1", "v2", "b" ],
                        [ "v1", "v3", "c" ],
                        [ "v2", "v4", "d" ],
                        [ "v3", "v4", "e" ] ] );
fq1 := PathAlgebra( Rationals, q1 );
fq2 := PathAlgebra( Rationals, q2 );
rels := [ fq2.b * fq2.d - fq2.c * fq2.e ];
quot := fq2 / rels;
t := TensorProductOfAlgebras( fq1, quot );
SimpleTensor( [ fq1.a, quot.b ], t );
t_dec := TensorProductDecomposition( t );
t_dec[ 1 ] = fq1;
