Q := Quiver( 1, [ [1,1,"a"], [1,1,"b"] ] );
Display(Q);
##########
VerticesOfQuiver(Q);
ArrowsOfQuiver(Q);
Q.a;
##########
kQ := PathAlgebra(Rationals, Q);
Display(kQ);
##########
gens := GeneratorsOfAlgebra(kQ);
v1 := gens[1];
a := gens[2];
b := gens[3];
id := One(kQ);
v1 = id;
##########
relations := [a^2,a*b-b*a, b*b];
A := kQ/relations;
##########
Q := Quiver( 3, [ [1,2,"a"], [2,3,"b"], [3,3,"c"] ]);
kQ := PathAlgebra(Rationals, Q);
relations := [kQ.c*kQ.c];
A := kQ/relations;
##########
projectives := IndecProjectiveModules(A);
proj1 := projectives[1];
Display(proj1);
##########
M := MatricesOfPathAlgebraModule(proj1);
M[1];
##########
injectives := IndecInjectiveModules(A);
simples := SimpleModules(A);
##########
s1 := simples[1];  
inj1 := injectives[1];
s1 = inj1;
IsomorphicModules(s1,inj1);
##########
M := RightModuleOverPathAlgebra( A, [0,1,1], [ ["b", [[1]] ] ] );
##########
N := RightModuleOverPathAlgebra( A, [1,2,2], [ ["a",[[1,1]] ], 
             ["b", [[1,0], [-1,0]] ], ["c", [[0,0],[1,0]] ] ] );
##########
f := RightModuleHomOverAlgebra(M,N, [ [[0]], [[1,1]], NullMat(1,2,Rationals)]);
Display(f);
##########
MatricesOfPathAlgebraMatModuleHomomorphism(f);
