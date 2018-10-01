Q := Quiver(2, [[1, 2, "a"], [2, 1, "b"],[1, 1, "c"]]);
P := PathAlgebra(Rationals, Q);
matrices := [["a", [[1,0,0],[0,1,0]]], 
 ["b", [[0,1],[1,0],[0,1]]],
 ["c", [[0,0],[1,0]]]];
M := RightModuleOverPathAlgebra(P,matrices);
mats := [ [[1,0,0], [0,1,0]], [[0,1],[1,0],[0,1]], 
          [[0,0],[1,0]] ];; 
N := RightModuleOverPathAlgebra(P,mats); 
arrows := ArrowsOfQuiver(Q);
mats := [[arrows[1], [[1,0,0],[0,1,0]]],
         [arrows[2], [[0,1],[1,0],[0,1]]], 
         [arrows[3], [[0,0],[1,0]]]];;
N := RightModuleOverPathAlgebra(P,mats); 
# Next we give the vertex simple associate to vertex 1. 
M := RightModuleOverPathAlgebra(P,[["a",[1,0]],["b",[0,1]],
             ["c",[[0]]]]);
# The zero module. 
M := RightModuleOverPathAlgebra(P,[["a",[0,0]],["b",[0,0]],
             ["c",[0,0]]]);
Dimension(M);
Basis(M);
matrices := [["a", [[1,0,0],[0,1,0]]], ["b",
 [[0,1],[1,0],[0,1]]], ["c", [[0,0],[1,0]]]];
M := RightModuleOverPathAlgebra(P,[2,3],matrices);
M := RightModuleOverPathAlgebra(P,[2,3],[]);  
A := P/[P.c^2 - P.a*P.b, P.a*P.b*P.c, P.b*P.c];         
Dimension(A);
Amod := RightAlgebraModule(A,\*,A);                       
RightAlgebraModuleToPathAlgebraMatModule(Amod);
##################
# Using the path algebra P from the above example. 
matrices := [["a", [[1,0,0],[0,1,0]]], 
["b", [[0,1],[1,0],[0,1]]], ["c", [[0,0],[1,0]]]];
M := RightModuleOverPathAlgebra(P,matrices);
B:=BasisVectors(Basis(M));
B[1] + B[3];
4*B[2];
m := 5*B[1] + 2*B[4]+B[5];
m^(P.a*P.b-P.c);
B[1]^P.a;
B[2]^P.b;
B[4]^(P.b*P.c);
##################
Q  := Quiver(3,[[1,2,"a"],[1,2,"b"],[2,2,"c"],[2,3,"d"],
[3,1,"e"]]);
KQ := PathAlgebra(Rationals, Q);
gens := GeneratorsOfAlgebra(KQ);
u := gens[1];; v := gens[2];;
w := gens[3];; a := gens[4];;
b := gens[5];; c := gens[6];;
d := gens[7];; e := gens[8];;
rels := [d*e,c^2,a*c*d-b*d,e*a];;
A := KQ/rels;
mat := [["a",[[1,2],[0,3],[1,5]]],["b",[[2,0],[3,0],[5,0]]],
["c",[[0,0],[1,0]]],["d",[[1,2],[0,1]]],["e",[[0,0,0],[0,0,0]]]];;
N := RightModuleOverPathAlgebra(A,mat);          
##################
N2 := DirectSumOfQPAModules([N,N]);
proj := DirectSumProjections(N2);
inc := DirectSumInclusions(N2);  
##################
F := GF(11);
Q := Quiver(["v","w", "x"],[["v","w","a"],["v","w","b"],
["w","x","c"]]);
A := PathAlgebra(F,Q);
P := RightProjectiveModule(A,[A.v,A.v,A.w]);
Dimension(P);
##################
p1 := Vectorize(P,[A.b*A.c,A.a*A.c,A.c]);
p2 := Vectorize(P,[A.a,A.b,A.w]);
2*p1 + p2;
S := SubAlgebraModule(P,[p1,p2]);
Dimension(S);
##################
p2^(A.c - A.w);
##################
p1 < p2;
##################
PS := P/S;
Basis(PS);
##################