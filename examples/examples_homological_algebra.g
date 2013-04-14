Q := Quiver(3,[[1,2,"a"],[1,2,"b"],[2,2,"c"],[2,3,"d"],
[3,1,"e"]]);;
KQ := PathAlgebra(Rationals, Q);;
gen := GeneratorsOfAlgebra(KQ);;
a := gen[4];;
b := gen[5];;
c := gen[6];;
d := gen[7];;
e := gen[8];;
rels := [d*e,c^2,a*c*d-b*d,e*a];;
A := KQ/rels;;
mat := [["a",[[1,2],[0,3],[1,5]]],["b",[[2,0],[3,0],[5,0]]],
["c",[[0,0],[1,0]]],["d",[[1,2],[0,1]]],["e",[[0,0,0],[0,0,0]]]];;
N := RightModuleOverPathAlgebra(A,mat);;
#####

# N should be defined!
B := BasisVectors(Basis(N));
g := SubRepresentationInclusion(N,[B[1],B[4]]);
f := SubRepresentationInclusion(N,[B[1],B[2]]);
LiftingInclusionMorphisms(f,g);
S := SimpleModules(A); 
homNS := HomOverAlgebra(N,S[1]);
f := homNS[1];
p := ProjectiveCover(S[1]);
LiftingMorphismFromProjective(f,p);
################
hom := HomOverAlgebra(N,N);
g := MorphismOnKernel(hom[1],hom[2],hom[1],hom[2]);
IsomorphicModules(Source(g),Range(g));
p := ProjectiveCover(N);
N1 := Kernel(p);
pullback := PullBack(p,hom[1]);
Kernel(pullback[1]);
IsomorphicModules(N1,Kernel(pullback[1]));
t := LiftingMorphismFromProjective(p,p*hom[1]);
s := MorphismOnKernel(p,p,t,hom[1]);    
Source(s) = N1;
q := KernelInclusion(p);
pushout := PushOut(q,s);
U := CoKernel(pushout[1]);
IsomorphicModules(U,N);
###################
S := SimpleModules(A);
Ext := ExtOverAlgebra(S[2],S[2]);
Length(Ext[2]);
# i.e. Ext^1(S[2],S[2]) is 1-dimensional
pushout := PushOut(Ext[2][1],Ext[1]);   
f := CoKernelProjection(pushout[1]);
U := Range(pushout[1]);            
