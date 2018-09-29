Q := Quiver(3,[[1,2,"a"],[1,2,"b"],[2,2,"c"],[2,3,"d"],[3,1,"e"]]);;
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
###########
L := RightModuleOverPathAlgebra(A,[["a",[0,1]],["b",[0,1]],
["c",[[0]]],["d",[[1]]],["e",[1,0]]]);
DimensionVector(L);
f := RightModuleHomOverAlgebra(L,N,[[[0,0,0]], [[1,0]], 
[[1,2]]]);
IsPathAlgebraMatModuleHomomorphism(f);
###########
B := BasisVectors(Basis(N)); 
PreImagesRepresentative(f,B[4]);     
PreImagesRepresentative(f,B[5]);
BL := BasisVectors(Basis(L));
ImageElm(f,BL[1]);
ImagesSet(f,L);
ImagesSet(f,BL);
###########
z := Zero(f);
f = z;
Range(f) = Range(z);
y := ZeroMapping(L,N); 
y = z;            
id := IdentityMapping(N);  
f*id;
# #This causes an error!
# id*f;
# quit;;
2*f + z;
###########
L := RightModuleOverPathAlgebra(A,[["a",[0,1]],["b",[0,1]],
["c",[[0]]],["d",[[1]]],["e",[1,0]]]);;
f := RightModuleHomOverAlgebra(L,N,[[[0,0,0]], [[1,0]], 
[[1,2]]]);
g := CoKernelProjection(f);
CoKernelOfWhat(g) = f;
h := ImageProjection(f);
ImageOfWhat(h) = f;
IsInjective(f); IsSurjective(f); IsIsomorphism(f); 
IsIsomorphism(h);
###########
S := SimpleModules(A)[1];;
H := HomOverAlgebra(N,S);; 
IsSplitMonomorphism(H[1]);  
IsSplitEpimorphism(H[1]);
IsSplitMonomorphism(f);
###########
L := RightModuleOverPathAlgebra(A,[["a",[0,1]],["b",[0,1]],
["c",[[0]]],["d",[[1]]],["e",[1,0]]]);
f := RightModuleHomOverAlgebra(L,N,[[[0,0,0]], [[1,0]], 
[[1,2]]]);;
IsZero(0*f);
g := KernelInclusion(f);
KnownAttributesOfObject(g);
KernelOfWhat(g) = f;
###########
MatricesOfPathAlgebraMatModuleHomomorphism(f);
Range(f);
Source(f);
Source(f) = L;
###########
hom := HomOverAlgebra(N,N);
g := hom[1];
M := CoKernel(g);
f := CoKernelProjection(g);
Range(f) = M;
endo := EndOverAlgebra(N);
RadicalOfAlgebra(endo);
B := BasisVectors(Basis(N));
p := HomFromProjective(B[1],N);
U := Image(p);
projinc := ImageProjectionInclusion(p);
U = Range(projinc[1]);                                      
Kernel(p);
###########
H:= HomOverAlgebra(N,N);;
RightMinimalVersion(H[1]);   
LeftMinimalVersion(H[1]);             
S := SimpleModules(A)[1];;
MinimalRightApproximation(N,S);
S := SimpleModules(A)[3];;
MinimalLeftApproximation(S,N);    
###########
f := RadicalOfModuleInclusion(N);
radN := Source(f);
g := SocleOfModuleInclusion(N);
U := SubRepresentationInclusion(N,[B[5]+B[6],B[7]]);
h := TopOfModuleProjection(N);