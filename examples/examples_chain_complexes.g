# alg should be defined!
alg;
# L, M, and N are alg-modules
# f: L --> M and g: M --> N are non-zero morphisms
cat := CatOfRightAlgebraModules(alg);
cat.zeroObj;
cat.isZeroObj(M);
cat.zeroMap(M,N);
cat.composeMaps(g,f);
cat.ker(g);
cat.isExact(g,f);
#################
# L, M and N are modules over the same algebra A
# cat is the category mod A
# f: L --> M and g: M --> N maps
C := FiniteComplex(cat, 1, [g,f]);
##################
Ms := StalkComplex(cat, M, 3);
##################
ses := ShortExactSequence(cat, f, g);
##################
C;                     
# Want to compute the homology in degree 2
f := DifferentialOfComplex(C,3);
g := KernelInclusion(DifferentialOfComplex(C,2));
# We know that Im f is included in Ker g, so can find the
# lifting morphism h from C_3 to Ker g.
h := LiftingInclusionMorphisms(g,f);
# The cokernel of h is Ker g / Im f 
Homology := CoKernel(h);
###################
C;
IsFiniteComplex(C);
UpperBound(C);
LowerBound(C);
LengthOfComplex(C);
HighestKnownDegree(C);
LowestKnownDegree(C);
###################
C;                
IsExactSequence(C);
IsExactInDegree(C,1);
IsExactInDegree(C,2);
###################
C;
Shift(C,1);
D := Shift(C,-1);
dc := DifferentialOfComplex(C,3);
dd := DifferentialOfComplex(D,4);
MatricesOfPathAlgebraMatModuleHomomorphism(dc);
MatricesOfPathAlgebraMatModuleHomomorphism(dd);
####################
C2;
C3;
YonedaProduct(C2,C3);
####################

