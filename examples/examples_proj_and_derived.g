########## alg should be defined!
alg;
cat := CatOfRightAlgebraModules(alg);
C := StalkComplex(cat, IndecInjectiveModules(alg)[1], 0);
ProjC := ProjectiveResolutionOfComplex(C);
InjC := ProjectiveToInjectiveComplex(ProjC);
TauC := TauOfComplex(C);
IsProjectiveComplex(C);
IsInjectiveComplex(C);
C;
##########
##########
##########
##########
##########