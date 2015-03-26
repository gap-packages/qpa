LoadPackage("QPA");
Q := Quiver( 2, [[ 1, 2, "a" ] ] );
kQ := PathAlgebra( Rationals, Q );
M := RightModuleOverPathAlgebra( kQ, [ [ "a", [ [ 1, 1 ] ] ] ] );
N := RightModuleOverPathAlgebra( kQ, [ [ "a", [ [ 1, 0 ], [ 0, 1 ] ] ] ] );
#O := RightModuleOverPathAlgebra( kQ, [ [ "a", [ [ 1, 1 ] ] ] ] );
f := RightModuleHomOverAlgebra( M, N, [ [ [ 1, 0 ] ], [ [ 1, 0 ], [ 0, 0 ] ] ] );
#g := RightModuleHomOverAlgebra( N, O, [ [ [ 1 ], [ 1 ] ], [ [ 0, 1 ], [ 1, 0 ] ] ] );

S1 := SimpleModules(kQ)[1];
S2 := SimpleModules(kQ)[2];
g := RightModuleHomOverAlgebra( N, S1, [  [ [0],[1] ],  [ [0],[0] ]  ] );
h := RightModuleHomOverAlgebra( M, S1, [  [ [1] ],  [ [0],[0] ]  ] );
