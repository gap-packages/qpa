LoadPackage("QPA");
Q := Quiver( 2, [[ 1, 2, "a" ] ] );
kQ := PathAlgebra( Rationals, Q );
M := RightModuleOverPathAlgebra( kQ, [ [ "a", [ [ 1, 1 ] ] ] ] );
N := RightModuleOverPathAlgebra( kQ, [ [ "a", [ [ 1, 0 ], [ 0, 1 ] ] ] ] );
O := RightModuleOverPathAlgebra( kQ, [ [ "a", [ [ 1, 1 ] ] ] ] );
f := RightModuleHomOverAlgebra( M, N, [ [ [ 1, 0 ] ], [ [ 1, 0 ], [ 0, 0 ] ] ] );
g := RightModuleHomOverAlgebra( N, O, [ [ [ 1 ], [ 1 ] ], [ [ 1, 1 ], [ 1, 1 ] ] ] );

S1 := SimpleModules(kQ)[1];
S2 := SimpleModules(kQ)[2];
#g := RightModuleHomOverAlgebra( N, S1, [  [ [0],[1] ],  [ [0],[0] ]  ] );
h := RightModuleHomOverAlgebra( M, S1, [  [ [1] ],  [ [0],[0] ]  ] );

M1 := RightModuleOverPathAlgebra( kQ, [ [ "a", [ [ 1 ] ] ] ] );
M2 := RightModuleOverPathAlgebra( kQ, [ [ "a", [ [ 1 ] ] ] ] );
ff := RightModuleHomOverAlgebra( M1, M2, [ [ [ 1 ] ], [ [ 1 ] ] ] );

D := CapFunctor( "DualOfModule",
                 RightModuleCategory( kQ ),
                 RightModuleCategory( OppositePathAlgebra( kQ ) ) );
AddObjectFunction( D, DualOfModule );
AddMorphismFunction( D, DualOfModuleHomomorphism );

C := ComplexFromMorphismList( [ f, g ], 0 );
