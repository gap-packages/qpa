LoadPackage( "QPA" );
Q := Quiver( 2, [ [1, 2, "a"] ] );
kQ := PathAlgebra( Rationals, Q );
M1 := RightModuleOverPathAlgebra( kQ, [ 1, 1 ], [ ["a", [[1]]] ] );
M2 := RightModuleOverPathAlgebra( kQ, [ 3, 3 ], [ ["a", [[1,0,0],[1,0,0],[0,1,0]]] ] );
M3 := RightModuleOverPathAlgebra( kQ, [ 1, 2 ], [ ["a", [[0,0]]] ] );
f := RightModuleHomOverAlgebra( M1, M2,
                                [ [[1,0,0]],
                                  [[1,0,0]] ] );
g := RightModuleHomOverAlgebra( M2, M3,
                                [ [[0],[0],[1]],
                                  [[0,0],[0,0],[1,0]] ] );
C := ComplexFromMorphismList( [ g, f ] );
homology := HomologyFunctor( RightModuleCategory( kQ ), 1 );
h := ApplyFunctor( homology, C );
