Q := Quiver( 3, [ [1,2,"a"], [2,3,"b"] ] );
PA := PathAlgebra( Rationals, Q );
rels := [ PA.a*PA.b ];
gb := GBNPGroebnerBasis( rels, PA );
I := Ideal( PA, gb );
grb := GroebnerBasis( I, gb );
alg := PA/I;
IsGroebnerBasis(gb);
IsGroebnerBasis(grb);
####################