# GAP Implementation
# This file was generated from
# $Id: alggrp.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallGlobalFunction( GroupAlgebra, function( arg )
    if Length( arg ) = 2 and IsField( arg[1] ) and IsGroup( arg[2] ) then
        return GroupRing( arg[1], arg[2] );
    else
        Error("usage: GroupAlgebra( <F>, <G> )");
    fi;
end );
InstallMethod( PowerSubalgebraOp,
    "for group algebras",
    true,
    [IsGroupAlgebra, IsPosInt], 0,
    function( A, n )
        local F, p, m, powerMap, G, gens, permGroup, autoGroup, orbits, fam,
              one, zero, basis, q;

        F := LeftActingDomain(A);
        p := Characteristic(F);
        m := DegreeOverPrimeField(F);
        G := GroupOfGroupAlgebra(A);

        if (m mod n) <> 0 then
            Error("<d> must be a factor of field degree");
        fi;

        if (Size(G) mod p) = 0 then
            TryNextMethod();
        fi;

        q := p^n;
        powerMap := GroupHomomorphismByFunction(G, G, x -> x^q, x -> x^(-q));
        autoGroup := InnerAutomorphismsAutomorphismGroup(AutomorphismGroup(G));
        if IsList(autoGroup) then
            gens := Concatenation(autoGroup, [powerMap]);
        else
            gens := Concatenation(GeneratorsOfGroup(autoGroup), [powerMap]);
        fi;
        permGroup := Group(gens);
        orbits := Orbits(permGroup, G);

        fam := ElementsFamily(FamilyObj(A));
        one := One(F);
        zero := Zero(F);

        basis := List( orbits, x -> ElementOfMagmaRing( fam, zero, 
                        List(x, x -> one), x) );
        return SubalgebraNC( A, basis, "basis" );
        
    end );
InstallMethod( GroupOfGroupAlgebra,
    "for group algebras",
    true,
    [ IsGroupAlgebra ], 0,
    UnderlyingMagma
);
