# GAP Implementation
# This file was generated from 
# $Id: modright.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallMethod( RightModuleByGenerators,
  "for ring and collection",
  true,
  [ IsRing, IsCollection ], 0,
  function( R, gens )
    local V;

    V := Objectify( NewType( FamilyObj( gens ),
                             IsRightModule and IsAttributeStoringRep ),
                      rec() );
    SetRightActingDomain( V, R );
    SetGeneratorsOfRightModule( V, AsList( gens ) );

    return V;

  end
);


InstallOtherMethod( RightModuleByGenerators,
  "for ring, homogeneous list, and vector",
  true,
  [ IsRing, IsHomogeneousList, IsVector ], 0,
  function ( R, gens, zero )
    local V;

    V := Objectify( NewType( CollectionsFamily( FamilyObj( zero ) ),
                    IsRightModule and IsAttributeStoringRep ),
                    rec() );
    SetRightActingDomain( V, R );
    SetGeneratorsOfRightModule( V, AsList( gens ) );
    SetZero( V, zero );

    if IsEmpty( gens ) then
      SetIsTrivial( V, true );
    fi;

    return V;

  end
);
