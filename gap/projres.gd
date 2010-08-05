# Projective Resolutions File
# This file was generated from
# $Id: projres.gd,v 1.2 2010/08/05 17:39:37 uid414323 Exp $
DeclareInfoClass( "InfoProjectiveResolutionFpPathAlgebraModule" );

DeclareCategory("IsProjectiveResolutionFpPathAlgebraModule",IsObject);
DeclareAttribute("ParentAlgebra", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("Maps", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("Module", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("Projectives", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("RProjectives", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("ProjectivesVertexList", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("ProjectivesFList", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("ProjectivesFPrimeList", IsProjectiveResolutionFpPathAlgebraModule);

#############################################################################
##
#R  IsProjectiveResolutionFpPathAlgebraModule
##
##  The following two lines create the representation and family for
##  our Projective Presentation objects.
##
DeclareRepresentation("IsProjectiveResolutionFpPathAlgebraModuleDefaultRep",
  IsAttributeStoringRep and IsPositionalObjectRep
  and IsProjectiveResolutionFpPathAlgebraModule,[]);

BindGlobal("ProjectiveResolutionFpPathAlgebraModuleFamily",
  NewFamily("ProjectiveResolutionFpPathAlgebraModuleFamily",
  IsProjectiveResolutionFpPathAlgebraModule, IsProjectiveResolutionFpPathAlgebraModule) );


DeclareOperation("ProjectiveResolutionFpPathAlgebraModule",
  [ IsAlgebra,
    IsRing and HasGroebnerBasisOfIdeal,
    IsRingElementTable]);

DeclareOperation("ProjectiveResolutionFpPathAlgebraModule",
  [ IsAlgebra,
    IsRing and HasGroebnerBasisOfIdeal,
    IsRingElementTable,
    IsPosInt ]
);

#DeclareOperation("ProjectiveResolutionFpPathAlgebraModule", [IsAlgebra,IsMatrix]);

DeclareOperation("XSetOfPathAlgebraVector",
  [ IsVertexProjectiveModule,
    IsRing and HasGroebnerBasisOfIdeal,
    IsPathAlgebraVector]
);

DeclareOperation("FirstPart",
  [ IsHomogeneousList, IsPathAlgebraVector ]
);

DeclareOperation("TipReduce",
  [ IsHomogeneousList, IsPathAlgebraVector ]
);

DeclareOperation("TipReduce",
  [ IsHomogeneousList ]
);

DeclareOperation( "LeftDivision",
    [IsPathAlgebraVector, IsPathAlgebraVector]);
