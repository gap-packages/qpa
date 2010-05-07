# Projective Resolutions File
# This file was generated from
# $Id: projres.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
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
  [ IsHomogeneousList, IsObject ]
);

DeclareOperation("TipReduce",
  [ IsHomogeneousList ]
);

DeclareOperation( "LeftDivision",
    [IsPathAlgebraVector, IsPathAlgebraVector]);
