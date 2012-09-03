# Projective Resolutions File
# This file was generated from
# $Id: projres.gd,v 1.3 2012/09/03 07:06:45 sunnyquiver Exp $
DeclareInfoClass( "InfoProjectiveResolutionFpPathAlgebraModule" );

DeclareCategory("IsProjectiveResolutionFpPathAlgebraModule",IsObject);
DeclareAttribute("Module", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("ParentAlgebra", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("RingIdeal", IsProjectiveResolutionFpPathAlgebraModule);
DeclareAttribute("Maps", IsProjectiveResolutionFpPathAlgebraModule,"mutable");
DeclareAttribute("Projectives", IsProjectiveResolutionFpPathAlgebraModule,"mutable");
DeclareAttribute("RProjectives", IsProjectiveResolutionFpPathAlgebraModule,"mutable");
DeclareAttribute("RProjectivesVertexList", IsProjectiveResolutionFpPathAlgebraModule,"mutable");
DeclareAttribute("ProjectivesFList", IsProjectiveResolutionFpPathAlgebraModule,"mutable");
DeclareAttribute("ProjectivesFPrimeList", IsProjectiveResolutionFpPathAlgebraModule,"mutable");

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

DeclareOperation("FindNextRProjective",
  [ IsProjectiveResolutionFpPathAlgebraModule, IsRing, IsPosInt ]
);

DeclareOperation("FindNextSyzygy",
  [ IsProjectiveResolutionFpPathAlgebraModule, IsPosInt ]
);

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
