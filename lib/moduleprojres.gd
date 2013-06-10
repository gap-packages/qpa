# Projective Resolutions File
# This file was generated from
# $Id: projres.gd,v 1.5 2012/09/27 08:55:07 sunnyquiver Exp $
DeclareInfoClass( "InfoProjectiveResolutionFpPathAlgebraModule" );

DeclareCategory("IsProjectiveResolutionFpPathAlgebraModule",IsObject);
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

DeclareOperation("ProjectiveResolutionOfPathAlgebraModule",
  [ IsPathAlgebraMatModule, IsPosInt ]
);

DeclareOperation("FindNextRProjective",
  [ IsProjectiveResolutionFpPathAlgebraModule, IsRing, IsPosInt ]
);

DeclareOperation("FindNextSyzygy",
  [ IsProjectiveResolutionFpPathAlgebraModule, IsPosInt ]
);

DeclareOperation("XSetOfPathAlgebraVector",
  [ IsHomogeneousList, IsHomogeneousList,
    IsRing and HasGroebnerBasisOfIdeal,
    IsPathAlgebraVector ]
);

DeclareOperation("FirstPart",
  [ IsHomogeneousList, IsHomogeneousList, IsPathAlgebraVector ]
);

DeclareOperation("TipReduce",
  [ IsHomogeneousList, IsPathAlgebraVector ]
);

DeclareOperation("TipReduce",
  [ IsHomogeneousList ]
);

DeclareOperation( "LeftDivision",
    [IsPathAlgebraVector, IsPathAlgebraVector]);
