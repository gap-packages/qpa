# GAP Declarations
# This file was generated from
# $Id: present2.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $

DeclareCategory("IsProjectivePresentation",IsObject);
DeclareAttribute("ParentAlgebra", IsProjectivePresentation);
DeclareAttribute("MatrixRepresentation", IsProjectivePresentation);
DeclareAttribute("MatrixNumCols", IsProjectivePresentation);
DeclareAttribute("MatrixNumRows", IsProjectivePresentation);
#DeclareAttribute("P0", IsProjectivePresentation);
#DeclareAttribute("P1", IsProjectivePresentation);
DeclareAttribute("P0VertexList", IsProjectivePresentation);
DeclareAttribute("P1VertexList", IsProjectivePresentation);


#############################################################################
##
#R  IsProjectivePresentation
##
##  The following two lines create the representation and family for
##  our Projective Presentation objects.
##
DeclareRepresentation("IsProjectivePresentationDefaultRep",
  IsAttributeStoringRep and IsPositionalObjectRep
  and IsProjectivePresentation,[]);

BindGlobal("ProjectivePresentationFamily",
  NewFamily("ProjectivePresentationFamily",
  IsProjectivePresentation, IsProjectivePresentation) );


DeclareOperation("ProjectivePresentation", [IsAlgebra,IsMatrix]);


##
#
DeclareProperty( "IsPathAlgebraModule", IsAlgebraModule );
DeclareProperty( "IsFpPathAlgebraModule", IsAlgebraModule );

