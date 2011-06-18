# GAP Declarations
# This file was generated from
# $Id: algpath.gd,v 1.4 2011/06/18 13:51:46 sunnyquiver Exp $
DeclareProperty( "IsPathRing", IsMagmaRingModuloSpanOfZero );
DeclareProperty( "IsPathRing", IsAlgebra );
DeclareSynonym( "IsPathAlgebra", IsAlgebra and IsPathRing );
DeclareFilter("IsElementOfPathRing");
DeclareGlobalFunction( "PathRing" );
DeclareGlobalFunction( "PathAlgebra" );
DeclareAttribute( "QuiverOfPathRing", IsPathRing );
DeclareSynonymAttr( "QuiverOfPathAlgebra", QuiverOfPathRing );
DeclareGlobalFunction("FactorPathAlgebraByRelators");

DeclareOperation( "LeadingTerm", [IsRingElement] );
DeclareSynonym("Tip", LeadingTerm);
DeclareSynonym( "TipCoefficient", LeadingCoefficient );
DeclareSynonym( "TipMonomial", LeadingMonomial );

DeclareOperation( "IsLeftUniform", [IsRingElement]);
DeclareOperation( "IsLeftUniform", [IsList,IsPath]);
DeclareOperation( "IsRightUniform", [IsRingElement]);
DeclareOperation( "IsRightUniform", [IsList,IsPath]);
DeclareOperation( "IsUniform", [IsRingElement]);

DeclareCategory( "IsSubalgebraFpPathAlgebra", IsAlgebra );
DeclareAttribute( "OrderingOfAlgebra", IsAlgebra );
DeclareAttribute( "GroebnerBasisOfIdeal", IsRing );
DeclareAttribute( "GroebnerBasisOfLeftIdeal", IsRing );
DeclareAttribute( "GroebnerBasisOfRightIdeal", IsRing );

DeclareCategory( "IsElementOfFpPathAlgebra", IsRingElement );
DeclareCategoryCollections( "IsElementOfFpPathAlgebra" );
DeclareCategoryFamily( "IsElementOfFpPathAlgebra" );

DeclareInfoClass( "InfoElementOfFpPathAlgebra" );
DeclareProperty("IsFullFpPathAlgebra", 
    IsFLMLOR and IsElementOfFpPathAlgebraCollection );
DeclareAttribute( "RelatorsOfFpAlgebra",
    IsSubalgebraFpPathAlgebra and IsFullFpPathAlgebra);
DeclareAttribute( "NormalFormFunction", IsFamily );
DeclareOperation( "ElementOfFpPathAlgebra", 
    [ IsElementOfFpPathAlgebraFamily, IsRingElement, IsBool ] );

DeclareHandlingByNiceBasis( "IsFpPathAlgebraElementsSpace",
    "for spaces of f.p. path algebras" );  

DeclareOperation( "BasisOfDomain", [IsSubalgebraFpPathAlgebra]);


# (This should have been an operation with argument filters
# [ IsElementOfPathAlgebra ], but that filter does not exist).
DeclareGlobalFunction( "PathAlgebraContainingElement");

DeclareOperation( "OriginalPathAlgebra", [ IsAlgebra ] ); 
DeclareOperation( "MakeUniformOnRight", [ IsHomogeneousList ] );
DeclareOperation( "GeneratorsTimesArrowsOnRight", [ IsHomogeneousList ] );
DeclareOperation( "nthPowerOfArrowIdeal", [ IsPathAlgebra, IS_INT ] );
DeclareOperation( "TruncatedPathAlgebra", [ IsField, IsQuiver, IS_INT ] );
DeclareOperation( "AddNthPowerToRelations", [ IsPathAlgebra, IsHomogeneousList, IS_INT ] );
