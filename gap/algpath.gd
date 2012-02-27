# GAP Declarations
# This file was generated from
# $Id: algpath.gd,v 1.5 2012/02/27 12:26:34 sunnyquiver Exp $
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

DeclareCategory( "IsQuotientOfPathAlgebra", IsAlgebra );
DeclareAttribute( "OrderingOfAlgebra", IsAlgebra );
DeclareAttribute( "GroebnerBasisOfIdeal", IsRing );
DeclareAttribute( "GroebnerBasisOfLeftIdeal", IsRing );
DeclareAttribute( "GroebnerBasisOfRightIdeal", IsRing );

DeclareCategory( "IsElementOfQuotientOfPathAlgebra", IsRingElement );
DeclareCategoryCollections( "IsElementOfQuotientOfPathAlgebra" );
DeclareCategoryFamily( "IsElementOfQuotientOfPathAlgebra" );

DeclareInfoClass( "InfoElementOfQuotientOfPathAlgebra" );
DeclareProperty("IsFullFpPathAlgebra", 
    IsFLMLOR and IsElementOfQuotientOfPathAlgebraCollection );
DeclareAttribute( "RelatorsOfFpAlgebra",
    IsQuotientOfPathAlgebra and IsFullFpPathAlgebra);
DeclareAttribute( "NormalFormFunction", IsFamily );
DeclareOperation( "ElementOfQuotientOfPathAlgebra", 
    [ IsElementOfQuotientOfPathAlgebraFamily, IsRingElement, IsBool ] );

DeclareHandlingByNiceBasis( "IsFpPathAlgebraElementsSpace",
    "for spaces of f.p. path algebras" );  

DeclareOperation( "BasisOfDomain", [IsQuotientOfPathAlgebra]);


# (This should have been an operation with argument filters
# [ IsElementOfPathAlgebra ], but that filter does not exist).
DeclareGlobalFunction( "PathAlgebraContainingElement");

DeclareOperation( "OriginalPathAlgebra", [ IsAlgebra ] ); 
DeclareOperation( "MakeUniformOnRight", [ IsHomogeneousList ] );
DeclareOperation( "GeneratorsTimesArrowsOnRight", [ IsHomogeneousList ] );
DeclareOperation( "NthPowerOfArrowIdeal", [ IsPathAlgebra, IS_INT ] );
DeclareOperation( "TruncatedPathAlgebra", [ IsField, IsQuiver, IS_INT ] );
DeclareOperation( "AddNthPowerToRelations", [ IsPathAlgebra, IsHomogeneousList, IS_INT ] );
