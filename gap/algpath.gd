# GAP Declarations
# This file was generated from
# $Id: algpath.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
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
