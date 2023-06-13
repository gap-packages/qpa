# GAP Declarations
# This file was generated from
# $Id: algpath.gd,v 1.7 2012/08/01 16:01:10 sunnyquiver Exp $
DeclareProperty( "IsPathRing", IsMagmaRingModuloSpanOfZero );
DeclareProperty( "IsPathRing", IsAlgebra );
DeclareSynonym( "IsPathAlgebra", IsAlgebra and IsPathRing );
DeclareFilter("IsElementOfPathRing");
DeclareGlobalFunction( "PathRing" );
DeclareGlobalFunction( "PathAlgebra" );
DeclareAttribute( "QuiverOfPathRing", IsPathRing );
DeclareSynonymAttr( "QuiverOfPathAlgebra", QuiverOfPathRing );
DeclareGlobalFunction("FactorPathAlgebraByRelators");
DeclareCategory( "IsQuiverAlgebra", IsAlgebra );
InstallTrueMethod( IsQuiverAlgebra, IsPathAlgebra );

DeclareOperation( "LeadingTerm", [IsRingElement] );
DeclareSynonym("Tip", LeadingTerm);
DeclareSynonym( "TipCoefficient", LeadingCoefficient );
DeclareSynonym( "TipMonomial", LeadingMonomial );

DeclareOperation( "IsLeftUniform", [IsRingElement]);
DeclareOperation( "IsLeftUniform", [IsList,IsPath]);
DeclareOperation( "IsRightUniform", [IsRingElement]);
DeclareOperation( "IsRightUniform", [IsList,IsPath]);
DeclareOperation( "IsUniform", [IsRingElement]);

DeclareCategory( "IsQuotientOfPathAlgebra", IsQuiverAlgebra );
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
DeclareAttribute( "RelationsOfAlgebra", IsQuiverAlgebra);
DeclareAttribute( "IdealOfQuotient", IsQuiverAlgebra);
DeclareAttribute( "GroebnerBasisFunction", IsQuiverAlgebra );
DeclareAttribute( "NormalFormFunction", IsFamily );
DeclareOperation( "ElementOfQuotientOfPathAlgebra", 
    [ IsElementOfQuotientOfPathAlgebraFamily, IsRingElement, IsBool ] );

DeclareOperation( "ElementOfPathAlgebra", 
    [ IsPathAlgebra, IsPath ] ); 	
	
	
DeclareHandlingByNiceBasis( "IsFpPathAlgebraElementsSpace",
    "for spaces of f.p. path algebras" );  

DeclareOperation( "BasisOfDomain", [IsQuotientOfPathAlgebra]);


# (This should have been an operation with argument filters
# [ IsElementOfPathAlgebra ], but that filter does not exist).
DeclareGlobalFunction( "PathAlgebraContainingElement");

DeclareAttribute( "OriginalPathAlgebra", IsQuiverAlgebra ); 
DeclareOperation( "MakeUniformOnRight", [ IsHomogeneousList ] );
DeclareOperation( "GeneratorsTimesArrowsOnRight", [ IsHomogeneousList ] );
DeclareOperation( "NthPowerOfArrowIdeal", [ IsPathAlgebra, IS_INT ] );
DeclareOperation( "TruncatedPathAlgebra", [ IsField, IsQuiver, IS_INT ] );
DeclareOperation( "AddNthPowerToRelations", [ IsPathAlgebra, IsHomogeneousList, IS_INT ] );

DeclareAttribute( "OppositePathAlgebra", IsAlgebra );
DeclareGlobalFunction( "OppositePathAlgebraElement"); # should be operation with args [ IsElementOfPathAlgebra ]
DeclareOperation( "OppositeRelations", [ IsDenseList ] );
DeclareOperation( "VertexPosition", [ IsElementOfQuotientOfPathAlgebra ] );

DeclareProperty( "IsSelfinjectiveAlgebra",  IsQuiverAlgebra ); 
DeclareProperty( "IsSymmetricAlgebra", IsQuiverAlgebra );
DeclareProperty( "IsWeaklySymmetricAlgebra", IsQuiverAlgebra );
DeclareProperty( "IsSchurianAlgebra", IsQuiverAlgebra );
DeclareProperty( "IsSemicommutativeAlgebra", IsPathAlgebra);
DeclareProperty( "IsSemicommutativeAlgebra", IsQuotientOfPathAlgebra);
InstallTrueMethod( IsHereditaryAlgebra, IsPathAlgebra );
InstallTrueMethod( IsFiniteGlobalDimensionAlgebra, IsPathAlgebra );

DeclareAttribute( "CoxeterPolynomial",  IsQuiverAlgebra  ); 
DeclareAttribute( "CoxeterMatrix", IsQuiverAlgebra ); 
DeclareOperation( "TipMonomialandCoefficientOfVector", [ IsQuiverAlgebra, IsCollection ] );
DeclareOperation( "TipReduceVectors", [ IsQuiverAlgebra, IsCollection ] );
DeclareOperation( "CoefficientsOfVectors", [ IsAlgebra, IsCollection, IsList ] );
DeclareProperty( "IsDistributiveAlgebra", IsQuotientOfPathAlgebra );
DeclareAttribute( "NakayamaAutomorphism", IsQuotientOfPathAlgebra ); 
DeclareAttribute( "NakayamaPermutation", IsQuotientOfPathAlgebra );
DeclareAttribute( "FrobeniusForm", IsQuotientOfPathAlgebra );
DeclareAttribute( "FrobeniusLinearFunctional", IsQuotientOfPathAlgebra );
DeclareAttribute( "OrderOfNakayamaAutomorphism", IsQuotientOfPathAlgebra );
DeclareAttribute( "AssociatedMonomialAlgebra", IsQuiverAlgebra );
DeclareProperty( "IsMonomialAlgebra", IsQuiverAlgebra );
DeclareOperation( "SaveAlgebra", [ IsQuiverAlgebra, IsString, IsString ]);
DeclareOperation( "ReadAlgebra", [ IsString ]);
DeclareProperty( "IsTriangularReduced", IsQuiverAlgebra );
DeclareOperation( "PathRemoval", [ IsElementOfMagmaRingModuloRelations, IsList ] );
DeclareOperation( "QuiverAlgebraOfAmodAeA", [ IsQuiverAlgebra, IsList ] );
DeclareOperation( "QuiverAlgebraOfeAe", [ IsQuiverAlgebra, IsObject ] );
DeclareOperation( "SupportOfQuiverAlgebraElement", [ IsQuiverAlgebra, IsObject ] );
DeclareOperation( "LeftSupportOfQuiverAlgebraElement", [ IsQuiverAlgebra, IsObject ] );
DeclareOperation( "RightSupportOfQuiverAlgebraElement", [ IsQuiverAlgebra, IsObject ] );