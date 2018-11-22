DeclareCategory( "IsCat", IsRecord );
DeclareRepresentation( "IsCatDefaultRep", IsComponentObjectRep,
                       [ "identity", "properties" ] );

#DeclareOperation( "\.", [ IsCat, IsString ] );

DeclareGlobalFunction( "Cat" );

DeclareOperation( "CatOfRightAlgebraModules", [ IsAlgebra ] );



DeclareCategory( "IsQPAComplex", IsObject );
DeclareCategoryFamily( "IsQPAComplex" );

DeclareCategory( "IsZeroComplex", IsQPAComplex );

DeclareRepresentation( "IsQPAComplexDefaultRep",
                       IsComponentObjectRep and IsAttributeStoringRep,
                       [ "differentials" ] );

DeclareAttribute( "CatOfComplex", IsQPAComplex );
DeclareOperation( "ObjectOfComplex", [ IsQPAComplex, IsInt ] );
DeclareOperation( "DifferentialOfComplex", [ IsQPAComplex, IsInt ] );
DeclareOperation( "\^", [ IsQPAComplex, IsInt ] );
DeclareAttribute( "DifferentialsOfComplex", IsQPAComplex );
DeclareOperation( "CyclesOfComplex", [ IsQPAComplex, IsInt ] );
DeclareOperation( "BoundariesOfComplex", [ IsQPAComplex, IsInt ] );
DeclareOperation( "HomologyOfComplex", [ IsQPAComplex, IsInt ] );
DeclareOperation( "UpperBound", [ IsQPAComplex ] );
DeclareOperation( "LowerBound", [ IsQPAComplex ] );
DeclareOperation( "IsFiniteComplex", [ IsQPAComplex ] );
DeclareOperation( "LengthOfComplex", [ IsQPAComplex ] );
DeclareOperation( "HighestKnownDegree", [ IsQPAComplex ] );
DeclareOperation( "LowestKnownDegree", [ IsQPAComplex ] );
DeclareProperty( "IsExactSequence", IsQPAComplex );
DeclareOperation( "IsExactInDegree", [ IsQPAComplex, IsInt ] );
DeclareProperty( "IsShortExactSequence", IsQPAComplex );
DeclareOperation( "ForEveryDegree", [ IsQPAComplex, IsFunction ] );
DeclareOperation( "IsPositiveRepeating", [ IsQPAComplex ] );
DeclareOperation( "IsNegativeRepeating", [ IsQPAComplex ] );
DeclareOperation( "PositiveRepeatDegrees", [ IsQPAComplex ] );
DeclareOperation( "NegativeRepeatDegrees", [ IsQPAComplex ] );

DeclareOperation( "Shift", [ IsQPAComplex, IsInt ] );
DeclareOperation( "ShiftUnsigned", [ IsQPAComplex, IsInt ] );
DeclareOperation( "YonedaProduct", [ IsQPAComplex, IsQPAComplex ] );

DeclareOperation( "GoodTruncationBelow", [ IsQPAComplex, IsInt ] );
DeclareOperation( "GoodTruncationAbove", [ IsQPAComplex, IsInt ] );
DeclareOperation( "GoodTruncation", [ IsQPAComplex, IsInt, IsInt ] );
DeclareOperation( "BrutalTruncationBelow", [ IsQPAComplex, IsInt ] );
DeclareOperation( "BrutalTruncationAbove", [ IsQPAComplex, IsInt ] );
DeclareOperation( "BrutalTruncation", [ IsQPAComplex, IsInt, IsInt ] );
DeclareOperation( "SyzygyTruncation", [ IsQPAComplex, IsInt ] );
DeclareOperation( "CosyzygyTruncation", [ IsQPAComplex, IsInt ] );
DeclareOperation( "SyzygyCosyzygyTruncation", [ IsQPAComplex, IsInt, IsInt ] );
DeclareOperation( "CutComplexAbove", [ IsQPAComplex ] );
DeclareOperation( "CutComplexBelow", [ IsQPAComplex ] );

DeclareGlobalFunction( "Complex" );
# Complex( cat, basePosition, differentials, [ "repeat", [ f, g, h ] ], "zero" );
# Complex( cat, basePosition, differentials, [ "next", function( d ) ... end ], "zero" );
# Complex( cat, basePosition, differentials, [ "pos", function( C, i ) ... end,  ], "zero" );
# - If positive or negative is "zero", then middle must be nonempty.

DeclareGlobalFunction( "FiniteComplex" );
DeclareGlobalFunction( "ZeroComplex" );
DeclareGlobalFunction( "StalkComplex" );
DeclareGlobalFunction( "ShortExactSequence" );
DeclareOperation( "ComplexByDifferentialList", [ IsCat, IsInfList ] );

DeclareOperation( "ProjectiveResolution", [ IsAlgebraModule ] );
DeclareOperation( "InjectiveResolution", [ IsAlgebraModule ] );

DeclareCategory( "IsChainMap", IsObject );
DeclareCategoryFamily( "IsChainMap" );

DeclareRepresentation( "IsChainMapDefaultRep",
                       IsComponentObjectRep and IsAttributeStoringRep,
                       [] );

DeclareAttribute( "Source", IsChainMap );
DeclareAttribute( "Range", IsChainMap );
DeclareAttribute( "MorphismsOfChainMap", IsChainMap );
DeclareOperation( "MorphismOfChainMap", [ IsChainMap, IsInt ] );
DeclareOperation( "\^", [ IsChainMap, IsInt ] );

DeclareOperation( "HighestKnownDegree", [ IsChainMap ] );
DeclareOperation( "LowestKnownDegree", [ IsChainMap ] );

DeclareOperation( "ChainMap",
                  [ IsQPAComplex, IsQPAComplex, IsInt, IsList, IsList, IsList ] );
DeclareOperation( "FiniteChainMap", [ IsQPAComplex, IsQPAComplex, IsInt, IsList ] );
DeclareOperation( "ZeroChainMap" , [ IsQPAComplex, IsQPAComplex ] );

DeclareGlobalFunction( "ComplexAndChainMaps" );

# ComplexAndChainMaps( sourceComplexes, rangeComplexes,
#                      basePosition, middle, positive, negative );

# ComplexAndChainMaps( [ C ], [], 0, [ [ g, f ] ], [ "pos", func ], "zero" );

DeclareOperation( "ComparisonLifting", [ IsPathAlgebraMatModuleHomomorphism, IsQPAComplex, IsQPAComplex ] );
DeclareOperation( "ComparisonLiftingToProjectiveResolution", [ IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "MappingCone", [ IsChainMap ] );