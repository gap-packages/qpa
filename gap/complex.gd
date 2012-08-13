DeclareCategory( "IsCat", IsRecord );
DeclareRepresentation( "IsCatDefaultRep", IsComponentObjectRep,
                       [ "identity", "properties" ] );

#DeclareOperation( "\.", [ IsCat, IsString ] );

DeclareGlobalFunction( "Cat" );

DeclareOperation( "CatOfRightAlgebraModules", [ IsAlgebra ] );



DeclareCategory( "IsComplex", IsObject );
DeclareCategoryFamily( "IsComplex" );

DeclareCategory( "IsZeroComplex", IsComplex );

DeclareRepresentation( "IsComplexDefaultRep",
                       IsComponentObjectRep and IsAttributeStoringRep,
                       [ "differentials" ] );

DeclareAttribute( "CatOfComplex", IsComplex );
DeclareOperation( "ObjectOfComplex", [ IsComplex, IsInt ] );
DeclareOperation( "DifferentialOfComplex", [ IsComplex, IsInt ] );
DeclareOperation( "\^", [ IsComplex, IsInt ] );
DeclareAttribute( "DifferentialsOfComplex", IsComplex );
DeclareOperation( "CyclesOfComplex", [ IsComplex, IsInt ] );
DeclareOperation( "BoundariesOfComplex", [ IsComplex, IsInt ] );
DeclareOperation( "HomologyOfComplex", [ IsComplex, IsInt ] );
DeclareOperation( "UpperBound", [ IsComplex ] );
DeclareOperation( "LowerBound", [ IsComplex ] );
DeclareOperation( "IsFiniteComplex", [ IsComplex ] );
DeclareOperation( "LengthOfComplex", [ IsComplex ] );
DeclareOperation( "HighestKnownDegree", [ IsComplex ] );
DeclareOperation( "LowestKnownDegree", [ IsComplex ] );
DeclareProperty( "IsExactSequence", IsComplex );
DeclareOperation( "IsExactInDegree", [ IsComplex, IsInt ] );
DeclareProperty( "IsShortExactSequence", IsComplex );
DeclareOperation( "ForEveryDegree", [ IsComplex, IsFunction ] );
DeclareOperation( "IsPositiveRepeating", [ IsComplex ] );
DeclareOperation( "IsNegativeRepeating", [ IsComplex ] );
DeclareOperation( "PositiveRepeatDegrees", [ IsComplex ] );
DeclareOperation( "NegativeRepeatDegrees", [ IsComplex ] );

DeclareOperation( "Shift", [ IsComplex, IsInt ] );
DeclareOperation( "ShiftUnsigned", [ IsComplex, IsInt ] );
DeclareOperation( "YonedaProduct", [ IsComplex, IsComplex ] );

DeclareOperation( "GoodTruncationBelow", [ IsComplex, IsInt ] );
DeclareOperation( "GoodTruncationAbove", [ IsComplex, IsInt ] );
DeclareOperation( "GoodTruncation", [ IsComplex, IsInt, IsInt ] );
DeclareOperation( "BrutalTruncationBelow", [ IsComplex, IsInt ] );
DeclareOperation( "BrutalTruncationAbove", [ IsComplex, IsInt ] );
DeclareOperation( "BrutalTruncation", [ IsComplex, IsInt, IsInt ] );
DeclareOperation( "SyzygyTruncation", [ IsComplex, IsInt ] );
DeclareOperation( "CosyzygyTruncation", [ IsComplex, IsInt ] );
DeclareOperation( "SyzygyCosyzygyTruncation", [ IsComplex, IsInt, IsInt ] );
DeclareOperation( "CutComplexAbove", [ IsComplex ] );
DeclareOperation( "CutComplexBelow", [ IsComplex ] );

DeclareGlobalFunction( "Complex" );
# Complex( cat, basePosition, differentials, [ "repeat", [ f, g, h ] ], "zero" );
# Complex( cat, basePosition, differentials, [ "next", function( d ) ... end ], "zero" );
# Complex( cat, basePosition, differentials, [ "pos", function( C, i ) ... end,  ], "zero" );
# - If positive or negative is "zero", then middle must be nonempty.

DeclareGlobalFunction( "FiniteComplex" );
DeclareGlobalFunction( "ZeroComplex" );
DeclareGlobalFunction( "StalkComplex" );
DeclareGlobalFunction( "ShortExactSequence" );
DeclareGlobalFunction( "ComplexByDifferentialList" );

DeclareOperation( "ProjectiveResolution", [ IsAlgebraModule ] );

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

DeclareGlobalFunction( "ChainMap" );
DeclareGlobalFunction( "FiniteChainMap" );
DeclareGlobalFunction( "ZeroChainMap" );

DeclareGlobalFunction( "ComplexAndChainMaps" );

# ComplexAndChainMaps( sourceComplexes, rangeComplexes,
#                      basePosition, middle, positive, negative );

# ComplexAndChainMaps( [ C ], [], 0, [ [ g, f ] ], [ "pos", func ], "zero" );

