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
DeclareOperation( "YonedaProduct", [ IsComplex, IsComplex ] );

DeclareGlobalFunction( "Complex" );
# Complex( cat, basePosition, differentials, [ "repeat", [ f, g, h ] ], "zero" );
# Complex( cat, basePosition, differentials, [ "next", function( d ) ... end ], "zero" );
# Complex( cat, basePosition, differentials, [ "pos", function( C, i ) ... end,  ], "zero" );
# - If positive or negative is "zero", then middle must be nonempty.

DeclareGlobalFunction( "FiniteComplex" );
DeclareGlobalFunction( "ZeroComplex" );
DeclareGlobalFunction( "SingleObjectComplex" );
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

