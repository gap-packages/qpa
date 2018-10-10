DeclareOperation( "HighLevelGroebnerBasis", [ IsList, IsPathAlgebra ] );
DeclareOperation( "RemainderOfDivision", [ IsElementOfMagmaRingModuloRelations, IsList, IsPathAlgebra ] );
DeclareOperation( "ReducedListQPA", [ IsList, IsPathAlgebra ] );
DeclareOperation( "TipReducedList", [ IsList, IsPathAlgebra ] );
DeclareOperation( "LeftmostOccurrence", [ IsList, IsList ] );
DeclareSynonym( "TipWalk", x -> WalkOfPath(TipMonomial(x)) );
