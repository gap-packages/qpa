DeclareOperation( "GBNPGroebnerBasis",
    [ IsList, IsPathAlgebra ] );

DeclareOperation( "RemainderOfDivision",
    [ IsElementOfMagmaRingModuloRelations, IsList, IsPathAlgebra ] );

DeclareOperation( "ReducedList",
    [ IsList, IsPathAlgebra ] );

DeclareOperation( "TipReducedList",
    [ IsList, IsPathAlgebra ] );

DeclareOperation( "LeftmostOccurrence",
    [ IsList, IsList ] );

DeclareSynonym( "TipWalk",
    x -> WalkOfPath(TipMonomial(x)) );

DeclareOperation( "MakeUniform",
    [ IsElementOfMagmaRingModuloRelations ] );

DeclareOperation( "MakeUniform",
    [ IsList ] );

DeclareOperation( "QPA_InArrowIdeal",
    [ IsList, IsPathAlgebra ] );

DeclareOperation( "QPA_RelationsForPathAlgebra",
    [ IsPathAlgebra ] );

