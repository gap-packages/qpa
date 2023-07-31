BindGlobal("PositiveInfinity", infinity); 
BindGlobal("NegativeInfinity", -infinity);

DeclareCategory( "IsInfList", IsObject );
DeclareCategory( "IsHalfInfList", IsObject );

DeclareOperation( "\^", [ IsInfList, IsInt ] );
DeclareOperation( "\^", [ IsHalfInfList, IsInt ] );

DeclareOperation( "MiddleStart", [ IsInfList ] );
DeclareOperation( "MiddleEnd", [ IsInfList ] );
DeclareOperation( "MiddlePart", [ IsInfList ] );
DeclareOperation( "PositivePart", [ IsInfList ] );
DeclareOperation( "NegativePart", [ IsInfList ] );
DeclareOperation( "LowestKnownPosition", [ IsInfList ] );
DeclareOperation( "HighestKnownPosition", [ IsInfList ] );
DeclareOperation( "UpperBound", [ IsInfList ] );
DeclareOperation( "LowerBound", [ IsInfList ] );

DeclareOperation( "StartPosition", [ IsHalfInfList ] );
DeclareOperation( "Direction", [ IsHalfInfList ] );
DeclareOperation( "InfListType", [ IsHalfInfList ] );
DeclareOperation( "RepeatingList", [ IsHalfInfList ] );
DeclareOperation( "ElementFunction", [ IsHalfInfList ] );
DeclareOperation( "IsStoringValues", [ IsHalfInfList ] );
DeclareOperation( "NewValueCallback", [ IsHalfInfList ] );
DeclareOperation( "IsRepeating", [ IsHalfInfList ] );
DeclareOperation( "InitialValue", [ IsHalfInfList ] );
DeclareOperation( "LowestKnownPosition", [ IsHalfInfList ] );
DeclareOperation( "HighestKnownPosition", [ IsHalfInfList ] );
DeclareOperation( "MakeRepeatingList", [ IsHalfInfList, IsPosInt, IsPosInt ] );

DeclareOperation( "FinitePartAsList", [ IsInfList, IsInt, IsInt ] );
DeclareOperation( "PositivePartFrom", [ IsInfList, IsInt ] );
DeclareOperation( "NegativePartFrom", [ IsInfList, IsInt ] );
DeclareOperation( "Splice", [ IsInfList, IsInfList, IsInt ] );
DeclareOperation( "ShiftedSplice", [ IsInfList, IsInt, IsInfList, IsInt, IsInt ] );
DeclareOperation( "Shift", [ IsHalfInfList, IsInt ] );
DeclareOperation( "Shift", [ IsInfList, IsInt ] );
DeclareOperation( "Cut", [ IsHalfInfList, IsInt ] );
DeclareGlobalFunction( "InfConcatenation" );

DeclareOperation( "InfList", [ IsInfList, IsFunction ] );
DeclareOperation( "HalfInfList", [ IsHalfInfList, IsFunction ] );

DeclareRepresentation( "IsHalfInfListDefaultRep",
                       IsComponentObjectRep and IsAttributeStoringRep,
                       [ "values", "start", "direction",
                         "type", "func", "repeatingList",
                         "storingValues", "initialValue",
                         "callback" ] );

DeclareRepresentation( "IsInfListDefaultRep",
                       IsComponentObjectRep and IsAttributeStoringRep,
                       [ "basePosition", "middle", "positive", "negative" ] );

DeclareGlobalFunction( "MakeHalfInfList" );

# MakeHalfInfList( start, direction, typeWithArgs, callback )
# MakeHalfInfList( start, direction, [ "repeat", repeatList ], callback )
# MakeHalfInfList( start, direction, [ "next", nextFunction, initialValue ], callback )
# MakeHalfInfList( start, direction, [ "pos", posFunction ], callback )
# MakeHalfInfList( start, direction, [ "pos", posFunction, storeValues ], callback )

DeclareGlobalFunction( "MakeInfList" );

# MakeInfList( basePosition, middle, positive, negative, callback )
# MakeInfList( basePosition, middle, [ "repeat", repeatList ], [ "repeat", repeatList ],
#              callback )
# MakeInfList( basePosition, middle,
#              [ "next", nextFunction, initialValue ],
#              [ "repeat", repeatList ], callback )
# ...

DeclareGlobalFunction( "MakeInfListFromHalfInfLists" );

# MakeInfListFromHalfInfLists( middle, positive, negative )

DeclareGlobalFunction( "FunctionInfList" );
DeclareGlobalFunction( "ConstantInfList" );
DeclareGlobalFunction( "FiniteInfList" );
