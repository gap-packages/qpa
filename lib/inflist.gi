InstallValue( PositiveInfinity,
              Objectify( NewType( FamilyObj(0),
                                  IsInfiniteNumber and IsPositionalObjectRep ),
                         [] ) );
InstallValue( NegativeInfinity,
              Objectify( NewType( FamilyObj(0),
                                  IsInfiniteNumber and IsPositionalObjectRep ),
                         [] ) );

InstallMethod( \=,
    [ IsInfiniteNumber, IsInfiniteNumber ],
    IsIdenticalObj );
InstallMethod( \=,
    [ IsInfiniteNumber, IsObject ],
    ReturnFalse );
InstallMethod( \=,
    [ IsObject, IsInfiniteNumber ],
    ReturnFalse );

InstallMethod( \<,
[ IsInfiniteNumber, IsInt ],
function( infnum, num )
    return infnum = NegativeInfinity;
end );
InstallMethod( \<,
[ IsInt, IsInfiniteNumber ],
function( num, infnum )
        return infnum = PositiveInfinity;
    end );
InstallMethod( \<,
[ IsInfiniteNumber, IsInfiniteNumber ],
function( n1, n2 )
        return n1 = NegativeInfinity and n2 = PositiveInfinity;
    end );

InstallMethod( PrintObj,
[ IsInfiniteNumber ],
function( n )
        if n = PositiveInfinity then
            Print( "+inf" );
	else
            Print( "-inf" );
	fi;
    end );



InstallGlobalFunction( MakeHalfInfList,
function( start, direction, typeWithArgs, callback )
    local type, repeatingList, func, initialValue, storeValues, data, list;

    type := typeWithArgs[ 1 ];
    if type = "repeat" then
        repeatingList := typeWithArgs[ 2 ];
        func := fail;
        storeValues := false;
    else
        repeatingList := fail;
        func := typeWithArgs[ 2 ];
        storeValues := true;
    fi;
    initialValue := fail;
    if type = "next" then
        initialValue := typeWithArgs[ 3 ];
    elif type = "pos" and Length( typeWithArgs ) = 3 then
        storeValues := typeWithArgs[ 3 ];
    fi;

    data := rec( values := [],
                 start := start,
                 direction := direction,
                 type := type,
                 func := func,
                 repeatingList := repeatingList,
                 storingValues := storeValues,
                 initialValue := initialValue,
                 callback := callback );

    list := Objectify( NewType( NewFamily( "HalfInfListsFamily" ),
                                IsHalfInfList and IsHalfInfListDefaultRep ),
                       data );

    return list;

end );

InstallMethod( StartPosition,
[ IsHalfInfList ],
function( list )
    return list!.start;
end );

InstallMethod( Direction,
[ IsHalfInfList ],
function( list )
    return list!.direction;
end );

InstallMethod( InfListType,
[ IsHalfInfList ],
function( list )
    return list!.type;
end );

InstallMethod( RepeatingList,
[ IsHalfInfList ],
function( list )
    return list!.repeatingList;
end );

InstallMethod( ElementFunction,
[ IsHalfInfList ],
function( list )
    return list!.func;
end );

InstallMethod( IsStoringValues,
[ IsHalfInfList ],
function( list )
    return list!.storingValues;
end );

InstallMethod( NewValueCallback,
[ IsHalfInfList ],
function( list )
    return list!.callback;
end );

InstallMethod( IsRepeating,
[ IsHalfInfList ],
function( list )
    return list!.type = "repeat";
end );

InstallMethod( InitialValue,
[ IsHalfInfList ],
function( list )
    return list!.initialValue;
end );

InstallMethod( \^,
[ IsHalfInfList and IsHalfInfListDefaultRep, IsInt ],
function( list, pos )
    local index, i, callCallbackFunction;

    index := (pos - StartPosition( list )) * Direction( list );
    
    if index < 0 then
        Error( "illegal half inf list position ", pos, "\n" );
    fi;

    callCallbackFunction := function( index )
        local pos, cb;
        pos := StartPosition( list ) + index * Direction( list );
        cb := NewValueCallback( list );
        if IsFunction( cb ) then
            cb( pos, Direction( list ), InfListType( list ) );
        fi;
    end;

    if InfListType( list ) = "repeat" then
        index := index mod Length( RepeatingList( list ) );
        return RepeatingList( list )[ index + 1 ];
    elif InfListType( list ) = "next" then
        if Length( list!.values ) = 0 then
            list!.values[ 1 ] := ElementFunction( list )( InitialValue( list ) );
            callCallbackFunction( 0 );
        fi;
        for i in [ Length( list!.values ) .. index ] do
            list!.values[ i + 1 ] := ElementFunction( list )( list!.values[ i ] );
            callCallbackFunction( i );
        od;
        return list!.values[ index + 1 ];
    elif InfListType( list ) = "pos" then
        if IsStoringValues( list ) then
            for i in [ Length( list!.values ) .. index ] do
                list!.values[ i + 1 ]
                  := ElementFunction( list )( StartPosition( list )
                                              + i * Direction( list ) );
                callCallbackFunction( i );
            od;
            return list!.values[ index + 1 ];
        else
            return ElementFunction( list )( pos );
        fi;
    fi;
    
end );

InstallMethod( LowestKnownPosition,
[ IsHalfInfList ],
function( list )
    if Direction( list ) = 1 then
        if Length( list!.values ) = 0 and InfListType( list ) <> "repeat" then
            return "none";
        else
            return StartPosition( list );
        fi;
    else
        if InfListType( list ) = "repeat" then
            return NegativeInfinity;
        elif Length( list!.values ) = 0 and InfListType( list ) <> "repeat" then
            return "none";
        else
            return StartPosition( list ) - (Length( list!.values ) - 1);
        fi;
    fi;
end );

InstallMethod( HighestKnownPosition,
[ IsHalfInfList ],
function( list )
    if Direction( list ) = -1 then
        if Length( list!.values ) = 0 and InfListType( list ) <> "repeat" then
            return "none";
        else
            return StartPosition( list );
        fi;
    else
        if InfListType( list ) = "repeat" then
            return PositiveInfinity;
        elif Length( list!.values ) = 0 and InfListType( list ) <> "repeat" then
            return "none";
        else
            return StartPosition( list ) + (Length( list!.values ) - 1);
        fi;
    fi;
end );




InstallGlobalFunction( MakeInfList,
function( basePosition, middle, positive, negative, callback )
    local posList, negList;

    # TODO: automatic initialValue for "next" parts if middle is nonempty

    # TODO: allow empty positive/negative?

    posList := MakeHalfInfList( basePosition + Length( middle ),
                                1, positive, callback );
    negList := MakeHalfInfList( basePosition - 1,
                                -1, negative, callback );

    return MakeInfListFromHalfInfLists( basePosition, middle, posList, negList );

end );

InstallGlobalFunction( MakeInfListFromHalfInfLists,
function( basePosition, middle, positive, negative )
    local list;

    if negative <> fail and basePosition <> StartPosition( negative ) + 1 then
        Error( "negative list starts at incorrect position" );
    fi;
    if positive <> fail and basePosition + Length( middle ) <> StartPosition( positive ) then
        Error( "positive list starts at incorrect position" );
    fi;

    list := Objectify( NewType( NewFamily( "InfListsFamily" ),
                                IsInfList and IsInfListDefaultRep ),
                       rec( basePosition := basePosition,
                            middle := middle,
                            positive := positive,
                            negative := negative ) );

    return list;

end );

InstallGlobalFunction( FunctionInfList,
function( func )
    local positiveNegativeList;
    positiveNegativeList := [ "pos", func, false ];
    return MakeInfList( 0, [], positiveNegativeList, positiveNegativeList, false );
end );

InstallGlobalFunction( ConstantInfList,
function( value )
    return MakeInfList( 0, [], [ "repeat", [ value ] ], [ "repeat", [ value ] ], false );
end );

InstallGlobalFunction( FiniteInfList,
function( basePosition, list )
    return MakeInfListFromHalfInfLists( basePosition, list, fail, fail );
end );

InstallMethod( MiddleStart,
[ IsInfList ],
function( list )
    return list!.basePosition;
end );

InstallMethod( MiddleEnd,
[ IsInfList ],
function( list )
    return list!.basePosition + Length( list!.middle ) - 1;
end );

InstallMethod( MiddlePart,
[ IsInfList ],
function( list )
    return list!.middle;
end );

InstallMethod( PositivePart,
[ IsInfList ],
function( list )
    return list!.positive;
end );

InstallMethod( NegativePart,
[ IsInfList ],
function( list )
    return list!.negative;
end );

InstallMethod( LowestKnownPosition,
[ IsInfList ],
function( list )
    if NegativePart( list ) <> fail
       and LowestKnownPosition( NegativePart( list ) ) <> "none" then
        return LowestKnownPosition( NegativePart( list ) );
    elif Length( MiddlePart( list ) ) > 0 then
        return MiddleStart( list );
    else
        return LowestKnownPosition( PositivePart( list ) );
    fi;
end );

InstallMethod( HighestKnownPosition,
[ IsInfList ],
function( list )
    if PositivePart( list ) <> fail
       and HighestKnownPosition( PositivePart( list ) ) <> "none" then
        return HighestKnownPosition( PositivePart( list ) );
    elif Length( MiddlePart( list ) ) > 0 then
        return MiddleEnd( list );
    else
        return HighestKnownPosition( NegativePart( list ) );
    fi;
end );

InstallMethod( UpperBound,
[ IsInfList ],
function( list )
    if PositivePart( list ) <> fail then
        return PositiveInfinity;
    elif Length( MiddlePart( list ) ) > 0 then
        return MiddleEnd( list );
    elif NegativePart( list ) <> fail then
        return StartPosition( NegativePart( list ) );
    else
        return NegativeInfinity;
    fi;
end );

InstallMethod( LowerBound,
[ IsInfList ],
function( list )
    if NegativePart( list ) <> fail then
        return NegativeInfinity;
    elif Length( MiddlePart( list ) ) > 0 then
        return MiddleStart( list );
    elif PositivePart( list ) <> fail then
        return StartPosition( PositivePart( list ) );
    else
        return PositiveInfinity;
    fi;
end );

InstallMethod( \^,
[ IsInfList, IsInt ],
function( list, pos )
    local positive, negative;
    positive := PositivePart( list );
    negative := NegativePart( list );
    if pos >= MiddleStart( list ) and pos <= MiddleEnd( list ) then
        return MiddlePart( list )[ pos - MiddleStart( list ) + 1 ];
    elif positive <> fail and pos >= StartPosition( positive ) then
        return positive^pos;
    elif negative <> fail and pos <= StartPosition( negative ) then
        return negative^pos;
    else
        Error( "position ", pos, " outside index range of inflist ", list, "\n" );
    fi;
end );

InstallMethod( Shift,
[ IsHalfInfList, IsInt ],
function( list, shift )

    if IsRepeating( list ) then
        return MakeHalfInfList( StartPosition( list ) - shift, Direction( list ),
                                [ "repeat", RepeatingList( list ) ], false );
    else
        return MakeHalfInfList( StartPosition( list ) - shift, Direction( list ),
                                [ "pos", i -> list^(i + shift) ], false );
    fi;

end );

InstallMethod( Shift,
[ IsInfList and IsInfListDefaultRep, IsInt ],
function( list, shift )
    return MakeInfListFromHalfInfLists
           ( MiddleStart( list ) - shift,
             MiddlePart( list ),
             Shift( PositivePart( list ), shift ),
             Shift( NegativePart( list ), shift ) );
    end );

InstallMethod( FinitePartAsList,
[ IsInfList, IsInt, IsInt ],
function( list, startPos, endPos )
    return List( [ startPos .. endPos ], i -> list^i );
end );

InstallMethod( Cut,
[ IsHalfInfList, IsInt ],
function( list, pos )
    local middle, half, middleLength, middlePositions, newStartPos;

    if Direction( list ) = 1 and pos < StartPosition( list )
       or Direction( list ) = -1 and pos > StartPosition( list ) then
        Error( "position outside range of list" );
    fi;
    if pos = StartPosition( list ) then
        middle := [];
        half := list;
    elif InfListType( list ) = "repeat" then
        middleLength := ((StartPosition( list ) - pos) * Direction( list ))
                        mod Length( RepeatingList( list ) );
        if Direction( list ) = 1 then
            middlePositions := [ pos .. pos + middleLength - 1 ];
        else
            middlePositions := [ pos - middleLength + 1 .. pos ];
        fi;
        middle := List( middlePositions, i -> list^i );
        newStartPos := pos + middleLength * Direction( list );
        half := Shift( list, StartPosition( list ) - newStartPos );
    else
        middle := [];
        half := MakeHalfInfList( pos, Direction( list ),
                                 [ "pos", i -> list^i ],
                                 false );
    fi;

    if Direction( list ) = 1 then
        return MakeInfListFromHalfInfLists( pos, middle, half, fail );
    else
        return MakeInfListFromHalfInfLists( pos + 1 - Length( middle ),
                                            middle, fail, half );
    fi;

end );

InstallMethod( PositivePartFrom,
[ IsInfList, IsInt ],
function( list, pos )
    local middle, positiveStart;

    positiveStart := StartPosition( PositivePart( list ) );
    if pos >= positiveStart then
        return Cut( PositivePart( list ), pos );
    else
        middle := FinitePartAsList( list, pos, positiveStart - 1 );
        return MakeInfListFromHalfInfLists( pos, middle, PositivePart( list ), fail );
    fi;

end );

InstallMethod( NegativePartFrom,
[ IsInfList, IsInt ],
function( list, pos )
    local middle, negativeStart;

    negativeStart := StartPosition( NegativePart( list ) );
    if pos <= negativeStart then
        return Cut( NegativePart( list ), pos );
    else
        middle := FinitePartAsList( list, negativeStart + 1, pos );
        return MakeInfListFromHalfInfLists( negativeStart + 1, middle, fail,
                                            NegativePart( list ) );
    fi;

end );

InstallMethod( Splice,
[ IsInfList, IsInfList, IsInt ],
function( positiveList, negativeList, joinPosition )
    local positiveCut, negativeCut, middle;

    positiveCut := PositivePartFrom( positiveList, joinPosition + 1 );
    negativeCut := NegativePartFrom( negativeList, joinPosition );
    middle := Concatenation( MiddlePart( negativeCut ), MiddlePart( positiveCut ) );

    return MakeInfListFromHalfInfLists( MiddleStart( negativeCut ),
                                        middle,
                                        PositivePart( positiveCut ),
                                        NegativePart( negativeCut ) );

end );

InstallGlobalFunction( InfConcatenation,
function( arg )
    local middle, basePosition, positive, negative;

    # TODO: error if arguments are of wrong type

    if Length( arg ) = 0 then
        return FiniteInfList( 0, [] );
    elif Length( arg ) = 1 then
        return arg[ 1 ];
    fi;

    middle := CallFuncList( Concatenation, List( Reversed( arg ), MiddlePart ) );
    basePosition := MiddleEnd( arg[ 1 ] ) - Length( middle ) + 1;
    positive := PositivePart( arg[ 1 ] );
    negative := NegativePart( arg[ Length( arg ) ] );
    if negative <> fail then
        # Want to have StartPosition( negative ) = basePosition - 1
        negative := Shift( negative, StartPosition( negative ) - basePosition + 1 );
    fi;

    return MakeInfListFromHalfInfLists( basePosition, middle, positive, negative );

end );

InstallMethod( InfList,
[ IsInfList, IsFunction ],
function( list, func )

    return MakeInfListFromHalfInfLists
           ( MiddleStart( list ),
             List( MiddlePart( list ), func ),
             HalfInfList( PositivePart( list ), func ),
             HalfInfList( NegativePart( list ), func ) );

end );

InstallMethod( HalfInfList,
[ IsHalfInfList, IsFunction ],
function( list, func )

    if IsRepeating( list ) then
        return MakeHalfInfList( StartPosition( list ), Direction( list ),
                                [ "repeat", List( RepeatingList( list ), func ) ],
                                false );
    else
        return MakeHalfInfList( StartPosition( list ), Direction( list ),
                                [ "pos", i -> func( list^i ) ],
                                false );
    fi;

end );

InstallValue( IntegersList, FunctionInfList( IdFunc ) );

