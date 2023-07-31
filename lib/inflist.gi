#######################################################################
##
#F  MakeHalfInfList( <start>, <direction>, <typeWithArgs>, <callback>, <repeatifyCallback> )
##
##  Creates a IsHalfInfList with start index <start>, <direction> is
##  -1 for negative or 1 for positive, <typeWithArgs> is a list of either
##  two or three arguments describing the nature of the infinite list
##  (see the documentation), and <callback> is a function which is called
##  whenever a new list element is computed.
##  
InstallGlobalFunction( MakeHalfInfList,
function( start, direction, typeWithArgs, callback, repeatifyCallback )
    local type, repeatingList, func, initialValue, storeValues, data, list,
          repeatify;

    type := typeWithArgs[ 1 ];
    if type = "next/repeat" then
        type := "next";
        repeatify := true;
    else
        repeatify := false;
    fi;
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
                 callback := callback,
                 repeatify := repeatify,
                 repeatifyCallback := repeatifyCallback );

    list := Objectify( NewType( NewFamily( "HalfInfListsFamily" ),
                                IsHalfInfList and IsHalfInfListDefaultRep ),
                       data );

    return list;

end );

#######################################################################
##
#M  StartPosition( <list> )
##
##  Returns the start position of a IsHalfInfList <list>.
##  
InstallMethod( StartPosition,
[ IsHalfInfList ],
function( list )
    return list!.start;
end );

#######################################################################
##
#M  Direction( <list> )
##
##  Returns the direction of a IsHalfInfList <list>.
##  
InstallMethod( Direction,
[ IsHalfInfList ],
function( list )
    return list!.direction;
end );

#######################################################################
##
#M  InfListType( <list> )
##
##  Returns the type of a IsHalfInfList <list>.
##  
InstallMethod( InfListType,
[ IsHalfInfList ],
function( list )
    return list!.type;
end );

#######################################################################
##
#M  RepeatingList( <list> )
##
##  Returns the repeatingList of a IsHalfInfList <list>.
##  
InstallMethod( RepeatingList,
[ IsHalfInfList ],
function( list )
    return list!.repeatingList;
end );

#######################################################################
##
#M  ElementFunction( <list> )
##
##  Returns the element-function of a IsHalfInfList <list>.
##  
InstallMethod( ElementFunction,
[ IsHalfInfList ],
function( list )
    return list!.func;
end );

#######################################################################
##
#M  IsStoringValues( <list> )
##
##  Returns the value of the (boolean) storingValues of a 
##  IsHalfInfList <list>.
##  
InstallMethod( IsStoringValues,
[ IsHalfInfList ],
function( list )
    return list!.storingValues;
end );

#######################################################################
##
#M  NewValueCallback( <list> )
##
##  Returns the callback function of a IsHalfInfList <list>.
##  
InstallMethod( NewValueCallback,
[ IsHalfInfList ],
function( list )
    return list!.callback;
end );

#######################################################################
##
#M  IsRepeating( <list> )
##
##  Returns true if the type of the IsHalfInfList <list> is "repeat",
##  false otherwise.
##  
InstallMethod( IsRepeating,
[ IsHalfInfList ],
function( list )
    return list!.type = "repeat";
end );

#######################################################################
##
#M  InitialValue( <list> )
##
##  Returns the initial value of a IsHalfInfList <list>.
##  
InstallMethod( InitialValue,
[ IsHalfInfList ],
function( list )
    return list!.initialValue;
end );

#######################################################################
##
#M  \^( <list>, <pos> )
##
##  Returns the stored value at index <pos> in the IsHalfInfList
##  <list>.
##  
InstallMethod( \^,
[ IsHalfInfList and IsHalfInfListDefaultRep, IsInt ],
function( list, pos )
    local index, i, callCallbackFunction, r;

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
            if list!.repeatify and i mod 2 = 0 then
                if list!.values[ i / 2 + 1 ] = list!.values[ i + 1 ] then
                    r := MakeRepeatingList( list, i / 2 + 1, i + 1 );
                    list!.repeatifyCallback( r[ 1 ], r[ 2 ] );
                    list!.repeatify := false;
                    return r[ 2 ]^pos;
                fi;
            fi;
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

InstallMethod( MakeRepeatingList,
               [ IsHalfInfList, IsPosInt, IsPosInt ],
function( list, collisionIndex1, collisionIndex2 )
  local values, i1, i2, i, repeatStartIndex, repeatEndIndex,
        tail, repeatedList, newList;
  values := list!.values;
  i1 := collisionIndex1;
  i2 := collisionIndex2;
  while i1 > 0 and values[ i1 ] = values[ i2 ] do
    i1 := i1 - 1;
    i2 := i2 - 1;
  od;
  repeatStartIndex := i1 + 1;
  i := repeatStartIndex + 1;
  while values[ i ] <> values[ repeatStartIndex ] do
    i := i + 1;
  od;
  repeatEndIndex := i - 1;
  tail := values{ [ 1 .. ( repeatStartIndex - 1 ) ] };
  if Direction( list ) = -1 then
    tail := Reversed( tail );
  fi;
  repeatedList := values{ [ repeatStartIndex .. repeatEndIndex ] };
  newList := MakeHalfInfList( StartPosition( list ) + ( repeatStartIndex - 1 ) * Direction( list ),
                              Direction( list ),
                              [ "repeat", repeatedList ],
                              false, false );
  return [ tail, newList ];
end );

#######################################################################
##
#M  LowestKnownPosition( <list> )
##
##  Returns the lowest index of the list where the value is known
##  without computation.
##  
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
          return -infinity;
        elif Length( list!.values ) = 0 and InfListType( list ) <> "repeat" then
            return "none";
        else
            return StartPosition( list ) - (Length( list!.values ) - 1);
        fi;
    fi;
end );

#######################################################################
##
#M  HighestKnownPosition( <list> )
##
##  Returns the highest index of the list where the value is known
##  without computation.
##  
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
          return infinity;
        elif Length( list!.values ) = 0 and InfListType( list ) <> "repeat" then
            return "none";
        else
            return StartPosition( list ) + (Length( list!.values ) - 1);
        fi;
    fi;
end );


#######################################################################
##
#F  MakeInfList( <basePosition>, <middle>, <positive>, <negative>,
##               <callback> )
##
##  Creates an IsInfList object.  <basePostition> tells in which index
##  the values should be put, <middle> is a list of known values,
##  <positive> is a list describing the positive part of the infinite
##  list, and <negative> part is a list describing the negative part
##  of the list.  <callback> is a function to be called whenever a 
##  new value of the infinite list is computed.
##  
InstallGlobalFunction( MakeInfList,
function( basePosition, middle, positive, negative, callback )
    local list, posList, negList, posRepeatify, negRepeatify;

    # TODO: automatic initialValue for "next" parts if middle is nonempty

    # TODO: allow empty positive/negative?

    posRepeatify := function( tail, repeating )
      list!.middle := Concatenation( list!.middle, tail );
      list!.positive := repeating;
    end;
    posList := MakeHalfInfList( basePosition + Length( middle ),
                                1, positive, callback, posRepeatify );

    negRepeatify := function( tail, repeating )
      list!.middle := Concatenation( tail, list!.middle );
      list!.negative := repeating;
    end;
    negList := MakeHalfInfList( basePosition - 1,
                                -1, negative, callback, negRepeatify );

    list := MakeInfListFromHalfInfLists( basePosition, middle, posList, negList );
    return list;

end );

#######################################################################
##
#F  MakeInfListFromHalfInfLists( <basePosition>, <middle>, <positive>,
##                                <negative> )
##  
##  Creates an IsInfList from a middle part (a list) and two IsHalfInfLists
##  <positive> and <negative>.
##
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

#######################################################################
##
#F  FunctionInfList( <func> )
##  
##  Creates an IsInfList where all list entries are described by a 
##  function.
##
InstallGlobalFunction( FunctionInfList,
function( func )
    local positiveNegativeList;
    positiveNegativeList := [ "pos", func, false ];
    return MakeInfList( 0, [], positiveNegativeList, positiveNegativeList, false );
end );

#######################################################################
##
#F  ConstantInfList( <value> )
##  
##  Creates an IsInfList with <value> in each position.
##
InstallGlobalFunction( ConstantInfList,
function( value )
    return MakeInfList( 0, [], [ "repeat", [ value ] ], [ "repeat", [ value ] ], false );
end );

#######################################################################
##
#F  FiniteInfList( <basePosition>, <list> )
##  
##  Creates an IsInfList which is finite.  Only indexes in the interval
##  [basePostition, basePosition + length(<list>) - 1] is allowed.
##
InstallGlobalFunction( FiniteInfList,
function( basePosition, list )
    return MakeInfListFromHalfInfLists( basePosition, list, fail, fail );
end );

#######################################################################
##
#M  MiddleStart( <list> )
##  
##  Returns the first index of the middle part, or the basePosition,
##  of the IsInfList <list>.
##
InstallMethod( MiddleStart,
[ IsInfList ],
function( list )
    return list!.basePosition;
end );

#######################################################################
##
#M  MiddleEnd( <list> )
##  
##  Returns the last index of the middle part of the IsInfList <list>.
##
InstallMethod( MiddleEnd,
[ IsInfList ],
function( list )
    return list!.basePosition + Length( list!.middle ) - 1;
end );

#######################################################################
##
#M  MiddlePart( <list> )
##  
##  Returns the middle part of the IsInfList <list>, as a list.
##
InstallMethod( MiddlePart,
[ IsInfList ],
function( list )
    return list!.middle;
end );

#######################################################################
##
#M  PositivePart( <list> )
##  
##  Returns the positive part of the IsInfList <list>, as an IsHalfInfList.
##
InstallMethod( PositivePart,
[ IsInfList ],
function( list )
    return list!.positive;
end );

#######################################################################
##
#M  NegativePart( <list> )
##  
##  Returns the negative part of the IsInfList <list>, as an IsHalfInfList.
##
InstallMethod( NegativePart,
[ IsInfList ],
function( list )
    return list!.negative;
end );

#######################################################################
##
#M  LowestKnownPosition( <list> )
##  
##  Returns the lowest index where the value of the IsInfList <list>
##  is computed or otherwise known. 
##
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

#######################################################################
##
#M  HighestKnownPosition( <list> )
##  
##  Returns the highest index where the value of the IsInfList <list>
##  is computed or otherwise known. 
##
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

#######################################################################
##
#M  UpperBound( <list> )
##  
##  Returns the highest index where the value of the IsInfList <list>
##  exists (but is not necessarily known).
##
InstallMethod( UpperBound,
[ IsInfList ],
function( list )
    if PositivePart( list ) <> fail then
      return infinity;
    elif Length( MiddlePart( list ) ) > 0 then
        return MiddleEnd( list );
    elif NegativePart( list ) <> fail then
        return StartPosition( NegativePart( list ) );
    else
      return -infinity;
    fi;
end );

#######################################################################
##
#M  LowerBound( <list> )
##  
##  Returns the smallest index where the value of the IsInfList <list>
##  exists (but is not necessarily known).
##
InstallMethod( LowerBound,
[ IsInfList ],
function( list )
    if NegativePart( list ) <> fail then
      return -infinity;
    elif Length( MiddlePart( list ) ) > 0 then
        return MiddleStart( list );
    elif PositivePart( list ) <> fail then
        return StartPosition( PositivePart( list ) );
    else
      return infinity;
    fi;
end );

#######################################################################
##
#M  \^( <list>, <pos> )
##  
##  Returns the value of the IsInfList <list> at the index <pos>.
##
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

#######################################################################
##
#M  Shift( <list>, <shift> )
##  
##  <list> is an IsHalfInfList.  The method shifts the list <shift> 
##  positions; to the right if <shift> is positive or to the left if
##  <shift> is negative.
##
InstallMethod( Shift,
[ IsHalfInfList, IsInt ],
function( list, shift )

    if IsRepeating( list ) then
        return MakeHalfInfList( StartPosition( list ) - shift, Direction( list ),
                                [ "repeat", RepeatingList( list ) ], false, false );
    else
        return MakeHalfInfList( StartPosition( list ) - shift, Direction( list ),
                                [ "pos", i -> list^(i + shift) ], false, false );
    fi;

end );

#######################################################################
##
#M  Shift( <list>, <shift> )
##  
##  <list> is an IsInfList.  The method shifts the list <shift> 
##  positions; to the right if <shift> is positive or to the left if
##  <shift> is negative.
##
InstallMethod( Shift,
[ IsInfList and IsInfListDefaultRep, IsInt ],
function( list, shift )
    return MakeInfListFromHalfInfLists
           ( MiddleStart( list ) - shift,
             MiddlePart( list ),
             Shift( PositivePart( list ), shift ),
             Shift( NegativePart( list ), shift ) );
    end );

#######################################################################
##
#M  FinitePartAsList( <list>, <startPos>, <endPos> )
##  
##  <list> is an IsInfList.  The part of the list starting at index
##  <startPos> and ending at index <endPos> is returned as a list.
##  Note that <endPos> > <startPos>.
##
InstallMethod( FinitePartAsList,
[ IsInfList, IsInt, IsInt ],
function( list, startPos, endPos )
    return List( [ startPos .. endPos ], i -> list^i );
end );

#######################################################################
##
#M  Cut( <list>, <pos> )
##
##  Returns a new IsHalfInfList which the old IsHalfInfList except
##  that the positions with indices smaller than <pos> is removed.
##  
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
                                 false, false );
    fi;

    if Direction( list ) = 1 then
        return MakeInfListFromHalfInfLists( pos, middle, half, fail );
    else
        return MakeInfListFromHalfInfLists( pos + 1 - Length( middle ),
                                            middle, fail, half );
    fi;

end );

#######################################################################
##
#M  PositivePartFrom( <list>, <pos> )
##
##  Returns a new IsInfList with only the positions with indices 
##  greater than or equal to <pos>.
##  
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

#######################################################################
##
#M  NegativePartFrom( <list>, <pos> )
##
##  Returns a new IsInfList with only the positions with indices 
##  smaller than or equal to <pos>.
##  
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

#######################################################################
##
#M  Splice( <positiveList>, <negativeList>, <joinPosition> )
##
##  Returns a new IsInfList which is identical to <positiveList> for
##  indices greater than <joinPosition> and identical to <negativeList>
##  for indices smaller than or equal to <joinPosition>.
##  
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

#######################################################################
##
#F  InfConcatenation( <arg> )
##
##  <arg> is a list of IsInfLists. The method creates a new IsInfList
##  where the middle part is the middle parts of the respective lists
##  of <arg>, the positive part is the positive part of the first list
##  and the negative part is the negative part of the last list.
##  
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

#######################################################################
##
#M  InfList( <list>, <func> )
##
##  <list> is an InfList and <func> is a function which can take 
##  elements of <list> as argument.  The method returns a new list
##  where the elements are the elements of <list> with <func> applied
##  to them.
##  
InstallMethod( InfList,
[ IsInfList, IsFunction ],
function( list, func )

    return MakeInfListFromHalfInfLists
           ( MiddleStart( list ),
             List( MiddlePart( list ), func ),
             HalfInfList( PositivePart( list ), func ),
             HalfInfList( NegativePart( list ), func ) );

end );

#######################################################################
##
#M  HalfInfList( <list>, <func> )
##
##  <list> is a HalfInfList and <func> is a function which can take 
##  elements of <list> as argument.  The method returns a new list
##  where the elements are the elements of <list> with <func> applied
##  to them.
##  
InstallMethod( HalfInfList,
[ IsHalfInfList, IsFunction ],
function( list, func )

    if IsRepeating( list ) then
        return MakeHalfInfList( StartPosition( list ), Direction( list ),
                                [ "repeat", List( RepeatingList( list ), func ) ],
                                false, false );
    else
        return MakeHalfInfList( StartPosition( list ), Direction( list ),
                                [ "pos", i -> func( list^i ) ],
                                false, false );
    fi;

end );

#######################################################################
##
#V  IntegersList
##
##  An IsInfList with the integer i at position i.
##
BindGlobal( "IntegersList", FunctionInfList( IdFunc ) );
