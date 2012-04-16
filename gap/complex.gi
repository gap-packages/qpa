InstallGlobalFunction( Cat,
function( identity, properties )
    return Objectify( NewType( NewFamily( "CatFamily", IsCat ),
                               IsCat and IsCatDefaultRep ),
                      rec( identity := identity,
                           properties := properties ) );
end );

InstallMethod( \.,
[ IsCat, IsPosInt ],
function( cat, propName )
    return cat!.properties.( NameRNam( propName ) );
end );

InstallMethod( \=,
[ IsCat, IsCat ],
function( cat1, cat2 )
    return cat1!.identity = cat2!.identity;
end );

InstallMethod( PrintObj,
[ IsCat ],
function( cat )
    Print( "<cat: ", cat.name, ">" );
end );

InstallMethod( ViewObj,
[ IsCat ],
function( cat )
    Print( cat );
end );

InstallMethod( CatOfRightAlgebraModules,
[ IsAlgebra ],
function( A )
    local isZeroObj, composeMaps, isExact, objStr;
    isZeroObj := function( M )
        return Dimension( M ) = 0;
    end;
    composeMaps := function( g, f )
        return f * g;
    end;
    isExact := function( g, f )
        return Dimension( Kernel( g ) ) = Dimension( Image( f ) );
    end;
    objStr := function( M )
        local dims;
        dims := JoinStringsWithSeparator( DimensionVector( M ), "," );
        return Concatenation( "(", dims, ")" );
    end;
    return Cat( [ "right modules", A ],
                rec( name := "right modules over algebra",
                     zeroObj := ZeroModule( A ),
                     isZeroObj := isZeroObj,
                     zeroMap := ZeroMapping,
                     isZeroMapping := IsZero,
                     composeMaps := composeMaps,
                     ker := Kernel,
                     im := Image,
                     isExact := isExact,
                     objStr := objStr ) );
end );



InstallGlobalFunction( FiniteComplex,
function( cat, basePosition, differentials )
    return Complex( cat, basePosition, differentials, "zero", "zero" );
end );

InstallGlobalFunction( ZeroComplex,
function( cat )
    local fam, C;
    fam := NewFamily( "ComplexesFamily", IsComplex );
    fam!.cat := cat;
    C := Objectify( NewType( fam, IsZeroComplex and IsComplexDefaultRep ),
                    rec( ) );
    SetCatOfComplex( C, cat );
    SetDifferentialsOfComplex( C, ConstantInfList( cat.zeroMap( cat.zeroObj, cat.zeroObj ) ) );
    return C;
end );

InstallGlobalFunction( SingleObjectComplex,
function( cat, obj, degree )
    return FiniteComplex( cat, degree,
                          [ cat.zeroMap( obj, cat.zeroObj ),
                            cat.zeroMap( cat.zeroObj, obj ) ] );
end );

InstallGlobalFunction( ShortExactSequence,
function( cat, f, g )
    local SES;
    SES := FiniteComplex( cat, 0, [ g, f ] );
    if not IsShortExactSequence( SES ) then
        Error( "not exact\n" );
    fi;
    return SES;
end );

InstallGlobalFunction( Complex,
function( cat, basePosition, middle, positive, negative )
    local checkDifferentials, checkDifferentialList, checkDifferentialListWithRepeat,
          positiveRepeat, negativeRepeat,
          fam, C, basePositionL, middleL, positiveL, negativeL,
          firstMiddleObj, lastMiddleObj, checkNewDifferential;

    checkDifferentials := function( topDegree, indices, lists, listNames )
        local degrees, diffs, diffNames;
        degrees := [ topDegree, topDegree - 1 ];
        diffs := List( [1,2], i -> lists[ i ][ indices[ i ] ] );
        diffNames := List( [1,2], i -> Concatenation(
          "differential ", String( degrees[ i ] ),
          " (element ", String( indices[ i ] ), " in ", String( listNames[ i ] ), " list)" ) );
        if Range( diffs[ 1 ] ) <> Source( diffs[ 2 ] ) then
            Error( "range of ", diffNames[ 1 ], " is not the same as source of ",
                   diffNames[ 2 ], ".\n" );
        fi;
        if not cat.isZeroMapping( cat.composeMaps( diffs[ 2 ], diffs[ 1 ] ) ) then
            Error( "non-zero composition of ", diffNames[ 2 ],
                   " and ", diffNames[ 1 ], ".\n" );
        fi;
    end;

    checkDifferentialList := function( list, startDegree, listName )
        local i;
        for i in [ 2 .. Length( list ) ] do
            checkDifferentials( startDegree + i - 1,
                                [ i, i - 1 ],
                                [ list, list ],
                                [ listName, listName ] );
        od;
    end;

    checkDifferentialListWithRepeat := function( list, startDegree, listName,
                                                 repeatDirection )
                                       checkDifferentialList( list, startDegree, listName );
        checkDifferentials( startDegree + Length( list ) * repeatDirection,
                            [ 1, Length( list ) ],
                            [ list, list ], [ listName, listName ] );
    end;

    checkDifferentialList( middle, basePosition, "middle" );
    if positive[ 1 ] = "repeat" then
        positiveRepeat := positive[ 2 ];
        checkDifferentialListWithRepeat( positiveRepeat, basePosition + Length( middle ),
                                         "positive", 1 );
        if Length( middle ) > 0 then
            checkDifferentials( basePosition + Length( middle ),
                                [ 1, Length( middle ) ],
                                [ positiveRepeat, middle ],
                                [ "positive", "middle" ] );
        fi;
    fi;
    if negative[ 1 ] = "repeat" then
        negativeRepeat := Reversed( negative[ 2 ] );
        checkDifferentialListWithRepeat( negativeRepeat, basePosition + Length( middle ),
                                         "negative", 1 );
        if Length( middle ) > 0 then
            checkDifferentials( basePosition,
                                [ 1, Length( negativeRepeat ) ],
                                [ middle, negativeRepeat ],
                                [ "middle", "negative" ] );
        fi;
    fi;
    if positive[ 1 ] = "repeat" and negative[ 1 ] = "repeat" and Length( middle ) = 0 then
        checkDifferentials( basePosition,
                            [ 1, Length( negativeRepeat ) ],
                            [ positiveRepeat, negativeRepeat ],
                            [ "positive", "negative" ] );
    fi;

    # Create the complex object

    fam := NewFamily( "family of complexes", IsComplex );
    fam!.cat := cat;

    C := Objectify( NewType( fam, IsComplex and IsComplexDefaultRep ),
                    rec( ) );
    SetCatOfComplex( C, cat );

    # Normalize middle list if positive or negative is zero

    basePositionL := basePosition;
    middleL := middle;
    if positive = "zero" or negative = "zero" then
        firstMiddleObj := Range( middle[ 1 ] );
        lastMiddleObj := Source( middle[ Length( middle ) ] );
    fi;
    if positive = "zero" then
        # add zero object at the end if necessary:
        if not cat.isZeroObj( lastMiddleObj ) then
            middleL := Concatenation( middleL,
                                      [ cat.zeroMap( cat.zeroObj, lastMiddleObj ) ] );
        fi;
        # cut away superfluous zero objects:
        while cat.isZeroObj( Range( middleL[ Length( middleL ) ] ) ) do
            middleL := middleL{ [ 1 .. Length( middleL ) - 1 ] };
        od;
    fi;
    if negative = "zero" then
        # add zero object at the end if necessary:
        if not cat.isZeroObj( lastMiddleObj ) then
            middleL := Concatenation( [ cat.zeroMap( firstMiddleObj, cat.zeroObj ) ],
                                      middleL );
            basePositionL := basePositionL - 1;
        fi;
        # cut away superfluous zero objects:
        while cat.isZeroObj( Source( middleL[ 1 ] ) ) do
            middleL := middleL{ [ 2 .. Length( middleL ) ] };
            basePositionL := basePositionL + 1;
        od;
    fi;

    if positive = "zero" then
        positiveL := [ "repeat", [ cat.zeroMap( cat.zeroObj, cat.zeroObj ) ] ];
    elif positive[ 1 ] = "pos" then
        positiveL := ShallowCopy( positive );
        positiveL[ 2 ] := function( i )
            return positive[ 2 ]( C, i );
        end;
    else
        positiveL := positive;
    fi;
    if negative = "zero" then
        negativeL := [ "repeat", [ cat.zeroMap( cat.zeroObj, cat.zeroObj ) ] ];
    elif negative[ 1 ] = "pos" then
        negativeL := ShallowCopy( negative );
        negativeL[ 2 ] := function( i )
            return negative[ 2 ]( C, i );
        end;
    else
        negativeL := negative;
    fi;
            
    checkNewDifferential := function( i, dir, type )
        local degrees, diffs;
        if dir = 1 then
            degrees := [ i - 1, i ];
        else
            degrees := [ i, i + 1 ];
        fi;
        diffs := [ DifferentialOfComplex( C, degrees[ 1 ] ),
                   DifferentialOfComplex( C, degrees[ 2 ] ) ];
        if Source( diffs[ 1 ] ) <> Range( diffs[ 2 ] ) then
            Error( "source of differential ", degrees[ 1 ],
                   " is not the same as range of differential ",
                   degrees[ 2 ], " in complex\n   ", C, "\n" );
        fi;
        if not cat.isZeroMapping( cat.composeMaps( diffs[ 1 ], diffs[ 2 ] ) ) then
            Error( "nonzero composition of differentials ", degrees[ 1 ],
                   " and ", degrees[ 2 ], " in complex\n   ", C, "\n" );
        fi;
    end;

    SetDifferentialsOfComplex( C,
                               MakeInfList( basePositionL, middleL, positiveL, negativeL,
                                            checkNewDifferential ) );

    return C;

end );

InstallMethod( DifferentialOfComplex,
[ IsComplex, IsInt ],
function( C, i )
    return DifferentialsOfComplex( C )^i;
end );

InstallMethod( ObjectOfComplex,
[ IsComplex, IsInt ],
function( C, i )
    return Source( DifferentialOfComplex( C, i ) );
end );

InstallMethod( \^,
[ IsComplex, IsInt ],
function( C, i )
    return DifferentialOfComplex( C, i );
end );

InstallMethod( CyclesOfComplex,
[ IsComplex, IsInt ],
function( C, i )
    return Kernel( DifferentialOfComplex( C, i ) );
end );

InstallMethod( BoundariesOfComplex,
[ IsComplex, IsInt ],
function( C, i )
    return Image( DifferentialOfComplex( C, i + 1 ) );
end );

InstallMethod( HomologyOfComplex, # TODO: this does not work
[ IsComplex, IsInt ],
function( C, i )
    return CyclesOfComplex( C, i ) / BoundariesOfComplex( C, i );
end );

InstallMethod( IsFiniteComplex,
[ IsComplex ],
function( C )
    local upbound, lowbound;

    upbound := UpperBound( C );
    lowbound := LowerBound( C );

    if IsZeroComplex( C ) then
        return true;
    elif IsInt( upbound ) and IsInt( lowbound ) then
        return true;
    elif IsInfiniteNumber( upbound ) or IsInfiniteNumber( lowbound ) then
        return false;
    else
        return fail;
    fi;
end );

InstallMethod( LengthOfComplex,
[ IsComplex ],
function( C )
    local finiteness;
    finiteness := IsFiniteComplex( C );
    if IsZeroComplex( C ) then
        return 0;
    elif finiteness = true then
        return UpperBound( C ) - LowerBound( C ) + 1;
    elif finiteness = false then
        return PositiveInfinity;
    else
        return fail;
    fi;
end );

InstallMethod( HighestKnownDegree,
[ IsComplex ],
function( C )
    return HighestKnownPosition( DifferentialsOfComplex( C ) );
end );

InstallMethod( LowestKnownDegree,
[ IsComplex ],
function( C )
    return LowestKnownPosition( DifferentialsOfComplex( C ) );
end );

InstallMethod( IsExactSequence,
[ IsComplex ],
function( C )
    return ForEveryDegree( C, CatOfComplex( C ).isExact );
end );

InstallMethod( IsExactInDegree,
[ IsComplex, IsInt ],
function( C, i )
    return CatOfComplex( C ).isExact( DifferentialOfComplex( C, i ),
                                      DifferentialOfComplex( C, i + 1 ) );
end );

InstallMethod( IsShortExactSequence,
[ IsComplex ],
function( C )
    local length;
    length := LengthOfComplex( C );
    if length = fail then return fail; fi;
    if length = 3 then
        return IsExactSequence( C );
    else
        return false;
    fi;
end );

InstallMethod( ForEveryDegree, # TODO: misleading name?
[ IsComplex, IsFunction ],
function( C, func )
    local diffs, pos, neg, i;
    diffs := DifferentialsOfComplex( C );
    pos := PositivePart( diffs );
    neg := NegativePart( diffs );
    if not (IsRepeating( PositivePart( diffs ) ) and
            IsRepeating( NegativePart( diffs ) ) ) then
        return fail;
    fi;
    for i in [ StartPosition( neg ) - Length( RepeatingList( neg ) )
               .. StartPosition( pos ) + Length( RepeatingList( pos ) ) ] do
        if not func( diffs^i, diffs^(i+1) ) then
            return false;
        fi;
    od;
    return true;
end );

InstallMethod( UpperBound,
[ IsComplex ],
function( C )
    local cat, diffs, positive, i;

    cat := CatOfComplex( C );
    diffs := DifferentialsOfComplex( C );
    positive := PositivePart( diffs );

    if IsZeroComplex( C ) then
        return NegativeInfinity;
    elif IsRepeating( positive ) and Length( RepeatingList( positive ) ) = 1
      and cat.isZeroObj( ObjectOfComplex( C, StartPosition( positive ) ) ) then
        i := MiddleEnd( diffs ) - 1;
        while cat.isZeroObj( ObjectOfComplex( C, i ) ) do
            i := i - 1;
            if i < MiddleStart( diffs ) then
                return fail;
            fi;
        od;
        return i;
    elif IsRepeating( positive ) then
        return PositiveInfinity;
    else
        return fail;
    fi;
end );

InstallMethod( LowerBound,
[ IsComplex ],
function( C )
    local cat, diffs, negative, i;

    cat := CatOfComplex( C );
    diffs := DifferentialsOfComplex( C );
    negative := NegativePart( diffs );

    if IsZeroComplex( C ) then
        return PositiveInfinity;
    elif IsRepeating( negative ) and Length( RepeatingList( negative ) ) = 1
      and cat.isZeroObj( ObjectOfComplex( C, StartPosition( negative ) ) ) then
        i := MiddleStart( diffs );
        while cat.isZeroObj( ObjectOfComplex( C, i ) ) do
            i := i + 1;
            if i > MiddleEnd( diffs ) then
                return fail;
            fi;
        od;
        return i;
    elif IsRepeating( negative ) then
        return NegativeInfinity;
    else
        return fail;
    fi;
end );

InstallMethod( IsPositiveRepeating,
[ IsComplex ],
function( C )
    return IsRepeating( PositivePart( DifferentialsOfComplex( C ) ) );
end );

InstallMethod( IsNegativeRepeating,
[ IsComplex ],
function( C )
    return IsRepeating( NegativePart( DifferentialsOfComplex( C ) ) );
end );

InstallMethod( PositiveRepeatDegrees,
[ IsComplex ],
function( C )
    local positive, first, last;
    positive := PositivePart( DifferentialsOfComplex( C ) );
    if not IsRepeating( positive ) then
        return fail;
    fi;
    first := StartPosition( positive );
    last := first + Length( RepeatingList( positive ) ) - 1;
    return [ first .. last ];
end );

InstallMethod( NegativeRepeatDegrees,
[ IsComplex ],
function( C )
    local negative, first, last;
    negative := NegativePart( DifferentialsOfComplex( C ) );
    if not IsRepeating( negative ) then
        return fail;
    fi;
    first := StartPosition( negative );
    last := first - Length( RepeatingList( negative ) ) + 1;
    return [ last .. first ];
end );

InstallMethod( Shift,
[ IsComplex, IsInt ],
function( C, shift )
    local newDifferentials;

    newDifferentials := Shift( DifferentialsOfComplex( C ), shift );
    if shift mod 2 = 1 then
        newDifferentials := InfList( newDifferentials, d -> -d );
    fi;

    return ComplexByDifferentialList( CatOfComplex( C ), newDifferentials );

end );

InstallMethod( YonedaProduct,
[ IsComplex, IsComplex ],
function( C1, C2 )
    local cat, lowbound1, upbound2, diff1, diff2, connection, diffs;

    cat := CatOfComplex( C1 );

    lowbound1 := LowerBound( C1 );
    upbound2 := UpperBound( C2 );

    if not IsInt( lowbound1 ) then
        Error( "first argument in Yoneda product must be bounded below" );
    fi;
    if not IsInt( upbound2 ) then
        Error( "second argument in Yoneda product must be bounded above" );
    fi;
    if ObjectOfComplex( C1, lowbound1 ) <> ObjectOfComplex( C2, upbound2 ) then
        Error( "non-compatible complexes for Yoneda product" );
    fi;

    diff1 := Shift( DifferentialsOfComplex( C1 ), lowbound1 - upbound2 + 1 );
    diff2 := DifferentialsOfComplex( C2 );
    connection := cat.composeMaps( C2^upbound2, C1^(lowbound1+1) );

    diffs := InfConcatenation( PositivePartFrom( diff1, upbound2 + 1 ),
                               FiniteInfList( 0, [ connection ] ),
                               NegativePartFrom( diff2, upbound2 - 1 ) );

    return ComplexByDifferentialList( cat, diffs );

end );

InstallGlobalFunction( ComplexByDifferentialList,
[ IsCat, IsInfList ],
function( cat, differentials )
    local C, fam;

    fam := NewFamily( "ComplexesFamily", IsComplex );
    fam!.cat := cat;
    C := Objectify( NewType( fam, IsComplex and IsComplexDefaultRep ),
                    rec( ) );
    SetCatOfComplex( C, cat );
    SetDifferentialsOfComplex( C, differentials );
    return C;

end );

InstallMethod( PrintObj,
[ IsZeroComplex ],
function( C )
    Print( "0 -> 0" );
end );

InstallMethod( PrintObj,
[ IsComplex ],
function( C )
    local cat, diffs, i, upbound, lowbound, top, bottom;

    cat := CatOfComplex( C );

    diffs := DifferentialsOfComplex( C );
    upbound := UpperBound( C );
    lowbound := LowerBound( C );
    if IsInt( upbound ) then
        top := MiddleEnd( diffs ) - 1;
    elif IsPositiveRepeating( C ) then
        top := MiddleEnd( diffs );
    else
        top := HighestKnownDegree( C );
    fi;
    if IsNegativeRepeating( C ) then
        bottom := MiddleStart( diffs );
    else
        bottom := LowestKnownDegree( C );
    fi;
    
    if IsInt( upbound ) then
        Print( "0 -> " );
    else
        Print( "--- -> " );
    fi;

    if IsPositiveRepeating( C ) and upbound = PositiveInfinity then
        Print( "[ " );
        for i in Reversed( PositiveRepeatDegrees( C ) ) do
            Print( i, ":", cat.objStr( ObjectOfComplex( C, i ) ), " -> " );
        od;
        Print( "] " );
    fi;

    for i in [ top, top - 1 .. bottom ] do
        Print( i, ":", cat.objStr( ObjectOfComplex( C, i ) ), " -> " );
    od;

    if IsNegativeRepeating( C ) and lowbound = NegativeInfinity then
        Print( "[ " );
        for i in Reversed( NegativeRepeatDegrees( C ) ) do
            Print( i, ":", cat.objStr( ObjectOfComplex( C, i ) ), " -> " );
        od;
        Print( "] " );
    fi;

    if IsInt( lowbound ) then
        Print( "0" );
    else
        Print( "---" );
    fi;
end );


InstallMethod( ProjectiveResolution,
[ IsPathAlgebraMatModule ],
function( M )
    local nextDifferential, cover;

    nextDifferential := function( d )
        return ProjectiveCover( Kernel( d ) ) * KernelInclusion( d );
    end;

    cover := ProjectiveCover( M );

    return Complex( CatOfRightAlgebraModules( ActingAlgebra( M ) ),
                    0,
                    [ cover ],
                    [ "next", nextDifferential, cover ],
                    "zero" );
end );



InstallGlobalFunction( "ChainMap",
function( source, range, basePosition, middle, positive, negative )
    local cat, fam, map, positiveL, negativeL, numZeroMaps, i,
          correctDomainAt, correctCodomainAt, commutesAt,
          checkDomainAndCodomain, checkCommutes, checkNewMorphism,
          morphisms, firstCheckDegree, lastCheckDegree;

    cat := CatOfComplex( source );
    if cat <> CatOfComplex( range ) then
        Error( "source and range of chain map must be complexes over the same cat" );
    fi;

    fam := NewFamily( "ChainMapsFamily", IsChainMap );
    map := Objectify( NewType( fam, IsChainMap and IsChainMapDefaultRep ),
                      rec( ) );

    SetSource( map, source );
    SetRange( map, range );

    if positive = "zero" then
        positiveL := [ "pos",
                       i -> cat.zeroMap( ObjectOfComplex( source, i ),
                                         ObjectOfComplex( range, i ) ),
                       false ];
    elif positive[ 1 ] = "pos" then
        positiveL := ShallowCopy( positive );
        positiveL[ 2 ] := function( i )
            return positive[ 2 ]( map, i );
        end;
    else
        positiveL := positive;
    fi;
    if negative = "zero" then
        negativeL := [ "pos",
                       i -> cat.zeroMap( ObjectOfComplex( source, i ),
                                         ObjectOfComplex( range, i ) ),
                       false ];
    elif negative[ 1 ] = "pos" then
        negativeL := ShallowCopy( negative );
        negativeL[ 2 ] := function( i )
            return negative[ 2 ]( map, i );
        end;
    else
        negativeL := negative;
    fi;

    # if positive = "zero" then
    #     numZeroMaps := 0;
    #     for i in [ Length( middle ), Length( middle ) - 1 .. 1 ] do
    #         if cat.isZeroMapping( middle[ i ] ) then
    #             numZeroMaps := numZeroMaps + 1;
    #         else
    #             break;
    #         fi;
    #     od;
    #     SetUpperBound( map, basePosition + Length( middle ) - 1 - numZeroMaps );
    # fi;
    # if negative = "zero" then
    #     numZeroMaps := 0;
    #     for i in [ 1 .. Length( middle ) ] do
    #         if cat.isZeroMapping( middle[ i ] ) then
    #             numZeroMaps := numZeroMaps + 1;
    #         else
    #             break;
    #         fi;
    #     od;
    #     SetLowerBound( map, basePosition + numZeroMaps );
    # fi;

    correctDomainAt := function( i )
        return Source( map^i ) = ObjectOfComplex( source, i );
    end;
    correctCodomainAt := function( i )
        return Range( map^i ) = ObjectOfComplex( range, i );
    end;

    commutesAt := function( i )
        return cat.composeMaps( range^i, map^i ) = cat.composeMaps( map^(i-1), source^i );
    end;

    checkDomainAndCodomain := function( i )
        if not correctDomainAt( i ) then
            Error( "incorrect source of chain map morphism in degree ", i, "\n" );
        fi;
        if not correctCodomainAt( i ) then
            Error( "incorrect range of chain map morphism in degree ", i, "\n" );
        fi;
    end;

    checkCommutes := function( i )
        if not commutesAt( i ) then
            Error( "chain map morphisms at degrees (", i, ",", i-1, ") ",
                   "do not commute with complex differentials\n" );
        fi;
    end;

    checkNewMorphism := function( i, dir, type )
        local topDegreeOfNewSquare;
        if dir = 1 then
            topDegreeOfNewSquare := i;
        else
            topDegreeOfNewSquare := i + 1;
        fi;
        checkDomainAndCodomain( i );
        checkCommutes( topDegreeOfNewSquare );
    end;

    morphisms := MakeInfList( basePosition, middle, positiveL, negativeL,
                              checkNewMorphism );
    SetMorphismsOfChainMap( map, morphisms );

    # for i in [ LowestKnownDegree( map ) .. HighestKnownDegree( map ) ] do
    #     checkDomainAndCodomain( i );
    # od;
    # for i in [ LowestKnownDegree( map ) + 1 .. HighestKnownDegree( map ) ] do
    #     checkCommutes( i );
    # od;

    # TODO: hvor mye som må sjekkes (må se på kompleksene)
    if IsRepeating( NegativePart( morphisms ) ) then
        firstCheckDegree := StartPosition( NegativePart( morphisms ) )
                            - Length( RepeatingList( NegativePart( morphisms ) ) );
    else
        firstCheckDegree := MiddleStart( morphisms );
    fi;
    if IsRepeating( PositivePart( morphisms ) ) then
        lastCheckDegree := StartPosition( PositivePart( morphisms ) )
                           + Length( RepeatingList( PositivePart( morphisms ) ) );
    else
        lastCheckDegree := MiddleEnd( morphisms );
    fi;
    for i in [ firstCheckDegree .. lastCheckDegree ] do
        checkDomainAndCodomain( i );
    od;
    for i in [ firstCheckDegree + 1 .. lastCheckDegree ] do
        checkCommutes( i );
    od;

    return map;

end );

InstallMethod( HighestKnownDegree,
[ IsChainMap ],
function( map )
    return HighestKnownPosition( MorphismsOfChainMap( map ) );
end );

InstallMethod( LowestKnownDegree,
[ IsChainMap ],
function( map )
    return LowestKnownPosition( MorphismsOfChainMap( map ) );
end );

InstallMethod( MorphismOfChainMap,
[ IsChainMap, IsInt ],
function( map, i )
    return MorphismsOfChainMap( map )^i;
end );

InstallMethod( \^,
[ IsChainMap, IsInt ],
function( map, i )
    return MorphismsOfChainMap( map )^i;
end );

InstallMethod( PrintObj,
[ IsChainMap ],
function( map )
    Print( "<chain map>" );
end );

InstallGlobalFunction( ComplexAndChainMaps,
function( sourceComplexes, rangeComplexes,
          basePosition, middle, positive, negative )

    local cat, C, inMaps, outMaps, numInMaps, numOutMaps,
          positiveL, negativeL, list,
          positiveDiffs, positiveInMaps, positiveOutMaps,
          negativeDiffs, negativeInMaps, negativeOutMaps,
          middleDiffs, middleInMaps, middleOutMaps;

    cat := CatOfComplex( Concatenation( sourceComplexes, rangeComplexes )[ 1 ] );

    if positive = "zero" then
        positiveL := [ "repeat", [ fail ] ];
    elif positive[ 1 ] = "pos" then
        positiveL := ShallowCopy( positive );
        positiveL[ 2 ] := function( i )
            return CallFuncList( positive[ 2 ],
                                 Concatenation( [ C ], inMaps, outMaps, [ i ] ) );
        end;
    else
        positiveL := positive;
    fi;
    if negative = "zero" then
        negativeL := [ "repeat", [ fail ] ];
    elif negative[ 1 ] = "pos" then
        negativeL := ShallowCopy( negative );
        negativeL[ 2 ] := function( i )
            return CallFuncList( negative[ 2 ],
                                 Concatenation( [ C ], inMaps, outMaps, [ i ] ) );
        end;
    else
        negativeL := negative;
    fi;

    list := MakeInfList( basePosition, middle, positiveL, negativeL, false );

    numInMaps := Length( sourceComplexes );
    numOutMaps := Length( rangeComplexes );

    if positive = "zero" then
        positiveDiffs := "zero";
        positiveInMaps := List( [ 1 .. numInMaps ], i -> "zero" );
        positiveOutMaps := List( [ 1 .. numOutMaps ], i -> "zero" );
    else
        positiveDiffs := [ "pos",
                           function( C, i ) return \^(list, i)[ 1 ]; end ];
        positiveInMaps :=
          List( [ 1 .. numInMaps ],
                i -> [ "pos", function( map, pos ) return \^(list, pos)[ 1 + i ]; end ] );
        positiveOutMaps :=
          List( [ 1 .. numOutMaps ],
                i -> [ "pos", function( map, pos )
                                 return  \^(list, pos)[ 1 + numInMaps + i ]; end ] );
    fi;
    if negative = "zero" then
        negativeDiffs := "zero";
        negativeInMaps := List( [ 1 .. numInMaps ], i -> "zero" );
        negativeOutMaps := List( [ 1 .. numOutMaps ], i -> "zero" );
    else
        negativeDiffs := [ "pos", function( C, i ) return \^(list, i)[ 1 ]; end ];
        negativeInMaps :=
          List( [ 1 .. numInMaps ],
                i -> [ "pos", function( map, pos ) return \^(list, pos)[ 1 + i ]; end ] );
        negativeOutMaps :=
          List( [ 1 .. numOutMaps ],
                i -> [ "pos", function( map, pos )
                                 return \^(list, pos)[ 1 + numInMaps + i ]; end ] );
    fi;

    middleDiffs := List( middle, L -> L[ 1 ] );
    middleInMaps := List( [ 1 .. numInMaps ],
                          i -> List( middle, L -> L[ 1 + i ] ) );
    middleOutMaps := List( [ 1 .. numOutMaps ],
                           i -> List( middle, L -> L[ 1 + numInMaps + i ] ) );

    C := Complex( cat, basePosition, middleDiffs, positiveDiffs, negativeDiffs );
    inMaps := List( [ 1 .. numInMaps ],
                    i -> ChainMap( sourceComplexes[ i ],
                                   C,
                                   basePosition,
                                   middleInMaps[ i ],
                                   positiveInMaps[ i ],
                                   negativeInMaps[ i ] ) );
    outMaps := List( [ 1 .. numOutMaps ],
                     i -> ChainMap( C,
                                    rangeComplexes[ i ],
                                    basePosition,
                                    middleOutMaps[ i ],
                                    positiveOutMaps[ i ],
                                    negativeOutMaps[ i ] ) );

    return Concatenation( [ C ], inMaps, outMaps );
                       
end );

