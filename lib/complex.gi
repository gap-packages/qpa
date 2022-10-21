#######################################################################
##
#F  Cat( <identity>, <properties> ) 
##
##  Constructs a category object from two pieces of data:
##  -- <identity>, a list with information identifying the category, 
##     f.ex. a name and some object the cat is constructed from.
##  -- <properties>, a record with the following entries:
##     name (will be used for printing a categeory), zeroObj (the zero
##     object of the cat), isZeroObj (a test function), zeroMap (an 
##     operation on two objects of the cat), isZeroMapping (a test
##     function), composeMaps (an operation on two maps of the cat), 
##     ker (an operation on the maps of the cat), im (an operation on
##     the maps of the cat), isExact (an operation on two maps of the 
##     cat) and objStr (a function printing the objects of the cat).
##  
InstallGlobalFunction( Cat,
function( identity, properties )
    return Objectify( NewType( NewFamily( "CatFamily", IsCat ),
                               IsCat and IsCatDefaultRep ),
                      rec( identity := identity,
                           properties := properties ) );
end );

#######################################################################
##
#M  \.( <cat>, <propName> ) 
##
##  Method for retrieving properties from a category, where <propName>
##  is one of the properties listed in the Cat() documentation above.
##  
InstallMethod( \.,
[ IsCat, IsPosInt ],
function( cat, propName )
    return cat!.properties.( NameRNam( propName ) );
end );

#######################################################################
##
#M  \=( <cat1>, <cat2> ) 
##
##  Equality test for two categories. Two categories are equal if
##  and only if their identity are the same.
##  
InstallMethod( \=,
[ IsCat, IsCat ],
function( cat1, cat2 )
    return cat1!.identity = cat2!.identity;
end );

#######################################################################
##
#M  PrintObj( <cat> ) 
##
##  Printing a category, using its name.
##  
InstallMethod( PrintObj,
[ IsCat ],
function( cat )
    Print( "<cat: ", cat.name, ">" );
end );

#######################################################################
##
#M  ViewObj( <cat> ) 
##
##  Viewing a cateory, adding no more information than the PrintObj.
##  
InstallMethod( ViewObj,
[ IsCat ],
function( cat )
    Print( cat );
end );

#######################################################################
##
#M  CatOfRightAlgebraModules( <A> ) 
##
##  Returns the category of right modules over some (quotient of a) 
##  path algebra <A>.
##  
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


#######################################################################
##
#F  FiniteComplex( <cat>, <basePosition>, <differentials> ) 
##
##  This function returns a finite complex with objects in <cat>.  The
##  differentials are given in the list <differentials> = [d1, ..., dN],
##  an <basePosition> is some integer i. The returned complex has the 
##  map d1 from degree (i) to degree (i-1).
##  
InstallGlobalFunction( FiniteComplex,
function( cat, basePosition, differentials )
    return Complex( cat, basePosition, differentials, "zero", "zero" );
end );

#######################################################################
##
#F  ZeroComplex( <cat> ) 
##
##  Returns the complex in which all objects are the zero object in
##  <cat>.
##  
InstallGlobalFunction( ZeroComplex,
function( cat )
    local fam, C;
    fam := NewFamily( "ComplexesFamily", IsQPAComplex );
    fam!.cat := cat;
    C := Objectify( NewType( fam, IsZeroComplex and IsQPAComplexDefaultRep ),
                    rec( ) );
    SetCatOfComplex( C, cat );
    SetDifferentialsOfComplex( C, ConstantInfList( cat.zeroMap( cat.zeroObj, cat.zeroObj ) ) );
    return C;
end );

#######################################################################
##
#F  StalkComplex( <cat>, <obj>, <degree> ) 
##
##  Returns the stalk complex with the object <obj> from <cat> in
##  degree <degree>.
##  
InstallGlobalFunction( StalkComplex,
function( cat, obj, degree )
    return FiniteComplex( cat, degree,
                          [ cat.zeroMap( obj, cat.zeroObj ),
                            cat.zeroMap( cat.zeroObj, obj ) ] );
end );

#######################################################################
##
#F  ShortExactSequence( <cat>, <f>, <g> ) 
##
##  Returns a complex with three non-zero consecutive objects, and zero
##  objects elsewhere, such that the complex is exact: The image of <f>
##  should be the kernel of <g>, and <f> should be injective, and <g>
##  should be surjective. The function checks that this is the case, 
##  and returns an error otherwise.
##  
InstallGlobalFunction( ShortExactSequence,
function( cat, f, g )
    local SES;
    SES := FiniteComplex( cat, 0, [ g, f ] );
    if not IsShortExactSequence( SES ) then
        Error( "not exact\n" );
    fi;
    return SES;
end );

#######################################################################
##
#F  Complex( <cat>, <basePosition>, <middle>, <positive>, <negative> ) 
##
##  Constructs a complex, not necessarily finite, from the given data.
##  See the QPA manual for detailed information on the input data.
##  
InstallGlobalFunction( Complex,
function( cat, basePosition, middle, positive, negative )
    local checkDifferentials, checkDifferentialList, checkDifferentialListWithRepeat,
          positiveRepeat, negativeRepeat,
          fam, C, basePositionL, middleL, positiveL, negativeL,
          firstMiddleObj, lastMiddleObj, checkNewDifferential,
          firstMap, next, secondMap;

    # check that all consecutive differentials compose to zero
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

    fam := NewFamily( "family of complexes", IsQPAComplex );
    fam!.cat := cat;

    C := Objectify( NewType( fam, IsQPAComplex and IsQPAComplexDefaultRep ),
                    rec( ) );
    SetCatOfComplex( C, cat );

    # Normalize middle list if positive or negative is zero

    basePositionL := basePosition;
    middleL := middle;
    if positive = "zero" then
        lastMiddleObj := Source( middle[ Length( middle ) ] );
        # add zero object at the end if necessary:
        if not cat.isZeroObj( lastMiddleObj ) then
            middleL := Concatenation( middleL,
                                      [ cat.zeroMap( cat.zeroObj, lastMiddleObj ) ] );
        fi;
        # cut away superfluous zero objects:
        while Length( middleL ) > 0 and cat.isZeroObj( Range( middleL[ Length( middleL ) ] ) ) do
            middleL := middleL{ [ 1 .. Length( middleL ) - 1 ] };
        od;
    fi;
    if negative = "zero" then
        firstMiddleObj := Range( middle[ 1 ] );
        # add zero object at the end if necessary:
        if not cat.isZeroObj( firstMiddleObj ) then
            middleL := Concatenation( [ cat.zeroMap( firstMiddleObj, cat.zeroObj ) ],
                                      middleL );
            basePositionL := basePositionL - 1;
        fi;
        # cut away superfluous zero objects:
        while Length( middleL ) > 0 and cat.isZeroObj( Source( middleL[ 1 ] ) ) do
            middleL := middleL{ [ 2 .. Length( middleL ) ] };
            basePositionL := basePositionL + 1;
        od;
    fi;

    if positive = "zero" and Length( middleL ) = 0 and negative[ 1 ] = "next/repeat" then
      firstMap := negative[ 3 ];
      next := negative[ 2 ];
      secondMap := next( firstMap );
      if cat.isZeroMapping( firstMap ) and cat.isZeroMapping( secondMap ) then
        return ZeroComplex( cat );
      fi;
    fi;
    if negative = "zero" and Length( middleL ) = 0 and positive[ 1 ] = "next/repeat" then
      firstMap := positive[ 3 ];
      next := positive[ 2 ];
      secondMap := next( firstMap );
      if cat.isZeroMapping( firstMap ) and cat.isZeroMapping( secondMap ) then
        return ZeroComplex( cat );
      fi;
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

#######################################################################
##
#M  DifferentialOfComplex( <C>, <i> ) 
##
##  Returns the differential in degree <i> of the complex <C>.
##  
InstallMethod( DifferentialOfComplex,
[ IsQPAComplex, IsInt ],
function( C, i )
    return DifferentialsOfComplex( C )^i;
end );

#######################################################################
##
#M  ObjectOfComplex( <C>, <i> ) 
##
##  Returns the object in degree <i> of the complex <C>.
##  
InstallMethod( ObjectOfComplex,
[ IsQPAComplex, IsInt ],
function( C, i )
    return Source( DifferentialOfComplex( C, i ) );
end );

#######################################################################
##
#M  \^( <C>, <i> ) 
##
##  Is this in use??
##  
InstallMethod( \^,
[ IsQPAComplex, IsInt ],
function( C, i )
    return DifferentialOfComplex( C, i );
end );

#######################################################################
##
#M  CyclesOfComplex( <C>, <i> ) 
##
##  For a complex <C> and an integer <i>. Returns the i-cycle of the
##  complex, that is the subobject Ker(d_i) of the object in degree i.
##  
InstallMethod( CyclesOfComplex,
[ IsQPAComplex, IsInt ],
function( C, i )
    return Kernel( DifferentialOfComplex( C, i ) );
end );

#######################################################################
##
#M  BoundariesOfComplex( <C>, <i> ) 
##
##  For a complex <C> and an integer <i>. Returns the i-boundary of the
##  complex, that is the subobject Im(d_{i+1}) of the object in degree i.
##  
InstallMethod( BoundariesOfComplex,
[ IsQPAComplex, IsInt ],
function( C, i )
    return Image( DifferentialOfComplex( C, i + 1 ) );
end );

#######################################################################
##
#M  HomologyOfComplex( <C>, <i> ) 
##
##  For a complex <C> and an integer <i>. Returns the ith homology of 
##  the complex.
##
##  TODO: Does not currently work (see the documentation).
##  
InstallMethod( HomologyOfComplex, # TODO: this does not work
[ IsQPAComplex, IsInt ],
function( C, i )
    return CyclesOfComplex( C, i ) / BoundariesOfComplex( C, i );
end );

#######################################################################
##
#M  IsFiniteComplex( <C> ) 
##
##  Returns true if the complex <C> is a finite complex, false otherwise.
##  
InstallMethod( IsFiniteComplex,
[ IsQPAComplex ],
function( C )
    local upbound, lowbound;

    upbound := UpperBound( C );
    lowbound := LowerBound( C );

    if IsZeroComplex( C ) then
        return true;
    elif IsInt( upbound ) and IsInt( lowbound ) then
        return true;
      elif ( upbound in [ -infinity, infinity ] ) or ( lowbound in [ -infinity, infinity ] ) then
        return false;
    else
        return fail;
    fi;
end );

#######################################################################
##
#M  LengthOfComplex( <C> ) 
##
##  Returns the length of the complex <C>. If C is a zero complex, then
##  the length is zero. If C is a finite complex, the length is the
##  upper bound - the lower bound + 1. If C is an infinite complex, the
##  length is infinity.
##  
InstallMethod( LengthOfComplex,
[ IsQPAComplex ],
function( C )
    local finiteness;
    finiteness := IsFiniteComplex( C );
    if IsZeroComplex( C ) then
        return 0;
    elif finiteness = true then
        return UpperBound( C ) - LowerBound( C ) + 1;
    elif finiteness = false then
      return infinity;
    else
        return fail;
    fi;
end );

#######################################################################
##
#M  HighestKnownDegree( <C> ) 
##
##  Returns the greatest integer i such that the object at position i 
##  is known (or computed). For a finite complex, this will be infinity.
##  
InstallMethod( HighestKnownDegree,
[ IsQPAComplex ],
function( C )
    return HighestKnownPosition( DifferentialsOfComplex( C ) );
end );

#######################################################################
##
#M  LowestKnownDegree( <C> ) 
##
##  Returns the smallest integer i such that the object at position i 
##  is known (or computed). For a finite complex, this will be negative
##  infinity.
##  
InstallMethod( LowestKnownDegree,
[ IsQPAComplex ],
function( C )
    return LowestKnownPosition( DifferentialsOfComplex( C ) );
end );

#######################################################################
##
#M  IsExactSequence( <C> ) 
##
##  True if the complex <C> is exact in every degree. If the complex
##  is not finite and not repeating, the function fails.
##  
InstallMethod( IsExactSequence,
[ IsQPAComplex ],
function( C )
    return ForEveryDegree( C, CatOfComplex( C ).isExact );
end );

#######################################################################
##
#M  IsExactInDegree( <C>, <i> ) 
##
##  Returns true if the complex <C> is exact in degree <i>.
##  
InstallMethod( IsExactInDegree,
[ IsQPAComplex, IsInt ],
function( C, i )
    return CatOfComplex( C ).isExact( DifferentialOfComplex( C, i ),
                                      DifferentialOfComplex( C, i + 1 ) );
end );

#######################################################################
##
#M  IsShortExactSequence( <C> ) 
##
##  Returns true if the complex <C> is exact and has only three non-zero
##  objects, which are consecutive.
##  
InstallMethod( IsShortExactSequence,
[ IsQPAComplex ],
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

#######################################################################
##
#M  ForEveryDegree( <C>, <func> ) 
##
##  <C> is a complex, and <func> is a function operating on two conse-
##  cutive differentials, returning either true or false.
##
##  Returns true if func(d_i, d_{i+1}) is true for all i. Fails if this
##  is unknown, i.e. if the complex is infinite and not repeating.
##  
InstallMethod( ForEveryDegree, # TODO: misleading name?
[ IsQPAComplex, IsFunction ],
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

#######################################################################
##
#M  UpperBound( <C> ) 
##
##  Returns: 
##  If it exists: The smallest integer i such that the object at 
##  position i is non-zero, but for all j > i the object at position j 
##  is zero. If C is not a finite complex, the operation will return 
##  fail or infinity, depending on how C was defined.
##  
InstallMethod( UpperBound,
[ IsQPAComplex ],
function( C )
    local cat, diffs, positive, i;

    cat := CatOfComplex( C );
    diffs := DifferentialsOfComplex( C );
    positive := PositivePart( diffs );

    if IsZeroComplex( C ) then
      return -infinity;
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
      return infinity;
    else
        return fail;
    fi;
end );

#######################################################################
##
#M  LowerBound( <C> ) 
##
##  Returns:
##  If it exists: The greatest integer i such that the object at 
##  position i is non-zero, but for all j < i the object at position j 
##  is zero. If C is not a finite complex, the operation will return 
##  fail or negative infinity, depending on how C was defined.
##
InstallMethod( LowerBound,
[ IsQPAComplex ],
function( C )
    local cat, diffs, negative, i;

    cat := CatOfComplex( C );
    diffs := DifferentialsOfComplex( C );
    negative := NegativePart( diffs );

    if IsZeroComplex( C ) then
      return infinity;
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
      return -infinity;
    else
        return fail;
    fi;
end );

#######################################################################
##
#M  IsPositiveRepeating( <C> ) 
##
##  Returns true if the positive part of the differential list of the
##  complex <C> is repeating.
##
InstallMethod( IsPositiveRepeating,
[ IsQPAComplex ],
function( C )
    return IsRepeating( PositivePart( DifferentialsOfComplex( C ) ) );
end );

#######################################################################
##
#M  IsNegativeRepeating( <C> ) 
##
##  Returns true if the negative part of the differential list of the
##  complex <C> is repeating.
##
InstallMethod( IsNegativeRepeating,
[ IsQPAComplex ],
function( C )
    return IsRepeating( NegativePart( DifferentialsOfComplex( C ) ) );
end );

#######################################################################
##
#M  PositiveRepeatDegrees( <C> ) 
##
##  Returns a list describing the positions of the positive repeating
##  differentials.  The returned list is of the form [first .. last], 
##  where 'first' is the smallest degree such that the differential
##  ending there is in the positive repeating part. After the degree 
##  'last', the same differentials start repeating. Fails if
##  IsPositiveRepeating(C) is false.
##
InstallMethod( PositiveRepeatDegrees,
[ IsQPAComplex ],
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

#######################################################################
##
#M  NegativeRepeatDegrees( <C> ) 
##
##  Returns a list describing the positions of the negative repeating
##  differentials.  The returned list is of the form [last .. first], 
##  where 'last' is the greatest degree such that the differential
##  starting there is in the negative repeating part. After the degree 
##  'first', the same differentials start repeating. Fails if
##  IsNegativeRepeating(C) is false.
##
InstallMethod( NegativeRepeatDegrees,
[ IsQPAComplex ],
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

#######################################################################
##
#M  Shift( <C>, <i> ) 
##
##  Returns the complex C[i].  Note that the signs of the differentials
##  change if i is odd.
##
InstallMethod( Shift,
[ IsQPAComplex, IsInt ],
function( C, shift )
    local newDifferentials;

    newDifferentials := Shift( DifferentialsOfComplex( C ), shift );
    if shift mod 2 = 1 then
        newDifferentials := InfList( newDifferentials, d -> -d );
    fi;

    return ComplexByDifferentialList( CatOfComplex( C ), newDifferentials );

end );

#######################################################################
##
#M  ShiftUnsigned( <C>, <i> ) 
##
##  Returns a complex with the same objects as C[i], but with the 
##  differentials of C shifted without changing the signs.  A 'non-
##  algebraic' operation, but useful for manipulating complexes.
##
InstallMethod( ShiftUnsigned,
[ IsQPAComplex, IsInt ],
function( C, shift )
    local newDifferentials;
    
    newDifferentials := Shift( DifferentialsOfComplex( C ), shift );
    return ComplexByDifferentialList( CatOfComplex( C ), newDifferentials );

end );

#######################################################################
##
#M  YonedaProduct( <C1>, <C2> ) 
##
##  <C1> and <C2> are complexes such that the object in degree LowerBound(C)
##  equals the object in degree UpperBound(D).  The returned complex
##  is the Yoneda product of <C1> and <C2> (see the QPA documentation
##  for a more precise definition).
##  
InstallMethod( YonedaProduct,
[ IsQPAComplex, IsQPAComplex ],
function( C1, C2 )
    local cat, lowbound1, upbound2, diff1, diff2, connection, diffs;

    if IsZeroComplex( C1 ) then
      return C2;
    fi;

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

#######################################################################
##
#M  GoodTruncationBelow( <C>, <i> ) 
##
##  Not working at the moment.  Suppose that C is a complex
##    ... --> C_{i+1} --> C_i --> C_{i-1} --> ...
##
##  then the function should return the complex
##    ... --> C_{i+1} --> Z_i --> 0 --> 0 --> ...
##
##  where Z_i is the i-cycle of C.
##  
InstallMethod( GoodTruncationBelow,
[ IsQPAComplex, IsInt ],
function( C, i )
    local cat, difflist, truncpart, newpart, zeropart, newdifflist, kerinc;

    cat := CatOfComplex( C );
    difflist := DifferentialsOfComplex( C );
    truncpart := PositivePartFrom( difflist, i+2 );
    kerinc := KernelInclusion( DifferentialOfComplex( C, i ) );
    newpart := FiniteInfList( i, [ cat.zeroMap( Source(kerinc), cat.zeroObj ),
                                   LiftingInclusionMorphisms( kerinc, DifferentialOfComplex( C, i+1 ) ) ] );
    zeropart := NegativePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  i-1 );
    newdifflist := InfConcatenation( truncpart, newpart, zeropart );
    
    return ComplexByDifferentialList( cat, newdifflist );

end );

#######################################################################
##
#M  GoodTruncationAbove( <C>, <i> ) 
##
##  Not working at the moment.  Suppose that C is a complex
##    ... --> C_{i+1} --> C_i --> C_{i-1} --> ...
##
##  then the function should return the complex
##    ... --> 0 --> C_i/Z_i --> C_{i-1} --> 0 --> ...
##
##  where Z_i is the i-cycle of C.
##  
InstallMethod( GoodTruncationAbove,
 [ IsQPAComplex, IsInt ],
function( C, i )
    local cat, difflist, truncpart, newpart, zeropart, newdifflist, factor, factorinclusion,
          kerinc, factorproj;

    cat := CatOfComplex( C );
    difflist := DifferentialsOfComplex( C );
    truncpart := NegativePartFrom( difflist, i-1 );
    kerinc := KernelInclusion( DifferentialOfComplex( C, i ) );
    factorproj := CoKernelProjection( kerinc );
    factor := Range( factorproj );
    factorinclusion := cat.zeroMap( factor, ObjectOfComplex( C, i-1 ) ); #TODO
    newpart := FiniteInfList( i, [ factorinclusion, cat.zeroMap( cat.zeroObj, factor ) ] );

    zeropart := NegativePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  i+2 );
    newdifflist := InfConcatenation( zeropart, newpart, truncpart );
    
    return ComplexByDifferentialList( cat, newdifflist );

end );

# TODO!
# InstallMethod( GoodTruncation, [ IsQPAComplex, IsInt, IsInt ] );

#######################################################################
##
#M  BrutalTruncationBelow( <C>, <i> ) 
##
##  Suppose that C is a complex
##    ... --> C_{i+1} --> C_i --> C_{i-1} --> ...
##
##  then the function returns the complex
##    ... --> C_{i+1} --> C_i --> 0 --> 0 --> ...
##
InstallMethod( BrutalTruncationBelow,
[ IsQPAComplex, IsInt ],
function( C, i )
    local cat, difflist, truncpart, newpart, zeropart, newdifflist;
    
    cat := CatOfComplex( C );
    difflist := DifferentialsOfComplex( C );
    truncpart := PositivePartFrom( difflist, i+1 );
    newpart := FiniteInfList( i, [ cat.zeroMap( ObjectOfComplex( C, i), 
                                                cat.zeroObj) ] );
    zeropart := NegativePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  i-1 );
    newdifflist := InfConcatenation( truncpart, newpart, zeropart );
    
    return ComplexByDifferentialList( cat, newdifflist );

end );

#######################################################################
##
#M  BrutalTruncationAbove( <C>, <i> ) 
##
##  Suppose that C is a complex
##    ... --> C_{i+1} --> C_i --> C_{i-1} --> ...
##
##  then the function returns the complex
##    ... --> 0 --> C_i --> C_{i-1} --> 
##
InstallMethod( BrutalTruncationAbove,
[ IsQPAComplex, IsInt ],
function( C, i )
    local cat, difflist, truncpart, newpart, zeropart, newdifflist;
    
    cat := CatOfComplex( C );
    difflist := DifferentialsOfComplex( C );
    truncpart := NegativePartFrom( difflist, i );
    newpart := FiniteInfList( i+1, [ cat.zeroMap( cat.zeroObj,
                                                ObjectOfComplex( C, i )) ] );
    zeropart := PositivePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  i+2 );
    newdifflist := InfConcatenation( zeropart, newpart, truncpart );
    
    return ComplexByDifferentialList( cat, newdifflist );

end );

#######################################################################
##
#M  BrutalTruncation( <C>, <i>, <j> ) 
##
##  Suppose that C is a complex
##    ... --> C_{i+1} --> C_i --> C_{i-1} --> ...
##
##  then the function returns the complex
##    ... --> 0 --> C_i --> C_{i-1} --> ... --> C_j --> 0 --> ...
##
InstallMethod( BrutalTruncation, 
[ IsQPAComplex, IsInt, IsInt ],
function( C, i, j )
    local cat, difflist, middlediffs, truncpart, newpart1, zeropart1, 
          newpart2, zeropart2, newdifflist;
    
    if( j > i ) then
        Error( "First input integer must be greater than or equal to the second" );
    fi;
    
    return BrutalTruncationAbove( BrutalTruncationBelow( C, j ), i );
end );

#######################################################################
##
#O  SyzygyTruncation( <C>, <i> ) 
##
##  Suppose that C is a complex
##    ... --> C_{i+1} --> C_i --> C_{i-1} --> ...
##
##  then the function returns the complex
##    ... --> 0 --> ker(d_i) --> C_i --> C_{i-1} --> ...
##
InstallMethod( SyzygyTruncation, 
[ IsQPAComplex, IsInt ],
function( C, i )
    local cat, difflist, truncpart, kernelinc, newpart, kernel,
          zeropart, newdifflist;
    
    cat := CatOfComplex( C );
    difflist := DifferentialsOfComplex( C );
    truncpart := NegativePartFrom( difflist, i );

    kernelinc := KernelInclusion( DifferentialOfComplex( C, i ) );
    kernel := Source( kernelinc );
    newpart := FiniteInfList( i+1, [ kernelinc, 
                                     cat.zeroMap( cat.zeroObj, kernel ) ] );
    zeropart := PositivePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  i+3 );

    newdifflist := InfConcatenation( zeropart, newpart, truncpart );

    return ComplexByDifferentialList( cat, newdifflist );

end );

#######################################################################
##
#O  CosyzygyTruncation( <C>, <i> ) 
##
##  Suppose that C is a complex
##    ... --> C_{i+1} --> C_i --> C_{i-1} --> ...
##
##  then the function returns the complex
##    ... --> C_i --> C_{i-1} --> cok(d_i) --> 0 --> ...
##
InstallMethod( CosyzygyTruncation, 
[ IsQPAComplex, IsInt ],
function( C, i )
    local cat, difflist, truncpart, newpart,
          zeropart, newdifflist, cokerproj, coker;
    
    cat := CatOfComplex( C );
    difflist := DifferentialsOfComplex( C );
    truncpart := PositivePartFrom( difflist, i );

    cokerproj := CoKernelProjection( DifferentialOfComplex( C, i ) );
    coker := Range( cokerproj );
    newpart := FiniteInfList( i-2, [ cat.zeroMap( coker, cat.zeroObj ),
                                     cokerproj ] );

    zeropart := NegativePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  i-3 );

    newdifflist := InfConcatenation( truncpart, newpart, zeropart );

    return ComplexByDifferentialList( cat, newdifflist );

end );

#######################################################################
##
#O  SyzygyCosyzygyTruncation( <C>, <i>, <j> ) 
##
##  Suppose that C is a complex
##    ... --> C_{i+1} --> C_i --> C_{i-1} --> ...
##
##  then the function returns the complex
##    ... --> 0 --> ker(d_i) --> C_i --> ... --> C_{j+1} --> cok(d_j) --> 0 --> ...
##
InstallMethod( SyzygyCosyzygyTruncation, 
[ IsQPAComplex, IsInt, IsInt ],
function( C, i, j )
    local cat, difflist, truncpart, newdifflist, cokerproj, coker, 
          kernelinc, kernel, newpart1, zeropart1, newpart2, zeropart2, middlediffs;
    
    if( j > i ) then
        Error( "First input integer must be greater than or equal to the second" );
    fi;

    cat := CatOfComplex( C );

    difflist := DifferentialsOfComplex( C );
    middlediffs := FinitePartAsList( difflist, j, i );
    truncpart := FiniteInfList( j, middlediffs );

    kernelinc := KernelInclusion( DifferentialOfComplex( C, i ) );
    kernel := Source( kernelinc );
    newpart1 := FiniteInfList( i+1, [ kernelinc, 
                                     cat.zeroMap( cat.zeroObj, kernel ) ] );
    zeropart1 := PositivePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  i+3 );


    cokerproj := CoKernelProjection( DifferentialOfComplex( C, j ) );
    coker := Range( cokerproj );
    newpart2 := FiniteInfList( j-2, [ cat.zeroMap( coker, cat.zeroObj ),
                                     cokerproj ] );

    zeropart2 := NegativePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  j-3 );

    newdifflist := InfConcatenation( zeropart1, newpart1, truncpart, 
                                     newpart2, zeropart2 );

    return ComplexByDifferentialList( cat, newdifflist );
  
end );

#######################################################################
##
#O  CutComplexAbove( <C> ) 
##
##  For a bounded below complex C which is stored as an infinite complex,
##  but is known to be finite. Returns the same complex, but represented
##  as a finite complex.
##
InstallMethod( CutComplexAbove,
[ IsQPAComplex ],
function( C )
    local i, obj;
    if (IsInt(UpperBound(C))) then
        return C;
    else
       i := LowerBound(C);
       while true do
           obj := ObjectOfComplex(C, i);
           if (IsZero(DimensionVector(obj))) then
               return BrutalTruncationAbove(C, i-1);
           fi;
           i := i+1;
       od;
   fi;

end );               

#######################################################################
##
#O  CutComplexBelow( <C> ) 
##
##  For a bounded above complex C which is stored as an infinite complex,
##  but is known to be finite. Returns the same complex, but represented
##  as a finite complex.
##
InstallMethod( CutComplexBelow,
[ IsQPAComplex ],
function( C )
    local i, obj;
    if (IsInt(LowerBound(C))) then
        return C;
    else
       i := UpperBound(C);
       while true do
           obj := ObjectOfComplex(C, i);
           if (IsZero(DimensionVector(obj))) then
               return BrutalTruncationBelow(C, i+1);
           fi;
           i := i-1;
       od;
   fi;

end );               

#######################################################################
##
#O  ComplexByDifferentialList( <cat>, <differentials> ) 
##
##  <cat> is a category, and <differentials> is an InfList of
##  differentials.
##  
InstallMethod( ComplexByDifferentialList,
[ IsCat, IsInfList ],
function( cat, differentials )
    local C, fam;

    fam := NewFamily( "ComplexesFamily", IsQPAComplex );
    fam!.cat := cat;
    C := Objectify( NewType( fam, IsQPAComplex and IsQPAComplexDefaultRep ),
                    rec( ) );
    SetCatOfComplex( C, cat );
    SetDifferentialsOfComplex( C, differentials );
    return C;

end );

#######################################################################
##
#O  PrintObj( <C> ) 
##
##  Prints the zero complex
##  
InstallMethod( PrintObj,
[ IsZeroComplex ],
function( C )
    Print( "0 -> 0" );
end );

#######################################################################
##
#O  PrintObj( <C> ) 
##
##  Prints a non-zero complex
##  
InstallMethod( PrintObj,
[ IsQPAComplex ],
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

    if IsPositiveRepeating( C ) and upbound = infinity then
        Print( "[ " );
        for i in Reversed( PositiveRepeatDegrees( C ) ) do
            Print( i, ":", cat.objStr( ObjectOfComplex( C, i ) ), " -> " );
        od;
        Print( "] " );
    fi;

    if IsInt( top ) and IsInt( bottom ) then
      for i in [ top, top - 1 .. bottom ] do
        Print( i, ":", cat.objStr( ObjectOfComplex( C, i ) ), " -> " );
      od;
    fi;

    if IsNegativeRepeating( C ) and lowbound = -infinity then
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

#######################################################################
##
#O  ProjectiveResolution( <M> ) 
##
##  Returns the projective resolution of M, including M in degree -1.
##  
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
                    [ "next/repeat", nextDifferential, cover ],
                    "zero" );
end );

#######################################################################
##
#O  InjectiveResolution( < M > ) 
##
##  Returns the minimal injective resolution of M, including M in degree -1.
##  
InstallMethod( InjectiveResolution,
    [ IsPathAlgebraMatModule ],
    function( M )
    local nextDifferential, envelope;

    nextDifferential := function( d )
        return  CoKernelProjection( d )*InjectiveEnvelope( CoKernel( d ) );
    end;

    envelope := InjectiveEnvelope( M );

    return Complex( CatOfRightAlgebraModules( ActingAlgebra( M ) ),
                    0,
                    [ envelope ],
                    "zero",
                    [ "next/repeat", nextDifferential, envelope ]
                    );
end
);

#######################################################################
##
#F  ChainMap( <M>, <v> ) 
##
##  TODO: documentation for all chain map-functions
##  
InstallMethod( ChainMap,
               "for complexes, int and lists",
               [ IsQPAComplex, IsQPAComplex, IsInt, IsList, IsList, IsList ],
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

#######################################################################
##
#O  HighestKnownDegree( <map> ) 
##
##  Returns the highest degree that has been computed (or is otherwise
##  known) of the chain map <map>.  
##  
InstallMethod( HighestKnownDegree,
[ IsChainMap ],
function( map )
    return HighestKnownPosition( MorphismsOfChainMap( map ) );
end );

#######################################################################
##
#O  LowestKnownDegree( <map> ) 
##
##  Returns the lowest degree that has been computed (or is otherwise)
##  known) of the chain map <map>.
##  
InstallMethod( LowestKnownDegree,
[ IsChainMap ],
function( map )
    return LowestKnownPosition( MorphismsOfChainMap( map ) );
end );

#######################################################################
##
#O  MorphismOfChainMap( <map>, <i> ) 
##
##  Returns the morhpism in degree <i> of the map <map>.
##  
InstallMethod( MorphismOfChainMap,
[ IsChainMap, IsInt ],
function( map, i )
    return MorphismsOfChainMap( map )^i;
end );

#######################################################################
##
#O  \^( <map>, <i> ) 
##
##  Returns the morphism in degree <i> of the map <map>.
##  
InstallMethod( \^,
[ IsChainMap, IsInt ],
function( map, i )
    return MorphismsOfChainMap( map )^i;
end );

#######################################################################
##
#O  PrintObj( <map> ) 
##
##  Prints the chain map <map>.
##  
InstallMethod( PrintObj,
[ IsChainMap ],
function( map )
    Print( "<chain map>" );
end );

#######################################################################
##
#F  ComplexAndChainMAps( <sourceComplexes>, <rangeComplexes>,
##                       <basePosition>, <middle>, <positive>,
##                       <negative> ) 
##
##  Returns a list consisting of a newly created complex, togeher with
##  one or more newly created chain maps.  The new complex is either
##  source or range for the new chain map(s).
##
##  <sourceComplexes> is a list of the complexes to be sources of the
##  chain maps which will have the new complex as range.  Similarly,
##  <rangeComplexes> is a list of the complexes to be ranges of the new
##  chain maps which will have the new complex as source.
##  
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

#######################################################################
##
#F  FiniteChainMap( <source>, <range>, <baseDegree>, <morphisms> ) 
##
##  This function returns a finite chain map between two finite
##  complexes <source> and <range>.
##  
InstallMethod( FiniteChainMap,
               [ IsQPAComplex, IsQPAComplex, IsInt, IsList ],
function( source, range, baseDegree, morphisms )
    return ChainMap( source, range, baseDegree, morphisms, "zero", "zero" );
end );

#######################################################################
##
#F  ZeroChainMap( <source>, <range> )
##
##  This function returns the zero chain map between two
##  complexes <source> and <range>.
##  
InstallMethod( ZeroChainMap,
[ IsQPAComplex, IsQPAComplex ],
function( source, range )
    return ChainMap( source, range, 0, [], "zero", "zero" );
end );

#######################################################################
##
#O  ComparisonLifting( <f>, <PC>, <EC> )
##  
##  <f> is a morphism between two modules M and N.  <PC> is a complex
##  with M in degree i, zero in degrees j < i and projective modules
##  in degrees j > i.  <EC> is an exact complex with N in degree i
##  and zero in degrees j < i.  The function returns a chain map from
##  <PC> to <EC> which lifts the map <f> (it has <f> in degree i).
##  
InstallMethod( ComparisonLifting,
               [ IsPathAlgebraMatModuleHomomorphism, IsQPAComplex, IsQPAComplex ],
               function( f, PC, EC )
    local lbound, surjection, middle, nextLifting, nextLiftingFunction, i, j, chainmap;

    if ( LowerBound(PC) <> LowerBound(EC) ) then
       Error("The complexes must have the same lower bound");
    fi;

    lbound := LowerBound(PC);

    if (ObjectOfComplex(PC,lbound) <> Source( f ) or
        ObjectOfComplex(EC,lbound) <> Range( f ) ) then
       Error("The map has wrong source and/or range");
    fi;

    surjection := DifferentialOfComplex( PC , lbound + 1 )*f;

    middle := [ f, LiftingMorphismFromProjective( DifferentialOfComplex( EC, lbound + 1 ), surjection ) ];

    nextLifting := function( prev, i )
        local maps, kerinc, map, maps2, surj;

        if IsZero(DimensionVector( ObjectOfComplex( PC, i ))) then
            return ZeroMapping( ObjectOfComplex( PC, i ), ObjectOfComplex( EC, i ));
        fi;
        
        maps := ImageProjectionInclusion( DifferentialOfComplex( PC, i )*prev );
        kerinc := LiftingInclusionMorphisms( KernelInclusion( DifferentialOfComplex( EC, i-1 ) ), maps[2] );
        map := maps[1] * kerinc;
        maps2 := ImageProjectionInclusion( DifferentialOfComplex( EC, i ) );
        surj := maps2[1]*LiftingInclusionMorphisms( KernelInclusion( DifferentialOfComplex( EC, i-1 )), maps2[2]);
        
        return LiftingMorphismFromProjective(surj, map);
    end;

    if  IsInt( UpperBound( PC ) ) and IsInt( UpperBound ( EC )) then
        for i in [ lbound+2..UpperBound( PC ) ] do
            Append( middle, [ nextLifting( middle[i-lbound], i ) ] );
        od;
        for j in [ i+1..UpperBound( EC ) ] do
            Append( middle, [ ZeroMapping( ObjectOfComplex( PC, j ), ObjectOfComplex( EC, j ) ) ] );
        od;
        return FiniteChainMap( PC, EC, lbound, middle );
    else
        nextLiftingFunction := function( M, i ) return nextLifting( MorphismOfChainMap( M, i-1 ), i ); end; 
        return ChainMap( PC, EC, lbound, middle, [ "pos", nextLiftingFunction, true ], "zero" );
    fi;
end );


#######################################################################
##
#O  ComparisonLiftingToProjectiveResolution( <f> )
##  
##  <f> is a morphism between two modules M and N.  The function
##  returns the lifting of <f> to a chain map between the projective
##  resolutions of M and N.  The modules M and N are in degree -1 of
##  the resolutions, and the map in degree -1 is <f>.
##  
InstallMethod( ComparisonLiftingToProjectiveResolution,
[ IsPathAlgebraMatModuleHomomorphism ],
function( f )
    return ComparisonLifting( f, ProjectiveResolution( Source(f) ), ProjectiveResolution( Range(f) ) );
end );

#######################################################################
##
#O  MappingCone( <f> )
##  
##  <f> is a chain map between two complexes A and B.
##  This method returns a list [ C, i, p ] where C is the mapping cone
##  of <f>, i is the (chain map) inclusion of B into the cone and
##  p is the projection from the cone onto A[-1].
##  
InstallMethod( MappingCone,
[ IsChainMap ],
function( f )
    local C, i, A, B, dirsum, dirsum2, diff, middle, positiveFunction, negativeFunction;

    A := Source(f);
    B := Range(f);

#
#  Consider where to "start" the cone complex.
#
    if (IsInt(LowerBound(A))) then
        i := LowerBound(A);
    elif (IsInt(LowerBound(B))) then
        i := LowerBound(B);
    elif (IsInt(UpperBound(A))) then
        i := UpperBound(A);
    elif (IsInt(UpperBound(B))) then
        i := UpperBound(B);
    else
        i := 0;
    fi;

#
#  Construct the first differential of the cone, and the first projection/inclusion morphisms
#

    dirsum := DirectSumOfQPAModules( [ ObjectOfComplex(A,i-1), ObjectOfComplex(B,i) ] );
    dirsum2 := DirectSumOfQPAModules( [ ObjectOfComplex(A,i-2), ObjectOfComplex(B,i-1) ] );
    diff := MultiplyListsOfMaps( DirectSumProjections(dirsum),
                                 [[ -DifferentialOfComplex(A,i-1),
                                    ZeroMapping(ObjectOfComplex(B,i), ObjectOfComplex(A,i-2)) ],
                                  [ -MorphismOfChainMap(f,i-1), DifferentialOfComplex(B,i)]],
                                 DirectSumInclusions(dirsum2) );
    middle := [ [ diff, DirectSumInclusions(dirsum)[2], DirectSumProjections(dirsum)[1] ] ];

#
#  Positive degrees of the cone, the projection and the inclusion
#
    positiveFunction := function(C,inmap,outmap,i)
        local nextObj, prevObj, nextDiff;
        nextObj := DirectSumOfQPAModules( [ ObjectOfComplex(A,i-1), ObjectOfComplex(B,i) ] );
        prevObj := Source(DifferentialOfComplex(C,i-1));

        nextDiff :=  MultiplyListsOfMaps( DirectSumProjections(nextObj),
                                          [[ -DifferentialOfComplex(A,i-1),
                                             ZeroMapping(ObjectOfComplex(B,i), ObjectOfComplex(A,i-2)) ],
                                           [ -MorphismOfChainMap(f,i-1), DifferentialOfComplex(B,i)]],
                                          DirectSumInclusions(prevObj) );
        return [ nextDiff, DirectSumInclusions(nextObj)[2], DirectSumProjections(nextObj)[1] ];

    end ;

#
#  Negative degrees of the cone, the projection and the inclusion
#
    negativeFunction := function(C,inmap,outmap,i)
        local nextObj, prevObj, nextDiff;
        nextObj := DirectSumOfQPAModules( [ ObjectOfComplex(A,i-2), ObjectOfComplex(B,i-1) ] );
        prevObj := Range(DifferentialOfComplex(C,i+1));

        nextDiff := MultiplyListsOfMaps( DirectSumProjections(prevObj),
                                         [[ -DifferentialOfComplex(A,i-1),
                                            ZeroMapping(ObjectOfComplex(B,i), ObjectOfComplex(A,i-2)) ],
                                          [ -MorphismOfChainMap(f,i-1), DifferentialOfComplex(B,i)]],
                                         DirectSumInclusions(nextObj) );
        return [ nextDiff, DirectSumInclusions(prevObj)[2], DirectSumProjections(prevObj)[1] ];
    end ;

#
#  Returns the cone (as a complex) together with the projection and the inclusion
#
    C := ComplexAndChainMaps( [B], [Shift(A,-1)], i, middle,
                              ["pos", positiveFunction, true],
                              ["pos", negativeFunction, true] );
    return C;
end );

