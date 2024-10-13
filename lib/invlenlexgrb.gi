InstallMethod( TipInvLenInvLex,
  "for a path algebra element",
  IsElementOfPathRing,
  [ IsElementOfMagmaRingModuloRelations ], 0,
  function( e )
    local termlist, n, tip, fam, P, embedding;

    termlist := CoefficientsAndMagmaElements( e );
    n := Length( termlist ) - 1;
    if n < 1 then
       return Zero( e );
    else
       fam := FamilyObj( e );
       tip :=  ElementOfMagmaRing( fam,
                                   fam!.zeroRing,
                                  [ termlist[ 2 ] ], 
                                  [ termlist[ 1 ] ] );
       return tip;
    fi;
  end
);

InstallMethod( TipCoefficientInvLenInvLex,
  "for a path algebra element",
  IsElementOfPathRing,
  [ IsElementOfMagmaRingModuloRelations ], 0,
  function( e )
    local termlist, n, tip, fam, P, embedding;

    termlist := CoefficientsAndMagmaElements( e );
    n := Length( termlist ) - 1;
    fam := FamilyObj( e );
    if n < 1 then
       return Zero( LeftActingDomain( fam!.pathRing ) );
    else
       return termlist[ 2 ];
    fi;
  end
);

InstallMethod( TipMonomialInvLenInvLex,
  "for a path algebra element",
  IsElementOfPathRing,
  [ IsElementOfMagmaRingModuloRelations ], 0,
  function( e )
    local termlist, n, tip, fam, P, embedding;

    termlist := CoefficientsAndMagmaElements( e );
    n := Length( termlist ) - 1;
    if n < 1 then
        return Zero( e );
    else
	return termlist[ 1 ];
    fi;
  end
);

InstallMethod( TipWalkInvLenInvLex,
  "for a path algebra element",
  IsElementOfPathRing,
  [ IsElementOfMagmaRingModuloRelations ], 0,
  function( e )
    return WalkOfPath( TipMonomialInvLenInvLex( e ) );
  end
);

InstallMethod( TruncateElement,
  "truncate element modulo <arrows>^s",
  [ IsRingElement, IsPosInt ],
  function( elem, s )
  local temp, n, positions, paths, limit, truncation, truncated, i;

  temp := CoefficientsAndMagmaElements( elem );
  n := Length( temp );
  positions := Filtered( [ 1..n ], i -> i mod 2 <> 0 );
  paths := temp{ positions };
  limit := Length( Filtered( paths, p -> Length( WalkOfPath( p ) ) < s ) );
  truncated := Zero( elem );
  if limit = 0 then
    return truncated;
  fi;
  truncation := temp{ [ 1.. 2 * limit ] };
  for i in [ 1..limit ] do
    truncated := truncated + truncation[ 2 * i - 1 ] * One( elem ) * truncation[ 2 * i ];
  od;
  
  return truncated;
end
);


InstallMethod( GroebnerBasisInvLenInvLex,
  "compute the complete reduced Groebner Basis",
  [ IsList, IsPathAlgebra ],
  function( els, A )
  local B, t, Q, gb, el, el_tip, n, i, x, k, j, y, l, r, b, c, 
        overlap, remainder, min, max, num_paths, g, len, maxlenpaths;

  if not QPA_InArrowIdeal( els, A ) then
    Error( "elements do not belong to the arrow ideal of the path algebra,\n" );
  fi;
  B := A / els;
  t := LoewyLength( B ) + 1;
  Q := QuiverOfPathAlgebra( A );

  els := ReducedListInvLenInvLexQPA( MakeUniform( els ), A );
  gb := [];
  while Length( els ) > 0 do
    for el in els do
      el_tip := TipCoefficientInvLenInvLex( el );
      Add( gb, el * el_tip^(-1) );
    od;
    
    n := Length( gb );
    els := [];
    
    for i in [ 1..n ] do
      x := TipWalkInvLenInvLex( gb[ i ] );
      k := Length( x );

      for j in [ 1..n ] do
        y := TipWalkInvLenInvLex( gb[ j ] );
        l := Length( y );
        
        for r in [ Maximum( 0, k - l )..k - 1 ] do
          if x{ [ r + 1..k ] } = y{ [ 1..k - r ] } then
            b := x{ [ 1..r ] };
            c := y{ [ k - r + 1..l ] };
            
            overlap := gb[ i ] * Product( c, One( A ) ) - Product( b, One( A ) ) * gb[ j ];
            remainder := RemainderOfDivisionInvLenInvLexGroebnerBasis( overlap, gb, A );
            remainder := TruncateElement( remainder, t );
            
            if not IsZero( remainder ) then
              AddSet( els, remainder );
            fi;
          fi;
        od;
      od;
    od;
    gb := TipReducedListInvLenInvLex( gb, A );
    gb := ReducedListInvLenInvLexQPA( gb, A );
  od;
  
  return gb;
  end
);

InstallMethod( ReducedListInvLenInvLexQPA,
  "for a list of path-algebra elements",
  [ IsList, IsPathAlgebra ],
  function( els, A )
    local res, i, r;

    res := Filtered( els, el -> not IsZero( el ) );

    i := Length( res );
    while i > 0 do
      r := RemainderOfDivisionInvLenInvLexGroebnerBasis( res[ i ], 
                         res{ Concatenation( [ 1..i - 1 ], [ i + 1..Length( res ) ] ) }, A );

      if IsZero( r ) then
        Remove( res, i );
      else
        res[ i ] := r;
      fi;

      i := i - 1;
    od;

    return res;
  end
);

InstallMethod( TipReducedListInvLenInvLex,
  "for a list of path-algebra elements",
  [ IsList, IsPathAlgebra ],
  function( els, A )
    local res, el, i;

    res := [];

    for el in els do
      if not IsZero( el ) then
        AddSet( res, el );
      fi;
    od;

    i := Length( res );
    while i > 0 do
      if ForAny( [ 1..i - 1 ], j -> 
        LeftmostOccurrence( TipWalkInvLenInvLex( res[ i ] ), TipWalkInvLenInvLex( res[ j ] ) ) <> fail ) then
        Remove( res, i );
      fi;
      i := i - 1;
    od;

    return res;
  end
);

InstallMethod( RemainderOfDivisionInvLenInvLexGroebnerBasis,
  "for a path-algebra element and a list of path-algebra elements",
  [ IsElementOfMagmaRingModuloRelations, IsList, IsPathAlgebra ],
  function( y, X, A )
    local r, n, y_tip, y_wtip, divided, i, p, u, v;

    r := Zero( A );
    n := Length( X );

    while not IsZero( y ) do
      y_tip := TipInvLenInvLex( y );
      y_wtip := TipWalkInvLenInvLex( y_tip );

      divided := false;

      for i in [ 1..n ] do
        p := LeftmostOccurrence (y_wtip, TipWalkInvLenInvLex( X[ i ] ) );

        if p <> fail then
          u := Product( y_wtip{ [ 1..p[ 1 ] - 1 ] }, One( A ) );
          v := Product( y_wtip{ [ p[ 2 ] + 1..Length( y_wtip ) ] }, One( A ) );

          y := y - TipCoefficientInvLenInvLex( y_tip ) / TipCoefficientInvLenInvLex( X[ i ] ) * u * X[ i ] * v;

          divided := true;
          break;
        fi;
      od;

      if not divided then
        r := r + y_tip;
        y := y - y_tip;
      fi;
    od;

    return r;
  end
);

InstallMethod( LeftmostOccurrenceInvLenInvLex,
  "find second list as sublist of first list",
  [ IsList, IsList ],
  function(b, c)
    local lb, lc, i;

    lb := Length( b );
    lc := Length( c );

    for i in [ 1..lb - lc + 1 ] do
      if b{ [ i..i + lc - 1 ] } = c then
        return [ i, i + lc - 1 ];
      fi;
    od;

    return fail;
  end
);

InstallMethod( BiggestHomogeneousPartInvLenInvLex,
  "find second list as sublist of first list",
  [ IsElementOfMagmaRingModuloRelations ],
  function( e )
  local termlist, n, fam, d, i, bighom;

  if e = Zero( e ) then 
    return e;
  fi;
  termlist := CoefficientsAndMagmaElements( e );
  n := Length( termlist );
  if n < 1 then
    return Zero( e );
  else
    fam := FamilyObj( e );
    d := Length( WalkOfPath( termlist[ 1 ] ) );
    i := 1;
    bighom := Zero( e );
    while ( i <= n / 2 ) do
      if ( Length( WalkOfPath( termlist[ 2 * i - 1 ] ) )  = d ) then
        bighom := bighom + ElementOfMagmaRing( fam,
                                             fam!.zeroRing,
                                             [ termlist[ 2 * i ] ], 
                                             [ termlist[ 2 * i - 1 ] ] );
      fi;
      i := i + 1;
    od;
  fi;
  
  return bighom;
end
);

InstallMethod( BiggestHomogeneousPartInvLenInvLex,
  "find second list as sublist of first list",
  [ IsHomogeneousList ],
  function( elems )
  local fam;

    if elems = [] then 
      return elems;
    fi;
    if not ForAll( elems, IsElementOfMagmaRingModuloRelations ) then 
      Error( "Not all the entered elements are not elements of a path algebra,\n" );
    fi;
    fam := FamilyObj( elems[ 1 ] );    
    if not ForAll( elems, e -> FamilyObj( e ) = fam ) then
      Error ( "Not all the entered elements are from the same path algebra,\n" );
    fi;
    
    return List( elems, BiggestHomogeneousPartInvLenInvLex );
  end
);
