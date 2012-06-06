# GAP Implementation
# This file was generated from
# $Id: algpath.gi,v 1.6 2012/06/06 17:11:53 andrzejmroz Exp $


InstallMethod( IsPathRing,
    "for magma ring modulo the span of zero",
    true,
    [ IsMagmaRingModuloSpanOfZero ], 0,

   # Check to see underlying magma is a proper Quiver object:
    RM -> IsQuiver( UnderlyingMagma( RM ))
);


InstallMethod( IsPathRing,
    "for algebras that are not path rings",
    true,
    [ IsAlgebra ], 0,
    false
);


InstallGlobalFunction( PathRing,
  function(R, Q)
    local v, one, eFam;

    if not IsQuiver(Q) then
        Error( "<Q> must be a quiver." );
    fi;

    R := MagmaRingModuloSpanOfZero( R, Q, Zero(Q) );

    SetIsPathRing( R, true );
    eFam := ElementsFamily(FamilyObj(R));
    eFam!.pathRing := R;
    SetFilterObj(eFam,IsElementOfPathRing);

    # Create multiplicative identity element for this
    #  path ring; which is the sum of all vertices:
    one := Zero(R);
    for v in VerticesOfQuiver(Q) do
        one := one + ElementOfMagmaRing(eFam, eFam!.zeroRing, 
                                        [One(eFam!.zeroRing)], [v]);
    od;

    SetIsAssociative(R, true);
    SetOne(eFam, one);
    SetGeneratorsOfLeftOperatorRingWithOne(R, GeneratorsOfLeftOperatorRing(R));
    SetFilterObj(R, IsRingWithOne);

    return R;
  end
);


InstallGlobalFunction( PathAlgebra,
  function(F, Q)
    if not IsDivisionRing(F) then
        Error( "<F> must be a division ring or field." );
    fi;
    F := PathRing( F, Q );
    SetOrderingOfAlgebra(F, OrderingOfQuiver(Q));
    return F;
  end
);


InstallMethod( QuiverOfPathRing,
    "for a path ring",
    true,
    [ IsPathRing ], 0,
    RQ -> UnderlyingMagma(RQ)
);


InstallMethod( CanonicalBasis,
  "for path algebras",
  true,
  [ IsPathAlgebra ],
  NICE_FLAGS+1,
  function( A )
    local q, bv, i;

    q := QuiverOfPathAlgebra(A);
    if IsFinite(q) then
        bv := [];
        for i in Iterator(q) do
            Add(bv, i * One(A));
        od;
        return Basis(A, bv);
    else
       return fail;
    fi;
  end
);


# What's this for?
InstallMethod( NormalizedElementOfMagmaRingModuloRelations,
  "for elements of path algebras",
  true,
  [ IsElementOfMagmaRingModuloSpanOfZeroFamily 
    and IsElementOfPathRing, IsList ],
  0,
  function( Fam, descr )
        local zeromagma, len;

        zeromagma:= Fam!.zeroOfMagma;
        len := Length( descr[2] );
        if len > 0 and descr[2][1] = zeromagma then
            descr := [descr[1], descr[2]{[3..len]}];
            MakeImmutable(descr);
        fi;
        return descr;
  end
);


InstallMethod( \*,
  "for path algebra elements (faster)",
  IsIdenticalObj,
  [ IsZero and IsElementOfMagmaRingModuloRelations,
    IsElementOfMagmaRingModuloRelations ],
  0,
  function( x, y )
    return x;
  end
);


InstallMethod( \*,
  "for path algebra elements (faster)",
  IsIdenticalObj,
  [ IsElementOfMagmaRingModuloRelations,
    IsZero and IsElementOfMagmaRingModuloRelations ],
  0,
  function( x, y )
    return y;
  end
);


InstallMethod( \+,
  "for path algebra elements (faster)",
  IsIdenticalObj,
  [ IsZero and IsElementOfMagmaRingModuloRelations,
    IsElementOfMagmaRingModuloRelations ],
  0,
  function( x, y )
    return y;
  end
);


InstallMethod( \+,
  "for path algebra elements (faster)",
  IsIdenticalObj,
  [ IsElementOfMagmaRingModuloRelations,
    IsZero and IsElementOfMagmaRingModuloRelations ],
  0,
  function( x, y )
    return x;
  end
);


InstallGlobalFunction( FactorPathAlgebraByRelators,
  function( P, I, O )
    local A, fam, R, elementFam, relators, gb;
   
    relators := GeneratorsOfIdeal( I );

    if not IsIdenticalObj(OrderingOfQuiver(QuiverOfPathAlgebra(P)), O) then
        # Create a path algebra with a newly ordered quiver
        R := PathAlgebra(LeftActingDomain(P),
                     OrderedBy(QuiverOfPathAlgebra(P), O));
    
        # Create relators in $R$
        elementFam := ElementsFamily(FamilyObj(R));
        relators := List(relators, function(rel)
           local newRelator, extRep, zero, terms, i;

           newRelator := Zero(R);
           extRep := ExtRepOfObj(rel);
           zero := extRep[1];
           terms := extRep[2];

           for i in [1,3..Length(terms) - 1] do
                newRelator := newRelator + 
                              ObjByExtRep(elementFam, [zero, [terms[i], terms[i+1]]]);
           od;
           return newRelator;
        end);
    else
        R := P;
    fi;
    
    # Create a new family
    fam := NewFamily( "FamilyElementsFpPathAlgebra",
                      IsElementOfQuotientOfPathAlgebra );


    # Create the default type of elements
    fam!.defaultType := NewType( fam, IsElementOfQuotientOfPathAlgebra 
                                      and IsPackedElementDefaultRep );


    # If the ideal has a Groebner basis, create a function for
    # finding normal forms of the elements
    if HasGroebnerBasisOfIdeal( I ) 
#       and IsCompleteGroebnerBasis( GroebnerBasisOfIdeal( I ) )
    then
       fam!.normalizedType := NewType( fam, IsElementOfQuotientOfPathAlgebra
                                            and IsNormalForm
                                            and IsPackedElementDefaultRep );
       gb := GroebnerBasisOfIdeal( I );
       # Make sure the groebner basis is tip reduced now.
       TipReduceGroebnerBasis(gb);
       SetNormalFormFunction( fam, function(fam, x)
           local y;
           y := CompletelyReduce(gb, x);
           return Objectify( fam!.normalizedType, [y] );
       end );
    fi;

    fam!.pathAlgebra := R;
    fam!.ideal := I;
    fam!.familyRing := FamilyObj(LeftActingDomain(R));

    # Set the characteristic.
    if HasCharacteristic( R ) or HasCharacteristic( FamilyObj( R ) ) then
      SetCharacteristic( fam, Characteristic( R ) );
    fi;

    # Path algebras are always algebras with one
    A := Objectify(
        NewType( CollectionsFamily( fam ),
                IsQuotientOfPathAlgebra
            and IsAlgebraWithOne
            and IsWholeFamily
            and IsAttributeStoringRep ),
        rec() );

    SetLeftActingDomain( A, LeftActingDomain( R ) );
    SetGeneratorsOfAlgebraWithOne( A, 
        List( GeneratorsOfAlgebra( R ), 
            a -> ElementOfQuotientOfPathAlgebra( fam, a, false ) ) );

    SetZero( fam, ElementOfQuotientOfPathAlgebra( fam, Zero( R ), true ) );
    SetOne( fam, ElementOfQuotientOfPathAlgebra( fam, One( R ), true ) );
    UseFactorRelation( R, relators, A );

    SetOrderingOfAlgebra( A, O );
    SetQuiverOfPathAlgebra(A, QuiverOfPathAlgebra(R));
    SetIsFullFpPathAlgebra(A, true);

    fam!.wholeAlgebra := A;

    return A;
  end
);


InstallMethod( NaturalHomomorphismByIdeal,
  "for a path algebra and ideal",
  IsIdenticalObj,
  [ IsPathRing, IsFLMLOR ],
  0,
  function( A, I )
    local image, hom;

    image := FactorPathAlgebraByRelators( A, I, OrderingOfAlgebra(A) );

    if IsMagmaWithOne( A ) then
        hom := AlgebraWithOneHomomorphismByImagesNC( A, image,
                   GeneratorsOfAlgebraWithOne( A ),
                   GeneratorsOfAlgebraWithOne( image ) );
    else
        hom := AlgebraHomomorphismByImagesNC( A, image,
                   GeneratorsOfAlgebra( A ),
                   GeneratorsOfAlgebra( image ) );
    fi;

    SetIsSurjective( hom, true );

    return hom;
  end
);


# How is this used?
InstallMethod( \.,
  "for path algebras",
  true,
  [IsPathAlgebra, IsPosInt],
  0,
  function(A, name)
    local quiver, family;

    family := ElementsFamily(FamilyObj(A));
    quiver := QuiverOfPathAlgebra(A);
    return ElementOfMagmaRing(family, Zero(LeftActingDomain(A)),
        [One(LeftActingDomain(A))], [quiver.(NameRNam(name))] );
  end
);


InstallMethod( \<,
  "for path algebra elements",
  IsIdenticalObj,
  [IsElementOfMagmaRingModuloRelations, 
   IsElementOfMagmaRingModuloRelations], 10, # must be higher rank
  function( e, f )
    local i, swap;

    if not IsElementOfPathRing(FamilyObj(e)) then
        TryNextMethod();
    fi;
        
    e := Reversed(CoefficientsAndMagmaElements(e));
    f := Reversed(CoefficientsAndMagmaElements(f));

    # The reversal causes coefficients to come before
    # monomials. We need to swap these in order to
    # compare elements properly.
    for i in [ 1,3  .. Length( e )-1 ] do
        swap := e[i];
        e[i] := e[i+1];
        e[i+1] := swap;
    od;

    for i in [ 1,3  .. Length( f )-1 ] do
        swap := f[i];
        f[i] := f[i+1];
        f[i+1] := swap;
    od;

    for i in [ 1 .. Minimum(Length( e ), Length( f )) ] do
      if   e[i] < f[i] then
        return true;
      elif e[i] > f[i] then
        return false;
      fi;
    od;

    return Length( e ) < Length( f );

  end
);


InstallOtherMethod( LeadingTerm,
  "for a path algebra element",
  IsElementOfPathRing,
  [IsElementOfMagmaRingModuloRelations], 0,
  function(e)
    local termlist, n, tip, fam, P, embedding;

    termlist := CoefficientsAndMagmaElements(e);
    n := Length(termlist) - 1;
    if n < 1 then
       return Zero(e);
    else
       fam := FamilyObj(e);
       tip :=  ElementOfMagmaRing(fam,
                                  fam!.zeroRing,
                                  [termlist[n+1]], 
                                  [termlist[n]]);
       return tip;
    fi;

  end
);


InstallMethod( LeadingTerm,
  "for a quotient of path algebra element",
  true,
  [IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep],
  0,
  function(e)
    local lt, fam;

    fam := FamilyObj(e);
    lt :=  LeadingTerm(e![1]);

    return ElementOfQuotientOfPathAlgebra(fam, lt, IsNormalForm(e));

  end
);


InstallOtherMethod( LeadingCoefficient, 
  "for elements of path algebras",
  IsElementOfPathRing,
  [IsElementOfMagmaRingModuloRelations],
  0,
  function(e)
    local termlist, n, fam;

    termlist := CoefficientsAndMagmaElements(e);
    n := Length(termlist);
    if n < 1 then
        fam := FamilyObj(e);
        return fam!.zeroRing;
    else
        return termlist[n];
    fi;
  end
);


InstallOtherMethod( LeadingCoefficient, 
  "for elements of quotients of path algebras",
  true,
  [IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep], 0,
  function(e)

    return LeadingCoefficient(e![1]);

  end
);


InstallOtherMethod( LeadingMonomial,
  "for elements of path algebras",
  IsElementOfPathRing,
  [ IsElementOfMagmaRingModuloRelations ],
  0,
  function(e)
    local termlist, n, fam;

    termlist := CoefficientsAndMagmaElements(e);

    n := Length(termlist);

    if n < 1 then
        fam := FamilyObj(e);
        return fam!.zeroOfMagma;
    else
        return termlist[n-1];
    fi;

  end
);


InstallOtherMethod( LeadingMonomial, 
  "for elements of quotients of path algebras",
  true,
  [IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep],
  0,
  function(e)
      return LeadingMonomial(e![1]);
  end
);


InstallMethod( IsLeftUniform,
  "for path ring elements",
  IsElementOfPathRing,
  [ IsElementOfMagmaRingModuloRelations ],
  0,
  function( e )
      local terms, v;

      if IsZero(e) then
          return true;
      fi;

      terms := CoefficientsAndMagmaElements(e);
      v := SourceOfPath(terms[1]);
      return ForAll(terms{[1,3..Length(terms)-1]},
                    x -> SourceOfPath(x) = v);
  end
);


InstallMethod( IsLeftUniform,
  "for quotient of path algebra elements",
  true,
  [ IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep],
  0,
  function( e )
    return IsLeftUniform(e![1]);
  end
);


# Note: This will return true for zero entries in list
InstallMethod( IsLeftUniform,
  "for list of path ring elements",
  true,
  [ IsList, IsPath ],
  0,
  function( l, e )
      local terms, v, len, retflag, i;

      retflag := true;
      len := Length(l);

      # Strip off given vertex's coefficient:
#      terms := CoefficientsAndMagmaElements(e);
#      e := terms[1];
      terms := CoefficientsAndMagmaElements(e);
      e := terms[1];

      i := 1;
      while ( i <= len ) and ( retflag = true ) do
        terms := CoefficientsAndMagmaElements(l[i]);
        retflag :=  ForAll(terms{[1,3..Length(terms)-1]},
                      x -> SourceOfPath(x) = e);
        i := i + 1;
      od;

      return retflag;
  end
);



InstallMethod( IsRightUniform,
  "for path ring elements",
  IsElementOfPathRing,
  [ IsElementOfMagmaRingModuloRelations ],
  0,
  function( e )
    local terms, v;

    if IsZero(e) then
        return true;
    fi;

    terms := CoefficientsAndMagmaElements(e);
    v := TargetOfPath(terms[1]);
    return ForAll(terms{[1,3..Length(terms)-1]},
                  x -> TargetOfPath(x) = v);
  end
);


InstallMethod( IsRightUniform,
  "for quotient of path algebra elements",
  true,
  [ IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep], 0,
  function( e )
    return IsRightUniform(e![1]);
  end
);


# Note: This will return true for zero entries in list
InstallMethod( IsRightUniform,
  "for list of path ring elements",
  true,
  [ IsList, IsPath ],
  0,
  function( l, e )
      local terms, v, len, retflag, i;

      retflag := true;
      len := Length(l);

      # Strip off given vertex's coefficient:
#      terms := CoefficientsAndMagmaElements(e);
#      e := terms[1];
      terms := CoefficientsAndMagmaElements(e);
      e := terms[1];

      i := 1;
      while ( i <= len ) and ( retflag = true ) do
        terms := CoefficientsAndMagmaElements(l[i]);
        retflag :=  ForAll(terms{[1,3..Length(terms)-1]},
                      x -> TargetOfPath(x) = e);
        i := i + 1;
      od;

      return retflag;
  end
);



InstallMethod( IsUniform,
  "for path ring elements",
  IsElementOfPathRing,
  [ IsElementOfMagmaRingModuloRelations ], 0,
  function( e )

    if IsZero(e) then
        return true;
    fi;

    return IsRightUniform(e) and IsLeftUniform(e);
  end
);


InstallMethod( IsUniform,
  "for quotient of path algebra elements",
  true,
  [ IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep], 0,
  function( e )
    return IsUniform(e![1]);
  end
);


InstallOtherMethod( MappedExpression,
  "for a path algebra and two homogeneous lists",
  function( x, y, z)
    return IsElementOfPathRing(x)
           and IsElmsCollsX( x, y, z);
  end,
  [ IsElementOfMagmaRingModuloRelations, IsHomogeneousList,
    IsHomogeneousList ], 0,
  function( expr, gens1, gens2 )
    local MappedWord, genTable, i, mapped, one, coeff, mon;

    if IsZero(expr) then
        return Zero(gens2[1]);
    fi;

    expr := ExtRepOfObj( expr )[2];

    genTable := [];
    for i in [1..Length(gens1)] do
        coeff := LeadingCoefficient(gens1[i]);
        mon := LeadingMonomial(gens1[i]);
        if IsOne(coeff) and (IsVertex(mon) or IsArrow(mon)) then
            genTable[ExtRepOfObj(mon)[1]] := gens2[i];
        else
            Error( "gens1 must be a list of generators" );
        fi;
    od;

    one := One(expr[2]);

    MappedWord := function( word )
        local mapped, i, exponent;

        i := 2;
        exponent := 1;
        while i <= Length(word) and word[i] = word[i-1] do
            exponent := exponent + 1;
            i := i + 1;
        od;
        mapped := genTable[word[i-1]]^exponent;
        exponent := 1;

        while i <= Length(word) do
            i := i + 1;
            while i <= Length(word) and word[i] = word[i-1] do
                exponent := exponent + 1;
                i := i + 1;
            od;
            if exponent > 1 then
                mapped := mapped * genTable[word[i-1]]^exponent;
                exponent := 1;
            else
                mapped := mapped * genTable[word[i-1]];
            fi;
        od;

        return mapped;
    end;

    if expr[2] <> one then
        mapped := expr[2] * MappedWord(expr[1]);
    else
        mapped := MappedWord(expr[1]);
    fi;

    for i in [4,6 .. Length(expr)] do
        if expr[i] = one then
            mapped := mapped + MappedWord( expr[i-1] );
        else
            mapped := mapped + expr[i] * MappedWord( expr[ i-1 ] );
        fi;
    od;

    return mapped;
  end
);


InstallMethod( ImagesRepresentative,
  "for an alg. hom. from f. p. algebra, and an element",
  true,
  [ IsAlgebraHomomorphism
    and IsAlgebraGeneralMappingByImagesDefaultRep,
    IsRingElement ], 0,
  function( alghom, elem )

    local A;

    A := Source(alghom);
    if not (IsPathAlgebra(A) or IsQuotientOfPathAlgebra(A)) then
        TryNextMethod();
    fi;

    return MappedExpression( elem, alghom!.generators, alghom!.genimages );

  end
);


InstallOtherMethod( OrderedBy,
  "for a quotient of a path algebra",
  true,
  [IsQuotientOfPathAlgebra, IsQuiverOrdering], 0,
  function(A, O)
    local fam;
    fam := ElementsFamily(FamilyObj(A));
    return FactorPathAlgebraByRelators( fam!.pathAlgebra,
                                        fam!.relators,
                                        O);
  end
);


InstallMethod( \.,
  "for quotients of path algebras",
  true,
  [IsQuotientOfPathAlgebra, IsPosInt], 0,
  function(A, name)
    local parent, family;

    family := ElementsFamily(FamilyObj(A));
    parent := family!.pathAlgebra;

    return ElementOfQuotientOfPathAlgebra(family, parent.(NameRNam(name)), false );

  end
);


InstallMethod( RelatorsOfFpAlgebra,
  "for a quotient of a path algebra",
  true,
  [IsQuotientOfPathAlgebra and IsFullFpPathAlgebra], 0,
  A -> GeneratorsOfIdeal(ElementsFamily(FamilyObj(A))!.ideal)
);


InstallMethod( IsFiniteDimensional,
  "for quotients of path algebras",
  true,
  [IsQuotientOfPathAlgebra and IsFullFpPathAlgebra], 0,
  function( A )
    local gb, fam;

    fam := ElementsFamily(FamilyObj(A));
    gb := GroebnerBasisOfIdeal(fam!.ideal);
    if IsCompleteGroebnerBasis(gb) then
        return AdmitsFinitelyManyNontips(gb);
    elif IsFiniteDimensional(fam!.pathRing) then
        return true;
    else
        TryNextMethod();
    fi;
  end
);


InstallMethod( Dimension,
  "for quotients of path algebras",
  true,
  [IsQuotientOfPathAlgebra and IsFullFpPathAlgebra], 0,
  function( A )
    local gb, fam;

    fam := ElementsFamily(FamilyObj(A));
    gb := GroebnerBasisOfIdeal(fam!.ideal);
Print(gb,"\n");
    if IsCompleteGroebnerBasis(gb) then
        return NontipSize(gb);
    else
        TryNextMethod();
    fi;
  end
);


InstallMethod( CanonicalBasis,
  "for quotients of path algebras",
  true,
  [IsQuotientOfPathAlgebra], 0,
  function( A )
    local B, fam, zero, nontips, parent, parentFam, parentOne;

    fam := ElementsFamily( FamilyObj( A ) );
    zero := Zero(LeftActingDomain(A));
    parent := fam!.pathAlgebra;
    parentFam := ElementsFamily( FamilyObj( parent ) );
    parentOne := One(parentFam!.zeroRing);

    B := Objectify( NewType( FamilyObj( A ),
                        IsBasis and IsCanonicalBasisFreeMagmaRingRep ),
                    rec() );
    SetUnderlyingLeftModule( B, A );
    if HasGroebnerBasisOfIdeal( fam!.ideal )
       and IsCompleteGroebnerBasis( GroebnerBasisOfIdeal( fam!.ideal ) )
       and IsFiniteDimensional( A )
    then
        nontips := Nontips(GroebnerBasisOfIdeal( fam!.ideal ));
        nontips := List( nontips, 
                         x -> ElementOfMagmaRing( parentFam,
                                                  parentFam!.zeroRing,
                                                  [parentOne],
                                                  [x] ) );
        SetBasisVectors( B,
            List( EnumeratorSorted( nontips ), 
                  x -> ElementOfQuotientOfPathAlgebra( fam, x, true ) ) );
        B!.zerovector := List( BasisVectors( B ), x -> zero );
    fi;
    SetIsCanonicalBasis( B, true );
    return B;
  end
);


InstallMethod( BasisOfDomain,
  "for quotients of path algebras (CanonicalBasis)",
  true,
  [IsQuotientOfPathAlgebra], 10,
  CanonicalBasis
);


InstallMethod( Coefficients,
  "for canonical bases of quotients of path algebras",
  IsCollsElms,
  [IsCanonicalBasisFreeMagmaRingRep, 
   IsElementOfQuotientOfPathAlgebra and IsNormalForm], 0,
  function( B, e )
    local coeffs, data, elms, i, fam;

    data := CoefficientsAndMagmaElements( e![1] );
    coeffs := ShallowCopy( B!.zerovector );
    fam := ElementsFamily(FamilyObj( UnderlyingLeftModule( B ) ));
    elms := EnumeratorSorted( Nontips( GroebnerBasisOfIdeal(fam!.ideal)));
    for i in [1, 3 .. Length( data )-1 ] do
        coeffs[ PositionSet( elms, data[i] ) ] := data[i+1];
    od;
    return coeffs;
  end
);


InstallMethod( ElementOfQuotientOfPathAlgebra, 
  "for family of quotient path algebra elements and a ring element",
  true,
  [ IsElementOfQuotientOfPathAlgebraFamily, IsRingElement, IsBool ], 0,
  function( fam, elm, normal )
      return Objectify( fam!.defaultType, [ Immutable(elm) ]);
  end
);


InstallMethod( ElementOfQuotientOfPathAlgebra,
  "for family of quotient path algebra elements and a ring element (n.f.)",
  true,
  [ IsElementOfQuotientOfPathAlgebraFamily and HasNormalFormFunction,
    IsRingElement, IsBool ], 0,
  function( fam, elm, normal )
    if normal then
        return Objectify( fam!.normalizedType, [ Immutable(elm) ] );
    else
        return NormalFormFunction(fam)(fam, elm);
    fi;
  end
);


# Embeds a path (element of a quiver) into a path algebra (of the same quiver)
InstallMethod( ElementOfPathAlgebra, 
  "for a path algebra  and a path = element of a quiver",
  true,
  [ IsPathAlgebra, IsPath  ], 0,
  function( PA, path )
      local F;
	  
	  if not path in QuiverOfPathAlgebra(PA) then
	    Print("<path> should be an element of the quiver of the algebra <PA>\n");
		return false;
	  fi;
	  
	  F := LeftActingDomain(PA);
      return ElementOfMagmaRing(FamilyObj(Zero(PA)),Zero(F),[One(F)],[path]);
  end
);	



InstallMethod( ExtRepOfObj,
  "for element of a quotient of a path algebra",
  true,
  [ IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep ], 0,
  elm -> ExtRepOfObj( elm![1] )
);


InstallMethod( IsNormalForm,
  "for f.p. algebra elements",
  true,
  [IsElementOfQuotientOfPathAlgebra], 0,
  ReturnFalse
);


InstallMethod( \=,
   "for normal forms of elements of quotients of path algebras",
  IsIdenticalObj,
  [IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep and IsNormalForm,
   IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep and IsNormalForm],
  0,
  function(e, f)
    return ExtRepOfObj(e![1]) = ExtRepOfObj(f![1]);
  end
);


InstallMethod( \<,
  "for elements of quotients of path algebras",
  IsIdenticalObj,
  [IsElementOfQuotientOfPathAlgebra and IsNormalForm and IsPackedElementDefaultRep,
   IsElementOfQuotientOfPathAlgebra and IsNormalForm and IsPackedElementDefaultRep],
  0,
  function(e, f)
    return e![1] < f![1];
  end
);


InstallMethod(\+,
  "quotient path algebra elements",
  IsIdenticalObj,
  [IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra,
   IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra], 0,
  function(e, f)
    return ElementOfQuotientOfPathAlgebra(FamilyObj(e),e![1]+f![1], 
                                  IsNormalForm(e) and IsNormalForm(f));
  end
);


InstallMethod(\-,
  "quotient path algebra elements",
  IsIdenticalObj,
  [IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra,
   IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra], 0,
  function(e, f)
    return ElementOfQuotientOfPathAlgebra(FamilyObj(e),e![1]-f![1], 
                                  IsNormalForm(e) and IsNormalForm(f));
  end
);


InstallMethod(\*,
  "quotient path algebra elements",
  IsIdenticalObj,
  [IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra,
   IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra], 0,
  function(e, f)
    return ElementOfQuotientOfPathAlgebra(FamilyObj(e),e![1]*f![1], false);
  end
);


InstallMethod(\*,
  "ring el * quot path algebra el",
  IsRingsMagmaRings,
  [IsRingElement,
   IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra],0,
  function(e,f)
    return ElementOfQuotientOfPathAlgebra(FamilyObj(f),e*f![1], IsNormalForm(f));
  end
);


InstallMethod(\*,
  "quot path algebra el*ring el",
  IsMagmaRingsRings,
    [IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra,
     IsRingElement],0,
    function(e,f)
        return ElementOfQuotientOfPathAlgebra(FamilyObj(e),e![1]*f, IsNormalForm(e));
end);


InstallMethod(\*,
  "quiver element and quotient of path algebra element",
  true, # should check to see that p is in correct quiver
  [IsPath, 
   IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra], 0,
  function( p, e )
    return ElementOfQuotientOfPathAlgebra(FamilyObj(e), p*e![1], false);
  end
);


InstallMethod(\*,
  "quotient of path algebra element and quiver element",
  true, # should check to see that p is in correct quiver
  [IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra,
   IsPath], 0,
  function( e, p )
    return ElementOfQuotientOfPathAlgebra(FamilyObj(e), e![1] * p, false);
  end
);


InstallMethod(AdditiveInverseOp,
  "quotient path algebra elements",
  true,
  [IsPackedElementDefaultRep and IsElementOfQuotientOfPathAlgebra], 0,
  function(e)
    return
        ElementOfQuotientOfPathAlgebra(FamilyObj(e),AdditiveInverse(e![1]),
                               IsNormalForm(e));
  end
);


InstallOtherMethod( OneOp,
  "for quotient path algebra algebra element",
  true,
  [ IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep ], 0,
  function( e )
  local one;
  one:= One( e![1] );
  if one <> fail then
    one:= ElementOfQuotientOfPathAlgebra( FamilyObj( e ), one, true );
  fi;

  return one;

  end
);


InstallMethod( ZeroOp,
  "for a quotient path algebra element",
  true,
  [ IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep ], 0,
  e -> ElementOfQuotientOfPathAlgebra( FamilyObj( e ), Zero( e![1] ), true )
);


InstallOtherMethod( MappedExpression,
  "for f.p. path algebra, and two lists of generators",
  IsElmsCollsX,
  [ IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep,
    IsHomogeneousList, IsHomogeneousList ], 0,
  function( expr, gens1, gens2 )
    return MappedExpression( expr![1], List(gens1, x -> x![1]), gens2 );
  end
);


InstallMethod( PrintObj,
  "quotient of path algebra elements",
  true,
  [ IsElementOfQuotientOfPathAlgebra and IsPackedElementDefaultRep ], 0,
  function( e )
    Print( "[", e![1], "]" );
  end
);


InstallMethod( ObjByExtRep,
  "for family of f.p. algebra elements with normal form",
  true,
  [ IsElementOfQuotientOfPathAlgebraFamily,
    IsList ], 0,
  function( Fam, descr )
    local pathAlgFam;

    pathAlgFam := ElementsFamily(FamilyObj(Fam!.pathAlgebra));

    return ElementOfQuotientOfPathAlgebra(Fam, 
                  ObjByExtRep( pathAlgFam, descr ), false);
  end
);


InstallHandlingByNiceBasis( "IsFpPathAlgebraElementsSpace",
  rec(
    detect := function( F, gens, V, zero )
      return IsElementOfQuotientOfPathAlgebraCollection( V );
    end,

    NiceFreeLeftModuleInfo := function( V )
      local gens, monomials, gen, list, i, zero, info;

      gens := GeneratorsOfLeftModule( V );
      monomials := [];

      for gen in gens do
        list := CoefficientsAndMagmaElements( gen![1] );
        for i in [1, 3 .. Length(list) - 1] do
            AddSet( monomials, list[i] );
        od;
      od;

      V!.monomials := monomials;

      zero := Zero( V )![1]![1];

      info := rec(monomials := monomials,
                  zerocoeff := zero,
                  family := ElementsFamily( FamilyObj( V ) ) );

      if IsEmpty( monomials ) then
        info.zerovector := [ Zero( LeftActingDomain( V ) ) ];
      else
        info.zerovector := ListWithIdenticalEntries( Length( monomials ),
                                                         zero );
      fi;

      return info;
    end,

    NiceVector := function(V, v)
      local info, c, monomials, i, pos;

      info := NiceFreeLeftModuleInfo( V );
      c := ShallowCopy( info.zerovector );
      v := CoefficientsAndMagmaElements( v![1] );
      monomials := info.monomials;

      for i in [2, 4 .. Length(v)] do
        pos := PositionSet( monomials, v[ i-1 ] );
        if pos = fail then
            return fail;
        fi;
        c[ pos ] := v[i];
      od;

      return c;
    end,

    UglyVector := function( V, r )
      local elem, info, parentFam;
      info := NiceFreeLeftModuleInfo( V );
      if Length( r ) <> Length( info.zerovector ) then
        return fail;
      elif IsEmpty( info.monomials ) then
        if IsZero( r ) then
            return Zero( V );
        else
            return fail;
        fi;
      fi;
      parentFam := ElementsFamily( FamilyObj( info.family!.pathAlgebra ) );
      elem := ElementOfMagmaRing( parentFam, parentFam!.zeroRing,
                                  r, info.monomials );
      return ElementOfQuotientOfPathAlgebra( info.family, elem, true );
    end
  )
);
    
    
InstallGlobalFunction( PathAlgebraContainingElement,
        function( elem )
    if IsElementOfQuotientOfPathAlgebra( elem ) then
        return FamilyObj( elem )!.wholeAlgebra;
    else
        return FamilyObj( elem )!.pathRing;
    fi;
end );


# If there was an IsElementOfPathAlgebra filter, it would be
# better to write PathAlgebraContainingElement as an operation
# with the following two methods:
#
# InstallMethod( PathAlgebraContainingElement,
#         "for an element of a path algebra",
#         [ IsElementOfPathAlgebra ],
#         function( elem )
#     return FamilyObj( elem )!.pathRing;
# end );
#
# InstallMethod( PathAlgebraContainingElement,
#         "for an element of a quotient of a path algebra",
#         [ IsElementOfQuotientOfPathAlgebra ],
#         function( elem )
#     return FamilyObj( elem )!.wholeAlgebra;
# end );


InstallMethod( OriginalPathAlgebra,
   "for a quotient of a path algebra",
   [ IsAlgebra ],
   function( quot )
      if IsPathAlgebra( quot ) then
         return quot;
      elif IsQuotientOfPathAlgebra( quot ) then
         return ElementsFamily( FamilyObj( quot ) )!.pathAlgebra;
      else 
         Error("the algebra entered was not a quotient of a path algebra.");
      fi;
end );

InstallMethod( MakeUniformOnRight,
  "return array of right uniform path algebra elements",
  true,
  [ IsHomogeneousList ],
  0,
  function( gens )
    local l, t, v2, terms, uterms, newg, g, test;

    uterms := [];

    # get coefficients and monomials for all
    #  terms for gens:

    if Length(gens) = 0 then 
       return gens;
    else 
       test := IsPathAlgebra(PathAlgebraContainingElement(gens[1]));
       for g in gens do
          if test then 
             terms := CoefficientsAndMagmaElements(g);
          else
             terms := CoefficientsAndMagmaElements(g![1]);
          fi;    
          # Get terminus vertices for all terms:
          l := List(terms{[1,3..Length(terms)-1]}, x -> TargetOfPath(x));

          # Make the list unique:
          t := Unique(l);
 
          # Create uniformized array based on g:
          for v2 in t do
             newg := g*v2;
             if not IsZero(newg) then
                Add(uterms,newg);
             fi;
          od;
       od;
    fi;
    return Unique(uterms);
end
);

InstallMethod( GeneratorsTimesArrowsOnRight, 
   "for a path algebra",
   [ IsHomogeneousList ], 0,
   function( x ) 

   local A, fam, i, n, longer_list, extending_arrows, a, y;

   if Length(x) = 0 then
      Print("The entered list is empty.\n");
      return x;
   else
      A   := PathAlgebraContainingElement(x[1]);
      fam := ElementsFamily(FamilyObj(A));
      y   := MakeUniformOnRight(x);
      n   := Length(y); 

      longer_list := [];
      for i in [1..n] do
         extending_arrows := OutgoingArrowsOfVertex(TargetOfPath(TipMonomial(y[i]))); 
         for a in extending_arrows do
            if not IsZero(y[i]*a) then  
               Add(longer_list,y[i]*a);
            fi;
         od;      
      od; 

      return longer_list;
   fi;
end
);

InstallMethod( NthPowerOfArrowIdeal, 
   "for a path algebra",
   [ IsPathAlgebra, IS_INT ], 0,
   function( A, n ) 

   local num_vert, num_arrows, list, i;

   num_vert   := NumberOfVertices(QuiverOfPathAlgebra(A));
   num_arrows := NumberOfArrows(QuiverOfPathAlgebra(A));
   list := GeneratorsOfAlgebra(A){[1+num_vert..num_arrows+num_vert]};
   for i in [1..n-1] do 
      list := GeneratorsTimesArrowsOnRight(list);
   od;

   return list;
end
);

InstallMethod( TruncatedPathAlgebra, 
   "for a path algebra",
   [ IsField, IsQuiver, IS_INT ], 0,
   function( K, Q, n ) 

   local KQ, rels, I, gb, gbb; 

   KQ   := PathAlgebra(K,Q);
   rels := NthPowerOfArrowIdeal(KQ,n);
   I    := Ideal(KQ,rels);
   gb   := GBNPGroebnerBasis(rels,KQ);
   gbb  := GroebnerBasis(I,gb);

   return KQ/I;
end
);

InstallMethod( AddNthPowerToRelations, 
   "for a set of relations in a path algebra",
   [ IsPathAlgebra, IsHomogeneousList, IS_INT ], 0, 
   function ( pa, rels, n );
   
   Append(rels,NthPowerOfArrowIdeal(pa,n));   

   return rels; 

end
);