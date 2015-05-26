# GAP Implementation
# This file was generated from
# $Id: algpath.gi,v 1.8 2012/06/18 16:00:57 andrzejmroz Exp $

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
    if NumberOfArrows(Q) > 0 then 
        SetGlobalDimension(F, 1);
        SetFilterObj( F, IsHereditaryAlgebra ); 
    else
        SetGlobalDimension(F, 0);
        SetFilterObj( F, IsSemisimpleAlgebra ); 
    fi;
    if IsAcyclicQuiver(Q) then
        SetFilterObj( F, IsAdmissibleQuotientOfPathAlgebra);
    fi;
    return F;
  end
);

InstallMethod( QuiverOfPathRing,
    "for a path ring",
    true,
    [ IsPathRing ], 0,
    RQ -> UnderlyingMagma(RQ)
);

InstallMethod( ViewObj,
    "for a path algebra",
    true,
    [ IsPathAlgebra ], NICE_FLAGS + 1,
    function( A )

    Print("<");
    View(LeftActingDomain(A));
    Print("[");
    View(QuiverOfPathAlgebra(A));
    Print("]>");
end
);

InstallMethod( PrintObj,
    "for a path algebra",
    true,
    [ IsPathAlgebra ], NICE_FLAGS + 1,
    function( A )

    local Q;

    Print("<Path algebra of the quiver ");
    if HasName(QuiverOfPathAlgebra(A)) then
        Print(Name(QuiverOfPathAlgebra(A)));
    else
        View(QuiverOfPathAlgebra(A));
    fi;
    Print(" over the field ");
    View(LeftActingDomain(A));
    Print(">");
end
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

InstallMethod( ViewObj,
    "for a quotient of a path algebra",
    true,
    [ IsQuotientOfPathAlgebra ], NICE_FLAGS + 1,
    function( A )

    local fam;

    fam := ElementsFamily(FamilyObj(A));
    Print("<");
    View(LeftActingDomain(A));
    Print("[");
    View(QuiverOfPathAlgebra(A));
    Print("]/");
    View(fam!.ideal);
    Print(">");
end
);

InstallMethod( PrintObj,
    "for a quotient of a path algebra",
    true,
    [ IsQuotientOfPathAlgebra ], NICE_FLAGS + 1,
    function( A )

    local fam;

    fam := ElementsFamily(FamilyObj(A));
    Print("<A quotient of the path algebra ");
    View(OriginalPathAlgebra(A));
    Print(" modulo the ideal ");
    View(fam!.ideal);
    Print(">");
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
        if IsOne(coeff) and (IsQuiverVertex(mon) or IsArrow(mon)) then
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



################################################################
# Property IsFiniteDimensional for quotients of path algebras
# (analogue for path algebras uses standard GAP method)
# It uses Groebner bases machinery (computes G.b. only in case
# it has not been computed).

InstallMethod( IsFiniteDimensional,
  "for quotients of path algebras",
  true,
  [IsQuotientOfPathAlgebra and IsFullFpPathAlgebra], 0,
  function( A )
    local gb, fam, I, rels;
  
    fam := ElementsFamily(FamilyObj(A));
    I := fam!.ideal;

    if HasGroebnerBasisOfIdeal(I) then
      gb := GroebnerBasisOfIdeal(I);
    else
      rels := GeneratorsOfIdeal(I);     
      gb := GBNPGroebnerBasis(rels, fam!.pathAlgebra);
      gb := GroebnerBasis(I, gb);
    fi;
    
    if IsCompleteGroebnerBasis(gb) then
      return AdmitsFinitelyManyNontips(gb);
    elif IsFiniteDimensional(fam!.pathAlgebra) then
      return true;
    else
      TryNextMethod();
    fi;
  end
); # IsFiniteDimensional


################################################################
# Attribute Dimension for quotients of path algebras
# (analogue for path algebras uses standard GAP method,
# note that for infinite dimensional path algebras can fail!)
# It uses Groebner bases machinery (computes G.b. only in case
# it has not been computed).
# It returns "infinity" for infinite dimensional algebra.

InstallMethod( Dimension,
  "for quotients of path algebras",
  true,
  [IsQuotientOfPathAlgebra and IsFullFpPathAlgebra], 0,
  function( A )
    local gb, fam, I, rels;

    fam := ElementsFamily(FamilyObj(A));
    I := fam!.ideal;
    
    if HasGroebnerBasisOfIdeal(I) then
      gb := GroebnerBasisOfIdeal(I);
    else
      rels := GeneratorsOfIdeal(I);     
      gb := GBNPGroebnerBasis(rels, fam!.pathAlgebra);
      gb := GroebnerBasis(I, gb);
    fi;

    if IsCompleteGroebnerBasis(gb) then
      return NontipSize(gb);
    else
      TryNextMethod();
    fi;
  end
); # Dimension

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
   [ IsQuiverAlgebra ],
   function( quot )
      if IsPathAlgebra( quot ) then
         return quot;
      elif IsQuotientOfPathAlgebra( quot ) then
         return ElementsFamily( FamilyObj( quot ) )!.pathAlgebra;
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

    local KQ, rels; 

    KQ   := PathAlgebra(K,Q);
    rels := NthPowerOfArrowIdeal(KQ,n);

    return KQ/rels;
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

InstallMethod( OppositePathAlgebra,
        "for a path algebra",
        [ IsPathAlgebra ],
        function( PA )
    local PA_op, field, quiver_op;

    field := LeftActingDomain( PA );
    quiver_op := OppositeQuiver( QuiverOfPathAlgebra( PA ) );
    PA_op := PathAlgebra( field, quiver_op );

    SetOppositePathAlgebra( PA_op, PA );

    return PA_op;
end );


InstallMethod( OppositeAlgebra,
        "for a path algebra",
        [ IsPathAlgebra ],
        OppositePathAlgebra );


InstallGlobalFunction( OppositePathAlgebraElement,
        function( elem )
    local PA, PA_op, terms, op_term;

    PA := PathAlgebraContainingElement( elem );
    PA_op := OppositePathAlgebra( PA );

    if elem = Zero(PA) then 
       return Zero(PA_op);
    else 
       op_term :=
       function( t )
          return t.coeff * One( PA_op ) * OppositePath( t.monom );
       end;

       terms := PathAlgebraElementTerms( elem );
       return Sum( terms, op_term );
    fi;
end );


InstallMethod( OppositeRelations,
        "for a list of relations",
        [ IsDenseList ],
        function( rels )
    return List( rels, OppositePathAlgebraElement );
end );
        
        
InstallMethod( OppositePathAlgebra,
        "for a quotient of a path algebra",
        [ IsQuotientOfPathAlgebra ],
        function( quot )
    local PA, PA_op, rels, rels_op, I_op, gb, gbb, quot_op;

    PA := OriginalPathAlgebra( quot );
    PA_op := OppositePathAlgebra( PA );
    rels  := RelatorsOfFpAlgebra( quot );
    rels_op := OppositeRelations( rels );
    I_op := Ideal( PA_op, rels_op );
    gb   := GBNPGroebnerBasis(rels_op,PA_op);
    gbb  := GroebnerBasis(I_op,gb);
    quot_op := PA_op / I_op;
    
#    if HasIsAdmissibleQuotientOfPathAlgebra(quot) and IsAdmissibleQuotientOfPathAlgebra(quot) then
#        SetIsAdmissibleQuotientOfPathAlgebra(quot_op, true);
#    fi; 
    if IsAdmissibleQuotientOfPathAlgebra(quot) then
        SetFilterObj(quot_op, IsAdmissibleQuotientOfPathAlgebra );
    fi; 

    SetOppositePathAlgebra( quot_op, quot );

    return quot_op;
end );


InstallMethod( OppositeAlgebra,
        "for a quotient of a path algebra",
        [ IsQuotientOfPathAlgebra ],
        OppositePathAlgebra );

InstallMethod( Centre,
   "for a path algebra",
   [ IsQuotientOfPathAlgebra ], 0,
   function( A ) 

   local Q, K, num_vert, vertices, arrows, B, cycle_list, i, j, b, 
         pos, commutators, center, zero_count, a, c, x, cycles, matrix,
         data, coeffs, fam, elms, solutions, gbb;
   
   if not IsFiniteDimensional(A) then
       Error("the entered algebra is not finite dimensional,\n");
   fi;
   B := CanonicalBasis(A);
   Q := QuiverOfPathAlgebra(A); 
   K := LeftActingDomain(A);
   num_vert := NumberOfVertices(Q);
   vertices := VerticesOfQuiver(Q);
   arrows   := GeneratorsOfAlgebra(A){[2+num_vert..1+num_vert+NumberOfArrows(Q)]};
  
   cycle_list := [];
   for b in B do
      if SourceOfPath(TipMonomial(b)) = TargetOfPath(TipMonomial(b)) then 
         Add(cycle_list,b);
      fi;
   od;
   commutators := [];
   center := [];
   cycles := Length(cycle_list);
   for c in cycle_list do
      for a in arrows do
         Add(commutators,c*a-a*c);
      od;
   od;

   matrix := NullMat(cycles,Length(B)*NumberOfArrows(Q),K);

   for j in [1..cycles] do
      for i in [1..NumberOfArrows(Q)] do
         matrix[j]{[Length(B)*(i-1)+1..Length(B)*i]} := Coefficients(B,commutators[i+(j-1)*NumberOfArrows(Q)]); 
      od;
   od;

   solutions := NullspaceMat(matrix);

   center := [];
   for i in [1..Length(solutions)] do
      x := Zero(A);
      for j in [1..cycles] do 
         if solutions[i][j] <> Zero(K) then 
            x := x + solutions[i][j]*cycle_list[j];
         fi;
      od;
      center[i] := x;
   od;
   return SubalgebraWithOne(A,center,"basis");
end
);



#######################################################################
##
#O  VertexPosition(<elm>)
##
##  This function assumes that the input is a residue class of a trivial
##  path in finite dimensional quotient of a path algebra, and it finds  
##  the position of this trivial path/vertex in the list of vertices for
##  the quiver used to define the algebra.
##
InstallMethod ( VertexPosition, 
   "for an element in a quotient of a path algebra",
   true,
   [ IsElementOfQuotientOfPathAlgebra ],
   0,
   function( elm )

   return elm![1]![2][1]!.gen_pos;
end
);

#######################################################################
##
#O  VertexPosition(<elm>)
##
##  This function assumes that the input is a trivial path in finite 
##  dimensional a path algebra, and it finds the position of this 
##  trivial path/vertex in the list of vertices for the quiver used 
##  to define the algebra.
##
InstallOtherMethod ( VertexPosition, 
   "for an element in a PathAlgebra",
   true,
   [ IsElementOfMagmaRingModuloRelations ],
   0,
   function( elm )

   if "pathRing" in NamesOfComponents(FamilyObj(elm)) then 
      return elm![2][1]!.gen_pos;
   else
      return false;
   fi;
end
);

#######################################################################
##
#O  IsSelfinjectiveAlgebra ( <A> )
##
##  This function returns true or false depending on whether or not 
##  the algebra  <A>  is selfinjective.
##
##  This test is based on the paper by T. Nakayama: On Frobeniusean 
##  algebras I. Annals of Mathematics, Vol. 40, No.3, July, 1939. It 
##  assumes that the input is a basic algebra, which is always the 
##  case for a quotient of a path algebra modulo an admissible ideal.
##
##
InstallMethod( IsSelfinjectiveAlgebra, 
    "for a finite dimension (quotient of) a path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A ) 

    local Q, soclesofproj, test; 
   
    if not IsFiniteDimensional(A) then 
        return false;
    fi;
    if not IsPathAlgebra(A) and not IsAdmissibleQuotientOfPathAlgebra(A) then 
        TryNextMethod();
    fi;
    #
    # If the algebra is a path algebra.
    #
    if IsPathAlgebra(A) then 
        Q := QuiverOfPathAlgebra(A);
        if NumberOfArrows(Q) > 0 then 
            return false;
        else
            return true;
        fi;
    fi;
    #
    # By now we know that the algebra is a quotient of a path algebra.
    # First checking if all indecomposable projective modules has a simpel
    # socle, and checking if all simple modules occur in the socle of the
    # indecomposable projective modules. 
    #
    soclesofproj := List(IndecProjectiveModules(A), p -> SocleOfModule(p));
    if not ForAll(soclesofproj, x -> Dimension(x) = 1) then
        return false;
    fi;
    test := Sum(List(soclesofproj, x -> DimensionVector(x)));
    
    return ForAll(test, t -> t = 1); 
end
);



InstallOtherMethod( CartanMatrix, 
    "for a finite dimensional (quotient of) path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A ) 

    local P, C, i;
   
    if not IsFiniteDimensional(A) then 
	Print("Algebra is not finite dimensional!\n");
        return fail;
    fi;
    if not IsPathAlgebra(A) and not IsAdmissibleQuotientOfPathAlgebra(A) then 
	Print("Not a bound quiver algebra!\n");
        return fail;
    fi;
    P := IndecProjectiveModules(A);
    C := [];
    for i in [1..Length(P)] do
        Add(C,DimensionVector(P[i]));
    od;

    return C;
end
); # CartanMatrix

InstallMethod( CoxeterMatrix, 
    "for a finite dimensional (quotient of) path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A ) 

    local P, C, i;

    C := CartanMatrix(A);
    if C = fail then
	Print("Unable to determine the Cartan matrix!\n");
        return fail;
    fi;
    if DeterminantMat(C) <> 0 then 
        return (-1)*C^(-1)*TransposedMat(C);
    else
        Print("The Cartan matrix is not invertible!\n");
	return fail;
    fi;
end
); # CoxeterMatrix

InstallMethod( CoxeterPolynomial, 
    "for a finite dimensional (quotient of) path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A ) 

    local P, C, i;

    C := CoxeterMatrix(A);
    if C <> fail then 
        return CharacteristicPolynomial(C);
    else
        return fail;
    fi;
end
); # CoxeterPolynomial


InstallMethod( TipMonomialandCoefficientOfVector, 
   "for a path algebra",
   [ IsQuiverAlgebra, IsCollection ], 0,
   function( A, x ) 

   local pos, tipmonomials, sortedtipmonomials, tippath, i, n;

   pos := 1;
   tipmonomials := List(x,TipMonomial);
   sortedtipmonomials := ShallowCopy(tipmonomials);
   Sort(sortedtipmonomials,\<);

   if Length(tipmonomials) > 0 then
      tippath := sortedtipmonomials[Length(sortedtipmonomials)];
      pos := Minimum(Positions(tipmonomials,tippath));
   fi;

   return [pos,TipMonomial(x[pos]),TipCoefficient(x[pos])];
end
);


InstallMethod( TipReduceVectors, 
   "for a path algebra",
   [ IsQuiverAlgebra, IsCollection ], 0,
   function( A, x ) 

   local i, j, k, n, y, s, t, pos_m, z, stop;

   n := Length(x);
   if n > 0 then 
      for k in [1..n] do
      	 for j in [1..n] do
            if j <> k then  
               s := TipMonomialandCoefficientOfVector(A,x[k]);
               t := TipMonomialandCoefficientOfVector(A,x[j]); 
               if ( s <> fail ) and ( t <> fail ) then
                  if ( s[1] = t[1] ) and ( s[2] = t[2] ) and 
                     ( s[3] <> Zero(LeftActingDomain(A)) ) then 
               	     x[j] := x[j] - (t[3]/s[3])*x[k];
                  fi;
               fi;
            fi;
      	 od;
      od;
      return x;
   fi;
end
);

#
# This has a bug !!!!!!!!!!!!!!!!!!!
#
InstallMethod( CoefficientsOfVectors, 
   "for a path algebra",
   [ IsAlgebra, IsCollection, IsList ], 0,
   function( A, x, y ) 

   local i, j, n, s, t, tiplist, vector, K, zero_vector, z;

   tiplist := [];
   for i in [1..Length(x)] do 
      Add(tiplist,TipMonomialandCoefficientOfVector(A,x[i]){[1..2]});
   od;

   vector := []; 
   K := LeftActingDomain(A);
   for i in [1..Length(x)] do
      Add(vector,Zero(K));
   od;
   zero_vector := [];
   for i in [1..Length(y)] do
      Add(zero_vector,Zero(A));
   od;

   z := y;
   if z = zero_vector then
      return vector;
   else 
      repeat
         s := TipMonomialandCoefficientOfVector(A,z);
         j := Position(tiplist,s{[1..2]});
         if j = fail then 
            return [x,y];
         fi;
         t := TipMonomialandCoefficientOfVector(A,x[j]); 
         z := z - (s[3]/t[3])*x[j];
         vector[j] := vector[j] + (s[3]/t[3]);
      until 
         z = zero_vector;
      return vector;
   fi;
end
);




#######################################################################
##
#P  IsSymmetricAlgebra( <A> )
##
##  This function determines if the algebra  A  is a symmetric algebra,
##  if it is a (quotient of a) path algebra. 
##
InstallMethod( IsSymmetricAlgebra, 
    "for a quotient of a path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A ) 

    local Q, fam, Aenv, M, DM, mats, op_name, de_op_name, vertices, arrows, 
          new_vertices, new_arrows, stringvertices, stringarrows, 
          vertex_positions, arrow_positions, newdimvector, newmats, 
          matrices, a, arrowentry, MM;
    
    if not IsFiniteDimensional(A) then 
        return false;
    fi;
    if IsPathAlgebra(A) then
        Q := QuiverOfPathAlgebra(A);
        return Length(ArrowsOfQuiver(QuiverOfPathAlgebra(A))) = 0;
    fi;    
    #
    # By now we know that the algebra is a finite dimensional quotient of a path algebra.
    #
    Aenv := EnvelopingAlgebra(A);
    M    := AlgebraAsModuleOverEnvelopingAlgebra(A);
    DM   := DualOfModule(M);
    #
    #   Finding DM as a module over Aenv.
    #
    mats    := MatricesOfPathAlgebraModule(DM);
    op_name := OppositeQuiverNameMap(QuiverOfPathAlgebra(A));        
    de_op_name := OppositeQuiverNameMap(OppositeQuiver(QuiverOfPathAlgebra(A)));
    vertices   := VerticesOfQuiver(QuiverOfPathAlgebra(Aenv));
    arrows     := ArrowsOfQuiver(QuiverOfPathAlgebra(Aenv));
    new_vertices := List(vertices, x -> Concatenation(op_name(String(ProjectFromProductQuiver(2,x))),"_",de_op_name(String(ProjectFromProductQuiver(1,x)))));
    new_arrows   := List(arrows, x -> Concatenation(op_name(String(ProjectFromProductQuiver(2,x))),"_",de_op_name(String(ProjectFromProductQuiver(1,x)))));
    stringvertices := List(vertices, x -> String(x));
    stringarrows   := List(arrows, x -> String(x));
    #
    #   Finding the permutation of the vertices and the arrows.
    #
    vertex_positions := List(new_vertices, x -> Position(stringvertices, x));
    arrow_positions  := List(new_arrows, x -> Position(stringarrows, x));
    #
    #   Finding the new dimension vector and the matrices of  DM  as a module over Aenv.
    #
    newdimvector := List([1..Length(vertices)], x -> DimensionVector(DM)[vertex_positions[x]]);
    newmats  := List([1..Length(mats)], x -> mats[arrow_positions[x]]);
    #
    #   Creating the input for construction  DM  as a module over Aenv.
    #
    matrices := [];
    for a in arrows do
        if newdimvector[Position(vertices,SourceOfPath(a))] <> 0 and 
           newdimvector[Position(vertices,TargetOfPath(a))] <> 0 then 
            arrowentry := List([1..2], x -> []);
            arrowentry[1] := String(a);
            arrowentry[2] := newmats[Position(arrows,a)];
            Add(matrices, arrowentry);
        fi;
    od;
    MM := RightModuleOverPathAlgebra(Aenv, newdimvector, matrices);
    
    return IsomorphicModules(M,MM);
end
);

#######################################################################
##
#P  IsWeaklySymmetricAlgebra( <A> )
##
##  This function determines if the algebra  A  is a weakly symmetric 
##  algebra, if it is a (quotient of a) path algebra. 
##
InstallMethod( IsWeaklySymmetricAlgebra, 
   "for a quotient of a path algebra",
   [ IsQuiverAlgebra ], 0,
   function( A ) 

   local P;
   
   if IsSelfinjectiveAlgebra(A) then
       P := IndecProjectiveModules(A);
       return ForAll(P, x -> DimensionVector(SocleOfModule(x)) = DimensionVector(TopOfModule(x)));
   else   
       return false; 
   fi;
end
);

InstallOtherMethod( \/,
    "for PathAlgebra and a list of generators",
    IsIdenticalObj,
    [ IsPathAlgebra, IsList ],
    function( KQ, relators )
    
    local gb, I, A;
    
    if not ForAll(relators, x -> ( x in KQ ) ) then
        Error("the entered list of elements should be in the path algebra!");
    fi;
    gb := GBNPGroebnerBasis(relators, KQ);
    I := Ideal(KQ,gb);
    GroebnerBasis(I,gb);
    A := KQ/I; 
#    if IsAdmissibleIdeal(I) then 
#        SetIsAdmissibleQuotientOfPathAlgebra(A, true);
#    fi;
    if IsAdmissibleIdeal(I) then 
        SetFilterObj(A, IsAdmissibleQuotientOfPathAlgebra );
    fi;

    return A;
end 
); 

########################################################################
##
#P  IsSchurianAlgebra( <A> ) 
##  <A> = a path algebra or a quotient of a path algebra
##
##  It tests if an algebra <A> is a schurian  algebra.
##  By definition it means that:
##  for all x,y\in Q_0 dim A(x,y)<=1.
##  Note: This method fail when a Groebner basis
##  for ideal has not been computed before creating q quotient!
##
InstallMethod( IsSchurianAlgebra,
    "for PathAlgebra or QuotientOfPathAlgebra",
    true,
    [ IsQuiverAlgebra ],
    function( A )
    
    local Q, test;
    
    if not IsPathAlgebra(A) and not IsAdmissibleQuotientOfPathAlgebra(A) then 
       TryNextMethod();
    fi;
    if IsFiniteDimensional(A) then
        test := Flat(List(IndecProjectiveModules(A), x -> DimensionVector(x)));
        return ForAll(test, x -> ( x <= 1 ) );
    else
        Error("the entered algebra is not finite dimensional,\n");
    fi;    
end 
); 

#################################################################
##
#P  IsSemicommutativeAlgebra( <A> ) 
##  <A> = a path algebra
##
##  It tests if a path algebra <A> is a semicommutative  algebra. 
##  Note that for a path algebra it is an empty condition
##  (when <A> is schurian + acyclic, cf. description for quotients below).
##		
InstallMethod( IsSemicommutativeAlgebra,
    "for path algebras",
    true,
    [ IsPathAlgebra ], 0,
    function ( A )
    
    local Q; 
    
    if not IsSchurianAlgebra(A) then
        return false;
    fi;
    
    Q := QuiverOfPathAlgebra(A);
    return IsAcyclicQuiver(Q);  
end
); # IsSemicommutativeAlgebra for IsPathAlgebra


########################################################################
##
#P  IsSemicommutativeAlgebra( <A> ) 
##  <A> = a quotient of a path algebra
##
##  It tests if a quotient of a path algebra <A>=kQ/I is a semicommutative  algebra.
##  By definition it means that:
##  1. A is schurian (i.e. for all x,y\in Q_0 dim A(x,y)<=1).
##  2. Quiver Q of A is acyclic.
##  3. For all pairs of vertices (x,y) the following condition is satisfied:
##     for every two paths P,P' from x to y:
##     P\in I <=> P'\in I
##
  
InstallMethod( IsSemicommutativeAlgebra,
    "for quotients of path algebras",
    true,
    [ IsQuotientOfPathAlgebra ], 0,
    function ( A )
    
    local Q, PA, I, vertices, vertex, v, vt,
          paths, path, p,
          noofverts, noofpaths,
          eA, pp,
          inI, notinI;
    
    PA := OriginalPathAlgebra(A);
    Q := QuiverOfPathAlgebra(PA);
    if (not IsAcyclicQuiver(Q)) 
      or (not IsSchurianAlgebra(A))
      then
        return false;
    fi;
    
    I := ElementsFamily(FamilyObj(A))!.ideal;
    
    vertices := VerticesOfQuiver(Q);
    paths := [];
    for path in Q do
      if LengthOfPath(path) > 0 then
        Add(paths, path);
      fi;
    od;
    noofverts := Length(vertices);
    noofpaths := Length(paths);
    eA := [];
    for v in [1..noofverts] do
      eA[v] := []; # here will be all paths starting from v
      for p in [1..noofpaths] do
        pp := vertices[v]*paths[p];
        if pp <> Zero(Q) then
          Add(eA[v], pp);
        fi;  
      od;
    od;
    
    for v in [1..noofverts] do
      for vt in [1..noofverts] do
        # checking all paths from v to vt if they belong to I
        inI := 0;
        notinI := 0;
        for pp in eA[v] do # now pp is a path starting from v
          pp := pp*vertices[vt];
          if pp <> Zero(Q) then
            if ElementOfPathAlgebra(PA, pp) in I then
              inI := inI + 1;
            else
              notinI := notinI + 1;
            fi;
            if (inI > 0) and (notinI > 0) then
              return false;
            fi;
          fi;
        od;
      od;
    od;
    
    return true;
end
); # IsSemicommutativeAlgebra for IsQuotientOfPathAlgebra

#######################################################################
##
#P  IsDistributiveAlgebra( <A> ) 
##  
##  This function returns true if the algebra  <A>  is finite 
##  dimensional and distributive. Otherwise it returns false.
##
InstallMethod ( IsDistributiveAlgebra, 
    "for an QuotientOfPathAlgebra",
    true,
    [ IsQuotientOfPathAlgebra ],
    0,
    function( A )

    local pids, localrings, radicalseries, uniserialtest, radlocalrings, 
          i, j, module, testspace, flag, radtestspace;
    
    if not IsFiniteDimensional(A) then 
        Error("the entered algebra is not finite dimensional,\n");
    fi;
    #  Finding a complete set of primitive idempotents in  A.
    pids := List(VerticesOfQuiver(QuiverOfPathAlgebra(A)), v -> v*One(A));
    #  
    #  For each primitive idempotent e, compute eAe and check if 
    #  eAe is a uniserial algebra for all e.
    localrings := List(pids, e -> Subalgebra(A,e*BasisVectors(Basis(A))*e));  
    # 
    #  Check if all algebras in  localrings  are unisersial.
    #
    radicalseries := List(localrings, R -> RadicalSeriesOfAlgebra(R));
    uniserialtest := Flat(List(radicalseries, series -> List([1..Length(series) - 1], i -> Dimension(series[i]) - Dimension(series[i+1]))));
    if not ForAll(uniserialtest, x -> x = 1) then 
        Error("some algebra  eAe  is not a uniserial algebra,\n");
    fi;
    #  Check if the eAe-module eAf is uniserial or 
    #  the fAf-module eAf is uniserial for all pair
    #  of primitive idempotents e and f. 
    radlocalrings := List(localrings, R -> BasisVectors(Basis(RadicalOfAlgebra(R))));
    for i in [1..Length(pids)] do
        for j in [1..Length(pids)] do 
            if i <> j then 
                module := LeftAlgebraModule(localrings[i],\*,Subspace(A,pids[i]*BasisVectors(Basis(A))*pids[j]));
                # compute the radical series of this module over localrings[i] and 
                # check if all layers are one dimensional.
                testspace := ShallowCopy(BasisVectors(Basis(module))); 
                flag := true;
                while Length(testspace) <> 0 and flag do
                    if Length(radlocalrings[i]) = 0 then
                        radtestspace := [];
                    else
                        radtestspace := Filtered(Flat(List(testspace, t -> radlocalrings[i]*t)), x -> x <> Zero(x));
                    fi;
                    if Length(radtestspace) = 0 then
                        if Length(testspace) <> 1 then 
                            flag := false;
                        else
                            testspace := radtestspace;
                        fi;
                    else
                        radtestspace := Subspace(module, radtestspace); 
                        if Length(testspace) - Dimension(radtestspace) <> 1 then
                            flag := false;
                        else
                            testspace := ShallowCopy(BasisVectors(Basis(radtestspace)));
                        fi;
                    fi;
                od;
                if not flag then 
                    module := RightAlgebraModule(localrings[j],\*,Subspace(A,pids[i]*BasisVectors(Basis(A))*pids[j]));
                    # compute the radical series of this module over localrings[j] and
                    # check if all layers are one dimensional.
                    testspace := ShallowCopy(BasisVectors(Basis(module))); 
                    while Length(testspace) <> 0 do
                        if Length(radlocalrings[i]) = 0 then
                            radtestspace := [];
                        else
                            radtestspace := Filtered(Flat(List(testspace, t -> t*radlocalrings[j])), x -> x <> Zero(x));
                        fi;
                        if Length(radtestspace) = 0 then
                            if Length(testspace) <> 1 then 
                                return false;
                            else
                                testspace := radtestspace;
                            fi;
                        else
                            radtestspace := Subspace(module, radtestspace); 
                            if Length(testspace) - Dimension(radtestspace) <> 1 then
                                return false;
                            else
                                testspace := ShallowCopy(BasisVectors(Basis(radtestspace)));
                            fi;
                        fi;
                    od;
                fi;
            fi;
        od;
    od;
    
    return true;
end 
);          

InstallOtherMethod ( IsDistributiveAlgebra, 
    "for an PathAlgebra",
    true,
    [ IsPathAlgebra ],
    0,
    function( A );
    
    return IsTreeQuiver(QuiverOfPathAlgebra(A));
end
);

#######################################################################
##
#A  NakayamaPermutation( <A> )
##
##  Checks if the entered algebra is selfinjective, and returns false
##  otherwise. When the algebra is selfinjective, then it returns a 
##  list of two elements, where the first is the Nakayama permutation 
##  on the simple modules, while the second is the Nakayama permutation
##  on the indexing set of the simple modules. 
## 
InstallMethod( NakayamaPermutation, 
    "for an algebra",
    [ IsQuotientOfPathAlgebra ], 0,
    function( A )

    local perm, nakayamaperm, nakayamaperm_index;
    
    if not IsSelfinjectiveAlgebra(A) then 
        return false;
    fi;
    perm := List(IndecProjectiveModules(A), p -> SocleOfModule(p));
    nakayamaperm := function( x )
        local dimvector, pos;
        
        dimvector := DimensionVector(x);
        if Dimension(x) <> 1 then
            Error("not a simple module was entered as an argument for the Nakayama permutation,\n");
        fi;
        pos := Position(dimvector,1);
        return perm[pos];
    end;
    nakayamaperm_index := function( x )
        local pos;
        
        if not x in [1..Length(perm)] then
            Error("an incorrect value was entered as an argument for the Nakayama permutation,\n");
        fi;
        
        return Position(DimensionVector(perm[x]),1);
    end;
    return [nakayamaperm, nakayamaperm_index];
end
);

InstallOtherMethod( NakayamaPermutation, 
    "for an algebra",
    [ IsPathAlgebra ], 0,
    function( A )

    local n, nakayamaperm, nakayamaperm_index;
    
    if not IsSelfinjectiveAlgebra(A) then 
        return false;
    fi;
    n := NumberOfVertices(QuiverOfPathAlgebra(A));
    nakayamaperm := function( x )
        local dimvector, pos;
        
        dimvector := DimensionVector(x);
        if Dimension(x) <> 1 then
            Error("not a simple module was entered as an argument for the Nakayama permutation,\n");
        fi;
        return x;
    end;
    nakayamaperm_index := function( x )
        local pos;
        
        if not x in [1..n] then
            Error("an incorrect value was entered as an argument for the Nakayama permutation,\n");
        fi;
        
        return x;
    end;
    
    return [nakayamaperm, nakayamaperm_index];
end
);
#######################################################################
##
#A  NakayamaAutomorphism( <A> )
##
##  Checks if the entered algebra is selfinjective, and returns false
##  otherwise. When the algebra is selfinjective, then it returns the 
##  Nakayama automorphism of  <A>. 
##
InstallMethod( NakayamaAutomorphism, 
    "for an algebra",
    [ IsQuotientOfPathAlgebra ], 0,
    function( A )

    local ABasis, arrows, socleequations, soclesolutions, ABasisVector, temp,
          newspan, newABasis, V, pi_map, P, beta, nakaauto;
    
    if not IsSelfinjectiveAlgebra(A) then 
        return false;
    fi;
    ABasis := CanonicalBasis(A); 
    arrows := List(ArrowsOfQuiver(QuiverOfPathAlgebra(A)), x -> x*One(A));
    #
    # Finding the solce of the algebra A
    #
    socleequations := List(ABasis, b -> List(arrows, a -> Coefficients(ABasis, b*a)));
    socleequations := List(socleequations, s -> Flat(s));
    soclesolutions := NullspaceMat(socleequations);
    #
    # Finding a new basis for the algebra  A, where the first basis vectors are 
    # a basis for the socle of the algebra  A.
    #
    ABasisVector := List(ABasis, b -> Coefficients(ABasis, b));
    temp := BaseSteinitzVectors(ABasisVector, soclesolutions);     
    newspan := ShallowCopy(temp!.subspace);
    Append(newspan, temp!.factorspace);
    newABasis := List(newspan, b -> LinearCombination(ABasis, b));
    #
    # Redefining the vector space  A  to have the basis found above. 
    #
    V := Subspace(A, newABasis, "basis"); 
    #
    # Finding the linear map  \pi\colon A \to k being the sum of the dual 
    # basis elements corresponding to a basis of the socle of  A. 
    #
    pi_map := function(x) 
        return Sum(Coefficients(Basis(V), x){[1..Length(soclesolutions)]});
    end;
    #
    # Findind a matrix  P  such that  "x"*P*"y"^T = \pi(x*y), where "x" and
    # "y" means x and y in terms of the basis Basis(V). 
    #
    P := [];
    for beta in newABasis do
        Add(P,List(newABasis, b -> pi_map(beta*b)));
    od;
    # 
    # Since "a"*P*"y"^T = y*P*"\mu(a)"^T = "\mu(a)"*P^T*"y"^T, it follows 
    # that "a"*P = "\mu(a)"*P^T and therefore "\mu(a)" = "a"*P*P^(-T), where
    # \mu is the Nakayama automorphism. Hence, it is given by the following:
    # 
    nakaauto := function(x) 
        local tempo;
        
        if not x in A then 
            Error("the entered argument is not in the algebra,\n");
        fi;
        tempo := Coefficients(Basis(V), x)*P*(TransposedMat(P)^(-1));
        
        return LinearCombination(Basis(V), tempo);
    end;
        
    return nakaauto;
end
);

InstallOtherMethod( NakayamaAutomorphism, 
    "for an algebra",
    [ IsPathAlgebra ], 0,
    function( A )

    local nakaauto;
    
    if not IsSelfinjectiveAlgebra(A) then 
        return false;
    fi;
    nakaauto := function(x);
        if not x in A then 
            Error("argument not in the algebra,\n");
        fi;
        
        return x;
    end;
        
    return nakaauto; 
end
);

#######################################################################
##
#A  OrderOfNakayamaAutomorphism( <A> )
##
##  Checks if the entered algebra is selfinjective, and returns false
##  otherwise. When the algebra is selfinjective, then it returns a 
##  list of two elements, where the first is the Nakayama permutation 
##  on the simple modules, while the second is the Nakayama permutation
##  on the indexing set of the simple modules. 
##
InstallMethod( OrderOfNakayamaAutomorphism, 
    "for an algebra",
    [ IsQuotientOfPathAlgebra ], 0,
    function( A )

    local nakaauto, matrixofnakaauto; 
    
    nakaauto := NakayamaAutomorphism(A); 
    matrixofnakaauto := List(CanonicalBasis(A), b -> 
                             Coefficients(CanonicalBasis(A), nakaauto(b)));
    return Order(matrixofnakaauto);
end
);

InstallOtherMethod( OrderOfNakayamaAutomorphism, 
    "for a path algebra",
    [ IsPathAlgebra ], 0,
    function( A );

    if not IsSelfinjectiveAlgebra(A) then 
        return false;
    fi;

    return 1;
end
);

#######################################################################
##
#O  AssignGeneratorVariables( <A> )
##
##  Takes a quiver algebra  <A>  as an argument and creates variables, 
##  say v_1,...,v_n for the vertices, and a_1,...,a_t for the arrows
##  for the corresponding elements in  <A>, whenever the quiver for 
##  the quiver algebra  <A>  is was constructed with the vertices being 
##  named  v_1,...,v_n and the arrows being named  a_1,...,a_t.  
##  
InstallMethod ( AssignGeneratorVariables, 
    "for a quiver algebra",
    true,
    [ IsQuiverAlgebra ], 
    0,
    function( A )

    local Q, gens, g, s;
    
    Q := QuiverOfPathAlgebra(A); 
    gens := GeneratorsOfQuiver(Q); 
    
    for g in gens do
        s := String(g);
        if not IsValidIdentifier(s) then
            Error("Variable `", s, "' would not be a proper identifier");
        fi;
        if IS_READ_ONLY_GLOBAL(s) then
            Error("Variable `", s, "' is write protected.");
        fi;
    od;
    
    for g in gens do
        s := String(g);
        if ISBOUND_GLOBAL(s) then
            Info(InfoWarning + InfoGlobal, 1, "Global variable `", s,
                 "' is already defined and will be overwritten");
        fi;
        UNBIND_GLOBAL(s);
        ASS_GVAR(s, One(A)*g);
    od;
    Info(InfoWarning + InfoGlobal, 1, "Assigned the global variables ", gens);        
end 
);

#######################################################################
##
#O  AssociatedMonomialAlgebra( <M> )
##
##  Takes as an argument a quiver algebra  <A>  and returns the 
##  associated monomial algebra by using the Groebner basis  <A>  is 
##  endoved with and in particular the ordering of the vertices and the
##  arrows.  Taking another ordering of the vertices and the arrows
##  might changethe associated algebra. 
##  
InstallMethod(AssociatedMonomialAlgebra, 
    "for a finite dimensional quiver algebra",
    true,
    [ IsQuiverAlgebra ], 
    0,
    function( A )

    local fam, gb, relations; 

    if IsPathAlgebra(A) then 
        return A;
    fi;
    fam := ElementsFamily(FamilyObj(A));
    if IsMonomialIdeal(fam!.ideal) then
        return A;
    fi;
    gb := GroebnerBasisOfIdeal(fam!.ideal);
    relations := List(gb!.relations, r -> Tip(r));
    
    return OriginalPathAlgebra(A)/relations;
end
  );

#######################################################################
##
#A  IsMonomialAlgebra( <A> )
##
##  Returns true if the quiver algebra  <A>  is a monomial algebra and false otherwise.
##  
InstallMethod ( IsMonomialAlgebra, 
      "for a PathAlgebraMatModule",
      true,
      [ IsQuiverAlgebra ], 
      0,
      function( A )

      local I;
      
      if IsPathAlgebra(A) then 
          return true;
      else 
          I := ElementsFamily(FamilyObj(A))!.ideal;
          return IsMonomialIdeal(I);
      fi;
end
  ); 

##########################################################################
##
#P DeclareProperty( "IsSemisimpleAlgebra", [ A ] )
##
## The function is defined for an admissible quotient  <A>  of a path
## algebra, and it returns true if  <A>  is a semisimple algebra and 
## false otherwise.
## 
InstallMethod( IsSemisimpleAlgebra,
    "for an algebra",
    [ IsAdmissibleQuotientOfPathAlgebra ], 5, 
        
    function( A )
    local Q;
    
    Q := QuiverOfPathAlgebra(A);
    if Length(ArrowsOfQuiver(Q)) > 0 then
        return false;
    else
        return true;
    fi;
end 
  );

InstallOtherMethod( IsSemisimpleAlgebra,
    "for an algebra",
    [ IsAlgebra ], 0, 
        
    function( A );
    
    if IsFiniteDimensional(A) then 
        if Dimension(RadicalOfAlgebra(A)) = 0 then
            return true;
        else
            return false;
        fi;
    else
        TryNextMethod();
    fi;
end 
  );