# GAP Implementation
# This file was generated from 
# $Id: algebra.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallMethod( PowerSubalgebraOp,
    "for associative algebras",
    true,
    [IsAlgebra and IsCommutative and IsAssociative, IsPosInt], 0,
    function( A, n )
        local m, p, F, b, d, i, j, q, v, gamma, beta, fsBasis;

        F := LeftActingDomain(A);
        p := Characteristic(F);
        m := DegreeOverPrimeField(F);

        if not IsPosInt(m/n) then
            Error("<d> must be a factor of field degree");
        fi;

        q := p^n;
        d := Dimension(A);
        b := CanonicalBasis(A);
        v := BasisVectors(b);

        gamma := [];
        for j in [1..d] do
            Add(gamma, Coefficients(b, v[j] - v[j]^q));
        od;
        beta := NullspaceMat(gamma);

        fsBasis := [];
        for i in beta do
            Add(fsBasis, LinearCombination(b, i));
        od;

        return SubalgebraNC( A, fsBasis, "basis" );
    end );
InstallMethod( PowerSubalgebraOp,
    "for associative algebras",
    true,
    [IsAlgebra, IsPosInt], 0,
    function( A, n )
        local center;

        if not IsAssociative(A) then
            TryNextMethod();
        fi;

        center := Center(A);
        # Below is to work around problems.
        IsCommutative(center); IsAssociative(center);
        return PowerSubalgebraOp( center, n );
    end );
InstallMethod( OppositeAlgebra,
    "for algebras",
    true,
    [ IsAlgebra ], 0,
    function(A)
        local type, K, Aop, gens, fam;

        if HasIsCommutative(A) and IsCommutative(A) then
            return A;
        fi;

        type := NewType( NewFamily( "OppositeAlgebraElementsFamily",
                                    IsOppositeAlgebraElement ),
                         IsPackedElementDefaultRep );
        K := LeftActingDomain(A);

        if IsAlgebraWithOne(A) then
            gens := List( GeneratorsOfAlgebraWithOne(A), 
                          x -> Objectify( type, [x] ) );
            Aop := AlgebraWithOneByGenerators(K, gens);
        else
            gens := List( GeneratorsOfAlgebra(A),
                          x -> Objectify( type, [x] ) );
            Aop := AlgebraByGenerators(K, gens);
        fi;
        SetIsOppositeAlgebra(Aop, true);
        SetUnderlyingAlgebra(Aop, A);

        fam := ElementsFamily( FamilyObj( Aop ) );
        fam!.packedType := type;
        fam!.underlyingAlgebraEltsFam := ElementsFamily( FamilyObj( A ) );
        if IsAlgebraWithOne(Aop) then
            SetOne(fam, Objectify(type, [One(A)]));
        fi;

        return Aop;
    end );
InstallMethod(\in,
    "for opposite algebra elements and opposite algebras",
    IsElmsColls,
    [IsOppositeAlgebraElement, IsOppositeAlgebra], 0,
    function( a, A )
        return a![1] in UnderlyingAlgebra(A);
    end );
InstallMethod( PrintObj,
    "for opposite algebras",
    true,
    [ IsOppositeAlgebra ], 0,
    function(A)
        Print( "OppositeAlgebra( ", A, " )" );
    end );
InstallMethod(Dimension,
    "for opposite algebras",
    true,
    [ IsOppositeAlgebra ], 0,
    function(A)
        return Dimension(UnderlyingAlgebra(A));
    end );
InstallMethod(IsFiniteDimensional,
    "for opposite algebras",
    true,
    [ IsOppositeAlgebra ], 0,
    function(A)
        return IsFiniteDimensional(UnderlyingAlgebra(A));
    end );
DeclareRepresentation("IsBasisOfOppositeAlgebraDefaultRep",
    IsComponentObjectRep, ["underlyingBasis"] );
InstallMethod(Basis,
    "for opposite algebras",
    true,
    [ IsOppositeAlgebra ], 0,
    function(A)
        local B;
        
        B := Objectify( NewType( FamilyObj( A ),
                                 IsBasis and
                                 IsBasisOfOppositeAlgebraDefaultRep and 
                                 IsAttributeStoringRep ),
                        rec() );
        SetUnderlyingLeftModule(B, A);
        B!.underlyingBasis := Basis(UnderlyingAlgebra(A));
        return B;
    end);
InstallMethod(CanonicalBasis,
    "for opposite algebras",
    true,
    [ IsOppositeAlgebra ], 0,
    function(A)
        local B;
        
        B := Objectify( NewType( FamilyObj( A ),
                                 IsBasis and
                                 IsBasisOfOppositeAlgebraDefaultRep and 
                                 IsAttributeStoringRep ),
                        rec() );
        SetUnderlyingLeftModule(B, A);
        SetIsCanonicalBasis(B, true);
        B!.underlyingBasis := CanonicalBasis(UnderlyingAlgebra(A));
        if B!.underlyingBasis = fail then
            return fail;
        fi;
        return B;
    end);

InstallMethod(Coefficients,
    "for bases of opposite algebras and opposite algebra element",
    true,
    [ IsBasisOfOppositeAlgebraDefaultRep and IsBasis,
      IsOppositeAlgebraElement and IsPackedElementDefaultRep ], 0,
    function(B, a)
        return Coefficients(B!.underlyingBasis, a![1]);
    end );

InstallMethod(BasisVectors,
    "for bases of opposite algebras",
    true,
    [ IsBasisOfOppositeAlgebraDefaultRep and IsBasis ], 0,
    function(B)
        local vectors, type;

        vectors := BasisVectors( B!.underlyingBasis );
        type := ElementsFamily(FamilyObj(B))!.packedType;
        return List(vectors, x -> Objectify(type, [x]));
    end );
InstallMethod( \.,
    "for path algebras",
    true,
    [IsOppositeAlgebra, IsPosInt], 0,
    function(A, name)
        local family;
 
        family := ElementsFamily(FamilyObj(A));
        return Objectify(family!.packedType, 
                         [ UnderlyingAlgebra(A).(NameRNam(name)) ] );
    end );
InstallMethod( \=,
    "for two elements of an opposite algebra",
    IsIdenticalObj,
    [IsOppositeAlgebraElement and IsPackedElementDefaultRep,
     IsOppositeAlgebraElement and IsPackedElementDefaultRep], 0,
    function( a, b )
        return a![1] = b![1];
    end );
InstallMethod( \<,
    "for two elements of an opposite algebra",
    IsIdenticalObj,
    [IsOppositeAlgebraElement and IsPackedElementDefaultRep,
     IsOppositeAlgebraElement and IsPackedElementDefaultRep], 0,
    function( a, b )
        return a![1] < b![1];
    end );
InstallMethod( \+,
    "for two elements of an opposite algebra",
    IsIdenticalObj,
    [IsOppositeAlgebraElement and IsPackedElementDefaultRep,
     IsOppositeAlgebraElement and IsPackedElementDefaultRep], 0,
    function( a, b )
        return Objectify(TypeObj(a), [a![1] + b![1]]);
    end );
InstallMethod( \*,
    "for two elements of an opposite algebra",
    IsIdenticalObj,
    [IsOppositeAlgebraElement and IsPackedElementDefaultRep,
     IsOppositeAlgebraElement and IsPackedElementDefaultRep], 0,
    function( a, b )
        return Objectify(TypeObj(a), [b![1] * a![1]]);
    end );
InstallMethod( \*,
    "for scalar and an opposite algebra element",
    true,
    [IsScalar,
     IsOppositeAlgebraElement and IsPackedElementDefaultRep], 0,
    function( a, b )
        return Objectify(TypeObj(b), [b![1]*a]);
    end );
InstallMethod( \*,
    "for an opposite algebra element and a scalar",
    true,
    [ IsOppositeAlgebraElement and IsPackedElementDefaultRep,
      IsScalar ], 0,
    function( a, b )
        return Objectify(TypeObj(a), [b*a![1]]);
    end );
InstallMethod( AINV,
    "for an opposite algebra element",
    true,
    [ IsOppositeAlgebraElement and IsPackedElementDefaultRep ], 0,
    function(a)
        return Objectify(TypeObj(a), [-a![1]]);
    end );
InstallMethod( ZeroOp,
    "for an opposite algebra element",
    true,
    [ IsOppositeAlgebraElement and IsPackedElementDefaultRep ], 0,
    function(a)
        return Objectify(TypeObj(a), [0*a![1]]);
    end );
InstallMethod(PrintObj,
    "for an opposite algebra element",
    true,
    [ IsOppositeAlgebraElement and IsPackedElementDefaultRep ], 0,
    function(a)
        Print( "op( ", a![1], " )" );
    end );
InstallMethod(ExtRepOfObj,
    "for opposite algebra elements",
    true,
    [ IsOppositeAlgebraElement and IsPackedElementDefaultRep ], 0,
    function(a)
        return ExtRepOfObj(a![1]);
    end );
InstallMethod(ObjByExtRep,
    "for opposite algebra family, object",
    true,
    [ IsOppositeAlgebraElementFamily, IsObject ], 0,
    function(fam, obj)
        if fam!.underlyingAlgebraEltsFam <> FamilyObj(obj) then
            TryNextMethod();
        fi;
        return Objectify( fam!.packedType, [obj] );
    end );
#############################################################################
##
#M  AsLeftModuleGeneralMappingByImages( <alg_gen_map> )
##
##  If necessary then we compute a basis of the preimage,
##  and images of its basis vectors.
##
##  Note that we must prescribe also the products of basis vectors and
##  their images if <alg_gen_map> is not known to be a mapping.
##
InstallMethod( AsLeftModuleGeneralMappingByImages,
    "for an algebra general mapping by images",
    true,
    [     IsAlgebraGeneralMapping
      and IsAlgebraGeneralMappingByImagesDefaultRep ], 0,
    function( alg_gen_map )

    local origgenerators,  # list of algebra generators of the preimage
          origgenimages,   # list of images of `origgenerators'
          generators,      # list of left module generators of the preimage
          genimages,       # list of images of `generators'
          A,               # source of the general mapping
          left,            # is it necessary to multiply also from the left?
                           # (not if `A' is associative or a Lie algebra)
          maxdim,          # upper bound on the dimension
          MB,              # mutable basis of the preimage
          dim,             # dimension of the actual left module
          len,             # number of algebra generators
          i, j,            # loop variables
          gen,             # loop over generators
          prod;            #

    origgenerators := alg_gen_map!.generators;
    origgenimages  := alg_gen_map!.genimages;

    if IsBasis( origgenerators ) then

      generators := origgenerators;
      genimages  := origgenimages;

    else

      generators := ShallowCopy( origgenerators );
      genimages  := ShallowCopy( origgenimages );

      A:= Source( alg_gen_map );

      left:= not (    ( HasIsAssociative( A ) and IsAssociative( A ) )
                   or ( HasIsLieAlgebra( A ) and IsLieAlgebra( A ) ) );

      if HasDimension( A ) then
        maxdim:= Dimension( A );
      else
        maxdim:= infinity;
      fi;

      # $A_1$
      MB:= MutableBasis( LeftActingDomain( A ), generators,
                                     Zero( A ) );
      dim:= 0;
      len:= Length( origgenerators );

      while dim < NrBasisVectors( MB ) and dim < maxdim do

        # `MB' is a mutable basis of $A_i$.
        dim:= NrBasisVectors( MB );

        # Compute $\bigcup_{g \in S} ( A_i g \cup A_i g )$.
        for i in [ 1 .. len ] do
          gen:= origgenerators[i];
          for j in [ 1 .. Length( generators ) ] do
            prod:= generators[j] * gen;
            if not IsContainedInSpan( MB, prod ) then
              Add( generators, prod );
              Add( genimages, genimages[j] * origgenimages[i] );
              CloseMutableBasis( MB, prod );
            fi;
          od;
        od;

        if left then

          # Compute $\bigcup_{g \in S} ( A_i g \cup g A_i )$.
          for i in [ 1 .. len ] do
            gen:= origgenerators[i];
            for j in [ 1 .. Length( generators ) ] do
              prod:= gen * generators[j];
              if not IsContainedInSpan( MB, prod ) then
                Add( generators, prod );
                Add( genimages, origgenimages[i] * genimages[j] );
                CloseMutableBasis( MB, prod );
              fi;
            od;
          od;

        fi;

      od;

    fi;

    # Is this the proper fix?
    len := Length(generators);

    # Add the products of basis vectors to `generators',
    # and the products of their images to `genimages'.
    for i in [ 1 .. len ] do
      for j in [ 1 .. len ] do
        Add( generators, generators[i] * generators[j] );
        Add( genimages, genimages[i] * genimages[j] );
      od;
    od;

    # Construct and return the left module general mapping.
    return LeftModuleGeneralMappingByImages( A, Range( alg_gen_map ),
               generators, genimages );
    end );
InstallMethod( RadicalOfAlgebra,
    "for an algebra of l.m.b.m rep.",
    true,
    [ IsAlgebraWithOne and IsGeneralMappingCollection ], 0,
    function( A )
        local gens, uA, F, Igens, R, S, I;

        gens := GeneratorsOfAlgebraWithOne(A);
        if not IsLinearMappingByMatrixDefaultRep(gens[1]) then
            TryNextMethod();
        fi;

        S := Source(gens[1]);
        R := Range(gens[1]);
        F := LeftActingDomain(A);
        if HasBasis(A) and HasBasisVectors(Basis(A)) then
            gens := List(BasisVectors(Basis(A)), x -> TransposedMat(x!.matrix));
            uA := AlgebraWithOneByGenerators(F, gens, "basis");
        else
            gens := List(gens, x -> TransposedMat(x!.matrix));
            uA := AlgebraWithOneByGenerators(F, gens);
        fi;
        I := RadicalOfAlgebra(uA);

        Igens := BasisVectors(Basis(I));
        Igens := List(Igens, x -> TransposedMat(x));
        Igens := List(Igens, 
                      x->LeftModuleHomomorphismByMatrix(Basis(S),x,Basis(R)));

        I := TwoSidedIdealNC(A, Igens, "basis");
        SetNiceFreeLeftModuleInfo( I, NiceFreeLeftModuleInfo(A) ); 
        return I;
    end );
    
#######################################################################
##
#A  RadicalSeriesOfAlgebra( <A> ) 
##  
##  This function returns the radical series of the algebra  <A> in a 
##  list, where the first element is the algebra  <A> itself, then 
##  radical of <A>, radical square of <A>, and so on.
##
InstallMethod( RadicalSeriesOfAlgebra,
    "for an algebra",
    [ IsAlgebra ],
    function ( A )
    local   radical,    # radical of the algebra <A>
            S,          # radical series of the algebra <A>, result
            D;          # power of the radical

    radical := RadicalOfAlgebra(A);
    # Compute the series by repeated calling of `ProductSpace'.
    S := [ A ];
    D := radical;
    while Dimension(D) <> 0  do
      Add( S, D );
      D := ProductSpace( D, radical );
    od;
    Add(S,D); 
    
    # Return the series when it becomes zero.
    return S;
end 
);

