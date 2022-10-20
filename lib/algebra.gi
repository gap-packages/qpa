# GAP Implementation
# This file was generated from 
# $Id: algebra.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallImmediateMethod( GlobalDimension, IsSelfinjectiveAlgebra, 0, 
        function(A) 
    if HasIsSemisimpleAlgebra(A) and not IsSemisimpleAlgebra(A) then 
        return infinity;
    fi;
    TryNextMethod();
end);

InstallImmediateMethod( GlobalDimension, IsSemisimpleAlgebra, 0, A -> 0); 

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
end 
);
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
InstallMethod( AdditiveInverseSameMutability,
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

#######################################################################
##
#P  IsRadicalSquareZeroAlgebra( <A> ) 
##  
##  This function returns true if the algebra has the property that 
##  the radical squares to zero. Otherwise it returns false.
##
InstallMethod( IsRadicalSquareZeroAlgebra,
    "for an algebra",
    [ IsAlgebra ],
    function ( A )
    
    local radical; # radical of the algebra <A>
    
    radical := RadicalOfAlgebra(A);
    
    return Dimension(ProductSpace(radical,radical)) = 0;
end
  );

#######################################################################
##
#P  IsBasicAlgebra( <A> )
##
##  This function returns true if the entered algebra  <A>  is a (finite
##  dimensional) basic algebra and false otherwise. This method applies 
##  to algebras over finite fields. 
##
InstallMethod( IsBasicAlgebra,
    "for an algebra",
    [ IsAlgebra ], 0,
    function( A )

    local K, J, L;
    #
    # Only finite dimensional algebras are regarded as basic.
    # 
    if not IsFiniteDimensional(A) then
        Error("the entered algebra is not finite dimensional,\n");
    fi;
    K := LeftActingDomain(A);
    #
    # Here we only can deal with algebras over finite fields.
    # First we find a decomposition of the algebra modulo the
    # radical, and then we find all the primitive idempotents 
    # in each block. If each block only contains one primitive
    # idempotent, then the algebra is basic. 
    #
    if IsFinite(K) then 
        J := RadicalOfAlgebra(A);
        L := List(DirectSumDecomposition(A/J),PrimitiveIdempotents);
        L := List(L,Length);
        if Sum(L) <> Length(L) then
            return false;
        else
            return true;
        fi;
    else
        TryNextMethod();
    fi;
end
);

#######################################################################
##
#P  IsElementaryAlgebra( <A> )
##
##  This function returns true if the entered algebra  <A>  is a (finite
##  dimensional) elementary algebra and false otherwise. This method 
##  applies to algebras over finite fields. 
##
InstallMethod( IsElementaryAlgebra,
    "for an algebra",
    [ IsAlgebra ], 0,
    function( A )

    local K, J, D, L;
    
    K := LeftActingDomain(A);
    #
    # Here we only can deal with algebras over finite fields.
    # First we find a decomposition of the algebra modulo the
    # radical, and then we find all the primitive idempotents 
    # in each block. If each block only contains one primitive
    # idempotent, then the algebra is basic. 
    #
    if IsFinite(K) and IsBasicAlgebra(A) then 
        J := RadicalOfAlgebra(A);
        D := DirectSumDecomposition(A/J);
                # Since  A  is basic, we know that  D  only consists
                # of simple blocks of division algebras. Hence  A  is 
                # elementary if all the blocks have the same dimension 
                # over  K, as then they are all the same finite field.
                #
        L := List(D, Dimension); 
        if ForAll(L, x -> x = L[1]) then
            return true;
        else
            return false;
        fi;
    else
        TryNextMethod();
    fi;
end
  );

#######################################################################
##
#O  LiftingIdempotent( <f>, <e> )
##
##  When the domain of  <f>  is finite dimensional, then it checks if
##  the element  <e>  is an idempotent.  If so and  <e>  has a preimage,
##  then this operation returns a lifting of  <e>  to the domain of  <f>
##  whenever possible.  It returns  fail  if  the element  <e>  has no 
##  preimage. If the algorithm is unable to construct a lifting, it 
##  returns  false.  Using the algorithm described in the proof of 
##  Proposition 27.1 in Anderson and Fuller, Rings and categories of 
##  modules, second edition, GMT, Springer-Verlag.
##  
InstallMethod( LiftingIdempotent,  
    "for two PathAlgebraMatModule's",
    [ IsAlgebraGeneralMapping, IsObject ], 0,
    function( f, e )

    local elift, h, A, i, t, k;
    #
    # Is the input OK?
    # Checking if the domain of the map  <f>  is a finite dimensional algebra.
    #
    if not IsFiniteDimensional(Source(f)) then
        Error("the domain of <f> is not finite dimensional, ");
    fi;
    #
    # Checking if the entered element  <e>  is an idempotent.
    #
    if e <> e^2 then 
        Error("the entered element  <e>  is not an idempotent, ");
    fi; 
    #
    # Checking if the element  <e>  is in the range of the map  <f>, and next
    # if the element  <e>  is in the image of the map  <f>.
    #
    if not ( e in Range(f) ) then
        Error("the entered element <e> is not in the range of the entered homomorphism <f>, ");
    fi;
    elift := PreImagesRepresentative(f, e);
    if elift = fail then
        Error("the enter element <e> has not preimage in the domain of <f>, ");
    fi;
    #
    # Starting the algorithm described in the proof of Proposition 27.1 in Anderson & Fuller.
    #
    h := elift - elift^2;
    A := Source(f);
    if h = Zero(A) then
        return elift;
    fi;
    i := 0;
    repeat
        i := i + 1;
    until
        Dimension(Subspace(A, BasisVectors(Basis(A))*h^i)) = Dimension(Subspace(A, BasisVectors(Basis(A))*h^(i+1)));
    if h^i = Zero(A) then
        t := Binomial(i,1)*MultiplicativeNeutralElement(A);
        for k in [2..i] do 
            t := t + (-1)^(k-1)*Binomial(i,k)*elift^(k-1);
        od;
        return elift^i*t^i;
    else
        return false;
    fi;
end
  );

#######################################################################
##
#O  LiftingCompleteSetOfOrthogonalIdempotents( <map>, <idempotents> )
##
##  Given a map  <map> :  A --> B  and a complete set  <idempotents>
##  of orthogonal idempotents in  B, which all are in the image
##  of  <map>,  this function computes (when possible) a complete set of
##  orthogonal idempotents of preimages of the idempotents in  B. 
##
InstallMethod(LiftingCompleteSetOfOrthogonalIdempotents,
    "for an algebra mapping and list of idempotents",
    true,
    [IsAlgebraGeneralMapping,IsHomogeneousList],
    0,
    function(map, idempotents)

    local n, i, j, idems, A, liftidem, temp;
    
    #
    # Input OK?
    # Checking if the domain of  <map>  is a finite dimensional algebra.
    #
    if not IsFiniteDimensional(Source(map)) then
        Error("the domain of the entered  <map>  is not finite dimensional, \n");
    fi;
    #
    # Checking if the entered elements are a complete set of orthogonal
    # idempotents. 
    #
    n := Length(idempotents); 
    for i in [1..n] do
        for j in [1..n] do
            if i = j then
                if idempotents[i]^2 <> idempotents[i] then
                    Error("one of the entered elements is not an idempotent, ");
                fi;
            else
                if idempotents[i]*idempotents[j] <> Zero(Range(map)) then
                    Error("the enter elements are not orthogonal, ");
                fi;
            fi;
        od;
    od;
    if Sum(idempotents) <> MultiplicativeNeutralElement(Range(map)) then
        Error("the entered elements are not a complete set of idempotents, ");
    fi;
    #
    # Lift the idempotents in the range of  <map>  to idempotents in the 
    # domain of  <map>.
    #
    liftidem := [];
    Add(liftidem, LiftIdempotent(map, idempotents[1]));
    for i in [1..Length(idempotents) - 1] do
        temp := LiftTwoOrthogonalIdempotents(map, Sum(liftidem{[1..i]}), idempotents[i+1]);
        Add(liftidem, temp[2]);
    od;

    return liftidem;
end
);

#######################################################################
##
#O  FindMultiplicativeIdentity( <A> )
##
##  This function finds the injective envelope I(M) of the module  <M>  
##  in that it returns the map from M ---> I(M). 
##
InstallMethod ( FindMultiplicativeIdentity, 
    "for a finite dimensional algebra",
    true,
    [ IsAlgebra ],
    0,
    function( A )
    
    local   B,  constantterm,  matrix,  i,  b1,  partialmatrix,  b2,  
            totalmatrix,  totalconstantterm,  identity;
    
    if not IsFiniteDimensional( A ) then
        Error("The entered algebra is not finite dimensional,\n");
    fi;
    B := Basis( A ); 
    
    constantterm := [];
    matrix := [];
    i := 1;
    for b1 in B do
        partialmatrix := [];
        for b2 in B do 
            Add( partialmatrix, Coefficients( B, b2*b1 ) );
        od;
        Add( constantterm, Coefficients( B, b1 ) );
        Add( matrix, [1, i, partialmatrix ] );
        i := i + 1;
    od;
    totalmatrix := BlockMatrix( matrix, 1, i - 1 );
    totalconstantterm := Flat(constantterm);
    identity := SolutionMat( totalmatrix, totalconstantterm );
    
    if identity = fail then
        return fail;
    else
        return LinearCombination( B, identity );
    fi;
end
);

InstallMethod ( AlgebraAsQuiverAlgebra, 
    "for a finite dimensional algebra",
    true,
    [ IsAlgebra ],
    0,
    function( A )

  local F, idA, C, id, centralidempotentsinA, radA, g, centralidem, 
        factoralgebradecomp, c, temp, dimcomponents, pi, ppowerr, 
        gens, D, idD, order, generatorinD, K, i, alghom, vertices, 
        radAsquare, h, radAmodsquare, arrows, adjacencymatrix, j, t, 
        Q, KQ, Jt, n, images, AA, AAgens, AAvertices, AAarrows, f, B, 
        matrix, fam, b, tempx, walk, length, image, solutions, 
        idealgens, I;
    #
    # Checking if  <A>  is a finite dimensional algebra over a finite field.
    #
    F := LeftActingDomain(A);
    if not IsFinite(F) then
        Error("the entered algebra is not an algebra over a finite field, \n");
    fi;
    if not IsFiniteDimensional(A) then
        Error("the entered algebra is not finite dimensional, \n");
    fi; 
    if One( A ) = fail then
        idA := FindMultiplicativeIdentity( A );
        SetOne( A, idA );
    fi;
    #
    # Checking if the algebra is indecomposable.
    #
    C := Center(A);
    id := FindMultiplicativeIdentity( C );
    if id = fail then
        Error("The center of the entered algebra does not have a multiplicative identity,\n");
    else
        SetOne( C, id );
    fi;
    centralidempotentsinA := CentralIdempotentsOfAlgebra(C);
    if Length(centralidempotentsinA) > 1 then
        Error("the entered algebra is not indecomposable, \n");
    fi;
    #
    # Checking if the algebra is basic.  Remember all finite division rings are fields.
    #
    radA := RadicalOfAlgebra(A);
    g := NaturalHomomorphismByIdeal(A,radA);
    centralidem := CentralIdempotentsOfAlgebra(Range(g)); 
    factoralgebradecomp := [];
    for c in centralidem do
        temp := FLMLORByGenerators( LeftActingDomain(A), Filtered(c*BasisVectors(Basis(Range(g))), b -> b <> Zero(A)));
        SetParent(temp, A); 
        SetOne(temp, c);
        SetMultiplicativeNeutralElement(temp, c);        
        Add(factoralgebradecomp, temp);
    od;
    if not ForAll(factoralgebradecomp, IsCommutative) then
        Error("the entered algebra is not basic, \n");
    fi;
    #
    # Checking if the algebra is elementary.
    #
    dimcomponents := List(factoralgebradecomp, Dimension);
    if not ForAll(dimcomponents, d -> d = dimcomponents[1]) then 
        Error("the entered algebra is not elementary, hence not a quotient of a path algebra, \n");
    fi;
    #
    # Finding the field  K  over which the algebra is K-elementary. Recall C = Center(A).
    # 
    pi := NaturalHomomorphismByIdeal(C, RadicalOfAlgebra(C)); 
    #
    # Known by the Wedderburn-Malcev «Principal» theorem that A is an algebra over 
    # a field isomorphic to  C/radC. 
    #
    ppowerr := Size(Range(pi));
    gens := Filtered(BasisVectors(Basis(C)), b -> ForAll(centralidem, c -> c*ImageElm(g, b) <> Zero(One(Range(g))))); 
    D := Subalgebra(C, gens);
    idD := FindMultiplicativeIdentity( D );
    if idD = fail then
        Error("The center of the entered algebra does not have a multiplicative identity,\n");
    else
        SetOne( D, idD );
    fi;    
    order := function( obj )
        
        local one, pow, ord;
        
        if obj = Zero(obj) then
            return -1;
        fi;
        one := One( D );
        pow:= obj;
        ord:= 1;
        while pow <> one do
            ord:= ord + 1;
            pow:= pow * obj;
            if pow = Zero(obj) then
                return -1;
            fi;
        od;

        return ord;
    end;
    generatorinD := Filtered(Elements(D), d -> order(d) = ppowerr - 1);     
    K := GF(ppowerr);
    D := AsAlgebra(PrimeField(K), D);
    i := 1;
    repeat 
      alghom := AlgebraHomomorphismByImages( K, D, GeneratorsOfDivisionRing( K ), [ generatorinD[ i ] ] );
      i := i + 1;
    until 
      alghom <> fail or i = Length( generatorinD ) + 1;
    #
    # Finding representatives for the vertices in  A.
    #
    vertices := LiftingCompleteSetOfOrthogonalIdempotents(g, centralidem); 
    #
    # Finding the radical square in  <A> and storing it in  <radAsquare>. 
    #
    radAsquare := ProductSpace(radA,radA);
    if Dimension(radAsquare) = 0 then
        radAsquare := Ideal(A, []);
    else
        radAsquare := Ideal(A, BasisVectors(Basis(radAsquare)));
    fi;
    #
    # Finding the natural homomorphism  <A> ---> <A>/rad^2 <A> and 
    # finding the image of  rad <A>  in  <A>/rad^2 <A>  and storing it in  <radmodsquare>. 
    #
    h := NaturalHomomorphismByIdeal(A, radAsquare);
    radAmodsquare := Ideal(Range(h), List(BasisVectors(Basis(radA)), b -> ImageElm(h,b)));
    #
    # Finding a basis for the arrows for the algebra  <A>  inside  <A>. 
    # At the same time finding the adjacency matrix for the quiver of  <A>. 
    #
    arrows := List([1..Length(centralidem)], x -> List([1..Length(centralidem)], y -> []));
    adjacencymatrix := NullMat(Length(centralidem),Length(centralidem));
    for i in [1..Length(centralidem)] do
        for j in [1..Length(centralidem)] do
            arrows[i][j] := Filtered(ImageElm(h,vertices[i])*BasisVectors(Basis(radAmodsquare))*ImageElm(h,vertices[j]), y -> y <> Zero(y));
            arrows[i][j] := BasisVectors(Basis(Subspace(Range(h),arrows[i][j])));
            arrows[i][j] := List(arrows[i][j], x -> vertices[i]*PreImagesRepresentative(h,x)*vertices[j]);
            adjacencymatrix[i][j] := Length(arrows[i][j]);
        od; 
    od;
    #
    # Defining the quiver of the algebra  <A>  and storing it in  <Q>. 
    #
    t := Length( RadicalSeriesOfAlgebra( A ) ) - 1; # then (rad A)^t = (0) 
    Q := Quiver( adjacencymatrix );
    KQ := PathAlgebra( K, Q );
    Jt := NthPowerOfArrowIdeal( KQ, t );
    n := NumberOfVertices( Q );
    images := ShallowCopy( vertices );   #  images of the vertices/trivial paths
    for i in [ 1 .. n ] do
        for j in [ 1 .. n ] do
            Append( images,arrows[ i ][ j ] );
        od;
    od;
    #
    #  Define  AA := KQ/J^t, where t = the Loewy length of  <A>, and 
    #  in addition define f : AA ---> A. Find this as a linear map and find 
    #  the kernel, and construct the relations from this.
    #
    if Length( Jt ) = 0 then 
        AA := KQ;
        AAgens := GeneratorsOfAlgebra( AA );
        AAvertices := AAgens{ [ 1..n ] };
        AAarrows := AAgens{ [ n + 1..Length( AAgens ) ] };
        f := [ AA, A, AAgens{ [ 1..Length( AAgens ) ] }, images ];
    else
        AA := KQ/Jt;
        AAgens := GeneratorsOfAlgebra( AA );
        AAvertices := AAgens{ [ 2..n + 1 ] };
        AAarrows := AAgens{ [ n + 2..Length( AAgens ) ] };
        f := [ AA, A, AAgens{ [ 2..Length( AAgens ) ] }, images ];
    fi;
    #
    #  First giving the ring surjection  AA ---> A  as a linear map.  Stored
    #  in the matrix called  <matrix>.
    #
    B := BasisVectors( Basis( AA ) ); 
    matrix := []; 
    fam := ElementsFamily( FamilyObj( AA ) ); 
    for b in B do
        if IsPathAlgebra( AA ) then 
            temp := CoefficientsAndMagmaElements( b );
        else 
            temp := CoefficientsAndMagmaElements( b![ 1 ] );
        fi;
        n := Length( temp )/2;
        tempx := Zero( f[ 2 ] );
        for i in [ 0..n - 1 ] do
            # for each term compute the image. 
            walk := WalkOfPath( temp[ 2 * i + 1 ] ); 
            length := Length( walk ); 
            image := One( A ); 
            if length = 0 then 
                if IsPathAlgebra( AA ) then
                    image := image * f[ 4 ][ Position( f[ 3 ], temp[ 2 * i + 1 ] * One( AA ) ) ];
                else
                    image := image*f[4][Position(f[3],ElementOfQuotientOfPathAlgebra(fam, temp[2*i + 1]*One(KQ), true))];
                fi;
            else 
                for j in [ 1..length ] do
                    image := image * f[ 4 ][ Position( f[ 3 ], walk[ j ] * One( AA ) ) ];
                od;
            fi;
            tempx := tempx + ImageElm( alghom, temp[ 2 * i + 2 ] ) * image;
        od;
        Add( matrix, Coefficients( Basis( f[ 2 ] ), tempx ) );
    od;    
    #
    #  Finding a vector space basis for the kernel of the ring surjection  AA ---> A.
    #
    solutions := NullspaceMat( matrix );
    idealgens := List( solutions, x -> LinearCombination( B, x ) );  # solutions as elements in  AA.
    # 
    # Lifting the relations back to KQ.
    #
    if not IsPathAlgebra( AA ) then
        idealgens := List( idealgens, x -> x![ 1 ] );
    fi;    
    if Length( idealgens ) = 0 then 
        return [ AA, images ];
    else
        return [ KQ/idealgens, images ];
    fi;
end
  );

#######################################################################
##
#O  BlocksOfAlgebra
##
##  This function finds the block decomposition of a finite dimensional
##  algebra, and returns it as a list of algebras. 
##
InstallMethod ( BlocksOfAlgebra, 
    "for a finite dimensional algebra",
    true,
    [ IsAlgebra ],
    0,
    function( A )

  local cents, Alist, n, i;
  
  cents := CentralIdempotentsOfAlgebra( A );
  Alist := [];
  n := Length( cents );
  for i in [ 1..n ] do
    Add( Alist, FLMLORByGenerators( LeftActingDomain( A ), BasisVectors( Basis(A)) * cents[ i ] ) );
    SetParent( Alist[ i ], A );  
    SetOne( Alist[ i ], cents[ i ] );
    SetMultiplicativeNeutralElement( Alist[ i ], cents[ i ] );
    SetFilterObj( Alist[ i ], IsAlgebraWithOne );  
  od;

  return Alist;
end
  );

#######################################################################
##
#A  OppositeAlgebraHomomorphism( <f> )
##
##  Constructs the opposite of an algebra homomorphism. 
## 
InstallMethod( OppositeAlgebraHomomorphism, 
    "for an algebra",
    [ IsAlgebraHomomorphism ], 0,
    function( f )
    
    local   A,  B,  Aop,  Bop,  gensandimages,  gensAop,  imagesBop,  
            g;
    
    A := Source( f ); 
    B := Range( f );
    Aop := OppositeAlgebra( A );
    Bop := OppositeAlgebra( B );
    gensandimages := MappingGeneratorsImages( f );
    gensAop := List( gensandimages[ 1 ], OppositePathAlgebraElement );
    imagesBop := List( gensandimages[ 2 ], OppositePathAlgebraElement );
    g := AlgebraHomomorphismByImages( Aop, Bop, gensAop, imagesBop );
    g!.generators := gensAop;
    g!.genimages := imagesBop;
    
    return g;
end
  );

