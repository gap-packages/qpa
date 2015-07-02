# GAP Implementation
# This file was generated from 
# $Id: idempotent.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
#

#######################################################################
##
#A  PrimitiveIdempotents( <A> )
##
##  This attribute takes as an argument a finite dimensional simple 
##  algebra  <A>  for a finite field and returns a complete set of 
##  primitive idempotents  { e_i }  such that  
##                     A \simeq  Ae_1 + ... + Ae_n.
##
##  This function is based on Eberly, W., "Decomposition of algebras 
##  over finite fields and number fields", Computational Complexity 1 
##  (1991), 183–210. See page 84 in Craig A. Struble, "Analysis and 
##  implementation of algorithms for noncommutative algebra", PhD-
##  thesis, Virginia Tech, 2000. 
##
InstallMethod( PrimitiveIdempotents, 
    "for semisimple algebras",
    true,
    [ IsFiniteDimensional ], 0,
    function(A)
    local F, d, e, eA, r, m, E, b, one, i, s, j, x, e_hat;
    
    
    #
    # Input a (semi-)simple algebra?
    # 
    if not IsSemisimpleAlgebra(A) then
        Error("the entered algebra is not a semisimple algebra, \n");
    fi;
    if Length(CentralIdempotentsOfAlgebra(A)) > 1 then
        Error("the entered algebra is not a simple algebra, \n");
    fi;
    
    F := LeftActingDomain(A);
    if not IsFinite(F) then
        Error("the entered algebra is not over a finite field, \n");
    fi;
    d := Dimension(A);
    e := HOPF_SingularIdempotent(A);
    eA := e*A;

    r := Dimension(eA);  # Dimension of a copy of the center of A. 
    m := d / r;          # Matrix dimension (rows and cols)
    E := List([1..m], x -> Zero(A));
    E[1] := e;
    s := Zero(A);
    one := MultiplicativeNeutralElement(A);
    b := BasisVectors(Basis(A));

    for i in [2..m] do
        s := s + E[i-1];
        j := 0;
        repeat
            j := j+1;
            x := (one - s)*b[j]*e;
        until x <> Zero(A);
        e_hat := HOPF_LeftIdentity(x*A);
        E[i] := e_hat*(one - s);
    od;
    return E;
end
);

#######################################################################
##
#F  HOPF_MinPoly( <A>, <a> )
##
##  This function takes as an argument a finite dimensional algebra  <A>
##  over a field  K  and an element  <a>  in this algebra and returns 
##  the monic polynomial  f  of smallest degree such that  f(a) = 0  in
##  <A>. 
##
InstallGlobalFunction(HOPF_MinPoly,
    function( A, a )
    local F, vv, sp, x, cf, f;

    F := LeftActingDomain(A);
    vv := [ MultiplicativeNeutralElement( A ) ];
    sp := MutableBasis(F, vv);
    x := ShallowCopy(a);
    while not IsContainedInSpan( sp, x ) do
        Add(vv,x);
        CloseMutableBasis(sp, x);
        x := x*a;
    od;
    sp := UnderlyingLeftModule(ImmutableBasis(sp));
    cf := ShallowCopy( - Coefficients(BasisNC(sp,vv), x));
    Add(cf, One(F));
    f := ElementsFamily(FamilyObj(F));
    f := LaurentPolynomialByCoefficients( f, cf, 0 );
    return f;
end
);

#######################################################################
##
#F  HOPF_ZeroDivisor( <A> )
##
##  This function takes as an argument a finite dimensional simple 
##  algebra  <A>  over a finite field and returns a list of two elements  
##  [a, b]  in  <A>  such that  a*b = 0  in <A>, if  <A> is not 
##  commutative. If  <A>  is commutative, it returns an empty list.
##
##  This function is based on Ronyai, L., "Simple algebras are difficult", 
##  In 19th ACM Symposium on Theory of Computing (1987), pp. 398–408, and 
##  Ro ́nyai, L., "Computing the structure of finite algebras", Journal of 
##  Symbolic Computation 9 (1990), 355–373. See page 83 in Craig A. Struble,
##  "Analysis and implementation of algorithms for noncommutative algebra", 
##  PhD-thesis, Virginia Tech, 2000. 
##   
InstallGlobalFunction(HOPF_ZeroDivisor,
  function(A)
    local F, d, b, bv, cf, x, m, facts, one, f, g;
    if IsCommutative(A) then
        return [];
    fi;

    F := LeftActingDomain(A);
    if not (IsFinite(F) and IsFiniteDimensional(A)) then
        Error("Algebra A must be finite.");
    fi;
    d := Dimension(A);
    b := CanonicalBasis(A);
    bv := BasisVectors(b);
    repeat
        cf := List([1..d], x -> Random(F));
        x := LinearCombination(bv, cf);
        m := HOPF_MinPoly(A, x);
        facts := Factors( PolynomialRing(F), m);
    until Length(facts) > 1;
    one := MultiplicativeNeutralElement(A);
    f := facts[1];
    g := Product(facts{[2..Length(facts)]});
    return [Value(f, x, one), Value(g, x, one)];
end
);

#######################################################################
##
#F  HOPF_LeftIdentity( <A> )
##
##  This function takes as an argument a finite dimensional algebra 
##  <A>  and returns a left identity for  <A>, if such a thing exists.   
##  Otherwise it returns fail.
##
InstallGlobalFunction(HOPF_LeftIdentity,
  function(A)
    local F, b, bv, n, equ, zero, one, zerovec, vec, row, p, sol,
          i, j;

    F := LeftActingDomain(A);
    b := CanonicalBasis(A);
    bv := BasisVectors(b);
    n := Dimension(A);

    equ := [];
    zero := Zero(F);
    one := One(F);
    zerovec := ListWithIdenticalEntries(n^2, zero);
    vec := ShallowCopy(zerovec);

    for i in [1..n] do
        row := ShallowCopy(zerovec);
        for j in [1..n] do
            p := (j-1)*n;
            row{[p+1..p+n]} := Coefficients(b, b[i]*b[j]);
        od;
        Add(equ, row);
        vec[ (i-1)*n + i ] := one;
    od;
    sol := SolutionMat(equ,vec);
    if sol <> fail then
        sol := LinearCombination(bv, sol);
    fi;
    return sol;
end
);

#######################################################################
##
#F  HOPF_SingularIdempotent( <A> )
##
##  This function takes as an argument a finite dimensional simple 
##  algebra  <A>  and returns a primitive idempotent for  <A>?   
## 
##  This function is based on Ronyai, L., "Simple algebras are difficult", 
##  In 19th ACM Symposium on Theory of Computing (1987), pp. 398–408, and 
##  Ro ́nyai, L., "Computing the structure of finite algebras", Journal of 
##  Symbolic Computation 9 (1990), 355–373. See page 82 in Craig A. Struble,
##  "Analysis and implementation of algorithms for noncommutative algebra", 
##  PhD-thesis, Virginia Tech, 2000. 
##
InstallGlobalFunction(HOPF_SingularIdempotent,
  function(A)
    local z, b, bv, x, li, e;

    while not IsCommutative(A) do
        z := HOPF_ZeroDivisor(A);
        b := Basis(A);
        bv := BasisVectors(b);
        x := z[1];
        li := x*A;
        e := HOPF_LeftIdentity(li);

        A := Algebra(LeftActingDomain(A), List(bv, x->e*x*e));
    od;

    return MultiplicativeNeutralElement(A);
  end
);