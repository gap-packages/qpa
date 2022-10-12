# GAP Implementation
# This file was generated from 
# $Id: decomp.gi,v 1.6 2012/08/15 07:34:50 sunnyquiver Exp $

#######################################################################
##
#O  ComplementInFullRowSpace( <V> )
##
##  Computes the complement of a vector space in a row space. 
##
InstallMethod(ComplementInFullRowSpace, 
    "Compute the complement vector space of a row space", 
    true, 
    [IsRowSpace], 0,
    function(V)

    local M, gens, zero, n, i, needed, dims, t, K;
    
    K := LeftActingDomain(V);

    if Dimension(V) = 0 then
        Error("Cannot find complement of trivial space");
    fi;

    M := ShallowCopy(BasisVectors(Basis(V)));
    TriangulizeMat(M);
    dims := DimensionsMat(M);
    n := dims[2];
    needed := [1 .. n];

    zero := List([1 .. n], x -> Zero(K));

    if dims[1] = n then
        return TrivialSubspace(V);
    fi;

    gens := [];
    for i in [1 .. dims[1]] do
        t := Position(M[i], One(K));
        RemoveSet(needed, t);
    od;

    for i in [1 .. Length(needed)] do
        zero[needed[i]] := One(K);
        gens[i] := ShallowCopy(zero);
        zero[needed[i]] := Zero(K);
    od;

    return VectorSpace(K, gens);
end
);

#######################################################################
##
#O  LiftIdempotentsForDecomposition( <map>, <ids> )
##
##  Given a map  <map> :  A --> B (surjective?) and a set of 
##  (orthogonal) idempotents in  B, this function computes 
##  a (orthogonal) set of idempotents of preimages of the 
##  idempotents in  B? 
##
InstallMethod(LiftIdempotentsForDecomposition,
    "for an algebra mapping and list of idempotents",
    true,
    [IsAlgebraGeneralMapping,IsList],
    0,
    function(map,ids)

    local a, e, i, zero;

    ids := List(ids, x -> PreImagesRepresentative(map,x));
    a := Source(map);
    zero := Zero(a);
    e := Zero(LeftActingDomain(a))*ids[1];

    for i in [1..Length(ids)] do
        ids[i] := ids[i] - e*ids[i] - ids[i]*e + e*ids[i]*e;
        while not ids[i]^2 - ids[i] = zero do
            ids[i] := 3*ids[i]^2 - 2*ids[i]^3;
        od;
        e := e + ids[i];
    od;

    return AsSSortedList(ids);
end
);

#######################################################################
##
#O  IdempotentsForDecomposition( <a> )
##
##  Computes a complete of set primitive idempotents of the algebra  <a>. 
##
InstallMethod(IdempotentsForDecomposition,
    "for an algebra",
    true,
    [IsAlgebra],
    0,
    function(a)

    local map, semi, simp, cents, prims, i;

    map := NaturalHomomorphismByIdeal(a,RadicalOfAlgebra(a));
    semi := Range(map);
    cents := CentralIdempotentsOfAlgebra(semi);
    prims := [];

    for i in cents do
        simp := FLMLORByGenerators( LeftActingDomain(a), BasisVectors(Basis(semi*i)));
        SetParent(simp, a); 
        SetOne(simp, i);
        SetMultiplicativeNeutralElement(simp, i);   
        SetFilterObj(simp, IsAlgebraWithOne);
        prims := Concatenation(prims,PrimitiveIdempotents(simp));
    od;
   
    return LiftIdempotentsForDecomposition(map,prims);
end
);

#######################################################################
##
#O  DecomposeModule( <M> )
##
##  Given a module  <M>  this function computes a list of modules  L
##  such that  <M>  is isomorphic to the direct sum of the modules on 
##  the list  L. 
##
InstallMethod(DecomposeModule, 
    "for a path algebra matrix module", 
    true, 
    [IsPathAlgebraMatModule], 0, 
    function(M)

    local genmats, genmaps, basis, endo, idemmats, idemmaps, x;

    basis := CanonicalBasis(M);
    endo := EndOverAlgebra(M);
    genmaps := BasisVectors(Basis(endo));
    genmats := List(genmaps, x -> TransposedMat(x));
    idemmats := IdempotentsForDecomposition(AlgebraWithOne(LeftActingDomain(M), genmats));
    idemmaps := List(idemmats, x -> LeftModuleHomomorphismByMatrix(basis, TransposedMat(x), basis));
    for x in idemmaps do
        SetFilterObj(x, IsAlgebraModuleHomomorphism);
    od;
    idemmaps := List(idemmaps, x -> FromEndMToHomMM(M,x!.matrix));
    
    return List(idemmaps, x -> Image(x));
end
);

#######################################################################
##
#O  DecomposeModuleWithInclusions( <M> )
##
##  Given a module  <M>  this function computes a list of inclusions  L
##  such that  <M>  is isomorphic to the direct sum of the images of the
##  the inclusions in the list  L. 
##
InstallMethod(DecomposeModuleWithInclusions, 
    "for a path algebra matrix module", 
    true, 
    [IsPathAlgebraMatModule], 0, 
    function( M )

    local genmats, genmaps, basis, endo, idemmats, idemmaps, x;

    basis := CanonicalBasis( M );
    endo := EndOverAlgebra( M );
    genmaps := BasisVectors( Basis( endo ) );
    genmats := List( genmaps, x -> TransposedMat( x ) );
    idemmats := IdempotentsForDecomposition( AlgebraWithOne(LeftActingDomain( M ), genmats ) );
    idemmaps := List( idemmats, x -> LeftModuleHomomorphismByMatrix( basis, TransposedMat( x ), basis ) );
    for x in idemmaps do
        SetFilterObj( x, IsAlgebraModuleHomomorphism );
    od;
    idemmaps := List( idemmaps, x -> FromEndMToHomMM( M, x!.matrix ) );
    
    return List( idemmaps, x -> ImageInclusion( x ) );
end
);

#######################################################################
##
#O  DecomposeModule( <M> )
##
##  This function is an extension of above version of DecomposeModule
##  as this version takes adavantage of the fact that an argument might
##  have an direct sum decomposition already through being created by
##  the command  DirectSumOfQPAModules.
##
InstallOtherMethod( DecomposeModule,
    "for a path algebra",
    [ IsPathAlgebraMatModule and IsDirectSumOfModules ], 0,
    function( M ) 

    local directsummands, decomposition;

    if Dimension(M) = 0 then 
        return [];
    fi;
    directsummands := List(DirectSumInclusions(M), x -> Source(x));
    decomposition := List(directsummands, x -> DecomposeModule(x));
    
    return Flat(decomposition);
end
);



#InstallMethod(DecomposeModule, 
#  "for f. p. path algebra modules",
#  true, 
#  [IsFpPathAlgebraModule], 0,
#  function(V)
#    local endo, idemmats, idemmaps, dbasis, rbasis, x, genmats, genmaps;
#
#    endo:=End(ActingAlgebra(V), V);
#    genmaps:=BasisVectors(Basis(endo));
#    genmats:=List(genmaps, x->TransposedMat(x!.matrix));
#    dbasis:=genmaps[1]!.basissource;
#    rbasis:=genmaps[1]!.basisrange;
#    idemmats:=IdempotentsForDecomposition(AlgebraWithOne(LeftActingDomain(V),genmats));
#    idemmaps:=List(idemmats, x->LeftModuleHomomorphismByMatrix(dbasis, TransposedMat(x), rbasis));
#
#    for x in idemmaps do
#      SetFilterObj(x, IsAlgebraModuleHomomorphism);
#    od;
#
#    return List(idemmaps, x->Image(x,V));
#
#  end
#);
  
#######################################################################
##
#O  DecomposeModuleWithMultiplicities( <M> )
##
##  Given a PathAlgebraMatModule this function decomposes the module 
##  M  into indecomposable modules with multiplicities. First 
##  decomposing the module  M = M_1 + M_2 + ... + M_t, then checking
##  if  M_i  is isomorphic to M_1 for  i in [2..t], removing those
##  indices from [2..t] which are isomorphic to M_1, call this set 
##  rest, take the minimum from this, call it current, remove it from
##  rest, and continue as above until rest is empty. 
##
InstallMethod ( DecomposeModuleWithMultiplicities, 
    "for a IsPathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule ], 
    0,
    function( M )

    local L, current, non_iso_summands, rest, basic_summands, 
          multiplicities, temprest, i;

    if Dimension(M) = 0 then
        return fail; 
    fi;
#
#   Decompose the module  M.
#
    L := DecomposeModule(M);
#
#   Do the initial setup.
#
    current := 1;
    non_iso_summands := 1;
    rest := [2..Length(L)];
    basic_summands := [];
    multiplicities := [];
    Add(basic_summands,L[1]);
    Add(multiplicities,1);
# 
#   Go through the algorithm above to find the multiplicities.
#   
    while Length(rest) > 0 do
        temprest := ShallowCopy(rest);
        for i in temprest do
            if DimensionVector(L[current]) = DimensionVector(L[i]) then
                if CommonDirectSummand(L[current],L[i]) <> false then
                    multiplicities[non_iso_summands] := 
                        multiplicities[non_iso_summands] + 1;
                    RemoveSet(rest,i);
                fi;
            fi;
        od;
        if Length(rest) > 0 then 
            current := Minimum(rest);
            non_iso_summands := non_iso_summands + 1;
            Add(basic_summands,L[current]);
            RemoveSet(rest,current);
            Add(multiplicities,1);
        fi;
    od;

    return [basic_summands,multiplicities];
end
);

#######################################################################
##
#O  LiftIdempotent( <f>, <v> )
##
##  Given an onto homomorphism  <f>  of finite dimensional algebras and
##  an idempotent  <v>  in the range of  <f>, such that the kernel of  
##  <f>  is nilpotent, this function returns an idempotent  e  in the 
##  source of  <f>, such that  f(e) = v.
##
##  See for instance Andersen & Fuller: Rings and categories of modules, 
##  Proposition 27.1 for the algorithm used here.
##
InstallMethod( LiftIdempotent, 
    "for a morphism of algebras and one idempotent in the range",
    [ IsAlgebraGeneralMapping, IsRingElement ], 0,
    function( f, v )

    local g, x, y, nilindex, t, k;
    
    if v in Range(f) and v^2 = v then 
        g := PreImagesRepresentative(f, v);
        x := g - g*g;
        y := x;
        nilindex := 1;
        if y <> Zero(y) then 
            repeat 
                nilindex := nilindex + 1;
                y := x*y;
            until y = Zero(y);
        fi;
        if nilindex = 1 then
            return g;
        else
            t := Zero(y);
            for k in [1..nilindex] do
                t := t + (-1)^(k-1)*Binomial(nilindex,k)*g^(k-1);
            od;
	    return g^nilindex*t^nilindex;
        fi;
    else
        Error("entered map is not a IsAlgebraGeneralMapping or the entered element is not in the range of the map, ");
    fi;
end
);

#######################################################################
##
#O  LiftTwoOrtogonalIdempotents( <f>, <v>, <w> )
##
##  Given an onto homomorphism  <f>  of finite dimensional algebras, an 
##  idempotent  <v>  in the source of  <f>  and an idempotent  <w>  in 
##  the range of  <f>, such that the kernel of  <f>  is nilpotent and 
##  f(v)  and  w  are two orthogonal idempotents in the range of  <f>, 
##  this function returns an idempotent  e  in the source of  <f>, such
##  that  \{v, e\}  is a pair of orthogonal idempotents in the source of 
##  <f>. 
##
##  See for instance in Lam: A first course in noncommutative rings, 
##  Proposition 21.25 for the algorithm used here.
##
InstallMethod( LiftTwoOrthogonalIdempotents, 
    "for a morphism of algebras and one idempotent in the range",
    [ IsAlgebraGeneralMapping, IsRingElement, IsRingElement ], 0,
    function( f, v, w )

  local g, x, y, series, nilindex, temp;
    
    if  ImageElm( f, v ) * w <> Zero( w ) or w * ImageElm( f, v ) <> Zero( w ) then 
        Error("the entered idempotents are not orthogonal in the range of the algebra homomorphism,");
    fi;
    g := LiftIdempotent( f, w );
    x := g * v;
    y := x;
    series := One( g );
    nilindex := 1;       
    while not IsZero( y ) do 
      series := series + y;
      nilindex := nilindex + 1;
      y := x * y;
    od;
    temp := ( One( g ) - v) * series * g * ( One( g ) - x );
    return [ v, temp ];
  end
    );

#######################################################################
##
#O  BlockSplittingIdempotents( <M> )
##
##  Given a PathAlgebraMatModule  <M>  this function returns a set 
##  \{e_1,..., e_t\} of idempotents in the endomorphism of  <M>  such 
##  that  M \simeq  Im e_1\oplus \cdots \oplus Im e_t,
##  where each  Im e_i  is isomorphic to  X_i^{n_i}  for some 
##  decomposable module  X_i  and positive integer  n_i  for all i.
##  
InstallMethod( BlockSplittingIdempotents, 
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M )

    local   EndM,  map,  top,  orthogonalidempotents,  i,  etemp,  
            temp;
    
    if Dimension(M) = 0 then
        return [];
    fi;
    EndM := EndOverAlgebra(M);
    map := NaturalHomomorphismByIdeal( EndM, RadicalOfAlgebra( EndM ) ); 
    top := CentralIdempotentsOfAlgebra( Range( map ) );
    orthogonalidempotents := [ ];
    Add( orthogonalidempotents, LiftIdempotent( map, top[ 1 ] ) );
    for i in [ 1..Length( top ) - 1 ] do
        etemp := Sum(orthogonalidempotents{ [ 1..i ] } );
        temp := LiftTwoOrthogonalIdempotents( map, etemp, top[ i + 1 ] );
        Add( orthogonalidempotents, temp[ 2 ] );
    od;
    orthogonalidempotents := List( orthogonalidempotents, x -> FromEndMToHomMM( M, x ) );
    
    return orthogonalidempotents;
end
);

#######################################################################
##
#O  BlockDecompositionOfModule( <M> )
##
##  Given a PathAlgebraMatModule  <M>  this function returns a set of 
##  modules \{M_1,..., M_t\} such that  
##        M \simeq  M_1\oplus \cdots \oplus M_t,
##  where each  M_i  is isomorphic to  X_i^{n_i}  for some 
##  decomposable module  X_i  and positive integer  n_i  for all i.
##  
InstallMethod( BlockDecompositionOfModule, 
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M );
    
    return List(BlockSplittingIdempotents(M), e -> Image(e));
end
);

#######################################################################
##
#O  BasicVersionOfModule( <M> )
##
##  This function returns a basic version of the entered module  <M>, 
##  that is, if  <M> \simeq M_1^{n_1} \oplus \cdots \oplus M_t^{n_t}
##  where  M_i  is indecomposable, then  M_1\oplus \cdots \oplus M_t
##  is returned. At present, this function only work at best for 
##  finite dimensional (quotients of a) path algebra over a finite
##  field. If  <M>  is zero, then  <M> is returned.
##  
InstallMethod( BasicVersionOfModule, 
    "for a module over a (quotient of a) path algebra",
    [ IsPathAlgebraMatModule ], 0,
    function( M ) 

    local L;
    
    if Dimension(M) = 0 then
        return M;
    else
        L := DecomposeModuleWithMultiplicities(M);
        return DirectSumOfQPAModules(L[1]);
    fi;
end
  );

#######################################################################
##
#O  DecomposeModuleProbabilistic( <HomMM>, <M> )
##
##  Given a module  M  over a finite dimensional quotient of a path 
##  algebra over a finite field, this function tries to decompose the 
##  entered module  M  by choosing random elements in the endomorphism
##  ring of  M  which are non-nilpotent and non-invertible.  Such
##  elements splits the module in two direct summands, and the procedure
##  does this as long as it finds such elements.  The output is not 
##  guaranteed to be a list of indecomposable modules, but their direct
##  sum is isomorphic to the entered module  M. This was constructed as
##  joint effort by the participants at the workshop "Persistence, 
##  Representations, and Computation", February 26th - March 2nd, 2018".
##  This is an experimental function, so use with caution.
##  
InstallMethod( DecomposeModuleProbabilistic, 
   "for a path algebra matmodule",
   [ IsHomogeneousList, IsPathAlgebraMatModule ], 0,
   function( HomMM, M ) 
    
  local   k,  d,  l,  n,  V,  i,  p,  f,  g,  j,  h,  pi,  nu,  
          nuprime,  piprime,  phi,  HomImhImh,  gensImhImh,  t,  s,  
          VImhImh,  basisImhImh,  HomKerhKerh,  gensKerhKerh,  
          VKerhKerh,  basisKerhKerh;
    
    k := LeftActingDomain( M );
    if not IsFinite( k ) then 
        Error( "The entered module is not over a finite field.\n" );
    fi;
    d := Maximum( DimensionVector( M ) );
    l := Int( Ceil( Log2( 1.0*( d ) ) ) );
    #    n := Length( B );
    n := Length( HomMM );
    V := FullRowSpace( k, n );
    i := 1; 
    repeat
        i := i + 1; 
        p := Random( V ); 
        #        f := LinearCombination( B, p ); 
        f := LinearCombination( HomMM, p ); 
        if DeterminantMat( FromHomMMToEndM( f ) ) = Zero( k ) then
            g := f;
            for j in [ 1..l + 1 ] do
                g := g * g;
            od;
            if g <> Zero( g ) then
                h := g;
                pi := ImageProjection( h );
                nu := ImageInclusion( h );
                nuprime := KernelInclusion( h );
                piprime := CoKernelProjection( h );
                phi := IsomorphismOfModules( Range( piprime ), Source( nuprime ) );
                piprime := piprime * phi;
                #
                #
                HomImhImh := Unique( List( HomMM, x -> nu * x * pi ) );
                gensImhImh := List( HomImhImh, h -> Flat( FromHomMMToEndM( h ) ) );
                t := Length( gensImhImh[ 1 ] );
                s := Dimension( Range( pi ) );
                VImhImh := Subspace( k^t, gensImhImh );
                basisImhImh := List( BasisVectors( Basis( VImhImh ) ), b -> List( [ 1..s ], x -> b{ [ 1 + ( x - 1 ) * s..s * x ] } ) );
                HomImhImh := List( basisImhImh, b -> FromEndMToHomMM( Range( pi ), b ) );
                #
                #
                HomKerhKerh := Unique( List( HomMM, x -> nuprime * x * piprime ) );
                gensKerhKerh := List( HomKerhKerh, h -> Flat( FromHomMMToEndM( h ) ) );
                t := Length( gensKerhKerh[ 1 ] );
                s := Dimension( Range( piprime ) );
                VKerhKerh := Subspace( k^t, gensKerhKerh );
                basisKerhKerh := List( BasisVectors( Basis( VKerhKerh ) ), b -> List( [ 1..s ], x -> b{ [ 1 + ( x - 1 ) * s..s * x ] } ) );
                HomKerhKerh := List( basisKerhKerh, b -> FromEndMToHomMM( Range( piprime ), b ) );
                return Flat( [ DecomposeModuleProbabilistic( HomImhImh, Range( pi ) ), 
                               DecomposeModuleProbabilistic( HomKerhKerh, Range( piprime ) ) ] );
            fi;
        fi;
    until i = 2 * Size( k ) * n;
      
    return M; 
end
  );

InstallMethod ( DecomposeModuleViaTop, 
"for a PathAlgebraMatModule",
[ IsPathAlgebraMatModule ],
function( M )
    
  local K, HomMM, HomTopMM, EndTopMM, endTopMM, idempotents, EndM, V, 
        W, EndTopM, g;
    
  K := LeftActingDomain( M );
  if not IsFinite( K ) then
    Error( "The entered module is not a module over a finite field.\n" );
  fi;
  HomMM := HomOverAlgebra( M, M );
  HomTopMM := List( HomMM, TopOfModule );
  EndTopMM := List( HomTopMM, FromHomMMToEndM );
  endTopMM := Algebra( K, EndTopMM );
  idempotents := IdempotentsForDecomposition( endTopMM );
  EndM := List( HomMM, FromHomMMToEndM );
  V := Algebra( K, EndM, "basis" );
  SetOne( V, MultiplicativeNeutralElement( V ) );
  W := Algebra( K, EndTopMM );
  SetOne( W, MultiplicativeNeutralElement( W ) );
  EndTopM := List( HomMM, h -> FromHomMMToEndM( TopOfModule( h ) ) );
  g := AlgebraHomomorphismByImages( V, W, EndM, EndTopM );
  idempotents := LiftingCompleteSetOfOrthogonalIdempotents( g, idempotents );
  idempotents := List( idempotents, e -> FromEndMToHomMM( M, e ) );
  
  return List( idempotents, e -> Image( e ) );
end
  );

#######################################################################
##
#O  DecomposeModuleViaCharPoly( <M> )
##
##  Given a module  M  over a finite dimensional quotient of a path 
##  algebra over a finite field, this function tries to decompose the 
##  entered module  M  by choosing random elements in the endomorphism
##  ring of  M computing a factorization of their characteristic 
##  polynomials. 
##  
InstallMethod( DecomposeModuleViaCharPoly, 
   "for a path algebra matmodule",
   [ IsPathAlgebraMatModule ],
   function( M ) 
    
  local K, homMM, dimM, t, V, num_vert, num_repeats, p, h, matrices, 
        charpolys, nonzero_vert, m, factorlist, numfactors, 
        occurringprimes, multiplicities, maxmultiplicities, polys, 
        mats, f, mat, i, homs;
    
    K := LeftActingDomain( M );
    if not IsFinite( K ) then 
        Error( "The entered module is not over a finite field.\n" );
    fi;
    homMM := HomOverAlgebra( M, M );
    dimM := DimensionVector( M );
    t := Length( homMM );
    V := FullRowSpace( K, t );
    num_vert := Length( DimensionVector( M ) );
    num_repeats := 0;
    repeat
      num_repeats := num_repeats + 1;
      p := Random( V );
      h := LinearCombination( homMM, p );
      matrices := MatricesOfPathAlgebraMatModuleHomomorphism( h );
      charpolys := [ ];
      nonzero_vert := 0;
      for m in matrices do
        if not IsZero( m ) then 
          nonzero_vert := nonzero_vert + 1;
          Add(charpolys, CharacteristicPolynomial( m ) );
        fi;
      od;
      factorlist := List( charpolys, Factors );
      numfactors := List( factorlist, Length );
      occurringprimes := Unique( Flat( factorlist ) );
      if Length( occurringprimes ) > 1 then 
        multiplicities := List( [ 1..Length( occurringprimes ) ], 
                                j -> List( [1..nonzero_vert], i -> Length( Positions( factorlist[i], occurringprimes[ j ] ) ) ) );
        maxmultiplicities := List( multiplicities, m -> Maximum( m ) );
        polys := List( [ 1..Length( occurringprimes ) ], i -> occurringprimes[ i ]^maxmultiplicities[ i ] );
        
        mats := [];
        for f in polys do
          mat := [ ];
          for i in [ 1..num_vert ] do
            if dimM[ i ] = 0 then
              Add( mat, NullMat( 1, 1, K ) );
            else
              Add( mat, Value( f, matrices[ i ] ) );
            fi;
          od;
          Add( mats, mat );
        od;
        homs := [ ];
        for m in mats do
          Add( homs, RightModuleHomOverAlgebra( M, M, m ) );
        od;

        return Flat( List( List( homs, Kernel ), m -> DecomposeModuleViaCharPoly( m ) ) );
      fi;
    until 
      num_repeats = 2 * Size( K ) * Length( homMM );

    return [ M ];
end
  );