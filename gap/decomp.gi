# GAP Implementation
# This file was generated from 
# $Id: decomp.gi,v 1.2 2012/02/27 12:26:34 sunnyquiver Exp $
InstallMethod(ComplementInFullRowSpace, 
  "Compute the complement vector space of a row space", 
  true, 
  [IsRowSpace], 0,
  function(V)

    local M, gens, zero, n, i, needed, dims, t, K;
    K:=LeftActingDomain(V);

    if Dimension(V)=0 then
      Error("Cannot find complement of trivial space");
    fi;

    M:=ShallowCopy(BasisVectors(Basis(V)));
    TriangulizeMat(M);
    dims:=DimensionsMat(M);
    n:=dims[2];
    needed:=[1 .. n];

    zero:=List([1 .. n], x->Zero(K));

    if dims[1]=n then
      return TrivialSubspace(V);
    fi;

    gens:=[];
    for i in [1 .. dims[1]] do
      t:=Position(M[i], One(K));
      RemoveSet(needed, t);
    od;

    for i in [1 .. Length(needed)] do
      zero[needed[i]]:=One(K);
      gens[i]:=ShallowCopy(zero);
      zero[needed[i]]:=Zero(K);
    od;

    return VectorSpace(K, gens);

  end
);


InstallMethod(LiftIdempotentsForDecomposition,
  "for an algebra mapping and list of idempotents",
  true,
  [IsAlgebraGeneralMapping,IsList],
  0,
  function(map,ids)

    local a,e,i,zero;

    ids:=List(ids,x->PreImagesRepresentative(map,x));
    a:=Source(map);
    zero:=Zero(a);
    e:=Zero(LeftActingDomain(a))*ids[1];

    for i in [1..Length(ids)] do
      ids[i]:=ids[i]-e*ids[i]-ids[i]*e+e*ids[i]*e;

      while not ids[i]^2-ids[i]=zero do
       ids[i]:=3*ids[i]^2-2*ids[i]^3;
      od;

      e:=e+ids[i];

    od;

    return AsSSortedList(ids);

  end
);


InstallMethod(IdempotentsForDecomposition,
  "for an algebra",
  true,
  [IsAlgebra],
  0,
  function(a)

    local map,semi, simp, cents, prims,i;

    map:=NaturalHomomorphismByIdeal(a,RadicalOfAlgebra(a));
    semi:=Range(map);
    cents:=CentralIdempotentsOfAlgebra(semi);
    prims:=[];

    for i in cents do
      simp:=semi*i;
      prims:=Concatenation(prims,PrimitiveIdempotents(simp));
    od;
   
    return LiftIdempotentsForDecomposition(map,prims);

  end
);


InstallMethod(DecomposeModule, 
  "for a path algebra matrix module", 
  true, 
  [IsPathAlgebraModule], 0, 
  function(M)

    local genmats, genmaps, basis, endo, idemmats, idemmaps, x;

    basis:=CanonicalBasis(M);
    endo:=End(ActingAlgebra(M), M);
    genmaps:=BasisVectors(Basis(endo));
    genmats:=List(genmaps, x->TransposedMat(x!.matrix));
    idemmats:=IdempotentsForDecomposition(AlgebraWithOne(LeftActingDomain(M), genmats));
    idemmaps:=List(idemmats, x->LeftModuleHomomorphismByMatrix(basis, TransposedMat(x), basis));

    for x in idemmaps do
      SetFilterObj(x, IsAlgebraModuleHomomorphism);
    od;

    return List(idemmaps, x->Image(x,M));

  end
);


InstallMethod(DecomposeModule, 
  "for f. p. path algebra modules",
  true, 
  [IsFpPathAlgebraModule], 0,
  function(V)
    local endo, idemmats, idemmaps, dbasis, rbasis, x, genmats, genmaps;

    endo:=End(ActingAlgebra(V), V);
    genmaps:=BasisVectors(Basis(endo));
    genmats:=List(genmaps, x->TransposedMat(x!.matrix));
    dbasis:=genmaps[1]!.basissource;
    rbasis:=genmaps[1]!.basisrange;
    idemmats:=IdempotentsForDecomposition(AlgebraWithOne(LeftActingDomain(V),genmats));
    idemmaps:=List(idemmats, x->LeftModuleHomomorphismByMatrix(dbasis, TransposedMat(x), rbasis));

    for x in idemmaps do
      SetFilterObj(x, IsAlgebraModuleHomomorphism);
    od;

    return List(idemmaps, x->Image(x,V));

  end
);
