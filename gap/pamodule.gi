# GAP Implementation
# This file was generated from 
# $Id: pamodule.gi,v 1.4 2010/09/24 08:59:12 sunnyquiver Exp $

ZeroModElement:=function(fam,zero)
  local result,i;
  result:=[];
  for i in [1..Length(fam!.vertices)] do
    result[i]:=List([1..fam!.vertices[i]],x->zero);
    result[i][1]:=zero;
  od;

  return result;
end;


InstallGlobalFunction(PathModuleElem,
  function( arg )
    local i, fam, basis, vertices, elem;
   
    if Length( arg ) <> 2 then
      Error("Function PathModuleElem must take 2 arguments");
    elif not IsFamily( arg[1] ) and not IsCollection( arg[2] ) then
      Error("Arguments are of wrong type");
    else
      fam:=arg[1];
      basis:=arg[2];
      vertices:=fam!.vertices;

      if Length( vertices ) <> Length( basis ) then
         Error("Incorrect amount of vertices");
      fi;

      if not ForAll([1..Length(vertices)],i->vertices[i]=0 or vertices[i]=Length(basis[i])) then
         Error("Dimensions of element are incorrect");
      fi;

      elem:=Objectify( NewType( fam, IsPositionalObjectRep), []);
      elem![1]:=[];
      for i in [1 .. Length(vertices)] do
        elem![1][i]:=ShallowCopy(basis[i]);
      od;

      return elem;
    fi;

  end
);


InstallMethod(\+, 
  "for two module elements",
  IsIdenticalObj,
  [IsPathModuleElem, IsPathModuleElem], 0,
  function(elem1, elem2)

    local fam;
    fam:=FamilyObj(elem1);   

    return PathModuleElem( fam, elem1![1]+elem2![1] );
  end
);


InstallMethod(\-, 
  "for two module elements",
  IsIdenticalObj,
  [IsPathModuleElem, IsPathModuleElem], 0,
  function(elem1, elem2)

    local fam;
    fam:=FamilyObj(elem1);   

    return PathModuleElem( fam, elem1![1]-elem2![1] );
  end
);


InstallMethod(\*,
  "for a module element and a scalar",
  true,
  [IsScalar, IsPathModuleElem], 0,
  function(s,m)
    local fam;
    fam := FamilyObj(m);
    return PathModuleElem( fam, List(m![1],x->s*x));
  end
);


InstallOtherMethod(\*,
  "for a module element and a scalar",
  true,
  [IsPathModuleElem, IsScalar], 0,
  function(m,s)
    local fam;
    fam := FamilyObj(m);
    return PathModuleElem( fam, List(m![1],x->x*s));
  end
);


InstallOtherMethod( \^,
 "for a module element and a path",
 true,
 [IsPathModuleElem, IsRingElement ] , 0 ,
   function(elem, path)
     local i, x, result, fam, m, temp, matrix, source, target, K, zero,
         modterm, walk,  terms;
     fam:=FamilyObj(elem);

     if not IsIdenticalObj(ElementsFamily(FamilyObj(fam!.pathAlgebra)),
                           FamilyObj(path)) then
       Error("Arguments are not compatible");
     else
       terms:=CoefficientsAndMagmaElements(path);
       zero:=Zero(LeftActingDomain(fam!.pathAlgebra));
       result:=ZeroModElement(fam, zero);

       for i in [1,3 .. Length( terms ) -1] do
          modterm:=ZeroModElement(fam, zero);
          walk:=WalkOfPath(terms[i]);
          matrix:= One(zero);

          for x in walk do
             matrix:= matrix * fam!.matrices[x!.gen_pos];
          od;

          source:= SourceOfPath(terms[i]);
          target:= TargetOfPath(terms[i]);
          matrix:= terms[i+1] * matrix;
          modterm[target!.gen_pos]:= elem![1][source!.gen_pos] * matrix;
          result:= result + modterm;
       od;

       return PathModuleElem(fam, result);  

     fi; 
  end
);


InstallMethod(ZeroOp,
  "for elements of modules",
  true,
  [ IsPathModuleElem], 0,
  function( obj)
    return PathModuleElem(FamilyObj(obj), ZeroModElement(FamilyObj(obj),Zero(obj![1][1][1])));
  end
);


InstallMethod(AdditiveInverseOp,
  "for elements of modules",
  true,
  [IsPathModuleElem], 0,
  function(obj)
     return PathModuleElem(FamilyObj(obj), List(obj![1], x->-1*x));
  end
);


InstallMethod(\=,
  "for elements of modules",
  true,
  [ IsPathModuleElem, IsPathModuleElem], 0,
  function(m,n)
    return FamilyObj(m)=FamilyObj(n) and m![1]=n![1];
  end
);


InstallMethod( PrintObj,
  "for elements of modules",
  true,
  [ IsPathModuleElem ], 0, 
  function( obj ) 
     Print(obj![1]);
  end
);


InstallMethod( NewBasis,
  "for a space of path module elems",
  IsIdenticalObj,
  [ IsFreeLeftModule and IsPathModuleElemCollection,
    IsPathModuleElemCollection and IsList ], NICE_FLAGS+1,
  function( V, vectors )
    local B;  # The basis

    B := Objectify( NewType( FamilyObj( V ),
                             IsBasis and
                             IsBasisOfPathModuleElemVectorSpace and
                             IsAttributeStoringRep ),
                    rec() );
    SetUnderlyingLeftModule( B, V );
    SetBasisVectors( B, vectors );
    return B;
  end
);


InstallMethod( BasisNC,
  "for a space spanned by path module elems",
  IsIdenticalObj,
  [ IsFreeLeftModule and IsPathModuleElemCollection,
      IsPathModuleElemCollection and IsList ],
  NICE_FLAGS+1,
  function( V, vecs )
    local uvecs, uVecSpace, B;

    vecs:=Filtered(vecs, x->not IsZero(x));
    B:=NewBasis(V,vecs);
    uvecs:=List(vecs, x->Flat(x![1]));
    uVecSpace:=LeftModuleByGenerators(LeftActingDomain(V), uvecs);
    B!.delegateBasis:=BasisNC(uVecSpace,uvecs);

    return B;
  end
);


InstallOtherMethod( Basis,
  "for a space spanned by path module elems",
  IsIdenticalObj,
  [ IsFreeLeftModule and IsPathModuleElemCollection,
    IsPathModuleElemCollection and IsList ], NICE_FLAGS+1,
  function( V, vecs )
    local uvecs, uVecSpace, B;

    if not ForAll(vecs, x->x in V) then
       return fail;
    fi;

    vecs:=Filtered(vecs, x->not IsZero(x));
    B:=NewBasis(V,vecs);
    uvecs:=List(vecs, x->Flat(x![1]));
    uVecSpace:=LeftModuleByGenerators(LeftActingDomain(V), uvecs);
    B!.delegateBasis:=Basis(uVecSpace,uvecs);

    return B;
  end
);


InstallMethod(BasisOfDomain, 
  "for path module elem spaces",
  true,
  [IsFreeLeftModule and IsPathModuleElemCollection],
  NICE_FLAGS+1,
  function(V)
    return BasisNC(V, ShallowCopy(GeneratorsOfLeftModule(V)));
  end
);


InstallMethod(Basis,
  "for path module elem spaces",
  true,
  [IsFreeLeftModule and IsPathModuleElemCollection],
  NICE_FLAGS+1,
  function(V)
    return BasisNC(V, ShallowCopy(GeneratorsOfLeftModule(V)));
  end
);


InstallMethod( Coefficients, 
  "for a basis of a path module elem space and path module elem in the space",
  IsCollsElms,
  [IsBasisOfPathModuleElemVectorSpace, IsPathModuleElem],
  0,
  function(B,v)
    return Coefficients(B!.delegateBasis, Flat(v![1]));
  end
);


InstallMethod(CanonicalBasis,
  "for a path algebra matrix module",
  true,
  [IsFreeLeftModule and IsPathAlgebraMatModule], NICE_FLAGS+1,
  function(V)
    local basis, vgens, gens, x, i, fam, zero, MB, Vfam;

    gens:=BasisVectors(Basis(V));
    Vfam:=FamilyObj(gens[1]);
    gens:=List(gens, x->ExtRepOfObj(x));
    fam:=FamilyObj(gens[1]);
    vgens:=List(fam!.vertices, x->[]);

    for x in gens do
      for i in [1 .. Length(fam!.vertices)] do 
        if not Zero(x![1][i])=x![1][i] then
          Add(vgens[i], x![1][i]);
        fi;
      od;
    od;

    zero:=ZeroModElement(fam, Zero(gens[1]![1][1][1]));
    basis:=[];

    for x in [1 .. Length(fam!.vertices)] do
      for i in vgens[x] do
        zero[x]:=i;
        Add(basis, PathModuleElem(fam, ShallowCopy(zero)));
        zero[x]:=Zero(i);
      od;
    od;

    basis:=List(basis, x->ObjByExtRep(Vfam, x));
    MB:=BasisNC(V, basis);
    SetIsCanonicalBasis(MB,true);

    return MB;
  end
);


CreateModuleBasis:=function(fam)
  local vertices, track, base, K,i,x, basis;

  vertices:=fam!.vertices;
  K:=LeftActingDomain(fam!.pathAlgebra);
  base:=ZeroModElement(fam, Zero(K));
  track:=1;
  basis := [];

  for i in [1 .. Length(vertices)] do
    for x in [1 .. vertices[i]] do  
      base[i][x]:= One(K);
      basis[track]:= PathModuleElem(fam, base);
      track:=track + 1;
      base[i][x]:=Zero(K);
    od;
  od;

  return basis;
end;


InstallMethod(RightModuleOverPathAlgebra,
  "for a path algebra and list of matrices",
  true,
  [IsPathAlgebra, IsCollection], 0,
  function( R, gens )
    local a, dim, source, target, basis, i, x, Fam, 
          vertices, matrices, quiver, M, vlist, alist, K, dim_M;
    matrices:=[];
    quiver:=QuiverOfPathRing(R);
    vlist:=VerticesOfQuiver(quiver);
    K:=LeftActingDomain(R);
    alist:=ArrowsOfQuiver(quiver);
    vertices:=[];
#    
#  First checking if all arrows has been assigned some value.
#
    if Length(gens) < Length(alist) then 
       Error("Each arrow has not been assigned a matrix.");
    fi;
#
#  Setting the multiplication by the vertices.
#          
    for i in [1 .. Length(vlist)] do
      matrices[i]:= One(K);
    od;
#
#  Setting, partially, the multiplication by the arrows, taking into account
#  the possible different formats of the input. 
#
#  Input of the form ["a",[[..],...,[..]]], where "a" is the label of 
#  some arrow in the quiver
#
    if IsString(gens[1][1]) then                 
       for i in [1 .. Length ( gens )] do
          a:=gens[i][1];
          matrices[quiver.(a)!.gen_pos]:=gens[i][2];
       od;
#
#  Input of the form [[matrix_1],[matrix_2],...,[matrix_n]]
#
    elif IsMatrix(gens[1]) then
       for i in [1 .. Length ( gens )] do
          matrices[i + Length(vlist)]:=gens[i];
       od;
    else
#
#  Input of the form [[alist[1],[matrix_1]],...,[alist[n],[matrix_n]]] 
#  where alist is a list of the vertices in the quiver.
#  
       for i in [1 .. Length ( gens )] do
          a:=gens[i][1];
          matrices[a!.gen_pos]:=gens[i][2];
       od;
    fi;

    for i in [1 .. Length(vlist)] do
      vertices[i]:=-1;
    od; 
#
#  Setting dimensions and checking the usage of the zero space format.
#
    dim:=[];
    for x in alist do
      if IsMatrix(matrices[x!.gen_pos]) then
        dim:= DimensionsMat ( matrices[x!.gen_pos] );
      else
        dim:=matrices[x!.gen_pos];
        if ( dim[1] > 0 and dim[2] = 0 ) then
          matrices[x!.gen_pos]   := NullMat(dim[1],1,K);
        elif ( dim[1] = 0 and dim[2] > 0 ) then
            matrices[x!.gen_pos] := NullMat(1,dim[2],K);
	elif ( dim[1] = 0 and dim[2] = 0 ) then 
            matrices[x!.gen_pos] := NullMat(1,1,K);
        else
            Error("A vertex cannot have negative dimension or wrong usage of the zero space format.");
        fi;
      fi;

      source:=SourceOfPath(x);
      target:=TargetOfPath(x);
#
#  Checking if all the matrices entered are compatible with respect 
#  to the dimension of the vectorspaces at each vertex.
#
      if vertices[source!.gen_pos] = -1 then
        vertices[source!.gen_pos]:= dim[1];
      elif vertices[source!.gen_pos] <> dim[1] then
        Error("Dimensions of matrices do not match");
      fi;

      if vertices[target!.gen_pos] = -1 then
        vertices[target!.gen_pos]:= dim[2];
      elif vertices[target!.gen_pos] <> dim[2] then
        Error("Dimensions of matrices do not match");
      fi;
    od;

    dim_M := 0;
    for i in [1..Length(vlist)] do
    	dim_M := dim_M + vertices[i];
    od;
    Fam := NewFamily( "PathAlgModuleElementsFamily", IsPathModuleElem );
    SetFilterObj( Fam, IsPathModuleElemFamily );
    Fam!.vertices:=vertices;
    Fam!.matrices:=matrices;
    Fam!.pathAlgebra:=R;
    if dim_M > 0 then 
        basis := CreateModuleBasis(Fam);
    else 
    	basis := ZeroModElement(Fam, Zero(K));
    fi;
    M:= RightAlgebraModuleByGenerators(R, \^, basis);
    SetIsPathAlgebraMatModule(M,true);
    SetIsWholeFamily(M,true);

    return M;
  end
);


InstallMethod( ViewObj, 
  "for modules over path algebras",
  true,
  [ IsRightModule and IsPathModuleElemCollection ], 0,
  function (M)
    local fam, i;

    Print("<Module over path algebra ");
    fam:= ElementsFamily(FamilyObj(M));
    View(fam!.pathAlgebra);
    Print(" with ", Length(GeneratorsOfRightModule(M)), " generators ");
    Print(">");

  end
);


InstallMethod( \in, 
  "for a path module elem and a path algebra matrix module",
  IsElmsColls,
  [IsPathModuleElem, 
   IsFreeLeftModule and IsPathModuleElemCollection], NICE_FLAGS+1,
  function( v, V )
    local B, cf;

    B := BasisOfDomain(V);
    cf := Coefficients(B, v);

    return cf <> fail;
  end
);


InstallMethod( SubAlgebraModule,
  "for path algebra module and a list of submodule generators",
  IsIdenticalObj,
  [ IsFreeLeftModule and IsPathAlgebraMatModule,
    IsAlgebraModuleElementCollection and IsList], 0,
  function( V, gens )
    local sub, ngens, vlist, fam, alggens,x, i;

    fam:=FamilyObj(ExtRepOfObj(gens[1]));
    alggens:=GeneratorsOfAlgebra(fam!.pathAlgebra);
    ngens:=[];
    vlist:=List([1 .. Length(VerticesOfQuiver(QuiverOfPathAlgebra(fam!.pathAlgebra)))], x->alggens[x]);

    for x in gens do
       for i in vlist do 
          Add(ngens, x^i);
       od;
    od;

    ngens:=Filtered(ngens, x-> not IsZero(x));

    sub := Objectify( NewType( FamilyObj( V ),
                      IsLeftModule and IsAttributeStoringRep),
                      rec() );
    SetParent( sub, V );
    SetIsAlgebraModule( sub, true );
    SetIsPathAlgebraMatModule( sub, true );
    SetLeftActingDomain( sub, LeftActingDomain(V) );
    SetGeneratorsOfAlgebraModule( sub, ngens );

    if HasIsFiniteDimensional( V ) and IsFiniteDimensional( V ) then
      SetIsFiniteDimensional( sub, true );
    fi;

    if IsLeftAlgebraModuleElementCollection( V ) then
      SetLeftActingAlgebra( sub, LeftActingAlgebra( V ) );
    fi;

    if IsRightAlgebraModuleElementCollection( V ) then
      SetRightActingAlgebra( sub, RightActingAlgebra( V ) );
    fi;

    return sub;
  end
);


InstallMethod(SubmoduleAsModule, 
  "for submodules of path algebra matrix modules",
  true, 
  [IsAlgebraModule and IsPathAlgebraMatModule], 0,
  function( N )

    local gens, fam, R, F, alist, verts, sdims, vgens, x, i, 
          spaces, vertbases, mappings, oldmats, newmats, a, 
          apos, source, target, spos, tpos;

    gens:=List(BasisVectors(CanonicalBasis(N)), x->ExtRepOfObj(x));
    fam:=FamilyObj(gens[1]);
    R:=fam!.pathAlgebra;
    F:=LeftActingDomain(R);
    alist:=ArrowsOfQuiver(QuiverOfPathAlgebra(R));
    verts:=fam!.vertices;
    sdims:=List(verts, x->0);
    vgens:=List(verts, x->[]);

    for x in gens do
      for i in [1 .. Length(verts)] do 
        if not Zero(x![1][i])=x![1][i] then
          Add(vgens[i], x![1][i]);
          sdims[i]:=sdims[i]+1;
        fi;
      od;
    od;

    spaces:=[];
    vertbases:=[];

    for i in [1 .. Length(sdims)] do
      if sdims[i]<>0 then
        spaces[i]:=VectorSpace(F,vgens[i]);
        vertbases[i]:=CanonicalBasis(spaces[i]);
      fi;
    od;

    mappings:=List(vertbases,x->BasisVectors(x));
    oldmats:=ShallowCopy(fam!.matrices);
    newmats:=ShallowCopy(fam!.matrices);

    for a in alist do
      apos:=a!.gen_pos;
      source:=SourceOfPath(a);
      target:=TargetOfPath(a);
      spos:=source!.gen_pos;
      tpos:=target!.gen_pos;

      if sdims[spos]=0 or sdims[tpos]=0 then
        newmats[apos]:=[sdims[spos], sdims[tpos]];
      else    
        newmats[apos]:=List([1..Length(mappings[spos])],x->Coefficients(vertbases[tpos],mappings[spos][x]*oldmats[apos]));
      fi;
    od;

    return RightModuleOverPathAlgebra(R, List(alist, x->newmats[x!.gen_pos]));
  end
);


InstallOtherMethod(NaturalHomomorphismBySubAlgebraModule,
  "for modules over path algebras and their submodules",
  true,
  [IsPathAlgebraMatModule, IsPathAlgebraMatModule], 0,
  function(M,N)

    local pfam, Q, pgens, qgens, sgens, a, alist, R, pmats, qmats, pdims, sdims, qdims, hom,
          dualspaces, dualbasis, vgens, i, x, subspaces, F, Fam, vwspaces, vwbasis, vec, genmats,
          vwgens, apos, tpos, spos, image, basis, currmat, track, base, temp, images, qbasis;

    pgens:=List(BasisVectors(CanonicalBasis(M)), x->ExtRepOfObj(x));
    sgens:=List(BasisVectors(CanonicalBasis(N)), x->ExtRepOfObj(x));
    pfam:=FamilyObj(pgens[1]);

    if not IsSubspace(M,N) then
      Error("The second argument must be a submodule of the first argument");
    fi;

    pdims:=ShallowCopy(pfam!.vertices);
    sdims:=List(pdims, x->0);
    R:=pfam!.pathAlgebra;
    pmats:=ShallowCopy(pfam!.matrices);
    F:=LeftActingDomain(R);
    alist:=ArrowsOfQuiver(QuiverOfPathRing(R));
    vgens:=List(pdims, x->[]);

    for x in sgens do
      for i in [1 .. Length(pdims)] do 
        if not IsZero(x![1][i]) then
          Add(vgens[i], x![1][i]);
          sdims[i]:=sdims[i]+1;
        fi;
      od;
    od;

    qdims:=List([1..Length(pdims)], x->pdims[x]-sdims[x]);
    subspaces:=List(vgens, x->VectorSpace(F, x));
    vgens:=List(subspaces, x->BasisVectors(Basis(x)));
    dualspaces:=[];

    for i in [1 .. Length(pdims)] do
      if sdims[i]=0 then
        dualspaces[i]:=F^pdims[i];
      else 
        dualspaces[i]:=ComplementInFullRowSpace(subspaces[i]);
      fi;
    od;

    dualbasis:=List(dualspaces, x->BasisVectors(Basis(x)));
    vwgens:=List([1..Length(pdims)], x->Concatenation(vgens[x], dualbasis[x]));
    vwspaces:=List(vwgens, x->VectorSpace(F,x));
    vwbasis:=[];

    for i in [1 .. Length(pdims)] do
      if sdims[i]=0 then
        vwbasis[i]:=Basis(vwspaces[i]);
      else 
        vwbasis[i]:=Basis(vwspaces[i], vwgens[i]);
      fi;
    od;

    qmats:=ShallowCopy(pmats);

    for a in alist do
      currmat:=[];
      apos:=a!.gen_pos;
      spos:=SourceOfPath(a)!.gen_pos;
      tpos:=TargetOfPath(a)!.gen_pos;

      if qdims[spos]=0 then
        currmat:=[qdims[spos],qdims[tpos]];
      elif qdims[tpos]=0 then
        currmat:=[qdims[spos],qdims[tpos]];
      else    
        for i in [1..qdims[spos]] do
          image:=Coefficients(vwbasis[tpos], dualbasis[spos][i]*pmats[apos]);
          currmat[i]:=List([1..qdims[tpos]], x->image[x+sdims[tpos]]);
        od;
      fi;

      qmats[apos]:=currmat;
    od;

    genmats:=List(alist, x->ShallowCopy(qmats[x!.gen_pos]));
    Q:=RightModuleOverPathAlgebra(R, genmats);
    track:=1;
    images:=[];
    qbasis:=BasisVectors(Basis(Q));
    qgens:=List(qbasis, x->ExtRepOfObj(x));
    Fam:=FamilyObj(qgens[1]);
    base:=ZeroModElement(Fam, Zero(F));

    for i in [1 .. Length(qdims)] do
      for x in [1 .. pdims[i]] do
        if qdims[i]=0 then
          images[track]:=PathModuleElem(Fam, base);
          track:=track+1;
        else
          temp:=base;
          vec:=Coefficients(vwbasis[i],pgens[track]![1][i]);
          base[i]:=List([1 .. qdims[i]],x->vec[x+sdims[i]]);
          images[track]:=PathModuleElem(Fam, base);
          base:=temp;
          track:=track+1;
        fi;
      od;
    od; 

    images:=List(images, x->ObjByExtRep(FamilyObj(qbasis[1]), x));
    hom := LeftModuleHomomorphismByImagesNC( M, Q, BasisVectors(Basis(M)), images );
    SetIsSurjective(hom, true);

    return hom;

  end
);


HOPF_MatrixModuleHom:=function(R, M, N)
  local F, zero, Mfam, Nfam, quiver, Mdims, Ndims, Mmaps, Nmaps,
        l, arrows, i, arrowPos, sourcePos, targetPos,
        a, tMat, numColumns, numRows, r, c, row, col, startCol,
        startRow, equations, j, k, MBasis, NBasis, mat, x, y,
        Msum, Nsum, ns, basis, maps, b;

  if ActingAlgebra(M) <> R or ActingAlgebra(N) <> R then
    Error("<ring> must be right acting domain of <M> and <N>");
  fi;

  F := LeftActingDomain(R);
  zero := Zero(F);
  MBasis:=Basis(M);
  NBasis:=Basis(N);
  Mfam := FamilyObj(ExtRepOfObj(MBasis[1]));
  Nfam := FamilyObj(ExtRepOfObj(NBasis[1]));
  quiver := QuiverOfPathRing(R);
  Mdims := Mfam!.vertices;
  Ndims := Nfam!.vertices;
  Mmaps := Mfam!.matrices;
  Nmaps := Nfam!.matrices;
  l := Length(Mdims);
  arrows := ArrowsOfQuiver(quiver);
  
  # count the number of columns and rows
  numColumns := 0;
  startCol := [];

  for i in [1..l] do
      startCol[i] := numColumns + 1;
      numColumns := numColumns + Mdims[i]*Ndims[i];
  od;

  numRows := 0;

  for a in arrows do
    sourcePos := SourceOfPath(a)!.gen_pos;
    targetPos := TargetOfPath(a)!.gen_pos;
    numRows := numRows + Mdims[sourcePos]*Ndims[targetPos];
  od;
  
  equations := MutableNullMat(numRows, numColumns, F);
  
  # Create blocks
  startRow := 1;
  for a in arrows do
    arrowPos := a!.gen_pos;
    sourcePos := SourceOfPath(a)!.gen_pos;
    targetPos := TargetOfPath(a)!.gen_pos;

    # We always have a transposed version of the maps for $N$
    tMat := TransposedMat(Nmaps[arrowPos]);
    r := DimensionsMat(tMat)[1];
    c := DimensionsMat(tMat)[2];
    for i in [1..Mdims[sourcePos]] do
      row := startRow+(i-1)*r;
      col := startCol[sourcePos]+(i-1)*c;
      equations{[row..row+r-1]}{[col..col+c-1]} := tMat;
    od;

    # now subtract out appropriate copies of the $M$ map
    row := startRow;
    for i in [1..Mdims[sourcePos]] do
      for j in [1..Ndims[targetPos]] do
        col := startCol[targetPos] + j - 1;
        for k in [1..Mdims[targetPos]] do
          equations[row][col] := equations[row][col] - Mmaps[arrowPos][i][k];
          col := col + r;
        od;
        row := row + 1;
      od;
    od;
    startRow := startRow + Mdims[sourcePos]*Ndims[targetPos];
  od;
  
  ns := NullspaceMat(TransposedMat(equations));
  Msum:=Sum(Mdims);
  Nsum:=Sum(Ndims);
  basis := [];
  for b in ns do
    mat:=MutableNullMat(Msum, Nsum, F);
    k := 1;
    r:=0;
    c:=0;

    for i in [1..l] do
      for x in [1..Mdims[i]] do
        for y in [1..Ndims[i]] do
          mat[r+x][c+y] := b[k];
          k := k + 1;
        od;
      od;
      r:=r+Mdims[i];
      c:=c+Ndims[i];
    od;
    Add(basis, ShallowCopy(mat));
  od;

  maps:=List(basis, x->LeftModuleHomomorphismByMatrix(MBasis, x, NBasis));
  for x in maps do
    SetFilterObj(x, IsAlgebraModuleHomomorphism);
  od;

  return maps;
end;


InstallOtherMethod( Hom,
  "for a path ring, and two right modules over the path ring",
  true,
  [ IsPathRing, 
    IsPathAlgebraMatModule,
    IsPathAlgebraMatModule ], 0,
  function(R, M, N)
    local basis, Mbasis, Nbasis, gens, Msub, Nsub;

    if not IsIdenticalObj(R, ActingAlgebra(M)) or not IsIdenticalObj(R, ActingAlgebra(N))  then
      Error("<ring> must be right acting domain of <M> and <N>");
    fi;
        
    if not HasParent(M) then
      if not HasParent(N) then
        basis:=HOPF_MatrixModuleHom(R,M,N);
      else
        Mbasis:=CanonicalBasis(M);
        Nbasis:=CanonicalBasis(N);
        Nsub:=SubmoduleAsModule(N);
        gens:=HOPF_MatrixModuleHom(R, M, Nsub);
        basis:=List(gens, x->LeftModuleHomomorphismByMatrix(Mbasis, ShallowCopy(x!.matrix), Nbasis));
      fi;
    else
      if not HasParent(N) then
        Mbasis:=CanonicalBasis(M);
        Nbasis:=CanonicalBasis(N);
        Msub:=SubmoduleAsModule(M);
        gens:=HOPF_MatrixModuleHom(R, Msub, N);
        basis:=List(gens, x->LeftModuleHomomorphismByMatrix(Mbasis, ShallowCopy(x!.matrix), Nbasis));
      else
        Mbasis:=CanonicalBasis(M);
        Nbasis:=CanonicalBasis(N);
        Msub:=SubmoduleAsModule(M);
        Nsub:=SubmoduleAsModule(N);
        gens:=HOPF_MatrixModuleHom(R, Msub, Nsub);
        basis:=List(gens, x->LeftModuleHomomorphismByMatrix(Mbasis, ShallowCopy(x!.matrix), Nbasis));
      fi;
    fi;

    return VectorSpace(LeftActingDomain(R), basis, "basis");
  end
);


InstallOtherMethod( End,
  "for a path ring and a right module",
  true,
  [IsPathRing, IsPathAlgebraMatModule], 0,
  function( R, M )
    local basis, gens, Msub, Mbasis;

    if not IsIdenticalObj(R, ActingAlgebra(M)) then
      Error("<R> must be the right acting domain of <M>");
    fi;

    if not HasParent(M) then
      basis:=HOPF_MatrixModuleHom(R,M,M);
    else
      Mbasis:=CanonicalBasis(M);
      Msub:=SubmoduleAsModule(M);
      gens:=HOPF_MatrixModuleHom(R, Msub, Msub);
      basis:=List(gens, x->LeftModuleHomomorphismByMatrix(Mbasis, ShallowCopy(x!.matrix), Mbasis));
    fi;

    return AlgebraWithOne(LeftActingDomain(R), basis,"basis");
  end
);
