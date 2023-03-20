# GAP Implementation
# This file was generated from 
# $Id: pamodule.gi,v 1.20 2012/08/01 16:01:10 sunnyquiver Exp $

#######################################################################
##
#O  IsInFullMatrixRing( <M>, <R> )
##
##  This function returns true is the matrix  <M>  is in the full
##  matrix ring given by the field  <R>, and false otherwise.
##
InstallMethod ( IsInFullMatrixRing, 
    "for a matrix and a ring",
    true,
    [ IsMatrix, IsRing ], 
    0,
    function( M, R )

    local temp;
    
    temp := Flat( M ); 
    if ForAll( temp, x -> x in R ) then 
        return true;
    else
        return false;
    fi;
end
  );

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

InstallOtherMethod( \^,
 "for a module element and a path",
 true,
 [IsPathModuleElem, IsElementOfQuotientOfPathAlgebra ] , 0 ,
   function(elem, path);

   return elem^path![1];
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

InstallMethod( ExtRepOfObj, 	 
   "for elements of path algebra modules", 	 
   true, 	 
   [ IsPathModuleElem ], 0, 	 
   function( obj ) 	 
      return obj![1]; 	 
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
    if Length(gens) = 0 then
       MB := Objectify( NewType( FamilyObj( V ),
                             IsBasis and
                             IsBasisOfPathModuleElemVectorSpace and
                             IsAttributeStoringRep ),
                    rec() );
       SetUnderlyingLeftModule( MB, V );
       SetBasisVectors( MB, [] );
       SetIsEmpty(MB,true);
    else 
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
    fi;

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


###################################################
##
#O RightModuleOverPathAlgebra( <A>, <gens> )
##
## A=kQ/I is the algebra over which we want to construct
## the module. <gens> is a list of matrices, one for each 
## arrow of Q. The method checks that the relations
## are satisfied and if all matrices are over the correct
## field.
##
InstallMethod(RightModuleOverPathAlgebra,
    "for a (quotient of a) path algebra and list of matrices",
    true,
    [IsQuiverAlgebra, IsCollection], 0,
    function( R, gens )
    local tempgens, a, dim, source, target, basis, i, x, Fam, 
          vertices, matrices, quiver, M, vlist, alist, K, dim_M, 
          dim_vector, relationtest, I, terms, result, walk, matrix;

    matrices := [];
    quiver := QuiverOfPathRing(R);
    vlist  := VerticesOfQuiver(quiver);
    K      := LeftActingDomain(R);
    alist  := ArrowsOfQuiver(quiver);
    vertices := [];
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
        tempgens := Filtered(gens, g -> IsList(g[2][1]));
        if not ForAll(tempgens, g -> IsInFullMatrixRing(g[2],K)) then
            Error("not all matrices are over the correct field,\n");
        fi;
        for i in [ 1..Length( gens ) ] do
            a := gens[ i ][ 1 ];
            matrices[ quiver.(a)!.gen_pos ] := gens[ i ][ 2 ];
        od;
#
#  Input of the form [[matrix_1],[matrix_2],...,[matrix_n]]
#
    elif IsMatrix(gens[1]) then
        if not ForAll( gens, g -> IsInFullMatrixRing( g, K ) ) then
            Error("not all matrices are over the correct field,\n");
        fi;
        
        for i in [1 .. Length ( gens )] do
            matrices[i + Length(vlist)]:=gens[i];
        od;
    else
#
#  Input of the form [[alist[1],[matrix_1]],...,[alist[n],[matrix_n]]] 
#  where alist is a list of the arrows in the quiver.
#
        tempgens := Filtered(gens, g -> IsList(g[2][1]));
        if not ForAll( tempgens, g -> IsInFullMatrixRing( g[ 2 ], K) ) then
            Error("not all matrices are over the correct field,\n");
        fi;
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
  #
  # We might have isolated vertices, so no vectorspaces has been assign to these vertices.
  # In this variant of the function (not using dimension vector) we assume that the vector-
  # spaces are zero in these vertices.  These vertices are exactly those where the variable
  # vertices has value -1.
  #
  for i in [ 1..Length( vertices ) ] do
      if vertices[ i ] = -1 then
          vertices[ i ] := 0;
      fi;
  od;
#
#  Testing if the relations are satisfied, whenever we have a quotient of a path algebra.
#
    dim_M := 0;
    for i in [1..Length(vlist)] do
    	dim_M := dim_M + vertices[i];
    od;
   if dim_M = 0 or IsPathAlgebra(R) then 
      relationtest := true;
   else 
      relationtest := true;
      I := RelatorsOfFpAlgebra(R);
      for i in [1..Length(I)] do
         terms := CoefficientsAndMagmaElements(I[i]);
         result := [Zero(K)];
         for i in [1,3 .. Length( terms ) - 1] do
            walk := WalkOfPath(terms[i]);
            matrix := One(K);
            for x in walk do
               matrix := matrix * matrices[x!.gen_pos];
            od;
            matrix := terms[i+1] * matrix;
            result := result + matrix;
         od;
         dim := DimensionsMat(result);
         if result <> NullMat(dim[1],dim[2],K) then
            relationtest := false;
         fi;
      od;
   fi; 
#
# Creating the module if everything is OK, else give an error message.
#
   if relationtest then 
       Fam := NewFamily( "PathAlgModuleElementsFamily", IsPathModuleElem );
       SetFilterObj( Fam, IsPathModuleElemFamily );
       Fam!.vertices := vertices;
       Fam!.matrices := matrices;
       Fam!.pathAlgebra := OriginalPathAlgebra(R);
       if IsQuotientOfPathAlgebra(R) then 
          Fam!.quotientAlgebra := R;
       fi;
       if dim_M > 0 then 
          basis := CreateModuleBasis(Fam);
          dim_vector := vertices;
       else 
    	  basis := [ PathModuleElem( Fam, ZeroModElement( Fam, Zero( K ) ) ) ];
          dim_vector := List(vlist, x -> 0);
       fi;
       M := RightAlgebraModuleByGenerators(R, \^, basis);
       SetIsPathAlgebraMatModule(M,true);
       SetIsWholeFamily(M,true);
       SetDimensionVector(M,dim_vector);

       return M;
   else
      Print("The entered matrices for the module do not satisfy the relation(s).\n");
      return fail;
   fi;
end
);

###################################################
##
#O RightModuleOverPathAlgebraNC( <A>, <gens> )
##
## A=kQ/I is the algebra over which we want to construct
## the module. <gens> is a list of matrices, one for each 
## arrow of Q. The method checks that the relations
## are satisfied.
##
InstallMethod(RightModuleOverPathAlgebraNC,
  "for a (quotient of a) path algebra and list of matrices",
  true,
  [IsQuiverAlgebra, IsCollection], 0,
  function( R, gens )
    local a, dim, source, target, basis, i, x, Fam, 
          vertices, matrices, quiver, M, vlist, alist, K, dim_M, 
          dim_vector, relationtest, I, terms, result, walk, matrix;

    matrices := [];
    quiver := QuiverOfPathRing(R);
    vlist  := VerticesOfQuiver(quiver);
    K      := LeftActingDomain(R);
    alist  := ArrowsOfQuiver(quiver);
    vertices := [];
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
#
#  Testing if the relations are satisfied, whenever we have a quotient of a path algebra.
#
    dim_M := 0;
    for i in [1..Length(vlist)] do
    	dim_M := dim_M + vertices[i];
    od;
   if dim_M = 0 or IsPathAlgebra(R) then 
      relationtest := true;
   else 
      relationtest := true;
      I := RelatorsOfFpAlgebra(R);
      for i in [1..Length(I)] do
         terms := CoefficientsAndMagmaElements(I[i]);
         result := [Zero(K)];
         for i in [1,3 .. Length( terms ) - 1] do
            walk := WalkOfPath(terms[i]);
            matrix := One(K);
            for x in walk do
               matrix := matrix * matrices[x!.gen_pos];
            od;
            matrix := terms[i+1] * matrix;
            result := result + matrix;
         od;
         dim := DimensionsMat(result);
         if result <> NullMat(dim[1],dim[2],K) then
            relationtest := false;
         fi;
      od;
   fi; 
#
# Creating the module if everything is OK, else give an error message.
#
   if relationtest then 
       Fam := NewFamily( "PathAlgModuleElementsFamily", IsPathModuleElem );
       SetFilterObj( Fam, IsPathModuleElemFamily );
       Fam!.vertices := vertices;
       Fam!.matrices := matrices;
       Fam!.pathAlgebra := OriginalPathAlgebra(R);
       if IsQuotientOfPathAlgebra(R) then 
          Fam!.quotientAlgebra := R;
       fi;
       if dim_M > 0 then 
          basis := CreateModuleBasis(Fam);
          dim_vector := vertices;
       else 
    	  basis := [ PathModuleElem( Fam, ZeroModElement( Fam, Zero( K ) ) ) ];
          dim_vector := List(vlist, x -> 0);
       fi;
       M := RightAlgebraModuleByGenerators(R, \^, basis);
       SetIsPathAlgebraMatModule(M,true);
       SetIsWholeFamily(M,true);
       SetDimensionVector(M,dim_vector);

       return M;
   else
      Print("The entered matrices for the module do not satisfy the relation(s).\n");
      return fail;
   fi;
end
);


#######################################################################
##
#O  RightModuleOverPathAlgebra( <A>, <dim_vector>, <gens> )
##
##  This function constructs a right module over a (quotient of a) path
##  algebra  A  with dimension vector  dim_vector, and where the 
##  generators with a non-zero action is given in the list  gens. The 
##  format of the list gens is [["a",[matrix_a]],["b",[matrix_b]],...],
##  where "a" and "b" are labels of arrows used when the underlying quiver
##  was created and matrix_? is the action of the algebra element 
##  corresponding to the arrow with label "?". The action of the arrows
##  can be entered in any order. The function checks if the matrices are
##  over the correct field.
##
InstallOtherMethod(RightModuleOverPathAlgebra,
  "for a path algebra and list of matrices",
  true,
  [ IsQuiverAlgebra, IsList, IsList ], 0,
  function( A, dim_vector, gens )

  local quiver, vertices, K, arrows, num_vert, gens_of_quiver,
        matrices, i, entered_arrows, arrow_labels, matrices_set, a, b,
        origin, target, dim_mat, arrows_not_set, dim_M, relationtest,
        I, terms, result, matrix, walk, x, dim, Fam, basis, M;
#
#  Setting up the data we need.
#
   quiver   := QuiverOfPathRing(A);
   vertices := VerticesOfQuiver(quiver);
   K        := LeftActingDomain(A);
   arrows   := ArrowsOfQuiver(quiver);
   num_vert := NumberOfVertices(quiver);
   gens_of_quiver := GeneratorsOfQuiver(quiver);
#
#  Setting the multiplication by the vertices.
#          
   matrices := List([1 .. Length(vertices)], x -> One(K));
#
#  First finding the labels of the entered by the user,
#  then finding the labels used for the quiver.
#
   entered_arrows := List(gens, x -> x[1]);
   arrow_labels := List(arrows, x -> x!.String);    
   matrices_set := [];
#
#  If all labels entered by the user actually exist in the quiver in question,
#  start creating the matrices defining the multiplication by the arrows by 
#  first checking that they are in the full matrix ring over the correct field.
#
   if not ForAll(gens, g -> IsInFullMatrixRing(g[2],K)) then
       Error("not all matrices are over the correct field,\n");
   fi;
   if ForAll( entered_arrows, x -> x in arrow_labels )  then                 
      for i in [1 .. Length( gens )] do
         a := gens[i][1];
         b := quiver.(a);
         origin := Position(vertices,SourceOfPath(b));
         target := Position(vertices,TargetOfPath(b));
# 
#  If an entered map for an arrow ending or starting in a zero dimensional space is non-zero, 
#  give an error message, otherwise check for the right size and enter into list of acting 
#  matrices if OK.
# 
         if ( dim_vector[origin] = 0 or dim_vector[target] = 0 ) and gens[i][2] <> Zero(gens[i][2]) then  
            Error("a non-zero matrix was entered for the arrow  ",b,"  when it should be zero.\n"); 
         else
            dim_mat := DimensionsMat(gens[i][2]);
            if ( dim_vector[origin] <> dim_mat[1] ) or (dim_vector[target] <> dim_mat[2] ) then   
               Error("wrong size of the matrix for the arrow  ",b,".\n");
            else 
               matrices[quiver.(a)!.gen_pos] := gens[i][2];
               AddSet(matrices_set, quiver.(a)!.gen_pos);
            fi;
         fi;
      od;
#
#  Now set the matrices which are not set by the user, that is, all are going to be zero matrices. 
#
      arrows_not_set := [num_vert+1..num_vert+NumberOfArrows(quiver)];
      SubtractSet(arrows_not_set, matrices_set);
      for i in arrows_not_set do
         a := gens_of_quiver[i];
         origin := Position(vertices,SourceOfPath(a));
         target := Position(vertices,TargetOfPath(a)); 
         if dim_vector[origin] = 0 then 
            if dim_vector[target] = 0 then 
               matrices[i] := NullMat(1,1,K);
            else
               matrices[i] := NullMat(1,dim_vector[target],K);
            fi;
         else
            if dim_vector[target] = 0 then 
               matrices[i] := NullMat(dim_vector[origin],1,K);
            else
               matrices[i] := NullMat(dim_vector[origin],dim_vector[target],K);
            fi;            
         fi;
      od;
   else
      Print("The following arrow labels don't exist: ",Filtered(entered_arrows, x -> not x in arrow_labels),"\n");
      Error("input not on required form.\n");
   fi;

#
#  Testing if the relations are satisfied, whenever we have a quotient of a path algebra.
#
   dim_M := Sum(dim_vector);
   if dim_M = 0 or IsPathAlgebra(A) then 
      relationtest := true;
   else 
      relationtest := true;
      I := RelatorsOfFpAlgebra(A);
      for i in [1..Length(I)] do
         terms := CoefficientsAndMagmaElements(I[i]);
         result := [Zero(K)];
         for i in [1,3 .. Length( terms ) - 1] do
            walk := WalkOfPath(terms[i]);
            matrix := One(K);
            for x in walk do
               matrix := matrix * matrices[x!.gen_pos];
            od;
            matrix := terms[i+1] * matrix;
            result := result + matrix;
         od;
         dim := DimensionsMat(result);
         if result <> NullMat(dim[1],dim[2],K) then
            relationtest := false;
         fi;
      od;
   fi; 
#
# Creating the module if everything is OK, else give an error message.
#
   if relationtest then 
       Fam := NewFamily( "PathAlgModuleElementsFamily", IsPathModuleElem );
       SetFilterObj( Fam, IsPathModuleElemFamily );
       Fam!.vertices := dim_vector;
       Fam!.matrices := matrices;
       Fam!.pathAlgebra := OriginalPathAlgebra(A);
       if IsQuotientOfPathAlgebra(A) then 
          Fam!.quotientAlgebra := A;
       fi;       
       if dim_M > 0 then 
          basis := CreateModuleBasis(Fam);
       else 
    	  basis := [ PathModuleElem( Fam, ZeroModElement( Fam, Zero( K ) ) ) ];
       fi;
       M := RightAlgebraModuleByGenerators(A, \^, basis);
       SetIsPathAlgebraMatModule(M,true);
       SetIsWholeFamily(M,true);
       SetDimensionVector(M,dim_vector);

       return M;
   else
      Print("The entered matrices for the module do not satisfy the relation(s).\n");
      return fail;
   fi;
end
);

#######################################################################
##
#O  \=( <M>, <N> )
##
##  This function returns true if the modules  <M>  and  <N>  have
##  the same underlying vector spaces and the same defining matrices, 
##  and false otherwise.
##
InstallMethod ( \=, 
    "for a PathAlgebraMatModules",
    true,
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
    0,
    function( M, N );

    if ( RightActingAlgebra(M) = RightActingAlgebra(N) ) and 
       ( DimensionVector(M) = DimensionVector(N) ) and
       ( MatricesOfPathAlgebraModule(M) = MatricesOfPathAlgebraModule(N) ) then
        return true;
    else
        return false;
    fi;
end
);

InstallMethod( ViewObj, 
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule ], NICE_FLAGS + 1,
    function (M);

    Print("<",DimensionVector(M),">");
end
);
  
InstallMethod( PrintObj, 
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule ], NICE_FLAGS + 1,
    function (M)

    local A;

    A := RightActingAlgebra(M);
    Print("<Module over ");
    View(A);
    Print(" with dimension vector ", DimensionVector(M),">");
end
); 

InstallMethod( Display, 
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule ], NICE_FLAGS+1,
    function ( M )
    
    local Q, arrows, a;

    Print(M);
    Print(" and linear maps given by\n");
    Q := QuiverOfPathAlgebra(RightActingAlgebra(M));
    arrows := ArrowsOfQuiver(Q);
    for a in arrows do
        Print("for arrow ",a,":\n");
        PrintArray(MatricesOfPathAlgebraModule(M)[Position(arrows,a)]);
    od;
end
); 


#######################################################################
##
#O  \in ( <v>, <V> )
##
##  A test to check whether a given module element v is in the 
##  module V.
## 
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

#######################################################################
##
#O  SubAlgebraModule ( <V>, gens )
##
##  Given a module V and a list of module elements gens, this
##  method returns the module generated by the elements of gens.
## 
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

#######################################################################
##
#O  SubmoduleAsModule ( <N> )
##
##  Given a IsPathAlgebraMatModule N, this function returns a 
##  module objected constructed with the 
##  RightModuleOverPathAlgebra method.
##
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

#######################################################################
##
#O  NaturalHomomorphismBySubAlgebraModule ( <M>,<N> )
##
##  N is a submodule of M. Returns the natural homomorphism 
##  M --> N.  ?
##
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
  
  equations := NullMat(numRows, numColumns, F);
  
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
    mat:=NullMat(Msum, Nsum, F);
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

#######################################################################
##
#O  Hom( <R>, <M>, <N> )
##
##  R is an algebra, M and N two modules. Returns Hom_R(M,N) as
##  a vector space.
##
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

#######################################################################
##
#O  End( <R>, <M> )
##
##  For a k-algebra R and an R-module M. Returns End_R(M) as a k-algebra.
##  
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

#######################################################################
##
#O  MatricesOfPathAlgebraModule( <M> )
##
##  For an R-module M. Returns a list of the matrices giving the
##  linear transformations of M corresponding to arrows of Q, where
##  R = kQ/I.
##
InstallMethod ( MatricesOfPathAlgebraModule,
    "for a representation of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M )
#
#   M     = a representation of the quiver Q over K
#
        local n, i, A, Q, num_vert, num_arrow, m, fam, mat;

    if Dimension(M) = 0 then
       n := Length(ArrowsOfQuiver(QuiverOfPathAlgebra(RightActingAlgebra(M)))); 
       mat := [];
       for i in [1..n] do
          mat[i] := 0;
       od;

       return mat;      
    else
       A := RightActingAlgebra(M);
       Q := QuiverOfPathRing(A);
       num_vert  := Length(VerticesOfQuiver(Q));
       num_arrow := Length(ArrowsOfQuiver(Q));
       m   := ExtRepOfObj(Zero(M));
       fam := FamilyObj(m);

       return fam!.matrices{[num_vert+1..num_arrow+num_vert]};
    fi;
end
);

#######################################################################
##
#A  MinimalGeneratingSetOfModule( <M> )
##
##  This function computes a minimal generating set of the module M. 
##
InstallMethod( MinimalGeneratingSetOfModule,
  "for a path algebra module",
  true,
  [ IsPathAlgebraMatModule ], 0,
  function( M )

  local f;

    f := TopOfModuleProjection(M);
    return List(BasisVectors(Basis(Range(f))), x -> PreImagesRepresentative(f,x));
end
);

InstallMethod( DimensionVectorPartialOrder, 
   "for two path algebra matmodules",
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
   function( M, N) 

   local L1, L2;
 
   L1 := DimensionVector(M);
   L2 := DimensionVector(N);

   return ForAll(L2-L1,x->(x>=0));
end
);

#######################################################################
##
#A  AnnihilatorOfModule( <M> )
##
##  Given a module  M  over a (quotient of a) path algebra  A, this 
##  function computes a vectorspace basis for the annihilator of  M
##  in  A. 
##
InstallMethod ( AnnihilatorOfModule, 
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule ], 
    0,
    function( M )

    local A, B, fam, gb, nontips, matrix, n, temp, m,
          solutions, annihilator;

    A := RightActingAlgebra(M);
    #   
    #  If the module  M  is zero, return the whole algebra.
    #
    if Dimension(M) = 0 then 
        return BasisVectors(Basis(A));
    fi;
    B := BasisVectors(Basis(M));
    #
    #  Setting things up right according to the input.
    #
    if IsPathAlgebra(A) then 
        nontips := BasisVectors(Basis(A));
    elif IsQuotientOfPathAlgebra(A) then 
        fam := ElementsFamily(FamilyObj(A));
        gb := GroebnerBasisOfIdeal(fam!.ideal);
        nontips := List(Nontips(gb), x -> x*One(A));
    else
        Error("the module is not a module over a (quotient of a) path algebra,");
    fi;
    #
    #  Computing the linear system to solve in order to find the annihilator
    #  of the module  M.
    #
    matrix := [];
    for n in nontips do
        temp := [];
        for m in B do
            Add(temp,ExtRepOfObj(ExtRepOfObj(m^n)));
        od;
        Add(matrix,temp);
    od;

    matrix := List(matrix, x -> Flat(x));
    #
    #  Finding the solutions of the linear system, and creating the solutions
    #  as elements of the algebra  A.
    #
    solutions := NullspaceMat(matrix);
    annihilator := List(solutions, x -> LinearCombination(nontips,x));

    return annihilator;
end
);

#######################################################################
##
#O  LoewyLength ( <M> )
##
##  This function returns the Loewy length of the module  M, for a 
##  module over a (quotient of a) path algebra.
##
InstallMethod( LoewyLength, 
   "for a PathAlgebraMatModule",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) 

   local N, i;

   N := M;
   i := 0;
   if Dimension(M) = 0 then 
      return 0;
   else
      repeat 
         N := RadicalOfModule(N);
         i := i + 1;
      until
         Dimension(N) = 0;
      return i;
   fi;
end
);

#######################################################################
##
#O  LoewyLength ( <A> )
##
##  This function returns the Loewy length of the algebra  A, for a 
##  finite dimensional (quotient of a) path algebra (by an admissible
##  ideal).
##
InstallOtherMethod( LoewyLength, 
    "for (a quotient of) a path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A ) 

    local fam, N;
    
    if not IsFiniteDimensional(A) then
        Error("the entered algebra is not finite dimensional,\n");
    fi;
    if not IsPathAlgebra(A) and not IsAdmissibleQuotientOfPathAlgebra(A) then 
        TryNextMethod();
    fi;

    N := IndecProjectiveModules(A);
    N := List(N, x -> LoewyLength(x));
    return Maximum(N);
end
);

#######################################################################
##
#O  RadicalSeries ( <M> )
##
##  This function returns the radical series of the module  M, for a 
##  module over a (quotient of a) path algebra. It returns a list of 
##  dimension vectors of the modules: [ M/rad M, rad M/rad^2 M, 
##  rad^2 M/rad^3 M, .....].
##
InstallMethod( RadicalSeries, 
   "for a PathAlgebraMatModule",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) 

   local N, radlayers, i;

   if Dimension(M) = 0 then 
      return DimensionVector(M);
   else
      radlayers := [];
      i := 0;
      N := M;
      repeat     
         Add(radlayers,DimensionVector(TopOfModule(N)));
         N := RadicalOfModule(N);
         i := i + 1;
      until
         Dimension(N) = 0;
      return radlayers;
   fi;
end
);

#######################################################################
##
#O  SocleSeries ( <M> )
##
##  This function returns the socle series of the module  M, for a 
##  module over a (quotient of a) path algebra. It returns a list of 
##  dimension vectors of the modules: [..., soc(M/soc^3 M), 
##  soc(M/soc^2 M), soc(M/soc M), soc M].
##
InstallMethod( SocleSeries, 
   "for a PathAlgebraMatModule",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) 

   local N, series, i;

   if Dimension(M) = 0 then 
      return DimensionVector(M);
   else
      series := [];
      i := 0;
      N := DualOfModule(M);
      repeat     
         Add(series,DimensionVector(TopOfModule(N)));
         N := RadicalOfModule(N);
         i := i + 1;
      until
         Dimension(N) = 0;
      
      return Reversed(series);
   fi;
end
);

#######################################################################
##
#A  Dimension( <M> )
##
##  Returns the k-dimension of the module <M>.
##
InstallOtherMethod( Dimension,
   "for a PathAlgebraMatModule",
   [ IsPathAlgebraMatModule ], 0,
   function( M );

   return Sum(DimensionVector(M));
end
);

#######################################################################
##
#P  IsProjectiveModule( <M> )
##
##  Checks whether <M> is projective.
##
InstallMethod( IsProjectiveModule, 
   "for a module over a quotient of a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) 

   local top, dimension, i, P; 

   top := TopOfModule(M); 
   dimension := 0; 
   P := IndecProjectiveModules(RightActingAlgebra(M)); 
   for i in [1..Length(DimensionVector(M))] do 
      if DimensionVector(top)[i] <> 0 then 
         dimension := dimension + Dimension(P[i])*DimensionVector(top)[i];
      fi;
   od; 

   if dimension = Dimension(M) then 
      return true;
   else 
      return false;
   fi;
end
);

#######################################################################
##
#P  IsInjectiveModule( <M> )
##
##  Checks whether <M> is injective.
##
InstallMethod( IsInjectiveModule, 
   "for a module over a quotient of a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) ; 

   return IsProjectiveModule(DualOfModule(M));
end
);

#######################################################################
##
#P  IsSimpleQPAModule( <M> )
##
##  Checks whether <M> is simple.
##
InstallMethod( IsSimpleQPAModule, 
   "for a module over a quotient of a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) ; 

   if Dimension(M) = 1 then 
      return true;
   else
      return false;
   fi;
end
);

#######################################################################
##
#P  IsSemisimpleModule( <M> )
##
##  Checks whether <M> is semisimple.
##
InstallMethod( IsSemisimpleModule, 
   "for a module over a quotient of a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) ; 

   if Dimension(RadicalOfModule(M)) = 0 then 
      return true;
   else
      return false;
   fi;
end
);

#######################################################################
##
#M  DirectSumOfQPAModules( <L> )
##
##  <L> is a list of modules over a path algebra.  The function computes
##  and returns the direct sum of all modules in <L>.  The projections
##  and inclusions between the modules in <L> and the direct sum is stored
##  as attributes of the direct sum.
##
InstallMethod( DirectSumOfQPAModules,
   "for a list of modules over a path algebra",
   [ IsList ], 0,
   function( L ) 

    local   n,  A,  K,  Q,  arrows,  vertices,  dim_list,  dim_vect,  
            i,  temp,  j,  big_mat,  a,  origin,  target,  mat,  
            row_pos,  col_pos,  r,  direct_sum,  list_of_projs,  map,  
            maps,  dimension,  list_of_incls;

   n := Length( L );
   if n > 0 then 
      A := RightActingAlgebra( L[ 1 ] );
      K := LeftActingDomain( A );
      if ForAll( L, IsPathAlgebraMatModule ) and ForAll( L, x -> RightActingAlgebra( x ) = A ) then 
         Q := QuiverOfPathAlgebra( OriginalPathAlgebra( A ) );
         arrows := ArrowsOfQuiver( Q );
         vertices := VerticesOfQuiver( Q );
         dim_list := List( L, DimensionVector );
         dim_vect := Sum( dim_list ); 
         big_mat := [ ];
         i := 0;
         for a in arrows do
            i := i + 1;
            origin := Position( vertices, SourceOfPath( a ) );
            target := Position( vertices, TargetOfPath( a ) );
            mat := [ ];
            if dim_vect[ origin ] <> 0 and dim_vect[ target ] <> 0 then 
               mat := NullMat( dim_vect[ origin ], dim_vect[ target ], K );
               row_pos := 1;
               col_pos := 1;
               for r in [ 1..n ] do
                  if dim_list[ r ][ origin ] <> 0 and dim_list[ r ][ target ] <> 0 then
                     mat{ [ row_pos..( row_pos + dim_list[ r ][ origin ] - 1 ) ] }{ [ col_pos..( col_pos + dim_list[ r ][ target ] - 1 ) ] } := 
                       MatricesOfPathAlgebraModule( L[ r ] )[ i ];
                  fi;
                  row_pos := row_pos + dim_list[ r ][ origin ];
                  col_pos := col_pos + dim_list[ r ][ target ]; 
               od;
               mat := [ String( a ), mat ];
               Add( big_mat, mat );
            fi;
         od;
      fi;
      direct_sum := RightModuleOverPathAlgebra( A, dim_vect, big_mat );

      list_of_projs := [];
      for i in [ 1..n ] do 
         map := [ ];
         maps := [ ];
         for j in [ 1..Length( dim_vect ) ] do
            if dim_vect[ j ] = 0 then 
               if dim_list[ i ][ j ] = 0 then 
                  map := NullMat( 1, 1, K );
               else
                  map := NullMat( 1, dim_list[ i ][ j ], K );
               fi;
            else
               if dim_list[ i ][ j ] = 0 then 
                  map := NullMat( dim_vect[ j ], 1, K);
               else
                  dimension := 0;
                  if i > 1 then 
                     for r in [ 1..i-1 ] do 
                        dimension := dimension + dim_list[ r ][ j ];
                     od;
                  fi;
                  map := NullMat( dim_vect[ j ], dim_list[ i ][ j ], K );
                  map{ [ dimension + 1..dimension + dim_list[ i ][ j ] ] }{ [ 1..dim_list[ i ][ j ] ] } := IdentityMat( dim_list[ i ][ j ], K );
               fi;
            fi;
            Add( maps, map );
         od;
         Add( list_of_projs, RightModuleHomOverAlgebra( direct_sum, L[ i ], maps));
      od;

      list_of_incls := [ ];
      for i in [ 1..n ] do 
         map := [ ];
         maps := [ ];
         for j in [ 1..Length( dim_vect ) ] do
            if dim_vect[ j ] = 0 then 
               map := NullMat( 1, 1, K );
            else
               if dim_list[ i ][ j ] = 0 then 
                  map := NullMat( 1, dim_vect[ j ], K );
               else
                  dimension := 0;
                  if i > 1 then 
                     for r in [ 1..i - 1 ] do 
                        dimension := dimension + dim_list[ r ][ j ];
                     od;
                  fi;
                  map := NullMat( dim_list[ i ][ j ], dim_vect[ j ], K );
                  map{ [ 1..dim_list[ i ][ j ] ] }{ [ dimension + 1..dimension + dim_list[ i ][ j ] ] } := IdentityMat( dim_list[ i ][ j ], K );
               fi;
            fi;
            Add( maps, map );
         od;
         Add( list_of_incls, RightModuleHomOverAlgebra( L[ i ], direct_sum, maps ) );
      od;

      SetFilterObj( direct_sum, IsDirectSumOfModules );
      SetDirectSumProjections( direct_sum, list_of_projs );
      SetDirectSumInclusions( direct_sum, list_of_incls );

      return direct_sum;
   fi;

   return fail;
end
);


#######################################################################
##
#O  SupportModuleElement(<m>)
##
##  Given an element  <m>  in a PathAlgebraMatModule, this function
##  finds in which vertices this element is support, that is, has 
##  non-zero coordinates. 
##
InstallMethod ( SupportModuleElement, 
   "for an element in a PathAlgebraMatModule",
   true,
   [ IsRightAlgebraModuleElement ],
   0,
   function( m )

   local fam, A, Q, idempotents, support;
#
#  Finding the algebra over which the module  M, in which  m is an 
#  element, is a module over. 
# 
   if IsPathModuleElem(ExtRepOfObj(m)) then 
      fam := CollectionsFamily(FamilyObj(m))!.rightAlgebraElementsFam;   
      if "wholeAlgebra" in NamesOfComponents(fam) then 
         A := fam!.wholeAlgebra;
      elif "pathRing" in NamesOfComponents(fam) then 
         A := fam!.pathRing;
      fi;
   else 
      return fail;
   fi;
#
#  Finding the support of m.
#
   Q := QuiverOfPathAlgebra(OriginalPathAlgebra(A));
   idempotents := List(VerticesOfQuiver(Q), x -> x*One(A));
   support := Filtered(idempotents, x -> m^x <> Zero(m));

   return support;
end
);

#######################################################################
##
#O  RightAlgebraModuleToPathAlgebraMatModule( <M> ) 
##  
##  This function constructs a right module over a (quotient of a) path
##  algebra  A  from a RightAlgebraModule over the same algebra  A. The
##  function checks if  A  actually is a quotient of a path algebra and 
##  if the module  M  is finite dimensional. In either cases it returns
##  an error message. 
##
InstallMethod ( RightAlgebraModuleToPathAlgebraMatModule, 
    "for an RightAlgebraModule",
    true,
    [ IsRightAlgebraModuleElementCollection ],
    0,
    function( M )

    local A, vertices, num_vert, B, arrows, generators_in_vertices, 
          vertexwise_basis, i, j, vertices_Q, mat, arrows_as_paths, 
          a, partial_mat, source, target, source_index, target_index, 
          rows, cols, arrow, b, vector, dimvect;
    
    if not IsFiniteDimensional(M) then
        Error("the entered module is not finite dimensional,\n");
    fi;
    A := RightActingAlgebra(M); 
    if not ( IsPathAlgebra(A) or IsQuotientOfPathAlgebra(A) ) then 
        Error("the entered module is not a module over a quotient of a path algebra,\n");
    fi;
    vertices := VerticesOfQuiver(QuiverOfPathAlgebra(A));
    vertices := List(vertices, x -> x*One(A));
    num_vert := Length(vertices); 
    B := BasisVectors(Basis(M));
    arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(A));
    arrows := List(arrows, x -> x*One(A));    
    #
    #  Constructing a uniform basis of  M. 
    #
    generators_in_vertices := List(vertices, x -> []); 
    for i in [1..Length(B)] do
        for j in [1..num_vert] do
            if ( B[i]^vertices[j] <> Zero(M) ) then
                Add(generators_in_vertices[j],B[i]^vertices[j]);
            fi;
        od;
    od;    
    #
    #   Finding a K-basis of M in each vertex.
    #    
    vertexwise_basis := List(generators_in_vertices, x -> Basis(Subspace(M,x)));
    # 
    #   Finding the matrices defining the representation
    # 
    mat := [];
    vertices_Q:= VerticesOfQuiver(QuiverOfPathAlgebra(A)); 
    arrows_as_paths := List(ArrowsOfQuiver(QuiverOfPathAlgebra(A)), x -> One(A)*x);
    for a in arrows_as_paths do
        partial_mat := [];
        source := SourceOfPath(TipMonomial(a));
        target := TargetOfPath(TipMonomial(a));
        source_index := Position(vertices_Q,source);
        target_index := Position(vertices_Q,target);
        rows := Length(BasisVectors(vertexwise_basis[source_index]));
        cols := Length(BasisVectors(vertexwise_basis[target_index]));
        arrow := TipMonomial(a);
        if rows <> 0 and cols <> 0 then
            for b in BasisVectors(vertexwise_basis[source_index]) do
                vector := Coefficients(vertexwise_basis[target_index], b^a); 
                Add(partial_mat,vector);
            od;
            Add(mat,[String(arrow),partial_mat]);
        fi;
    od;
    
    dimvect := List(vertexwise_basis, b -> Length(b));
    
    return RightModuleOverPathAlgebra(A, dimvect, mat);
end
  );

#######################################################################
##
#P  IsRigidModule( <M> )
##
##  This function returns true if the entered module  <M>  is a rigid 
##  module, otherwise false.
##
InstallMethod( IsRigidModule,
    "for a PathAlgebraMatModule",
    [ IsPathAlgebraMatModule ], 0,
    function( M ); 
    
    return Length(ExtOverAlgebra(M,M)[2]) = 0;
end
  );

#######################################################################
##
#P  IsTauRigidModule( <M> )
##
##  This function returns true if the entered module  <M>  is a tau rigid 
##  module, otherwise false.
##
InstallMethod( IsTauRigidModule,
    "for a PathAlgebraMatModule",
    [ IsPathAlgebraMatModule ], 0,
    function( M ); 
    
    return Length(HomOverAlgebra(M,DTr(M))) = 0;
end
  );

#######################################################################
##
#P  IsIndecomposableModule( <M> )
##
##  This function returns true if the entered module  <M>  is an
##  indecomposable module, otherwise false if the field, over which
##  the algebra is defined, is finite. 
##
InstallMethod( IsIndecomposableModule,
    "for a PathAlgebraMatModule",
    [ IsPathAlgebraMatModule ], 0,
    function( M )
    
    local K; 

    if IsZero( M ) then
       return false;
    fi;
    K := LeftActingDomain(M);
    if IsFinite(K) then 
        return Length(DecomposeModule(M)) = 1;
    fi;
    if IsSimpleQPAModule( TopOfModule(M) ) then
        return true;
    fi;
    if IsSimpleQPAModule( SocleOfModule(M) ) then
        return true;
    fi;
    TryNextMethod();
end
  );

#######################################################################
##
#P  IsExceptionalModule( <M> )
##
##  This function returns true if the entered module  <M>  is an
##  exceptional module (ie. indecomposable and Ext^1(M,M)=(0), otherwise 
##  false, if the field, over which the algebra  <M>  is defined over, 
##  is finite.
##
InstallMethod( IsExceptionalModule,
    "for a PathAlgebraMatModule",
    [ IsPathAlgebraMatModule ], 0,
    function( M );
    
    return IsRigidModule(M) and IsIndecomposableModule(M);
end
  );

#######################################################################
##
#O  ComplexityOfModule( <M>, <n> )
##
##  Recall: If a function  f(x)  is a polynomial in  x, the degree of  
##  f(x)  is given by lim_{n --> infinity} log |f(n)|/log n. 
##  
##  Given a module  M  this function computes an estimate of the 
##  complexity of the module by approximating the complexity by 
##  considering the limit  lim_{m --> infinity} log dim(P(M)(m))/log m
##  where P(M)(m) is the m-th projective in a minimal projective 
##  resolution of  M  at stage  m.  This limit is estimated by 
##  log dim(P(M)(n))/log n.
##  
InstallMethod ( ComplexityOfModule, 
    "for a module of a fin. dim. quotient of a path algebra and a positive integer",
    true,
    [ IsPathAlgebraMatModule, IS_INT ], 
    0,
    function( M, n )

    local A, projres, temp;
    
    A := RightActingAlgebra(M);
    if HasIsFiniteGlobalDimensionAlgebra(A) and IsFiniteGlobalDimensionAlgebra(A) then
        return 0;
    fi;
    
    if IsAdmissibleQuotientOfPathAlgebra(A) then 
        projres := ProjectiveResolution(M); 
        temp := Dimension(ObjectOfComplex(projres,n))*1.0;
        if Dimension(ObjectOfComplex(projres,n)) = 0 then
            return 0;
        else
            temp := Log10(temp)/Log10(n*1.0);
            return Int(Round(temp + 1.0));
        fi;
    else
        Error("entered algebra is not finite dimensional or an admissible quotient of a path algebra,\n");
    fi;
end 
);

#######################################################################
##
#O  ComplexityOfAlgebra( <A>, <n> )
##
##  Recall: If a function  f(x)  is a polynomial in  x, the degree of  
##  f(x)  is given by lim_{n --> infinity} log |f(n)|/log n. 
##  
##  Given an algebra  A  this function computes an estimate of the 
##  maximal complexity of the simple modules by approximating the 
##  complexity of each simple  S  by  considering the limit  
##  lim_{m --> infinity} log dim(P(S)(m))/log m
##  where P(S)(m) is the m-th projective in a minimal projective 
##  resolution of a simple module  S  at stage  m.  This limit is 
##  estimated by  log dim(P(S)(n))/log n.
##  
InstallMethod ( ComplexityOfAlgebra, 
    "for a fin. dim. quotient of a path algebra and a positive integer",
    true,
    [ IsQuiverAlgebra, IS_INT ], 
    0,
    function( A, n )

    local S, temp;
    #
    #
    if HasIsFiniteGlobalDimensionAlgebra(A) and IsFiniteGlobalDimensionAlgebra(A) then
        return 0;
    fi;
    if IsAdmissibleQuotientOfPathAlgebra(A) then 
        S := SimpleModules(A);
        temp := List(S, s -> ComplexityOfModule(s,n));
        return Maximum(temp);
    else
        Error("entered algebra is not finite dimensional or an admissible quotient of a path algebra,\n");
    fi;
end 
  );

#######################################################################
##
#P  IsZero( <M> )
##
##  This function returns true if the entered module  <M>  is zero, 
##  otherwise it returns false. 
##
InstallOtherMethod( IsZero,
    "for a PathAlgebraMatModule",
    [ IsPathAlgebraMatModule ], 0,
    function( M );
    
    return Dimension(M) = 0;
end
  );

#######################################################################
##
#O  RestrictionViaAlgebraHomomorphism( < f, M > )
##
##  Given an algebra homomorphism  f : A ---> B and a module  M  over 
##  B, this function returns  M  as a module over  A.
##  
InstallMethod( RestrictionViaAlgebraHomomorphism, 
    "for a IsPathAlgebraMatModule and a IsAlgebraHomomorphism",
    [ IsAlgebraHomomorphism, IsPathAlgebraMatModule ], 0,
    function( f, M )

    local   K,  A,  B,  vertices,  V,  BV,  arrows,  mats,  a,  
            startpos,  endpos, dimvect;
    
    K := LeftActingDomain( M );
    A := Source( f );
    B := Range( f );
    if RightActingAlgebra( M ) <> B then 
        Error( "The entered module is not a module over the range of the algebra homomorphism.\n" );
    fi;
    vertices := One( A ) * VerticesOfQuiver( QuiverOfPathAlgebra( A ) );
    V := List( vertices, v -> List( BasisVectors( Basis( M ) ), m -> m ^ ImageElm( f, v ) ) ); 
    V := List( V, W -> Filtered( W, w -> w <> Zero( w ) ) );
    BV := List( V, function( W ) if IsEmpty(W) then return []; else return Basis( VectorSpace( K, W ) ); fi; end );
    arrows := ArrowsOfQuiver( QuiverOfPathAlgebra( A ) );
    mats := [ ]; 
    for a in arrows do
        startpos := Position( vertices, One( A ) * SourceOfPath ( a ) );
        endpos := Position( vertices, One( A ) * TargetOfPath ( a ) );        
        if not IsEmpty( BV[ startpos ] ) and not IsEmpty( BV[ endpos ] ) then
            Add( mats, [ String(a), List( BV[ startpos ], m -> Coefficients( BV[ endpos ], m ^ ImageElm( f, One( A ) * a ) ) ) ] );
        fi;
    od;
    
    dimvect := List(BV, b -> Length(b));
    
    return RightModuleOverPathAlgebra( A, dimvect, mats ); 
end
  );

