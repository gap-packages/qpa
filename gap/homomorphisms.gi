# GAP Implementation
# $Id: homomorphisms.gi,v 1.4 2010/11/19 13:24:48 sunnyquiver Exp $

InstallMethod( ImageElm, 
    "for a map between representations and an element in a representation.",
    [ IsPathAlgebraMatModuleMap, IsAlgebraModuleElement ], 0, 
    function( map, elem )

    local elt, i, n, fam, image, temp, zero, new_image;

    if elem in Source(map) then
  	    elt := ExtRepOfObj(elem);
        n := Basis(Range(map))[1];
        image := [];
        for i in [1..Length(elt![1])] do 
            image[i] := elt![1][i]*map!.maps[i];
        od;

        image := PathModuleElem(FamilyObj(Zero(Range(map))![1]),image);        
        return Objectify( TypeObj( n ), [ image ] );
    else
        Error("the element entered is not an element in the source of the map,");
    fi; 
end);

InstallMethod( ImagesSet, 
    "for a map between representations and a set of elements in a representation.",
    [ IsPathAlgebraMatModuleMap, IsCollection ], 0, 
    function( map, elts )

    local elt, B, images;

    images := [];
    if IsList(elts) then
        if IsFinite(elts) then
            B := elts;
        fi;
    else
        if IsPathAlgebraMatModule(elts) then 
            B := CanonicalBasis(elts);
        else
            Error("input of wrong type,");
        fi;
    fi;
    for elt in B do
       if ImageElm(map,elt) <> Zero(Range(map)) then 
          Add(images,ImageElm(map,elt));
       fi;
    od;
    return images;
end
);

InstallMethod( PreImagesElm, 
    "for a map between representations and an element in a representation.",
    [ IsPathAlgebraMatModuleMap, IsAlgebraModuleElement ], 0, 
    function( map, elem )

    local elt, i, m, fam, preimage;

    if elem in Range(map) then
        elt := ExtRepOfObj(elem)![1];
        m := Basis(Source(map))[1];
        preimage := [];
        for i in [1..Length(elt)] do 
            Add(preimage,SolutionMat(map!.maps[i],elt[i]));
        od;
        if ForAll(preimage,x-> x <> fail) then 
           preimage := PathModuleElem(FamilyObj(Zero(Source(map))![1]),preimage);        
           return Objectify( TypeObj( m ), [ preimage ] );
        else
	   return fail;
        fi;
    else
        Error("the element entered is not an element in the range of the map,");
    fi; 
end);


InstallMethod( RightModuleHomOverPathAlgebra,
    "for two representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule, IsList ], 0,
    function( M, N, linmaps)

  local A, K, mat_M, mat_N, map, Fam, dim_M, dim_N, Q, arrows, 
        vertices, a, i, origin, target, j;

  A := RightActingAlgebra(M); 
  if A <> RightActingAlgebra(N) then
     Error("the two modules are not over the same algebra, ");
  fi;
  dim_M := DimensionVector(M);
  dim_N := DimensionVector(N);
#
# Checking the number of matrices entered. 
#
  if Length( linmaps ) <> Length(dim_M) then 
     Error("the number of matrices entered is wrong,");
  fi;     
#
# Check the matrices of linmaps
#
  K := LeftActingDomain(A);
  for i in [1..Length(dim_M)] do
      if ( dim_M[i] = 0 ) then 
      	 if dim_N[i] = 0 then 
	        if linmaps[i] <> NullMat(1,1,K) then 
 	           Error("the dimension of matrix number ",i," is wrong,");
	        fi; 
         else
	        if linmaps[i] <> NullMat(1,dim_N[i],K) then
      	       Error("the dimension of matrix number ",i," is wrong,");
	        fi;
         fi;
      else
	     if dim_N[i] = 0 then
	        if linmaps[i] <> NullMat(dim_M[i],1,K) then 
      	       Error("the dimension of matrix number ",i," is wrong,");
            fi;
	     else 	        
            if DimensionsMat(linmaps[i])[1] <> dim_M[i] or DimensionsMat(linmaps[i])[2] <> dim_N[i] then
      	       Error("the dimension of matrix number ",i," is wrong,"); 
            fi;
	     fi;
      fi; 	 
  od;
# 
# Check commutativity relations with the matrices in M and N.
#
  mat_M := MatricesOfPathAlgebraMatModule(M);
  mat_N := MatricesOfPathAlgebraMatModule(N);
  Q := QuiverOfPathAlgebra(A);
  arrows   := ArrowsOfQuiver(Q);
  vertices := VerticesOfQuiver(Q);
  for a in arrows do
    i := Position(arrows,a); 
    origin := Position(vertices,SourceOfPath(a));
    target := Position(vertices,TargetOfPath(a)); 
    if mat_M[i]*linmaps[target] <> linmaps[origin]*mat_N[i] then 
        Error("entered map is not a module map between the representations entered, error for vertex number ",i," and arrow ",a,",");
    fi;
  od;
#
#
  map := Objectify( NewType( CollectionsFamily( GeneralMappingsFamily(
                                ElementsFamily( FamilyObj( M ) ),
                                ElementsFamily( FamilyObj( N ) ) ) ), 
#          IsPathAlgebraMatModuleMapRep ), rec( maps := linmaps ));
          IsPathAlgebraMatModuleMap and IsPathAlgebraMatModuleMapRep ), rec( maps := linmaps ));
  SetPathAlgebraOfMatModuleMap(map, A);
  SetSource(map, M);
  SetRange(map, N);  
  SetIsWholeFamily(map, true);

  return map;

end );

InstallMethod ( Zero, 
  "for a PathAlgebraMatModuleMap",
  true,
  [ IsPathAlgebraMatModuleMap ],
  0,
  function( f )
      local M, N, K, dim_M, dim_N, i, mats;

  M     := Source(f);
  N     := Range(f);
  K     := LeftActingDomain(Source(f));
  dim_M := DimensionVector(M);
  dim_N := DimensionVector(N);
  mats  := [];
  for i in [1..Length(dim_M)] do
     if dim_M[i] = 0 then
        if dim_N[i] = 0 then 
           Add(mats,NullMat(1,1,K));
        else
           Add(mats,NullMat(1,dim_N[i]));
        fi;
     else
        if dim_N[i] = 0 then 
           Add(mats,NullMat(dim_M[i],1,K));
        else
           Add(mats,NullMat(dim_M[i],dim_N[i]));
        fi;
     fi;
  od;
  
  return RightModuleHomOverPathAlgebra(Source(f),Range(f),mats);
end
);

InstallMethod ( ZeroMap, 
  " between two PathAlgebraMatModule's ",
  true,
  [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
  0,
  function( M, N )
      local A, K, dim_M, dim_N, i, mats;

  A := RightActingAlgebra(M);
  if ( A = RightActingAlgebra(N) ) then 
     K     := LeftActingDomain(M);
     dim_M := DimensionVector(M);
     dim_N := DimensionVector(N);
     mats  := [];
     for i in [1..Length(dim_M)] do
        if dim_M[i] = 0 then
           if dim_N[i] = 0 then 
              Add(mats,NullMat(1,1,K));
           else
              Add(mats,NullMat(1,dim_N[i]));
           fi;
        else
           if dim_N[i] = 0 then 
              Add(mats,NullMat(dim_M[i],1,K));
           else
              Add(mats,NullMat(dim_M[i],dim_N[i]));
           fi;
        fi;
     od;
  
     return RightModuleHomOverPathAlgebra(M,N,mats);
  else
     Error("the two modules entered are not modules of one and the same algebra, or they are not modules over the same (quotient of a) path algebra, ");
  fi;
end
);


InstallMethod ( \=, 
  "for a PathAlgebraMatModuleMap",
  true,
  [ IsPathAlgebraMatModuleMap, IsPathAlgebraMatModuleMap ],
  0,
  function( f, g )
      local a;

  if ( Source(f) = Source(g) ) and ( Range(f) = Range(g) ) 
      and ( f!.maps = g!.maps ) then
         return true;
  else
     return false;
  fi;
end
);

InstallMethod( SubRepInclusion,
  "for a path algebra module and list of its elements",
  true,
  [IsPathAlgebraMatModule, IsList], 0,
  function( M, gen )

  local A, q, K, num_vert, basis_M, vertices, arrows_as_path, arrows_of_quiver, 
  	newgen, g, v, submodspan, temp, new, V_list, m, V, basis_submod, submod_list, 
 	dim_vect_sub, i, cnt, j, dim_size, s, t, big_mat, a, mat, dom_a, im_a, pd, 
	pi, arrow, submodule, dim_vect_M, dim_size_M, inclusion, map;

    if not ForAll(gen, g -> g in M) then
       Error("entered elements are not in the module <M>, ");
    fi;

    A := RightActingAlgebra(M);
    q := QuiverOfPathAlgebra(A);
    K := LeftActingDomain(M);
    num_vert := Length(VerticesOfQuiver(q));
    basis_M  := Basis(M);
#
# vertices, vertices as elements of algebra
# arrows_as_path, arrows as elements of algebra
# arrows, as arrows in the quiver
#
if Length(gen) = 0 then 
   return ZeroMap(ZeroRepresentation(A),M);
else
    vertices := List(VerticesOfQuiver(q), x -> x*One(A));
    arrows_as_path := List(ArrowsOfQuiver(q), x -> x*One(A));
    arrows_of_quiver := GeneratorsOfQuiver(q){[1+num_vert..num_vert+Length(ArrowsOfQuiver(q))]};
#
# Ensuring uniform generators
#
    newgen := [];
    for g in gen do 
       for v in vertices do 
          if g^v <> Zero(M) then 
             Add(newgen,g^v);
          fi;
       od;
    od;
#
# Finding the linear span of the submodule generated by gen in M
#  
    submodspan := [];
    temp := [];
    new  := newgen;
    while new <> [] do
        for m in new do
            for a in arrows_as_path do
                if m^a <> Zero(M) then
                    Add(temp,m^a);            
                fi;
            od;
        od;
        Append(submodspan,new);
        new  := temp;
        temp := [];
    od;
    V_list := [];
    for m in submodspan do
        Add(V_list,Coefficients(basis_M,m));
    od;
    V := VectorSpace(K,V_list);
    basis_submod := CanonicalBasis(V);
#
# Converting elements in basis_submod to a list of elements in M 
#
    submod_list := [];
    for m in basis_submod do
        Add(submod_list,LinearCombination(basis_M,m));
    od;
#
# Finding the dimension vector of the submodule of M
#
    dim_vect_sub:=[];
    for i in [1.. num_vert] do
        cnt := 0;
        for j in [1..Length(submod_list)] do
            if submod_list[j]^vertices[i] <> Zero(M) then
                cnt:=cnt+1;   
            fi;
        od;
        Add(dim_vect_sub,cnt);    
    od;
    if Dimension(V) <> Sum(dim_vect_sub) then
       Error("Bug alert: Something is wrong in this code! Report this.\n");
    fi;
#    
# CanonicalBasis(V) gives the basis in "upper triangular form", so that 
# first comes the basis vectors for the vector space in vertex 1, in 
# vertex 2, ...., in vertex n. Next we find the intervals of the basis
# vector [[?..?],[?..?],...,[?..?]]. Is no basis vectors for a vertex, [] is entered.
#
    dim_size := [];
    s := 1;
    t := 0;
    for i in [1..Length(vertices)] do
        t := dim_vect_sub[i] + t;
	if t < s then 
	   Add(dim_size,[]);
        else
	   Add(dim_size,[s..t]);
        fi;
        s := t + 1;
    od;  
#
# Finding the submodule as a representation of the quiver
#  
    big_mat:=[];
    for a in arrows_as_path do
        mat := [];
        for v in vertices do
           if v*a <> Zero(A) then
              dom_a := v;
           fi;
        od;
        for v in vertices do
           if a*v <> Zero(A) then
              im_a := v;
           fi;
        od; 
        
        pd := Position(vertices,dom_a);
        pi := Position(vertices,im_a);
        
        arrow := arrows_of_quiver[Position(arrows_as_path,a)];
        if ( dim_vect_sub[pd] = 0 ) or ( dim_vect_sub[pi] = 0 ) then 
	       mat := [dim_vect_sub[pd],dim_vect_sub[pi]];
        else 
           for m in submod_list{dim_size[pd]} do
              Add(mat,Coefficients(basis_submod,Coefficients(basis_M,m^a)){dim_size[pi]}); 
           od;
        fi;
        Add(big_mat,[arrow,mat]);
    od;

    if IsPathAlgebra(A) then 
       submodule := RightModuleOverPathAlgebra(A,big_mat);
    else
       submodule := RightModuleOverQuotientOfPathAlgebra(A,big_mat); 
    fi;      

#
# Finding inclusion map of submodule into M
#
    mat := [];
    for i in [1..Length(basis_submod)] do
       Add(mat,basis_submod[i]);
    od;

    dim_vect_M := DimensionVector(M);
    dim_size_M := [];
    s := 1;
    t := 0;
    for i in [1..Length(vertices)] do
        t := dim_vect_M[i] + t;
        if t < s then 
           Add(dim_size_M,[]);
        else
           Add(dim_size_M,[s..t]);
        fi;
        s := t + 1;
    od;

    inclusion := [];
    for i in [1..Length(vertices)] do
       if dim_vect_sub[i] = 0 then
       	  if dim_vect_M[i] = 0 then 
       	     Add(inclusion,NullMat(1,1,K));
          else
             Add(inclusion,NullMat(1,dim_vect_M[i],K)); 
          fi;
       else 
          Add(inclusion,mat{dim_size[i]}{dim_size_M[i]});
       fi;
    od;
    return RightModuleHomOverPathAlgebra(submodule,M,inclusion);
fi;
end
);

InstallMethod( SubRep,
  "for a path algebra module and list of its elements",
  true,
  [IsPathAlgebraMatModule, IsList], 0,
  function( M, gen );

  return Source(SubRepInclusion(M,gen));
end
);


InstallMethod( RadicalOfRepInclusion,
  "for a path algebra module",
  true,
  [ IsPathAlgebraMatModule ], 0,
  function( M )

  local A, q, num_vert, arrows_as_path, basis_M, generators, a, b; 

    A := RightActingAlgebra(M);
    q := QuiverOfPathAlgebra(A);
    num_vert := Length(VerticesOfQuiver(q));
#
# Note arrows_as_path will change if A is a quotient of a path algebra !!!!
#
    arrows_as_path := List(ArrowsOfQuiver(q), x -> x*One(A));
#    arrows_as_path  := GeneratorsOfAlgebra(A){[1+num_vert..num_vert+Length(ArrowsOfQuiver(q))]};
    basis_M := Basis(M);
    generators := [];
    for a in arrows_as_path do
       for b in basis_M do 
	  if b^a <> Zero(M) then 
	     Add(generators,b^a);
          fi;
       od;
    od;
    
    return SubRepInclusion(M,generators);
end
);

InstallMethod( RadicalOfRep,
  "for a path algebra module",
  true,
  [ IsPathAlgebraMatModule ], 0,
  function( M )

  return Source(RadicalOfRepInclusion(M));
end
);

InstallMethod ( IsOneToOne, 
  "for a PathAlgebraMatModuleMap",
  true,
  [ IsPathAlgebraMatModuleMap ],
  0,
  function( f )
      local M, K, V_list, dim_K, dim_M, i, VS_list;

  M := Source(f);
  K := LeftActingDomain(M);
  V_list := [];
  dim_K  := 0;
  dim_M := DimensionVector(M);
#
# Computing the kernel of each f_i and making a vector space of Ker f_i with 
# basis given by the vectors supplied by NullspaceMat(f!.maps[i]).
#
  for i in [1..Length(dim_M)] do
      if dim_M[i] <> 0 then 
         dim_K := dim_K + Length(NullspaceMat(f!.maps[i]));
      fi;
  od;
  if dim_K = 0 then
     SetIsOneToOne(f,true);
     return true;
  else
     SetIsOneToOne(f,false);
     return false;
  fi;

end
);

InstallMethod ( KerInclusion, 
  "for a PathAlgebraMatModuleMap",
  true,
  [ IsPathAlgebraMatModuleMap ],
  0,
  function( f )
      local M, dim_M, V_list, V_dim, i, dim_K, VS_list, A, K, vertices, 
            arrows, mats, kermats, a, apos, j, matrix, k, Kerf, kermap;

  M := Source(f);
  K := LeftActingDomain(M);
  V_list := [];
  dim_K  := [];
  dim_M := DimensionVector(M);
  VS_list := [];
#
# Computing the kernel of each f_i and making a vector space of Ker f_i with 
# basis given by the vectors supplied by NullspaceMat(f!.maps[i]).
#
  for i in [1..Length(dim_M)] do
     if dim_M[i] <> 0 then
        Add(V_list,NullspaceMat(f!.maps[i]));
        Add(dim_K,Length(V_list[i]));
        if Length(V_list[i]) > 0 then 
           Add(VS_list,VectorSpace(K,V_list[i],"basis"));
        else
           Add(VS_list,TrivialSubmodule(VectorSpace(K,[[1]])));
        fi;
     else
        Add(V_list,[]);
        Add(dim_K,0);
        Add(VS_list,TrivialSubmodule(VectorSpace(K,[[1]])));
     fi;
  od;
  V_list := List(VS_list, V -> Basis(V));
  V_dim := Sum(dim_K);
  A := RightActingAlgebra(M);
if V_dim = 0 then 
   Print("The homomorphism is one-to-one.\n");
   SetIsOneToOne(f,true);
   return ZeroMap(ZeroRepresentation(A),Source(f)); 
else 
  vertices := VerticesOfQuiver(QuiverOfPathAlgebra(A));
  arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(A));
  mats := MatricesOfPathAlgebraMatModule(M);
  kermats := [];
#
# The matrices f_alpha of the representation M is used to compute
# the induced action on the basis of Ker f_i.
#
  for a in arrows do
     apos := Position(arrows,a);
     i := Position(vertices,SourceOfPath(a));
     j := Position(vertices,TargetOfPath(a));
     matrix := [];
     for k in [1..Length(V_list[i])] do
        Add(matrix,Coefficients(V_list[j],V_list[i][k]*mats[apos]));
     od;
     if dim_K[i] = 0 then 
     	matrix := [0,dim_K[j]];
     fi;
     if dim_K[j] = 0 then 
     	matrix := [dim_K[i],0];
     fi;
     Add(kermats,[a,matrix]);
  od;
#
# Extracting the basis vectors of each Ker f_i, as a subspace of M_i,
# which give the matrices of the inclusion of Ker f  into M.
#
  for i in [1..Length(dim_M)] do
     if Dimension(VS_list[i]) = 0 then
        V_list[i] := []; 
        for j in [1..dim_M[i]] do
	   Add(V_list[i],Zero(K));
	od;
        V_list[i] := [V_list[i]];
     else 
        V_list[i] := BasisVectors(V_list[i]);
     fi;
  od;
  if IsPathAlgebra(A) then 
     Kerf := RightModuleOverPathAlgebra(A,kermats);
  else 
     Kerf := RightModuleOverQuotientOfPathAlgebra(A,kermats);
  fi;
  kermap := RightModuleHomOverPathAlgebra(Kerf,M,V_list);
  SetKernelOfWhat(kermap,f);
  SetIsOneToOne(kermap,true);
  return kermap;
fi;
end
);

InstallMethod ( Ker, 
  "for a PathAlgebraMatModuleMap",
  true,
  [ IsPathAlgebraMatModuleMap ],
  0,
  function( f )

  return Source(KerInclusion(f));
end
);

InstallMethod ( ImProjectionInclusion, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleMap ],
   0,
   function( f )
      local M, N, image, images, A, K, num_vert, vertices, arrows, gen_list, B, i,fam, n, s, v,
      	    pos, dim_M, dim_N, V, W, projection, dim_image, basis_image, basis_M, basis_N, a, b, image_mat, 
            mat, mats, source_a, target_a, pos_a, bb, image_f, C, map, partmap,
            inclusion, image_inclusion, image_projection;

   M := Source(f);
   N := Range(f);
   A := RightActingAlgebra(M);
   K := LeftActingDomain(M); 
   num_vert := Length(VerticesOfQuiver(QuiverOfPathAlgebra(A)));
#
   vertices := List(VerticesOfQuiver(QuiverOfPathAlgebra(A)), x -> x*One(A));
#
# Finding a basis for the vector space in each vertex for the image
#
   image := ImagesSet(f,Source(f));
   gen_list := [];
   for i in [1..Length(vertices)] do
      Add(gen_list,[]);
   od;
   if Length(image) > 0 then 
      for n in image do
         for v in vertices do
            if n^v <> Zero(N) then
               pos := Position(vertices,v);
               Add(gen_list[pos],ExtRepOfObj(n)![1][pos]);
            fi;
         od;
      od;
      dim_N := DimensionVector(N);
      V := [];    
      W := [];
      basis_image := [];   
      basis_N := [];   
      for i in [1..Length(vertices)] do
         V[i] := K^dim_N[i];
         basis_N[i] := CanonicalBasis(V[i]);
         W[i] := Subspace(V[i],gen_list[i]);
         Add(basis_image,CanonicalBasis(W[i]));
      od;
#
# Finding the matrices of the representation image
#
      arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(A));
      mats := MatricesOfPathAlgebraMatModule(N);
      image_mat := [];
      for a in arrows do 
         mat := [];
         pos_a := Position(arrows,a);
         source_a := Position(vertices,SourceOfPath(a)*One(A));
         target_a := Position(vertices,TargetOfPath(a)*One(A));
         if ( Length(basis_image[source_a]) = 0 ) or ( Length(basis_image[target_a]) = 0 ) then 
            mat := [Length(basis_image[source_a]),Length(basis_image[target_a])];
         else 
            for b in basis_image[source_a] do
               Add(mat,Coefficients(basis_N[i],b*mats[pos_a]));
            od;
         fi;
#         Display(mat);
         Add(image_mat,[a,mat]);
      od;
      if IsPathAlgebra(A) then 
         image_f := RightModuleOverPathAlgebra(A,image_mat);
      else
         image_f := RightModuleOverQuotientOfPathAlgebra(A,image_mat);
      fi;
#
# Finding inclusion map from the image to Range(f)
#     
      inclusion := [];
      for i in [1..num_vert] do 
         mat := [];
         if Length(basis_image[i]) = 0 then 
            if dim_N[i] = 0 then 
               Add(inclusion,NullMat(1,1,K));
            else
               Add(inclusion,NullMat(1,dim_N[i],K));
            fi;
         else
            mat := [];
            for b in basis_image[i] do 
               Add(mat,b);
            od;
            Add(inclusion,mat);
         fi;
      od; 
#      Display(inclusion);
      image_inclusion := RightModuleHomOverPathAlgebra(image_f,Range(f),inclusion);
      SetImageOfWhat(image_inclusion,f);
      SetIsOneToOne(image_inclusion,true);
#
# Finding the projection for Source(f) to the image
#
      dim_M := DimensionVector(M);
      basis_M := [];
      for i in [1..num_vert] do
         Add(basis_M,CanonicalBasis(K^dim_M[i]));
      od;
      projection := [];
      for i in [1..num_vert] do
         mat := [];
         if Length(basis_image[i]) = 0 then
            if dim_M[i] = 0 then  
               mat := NullMat(1,1,K);
            else
               mat := NullMat(dim_M[i],1,K);
            fi;
            Add(projection,mat);
         else 
            for b in basis_M[i] do
               Add(mat,Coefficients(basis_image[i],b*f!.maps[i]));
            od;
            Add(projection,mat);
         fi;
      od;
      image_projection := RightModuleHomOverPathAlgebra(Source(f),image_f,projection);
      SetImageOfWhat(image_projection,f);
      SetIsOnto(image_projection,true);

      return [image_projection,image_inclusion];
   else
      Print("The homomorphism is the zero homomorphism.");
      return true;
   fi;
end
);

InstallMethod ( ImProjection, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleMap ],
   0,
   function( f );

   return ImProjectionInclusion(f)[1];
end
);

InstallMethod ( ImInclusion, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleMap ],
   0,
   function( f );

   return ImProjectionInclusion(f)[2];
end
);

InstallMethod ( Im, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleMap ],
   0,
   function( f );

   return Range(ImProjectionInclusion(f)[1]);
end
);

InstallMethod ( IsZeroMap, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleMap ],
   0,
   function( f );

   if Length(ImagesSet(f,Source(f))) = 0 then 
      return true;
   else
      return false;
   fi;
end
);

InstallMethod ( IsOnto, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleMap ],
   0,
   function( f )
      local image, dim_image, K, V;

   image := ImagesSet(f,Source(f));
   if Length(image) = 0 then 
      dim_image := 0;
   else
      image := List(image, x -> Flat(ExtRepOfObj(x)![1])); 
      K := LeftActingDomain(Source(f));
      V := K^Length(image[1]);
      dim_image :=  Dimension(Subspace(V,image));
   fi;
   if dim_image = Dimension(Range(f)) then 
      return true;
   else
      return false;
   fi;
end
);

InstallMethod ( IsIsom, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleMap ],
   0,
   function( f );

   return IsOneToOne(f) and IsOnto(f);
end
);

InstallMethod ( CokerProjection, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleMap ],
   0,
   function( f )
      local M, N, image, A, K, num_vert, vertices, arrows, basis_list,i, n, v,
      	    pos, dim_N, V, W, projection, coker, basis_coker, basis_N, a, b,
            mat, mats, source_a, target_a, pos_a, bb, cokermat, C, map, partmap,
            morph;

   M := Source(f);
   N := Range(f);
   A := RightActingAlgebra(M);
   K := LeftActingDomain(M); 
   num_vert := Length(VerticesOfQuiver(QuiverOfPathAlgebra(A)));
#
   vertices := List(VerticesOfQuiver(QuiverOfPathAlgebra(A)), x -> x*One(A));
#
# Finding a basis for the vector space in each vertex for the cokernel
#
   image := ImagesSet(f,Source(f));
   basis_list := [];
   for i in [1..Length(vertices)] do
      Add(basis_list,[]);
   od;
   for n in image do
      for v in vertices do
         if n^v <> Zero(N) then
            pos := Position(vertices,v);
            Add(basis_list[pos],ExtRepOfObj(n)![1][pos]);
         fi;
      od;
   od;
   dim_N := DimensionVector(N);
   V := [];    
   W := [];
   projection := [];
   basis_coker := [];   
   basis_N := [];   
   coker := [];
   for i in [1..Length(vertices)] do
      V[i] := K^dim_N[i];
      basis_N[i] := CanonicalBasis(V[i]);
      W[i] := Subspace(V[i],basis_list[i]);
      projection[i] := NaturalHomomorphismBySubspace(V[i],W[i]);
      Add(basis_coker,Basis(Range(projection[i])));
      coker[i] := Range(projection[i]);
   od;
#
# Finding the matrices of the representation of the cokernel 
#
   mats := MatricesOfPathAlgebraMatModule(N);
   arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(A));
   cokermat := [];
   for a in arrows do 
      mat := [];
      source_a := Position(vertices,SourceOfPath(a)*One(A));
      target_a := Position(vertices,TargetOfPath(a)*One(A));
      pos_a := Position(arrows,a);
      if ( Length(basis_coker[source_a]) = 0 ) or ( Length(basis_coker[target_a]) = 0 ) then
         mat := [Length(basis_coker[source_a]),Length(basis_coker[target_a])];
      else  
         for b in basis_coker[source_a] do
            bb := PreImagesRepresentative(projection[source_a],b);
            bb := bb*mats[pos_a]; # computing bb^a
            Add(mat,Coefficients(basis_coker[target_a],Image(projection[target_a],bb)));
         od;
      fi;
      Add(cokermat,[a,mat]);
   od;

   if IsPathAlgebra(A) then 
      C := RightModuleOverPathAlgebra(A,cokermat);
   else
      C := RightModuleOverQuotientOfPathAlgebra(A,cokermat);
   fi;
#
# Finding the map for Range(f) to the cokernel 
#
   map := [];
   for i in [1..Length(vertices)] do
      partmap := [];
      if dim_N[i] = 0 then
         partmap := NullMat(1,1,K);
      elif Length(basis_coker[i]) = 0 then 
         partmap := NullMat(dim_N[i],1,K);
      else 
         for b in basis_N[i] do 
            Add(partmap,Coefficients(basis_coker[i],Image(projection[i],b)));
         od;
      fi;
      Add(map,partmap);
   od;

   morph := RightModuleHomOverPathAlgebra(Range(f),C,map);
   SetCoKernelOfWhat(morph,f);
   SetIsOnto(morph,true);

   return morph;
end
);

InstallMethod ( Coker, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleMap ],
   0,
   function( f )

   return Range(CokerProjection(f));
end
);

InstallMethod( TopOfRepProjection,
  "for a path algebra module",
  true,
  [ IsPathAlgebraMatModule ], 0,
  function( M )

  local map; 

  map := RadicalOfRepInclusion(M);

  return CokerProjection(map);
end
);

InstallMethod( TopOfRep,
  "for a path algebra module",
  true,
  [ IsPathAlgebraMatModule ], 0,
  function( M );

  return Range(TopOfRepProjection(M));
end
);

InstallMethod( GeneratorsOfRep, 
  "for a path algebra module",
  true,
  [ IsPathAlgebraMatModule ], 0,
  function( M )

  local generators, map, B, i; 

  generators := [];
  if Dimension(M) <> 0 then 
     map := TopOfRepProjection(M);
     B := Basis(Range(map));
     generators := [];
     for i in [1..Length(B)] do
        Add(generators,PreImagesElm(map,B[i]));
     od;
  fi;
  
  return generators;
end
);


InstallMethod( \+,
  "for two PathAlgebraMatModuleMap's",
  true,
#  IsIdenticalObj,
  [ IsPathAlgebraMatModuleMap,
    IsPathAlgebraMatModuleMap ],
  0,
  function( f, g )
     local i, num_vert, x, Fam;

     if ( Source(f) = Source(g) ) and ( Range(f) = Range(g) ) then 
        x := [];
     	num_vert := Length(f!.maps);
     	for i in [1..num_vert] do
     	   x[i] := f!.maps[i] + g!.maps[i];
     	od;
	return RightModuleHomOverPathAlgebra(Source(f),Range(f),x);
     else
	Error("the two arguments entered do not live in the same homomorphism set, ");
     fi;
  end
);

InstallMethod( \*,
  "for two PathAlgebraMatModuleMap's",
  true,
  [ IsPathAlgebraMatModuleMap,
    IsPathAlgebraMatModuleMap ],
  0,
  function( f, g )
     local i, num_vert, x;

     if Range(f) = Source(g) then 
        x := [];
     	num_vert := Length(f!.maps);
     	for i in [1..num_vert] do
     	   x[i] := f!.maps[i]*g!.maps[i];
     	od;
	return RightModuleHomOverPathAlgebra(Source(f),Range(g),x);
     else
        Error("codomain of the first argument is not equal to the domain of the second argument, ");
     fi;
  end
);

InstallOtherMethod( \*,
  "for two PathAlgebraMatModuleMap's",
  true,
  [ IsScalar,
    IsPathAlgebraMatModuleMap ],
  0,
  function( a, g )
     local K, i, num_vert, x;

     K := LeftActingDomain(Source(g));
     if a in K then 
        x := [];
     	num_vert := Length(g!.maps);
     	for i in [1..num_vert] do
     	   x[i] := a*g!.maps[i];
     	od;
	return RightModuleHomOverPathAlgebra(Source(g),Range(g),x);
     else
	Error("the scalar is not in the same field as the algbra is over,");
     fi;
  end
);

InstallOtherMethod( \*,
  "for two PathAlgebraMatModuleMap's",
  true,
  [ IsPathAlgebraMatModuleMap, IsScalar ],
  0,
  function( f, a )
     local K, i, num_vert, x;

     K := LeftActingDomain(Source(f));
     if a in K then 
        x := [];
     	num_vert := Length(f!.maps);
     	for i in [1..num_vert] do
     	   x[i] := f!.maps[i]*a;
     	od;
	    return RightModuleHomOverPathAlgebra(Source(f),Range(f),x);
     else
	    Error("the scalar is not in the same field as the algbra is over,");
     fi;
  end
);

InstallMethod( HomOverPathAlgebra,
    "for two representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
  function( M, N )

  local R, F, zero, Mfam, Nfam, quiver, Mdims, Ndims, Mmaps, Nmaps,
        l, arrows, i, arrowPos, sourcePos, targetPos,
        a, tMat, numColumns, numRows, r, c, row, col, startCol,
        startRow, equations, j, k, MBasis, NBasis, mat, x, y,
        Msum, Nsum, ns, basis, maps, b, matrise;

  R := RightActingAlgebra(M); 
  if R <> RightActingAlgebra(N) then
     Error("the two representations are not for the same quiver.");
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
  Display(equations);
  
  # Create blocks
  startRow := 1;
  for a in arrows do
    arrowPos := a!.gen_pos;
    sourcePos := SourceOfPath(a)!.gen_pos;
    targetPos := TargetOfPath(a)!.gen_pos;

    # We always have a transposed version of the maps for $N$
    tMat := TransposedMat(Nmaps[arrowPos]);
    Display(tMat);
    r := DimensionsMat(tMat)[1];
    c := DimensionsMat(tMat)[2];
    for i in [1..Mdims[sourcePos]] do
      row := startRow+(i-1)*r;
      col := startCol[sourcePos]+(i-1)*c;
      Print("Intervals: ",[row..row+r-1],[col..col+c-1],"\n");
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

  Display(equations);
  
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

  maps := []; 
  for i in [1..Length(basis)] do
    matrise := [];
    r := 1;
    c := 1;
    for j in [1..l] do
        matrise[j] := basis[i]{[r..r+Mdims[j]-1]}{[c..c+Ndims[j]-1]};
        r := r + Mdims[j];
        c := c + Ndims[j];
    od;
    Add(maps,matrise);    
  od;

  for i in [1..Length(basis)] do
     maps[i] := Objectify( NewType( CollectionsFamily( GeneralMappingsFamily(
                                ElementsFamily( FamilyObj( M ) ),
                                ElementsFamily( FamilyObj( N ) ) ) ), 
                     IsPathAlgebraMatModuleMap and IsPathAlgebraMatModuleMapRep and IsAttributeStoringRep ), rec( maps := maps[i] ));
     SetPathAlgebraOfMatModuleMap(maps[i], R);
     SetSource(maps[i], M);
     SetRange(maps[i], N);
     SetIsWholeFamily(maps[i],true);
  od;

#
# LeftAlgebraModuleByGenerators(Rationals,\*,maps);
# gen:=GeneratorsOfAlgebraModule(last);
# ExtRepOfObj(gen[1])!.maps;
# CategoriesOfObject(ExtRepOfObj(gen[1]));
#
  return maps;
end
);

InstallMethod( EndOverPathAlgebra,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule ], 0,
  function( M )

  local EndM, R, F, dim_M, alglist, i, j, r, maps, A; 

  EndM := HomOverPathAlgebra(M,M);
  R := RightActingAlgebra(M); 
  F := LeftActingDomain(R);
  dim_M := DimensionVector(M);
  alglist := [];
  for i in [1..Length(dim_M)] do 
     Add(alglist, MatrixAlgebra(F,dim_M[i]));
  od;
  maps := [];
  for i in [1..Length(EndM)] do
     maps[i] := NullMat(Dimension(M),Dimension(M),F);
     r := 1; 
     for j in [1..Length(dim_M)] do 
        maps[i]{[r..r+dim_M[j]-1]}{[r..r+dim_M[j]-1]} := EndM[i]!.maps[j];
        r := r + dim_M[j];
     od; 
  od;
  A := DirectSumOfAlgebras(alglist); 

  return SubalgebraWithOne(A,maps,"basis"); 
end
);

InstallMethod( RightFacApproximation,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
  function( M, N )

  local homMN, i, generators; 

  homMN := HomOverPathAlgebra(M,N); 
  generators := [];
  for i in [1..Length(homMN)] do
     Add(generators,ImagesSet(homMN[i],Source(homMN[i])));
  od;

  return SubRepInclusion(N,generators);
end
);

InstallMethod( NumberOfNonIsoDirSummands,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule ], 0,
  function( M )

  local EndM, K, J, gens, I, A, top, AA, B, n,
  	i, j, genA, V, W, d;

  EndM := EndOverPathAlgebra(M);
  K    := LeftActingDomain(M);
  J    := RadicalOfAlgebra(EndM);
  gens := GeneratorsOfAlgebra(J);
  I    := Ideal(EndM,gens); 
  A    := EndM/I;
  top  := CentralIdempotentsOfAlgebra(A);
     
  return [Length(top),List(DirectSumDecomposition(EndM/I),Dimension)];
end
);

InstallMethod ( DualOfPathAlgebraMatModuleMap,
    "for a map between representations of a quiver",
    [ IsPathAlgebraMatModuleMap ], 0,
    function( f )

    local mats, M, N;
   
    mats := f!.maps;
    mats := List(mats, x -> TransposedMat(x));
    M := DualOfPathAlgebraMatModule(Source(f));
    N := DualOfPathAlgebraMatModule(Range(f));

    return RightModuleHomOverPathAlgebra(N,M,mats);

end
);

InstallMethod ( SocleOfPathAlgebraMatModuleInclusion,
    "for a map between representations of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M );

    return DualOfPathAlgebraMatModuleMap(TopOfRepProjection(DualOfPathAlgebraMatModule(M)));    
end
);

InstallMethod ( SocleOfPathAlgebraMatModule,
    "for a map between representations of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M );

    return Source(DualOfPathAlgebraMatModuleMap(TopOfRepProjection(DualOfPathAlgebraMatModule(M))));    
end
);

