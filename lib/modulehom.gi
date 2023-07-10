# GAP Implementation
# $Id: homomorphisms.gi,v 1.40 2012/09/28 12:57:10 sunnyquiver Exp $

#############################################################################
##
#M  ImageElm( <map>, <elm> )  . . . for PathAlgebraMatModuleMap and element
##
InstallMethod( ImageElm, 
    "for a map between representations and an element in a representation.",
    [ IsPathAlgebraMatModuleHomomorphism, IsAlgebraModuleElement ], 0, 
    function( map, elem )

    local elt, n, fam, image, temp, zero, new_image;

    if elem in Source(map) then
        if Dimension(Range(map)) = 0 then
            return Zero(Range(map));
        fi;
  	elt := ExtRepOfObj(elem);
        n := Zero(Range(map));
        image := List([1..Length(elt![1])], x -> elt![1][x]*map!.maps[x]);

        image := PathModuleElem(FamilyObj(Zero(Range(map))![1]),image);        
        return Objectify( TypeObj( n ), [ image ] );
    else
        Error("the element entered is not an element in the source of the map,");
    fi; 
end);

#############################################################################
##
#M  ImagesSet( <map>, <elms> ) . . for a PathAlgebraMatModuleMap and finite collection
##
InstallMethod( ImagesSet, 
    "for a map between representations and a set of elements in a representation.",
    [ IsPathAlgebraMatModuleHomomorphism, IsCollection ], 0, 
    function( map, elms )

    local elt, B, images;

    images := [];
    if IsList(elms) then
        if IsFinite(elms) then
            B := elms;
        fi;
    else
        if IsPathAlgebraMatModule(elms) then 
            B := BasisVectors(CanonicalBasis(elms));
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
#############################################################################
##
#M  PreImagesRepresentative( <map>, <elms> ) . . for a 
##                              PathAlgebraMatModuleMap and finite collection
##
InstallMethod( PreImagesRepresentative, 
    "for a map between representations and an element in a representation.",
    [ IsPathAlgebraMatModuleHomomorphism, IsAlgebraModuleElement ], 0, 
    function( map, elem )

    local elt, m, fam, preimage;

    if elem in Range(map) then
        elt := ExtRepOfObj(elem)![1];
        m := Basis(Source(map))[1];
        preimage := List([1..Length(elt)], x -> SolutionMat(map!.maps[x],elt[x]));
        if ForAll(preimage,x -> x <> fail) then 
           preimage := PathModuleElem(FamilyObj(Zero(Source(map))![1]),preimage);        
           return Objectify( TypeObj( m ), [ preimage ] );
        else
	   return fail;
        fi;
    else
        Error("the element entered is not an element in the range of the map,");
    fi; 
end);

#######################################################################
##
#O  RightModuleHomOverAlgebra( <M>, <N>, <linmaps> )
##
##  This function constructs a homomorphism  f  from the module  <M>  to
##  the module  <N>  from the linear maps given in  <linmaps>. The 
##  function checks if  <M>  and  <N>  are modules over the same algebra 
##  and checks if the linear maps  <linmaps>  defines a homomorphism from  
##  <M>  to  <N>. The source and the range of  f  can be recovered from 
##  f  via the function  Source(f)  and  Range(f). The linear maps can 
##  be recovered via  f!.maps  or  
##  MatricesOfPathAlgebraMatModuleHomomorphism(f).
##
InstallMethod( RightModuleHomOverAlgebra,
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
                     Error("the dimension of matrix number ",i," is wrong (A),");
                fi; 
            else
                if linmaps[i] <> NullMat(1,dim_N[i],K) then
                    Error("the dimension of matrix number ",i," is wrong (B),");
                fi;
            fi;
        else
            if dim_N[i] = 0 then
                if linmaps[i] <> NullMat(dim_M[i],1,K) then 
                    Error("the dimension of matrix number ",i," is wrong (C),");
                fi;
            else 	        
                if DimensionsMat(linmaps[i])[1] <> dim_M[i] or DimensionsMat(linmaps[i])[2] <> dim_N[i] then
                    Error("the dimension of matrix number ",i," is wrong (D),"); 
                fi;
            fi;
        fi; 	 
    od;
# 
# Check commutativity relations with the matrices in M and N.
#
    mat_M := MatricesOfPathAlgebraModule(M);
    mat_N := MatricesOfPathAlgebraModule(N);
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
                   IsPathAlgebraMatModuleHomomorphism and IsPathAlgebraMatModuleHomomorphismRep ), rec( maps := linmaps ));
    SetPathAlgebraOfMatModuleMap(map, A);
    SetSource(map, M);
    SetRange(map, N);  
    SetIsWholeFamily(map, true);
    
    return map;
end 
);

#######################################################################
##
#O  MatricesOfPathAlgebraMatModuleHomomorphism( <f> )
##
##  This function returns the matrices defining the homomorphism
##  <f>. 
##
InstallMethod( MatricesOfPathAlgebraMatModuleHomomorphism, 
    "for a PathAlgebraMatModuleHomomorphism",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )
    return f!.maps;  
end
);

#######################################################################
##
#M  ViewObj( <f> )
##
##  This function defines how View prints a 
##  PathAlgebraMatModuleHomomorphism  <f>.
##
InstallMethod( ViewObj, 
    "for a PathAlgebraMatModuleHomomorphism",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ], NICE_FLAGS+1,
    function ( f );

    Print("<");
    View(Source(f));
    Print(" ---> ");
    View(Range(f));
    Print(">");
end
); 

#######################################################################
##
#M  PrintObj( <f> )
##
##  This function defines how Print prints a 
##  PathAlgebraMatModuleHomomorphism  <f>.
##
InstallMethod( PrintObj, 
    "for a PathAlgebraMatModuleHomomorphism",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ], NICE_FLAGS + 1,
    function ( f );

    Print("<",Source(f)," ---> ",Range(f),">");
end
); 

#######################################################################
##
#M  Display( <f> )
##
##  This function defines how Display prints a 
##  PathAlgebraMatModuleHomomorphism  <f>.
##
InstallMethod( Display, 
    "for a PathAlgebraMatModuleHomomorphism",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ], NICE_FLAGS + 1,
    function ( f )
    
    local i;

    Print(f);
    Print("\nwith ");
    for i in [1..Length(DimensionVector(Source(f)))] do
        Print("linear map for vertex number ",i,":\n");
        PrintArray(f!.maps[i]);
    od;
end
); 

#######################################################################
##
#M  Zero( <f> )
##
##  This function returns the zero mapping with the same source and 
##  range as  <f>. 
##
InstallMethod ( Zero, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )
    
    local M, N, K, dim_M, dim_N, i, mats;

    M := Source(f);
    N := Range(f);
    K := LeftActingDomain(Source(f));
    dim_M := DimensionVector(M);
    dim_N := DimensionVector(N);
    mats := [];
    for i in [1..Length(dim_M)] do
        if dim_M[i] = 0 then
            if dim_N[i] = 0 then 
                Add(mats,NullMat(1,1,K));
            else
                Add(mats,NullMat(1,dim_N[i],K));
            fi;
        else
            if dim_N[i] = 0 then 
                Add(mats,NullMat(dim_M[i],1,K));
            else
                Add(mats,NullMat(dim_M[i],dim_N[i],K));
            fi;
        fi;
    od;
  
    return RightModuleHomOverAlgebra(Source(f),Range(f),mats);
end
);

#######################################################################
##
#M  ZeroMapping( <M>, <N> )
##
##  This function returns the zero mapping from the module  <M>  to the 
##  module  <N>. 
##
InstallMethod ( ZeroMapping, 
    " between two PathAlgebraMatModule's ",
    true,
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
    0,
    function( M, N )
    
    local A, K, dim_M, dim_N, i, mats;

    A := RightActingAlgebra(M);
    if ( A = RightActingAlgebra(N) ) then 
        K := LeftActingDomain(M);
        dim_M := DimensionVector(M);
        dim_N := DimensionVector(N);
        mats  := [];
        for i in [1..Length(dim_M)] do
            if dim_M[i] = 0 then
                if dim_N[i] = 0 then 
                    Add(mats,NullMat(1,1,K));
                else
                    Add(mats,NullMat(1,dim_N[i],K));
                fi;
            else
                if dim_N[i] = 0 then 
                    Add(mats,NullMat(dim_M[i],1,K));
                else
                    Add(mats,NullMat(dim_M[i],dim_N[i],K));
                fi;
            fi;
        od;
  
        return RightModuleHomOverAlgebra(M,N,mats);
    else
        Error("the two modules entered are not modules of one and the same algebra, or they are not modules over the same (quotient of a) path algebra, ");
    fi;
end
);

#######################################################################
##
#M  IdentityMapping( <M> )
##
##  This function returns the identity homomorphism from  <M>  to  <M>.
##
InstallMethod ( IdentityMapping, 
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule ],
    0,
    function( M )
    
    local K, dim_M, i, mats;
#
# Representing the identity map from M to M with identity matrices 
# of the right size, including if dim_M[i] = 0, then the identity 
# is represented by a one-by-one identity matrix.
#
    K := LeftActingDomain(M);
    dim_M := DimensionVector(M);
    mats := [];
    for i in [1..Length(dim_M)] do
        if dim_M[i] = 0 then
            Add(mats,NullMat(1,1,K));
        else
            Add(mats,IdentityMat(dim_M[i],K));
        fi;
    od;
    
    return RightModuleHomOverAlgebra(M,M,mats);
end
);

#######################################################################
##
#M  \=( <f>, <g> )
##
##  This function returns true if the homomorphisms  <f>  and  <g>  have
##  the same source, the same range and the matrices defining the 
##  homomorphisms are identitical.
##
InstallMethod ( \=, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ],
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

#######################################################################
##
#O  SubRepresentationInclusion( <M>, <gen> )
##
##  This function returns the inclusion from the submodule of  <M>  
##  generated by the elements  <gen>  to the module  <M>. The function
##  checks if all the elements on the list  <gen>  are elements of  the
##  module  <M>. 
##
InstallMethod( SubRepresentationInclusion,
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
    basis_M := Basis(M);
#
# vertices, vertices as elements of algebra
# arrows_as_path, arrows as elements of algebra
# arrows, as arrows in the quiver
#
    if Length(gen) = 0 then 
        return ZeroMapping(ZeroModule(A),M);
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
        new := newgen;
        while new <> [] do
            for m in new do
                for a in arrows_as_path do
                    if m^a <> Zero(M) then
                        Add(temp,m^a);            
                    fi;
                od;
            od;
            Append(submodspan,new);
            new := temp;
            temp := [];
        od;
        V_list := List(submodspan, x -> Coefficients(basis_M,x));
        V := VectorSpace(K,V_list);
        basis_submod := CanonicalBasis(V);
#
# Converting elements in basis_submod to a list of elements in M 
#
        submod_list := List(basis_submod, x -> LinearCombination(basis_M,x));
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
        big_mat := [];
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
            if dim_vect_sub[pd] <> 0 and dim_vect_sub[pi] <> 0 then 
                for m in submod_list{dim_size[pd]} do
                    Add(mat,Coefficients(basis_submod,Coefficients(basis_M,m^a)){dim_size[pi]}); 
                od;
                Add(big_mat,[String(arrow),mat]);
            fi;
        od;
        
        if IsPathAlgebra(A) then 
            submodule := RightModuleOverPathAlgebra(A, dim_vect_sub, big_mat);
        else
            submodule := RightModuleOverPathAlgebra(A, dim_vect_sub, big_mat); 
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
        return RightModuleHomOverAlgebra(submodule,M,inclusion);
    fi;
end
);

#######################################################################
##
#O  SubRepresentation( <M>, <gen> )
##
##  This function returns a module isomorphic to the submodule of  <M>  
##  generated by the elements  <gen>  to the module  <M>. The function
##  checks if all the elements on the list  <gen>  are elements of  the
##  module  <M>. 
##
InstallMethod( SubRepresentation,
    "for a path algebra module and list of its elements",
    true,
    [ IsPathAlgebraMatModule, IsList], 0,
    function( M, gen );

    return Source(SubRepresentationInclusion(M,gen));
end
);

#######################################################################
##
#O  RadicalOfModuleInclusion( <M> )
##
##  This function returns the inclusion from the radical of  <M>  
##  to the module  <M>. 
##
InstallMethod( RadicalOfModuleInclusion,
    "for a path algebra module",
    true,
    [ IsPathAlgebraMatModule ], 0,
    function( M )

    local A, q, num_vert, arrows_as_path, basis_M, generators, a, b, run_time; 
    
    #
    # Checking if the algebra  <A>  is finite dimensional.  If it is and is 
    # a path algebra, then the radical is the ideal generated by the arrows.  
    # If it is not a path algebra but an admissible quotient of a path algebra,
    # the radical is the ideal generated by the arrows, and this function 
    # applies. 
    #
    A := RightActingAlgebra(M);
    q := QuiverOfPathAlgebra(A);
    if not IsFiniteDimensional(A) then
        TryNextMethod();
    fi; 
    if not IsPathAlgebra(A) and not IsAdmissibleQuotientOfPathAlgebra(A) then
        TryNextMethod();
    fi;

    num_vert := Length(VerticesOfQuiver(q));
#
# Note arrows_as_path will change if A is a quotient of a path algebra !!!!
#
    arrows_as_path := List(ArrowsOfQuiver(q), x -> x*One(A));
    basis_M := Basis(M);
    generators := [];
    for a in arrows_as_path do
        for b in basis_M do 
            if b^a <> Zero(M) then 
                Add(generators,b^a);
            fi;
        od;
    od;
    
    generators := Unique(generators);
    return SubRepresentationInclusion(M,generators);
end
);

#######################################################################
##
#O  RadicalOfModule( <M> )
##
##  This function returns a module isomorphic to the radical of  <M>. 
##
InstallMethod( RadicalOfModule,
    "for a path algebra module",
    true,
    [ IsPathAlgebraMatModule ], 0,
    function( M )

    return Source(RadicalOfModuleInclusion(M));
end
);

#######################################################################
##
#M  IsInjective( <f> )
##
##  This function returns a module isomorphic to the radical of  <M>. 
##
InstallOtherMethod ( IsInjective, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )
    
    local M, K, V_list, dim_K, dim_M, i, VS_list;

    M := Source(f);
    K := LeftActingDomain(M);
    V_list := [];
    dim_K := 0;
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
        SetIsInjective(f,true);
        return true;
    else
        SetIsInjective(f,false);
        return false;
    fi;
end
);

#######################################################################
##
#O  KernelInclusion( <f> )
##
##  This function returns the inclusion from a module isomorphic to the 
##  kernel of  <f>  to the module  <M>. 
##
InstallMethod ( KernelInclusion, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
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
                Add(VS_list,TrivialSubmodule(VectorSpace(K,[[One(K)]])));
            fi;
        else
            Add(V_list,[]);
            Add(dim_K,0);
            Add(VS_list,TrivialSubmodule(VectorSpace(K,[[One(K)]])));
        fi;
    od;
    V_list := List(VS_list, V -> Basis(V));
    V_dim := Sum(dim_K);
    A := RightActingAlgebra(M);
    if V_dim = 0 then 
        SetIsInjective(f,true);
        kermap := ZeroMapping(ZeroModule(A),Source(f)); 
    else 
        vertices := VerticesOfQuiver(QuiverOfPathAlgebra(A));
        arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(A));
        mats := MatricesOfPathAlgebraModule(M);
        kermats := [];
#
# The matrices f_alpha of the representation M is used to compute
# the induced action on the basis of Ker f_i.
#
        for a in arrows do
            i := Position(vertices,SourceOfPath(a));
            j := Position(vertices,TargetOfPath(a));
            if dim_K[i] <> 0 and dim_K[j] <> 0 then
                apos := Position(arrows,a);
                matrix := [];
                for k in [1..Length(V_list[i])] do
                    Add(matrix,Coefficients(V_list[j],V_list[i][k]*mats[apos]));
                od;
                Add(kermats,[String(a),matrix]);
            fi;
        od;
#
# Extracting the basis vectors of each Ker f_i, as a subspace of M_i,
# which give the matrices of the inclusion of Ker f  into M.
#
        for i in [1..Length(dim_M)] do
            if Dimension(VS_list[i]) = 0 then
                V_list[i] := []; 
                if dim_M[i] = 0 then 
                    Add(V_list[i],Zero(K));
                else 
                    for j in [1..dim_M[i]] do
                        Add(V_list[i],Zero(K));
                    od;
                fi;
                V_list[i] := [V_list[i]];
            else 
                V_list[i] := BasisVectors(V_list[i]);
            fi;
        od;
        if IsPathAlgebra(A) then 
            Kerf := RightModuleOverPathAlgebra(A, dim_K, kermats);
        else 
            Kerf := RightModuleOverPathAlgebra(A, dim_K, kermats);
        fi;
        kermap := RightModuleHomOverAlgebra(Kerf,M,V_list);
    fi;  
    SetKernelOfWhat(kermap,f);
    SetIsInjective(kermap,true);
    return kermap;
end
);

#######################################################################
##
#M  Kernel( <f> )
##
##  This function returns a module isomorphic to the kernel of  <f>.
##
InstallOtherMethod ( Kernel, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    return Source(KernelInclusion(f));
end
);

#######################################################################
##
#O  ImageProjectionInclusion( <f> )
##
##  This function returns a projection from the source of  <f>  to a 
##  module isomorphic to the image of  <f>  and an inclusion from a 
##  module isomorphic to the image of  <f>  to the module  <M>. 
##
InstallMethod ( ImageProjectionInclusion, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )
    
    local M, N, image, images, A, K, num_vert, vertices, arrows, gen_list, B, i,fam, n, s, v,
          pos, dim_M, dim_N, V, W, projection, dim_image, basis_image, basis_M, basis_N, a, b, image_mat, 
          mat, mats, source_a, target_a, pos_a, bb, image_f, C, map, partmap,
          inclusion, image_inclusion, image_projection, dimvecimage_f;

    M := Source( f );
    N := Range( f );
    A := RightActingAlgebra( M );
    K := LeftActingDomain( M ); 
    num_vert := Length( VerticesOfQuiver( QuiverOfPathAlgebra( A ) ) );
#
    vertices := List( VerticesOfQuiver( QuiverOfPathAlgebra( A ) ), x -> x * One( A ) );
#
# Finding a basis for the vector space in each vertex for the image
#
    image := ImagesSet( f, Source( f ) );
    gen_list := [ ];
    for i in [ 1..Length( vertices ) ] do
        Add( gen_list,[ ] );
    od;
    if Length( image ) > 0 then 
        for n in image do
            for v in vertices do
                if n^v <> Zero( N ) then
                    pos := Position( vertices, v );
                    Add( gen_list[ pos ], ExtRepOfObj( n )![ 1 ][ pos ] );
                fi;
            od;
        od;
        dim_N := DimensionVector( N );
        V := [];    
        W := [];
        basis_image := [];   
        basis_N := [];   
        for i in [ 1..Length( vertices ) ] do
            V[ i ] := K^dim_N[ i ];
            basis_N[ i ] := CanonicalBasis( V[ i ] );
            W[ i ] := Subspace( V[ i ], gen_list[ i ] );
            Add( basis_image, CanonicalBasis( W[ i ] ) );
        od;
#
# Finding the matrices of the representation image
#
        arrows := ArrowsOfQuiver( QuiverOfPathAlgebra( A ) );
        mats := MatricesOfPathAlgebraModule( N );
        image_mat := [ ];
        for a in arrows do 
            mat := [ ];
            pos_a := Position( arrows, a );
            source_a := Position( vertices, SourceOfPath( a ) * One( A ) );
            target_a := Position( vertices, TargetOfPath( a ) * One( A ) );
            if ( Length( basis_image[ source_a ]) <> 0 ) and ( Length( basis_image[ target_a ] ) <> 0 ) then 
	        for b in basis_image[ source_a ] do
                    Add( mat, Coefficients( basis_image[ target_a ], b * mats[ pos_a ] ) );
                od;
                Add( image_mat, [ String( a ), mat ] );
            fi;
        od;

	dimvecimage_f := List( W, Dimension );
	image_f := RightModuleOverPathAlgebra( A, dimvecimage_f, image_mat );
#
# Finding inclusion map from the image to Range(f)
#     
        inclusion := [ ];
        for i in [ 1..num_vert ] do 
            mat := [ ];
            if Length( basis_image[ i ] ) = 0 then 
                if dim_N[ i ] = 0 then 
                    Add( inclusion, NullMat( 1, 1, K ) );
                else
                    Add( inclusion, NullMat( 1, dim_N[ i ], K ) );
                fi;
            else
                mat := [ ];
                for b in basis_image[ i ] do 
                    Add( mat, b );
                od;
                Add( inclusion, mat );
            fi;
        od; 
        image_inclusion := RightModuleHomOverAlgebra( image_f, Range( f ), inclusion );
        SetImageOfWhat( image_inclusion, f );
        SetIsInjective( image_inclusion, true );
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
        image_projection := RightModuleHomOverAlgebra(Source(f),image_f,projection);
        SetImageOfWhat(image_projection,f);
        SetIsSurjective(image_projection,true);
        
        return [image_projection,image_inclusion];
    else
        return [ZeroMapping(M,ZeroModule(A)),ZeroMapping(ZeroModule(A),N)];
    fi;
end
);

#######################################################################
##
#O  ImageProjection( <f> )
##
##  This function returns a projection from the source of  <f>  to a 
##  module isomorphic to the image of  <f>. 
##
InstallMethod ( ImageProjection, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f );

    return ImageProjectionInclusion(f)[1];
end
);

#######################################################################
##
#O  ImageInclusion( <f> )
##
##  This function returns an inclusion from a module isomorphic to the 
##  image of  <f>  to the module  <M>. 
##
InstallMethod ( ImageInclusion, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f );

    return ImageProjectionInclusion(f)[2];
end
);

#######################################################################
##
#M  ImagesSource( <f> )
##
##  This function returns a module isomorphic to the image of  <f>. 
##  TODO: Should this delegate to ImagesSet as for the general GAP
##  command?
##
InstallMethod ( ImagesSource, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f );

    return Range(ImageProjectionInclusion(f)[1]);
end
);

#######################################################################
##
#M  IsZero( <f> )
##
##  This function returns true if all the matrices of the homomorphism
##  <f>  are identically zero. 
##
InstallMethod ( IsZero, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f );

    return ForAll(f!.maps, IsZero);
end
);

#######################################################################
##
#M  IsSurjective( <f> )
##
##  This function returns true if the homomorphism  <f>  is surjective. 
##
InstallOtherMethod ( IsSurjective, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
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
        SetIsSurjective(f,true);
        return true;
    else
        SetIsSurjective(f,false);
        return false;
    fi;
end
);

#######################################################################
##
#M  IsIsomorphism( <f> )
##
##  This function returns true if the homomorphism  <f>  is an 
##  isomorphism. 
##
InstallMethod ( IsIsomorphism, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f );

    return IsInjective(f) and IsSurjective(f);
end
);

#######################################################################
##
#O  CoKernelProjection( <f> )
##
##  This function returns a projection from the range of  <f>  to a 
##  module isomorphic to the cokernel of  <f>. 
##
InstallMethod ( CoKernelProjection, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )
    
    local M, N, image, A, K, num_vert, vertices, arrows, basis_list,i, n, v,
          pos, dim_N, V, W, projection, coker, basis_coker, basis_N, a, b,
          mat, mats, source_a, target_a, pos_a, bb, cokermat, C, map, partmap,
          morph, dimvec_coker;

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
    mats := MatricesOfPathAlgebraModule(N);
    arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(A));
    cokermat := [];
    for a in arrows do
        source_a := Position(vertices,SourceOfPath(a)*One(A));
        target_a := Position(vertices,TargetOfPath(a)*One(A));
        if Length(basis_coker[source_a]) <> 0 and Length(basis_coker[target_a]) <> 0 then
            mat := [];
            pos_a := Position(arrows,a);
            for b in basis_coker[source_a] do
                bb := PreImagesRepresentative(projection[source_a],b);
                bb := bb*mats[pos_a]; # computing bb^a
                Add(mat,Coefficients(basis_coker[target_a],Image(projection[target_a],bb)));
            od;
            Add(cokermat,[String(a),mat]);
        fi;
    od;
    
    dimvec_coker := List(basis_coker, basis -> Length(basis));
    
    if IsPathAlgebra(A) then
        C := RightModuleOverPathAlgebra(A, dimvec_coker, cokermat);
    else
        C := RightModuleOverPathAlgebra(A, dimvec_coker, cokermat);
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
    
    morph := RightModuleHomOverAlgebra(Range(f),C,map);
    SetCoKernelOfWhat(morph,f);
    SetIsSurjective(morph,true);
    
    return morph;
end
);

#######################################################################
##
#O  CoKernelOfAdditiveGeneralMapping( <f> )
##
##  This function returns a module isomorphic to the cokernel of  <f>. 
##
InstallMethod ( CoKernelOfAdditiveGeneralMapping, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    SUM_FLAGS+1,
    function( f )

    return Range(CoKernelProjection(f));
end
);

#######################################################################
##
#A  TopOfModuleProjection( <M> )
##
##  This function computes the map from  M  to the top of  M, 
##  M ---> M/rad(M).
##
InstallMethod( TopOfModuleProjection, 
    "for a pathalgebramatmodule",
    true, 
    [ IsPathAlgebraMatModule ], 0,
    function( M ) 

    local K, A, Q, vertices, num_vert, incomingarrows, mats, arrows, subspaces,
          i, a, dim_M, Vspaces, Wspaces, naturalprojections, index,
          dim_top, matrices, topofmodule, topofmoduleprojection, W;

    A := RightActingAlgebra(M);
    if Dimension(M) = 0 then 
        return ZeroMapping(M,M);
    else
        K := LeftActingDomain(A);
        Q := QuiverOfPathAlgebra(A);
        vertices := VerticesOfQuiver(Q);
        num_vert := Length(vertices);
        incomingarrows := List([1..num_vert], x -> IncomingArrowsOfVertex(vertices[x]));
        mats := MatricesOfPathAlgebraModule(M);
        arrows := ArrowsOfQuiver(Q);
        subspaces := List([1..num_vert], x -> []);
        for i in [1..num_vert] do
            for a in incomingarrows[i] do
                Append(subspaces[i], StructuralCopy(mats[Position(arrows,a)]));
            od;
        od;
        dim_M := DimensionVector(M);
        Vspaces := List([1..num_vert], x -> FullRowSpace(K,dim_M[x]));
        Wspaces := List([1..num_vert], x -> []);
        for i in [1..num_vert] do
            if dim_M[i] <> 0 then 
                Wspaces[i] := Subspace(Vspaces[i],subspaces[i]);
            else
                Wspaces[i] := Subspace(Vspaces[i],[]);
            fi;
        od;
        naturalprojections := List([1..num_vert], x -> NaturalHomomorphismBySubspace(Vspaces[x],Wspaces[x]));
        dim_top := List([1..num_vert], x -> Dimension(Range(naturalprojections[x]))); 
        index := function( n )
            if n = 0 then 
                return 1;
            else
                return n;
            fi;
        end;
        matrices := [];
        for i in [1..num_vert] do
            if dim_top[i] <> 0 then
                Add(matrices, List(BasisVectors(Basis(Vspaces[i])), y -> ImageElm(naturalprojections[i],y)));
            else
                Add(matrices, NullMat(index(dim_M[i]),1,K));
            fi;
        od;
        topofmodule := RightModuleOverPathAlgebra(A,dim_top,[]);
        topofmoduleprojection := RightModuleHomOverAlgebra(M,topofmodule,matrices);

        SetTopOfModule(M,topofmodule);

        return topofmoduleprojection;
    fi;
end
);

#######################################################################
##
#A  TopOfModule( <M> )
##
##  This function computes the top  M/rad(M)  of the module M. 
##
InstallMethod( TopOfModule, 
    "for a pathalgebramatmodule",
    true, 
    [ IsPathAlgebraMatModule ], 0,
    function( M )

    return Range(TopOfModuleProjection(M));
end
);

InstallOtherMethod ( TopOfModule, 
"for a PathAlgebraMatModuleMap",
[ IsPathAlgebraMatModuleHomomorphism ],
function( f )
    
  local M, N, K, pi_M, pi_N, BTopM, dim_vector_TopM, dim_vector_TopN, 
        n, h, i, dim, matrix, j, b, btilde, bprime;
  
  M := Source( f );
  N := Range( f );
  K := LeftActingDomain( M );
  pi_M := TopOfModuleProjection( M );
  pi_N := TopOfModuleProjection( N );
  BTopM := BasisVectors( Basis( Range( pi_M ) ) );
  dim_vector_TopM := DimensionVector( Range( pi_M ) );
  dim_vector_TopN := DimensionVector( Range( pi_N ) );    
  n := Length( dim_vector_TopM ); 
  h := [ ];
  for i in [ 1..n ] do
    if dim_vector_TopM[ i ] > 0 then
      if dim_vector_TopN[ i ] = 0 then
        Add( h, NullMat( dim_vector_TopM[ i ], 1, K ) );
      else
#
# Assuming that the basisvectors of topM are listed as 
# first basisvectors for vertex 1, vertex 2, .....
#
        dim := Sum( dim_vector_TopM{ [ 1..i - 1 ] } );
        matrix := [ ];
        for j in [ 1..dim_vector_TopM[ i ] ] do
          b := BTopM[ dim + j ];
          btilde := PreImagesRepresentative( pi_M, b );
          bprime := ImageElm( pi_N, ImageElm( f, btilde ) );
          Add( matrix, ExtRepOfObj( ExtRepOfObj( bprime ) )[ i ] ); 
        od;
        Add( h, matrix );
      fi;
    else
      if dim_vector_TopN[ i ] = 0 then
        Add( h, NullMat( 1, 1, K ) );
      else
        Add( h, NullMat( 1, dim_vector_TopN[ i ], K ) );
      fi;
    fi;        
  od;
  
  return RightModuleHomOverAlgebra( Range( pi_M ), Range( pi_N ), h );
end
  );

#######################################################################
##
#M  \+( <f>, <g> )
##
##  This function returns the sum of two homomorphisms  <f>  and  <g>,
##  when the sum is defined, otherwise it returns an error message. 
##
InstallMethod( \+,
    "for two PathAlgebraMatModuleMap's",
    true,
#  IsIdenticalObj,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f, g )
    
    local i, num_vert, x, Fam;

    if ( Source(f) = Source(g) ) and ( Range(f) = Range(g) ) then 
     	num_vert := Length(f!.maps);
        x := List([1..num_vert], y -> f!.maps[y] + g!.maps[y]);
        return RightModuleHomOverAlgebra(Source(f),Range(f),x);
    else
	Error("the two arguments entered do not live in the same homomorphism set, ");
    fi;
end
);
  
#######################################################################
##
#M  \*( <f>, <g> )
##
##  This function returns the composition  <f*g>  of two homomorphisms  
##  <f>  and  <g>, that is, first the map  <f> then followed by  <g>,  
##  when the composition is defined, otherwise it returns an error 
##  message. 
##
InstallMethod( \*,
    "for two PathAlgebraMatModuleMap's",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f, g )
    
    local i, num_vert, x;

    if Range(f) = Source(g) then 
     	num_vert := Length(f!.maps);
        x := List([1..num_vert], y -> f!.maps[y]*g!.maps[y]);
	return RightModuleHomOverAlgebra(Source(f),Range(g),x);
    else
        Error("codomain of the first argument is not equal to the domain of the second argument, ");
    fi;
end
);
  
#######################################################################
##
#M  \*( <a>, <g> )
##
##  This function returns the scalar multiple  <a*g>  of a scalar  <a> 
##  with a homomorphism <g>, when this scalar multiplication is defined, 
##  otherwise it returns an error message. 
##
InstallOtherMethod( \*,
    "for two PathAlgebraMatModuleMap's",
    true,
    [ IsScalar, IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( a, g )
    
    local K, i, num_vert, x;

    K := LeftActingDomain(Source(g));
    if a in K then 
     	num_vert := Length(g!.maps);
        x := List([1..num_vert], y -> a*g!.maps[y]);
        return RightModuleHomOverAlgebra(Source(g),Range(g),x);
     else
         Error("the scalar is not in the same field as the algbra is over,");
     fi;
end
);
  
#######################################################################
##
#M  AdditiveInverseOp( <f> )
##
##  This function returns the additive inverse of the homomorphism 
##  <f>. 
##
InstallMethod( AdditiveInverseOp,
    "for a morphism in IsPathAlgebraMatModuleHomomorphism",
    [ IsPathAlgebraMatModuleHomomorphism ],

    function ( f ) 

    local i, num_vert, x;

    num_vert := Length(f!.maps);
    x := List([1..num_vert], y -> (-1)*f!.maps[y]);

    return RightModuleHomOverAlgebra(Source(f),Range(f),x);
end
);

#######################################################################
##
#M  \*( <f>, <a> )
##
##  This function returns the scalar multiple  <f*a>  of a homomorphism
##  <f>  and a scalar  <a>. 
##
InstallOtherMethod( \*,
    "for two PathAlgebraMatModuleMap's",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsScalar ],
    0,
    function( f, a )
    
    local K, i, num_vert, x;

    K := LeftActingDomain(Source(f));
    if a in K then 
     	num_vert := Length(f!.maps);
        x := List([1..num_vert], y -> f!.maps[y]*a);
	return RightModuleHomOverAlgebra(Source(f),Range(f),x);
    else
 	Error("the scalar is not in the same field as the algbra is over,");
    fi;
end
);
  
#######################################################################
##
#O  HomOverAlgebra( <M>, <N> )
##
##  This function computes a basis of the vector space of homomorphisms
##  from the module  <M>  to the module  <N>. The algorithm it uses is
##  based purely on linear algebra.
##
InstallMethod( HomOverAlgebra,
    "for two representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
    function( M, N )

    local A, F, dim_M, dim_N, num_vert, support_M, support_N, num_rows, num_cols, 
          block_rows, block_cols, block_intervals, 
          i, j, equations, arrows, vertices, v, a, source_arrow, target_arrow, 
          mats_M, mats_N, prev_col, prev_row, row_start_pos, col_start_pos, 
          row_end_pos, col_end_pos, l, m, n, hom_basis, map, mat, homs, x, y, k, b, 
          dim_hom, zero;

    A := RightActingAlgebra(M); 
    if A <> RightActingAlgebra(N) then
        Print("The two modules entered are not modules over the same algebra.");
        return fail;
    fi;
    F := LeftActingDomain(A);
    #
    # Finding the support of M and N
    # 
    vertices := VerticesOfQuiver(QuiverOfPathAlgebra(OriginalPathAlgebra(A)));
    dim_M := DimensionVector(M);
    dim_N := DimensionVector(N);
    num_vert := Length(dim_M);   
    support_M := [];
    support_N := [];
    for i in [1..num_vert] do
        if (dim_M[i] <> 0) then 
            AddSet(support_M,i);
        fi;
        if (dim_N[i] <> 0) then 
            AddSet(support_N,i);
        fi;
    od;
    #
    # Deciding the size of the equations, 
    # number of columns and rows
    #
    num_cols := 0;
    num_rows := 0;
    block_intervals := [];
    block_rows := [];
    block_cols := [];
    prev_col := 0;
    prev_row := 0;
    for i in support_M do
        num_rows := num_rows + dim_M[i]*dim_N[i];
        block_rows[i] := prev_row+1;
        prev_row := num_rows;
        for a in OutgoingArrowsOfVertex(vertices[i]) do
            source_arrow := Position(vertices,SourceOfPath(a));
            target_arrow := Position(vertices,TargetOfPath(a));
            if (target_arrow in support_N) and ( (source_arrow in support_N) or (target_arrow in support_M)) then 
                num_cols := num_cols + dim_M[source_arrow]*dim_N[target_arrow];
                Add(block_cols,[a,prev_col+1,num_cols]);
            fi;
            prev_col := num_cols; 
        od;
    od;
    #
    # Finding the linear equations for the maps between M and N
    #
    equations := NullMat(num_rows, num_cols, F);

    arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(OriginalPathAlgebra(A)));
    mats_M := MatricesOfPathAlgebraModule(M);
    mats_N := MatricesOfPathAlgebraModule(N);
    prev_col := 0;
    prev_row := 0;
    for i in support_M do
        for a in OutgoingArrowsOfVertex(vertices[i]) do
            source_arrow := Position(vertices,SourceOfPath(a));
            target_arrow := Position(vertices,TargetOfPath(a));
            if (target_arrow in support_N) and ( (source_arrow in support_N) or (target_arrow in support_M)) then
                for j in [1..dim_M[source_arrow]] do
                    row_start_pos := block_rows[source_arrow] + (j-1)*dim_N[source_arrow]; 
                    row_end_pos := block_rows[source_arrow] - 1 + j*dim_N[source_arrow];
                    col_start_pos := prev_col + 1 + (j-1)*dim_N[target_arrow];
                    col_end_pos := prev_col + j*dim_N[target_arrow];
                    if (source_arrow in support_N) then 
                        equations{[row_start_pos..row_end_pos]}{[col_start_pos..col_end_pos]} := mats_N[Position(arrows,a)];
                    fi;
                    if (target_arrow in support_M) then 
                        for m in [1..DimensionsMat(mats_M[Position(arrows,a)])[2]] do
                            for n in [1..dim_N[target_arrow]] do
                                b := block_rows[target_arrow]+(m-1)*dim_N[target_arrow];
                                equations[b+n-1][col_start_pos+n-1] := equations[b+n-1][col_start_pos+n-1]+(-1)*mats_M[Position(arrows,a)][j][m];
                            od;
                        od;
                    fi;
                od;
                prev_col := prev_col + dim_M[source_arrow]*dim_N[target_arrow];
            fi;
        od;
    od;
    #
    # Creating the maps between the module M and N
    #
    homs := [];
    if (num_rows <> 0) and (num_cols <> 0) then 
        dim_hom := 0; 
        hom_basis := NullspaceMat(equations);
        for b in hom_basis do
            map := [];
            dim_hom := dim_hom + 1;
            k := 1;
            for i in [1..num_vert] do 
                if dim_M[i] = 0 then 
                    if dim_N[i] = 0 then 
                        Add(map,NullMat(1,1,F));
                    else
                        Add(map,NullMat(1,dim_N[i],F));
                    fi;
                else
                    if dim_N[i] = 0 then 
                        Add(map,NullMat(dim_M[i],1,F));
                    else
                        mat := NullMat(dim_M[i],dim_N[i], F);
                        for y in [1..dim_M[i]] do 
                            for x in [1..dim_N[i]] do 
                                mat[y][x] := b[k];
                                k := k + 1;
                            od;
                        od;
                        Add(map,mat);
                    fi;
                fi;
            od;
            homs[dim_hom] := Objectify( NewType( CollectionsFamily( GeneralMappingsFamily(
                                     ElementsFamily( FamilyObj( M ) ),
                                     ElementsFamily( FamilyObj( N ) ) ) ), 
                                     IsPathAlgebraMatModuleHomomorphism and IsPathAlgebraMatModuleHomomorphismRep and IsAttributeStoringRep ), rec( maps := map ));
            SetPathAlgebraOfMatModuleMap(homs[dim_hom], A);
            SetSource(homs[dim_hom], M);
            SetRange(homs[dim_hom], N);
            SetIsWholeFamily(homs[dim_hom],true);
        od;
        return homs;
    else
        homs := [];
        if Dimension(M) = 0 or Dimension(N) = 0 then 
            return homs;
        else 
            dim_hom := 0;
            zero := [];
            for i in [1..num_vert] do
                if dim_M[i] = 0 then 
                    if dim_N[i] = 0 then 
                        Add(zero,NullMat(1,1,F));
                    else
                        Add(zero,NullMat(1,dim_N[i],F));
                    fi;
                else
                    if dim_N[i] = 0 then 
                        Add(zero,NullMat(dim_M[i],1,F));
                    else
                        Add(zero,NullMat(dim_M[i],dim_N[i],F));
                    fi;
                fi;
            od;      
            for i in [1..num_vert] do
                if (dim_M[i] <> 0) and (dim_N[i] <> 0) then 
                    for m in BasisVectors(Basis(FullMatrixSpace(F,dim_M[i],dim_N[i]))) do
                        dim_hom := dim_hom + 1;
                        homs[dim_hom] := ShallowCopy(zero);
                        homs[dim_hom][i] := m;
                    od;
                fi;
            od;
            for i in [1..dim_hom] do 
                homs[i] := Objectify( NewType( CollectionsFamily( GeneralMappingsFamily(
                                   ElementsFamily( FamilyObj( M ) ),
                                   ElementsFamily( FamilyObj( N ) ) ) ), 
                                   IsPathAlgebraMatModuleHomomorphism and IsPathAlgebraMatModuleHomomorphismRep and IsAttributeStoringRep ), rec( maps := homs[i] ));
                SetPathAlgebraOfMatModuleMap(homs[i], A);
                SetSource(homs[i], M);
                SetRange(homs[i], N);
                SetIsWholeFamily(homs[i],true);
            od;

            return homs;
        fi;
    fi;
end
);

#######################################################################
##
#O  HomOverAlgebraWithBasisFunction( <M>, <N> )
##
##  This function returns a list with two entries. The first entry is 
##  a basis of the vector space of homomorphisms from the module  <M>
##  to the module  <N>.  The second entry is a function from the space 
##  of homomorphisms from  <M>  to  <N>  to the vector space with the 
##  given by the first entry.  The algorithm it uses is based purely 
##  on linear algebra.
##
InstallMethod( HomOverAlgebraWithBasisFunction,
    "for two representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
    function( M, N )

    local homs, F, basis, V, W, HomB, transfermat, coefficients;


    homs := HomOverAlgebra( M, N );
    F := LeftActingDomain( RightActingAlgebra( M ) );
    
    if Length( homs ) > 0 then
      basis := List( homs, h -> Flat( h!.maps ) );
      V := FullRowSpace( F, Length( basis[ 1 ] ) );
      W := Subspace( V, basis, "basis" );
      HomB := Basis( W );
      transfermat := List( basis, b -> Coefficients( HomB, b ) );
      coefficients := function( h )
        local hflat, coeff;
               
        hflat := Flat( h!.maps );
        coeff := Coefficients( HomB, hflat );
        
        return coeff * transfermat^(-1);
      end;
    else
      coefficients := function( h );
        return [ Zero( F ) ];
      end;
    fi;
    
    return [ homs, coefficients ];
end
);

#######################################################################
##
#A  EndOverAlgebra( <M>, <N> )
##
##  This function computes endomorphism ring of the module  <M>  and
##  representing it as an general GAP algebra. The algorithm it uses is
##  based purely on linear algebra.
##
InstallMethod( EndOverAlgebra,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M )

    local EndM, R, F, dim_M, alglist, i, j, r, maps, A; 

    EndM := HomOverAlgebra(M,M);
    R := RightActingAlgebra(M); 
    F := LeftActingDomain(R);
    dim_M := DimensionVector(M);
    alglist := [];
    for i in [1..Length(dim_M)] do 
        if dim_M[i] <> 0 then 
            Add(alglist, MatrixAlgebra(F,dim_M[i]));
        fi;
    od;
    maps := [];
    for i in [1..Length(EndM)] do
        maps[i] := NullMat(Dimension(M),Dimension(M),F);
        r := 1; 
        for j in [1..Length(dim_M)] do 
            if dim_M[j] <> 0 then 
                maps[i]{[r..r+dim_M[j]-1]}{[r..r+dim_M[j]-1]} := EndM[i]!.maps[j];
            fi;
            r := r + dim_M[j];
        od; 
    od;
    A := DirectSumOfAlgebras(alglist); 
    
    return SubalgebraWithOne(A,maps,"basis"); 
end
);

#######################################################################
##
#A  RightFacApproximation( <M>, <C> )
##
##  This function computes a right Fac<M>-approximation of the module 
##  <C>. 
##
InstallMethod( RightFacApproximation,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
    function( M, C )

    local homMC, i, generators; 

    homMC := HomOverAlgebra(M,C); 
    generators := [];
    for i in [1..Length(homMC)] do
        Append(generators,ImagesSet(homMC[i],Source(homMC[i])));
    od;
    
    return SubRepresentationInclusion(C,generators);
end
);

#######################################################################
##
#O  NumberOfNonIsoDirSummands( <M> )
##
##  This function computes the number of non-isomorphic indecomposable 
##  direct summands of the module  <M>, and in addition returns the 
##  dimensions of the simple blocks of the semisimple ring  
##  End(M)/rad End(M). 
##
InstallMethod( NumberOfNonIsoDirSummands,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M )
 
    local EndM, K, J, gens, I, A, top, AA, B, n,
          i, j, genA, V, W, d;

    EndM := EndOverAlgebra(M);
    K := LeftActingDomain(M);
    J := RadicalOfAlgebra(EndM);
    gens := GeneratorsOfAlgebra(J);
    I := Ideal(EndM,gens); 
    A := EndM/I;
    top := CentralIdempotentsOfAlgebra(A);
    
    return [Length(top),List(DirectSumDecomposition(EndM/I),Dimension)];
end
);

#######################################################################
##
#A  DualOfModuleHomomorphism( <f> )
##
##  This function computes the dual of a homomorphism from the module  
##  <M>  to the module  <N>.
##
InstallMethod ( DualOfModuleHomomorphism,
    "for a map between representations of a quiver",
    [ IsPathAlgebraMatModuleHomomorphism ], 0,
    function( f )

    local mats, M, N;
   
    mats := f!.maps;
    mats := List(mats, x -> TransposedMat(x));
    M := DualOfModule(Source(f));
    N := DualOfModule(Range(f));

    return RightModuleHomOverAlgebra(N,M,mats);
end
);

#######################################################################
##
#A  SocleOfModuleInclusion( <M> )
##
##  This function computes the map from the socle of  M  to the module M, 
##  soc(M) ---> M.
##
InstallMethod( SocleOfModuleInclusion, 
    "for a pathalgebramatmodule",
    true, 
    [ IsPathAlgebraMatModule ], 0,
    function( M )

    local A, K, Q, vertices, num_vert, outgoingarrows, mats, arrows, dim_M, 
          subspaces, i, a, j, socle, matrixfunction, dim_socle, socleofmodule, 
          socleinclusion, V, temp;

    A := RightActingAlgebra(M);
    if Dimension(M) = 0 then 
        return ZeroMapping(M,M);
    else
        K := LeftActingDomain(A);
        Q := QuiverOfPathAlgebra(A);
        vertices := VerticesOfQuiver(Q);
        num_vert := Length(vertices);
        outgoingarrows := List([1..num_vert], x -> OutgoingArrowsOfVertex(vertices[x]));
        dim_M := DimensionVector(M);
        mats := MatricesOfPathAlgebraModule(M);
        arrows := ArrowsOfQuiver(Q);
        dim_socle := List([1..num_vert], x -> 0);
        socle := List([1..num_vert], x -> []);
        for i in [1..num_vert] do
            if ( Length(outgoingarrows[i]) = 0 ) or ( dim_M[i] = 0 ) then
                dim_socle[i] := dim_M[i];
                if  dim_M[i] = 0 then
                    socle[i] := NullMat(1,1,K);
                else
                    socle[i] := IdentityMat(dim_M[i],K);
                fi;
            else
                subspaces := List([1..dim_M[i]], y -> []);   
                for a in outgoingarrows[i] do
                    for j in [1..dim_M[i]] do
                        Append(subspaces[j], StructuralCopy(mats[Position(arrows,a)][j]));
                    od;
                od;
                V := FullRowSpace(K,dim_M[i]);
                temp := NullspaceMat(subspaces);
                dim_socle[i] := Dimension(Subspace(V,temp));
                if dim_socle[i] = 0 then
                    socle[i] := NullMat(1,dim_M[i],K);
                else
                    socle[i] := temp;
                fi;
            fi;
        od;
        socleofmodule := RightModuleOverPathAlgebra(A,dim_socle,[]);
        socleinclusion := RightModuleHomOverAlgebra(socleofmodule,M,socle);
#        SetSocleOfModule(M,socleofmodule);
        return socleinclusion;
    fi;
end
  );

#######################################################################
##
#A  SocleOfModule( <M> )
##
##  This function computes the socle  soc(M)  of the module M. 
##
InstallMethod( SocleOfModule, 
    "for a pathalgebramatmodule",
    true, 
    [ IsPathAlgebraMatModule ], 0,
    function( M )

    return Source(SocleOfModuleInclusion(M));
end
);

#######################################################################
##
#O  CommonDirectSummand( <M>, <N> )
##
##  This function is using the algorithm for finding a common direct 
##  summand presented in the paper "Gauss-Elimination und der groesste
##  gemeinsame direkte Summand von zwei endlichdimensionalen Moduln"
##  by K. Bongartz, Arch Math., vol. 53, 256-258, with the modification
##  done by Andrzej Mroz found in "On the computational complexity of Bongartz's
##  algorithm" (improving the complexity of the algorithm).
##
InstallMethod( CommonDirectSummand, 
    "for two path algebra matmodules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule  ], 0,
    function( M, N ) 

    local   HomMN,  HomNM,  mn,  nm,  m,  n,  l,  zero,  j,  i,  temp,  
          r,  f,  fnm,  nmf;

    if RightActingAlgebra( M ) <> RightActingAlgebra( N ) then 
        Print( "The two modules are not modules over the same algebra.\n" );
        return fail;
    else
        HomMN := HomOverAlgebra( M, N );
        HomNM := HomOverAlgebra( N, M );
        mn := Length( HomMN );
        nm := Length( HomNM );
      
        if mn = 0 or nm = 0 then 
            return false;
        fi;
      
        m := Maximum( DimensionVector( M ) );
        n := Maximum( DimensionVector( N ) );
        if n = m then
            l := n;
        else
            l := Minimum( [ n, m ] ) + 1;
        fi;
        
        n := Int( Ceil( Log2( 1.0*( l ) ) ) );
          
        zero := ZeroMapping( M, M );
      
        for j in [ 1..nm ] do
            for i in [ 1..mn ] do
                if l > 1 then  # because hom^0 * hom => error! 
                    temp := HomMN[ i ] * HomNM[ j ];
                    for r in [ 1..n ] do
                        temp := temp * temp;
                    od;
                    f := temp * HomMN[ i ];
                else 
                    f := HomMN[ i ];
                fi;
                
                fnm := f * HomNM[ j ];
              
                if fnm <> zero then
                    nmf := HomNM[ j ] * f; 
                    return [ Image( fnm ), Kernel( fnm ), Image( nmf ), Kernel( nmf ) ];
                fi;
            od;
        od;

        return false;
    fi; 
end
);

#######################################################################
##
#O  MaximalCommonDirectSummand( <M>, <N> )
##
##  This function is using the algorithm for finding a maximal common 
##  direct summand based on the algorithm presented in the paper 
##  "Gauss-Elimination und der groesste gemeinsame direkte Summand von 
##  zwei endlichdimensionalen Moduln" by K. Bongartz, Arch Math., 
##  vol. 53, 256-258, with the modification done by Andrzej Mroz found 
##  in "On the computational complexity of Bongartz's algorithm" 
##  (improving the complexity of the algorithm).
##
InstallMethod( MaximalCommonDirectSummand, 
    "for two path algebra matmodules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule  ], 0,
    function( M, N ) 

    local U, V, maxcommon, L;

    U := M;
    V := N;
    maxcommon := [];
    repeat 
        L := CommonDirectSummand(U,V);
        if L <> false and L <> fail then 
            Add(maxcommon,L[1]);
            U := L[2];
            V := L[4];
            if Dimension(L[2]) = 0 or Dimension(L[4]) = 0 then
                break;
            fi;
        fi;
    until  L = false or L = fail;
    
    if Length(maxcommon) = 0 then 
        return false;
    else 
        return [maxcommon,U,V];
    fi;     
end
);

#######################################################################
##
#O  IsomorphicModules( <M>, <N> )
##
##  This function returns true if the modules  <M>  and  <N>  are 
##  isomorphic, an error message if  <M>  and  <N>  are not modules over 
##  the same algebra and false otherwise.
##  
InstallMethod( IsomorphicModules, 
    "for two path algebra matmodules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule  ], 0,
    function( M, N ) 

    local L;

    if DimensionVector(M) <> DimensionVector(N) then 
        return false;
    elif Dimension(M) = 0 and Dimension(N) = 0 then 
        return true; 
    else 
        L := MaximalCommonDirectSummand(M,N);
        if L = false then 
            return false;
        else
            if Dimension(L[2]) = 0 and Dimension(L[3]) = 0 then 
                return true;
            else
                return false;
            fi;
        fi;
    fi;
end
); 

#######################################################################
##
#O  IsDirectSummand( <M>, <N> )
##
##  This function returns true if the module  <M>  is isomorphic to a 
##  direct of the module  <N>, an error message if  <M>  and  <N>  are 
##  not modules over the same algebra and false otherwise.
##  
InstallMethod( IsDirectSummand, 
    "for two path algebra matmodules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule  ], 0,
    function( M, N ) 

    local L;

    if not DimensionVectorPartialOrder(M,N) then 
        return false;
    else 
        L := MaximalCommonDirectSummand(M,N);
        if L = false then 
            return false;
        else 
            if Dimension(L[2]) = 0 then 
                return true;
            else
                return false;
            fi;
        fi;
    fi;
end
); 

#######################################################################
##
#O  IsInAdditiveClosure( <M>, <N> )
##
##  This function returns true if the module  <M>  is in the additive
##  closure of the module  <N>, an error message if  <M>  and  <N>  are 
##  not modules over the same algebra and false otherwise.
##
InstallMethod( IsInAdditiveClosure, 
    "for two path algebra matmodules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule  ], 0,
    function( M, N ) 

    local K, HomMN, HomNM, MM, i, j, HomMM, V_M;

    if RightActingAlgebra(M) <> RightActingAlgebra(N) then 
        return fail;
    else
#
# Computing Hom(M,N) and Hom(N,M), finding the subspace in Hom(M,M)
# spanned by Hom(M,N)*Hom(M,N), if they have the same dimension, then 
# the identity on  M  is in the linear span of Hom(M,N)*Hom(M,N) and  
# module M is in the additive closure of N. 
#
        K := LeftActingDomain(M);
        HomMN := HomOverAlgebra(M,N);
        HomNM := HomOverAlgebra(N,M);
        MM := [];
        for i in [1..Length(HomMN)] do
            for j in [1..Length(HomNM)] do
                Add(MM,HomMN[i]*HomNM[j]);
            od;
        od;
        MM := List(MM,x->x!.maps);
        for i in [1..Length(MM)] do
            MM[i] := List(MM[i],x->Flat(x)); 
            MM[i] := Flat(MM[i]);
        od;
        HomMM:=HomOverAlgebra(M,M);
        if Length(MM) = 0 then 
            V_M := TrivialSubspace(K);
        else 
            V_M := VectorSpace(K,MM);
        fi; 
        if Dimension(V_M) = Length(HomOverAlgebra(M,M)) then
            return true;
        else
            return false;
        fi;
    fi;
end
);

#######################################################################
##
#O  MorphismOnKernel( <f>, <g>, <beta>, <alpha> )
##
##  Given the following commutative diagram
##          B -- f --> C
##          |          |
##     beta |          | alpha
##          V          V 
##          B' - g --> C' 
##  this function finds the induced homomorphism on the kernels of 
##  f and g. 
##
InstallMethod ( MorphismOnKernel, 
    "for commutative diagram of PathAlgebraMatModuleMaps",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 
    0,
    function( f, g, beta, alpha )

    local K, kerf, kerg, B, dim_vec_kerf, dim_vec_kerg, map, j, i, p, mat;
#
#  Checking if a commutative diagram was entered
# 
    if f*alpha = beta*g then 
# 
#
        K := LeftActingDomain(Source(f));
        kerf := KernelInclusion(f);
        kerg := KernelInclusion(g);
#
#  Computing the information needed for lifting the morphism.
#
        B := BasisVectors(Basis(Source(kerf)));
        B := List(B, x -> ImageElm(kerf*beta,x));
        B := List(B, x -> PreImagesRepresentative(kerg,x));
#
#  Computing dimension vectors so that we can insert zero matrices of 
#  the right size. 
#
        dim_vec_kerf := DimensionVector(Source(kerf));
        dim_vec_kerg := DimensionVector(Source(kerg));
        map := [];
        j := 0;
        for i in [1..Length(dim_vec_kerf)] do
#
#  If the kernel of  f  is zero in vertex  i, then insert a zero matrix of
#  the right size, do not use any of the lifting information.
# 
            if dim_vec_kerf[i] = 0 then
                if  dim_vec_kerg[i] = 0 then 
                    Add(map,NullMat(1,1,K));
                else
                    Add(map,NullMat(1,dim_vec_kerg[i],K));
                fi;
            else
#
#  If the kernel of  f  is non-zero in vertex  i, then use the lifting 
#  information to compute the right matrix for the map from vertex  i  to  i.
# 
                mat := [];
                for p in [1..dim_vec_kerf[i]] do
                    j := j + 1;
                    Add(mat,ExtRepOfObj(ExtRepOfObj(B[j]))[i]);
                od;
                Add(map,mat);
            fi;
        od;
        return RightModuleHomOverAlgebra(Source(kerf),Source(kerg),map);
    else 
       return fail;
    fi;
end
);

#######################################################################
##
#O  MorphismOnImage( <f>, <g>, <beta>, <alpha> )
##
##  Given the following commutative diagram
##          B -- f --> C
##          |          |
##     beta |          | alpha
##          V          V 
##          B' - g --> C' 
##  this function finds the induced homomorphism on the images of 
##  f and g. 
##
InstallMethod ( MorphismOnImage, 
    "for commutative diagram of PathAlgebraMatModuleMaps",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 
    0,
    function( f, g, beta, alpha )

    local K, imagef, imageg, B, dim_vec_imagef, dim_vec_imageg, map, j, i, p, mat;
#
#  Checking if a commutative diagram was entered
# 
    if f*alpha = beta*g then 
# 
#
        K := LeftActingDomain(Source(f));
        imagef := ImageProjection(f);
        imageg := ImageProjection(g);
#
#  Computing the information needed for lifting the morphism.
#
        B := BasisVectors(Basis(Range(imagef)));
        B := List(B, x -> PreImagesRepresentative(imagef,x));
        B := List(B, x -> ImageElm(beta*imageg,x));
#
#  Computing dimension vectors so that we can insert zero matrices of 
#  the right size. 
#
        dim_vec_imagef := DimensionVector(Range(imagef));
        dim_vec_imageg := DimensionVector(Range(imageg));
        map := [];
        j := 0;
        for i in [1..Length(dim_vec_imagef)] do
#
#  If the image of  f  is zero in vertex  i, then insert a zero matrix of
#  the right size, do not use any of the lifting information.
# 
            if dim_vec_imagef[i] = 0 then
                if  dim_vec_imageg[i] = 0 then 
                    Add(map,NullMat(1,1,K));
                else
                    Add(map,NullMat(1,dim_vec_imageg[i],K));
                fi;
            else
#
#  If the cokernel of  f  is non-zero in vertex  i, then use the lifting 
#  information to compute the right matrix for the map from vertex  i  to  i.
# 
                mat := [];
                for p in [1..dim_vec_imagef[i]] do
                    j := j + 1;
                    Add(mat,ExtRepOfObj(ExtRepOfObj(B[j]))[i]);
                od;
                Add(map,mat);
            fi;
        od;
        return RightModuleHomOverAlgebra(Range(imagef),Range(imageg),map);
    else 
        return fail;
    fi;
end
);

#######################################################################
##
#O  MorphismOnCoKernel( <f>, <g>, <beta>, <alpha> )
##
##  Given the following commutative diagram
##          B -- f --> C
##          |          |
##     beta |          | alpha
##          V          V 
##          B' - g --> C' 
##  this function finds the induced homomorphism on the cokernels of 
##  f and g. 
##
InstallMethod ( MorphismOnCoKernel, 
    "for commutative diagram of PathAlgebraMatModuleMaps",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 
    0,
    function( f, g, beta, alpha )

    local K, cokerf, cokerg, B, dim_vec_cokerf, dim_vec_cokerg, map, j, i, p, mat;
#
#  Checking if a commutative diagram was entered
# 
    if f*alpha = beta*g then 
# 
#
        K := LeftActingDomain(Source(f));
        cokerf := CoKernelProjection(f);
        cokerg := CoKernelProjection(g);
#
#  Computing the information needed for lifting the morphism.
#
        B := BasisVectors(Basis(Range(cokerf)));
        B := List(B, x -> PreImagesRepresentative(cokerf,x));
        B := List(B, x -> ImageElm(alpha*cokerg,x));
#
#  Computing dimension vectors so that we can insert zero matrices of 
#  the right size. 
#
        dim_vec_cokerf := DimensionVector(Range(cokerf));
        dim_vec_cokerg := DimensionVector(Range(cokerg));
        map := [];
        j := 0;
        for i in [1..Length(dim_vec_cokerf)] do
#
#  If the cokernel of  f  is zero in vertex  i, then insert a zero matrix of
#  the right size, do not use any of the lifting information.
# 
            if dim_vec_cokerf[i] = 0 then
                if  dim_vec_cokerg[i] = 0 then 
                    Add(map,NullMat(1,1,K));
                else
                    Add(map,NullMat(1,dim_vec_cokerg[i],K));
                fi;
            else
#
#  If the cokernel of  f  is non-zero in vertex  i, then use the lifting 
#  information to compute the right matrix for the map from vertex  i  to  i.
# 
                mat := [];
                for p in [1..dim_vec_cokerf[i]] do
                    j := j + 1;
                    Add(mat,ExtRepOfObj(ExtRepOfObj(B[j]))[i]);
                od;
                Add(map,mat);
            fi;
        od;
        return RightModuleHomOverAlgebra(Range(cokerf),Range(cokerg),map);
    else 
        return fail;
    fi;
end
);


#######################################################################
##
#O  LiftingMorphismFromProjective( <f>, <g> )
##
##  Given the following diagram
##                P
##                |
##                | g
##            f   V 
##        B ----> C 
##  where  P  is a direct sum of indecomposable projective modules 
##  constructed via DirectSumOfQPAModules and  f  an epimorphism, 
##  this function finds a lifting of  g  to  B.
##
InstallMethod ( LiftingMorphismFromProjective, 
    "for two PathAlgebraMatModuleMaps",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 
    0,
    function( f, g )

    local P, K, B, inclusions, projections, i, m, hmap; 
#
#  Checking if the input is as required
# 
    P := Source(g);
    if not IsProjectiveModule(P) then
        Error( "The source of the second map is not projective.\n" );
    fi;
    if not IsSurjective(f) then
        Error( "The first map is not surjective.\n" );        
    fi;
    if not ( Range(f) = Range(g) ) then
        Error( "The entered maps do not have a common range.\n" );
    fi;
    if IsZero( P ) then
        return ZeroMapping( P, Source( f ) );
    fi;
    if IsDirectSumOfModules(P) then  	  
#   
#
        K := LeftActingDomain(P);
        inclusions  := DirectSumInclusions(P);
        projections := DirectSumProjections(P);
        hmap := [];
#
#  First construct the lifting from each indecomposable direct summand of  P.
# 
        for i in [1..Length(inclusions)] do
            m := MinimalGeneratingSetOfModule(Source(inclusions[i]))[1];
            if ImageElm(g,ImageElm(inclusions[i],m)) = Zero(Range(g)) then 
                Add(hmap,ZeroMapping(Source(inclusions[i]),Source(f)));
            else 
                m := PreImagesRepresentative(f,ImageElm(g,ImageElm(inclusions[i],m)));
                Add(hmap,HomFromProjective(m,Source(f)));
            fi;
        od;
#
#  Make sure that the partial liftings start in the right modules/variables. 
#
        hmap := List([1..Length(inclusions)], x -> RightModuleHomOverAlgebra(Source(inclusions[x]),Source(f),hmap[x]!.maps));
        return projections*hmap;
    else 
        return fail;
    fi;
end
);

#######################################################################
##
#O  LiftingInclusionMorphisms( <f>, <g> )
##
##  Given the following diagram
##                A
##                |
##                | g
##            f   V 
##        B ----> C 
##  where  f  and  g  are inclusions, this function constructs a 
##  morphism (inclusion) from  A  to  B, whenever the image of  g  is
##  contained in the image of  f.  Otherwise the function returns fail.
##
InstallMethod ( LiftingInclusionMorphisms, 
    "for two PathAlgebraMatModuleMaps",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 
    0,
    function( f, g )

    local pi, K, B, dim_vec_sourcef, dim_vec_sourceg, map, i, j, mat, p; 
#
#  Checking if the image of  g  is contained in the image of  f.
# 
    pi := CoKernelProjection(f);
    if ( IsInjective(f) and IsInjective(g) ) and ( Range(f) = Range(g) ) and ( g*pi = ZeroMapping(Source(g),Range(pi)) ) then 
# 
#
        K := LeftActingDomain(Source(g));
        B := BasisVectors(Basis(Source(g)));
        B := List(B, x -> PreImagesRepresentative(f,ImageElm(g,x)));
#
#  Computing dimension vectors so that we can insert zero matrices of 
#  the right size. 
#
        dim_vec_sourceg := DimensionVector(Source(g));
        dim_vec_sourcef := DimensionVector(Source(f));
        map := [];
        j := 0;
        for i in [1..Length(dim_vec_sourceg)] do
#
#  If the source of  g  is zero in vertex  i, then insert a zero matrix of
#  the right size, do not use any of the lifting information.
# 
            if dim_vec_sourceg[i] = 0 then
                if  dim_vec_sourcef[i] = 0 then 
                    Add(map,NullMat(1,1,K));
                else
                    Add(map,NullMat(1,dim_vec_sourcef[i],K));
                fi;
            else
#
#  If the source of  g  is non-zero in vertex  i, then use the lifting 
#  information to compute the right matrix for the map from vertex  i  to  i.
# 
                mat := [];
                for p in [1..dim_vec_sourceg[i]] do
                    j := j + 1;
                    Add(mat,ExtRepOfObj(ExtRepOfObj(B[j]))[i]);
                od;
                Add(map,mat);
            fi;
        od;
        return RightModuleHomOverAlgebra(Source(g),Source(f),map);
    else 
        return fail;
    fi;
end
);

#######################################################################
##
#O  IntersectionOfSubmodules( <f>, <g> )
##
##                                              f              g
##  Given two submodules/subrepresentations  M ---> X  and  N ---> X  of 
##  a module  X  by two monomorphism  f  and  g, this function computes 
##  the intersection of the images of  f  and  g  in the sense that it 
##  computes the pullback of the diagram
##             f'
##          E ---> N
##          |      |
##       g' |      | g 
##          |      |
##          V  f   V     
##          M ---> X
##  
##  and it returns the maps [ f', g', g'*f ]
##
InstallOtherMethod ( IntersectionOfSubmodules, 
    "for two IsPathAlgebraMatModuleHomomorphisms",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 
    0,
    function( f, g )

    local pullback;
#
#   Checking the input
#
    if not IsInjective(f) then 
        Error("the first argument is not a monomorphism,");
    fi;
    if not IsInjective(g) then 
        Error("the second argument is not a monomorphism,");
    fi;
    if not ( Range(f) = Range(g) ) then 
        Error("must have two submodules of the same module,");
    fi;
#
#   Doing the computations
#
    pullback := PullBack(f,g);
    
    return [pullback[1],pullback[2],pullback[2]*f];
end
);

#######################################################################
##
#O  IntersectionOfSubmodules( <args> )
##
##                                            f_i             
##  Given submodules/subrepresentations  M_i -----> X  for i = 1,...,n 
##  of a module  X  by  n  monomorphism  f_i, this function computes 
##  the intersection of the images of  f  and  g  in the sense that it 
##  computes the pullbacks of the diagrams
##             f_1'
##         E_1 ----> N           E_2 ---------> N
##          |        |            |             |
##     f_2' |        | f_2   f_3' |             | f_3  and so on.
##          |        |            |             |
##          V  f_1   V            V   f_2'*f_1  V
##          M -----> X           E_1 ---------> X
##  
##  It returns the map f_n'*f_(n-1)* ... * f_2'*f_1.
##
InstallMethod ( IntersectionOfSubmodules, 
    "for two IsPathAlgebraMatModuleHomomorphisms",
    true,
    [ IsDenseList ], 
    0,
    function( list )

    local A, f, i;
#
#   Checking the input
#
    if IsEmpty( list ) then 
        Error( "<list> must be non-empty" ); 
    fi;
    A := RightActingAlgebra(Source(list[1]));
    for f in list do
        if RightActingAlgebra( Source(f) ) <> A then
            Error( "all entries in <list> must be homomorphisms of modules over the same algebra");
        fi;
    od;
    if not ForAll(list, IsInjective) then 
        Error("not all the arguments are monomorphisms,");
    fi;
    A := Range(list[1]);
    if not ForAll(list, x -> Range(x) = A) then 
        Error("must have submodules of the same module,");
    fi;
#
#   Doing the computations.
#
    f := list[1];
    for i in [2..Length(list)] do
        f := IntersectionOfSubmodules(f,list[i])[3];
    od;

    return f;
end
);

#######################################################################
##
#O  SumOfSubmodules( <f>, <g> )
##
##                                              f              g
##  Given two submodules/subrepresentations  M ---> X  and  N ---> X  of 
##  a module  X  by two monomorphism  f  and  g, this function computes 
##  the sum of the images of  f  and  g  in the sense that it computes 
##  the submodule/subrepresentation  M + N  generated by the image of  f  
##  and the image of  g. It returns the maps M + N ---> X,  M ---> M + N,  
##  N ---> M + N, in that order.  
##
InstallOtherMethod ( SumOfSubmodules, 
    "for two IsPathAlgebraMatModuleHomomorphisms",
    true,
    [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 
    0,
    function( f, g )

    local sumMN, inclusions, projections, F;
#
#   Checking the input
#
    if not IsInjective(f) then 
        Error("the first argument is not a monomorphism,");
    fi;
    if not IsInjective(g) then 
        Error("the second argument is not a monomorphism,");
    fi;
    if not ( Range(f) = Range(g) ) then 
        Error("must have two submodules of the same module,");
    fi;
#
#   Doing the computations
#
    sumMN := DirectSumOfQPAModules([Source(f),Source(g)]);
    inclusions := DirectSumInclusions(sumMN);
    projections := DirectSumProjections(sumMN);

    F := ImageProjectionInclusion(projections[1]*f + projections[2]*g);

    return [F[2],inclusions[1]*F[1],inclusions[2]*F[1]];
end
);

#######################################################################
##
#O  SumOfSubmodules( <args> )
##
##                                            f_i             
##  Given submodules/subrepresentations  M_i -----> X  for i = 1,...,n 
##  of a module  X  by  n  monomorphism  f_i, this function computes 
##  the intersection of the images of  f  and  g  in the sense that it 
##  computes the pullbacks of the diagrams
##             f_1'
##         E_1 ----> N           E_2 ---------> N
##          |        |            |             |
##     f_2' |        | f_2   f_3' |             | f_3  and so on.
##          |        |            |             |
##          V  f_1   V            V   f_2'*f_1  V
##          M -----> X           E_1 ---------> X
##  
##  and it returns the map f_n'*f_(n-1)* ... * f_2'*f_1.
##
InstallMethod ( SumOfSubmodules, 
    "for two IsPathAlgebraMatModuleHomomorphisms",
    true,
    [ IsDenseList ], 
    0,
    function( list )

    local A, f, i;
#
#   Checking the input
#
    if IsEmpty( list ) then 
        Error( "<list> must be non-empty" ); 
    fi;
    A := RightActingAlgebra(Source(list[1]));
    for f in list do
        if RightActingAlgebra( Source(f) ) <> A then
            Error( "all entries in <list> must be homomorphisms of modules over the same algebra");
        fi;
    od;
    A := Range(list[1]);
    if not ForAll(list, x -> Range(x) = A) then 
        Error("must have submodules of the same module,");
    fi;   
    if not ForAll(list, IsInjective) then 
        Error("not all the arguments are monomorphisms,");
    fi;
#
#   Doing the computations
#
    f := list[1];
    for i in [2..Length(list)] do
        f := SumOfSubmodules(f,list[i])[1];
    od;

    return f;
end
);

#######################################################################
##
#O  HomFactoringThroughProjOverAlgebra( <M>, <N> )
##
##  Given two modules  M  and  N  over a finite dimensional quotient
##  of a path algebra, this function computes a basis of the homomorphisms
##  that factor through a projective module.
##
InstallMethod ( HomFactoringThroughProjOverAlgebra,
    "for two PathAlgebraMatModule's",
    true,
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
    0,
    function( M, N )

    local f, HomMP_N, homthroughproj, flathomthroughproj, K, V, B, dim_M,
          dim_N, homs, hom, i, j, matrices, dimension, d1, d2, matrix;
    #
    # Testing input
    #
    if RightActingAlgebra(M) <> RightActingAlgebra(N) then
        Error("the entered modules are not modules over the same algebra,");
    fi;
    #
    # if one module is zero, return the empty list
    #
    if Dimension(M) = 0 or Dimension(N) = 0 then
        return [];
    fi;
    #
    # The maps factoring through projectives are those factoring through
    # the projective cover of  N.
    #
    f := ProjectiveCover(N);
    HomMP_N := HomOverAlgebra(M,Source(f));
    homthroughproj := List(HomMP_N, x -> x*f);
    flathomthroughproj := List(homthroughproj, x -> Flat(x!.maps));
    flathomthroughproj := Filtered(flathomthroughproj, x -> x <> Zero(x));
    homthroughproj := [];
    if Length(flathomthroughproj) > 0 then
        K := LeftActingDomain(M);
        V := VectorSpace(K, flathomthroughproj);
        B := BasisVectors(Basis(V));
        dim_M := DimensionVector(M);
        dim_N := DimensionVector(N);
        homs := [];
        for hom in B do
            matrices := [];
            dimension := 0;
            for i in [1..Length(dim_M)] do
                if dim_M[i] = 0 then
                    d1 := 1;
                else
                    d1 := dim_M[i];
                fi;
                if dim_N[i] = 0 then
                    d2 := 1;
                else
                    d2 := dim_N[i];
                fi;
                matrix := [];
                for j in [1..d1] do
                    Add(matrix,hom{[dimension + 1..dimension + d2]});
                    dimension := dimension + d2;
                od;
                Add(matrices,matrix);
            od;
            Add(homs,matrices);
        od;
        homthroughproj := List(homs, x -> RightModuleHomOverAlgebra(M,N,x));
    fi;

    return homthroughproj;
end
);

#######################################################################
##
#O  EndModuloProjOverAlgebra( <M>, <N> )
##
##  Given a module  M  over a finite dimensional quotient of a path
##  algebra, this function computes the homomorphism from
##  End(M) ----> End(M)/P(M,M), where P(M,M) is the ideal in End(M)
##  generated by the morphisms factoring through projective modules.
##
InstallMethod ( EndModuloProjOverAlgebra,
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule ],
    0,
    function( M )

    local EndM, factorproj, I;

    if Dimension(M) = 0 then
        Error("zero module given as input for EndModuleProjOverAlgebra,");
    else
        EndM := EndOverAlgebra(M);
        factorproj := HomFactoringThroughProjOverAlgebra(M,M);
        factorproj := List(factorproj, x -> FromHomMMToEndM(x));
        I := Ideal(EndM,factorproj);

        return NaturalHomomorphismByIdeal(EndM,I);
    fi;
end
);

#######################################################################
##
#O  FromHomMMToEndM( <f> )
##
##  This function gives a translation from endomorphisms of a module  M  
##  to the corresponding enodomorphism represented in the algebra
##  EndOverAlgebra(M). 
##
InstallMethod ( FromHomMMToEndM, 
    "for an endomorphism to an element of EndOverAlgebra",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    local K, dim_vect, end_f, r, j;

    K := LeftActingDomain(Source(f));
    dim_vect := DimensionVector(Source(f));
    end_f := NullMat(Dimension(Source(f)),Dimension(Source(f)),K);
    r := 1; 
    for j in [1..Length(dim_vect)] do 
        if dim_vect[j] <> 0 then 
            end_f{[r..r+dim_vect[j]-1]}{[r..r+dim_vect[j]-1]} := f!.maps[j];
            r := r + dim_vect[j];
        fi;
    od; 

    return end_f;
end
);

#######################################################################
##
#O  FromEndMToHomMM( <M>, <mat> )
##
##  This function gives a translation from an element in 
##  EndOverAlgebra(<M>) to an endomorphism of the module  M. 
##
InstallMethod ( FromEndMToHomMM, 
    "for a subset of EndOverAlgebra to HomOverAlgebra",
    true,
    [ IsPathAlgebraMatModule, IsMatrix ],
    0,
    function( M, mat )

    local K, dim_vect, maps, i, r;

    K := LeftActingDomain(M); 
    dim_vect := DimensionVector(M);

    maps := [];
    r := 1;
    for i in [1..Length(dim_vect)] do
        if dim_vect[i] = 0 then 
            Add(maps,NullMat(1,1,K));
        else
            Add(maps,mat{[r..r+dim_vect[i]-1]}{[r..r+dim_vect[i]-1]});
            r := r + dim_vect[i];
        fi;
    od;

    return RightModuleHomOverAlgebra(M,M,maps);
end
);

#######################################################################
##
#P  IsRightMinimal( <f> )
##
##  This function returns true is the homomorphism  <f>  is right 
##  minimal. 
##
InstallMethod ( IsRightMinimal, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    local B, C, BB, mat, Ann_f, radEndB;

    B := Source(f);
    C := Range(f);
    BB := HomOverAlgebra(B,B);
    mat := List(BB, x -> x*f);
    mat := List(mat,x -> Flat(x!.maps));
    Ann_f := NullspaceMat(mat);
    Ann_f := List(Ann_f,x -> LinearCombination(BB,x));
    radEndB := RadicalOfAlgebra(EndOverAlgebra(B)); 
    Ann_f := List(Ann_f, x -> FromHomMMToEndM(x));

    if ForAll(Ann_f, x -> x in radEndB) then 
        return true;
    else
        return false;
    fi;
end
);

#######################################################################
##
#P  IsLeftMinimal( <f> )
##
##  This function returns true is the homomorphism  <f>  is left 
##  minimal. 
##
InstallMethod ( IsLeftMinimal, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    local A, B, BB, mat, Ann_f, radEndB;

    A := Source(f);
    B := Range(f);
    BB := HomOverAlgebra(B,B);
    mat := List(BB, x -> f*x );
    mat := List(mat,x -> Flat(x!.maps));
    Ann_f := NullspaceMat(mat);
    Ann_f := List(Ann_f,x -> LinearCombination(BB,x));
    radEndB := RadicalOfAlgebra(EndOverAlgebra(B)); 
    Ann_f := List(Ann_f, x -> FromHomMMToEndM(x));

    if ForAll(Ann_f, x -> x in radEndB) then 
        return true;
    else
        return false;
    fi;
end
);

#######################################################################
##
#A  RightInverseOfHomomorphism( <f> )
##
##  This function returns false if the homomorphism  <f>  is not a  
##  splittable monomorphism, otherwise it returns a splitting of the
##  split monomorphism  <f>. 
##
InstallMethod ( RightInverseOfHomomorphism, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    local B, C, CB, id_B, mat, flat_id_B, split_f;

    B := Source(f);
    C := Range(f);
    if IsInjective(f) then 
        CB := HomOverAlgebra(C,B);
        if Length(CB) = 0 then 
            return false;
        else 
            mat := [];
            mat := List(CB, x -> f*x);
            id_B := IdentityMapping(B); 
            mat := List(mat,x -> Flat(x!.maps));
            flat_id_B := Flat(id_B!.maps); 
            split_f := SolutionMat(mat,flat_id_B);
            
            if split_f <> fail then 
                split_f := LinearCombination(CB,split_f);
                SetIsSplitEpimorphism(split_f,true);
                SetIsSplitMonomorphism(f,true);
                SetLeftInverseOfHomomorphism(split_f,f); 
                return split_f;
            else
                SetIsSplitMonomorphism(f,false);
                return false;
            fi;
        fi;
    else
        return false;
    fi;
end
  );

#######################################################################
##
#A  LeftInverseOfHomomorphism( <f> )
##
##  This function returns false if the homomorphism  <f>  is not a  
##  splittable epimorphism, otherwise it returns a splitting of the
##  split epimorphism  <f>. 
##
InstallMethod ( LeftInverseOfHomomorphism, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    local B, C, CB, id_C, mat, flat_id_C, split_f;

    B := Source(f);
    C := Range(f);
    if IsSurjective(f) then 
        CB := HomOverAlgebra(C,B);
        if Length(CB) = 0 then 
            return false;
        else 
            mat := List(CB, x -> x*f );
            id_C := IdentityMapping(C); 
            mat := List(mat,x -> Flat(x!.maps));
            flat_id_C := Flat(id_C!.maps); 
            split_f := SolutionMat(mat,flat_id_C);
            
            if split_f <> fail then 
                split_f := LinearCombination(CB,split_f);
                SetIsSplitMonomorphism(split_f,true);
                SetIsSplitEpimorphism(f,true);
                SetRightInverseOfHomomorphism(split_f,f); 
                return split_f;
            else
                SetIsSplitEpimorphism(f,false);
                return false;
            fi;
        fi;
    else
        return false;
    fi;
end
  );

#######################################################################
##
#P  IsSplitMonomorphism( <f> )
##
##  This function returns false if the homomorphism  <f>  is not a  
##  splittable monomorphism, otherwise it returns a splitting of the
##  homomorphism  <f>. 
##
InstallMethod ( IsSplitMonomorphism, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    local B, C, CB, id_B, mat, flat_id_B, split_f;

    B := Source(f);
    C := Range(f);
    if IsInjective(f) then 
        CB := HomOverAlgebra(C,B);
        if Length(CB) = 0 then 
            return false;
        else 
            mat := [];
            mat := List(CB, x -> f*x);
            id_B := IdentityMapping(B); 
            mat := List(mat,x -> Flat(x!.maps));
            flat_id_B := Flat(id_B!.maps); 
            split_f := SolutionMat(mat,flat_id_B);
            
            if split_f <> fail then 
                split_f := LinearCombination(CB,split_f);
                SetLeftInverseOfHomomorphism(split_f,f);
                SetRightInverseOfHomomorphism(f,split_f);
                SetIsSplitEpimorphism(split_f,true);
                return true;
            else
                return false;
            fi;
        fi;
    else
        return false;
    fi;
end
  );

#######################################################################
##
#P  IsSplitEpimorphism( <f> )
##
##  This function returns false if the homomorphism  <f>  is not a  
##  splittable epimorphism, otherwise it returns true. 
##
InstallMethod ( IsSplitEpimorphism, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    local B, C, CB, id_C, mat, flat_id_C, split_f;

    B := Source(f);
    C := Range(f);
    if IsSurjective(f) then 
        CB := HomOverAlgebra(C,B);
        if Length(CB) = 0 then 
            return false;
        else 
            mat := List(CB, x -> x*f );
            id_C := IdentityMapping(C); 
            mat := List(mat,x -> Flat(x!.maps));
            flat_id_C := Flat(id_C!.maps); 
            split_f := SolutionMat(mat,flat_id_C);
            
            if split_f <> fail then 
                split_f := LinearCombination(CB,split_f);
                SetLeftInverseOfHomomorphism(f, split_f);
                SetRightInverseOfHomomorphism(split_f,f);
                SetIsSplitMonomorphism(split_f,true);
                return true;
            else
                return false;
            fi;
        fi;
    else
        return false;
    fi;
end
);

InstallMethod ( MoreRightMinimalVersion, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f ) 

    local B, g, A, HomBA, HomAA, n, i, j, gg, hh;

    B := Source(f);
    g := KernelInclusion(f);
    A := Source(g);
    HomBA := HomOverAlgebra(B,A);
    if Length(HomBA) = 0 then 
        return f;
    else 
        HomAA  := HomOverAlgebra(A,A);
        n := Maximum(Concatenation(DimensionVector(A),DimensionVector(B)));
        for i in [1..Length(HomBA)] do 
            for j in [1..Length(HomAA)] do
                gg := (g*HomBA[i]*HomAA[j])^n;
                if gg <> ZeroMapping(A,A) then
                    hh := (HomBA[i]*HomAA[j]*g)^n;
                    return [KernelInclusion(hh)*f,ImageInclusion(hh)*f];
                fi;
            od;
        od;
        return f;
    fi;
end
);

InstallMethod ( MoreLeftMinimalVersion, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f ) 

    local B, g, C, HomCB, HomCC, n, i, j, gg, hh, t;

    B := Range(f);
    g := CoKernelProjection(f);
    C := Range(g);
    HomCB := HomOverAlgebra(C,B);
    if Length(HomCB) = 0 then 
        return f;
    else 
        HomCC := HomOverAlgebra(C,C);
        n := Maximum(Concatenation(DimensionVector(B),DimensionVector(C)));
        for i in [1..Length(HomCC)] do 
            for j in [1..Length(HomCB)] do
                gg := (HomCC[i]*HomCB[j]*g)^n;
                if gg <> ZeroMapping(C,C) then
                    hh := (g*HomCC[i]*HomCB[j])^n;
#                    t := IsSplitMonomorphism(KernelInclusion(hh));
                    t := RightInverseOfHomomorphism(KernelInclusion(hh));
                    return [f*t,f*ImageProjection(hh)];
                fi;
            od;
        od;
        return f;
    fi;
end
);

#######################################################################
##
#A  RightMinimalVersion( <f> )
##
##  This function returns a right minimal version  f'  of the 
##  homomorphism  <f>  in addition to a list of modules  B such that 
##  Source(f') direct sum the modules on the list  B  is isomorphic to
##  Source(f).
##
InstallMethod ( RightMinimalVersion, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    local Bprime, g, L;

    Bprime := [];
    g := f;
    repeat
        L := MoreRightMinimalVersion(g);
        if L <> g then 
            g := L[1];
            Add(Bprime,Source(L[2]));
        fi;
    until 
        L = g;
    
    SetIsRightMinimal(g,true);

    return [g,Bprime];
end
);

#######################################################################
##
#A  LeftMinimalVersion( <f> )
##
##  This function returns a left minimal version  f'  of the 
##  homomorphism  <f>  in addition to a list of modules  B such that 
##  Range(f') direct sum the modules on the list  B  is isomorphic to
##  Range(f).
##
InstallMethod ( LeftMinimalVersion, 
    "for a PathAlgebraMatModuleMap",
    true,
    [ IsPathAlgebraMatModuleHomomorphism ],
    0,
    function( f )

    local Bprime, g, L;

    Bprime := [];
    g := f;
    repeat
        L := MoreLeftMinimalVersion(g);
        if L <> g then 
            g := L[1];
            Add(Bprime,Range(L[2]));
        fi;
    until 
        L = g;

    SetIsLeftMinimal(g,true);

    return [g,Bprime];
end
);

#######################################################################
##
#O  HomFromProjective(<m>,<M>)
##
##  This function checks if the element  <m>  is an element of  <M>  and
##  if the element  <m>  is uniform. If so, the function constructs the 
##  map from the indecomposable projective generated by the vertex of 
##  support to the module  <M>  given by sending the vertex to the 
##  element  <m>. 
##
InstallMethod ( HomFromProjective, 
    "for an element in a PathAlgebraMatModule and the PathAlgebraMatModule",
    IsElmsColls,
    [ IsRightAlgebraModuleElement, IsPathAlgebraMatModule ],
    0,
    function( m, M )

    local A, K, num_vert, support, pos, P, B, mats, zeros, i, f;
#
# Finding the algebra acting on M and the number of vertices.
#
    A := RightActingAlgebra(M);
    K := LeftActingDomain(M);
    num_vert := Length(m![1]![1]);
#   
# Finding the support of m.
#
    support := SupportModuleElement(m);
#
# If the element  m  is uniform, we proceed to find the map.
# 
    if Length(support) = 1 then 
#
# And we need the "number of" the vertex ...
#
        pos := VertexPosition(support[1]);
#
# And the basis for the correct projective module, and the proj. module itself
#
        B := BasisOfProjectives(A)[pos];
        P := ShallowCopy(IndecProjectiveModules(A)[pos]);
#
# Then we calculate the matrices for the homomorphism
#
        mats := List([1..num_vert],x -> List(B[x], y -> ExtRepOfObj(m^y)![1][x]));
        zeros := Filtered([1..num_vert], x -> DimensionVector(P)[x] = 0);
        for i in zeros do 
            if DimensionVector(M)[i] <> 0 then         
                mats[i] := NullMat(1,DimensionVector(M)[i],K);
            else 
                mats[i] := NullMat(1,1,K);
            fi;
        od;
#
# Construct the homomorphism
#
        f := RightModuleHomOverAlgebra(P,M,mats);
        return f;
    else
        return fail;
    fi;
end
);

######################################################
##
#O InverseOfIsomorphism( <f> )
##
## <f> is a map which is iso. 
##
## Output is the map f^(-1), or the inverse of f.
##
InstallMethod( InverseOfIsomorphism,
               [ IsPathAlgebraMatModuleHomomorphism ],
               function( f )
    local maps, invmaps, m;
    if not(IsIsomorphism(f)) then
        Error("entered map is not an isomorphism");
    fi;
    maps := f!.maps;
    invmaps := [];
    for m in maps do
        if IsZero(m) then
            Append(invmaps, [m]);
        else
            Append(invmaps, [Inverse(m)]);
        fi;
    od;

    return RightModuleHomOverAlgebra(Range(f), Source(f), invmaps);

end);

#######################################################################
##
#O  HomomorphismFromImages( <M>, <N>, <genImages> )
##
##  <M> and <N> are modules over the same algebra.  This operation
##  computes a homomorphism f: M ---> N such that the basis
##
##    B := BasisVectors( Basis( M ) );
##
##  is sent to the list <genImages>, that is, B[i] is sent to genImages[i].
##  If this does not give a well-defined homomorphism, the final 
##  call to RightModuleHomOverAlgebra will fail.
##
InstallMethod( HomomorphismFromImages,
[ IsPathAlgebraMatModule, IsPathAlgebraMatModule, IsList ],
    function( M, N, genImages )

    local B, K, dim_vec_N, dim_vec_M, map, i, j, mat, p; 

    K := LeftActingDomain(M);

    B := BasisVectors( Basis( M ) );

#
#  Checks that the element B[i] is sent to an element of N with support in 
#  the same vertex as B[i].
#
    for i in [1..Length(B)] do
        if ( SupportModuleElement( B[i] ) <> SupportModuleElement( genImages[i] ) ) and
           ( not IsZero( genImages[i] ) ) then
            Error("does not give a homomorphism -- wrong support,");
        fi;
    od;
    
#
#  Computing dimension vectors so that we can insert zero matrices of 
#  the right size. 
#
    dim_vec_M := DimensionVector(M);
    dim_vec_N := DimensionVector(N);
    map := [];
    j := 0;
    for i in [1..Length(dim_vec_M)] do
#
#  If M is zero in vertex  i, then insert a zero matrix of
#  the right size, do not use any of the lifting information.
# 
        if dim_vec_M[i] = 0 then
            if  dim_vec_N[i] = 0 then 
                Add(map,NullMat(1,1,K));
            else
                Add(map,NullMat(1,dim_vec_N[i],K));
            fi;
        else
#
#  If M is non-zero in vertex  i, then use genImages
#  to compute the right matrix for the map from vertex  i  to  i.
# 
            mat := [];
            for p in [1..dim_vec_M[i]] do
                j := j + 1;
                if not ( IsPathModuleElem(genImages[j]) ) then
                    if IsList( ExtRepOfObj( genImages[j] ) ) then
                        Add(mat, ExtRepOfObj( genImages[j] )[i]);
                    else
                        Add(mat, ExtRepOfObj(ExtRepOfObj(genImages[j]))[i]);
                    fi;
                else
                    Add(mat,ExtRepOfObj(ExtRepOfObj(genImages[j]))[i]);
                fi;
            od;
            Add(map,mat);
        fi;
    od;
    return RightModuleHomOverAlgebra(M, N, map);

end
);

######################################################
##
#O MultiplyListsOfMaps( <projections>, <matrix>, <inclusions> )
##
## <projections> is a list of m maps, <matrix> is a list of
## m lists of n maps, <inclusions> is a list of n maps. 
## Considering <projections> as a 1xm-matrix, <matrix> as an
## mxn-matrix and and <inclusions> as an nx1-matrix, the 
## matrix product is computed. Naturally, the maps in the 
## matrices must be composable.
##
## Output is a map (not a 1x1-matrix).
##
## Utility method not supposed to be here, but there seemed
## to be no existing GAP method for this.
##
InstallMethod( MultiplyListsOfMaps,
                    [ IsList, IsList, IsList ],
                    function( projections, matrix, inclusions )
    local sum, list, n, m, i, j;

    n := Length(projections);
    m := Length(inclusions);
    list := [1..m];
    sum := 0;

    for i in [1..m] do
       list[i] := projections*matrix[i];
    od;

    sum := list*inclusions;
    return sum;
end);

#######################################################################
##
#O  IsomorphismOfModules( <M>, <N> )
##
##  Given two modules  <M>  and  <N>  over a (quotient of a) path algebra
##  this function return an isomorphism from  <M>  to  <N>  is the two
##  modules are isomorphic, and false otherwise.
##
InstallMethod ( IsomorphismOfModules, 
    "for two PathAlgebraMatModules",
    true,
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 
    0,
    function( M, N )

    local K, dim_M, mats, i, maps, splittingmaps, fn, fsplit, gn, gsplit, 
          f_alpha, f_beta, g_alpha, g_beta, q, B, genImages, delta, 
          recursiveIOM;
    # 
    # If  M  and  N are not modules over the same algebra, 
    # return an error message.
    #
    if RightActingAlgebra(M) <> RightActingAlgebra(N) then 
        return false;
    fi;
    #
    # If the dimension vectors of  M  and  N  are different, they
    # are not isomorphic, so return  false.
    #
    if DimensionVector(M) <> DimensionVector(N) then
        return false;
    fi;
    #
    # If both  M  and  N  are the zero module, return the zero map.
    #
    if Dimension(M) = 0 and Dimension(N) = 0 then
        return ZeroMapping(M,N);
    fi;
    #
    # If  M  and  N  are the same representations (i.e. with the 
    # same vector spaces and defining matrices, then return the 
    # "identity homomorphism". 
    #
    if M = N then 
        K := LeftActingDomain(M);
        dim_M := DimensionVector(M);
        mats := [];
        for i in [1..Length(dim_M)] do
            if dim_M[i] = 0 then
                Add(mats,NullMat(1,1,K));
            else
                Add(mats,IdentityMat(dim_M[i],K));
            fi;
        od;
        
        return RightModuleHomOverAlgebra(M,N,mats);
    fi;
    #
    # Finds a homomorphism from  M  to  N  and corresponding homomorphism
    # from  N  to  M  identifying a common non-zero direct summand of  M  
    # and  N, if such exists. Otherwise it returns false. 
    #
    splittingmaps := function( M, N) 
        local K, HomMN, HomNM, mn, nm, n, m, i, j,  
              l, zero, f, fnm, nmf; 
            
        HomMN := HomOverAlgebra(M,N);
        HomNM := HomOverAlgebra(N,M);
        mn := Length(HomMN);
        nm := Length(HomNM);
      
        if mn = 0 or nm = 0 then 
            return false;
        fi;
      
        m := Maximum(DimensionVector(M));
        n := Maximum(DimensionVector(N));
        if n = m then
            l := n;
        else
            l := Minimum([n,m]) + 1;
        fi;
      
        zero := ZeroMapping(M,M);
      
        for j in [1..nm] do
            for i in [1..mn] do
                if l>1 then  # because hom^0*hom => error!       
                    f := (HomMN[i]*HomNM[j])^(l-1)*HomMN[i];
                else 
                    f := HomMN[i];
                fi;
                
                fnm := f*HomNM[j];
              
                if fnm <> zero then
                    nmf := HomNM[j]*f; 
                    return [fnm, nmf, HomMN[i]];
                fi;
            od;
        od;

        return false;
    end; 
    
    maps := splittingmaps(M,N);
    #
    # No common non-zero direct summand,  M  and  N  are not isomorphic,
    # return false.
    #
    if maps = false then
        return false;
    fi;
    #
    # M  and  N  has a common non-zero direct summand, find the splitting of
    # the inclusion of this common direct summand in  M  and in  N.
    #
    fn := ImageInclusion(maps[1]);
#    fsplit := IsSplitMonomorphism(fn);
    fsplit := RightInverseOfHomomorphism(fn);
    gn := ImageInclusion(maps[2]);
#    gsplit := IsSplitMonomorphism(gn);
#     gsplit := LeftInverseOfHomomorphism(gn);    
    gsplit := RightInverseOfHomomorphism(gn);    
    #
    # Find the idempotents corresponding to these direct summands.
    #
    f_alpha := ImageInclusion(fsplit*fn);
    f_beta := KernelInclusion(fsplit*fn);    
    g_alpha := ImageInclusion(gsplit*gn);
    g_beta := KernelInclusion(gsplit*gn);
    #
    # If this direct summand is all of  M  and  N, return the appropriate
    # isomorphism.
    #
    if Dimension(Source(f_beta)) = 0 and Dimension(Source(g_beta)) = 0 then
       return maps[1]*maps[3]; 
    fi;
    #
    # If this direct summand is not all of  M  and  N, recursively construct
    # a possible isomorphism between the complements of this common direct 
    # summand in  M  and in  N.
    #
    q := f_alpha*maps[1]*maps[3];
    
    B := BasisVectors(Basis(Source(q)));
    genImages := List(B, x -> PreImagesRepresentative(g_alpha, ImageElm(q,x)));
    
    delta := HomomorphismFromImages(Source(f_alpha), Source(g_alpha), genImages);
    
    recursiveIOM := IsomorphismOfModules(Source(f_beta),Source(g_beta));
    
    if recursiveIOM <> false then
#        return IsSplitMonomorphism(f_alpha) * (delta * g_alpha) + 
#               IsSplitMonomorphism(f_beta) * (recursiveIOM * g_beta);
        return RightInverseOfHomomorphism(f_alpha) * (delta * g_alpha) + 
               RightInverseOfHomomorphism(f_beta) * (recursiveIOM * g_beta);        
    else 
        return false;
    fi;
end
  ); # IsomorphismOfModules

#######################################################################
##
#O  EndOfModuleAsQuiverAlgebra( <M> )
##
##  The argument of this function is a PathAlgebraMatModule  <M>  over
##  a quiver algebra over a field  K. It checks if the endomorphism ring
##  of  <M>  is  K-elementary for some (finite) field K. If the endomorphism
##  ring of  <M>  is K-elementary, then it returns a list of three elements, 
##  (i)   the endomorphism ring of  <M>, 
##  (ii)  the adjacency matrix of the quiver of the endomorphism ring and 
##  (iii) the endomorphism ring as a quiver algebra.
## 
InstallMethod( EndOfModuleAsQuiverAlgebra, 
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M )

  local EndM, EndMAsQuiverAlg;

    EndM := EndOverAlgebra( M );
    EndMAsQuiverAlg := AlgebraAsQuiverAlgebra( EndM );
    
    if EndMAsQuiverAlg <> fail then 
      return [ EndM, AdjacencyMatrixOfQuiver( QuiverOfPathAlgebra( EndMAsQuiverAlg[ 1 ] ) ), EndMAsQuiverAlg[ 1 ] ];
    else
      return fail;
    fi;
end
  );

#######################################################################
##
#O  TraceOfModule( <M>, <N> )
##
##  This function computes trace of the module  <M>  in  <N> by doing 
##  the following: computes a basis for Hom_A(M,N) and then computes 
##  the sum of all the images of the elements in this basis. This is 
##  also a minimal right Fac(M)-approximation. 
##
InstallMethod( TraceOfModule,
    "for two PathAlgebraMatModules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
    
    function( M, N )
    local homMN, B, trace;
    #
    # Checking if the modules  <M>  and  <N>  are modules over the same algebra.
    #
    if RightActingAlgebra(M) <> RightActingAlgebra(N) then
        Error("the entered modules are not modules over the same algebra,\n");
    fi;
    #
    # If  <M>  or  <N>  are zero, then the trace is the zero module.
    #
    if M = ZeroModule(RightActingAlgebra(M)) or N = ZeroModule(RightActingAlgebra(N)) then
        return ZeroMapping(ZeroModule(RightActingAlgebra(N)),N);
    fi;
    #
    # Finding a basis for  Hom(M,N), and taking the sum of all the images of these
    # homomorphisms from  <M>  to  <N>.  This is the trace of  <M>  in  <N>. 
    #
    homMN := HomOverAlgebra(M,N);
    if Length(homMN) = 0 then
        return ZeroMapping(ZeroModule(RightActingAlgebra(N)),N);
    fi;
    B := BasisVectors(Basis(M));
    trace := Flat(List(B, b -> List(homMN, f -> ImageElm(f,b))));
    
    return SubRepresentationInclusion(N,trace);
end
  );

#######################################################################
##
#O  RejectOfModule( <N>, <M> )
##
##  This function computes the reject of the module  <M>  in  <N> by doing 
##  the following: computes a basis for Hom_A(N,M) and then computes 
##  the intersections of all the kernels of the elements in this basis. 
##
InstallMethod( RejectOfModule,
    "for two PathAlgebraMatModules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
    
    function( N, M )

    local  homNM;
    #
    # Checking if the modules  <M>  and  <N>  are modules over the same algebra.
    #
    if RightActingAlgebra(M) <> RightActingAlgebra(N) then
        Error("the entered modules are not modules over the same algebra,\n");
    fi;
    #
    # If  Hom( N, M )  is zero, then the reject is the identity homomorphism N ---> N.
    #
    homNM := HomOverAlgebra( N, M );
    if Length( homNM ) = 0 then 
        return IdentityMapping( N );
    fi;
    #
    # Using the basis we found for  Hom( N, M )  above, and taking the intersection of 
    # all the kernels of these homomorphisms in the basis.  This is the reject of  
    # <M>  in  <N>.
    #
    if Length( homNM ) = 1 then
       return KernelInclusion( homNM[ 1 ] );
    fi;
    if Length( homNM ) > 1 then
	return IntersectionOfSubmodules( List( homNM, f -> KernelInclusion( f ) ) );
    fi;
end
  );

#######################################################################
##
#O  AllSimpleSubmodulesOfModule( <M> )
##
##  Returns all the different simple submodules of a module given as
##  inclusions into the module  <M>. 
##
##  Algorithm used is based on one created by Didrik Fosse.
##
InstallMethod( AllSimpleSubmodulesOfModule, 
"for a representation",
[ IsPathAlgebraMatModule ], 
function( M )
    local   field,  simples,  socleinclusion,  soc,  blocks,  b,  s,  
            simp;
    
    field := LeftActingDomain( M );
    if not IsFinite( field ) then
        Error( "The module is not a module over an algebra over a finite field.\n" );
    fi;
    if IsZero( M ) or IsSimpleQPAModule( M ) then 
        return [ IdentityMapping( M ) ];
    fi;
    #
    # Finding all simple submodules of M. 
    #
    simples := [ ];
    socleinclusion := SocleOfModuleInclusion( M );
    soc := Source( socleinclusion ); 
    blocks := BlockSplittingIdempotents( soc ); 
    blocks := List( blocks, b -> ImageInclusion( b ) );
    for b in blocks do
        for s in Subspaces( Source( b ), 1 ) do 
            simp := SubRepresentationInclusion( Source( b ), BasisVectors( Basis( s ) ) );
            Add( simples, simp * b * socleinclusion );
        od;
    od;
    
    return simples; 
end
  );

#######################################################################
##
#O  AllSubmodulesOfModule( <M> )
##
##  Returns all the different submodules of a module given as inclusions
##  into the module  <M>. 
##  
##  Algorithm used is based on one created by Didrik Fosse.
##
InstallMethod( AllSubmodulesOfModule, 
"for a representation",
        [IsPathAlgebraMatModule], 
        function(M)
    local   field,  submodules,  length,  simples,  listofsubmodules,  
            previousstep,  dim,  dimsub,  Vspaces,  newsubmodules,  U,  
            MmodU,  allsimplesinU,  s,  newsubmodule,  V;
    
    field := LeftActingDomain( M );
    if not IsFinite( field ) then
        Error( "The module is not a module over an algebra over a finite field.\n" );
    fi;
    if IsZero( M ) then 
        return [ IdentityMapping( M ) ];
    fi;
    submodules := [ SubRepresentationInclusion( M, [ ] ), IdentityMapping( M ) ];
    if IsSimpleQPAModule( M ) then
        return submodules;
    fi;
    
    length := Dimension( M );
    #
    # Finding all simple submodules of M. 
    #
    simples := AllSimpleSubmodulesOfModule( M );
    listofsubmodules := [ ];
    listofsubmodules[ 1 ] := [ SubRepresentationInclusion( M, [ ] ) ];
    Add( listofsubmodules, simples );
    previousstep := simples; 
    dim := Dimension( M );
    dimsub := 1;
    while dimsub < dim do
        dimsub := dimsub + 1;
        Vspaces := [ ];
        newsubmodules := [ ];
        for U in previousstep do
            MmodU := CoKernelProjection( U );
            allsimplesinU := AllSimpleSubmodulesOfModule( Range( MmodU ) ); 
            for s in allsimplesinU do 
                newsubmodule := PullBack( MmodU, s )[ 2 ];
                V := Subspace( UnderlyingLeftModule( M ), 
                             List( BasisVectors( Basis( Source( newsubmodule ) ) ), b -> ExtRepOfObj( ImageElm( newsubmodule, b ) ) ) );
                if Length( newsubmodules ) = 0 then
                    Add( newsubmodules, newsubmodule );
                    Add( Vspaces, V );
                else 
                    if not ( V in Vspaces ) then
                        Add( newsubmodules, newsubmodule );
                        Add( Vspaces, V ); 
                    fi;
                fi;
            od;
        od;
        Add( listofsubmodules, newsubmodules );
        previousstep := newsubmodules;
    od; 
    
    return listofsubmodules;
end
  );

#######################################################################
##
#O  RestrictionViaAlgebraHomomorphismMap( < f, h > )
##
##  Given an algebra homomorphism  f : A ---> B and a homomorphism
##  h : M ---> N of modules over B, this function returns the induced
##  homomorphism h : M ---> N  as a homomorphism of modules over  A.
##  
InstallMethod( RestrictionViaAlgebraHomomorphismMap, 
    "for a IsAlgebraHomomorphism and a IsPathAlgebraMatModuleHomomorphism",
    [ IsAlgebraHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 0,
    function( f, h )

    local   M,  N,  K,  A,  B,  vertices,  VM,  BVM,  VN,  BVN,  mats,  
            mat,  i,  resM,  resN;

    M := Source( h );
    N := Range( h );
    K := LeftActingDomain( M );
    A := Source( f );
    B := Range( f );
    if RightActingAlgebra( M ) <> B then 
        Error( "The entered homomorphism is not over the range of the algebra homomorphism.\n" );
    fi;
    vertices := One( A ) * VerticesOfQuiver( QuiverOfPathAlgebra( A ) );
    VM := List( vertices, v -> List( BasisVectors( Basis( M ) ), m -> m ^ ImageElm( f, v ) ) ); 
    VM := List( VM, W -> Filtered( W, w -> w <> Zero( w ) ) );
    BVM := List( VM, function( W ) if IsEmpty(W) then return []; else return Basis( VectorSpace( K, W ) ); fi; end );
    VN := List( vertices, v -> List( BasisVectors( Basis( N ) ), m -> m ^ ImageElm( f, v ) ) ); 
    VN := List( VN, W -> Filtered( W, w -> w <> Zero( w ) ) );
    BVN := List( VN, function( W ) if IsEmpty(W) then return []; else return Basis( VectorSpace( K, W ) ); fi; end );
    mats := [ ];
    for i in [ 1..Length( vertices ) ] do
    	mat := List( BasisVectors( BVM[ i ] ), b -> Coefficients( BVN[ i ], ImageElm( h, b ) ) );	
	Add( mats, mat );
    od;	
    resM := RestrictionViaAlgebraHomomorphism( f, M );
    resN := RestrictionViaAlgebraHomomorphism( f, N );
    
    return RightModuleHomOverAlgebra( resM, resN, mats ); 
end
  );

#######################################################################
##
#O  AllModulesOfLengthPlusOne( <M> )
##
##  Given a list of all modules of length <n> for an algebra, this
##  function constructs all modules of length <n + 1>
##  
InstallMethod( AllModulesOfLengthPlusOne, 
"for a representation",
[ IsList ], 
function( list )
    local   A,  K,  d,  simples,  firstsyzygysimples,  extensions,  S,  
            temp,  M,  ext,  modules,  t,  i,  u,  j,  s,  coeff,  h,  
            firstsyzygylist,  nonisomodules,  noniso;
    
    if Length( list ) = 0 then 
        Error( "List entered into AllModulesOfLengthPlusOne is empty.\n" );
    fi;
    A := RightActingAlgebra( list[ 1 ] );
    K := LeftActingDomain( A );
    if not IsFinite( K ) then
        Error( "The algebra is not over a finite field.\n" );
    fi;
    if not ForAll( list, m -> RightActingAlgebra( m ) = A ) then
        Error( "List entered into AllModulesOfLengthPlusOne doesn't consist of modules over the same algebra.\n" );
    fi;    
    d := Dimension( list[ 1 ] );
    if not ForAll( list, m -> Dimension( m ) = d ) then
        Error( "List entered into AllModulesOfLengthPlusOne doesn't have the same length.\n" );
    fi;
    
    simples := SimpleModules( A );
    firstsyzygysimples := List( simples, s -> KernelInclusion( ProjectiveCover( s ) ) );
    
    extensions := [ ];
    for S in simples do
        temp := [ ];
        for M in list do
            ext := ExtOverAlgebra( S, M );
            Add( temp, ext[ 2 ] );
        od;
        Add( extensions, temp );
    od;
    
    modules := [ ];
    t := Length( extensions );  # number of simples
    for i in [ 1..t ] do
        u := Length( extensions[ i ] ); # number of elements in <list>
        for j in [ 1..u ] do
            s := Length( extensions[ i, j ] ); # dimension of Ext-space
            if s = 1 then
                Add( modules, DirectSumOfQPAModules( [ list[ u ], simples[ i ] ] ) );
                Add( modules, Range( PushOut( extensions[ i, j ][ 1 ], firstsyzygysimples[ i ] )[ 2 ] ) );
            elif s > 1 then
                for coeff in K^s do
                    h := LinearCombination( extensions[ i, j ], coeff );
                    Add( modules, Range( PushOut( h, firstsyzygysimples[ i ] )[ 2 ] ) );
                od;
            fi;
        od;
    od;
    modules := Unique( modules );

    firstsyzygylist := List( list, m -> KernelInclusion( ProjectiveCover( m ) ) );
    
    extensions := [ ];
    for M in list do
        temp := [ ];
        for S in simples do
            ext := ExtOverAlgebra( M, S );
            Add( temp, ext[ 2 ] );
        od;
        Add( extensions, temp );
    od;

    t := Length( extensions );  # number of elements in <list>
    for i in [ 1..t ] do
        u := Length( extensions[ i ] ); # number of simples
        for j in [ 1..u ] do
            s := Length( extensions[ i, j ] ); # dimension of Ext-space
            if s = 1 then
                Add( modules, DirectSumOfQPAModules( [ list[ i ], simples[ j ] ] ) );
                Add( modules, Range( PushOut( extensions[ i, j ][ 1 ], firstsyzygylist[ i ] )[ 2 ] ) );
            elif s > 1 then
                for coeff in K^s do
                    h := LinearCombination( extensions[ i, j ], coeff );
                    Add( modules, Range( PushOut( h, firstsyzygylist[ i ] )[ 2 ] ) );
                od;
            fi;
        od;
    od;
    modules := Unique( modules );

    s := Length( modules ); 
    nonisomodules := [ ];
    Add( nonisomodules, modules[ 1 ] );
    noniso := 1;
    for i in [ 2..s ] do
        if ForAll( [ 1..noniso ], j -> not IsomorphicModules( nonisomodules[ j ], modules[ i ] ) ) then
            Add( nonisomodules, modules[ i ] );
            noniso := noniso + 1;
        fi;
    od;
    
    return nonisomodules; 
end
  );

#######################################################################
##
#O  AllModulesOfLengthAtMost( <A>, n )
##
##  Returns all the different modules over an algebra of length at most
##  <n> over the algebra <A>.
## 
InstallMethod( AllModulesOfLengthAtMost, 
"for an algebra and an integer",
[ IsQuiverAlgebra, IS_INT ], 
function( A, n )
    local   field,  simples,  allmodules,  previousstep,  i;
    
    field := LeftActingDomain( A );
    if not IsFinite( field ) then
        Error( "The algebra is not over a finite field.\n" );
    fi;
    if n < 0 then
        return [ ];
    fi;
    if n = 0 then
        return [ ZeroModule( A ) ];
    fi;
    simples := SimpleModules( A );
    allmodules := Concatenation( [ ZeroModule( A ) ], simples );
    if n = 1 then 
        return allmodules;
    fi;
    
    previousstep := simples;
    for i in [ 2..n ] do
        previousstep := AllModulesOfLengthPlusOne( previousstep );
        allmodules := Concatenation( allmodules, previousstep );
    od;

    return allmodules;
end
  );

#######################################################################
##
#O  AllIndecModulesOfLengthAtMost( <A>, n )
##
##  Returns all the different indecomposable modules over an algebra 
##  of length at most <n> over the algebra <A>.
## 
InstallMethod( AllIndecModulesOfLengthAtMost, 
"for an algebra and an integer",
[ IsQuiverAlgebra, IS_INT ], 
function( A, n )
    local   allmodules;
    
    allmodules := AllModulesOfLengthAtMost( A, n );
    
    return Filtered( allmodules{ [ 2..Length( allmodules ) ] }, IsIndecomposableModule );
     
end
  );

#######################################################################
##
#A  FromIdentityToDoubleStarHomomorphism( <M> )
##
##  Returns the homomorphism from  M  to  M^{**}  for module  <M>.
##
InstallMethod( FromIdentityToDoubleStarHomomorphism, 
    "for a path algebra module",
    [ IsPathAlgebraMatModule ],
    function( M )
    
    local   A,  Aop,  K,  P,  num_vert,  V,  BV,  BasisVflat,  b,  
            Mstar,  dim_vect_star,  Pop,  Vop,  BVop,  BasisVflatop,  
            Mdstar,  dim_vect_dstar,  BM,  dimM,  matrices,  basisPop,  
            fam,  i,  temp,  start,  j,  m,  mats,  l,  mstar,  h,  
            hflat,  startpos,  endpos;
    # 
    # Setting up necessary infrastructure.
    # 
    A := RightActingAlgebra( M );
    Aop := OppositeAlgebra( A );
    K := LeftActingDomain( M ); 
    P := IndecProjectiveModules( A );
    num_vert := Length( P );
    #
    # Finding the vector spaces in each vertex in the representation of 
    # Hom_A(M,A)  and making a vector space of the flatten matrices in 
    # each basis vector in  Hom_A(M,e_iA).
    #
    V := List( P, p -> HomOverAlgebra( M, p ) );
    BV := List( V, v -> List( v, x -> 
                  Flat( MatricesOfPathAlgebraMatModuleHomomorphism( x ) ) ) );
    BasisVflat := [];
    for b in BV do
        if Length( b ) <> 0 then 
            Add( BasisVflat, 
                Basis( Subspace( FullRowSpace( K, Length( b[ 1 ] ) ), b, "basis"), b ) ); 
        else 
            Add( BasisVflat, Basis( TrivialSubspace( K^1 ) ) );
        fi;
    od;
    #
    # Computing M^*.
    #
    Mstar := StarOfModule( M );
    dim_vect_star := DimensionVector( Mstar );
    #
    # Finding the vector spaces in each vertex in the representation of 
    # Hom_Aop(M^*,Aop)  and making a vector space of the flatten matrices in 
    # each basis vector in  Hom_Aop(M^*,e^op_iAop).
    #
    Pop := IndecProjectiveModules( Aop );
    Vop := List( Pop, p -> HomOverAlgebra( Mstar, p ) );
    BVop := List( Vop, v -> List( v, x -> 
                  Flat( MatricesOfPathAlgebraMatModuleHomomorphism( x ) ) ) );
    BasisVflatop := [];
    for b in BVop do
        if Length( b ) <> 0 then 
            Add( BasisVflatop, 
                Basis( Subspace( FullRowSpace( K, Length( b[ 1 ] ) ), b, "basis"), b ) ); 
        else 
            Add( BasisVflatop, Basis( TrivialSubspace( K^1 ) ) );
        fi;
    od;
    #
    # Computing M^{**}
    #    
    Mdstar := StarOfModule( Mstar ); 
    dim_vect_dstar := DimensionVector( Mdstar );
    #
    # Computing the morphism from  M ----> M^**
    #    
    BM := BasisVectors( Basis( M ) );
    dimM := DimensionVector( M );
    matrices := [ ];
    basisPop := List( BasisOfProjectives( Aop ), b -> Flat( b ) );
    basisPop := List( basisPop, b -> Basis( Subspace( Aop, b, "basis" ) ) );
    fam := ElementsFamily( FamilyObj( A ) );
    for i in [ 1..num_vert ] do
        if dimM[ i ] = 0 then
            if dim_vect_dstar[ i ] = 0 then 
                Add( matrices, NullMat( 1, 1, K ) );
            else
                Add( matrices, NullMat( 1, dim_vect_dstar[ i ], K ) );
            fi;
        else
            if dim_vect_dstar[ i ] = 0 then
                Add( matrices, NullMat( dimM[ i ], 1, K ) );
            else
            #
            # For each m in M[ i ] construct the map from  M^* ---> A given by  g |---> g( m )
            #
                temp := [ ];
                start := Sum( dimM{ [ 1..i - 1 ] } );
                for j in [ 1..dimM[ i ] ] do
                    m := BM[ start + j ];
                    mats := [ ];
                    for l in [ 1..num_vert ] do
                        if Length( V[ l ] ) = 0 then
                            if DimensionVector( Pop[ i ] )[ l ] = 0 then
                                Add( mats, NullMat( 1, 1, K ) );
                            else
                                Add( mats, NullMat( 1, DimensionVector( Pop[ i ] )[ l ], K ) );
                            fi;
                        else
                            if DimensionVector( Pop[ i ] )[ l ] = 0 then
                                Add( mats, NullMat( Length( V[ l ] ), 1, K ) );
                            else
                                mstar := List( V[ l ], g -> ImageElm( g, m ) );
				if IsPathAlgebra( A ) then
				    mstar := List( mstar, w -> ElementInIndecProjective( A, w, l ) ); 
				else
				    mstar := List( mstar, w -> ElementOfQuotientOfPathAlgebra( fam, ElementInIndecProjective( A, w, l ), false ) );
				fi;
                                mstar := List( mstar, w -> OppositePathAlgebraElement( w ) );
                                mstar := List( mstar, w -> Coefficients( basisPop[ i ], w ) );
                                mstar := List( mstar, w -> LinearCombination( Basis( Pop[ i ] ), w ) );
                                mstar := List( mstar, w -> ExtRepOfObj( ExtRepOfObj( w ) )[ l ]  );
                                Add( mats, mstar );
                            fi;
                        fi;
                    od;
                    h := RightModuleHomOverAlgebra( Mstar, Pop[ i ], mats );
                    hflat := Flat( MatricesOfPathAlgebraMatModuleHomomorphism( h ) );
                    mstar := Coefficients( BasisVflatop[ i ], hflat );
                    startpos := Sum( dim_vect_dstar{ [ 1..i - 1 ] } ) + 1;
                    endpos := startpos + dim_vect_dstar[ i ] - 1;
                    mstar := LinearCombination( BasisVectors( Basis( Mdstar ) ){ [ startpos..endpos ] }, mstar );
                    mstar := ExtRepOfObj( ExtRepOfObj( mstar ) )[ i ];
                    Add( temp, mstar );
                od;
                Add( matrices, temp );
            fi;
        fi;
    od;
    
    return RightModuleHomOverAlgebra( M, Mdstar, matrices ); 
end
);

#######################################################################
##
#O  MatrixOfHomomorphismBetweenProjectives( <f> )
##
##  Given a homomorphism  <f> from <P = \oplus v_iA> to <P' = \oplus w_iA>
##  where <v_i> and <w_i> are vertices, this function finds the 
##  homomorphism as a matrix in <\oplus v_iAw_i>. 
##  
InstallMethod( MatrixOfHomomorphismBetweenProjectives, 
  "for a homomorphism between projectives",
  true,
  [ IsPathAlgebraMatModuleHomomorphism ], 0,
  function( f )
    
    local   zero,  projections,  inclusions,  A,  imageofgenerators,  
            l,  indexofprojectives;
    
    zero := Zero( OriginalPathAlgebra( RightActingAlgebra( Source( f ) ) ) );
    if IsZero( Range( f ) ) then 
        return [ [ zero ] ];
    fi;
    projections := DirectSumProjections( Range( f ) );
    if IsZero( Source( f ) ) then
        return [ List( [ 1..Length( projections ) ], i -> zero ) ];
    fi;
    inclusions := DirectSumInclusions( Source( f ) );
    if not IsDirectSumOfModules( Source( f ) ) or not IsDirectSumOfModules( Range( f ) ) then
        Error( "The enter projectives are not a direct sum indecomposable projectives.\n" );
    fi;
    
    A := RightActingAlgebra( Source( f ) );
    imageofgenerators := List( inclusions, x -> ImageElm( x,  MinimalGeneratingSetOfModule( Source( x ) )[ 1 ] ) );
    imageofgenerators := List( imageofgenerators, x -> ImageElm( f, x ) );
    l := Length( projections );
    imageofgenerators := List( imageofgenerators, x -> List( projections, p -> ImageElm( p, x ) ) );
    indexofprojectives := List( projections, x -> SupportModuleElement( MinimalGeneratingSetOfModule( Range( x ) )[ 1 ] )[ 
1 ] );
    indexofprojectives := List( indexofprojectives, x -> VertexPosition( x ) ); 
    imageofgenerators := List( imageofgenerators, x -> List( [ 1..l ], 
                                 y -> ElementInIndecProjective( A, x[ y ], indexofprojectives[ y ] ) ) ); 
    
    return TransposedMat( imageofgenerators );
end
  );

InstallMethod( FromMatrixToHomomorphismOfProjectives, 
    "for a matrix of elements in a quotient of a path algebra, and two list of vertices",
    true,
    [ IsQuiverAlgebra, IsMatrix, IsHomogeneousList, IsHomogeneousList ],
    0,
    function( A, mat, vert1, vert0 )               

  local m, P0vert, supprows, rowcount, row, tempsupp, r, test, P1vert, 
        suppcolumns, P1, P1projections, P0, P0inclusions, newmat, P, 
        dimmat, i, temp, j, fam, vfam, elem;
    
    for m in mat do
      if not ForAll( m, a -> a in A ) then
        Error( "The entered matrix is not a homogeneous matrix in one and the same algebra.\n" );
      fi;
    od;
    
    P0vert := [ ];
    supprows := List( mat, row -> List( row, r -> LeftSupportOfQuiverAlgebraElement( A, r ) ) );
    rowcount := 0;
    for row in supprows do
      rowcount := rowcount + 1;
      tempsupp := [];
      for r in row do
        if r <> [] then
          if Length( r ) > 1 then
            Error( "The entered matrix is not given by uniform elements.\n" );
          else
            Add( tempsupp, r[ 1 ] );
          fi;
        fi;
      od;
      if Length( tempsupp ) > 0 then 
        test := tempsupp[ 1 ];
        if not ForAll( tempsupp, t -> t = test ) then
          Error( "The entered matrix is not row homogeneous.\n" );
        fi;
      fi;
      if tempsupp = [] then 
        Add( P0vert, vert0[ rowcount ] );
      else
        Add( P0vert, test );
      fi;
    od;
    
    P1vert := [ ];
    suppcolumns := List( TransposedMat( mat ), row -> List( row, r -> RightSupportOfQuiverAlgebraElement( A, r ) ) );
    rowcount := 0;
    for row in suppcolumns do
      rowcount := rowcount + 1;
      tempsupp := [];
      for r in row do
        if r <> [] then
          if Length( r ) > 1 then
            Error( "The entered matrix is not given by uniform elements.\n" );
          else
            Add( tempsupp, r[ 1 ] );
          fi;
        fi;
      od;
      if Length( tempsupp ) > 0 and not ForAll( tempsupp, t -> t = tempsupp[ 1 ] ) then
        Error( "The entered matrix is not column homogeneous.\n" );
      fi;
      if tempsupp = [] then
        Add( P1vert, vert1[ rowcount ] );
      else
        Add( P1vert, tempsupp[ 1 ] );
      fi;
    od;

    if ( P1vert <> vert1 ) or ( P0vert <> vert0 ) then
      Error( "Discrepancy between entered list of vertices and support in the entered matrix.\n" );
    fi;
    
    P1 := DirectSumOfQPAModules( List( vert1, v -> IndecProjectiveModules( A )[ v ] ) );
    P1projections := DirectSumProjections( P1 );
    P0 := DirectSumOfQPAModules( List( vert0, v -> IndecProjectiveModules( A )[ v ] ) );
    P0inclusions := DirectSumInclusions( P0 );
    
    newmat := [];
    P := IndecProjectiveModules( A );
    dimmat := DimensionsMat( mat );
    for i in [ 1..dimmat[ 1 ] ] do
      temp := [];
      for j in [ 1..dimmat[ 2 ] ] do
        r := mat[ i ][ j ];
        if IsZero( r ) then
          Add( temp, ZeroMapping( Range( P1projections[ j ] ), Source( P0inclusions[ i ] ) ) );
        else
          fam := FamilyObj( Zero( Source( P0inclusions[ i ] ) )![ 1 ] );
          vfam := FamilyObj( Zero( Source( P0inclusions[ i ] ) ) );
          elem := ElementIn_vA_AsElementInIndecProj( A, r );
          elem := ExtRepOfObj( ExtRepOfObj( elem ) );
          elem := ObjByExtRep(vfam, PathModuleElem( fam, elem ) );
          Add( temp, HomFromProjective( elem, Source( P0inclusions[ i ] ) ) );
        fi;
      od;
      Add( newmat, temp );
    od;
    
    for i in [ 1..dimmat[ 1 ] ] do
      for j in [ 1..dimmat[ 2 ] ] do
        newmat[ i ][ j ] := P1projections[ j ] * newmat[ i ][ j ] * P0inclusions[ i ];
      od;
    od;
    
    return Sum( Flat( newmat ) );
end
  );
