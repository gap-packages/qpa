# GAP Implementation
# $Id: specialreps.gi,v 1.12 2012/05/22 14:54:46 sunnyquiver Exp $

#######################################################################
##
#A  IndecProjectiveModules( <A> )
##
##  This function constructs all the indecomposable projective modules
##  over a finite dimensional path algebra or quotient of path algebra. 
##  It returns all the indecomposable projective modules as a list of 
##  modules corresponding to the numbering of the vertices. 
##
InstallMethod ( IndecProjectiveModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A ) 
    local   Q,  num_vert,  fam,  K,  num_arrows,  vertices,  
            arrows_as_paths,  mat_list,  list_of_min_gen,  
            dimension_vector_list,  p,  P,  B,  length_B,  
            intervals_of_basis,  i,  j,  zero_vertices,  vertices_Q,  
            mat,  a,  partial_mat,  mat_index,  source,  target,  
            source_index,  target_index,  rows,  cols,  arrow,  
            vector,  indec_proj_list;
    #
    #    Testing input
    #   
    if not IsFiniteDimensional(A) then
        Error("the entered algebra is not finite dimensional,\n");
    fi;
    if not IsPathAlgebra(A) and not IsAdmissibleQuotientOfPathAlgebra(A) then 
        TryNextMethod();
    fi;
    Q := QuiverOfPathAlgebra(A);
    num_vert := Length(VerticesOfQuiver(Q));
    fam := ElementsFamily(FamilyObj(A));
    #
    #    Finding vertices and arrows as elements of the algebra  A  for later use. 
    #
    K := LeftActingDomain(A);
    num_arrows := Length(ArrowsOfQuiver(Q)); 
    vertices := List(VerticesOfQuiver(Q), v -> v*One(A));
    arrows_as_paths := List(ArrowsOfQuiver(Q), a -> a*One(A));
    # 
    #
    #
    mat_list := [];
    list_of_min_gen := [];
    dimension_vector_list := [ ];
    for p in [1..num_vert] do
        #
        #   Finding a K-basis for indecomposable projective associated with vertex  p.
        #
        P := RightIdeal(A,[vertices[p]]);
        B := CanonicalBasis(P);
        length_B := Length(B);
        intervals_of_basis := List([1..num_vert], x -> []);
        for i in [1..Length(B)] do
            for j in [1..num_vert] do
                if ( B[i]*vertices[j] <> Zero(A) ) then
                    Add(intervals_of_basis[j],i);
                fi;
            od;
        od;
        Add( dimension_vector_list, List( intervals_of_basis, Length ) );
        #  
        #   Finding where the indecomposable projective has no support.
        #
        zero_vertices := [];
        for i in [1..num_vert] do
            if ( Length(intervals_of_basis[i]) = 0 ) then
                Add(zero_vertices,i);
            fi;
        od;
        vertices_Q := VerticesOfQuiver(Q); 
        # 
        #   Finding the matrices defining the representation
        # 
        mat := [];
        for a in arrows_as_paths do
            partial_mat := [];
            mat_index := Position(arrows_as_paths,a);
            source := SourceOfPath(TipMonomial(a));
            target := TargetOfPath(TipMonomial(a));
            source_index := Position(vertices_Q,source);
            target_index := Position(vertices_Q,target);
            rows := Length(intervals_of_basis[source_index]);
            cols := Length(intervals_of_basis[target_index]);
            arrow := TipMonomial(a);
            if ( rows <> 0 and cols <> 0 ) then
                for i in intervals_of_basis[source_index] do
                    vector := Coefficients(Basis(P), Basis(P)[i]*a){intervals_of_basis[target_index]};
                    Add(partial_mat,vector);
                od;
                Add(mat,[ String( arrow ),partial_mat]);
            fi;
        od;
        Add(mat_list,mat);
        #
        #   Constructing a generator for each indecomposable projective
        #
        for i in [1..Length(intervals_of_basis)] do
            if intervals_of_basis[i] = [] then
                intervals_of_basis[i] := [Zero(K)];
            else 
                for j in [1..Length(intervals_of_basis[i])] do
                    if intervals_of_basis[i][j] > 1 then 
                        intervals_of_basis[i][j] := Zero(K);
                    else
                        intervals_of_basis[i][j] := One(K);
                    fi;
                od;
            fi;
        od;
        Add(list_of_min_gen,intervals_of_basis);
    od;
    #
    #   Construct the indecomposable projectives and set the generating set
    # 
    indec_proj_list := [];
    for i in [1..num_vert] do    
        Add(indec_proj_list,RightModuleOverPathAlgebra( A, dimension_vector_list[ i ], mat_list[ i ] ) );
        list_of_min_gen[i] := PathModuleElem(FamilyObj(Zero(indec_proj_list[i])![1]),list_of_min_gen[i]); 
        list_of_min_gen[i] := Objectify( TypeObj( Zero(indec_proj_list[i]) ), [ list_of_min_gen[i] ] );
        SetMinimalGeneratingSetOfModule(indec_proj_list[i],[list_of_min_gen[i]]);
        SetIsIndecomposableModule(indec_proj_list[i], true);
        SetIsProjectiveModule(indec_proj_list[i], true);
    od;
    
    return indec_proj_list;
end
);


#######################################################################
##
#A  IndecInjectiveModules( <A> )
##
##  This function constructs all the indecomposable injective modules
##  over a finite dimensional path algebra or quotient of path algebra. 
##  It returns all the indecomposable injective modules as a list of 
##  modules corresponding to the numbering of the vertices. 
##
InstallMethod ( IndecInjectiveModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A )
        local A_op, P_op; 

        A_op := OppositeAlgebra(A);
        P_op := IndecProjectiveModules(A_op);
        
        return List(P_op, x -> DualOfModule(x));
    end
);
    
#######################################################################
##
#A  SimpleModules( <A> )
##
##  This function constructs all the simple modules over a finite 
##  dimensional path algebra or quotient of path algebra. It returns 
##  all the simple modules as a list of modules corresponding to the 
##  numbering of the vertices. 
##
InstallMethod ( SimpleModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A )

    local KQ, num_vert, simple_rep, zero, v, temp, s;
#
    if not IsFiniteDimensional(A) then
        Error("argument entered is not a finite dimensional algebra,\n");
    fi;
    if ( not IsPathAlgebra(A) ) and ( not IsAdmissibleQuotientOfPathAlgebra(A) ) then
        Error("argument entered is not a quotient of a path algebra by an admissible ideal,\n");
    fi;
    KQ := OriginalPathAlgebra(A); 
    num_vert := NumberOfVertices(QuiverOfPathAlgebra(A));
    simple_rep := [];
    zero := List([1..num_vert], x -> 0);
    for v in [1..num_vert] do
    	temp := ShallowCopy(zero);
        temp[v] := 1; 
        Add(simple_rep,RightModuleOverPathAlgebra(A,temp,[]));
    od;
    for s in simple_rep do
        SetIsSimpleQPAModule(s, true);
    od;
    return simple_rep;
end
);

#######################################################################
##
#A  ZeroModule( <A> )
##
##  This function constructs the zero module over a path algebra or 
##  quotient of path algebra. 
##
InstallMethod ( ZeroModule, 
    "for a finite dimensional quotient of a path algebra",
    [ IsQuiverAlgebra ], 0,
    function( A )

    local KQ, num_vert, zero;

    KQ := OriginalPathAlgebra(A); 
    num_vert := NumberOfVertices(QuiverOfPathAlgebra(A));
    zero := List([1..num_vert], x -> 0);

    return RightModuleOverPathAlgebra(A,zero,[]);
end
);

#######################################################################
##
#O  BasisOfProjectives(<A>)
##
##  Given a finite dimensional qoutient of a path algebra  <A>, this 
##  function finds in the basis of each indecomposable projective 
##  in terms of paths (nontips of the ideal  I defining  A).
##
InstallMethod ( BasisOfProjectives, 
    "for a finite dimensional quotient of a path algebra",
    [ IsQuotientOfPathAlgebra ], 0,
    function( A )

    local Q, num_vert, fam, vertices, basis_of_projs, 
          v, basis_of_proj, P, B, i, j;
#
#    Testing input if finite dimensional.
#       
   fam := ElementsFamily(FamilyObj(A));
   if HasGroebnerBasisOfIdeal(fam!.ideal) and 
          AdmitsFinitelyManyNontips(GroebnerBasisOfIdeal(fam!.ideal)) then 
#
      Q := QuiverOfPathAlgebra(OriginalPathAlgebra(A));
      num_vert := NumberOfVertices(Q); 
      vertices := List(VerticesOfQuiver(Q), x -> x*One(A));
      basis_of_projs := [];
      for v in vertices do
         basis_of_proj := List([1..num_vert], x -> []);
         P := RightIdeal(A,[v]);
         B := CanonicalBasis(P);
         for i in [1..Length(B)] do
            for j in [1..num_vert] do
               if ( B[i]*vertices[j] <> Zero(A) ) then
                   Add(basis_of_proj[j],B[i]);
               fi;
            od;
         od;
         Add(basis_of_projs,basis_of_proj);
      od;
      return basis_of_projs;
   else
      Print("Need to have a finite dimensional quotient of a path algebra as argument.\n");
      return fail;      
   fi;
end
);

#######################################################################
##
#O  BasisOfProjectives(<A>)
##
##  Given a finite dimensional path algebra  <A>, this function finds 
##  in the basis of each indecomposable projective in terms of paths.
##
InstallOtherMethod ( BasisOfProjectives, 
   "for a finite dimensional quotient of a path algebra",
   [ IsPathAlgebra ], 0,
   function( A )
   
   local Q, num_vert, fam, vertices, basis_of_projs, 
          v, basis_of_proj, P, B, i, j;
#
#  Testing input if finite dimensional.
#       
   Q := QuiverOfPathAlgebra(A); 
   num_vert := NumberOfVertices(Q);
   if not IsAcyclicQuiver(Q)  then
      Print("Need to have a finite dimensional path algebra as argument.\n");
      return fail;
   else
      vertices := List(VerticesOfQuiver(Q), x -> x*One(A));
      basis_of_projs := [];
      for v in vertices do
         basis_of_proj := List([1..num_vert], x -> []);
         P := RightIdeal(A,[v]);
         B := CanonicalBasis(P);
         for i in [1..Length(B)] do
            for j in [1..num_vert] do
               if ( B[i]*vertices[j] <> Zero(A) ) then
                  Add(basis_of_proj[j],B[i]);
               fi;
            od;
         od;
         Add(basis_of_projs,basis_of_proj);
      od;
   fi;
  
   return basis_of_projs;
end
);

#######################################################################
##
#O  ElementInIndecProjective( <A>, <m>, <s> )
##
##  Given an admissible quotient A of a path algebra, an element
##  m  in the s-th indecomposable projective module  v_sA  over A, this 
##  function computes the element  m  in the module  v_sA  as a preimage 
##  in the ideal  v_sR. 
##
InstallMethod ( ElementInIndecProjective, 
    "for a fin. dim. quotient of a path algebra, an element in an indec. projective module over this algebra and a positive integer",
    true,
    [ IsQuotientOfPathAlgebra, IsAlgebraModuleElement, IS_INT ], 
    0,
    function( A, m, s )

    local   B, i, elem;
        
    if not IsAdmissibleQuotientOfPathAlgebra( A ) then
        Error("The entered algebra is not a finite dimensional algebra,\n");    
    fi;
        
    B := ShallowCopy( BasisOfProjectives( A )[ s ] );
    for i in [ 1..Length( B ) ] do
        if Length( B[ i ] ) = 0 then 
            B[ i ] := Zero( OriginalPathAlgebra( A ) );
        else
            B[ i ] := List( B[ i ], x -> x![ 1 ] );
        fi;
    od;
    B := Flat( B );
    elem := Flat( ExtRepOfObj( ExtRepOfObj( m ) ) );
    
    return LinearCombination( B, elem );
end 
);

#######################################################################
##
#O  ElementInIndecProjective( <A>, <m>, <s> )
##
##  Given a finite dimension path algebra A, an element m  in the s-th 
##  indecomposable projective module  v_sA  over A, this function 
##  computes the element  m  in the module  v_sA  as a preimage in the 
##  ideal  v_sA. 
##
InstallOtherMethod ( ElementInIndecProjective, 
 "for a finite dimensional path algebra, an element in an indec. projective module over this algebra and a positive integer",
    true,
    [ IsPathAlgebra, IsAlgebraModuleElement, IS_INT ], 
    0,
    function( A, m, s )

    local   B, i, elem;
        
    if not IsAcyclicQuiver( QuiverOfPathAlgebra( A ) ) then
        Error("The entered algebra is not a finite dimensional algebra,\n");    
    fi;
        
    B := ShallowCopy( BasisOfProjectives( A )[ s ] );
    for i in [ 1..Length( B ) ] do
        if Length( B[ i ] ) = 0 then 
            B[ i ] := Zero( OriginalPathAlgebra( A ) );
        fi;
    od;
    B := Flat( B );
    elem := Flat( ExtRepOfObj( ExtRepOfObj( m ) ) );
    
    return LinearCombination( B, elem );
end 
);

#######################################################################
##
#O  ElementIn_vA_AsElementInIndecProj( <A>, <m> )
##
##  Given a finite dimension path algebra A, an element m  in  vA  for 
##  for some vertex v, this function computes the element  m  
##  correspond to in the indecomposable projective representation 
##  corresponding to the vertex  v.
##
InstallMethod ( ElementIn_vA_AsElementInIndecProj, 
    "for an element in a quotient of a path algebra",
    true,
    [ IsQuiverAlgebra, IsObject ],
    0,
    function( A, m )
  
  local leftsupport, pos, P, dimP, B, rightsupport, Q, vertices, elem, 
        r, mpart, mtemp, mtip, posmtip;
  
  if not ( m in A ) then 
    Error( "The entered element is not in the given algebra.\n" );
  fi;
  
  if IsZero( m ) then
    Error( "The entered element is zero.\n" );
  fi;
  
  leftsupport := LeftSupportOfQuiverAlgebraElement( A, m );
  if Length( leftsupport ) > 1 then
    Error( "Entered element is not left uniform.\n" );
  fi;

  pos := leftsupport[ 1 ];
  P := IndecProjectiveModules( A )[ pos ];
  dimP := DimensionVector( P );
  B := BasisOfProjectives( A )[ pos ];
  rightsupport := RightSupportOfQuiverAlgebraElement( A, m );

  Q := QuiverOfPathAlgebra( OriginalPathAlgebra( A ) );
  vertices := List( VerticesOfQuiver( Q ), v -> One( A ) * v );

  elem := Zero( P );
  for r in rightsupport do
    mpart := m * vertices[ r ];
    mtemp := mpart;
    while not IsZero( mtemp ) do
      mtip := Tip( mtemp );
      mtemp := mtemp - TipCoefficient( mtip ) * mtip;
      posmtip := Position( B[ r ], mtip );
      elem := elem + TipCoefficient( mtip ) * BasisVectors( Basis( P ) )[ Sum( dimP{ [ 1..r - 1 ] } ) + posmtip ];
    od;
  od;
  
  return elem;
end
);
