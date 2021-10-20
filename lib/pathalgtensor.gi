# GAP Implementation
# $Id: patensor.gi,v 1.5 2012/04/13 08:12:39 kristink Exp $

DeclareRepresentation(
        "IsQuiverProductDecompositionRep",
        IsComponentObjectRep and IsAttributeStoringRep,
        [ "factors", "lookup_table", "name_composer" ] );


InstallMethod( QuiverProduct,
        "for two quivers",
        ReturnTrue,
        [ IsQuiver, IsQuiver ],
        function( q1, q2 )

    local vertex_lists, vertex_names,
          arrow_lists, arrow_names, arrow_defs,
          lookup_table, make_name, make_arrow_def,
          i,
          product, decomposition;

    make_name := function( factors )
        return Concatenation( String( factors[ 1 ] ), "_", String( factors[ 2 ] ) );
    end;
    make_arrow_def := function( factors )
        return [ make_name( List( factors, SourceOfPath ) ),
                 make_name( List( factors, TargetOfPath ) ),
                 make_name( factors ) ];
    end;

    vertex_lists := Cartesian( VerticesOfQuiver( q1 ), VerticesOfQuiver( q2 ) );
    vertex_names := List( vertex_lists, make_name );

    arrow_lists := Concatenation( Cartesian( VerticesOfQuiver( q1 ), ArrowsOfQuiver( q2 ) ),
                                  Cartesian( ArrowsOfQuiver( q1 ), VerticesOfQuiver( q2 ) ) );
    arrow_names := List( arrow_lists, make_name );
    arrow_defs := List( arrow_lists, make_arrow_def );

    lookup_table := rec( );
    for i in [ 1 .. Length( vertex_lists ) ] do
        lookup_table.( vertex_names[ i ] ) := vertex_lists[ i ];
    od;
    for i in [ 1 .. Length( arrow_lists ) ] do
        lookup_table.( arrow_names[ i ] ) := arrow_lists[ i ];
    od;

    decomposition := Objectify( NewType( ListsFamily,
                                         IsQuiverProductDecomposition and IsQuiverProductDecompositionRep ),
                                rec( factors := [ q1, q2 ],
                                     lookup_table := lookup_table,
                                     name_composer := make_name ) );
    SetLength( decomposition, 2 );
    SetIsFinite( decomposition, true );

    product := Quiver( vertex_names, arrow_defs );
    SetQuiverProductDecomposition( product, decomposition );

    return product;

end );

InstallMethod( IsBound\[\], "for quiver product decomposition and positive integer",
               [ IsQuiverProductDecomposition, IsPosInt ],
function( dec, i )
  return i <= Length( dec );
end );

InstallMethod( ViewObj, "for quiver product decomposition",
               [ IsQuiverProductDecomposition ],
               NICE_FLAGS + 1,
function( dec )
  Print( "<quiver product decomposition>" );
end );

InstallMethod( PrintObj, "for quiver product decomposition",
               [ IsQuiverProductDecomposition ],
               NICE_FLAGS + 1,
function( dec )
  Print( "<quiver product decomposition>" );
end );

InstallMethod( Display, "for quiver product decomposition",
               [ IsQuiverProductDecomposition ],
function( dec )
  Print( "<quiver product decomposition>" );
end );


InstallMethod( WalkOfPathOrVertex,
        "for a path",
        [ IsPath ],
        WalkOfPath );
InstallMethod( WalkOfPathOrVertex,
        "for a vertex",
        [ IsQuiverVertex ],
        function( v ) return [ v ]; end );


InstallMethod( ReversePath,
        "for a path",
        [ IsPath ],
        function( path )
    return Product( Reversed( WalkOfPath( path ) ) );
end );


InstallMethod( IncludeInProductQuiver,
        "for list and quiver with product decomposition",
        ReturnTrue,
        [ IsDenseList, IsQuiver and HasQuiverProductDecomposition ],
        function( factors, q )

    local dec, make_name, walks, join, product_walk;

    dec := QuiverProductDecomposition( q );
    make_name := dec!.name_composer;
    
    walks := List( factors, WalkOfPathOrVertex );
    join := [ TargetOfPath( factors[ 1 ] ),
              SourceOfPath( factors[ 2 ] ) ];
    
    product_walk :=
      List( Concatenation( Cartesian( walks[ 1 ], [ join[ 2 ] ] ),
                           Cartesian( [ join[ 1 ] ], walks[ 2 ] ) ),
            function( pair ) return q.( make_name( pair ) ); end);

    return Product( product_walk );
end );


InstallMethod( ProjectFromProductQuiver,
        "for positive int and path",
        ReturnTrue,
        [ IsPosInt, IsPath ],
        function( index, path )

    local quiver, dec, project_arrow_or_vertex, walk;

    quiver := FamilyObj( path )!.quiver;
    dec := QuiverProductDecomposition( quiver );

    project_arrow_or_vertex := function( p )
        return dec!.lookup_table.( String( p ) )[ index ];
    end;

    walk := WalkOfPathOrVertex( path );
    return Product( List( walk, project_arrow_or_vertex ) );

end );


InstallMethod( \[\],
        "for quiver product decomposition",
        ReturnTrue,
        [ IsQuiverProductDecomposition, IsPosInt ],
        function( dec, index )
    return dec!.factors[ index ];
end );


InstallMethod( IncludeInPathAlgebra,
        "for path and path algebra",
        ReturnTrue,
        [ IsPath, IsQuiverAlgebra ],
        function( path, pa )
    local walk;

    walk := WalkOfPathOrVertex( path );
    return Product( List( walk, function( p ) return pa.( String( p ) ); end ) );
end );


InstallMethod( VerticesOfPathAlgebra,
        "for (quotient of) path algebra",
        [ IsQuiverAlgebra ],
        function( pa )
    return List( VerticesOfQuiver( QuiverOfPathAlgebra( pa ) ),
                 function( v ) return IncludeInPathAlgebra( v, pa ); end );
end );


InstallGlobalFunction( PathAlgebraElementTerms,
        function( elem )

    local terms, tip, tip_monom, tip_coeff, rest;

    terms := [ ];

    rest := elem;
    while true do # TODO should be while not IsZero( rest ), but that doesn't work
        tip := Tip( rest );
        tip_monom := TipMonomial( rest );
        tip_coeff := TipCoefficient( rest );
        if IsZero( tip_coeff ) then
            break;
        fi;
        rest := rest - tip;
        terms[ Length(terms) + 1 ] := rec( monom := tip_monom,
                                           coeff := tip_coeff );
    od;

    return terms;
end );


InstallMethod( SimpleTensor,
        "for list and path algebra",
        ReturnTrue,
        [ IsDenseList, IsQuiverAlgebra ],
        function( factors, pa )

    local parts, pairs, inc_paths, inc_terms;

    if IsZero( factors[ 1 ] ) or IsZero( factors[ 2 ] ) then
       return Zero( pa );
    fi;
    inc_paths :=
      function( paths )
        return IncludeInPathAlgebra( IncludeInProductQuiver( paths, QuiverOfPathAlgebra( pa ) ),
                                     pa );
    end;

    inc_terms :=
      function( terms )
        return terms[ 1 ].coeff * terms[ 2 ].coeff
               * inc_paths( [ terms[ 1 ].monom, terms[ 2 ].monom ] );
    end;

    parts := List( factors, PathAlgebraElementTerms );
    pairs := Cartesian( parts );
    return Sum( pairs, inc_terms );
end );

InstallMethod( TensorProductOfAlgebras,
        "for two path algebras",
        ReturnTrue,
        [ IsQuiverAlgebra, IsQuiverAlgebra ],
        function( pa1, pa2 )
    return TensorProductOfPathAlgebras( [ pa1, pa2 ] );
end );


#####################################################################
#
# Helene Bjørshol proved in her bachelor thesis at IMF, NTNU, spring
# 2021 that a Groebner basis of the tensor product of two quiver
# algebras, are given by the Groebner basis elements induced from 
# each factor of the tensor product in addition to all the 
# commutativity relations  
#     ( a x o(b) ) * ( t(a) x b ) - ( o(a) x b ) * ( a x t(b) )
# for all arrows a and b from the first and the second factor of the
# tensor product.  This results in a huge speed up of the function
# compared to the previous implementation.
#
InstallGlobalFunction( TensorProductOfPathAlgebras,
        function( PAs )

    local   field,  Qs,  product_q,  product_pa,  inc_q,  inc_pa,  
            make_comm_rel,  paths,  comm_rels,  get_relators,  
            get_groebner_basis,  orig_rels,  induced_rels,  
            tensor_rels,  orig_grbs,  induced_grb,  tensor_grb,  I,  
            product_pa_grb,  tensor_product;

    field := LeftActingDomain( PAs[ 1 ] );
    if field <> LeftActingDomain( PAs[ 2 ] ) then
        Error( "path algebras must be over the same field" );
    fi;

    # Construct the product quiver and its path algebra:
    Qs := List( PAs, QuiverOfPathAlgebra );
    product_q := QuiverProduct( Qs[ 1 ], Qs[ 2 ] );
    product_pa := PathAlgebra( field, product_q );

    # Inclusion functions for product quiver and path algebra:
    inc_q := function( p )
        return IncludeInProductQuiver( p, product_q );
    end;
    inc_pa := function( p )
        return IncludeInPathAlgebra( p, product_pa );
    end;

    # Function for creating the commutativity relation of a
    # given pair [a,b] where a is an arrow in the first quiver
    # and b an arrow in the second:
    make_comm_rel :=
      function( arrows )
        local paths;
        paths := [ inc_q( [ arrows[ 1 ], SourceOfPath( arrows[ 2 ] ) ] ) *
                   inc_q( [ TargetOfPath( arrows[ 1 ] ), arrows[ 2 ] ] ),
                   inc_q( [ SourceOfPath( arrows[ 1 ] ), arrows[ 2 ] ] ) *
                   inc_q( [ arrows[ 1 ], TargetOfPath( arrows[ 2 ] ) ] ) ];
        return inc_pa( paths[ 1 ] ) - inc_pa( paths[ 2 ] );
    end;

    # Create all the commutativity relations:
    comm_rels := List( Cartesian( ArrowsOfQuiver( Qs[ 1 ] ),
                                  ArrowsOfQuiver( Qs[ 2 ] ) ),
                       make_comm_rel );

    # Function for finding the relations of a quotient of a path
    # algebra (or the empty list if the algebra is not a quotient):
    get_relators :=
      function( pa )
        if IsQuotientOfPathAlgebra( pa ) then
            return RelatorsOfFpAlgebra( pa );
        else
            return [ ];
        fi;
    end;
    get_groebner_basis := function( C )
        local fam;
        
        if IsQuotientOfPathAlgebra( C ) then
            fam := ElementsFamily( FamilyObj( C ) );
            return GroebnerBasisOfIdeal( fam!.ideal );
        else
            return [ ];
        fi;
    end;
    
    # The relations of the two original path algebras:
    orig_rels := List( PAs, get_relators );
    
    # The original relations included into the product quiver:
    induced_rels := List( Concatenation( Cartesian( orig_rels[ 1 ], VerticesOfPathAlgebra( PAs[ 2 ] ) ),
                                         Cartesian( VerticesOfPathAlgebra( PAs[ 1 ] ), orig_rels[ 2 ] ) ),
                          function( factors ) return SimpleTensor( factors, product_pa ); end );

    # All the relations for the tensor product:
    tensor_rels := Concatenation( comm_rels, induced_rels );
    
    # The Groebner basis of the two original path algebras:
    orig_grbs := List( PAs, get_groebner_basis );
    
    # The original Groebner bases included into the product quiver:
    induced_grb := List( Concatenation( Cartesian( orig_grbs[ 1 ], VerticesOfPathAlgebra( PAs[ 2 ] ) ),
                                         Cartesian( VerticesOfPathAlgebra( PAs[ 1 ] ), orig_grbs[ 2 ] ) ),
                          function( factors ) return SimpleTensor( factors, product_pa ); end );
                          
    # All the Groebner basis elements for the tensor product:
    tensor_grb := Concatenation( induced_grb, comm_rels );
    
    I := Ideal( product_pa, tensor_rels );
    product_pa_grb := GroebnerBasis( I, tensor_grb );
    tensor_product := product_pa / I;
    SetTensorProductDecomposition( tensor_product, PAs );

    return tensor_product;
end );

#######################################################################
##
#O  TensorAlgebraInclusion ( < T, n > )
##
##  Returns the inclusion  A ---> A \otimes B or the inclusion
##  B ---> A \otimes B if n = 1 or n = 2 respectively. 
##  
InstallMethod( TensorAlgebraInclusion, 
    "for a IsQuiverAlgebra",
    [ IsQuiverAlgebra, IS_INT ], 0,
    function( T, n )

    local   decomp,  A,  gens,  inclusion,  images,  f;
    
    if not HasTensorProductDecomposition( T ) then
        Error( "The entered algebra is not a tensor product of two algebras.\n" );
    fi;
    decomp := TensorProductDecomposition( T );
    if n = 1 then
        A := decomp[ 1 ];
    else
        A := decomp[ 2 ];
    fi;
    gens := GeneratorsOfAlgebra( A ); 
    
    inclusion := function( x ) 
        if n = 1 then
            return SimpleTensor( [ x, One( decomp[ 2 ] ) ], T );
        else
            return SimpleTensor( [ One( decomp[ 1 ] ), x ], T );
        fi;
    end;
    images := List( gens, g -> inclusion( g ) ); 
    f := AlgebraHomomorphismByImages( A, T, gens, images ); 
    f!.generators := gens;
    f!.genimages := images;
    
    return f; 
end
  );

InstallMethod( EnvelopingAlgebra,
        "for an algebra",
        [ IsQuiverAlgebra ],
        function( pa )
    local envalg;

    envalg := TensorProductOfAlgebras( OppositeAlgebra( pa ), pa );
    SetIsEnvelopingAlgebra( envalg, true );
    if  HasIsAdmissibleQuotientOfPathAlgebra(pa) and IsAdmissibleQuotientOfPathAlgebra(pa) then
        SetFilterObj(envalg, IsAdmissibleQuotientOfPathAlgebra );
    fi; 
    return envalg;
end );


InstallMethod( IsEnvelopingAlgebra,
               "for an algebra",
               [ IsAlgebra ],
               function( alg )
    return false;
end );

InstallMethod( AlgebraAsModuleOverEnvelopingAlgebra,
        "for an enveloping algebra of a path algebra",
        [ IsQuiverAlgebra ],
        function ( A )
    local   env,  Q,  QxQ,  basis,  basis_vectors,  vertices,  
            vertex_indices,  vector_spaces,  i,  j,  dimension_vector,  
            get_coefficients,  make_map,  arrows,  
            module_specification,  a;

    if not IsFiniteDimensional(A) then
       return fail;
    fi;
    env := EnvelopingAlgebra(A); 

    Q := QuiverOfPathAlgebra( A );
    QxQ := QuiverOfPathAlgebra( env );

    basis := CanonicalBasis( A );
    basis_vectors := BasisVectors( basis );

    vertices := VerticesOfQuiver( Q );
    vertex_indices := [ 1 .. Length( vertices ) ];

    vector_spaces := NullMat( Length( vertices ), Length( vertices ) );
    for i in vertex_indices do
        for j in vertex_indices do
            vector_spaces[ i ][ j ] := PositionsNonzero( vertices[ i ] * basis_vectors * vertices[ j ] );
        od;
    od;
    dimension_vector := List( Concatenation( vector_spaces ), Length );
    
    get_coefficients := function( elem )
        return Coefficients( basis, elem );
    end;

    make_map := function( a )
        local components, source, target, source_i, target_i, source_space, target_space,
              dims, map_on_basis, map_on_basis_coeffs;

        components := [ ProjectFromProductQuiver( 1, a ),
                        ProjectFromProductQuiver( 2, a ) ];

        source := List( components, SourceOfPath );
        target := List( components, TargetOfPath );
        source_i := List( source, VertexIndex );
        target_i := List( target, VertexIndex );
        source_space := vector_spaces[ source_i[ 1 ] ][ source_i[ 2 ] ];
        target_space := vector_spaces[ target_i[ 1 ] ][ target_i[ 2 ] ];

        dims := [ Length( source_space ), Length( target_space ) ];
        if dims[ 1 ] = 0 or dims[ 2 ] = 0 then
            return fail;
        fi;

        map_on_basis := OppositePath( components[ 1 ] ) * basis_vectors{ source_space } * components[ 2 ];
        map_on_basis_coeffs := List( map_on_basis, get_coefficients );
        return map_on_basis_coeffs
               { [ 1 .. Length( map_on_basis_coeffs ) ] }
               { target_space };
    end;

    arrows := ArrowsOfQuiver( QxQ );
    module_specification := [ ];
    for a in arrows do
        if make_map( a ) <> fail then 
            Add( module_specification, [ String( a ), make_map( a ) ] );
        fi;
    od;
    
    return RightModuleOverPathAlgebra( env, dimension_vector, module_specification );

end );

#######################################################################
##
#O  DualOfAlgebraAsModuleOverEnvelopingAlgebra( <A> )
##
##  Checks if the entered algebra is finite dimensional and returns 
##  false otherwise. When the algebra is finite dimensional, then it 
##  returns the dual of the algebra as module over the enveloping algebra
##  of  A. 
## 
InstallMethod( DualOfAlgebraAsModuleOverEnvelopingAlgebra, 
    "for an algebra",
    [ IsQuiverAlgebra ], 0,
    function( A )
    
    local Aenv, M, DM, mats, op_name, de_op_name, vertices, arrows, 
          new_vertices, new_arrows, stringvertices, stringarrows, 
          vertex_positions, arrow_positions, newdimvector, newmats, 
          matrices, a, arrowentry;
    
    if not IsFiniteDimensional(A) then
        return fail;
    fi;
    #
    # By now we know that the algebra is finite dimensional.
    #
    Aenv := EnvelopingAlgebra(A);
    M    := AlgebraAsModuleOverEnvelopingAlgebra(A);
    DM   := DualOfModule(M);
    #
    #   Finding DM as a module over Aenv.
    #
    mats := MatricesOfPathAlgebraModule(DM);
    op_name := OppositeQuiverNameMap(QuiverOfPathAlgebra(A));        
    de_op_name := OppositeQuiverNameMap(OppositeQuiver(QuiverOfPathAlgebra(A)));
    vertices := VerticesOfQuiver(QuiverOfPathAlgebra(Aenv));
    arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(Aenv));
    new_vertices := List(vertices, x -> Concatenation(op_name(String(ProjectFromProductQuiver(2,x))),"_",de_op_name(String(ProjectFromProductQuiver(1,x)))));
    new_arrows := List(arrows, x -> Concatenation(op_name(String(ProjectFromProductQuiver(2,x))),"_",de_op_name(String(ProjectFromProductQuiver(1,x)))));
    stringvertices := List(vertices, x -> String(x));
    stringarrows := List(arrows, x -> String(x));
    #
    #   Finding the permutation of the vertices and the arrows.
    #
    vertex_positions := List(new_vertices, x -> Position(stringvertices, x));
    arrow_positions := List(new_arrows, x -> Position(stringarrows, x));
    #
    #   Finding the new dimension vector and the matrices of  DM  as a module over Aenv.
    #
    newdimvector := List([1..Length(vertices)], x -> DimensionVector(DM)[vertex_positions[x]]);
    newmats := List([1..Length(mats)], x -> mats[arrow_positions[x]]);
    #
    #   Creating the input for construction  DM  as a module over Aenv.
    #
    matrices := [];
    for a in arrows do
        if newdimvector[Position(vertices,SourceOfPath(a))] <> 0 and 
           newdimvector[Position(vertices,TargetOfPath(a))] <> 0 then 
            arrowentry := [[],[]];
            arrowentry[1] := String(a);
            arrowentry[2] := newmats[Position(arrows,a)];
            Add(matrices, arrowentry);
        fi;
    od;
    
    return RightModuleOverPathAlgebra(Aenv, newdimvector, matrices);
end
  );

################################################################
# Attribute TrivialExtensionOfQuiverAlgebraLevel quiver algebras
# is created so that one is avoiding naming problems when taking
# iterated trivial extensions of a quiver algebra.  This is due
# to a bug reported by Rene Marczinzik.
# If this attribute is not set, then the level is 0.  The
# trivial extension of such an algebra, will have level 1, and 
# so on.
#
InstallMethod( TrivialExtensionOfQuiverAlgebraLevel,
  "for quiveralgebras",
  true,
  [ IsQuiverAlgebra ], 0,
  function( A );

    if not HasTrivialExtensionOfQuiverAlgebraLevel( A ) then
        return 0;
    else
        return TrivialExtensionOfQuiverAlgebraLevel( A );
    fi;
end );


#######################################################################
##
#O  TrivialExtensionOfQuiverAlgebra( <A> )
##
##  Constructs the trivial extension  T(A)  of the algebra  A, that is,  
##  T(A) = A + D(A) by finding the quiver and relations of  T(A). 
## 
InstallMethod( TrivialExtensionOfQuiverAlgebra, 
    "for an algebra",
    [ IsQuiverAlgebra ], 0,
    function( A )
    
    local   arrow_add_string,  Q,  vertices,  arrows,  new_vertices,  
            new_arrows,  DA,  TopOfDAProjection,  B,  temp,  
            de_op_name,  additional_arrows,  i,  t,  Q_TE,  te_arrows,  
            string_te_arrows,  K,  KQ_TE,  relations,  new_relations,  
            r,  coeffandelem,  n,  temprel,  templist,  num_arrows,  
            V,  nontipsofradA,  nontipsofradAinKQ_TE,  W,  VWV,  
            Jpower,  AA,  BAA,  te_arrow_rep,  Qarrows,  
            Qarrows_labels,  Q_TE_to_Q,  KQ,  matrix,  setofvectors,  
            b,  occurring_arrows,  positions,  initialterm,  
            terminalterm,  Aenv,  fam,  prod,  add_arrows_labels,  
            solutions,  tempBAA,  Solutions,  TA;
    
    if not IsFiniteDimensional(A) then
        return fail;
    fi;
    arrow_add_string := String( TrivialExtensionOfQuiverAlgebraLevel( A ) + 1 );
    Q := QuiverOfPathAlgebra(A);
    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);
    #
    # The vertices in the trivial extension are the "old" ones.
    #
    new_vertices := List(vertices, v -> String(v));
    #
    # Old arrows from the original quiver.
    #
    new_arrows := List(arrows, a -> [String(SourceOfPath(a)),String(TargetOfPath(a)),String(a)]);
    #
    # Finding the "new" arrows, te_arrows.
    #
    DA := DualOfAlgebraAsModuleOverEnvelopingAlgebra(A);
    TopOfDAProjection := TopOfModuleProjection(DA);
    B := BasisVectors(Basis(Range(TopOfDAProjection)));
    temp := Flat(List(B, b -> SupportModuleElement(b))); 
    temp := List(temp, t -> CoefficientsAndMagmaElements(t![1])[1]);
    de_op_name := OppositeQuiverNameMap(OppositeQuiver(QuiverOfPathAlgebra(A))); 
    temp := List(temp, t -> [de_op_name(String(ProjectFromProductQuiver(1,t))), String(ProjectFromProductQuiver(2,t))]);
    additional_arrows := [];
    i := 0;
    for t in temp do
        i := i + 1;
        Add(additional_arrows, [t[1],t[2],Concatenation("te_a",arrow_add_string,"_",String(Position(temp,t)),"_",String(i))]);
    od;
    Append(new_arrows, additional_arrows);
    Q_TE := Quiver(new_vertices, new_arrows);
    te_arrows := ArrowsOfQuiver(Q_TE);
    string_te_arrows := List(te_arrows, a -> String(a));
    K := LeftActingDomain(A);
    KQ_TE := PathAlgebra(K,Q_TE);
    #
    # The new relations 
    #
    relations := [];
    if not IsPathAlgebra(A) then 
        relations := RelatorsOfFpAlgebra(A);
    fi;
    #
    # Transferring the relations from  A  to  T(A).
    #
    new_relations := [];
    for r in relations do
        coeffandelem := CoefficientsAndMagmaElements(r);
        n := Length(coeffandelem)/2;
        temprel := Zero(KQ_TE);
        for i in [1..n] do
            templist := List(WalkOfPath(coeffandelem[2*i-1]), w -> te_arrows[Position(string_te_arrows, String(w))]);
            temprel := temprel + coeffandelem[2*i]*One(KQ_TE)*Product(templist);  
        od;
        Add(new_relations, temprel); 
    od;
    #
    # The products te_arrows * <paths in Qarrows> * te_arrows are zero. 
    # First adding te_arrows * te_arrows.
    #
    num_arrows := NumberOfArrows(Q);    
    V := Subspace(KQ_TE, List(ArrowsOfQuiver(Q_TE){[num_arrows + 1..NumberOfArrows(Q_TE)]}, x -> One(KQ_TE)*x));
    Append(new_relations, BasisVectors(Basis(ProductSpace(V,V))));
    #
    # Now adding te_arrows * <non-trival paths in Qarrows> * te_arrows.
    #
    if IsPathAlgebra(A) then 
        nontipsofradA := BasisVectors(Basis(A)){[NumberOfVertices(Q) + 1..Dimension(A)]};
    else
        nontipsofradA := List(BasisVectors(Basis(A)){[NumberOfVertices(Q) + 1..Dimension(A)]}, x -> x![1]);
    fi;
    #
    # Transferring the nontips of rad A from  A  to  T(A).
    #
    nontipsofradAinKQ_TE := [];
    for r in nontipsofradA do
        coeffandelem := CoefficientsAndMagmaElements(r);
        n := Length(coeffandelem)/2;
        temprel := Zero(KQ_TE);
        for i in [1..n] do
            templist := List(WalkOfPath(coeffandelem[2*i-1]), w -> te_arrows[Position(string_te_arrows, String(w))]);
            temprel := temprel + coeffandelem[2*i]*One(KQ_TE)*Product(templist);  
        od;
        Add(nontipsofradAinKQ_TE, temprel); 
    od;
    W := Subspace(KQ_TE, nontipsofradAinKQ_TE);
    W := ProductSpace(W,V);
    VWV := ProductSpace(V,W);
    Append(new_relations, BasisVectors(Basis(VWV)));
    #
    # Factoring out LoewyLength(A) + 2 power of the arrow ideal in KQ_TE, and 
    # calling it  AA, this is done so that we can find the remaining relations.
    #
    Jpower := NthPowerOfArrowIdeal(KQ_TE, LoewyLength(A) + 2);
    temprel := ShallowCopy(new_relations);
    Append(temprel,Jpower);
    AA := KQ_TE/temprel;
    BAA := BasisVectors(Basis(AA));
    
    te_arrow_rep := List(B, b -> PreImagesRepresentative(TopOfDAProjection, b));
    Qarrows := ArrowsOfQuiver(Q_TE){[1..num_arrows]};
    Qarrows := Flat(List(Qarrows, a -> WalkOfPath(a))); 
    #
    # Finding the relations involving the te_arrows.
    #
    Qarrows_labels := List(arrows, a -> String(a));
    Q_TE_to_Q := function( a );
        return arrows[Position(Qarrows_labels, String(a))];
    end;
    
    KQ := OriginalPathAlgebra(A);
    matrix := []; 
    setofvectors := [];
    for b in BAA{[Length(vertices) + 1..Length(BAA)]} do
        if IsPathAlgebra(AA) then 
            temp := CoefficientsAndMagmaElements(b);
        else 
            temp := CoefficientsAndMagmaElements(b![1]);
        fi;
        occurring_arrows := WalkOfPath(temp[1]);
        positions := PositionsProperty(occurring_arrows, x -> not x in Qarrows);
        if Length(positions) = 1 then
            Add(setofvectors, Position(BAA, b)); 
            initialterm := One(KQ);
            if positions[1] > 1 then 
                temp := List(occurring_arrows{[1..positions[1] - 1]}, a -> Q_TE_to_Q(a));
                initialterm := initialterm*Product(temp);
            fi;
            terminalterm := One(KQ);
            if positions[1] < Length(occurring_arrows) then
                temp := List(occurring_arrows{[positions[1] + 1..Length(occurring_arrows)]}, a -> Q_TE_to_Q(a));                
                terminalterm := terminalterm*Product(temp);
            fi;
            Add(matrix, [initialterm, occurring_arrows[positions[1]], terminalterm]);
        fi;
    od;
    setofvectors := Set(setofvectors);
    Aenv := RightActingAlgebra(DA);
    fam := ElementsFamily(FamilyObj(A)); 
    if IsPathAlgebra(A) then
        prod := List(matrix, m -> [OppositePathAlgebraElement(m[1]), m[2], m[3]]); 
    else
        prod := List(matrix, m -> [OppositePathAlgebraElement(ElementOfQuotientOfPathAlgebra(fam, m[1]*One(KQ), true)), m[2], ElementOfQuotientOfPathAlgebra(fam, m[3]*One(KQ), true)]); 
    fi;
    add_arrows_labels := List(additional_arrows, a -> a[3]);
    prod := List(prod, m -> te_arrow_rep[Position(add_arrows_labels,String(m[2]))]^SimpleTensor([m[1],m[3]], Aenv));
    matrix := List(prod, x -> Flat(x![1]![1]));
    solutions := NullspaceMat(matrix);
    tempBAA := BAA{setofvectors};
    tempBAA := List(tempBAA, t -> t![1]);
    #
    # Getting the relations involving the arrows  te_arrows
    #
    Solutions := List(solutions, s -> LinearCombination(tempBAA, s));
    
    Append(new_relations, Solutions); 
    TA := KQ_TE/new_relations;
    SetTrivialExtensionOfQuiverAlgebraLevel( TA, TrivialExtensionOfQuiverAlgebraLevel( A ) + 1 );
    
    return TA;
end
  );

#######################################################################
##
#A  EnvelopingAlgebraHomomorphism( <f> )
##
##  Given an algebra homomorphism, this constructs corresponding algebra
##  homomorphism between the enveloping algebras.
## 
InstallMethod( EnvelopingAlgebraHomomorphism, 
    "for an algebra",
    [ IsAlgebraHomomorphism ], 0,
    function( f )
    
    local   A,  Aop,  Aenv,  B,  Benv,  fop,  gensAenv,  
            coeffgensAenv,  firstcoordinates,  secondcoordinates,  
            coeffisients,  n,  images,  g;
    
    A := Source( f );
    Aop := OppositeAlgebra( A );
    Aenv := EnvelopingAlgebra( A );
    B := Range( f );
    Benv := EnvelopingAlgebra( B );
    fop := OppositeAlgebraHomomorphism( f );
    gensAenv := GeneratorsOfAlgebraWithOne( Aenv );
    coeffgensAenv := List( gensAenv, g -> CoefficientsAndMagmaElements( g![ 1 ] ) );
    
    firstcoordinates := List( coeffgensAenv, m -> ProjectFromProductQuiver( 1, m[ 1 ] ) );
    secondcoordinates := List( coeffgensAenv, m -> ProjectFromProductQuiver( 2, m[ 1 ] ) );
    coeffisients := List( coeffgensAenv, m -> m[ 2 ] );
    firstcoordinates := List( firstcoordinates, m -> ImageElm( fop, One( Aop ) * m ) );
    secondcoordinates := List( secondcoordinates, m -> ImageElm( f, One( A ) * m ) );
    n := Length( coeffgensAenv );    
    images := List( [ 1..n ], i -> coeffisients[ i ] * SimpleTensor( [ firstcoordinates[ i ], secondcoordinates[ i ] ] , Benv ) );
    
    g := AlgebraHomomorphismByImages( Aenv, Benv, gensAenv, images );
    g!.generators := gensAenv;
    g!.genimages := images;
    
    return g;
end
  );

#######################################################################
##
#O  TrivialExtensionOfQuiverAlgebraProjection( <A> )
##
##  Constructs the natural algebra projection form trivial extension  
##  T(A)  to  A.
## 
InstallMethod( TrivialExtensionOfQuiverAlgebraProjection, 
    "for an algebra",
    [ IsQuiverAlgebra ], 0,
    function( A )
    
    local   TA,  gensTA,  gensA,  newarrows,  i,  f;
    
    if not IsFiniteDimensional( A ) then
        return fail;
    fi;
    TA := TrivialExtensionOfQuiverAlgebra( A );
    gensTA := GeneratorsOfAlgebraWithOne( TA );
    gensA := ShallowCopy( GeneratorsOfAlgebraWithOne( A ) );
    newarrows := Length( gensTA ) - Length( gensA );
    for i in [ 1..newarrows ] do
        Add( gensA, Zero( A ) );
    od;
    f := AlgebraHomomorphismByImages( TA, A, gensTA, gensA );
    f!.generators := gensTA;
    f!.genimages := gensA;
    
    return f;
end
  );