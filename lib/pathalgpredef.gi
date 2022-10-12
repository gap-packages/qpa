# GAP Implementation
# This file was generated from 
# $Id: predefalgs.gi,v 1.6 2012/08/01 16:01:10 sunnyquiver Exp $
InstallImmediateMethod( IsFiniteTypeAlgebra, IsNakayamaAlgebra, 0, A -> true);
        
#######################################################################
##
#O  NakayamaAlgebra( K, <admiss_seq> )
##
##  Given a field  <K>  and an admissible sequence  <admiss_seq>  this 
##  function constructs a Nakayama algebra with this data. 
##
InstallMethod ( NakayamaAlgebra,
    "for a field and an admissible sequence",
    [IsField, IsList], 0,
    function( K, admiss_seq )
    
    local num_vert, # number of vertices in the quiver we are creating.
          test,     # the result of the test of being an admissible seq.
          i, j,     # integer variables for for-loops.
          q, r,     # for writing a number n: n = q*num_vert + r.
          Q, KQ,    # the quiver and the path algebra for the 
          # Nakayama algebra we are creating.
          rels,     # relations of the Nakayama algebra we are creating. 
          new_rel,  # partial relation. 
          arrows,   # arrows of the quiver we are creating.
          cycle,    # if we are in the A-tilde case and the relations.
                    # are long, this is used to store the cycles involved
                    # in the relations. 
          I,        # the ideal of relations.
          gb, gbb,  # Groebner basis for the ideal I.
          A;
#
# Testing if we have an admissible sequence
#
    num_vert := Length(admiss_seq);
    test := true;
    i := 1;
    while ( test and ( i < num_vert ) ) do
        test := ( ( admiss_seq[i+1] >= admiss_seq[i] - 1 ) and
                  ( admiss_seq[i] - 1 >= 1 ) );  
        i := i + 1;
    od;
    if test then 
        test := ( admiss_seq[1] >= admiss_seq[num_vert] - 1 );
    fi;
#
# If an admissible sequence, the case admiss_seq[num_vert] = 1, is the
# A_n case, and the other case is A_n-tilde case. 
#
    if test then 
        if ( admiss_seq[num_vert] = 1 ) then
#
# This is the A_n case.
#
            arrows := [];
            for i in [1..num_vert - 1] do
                Add(arrows,[i,i+1]);
            od;
            Q := Quiver(num_vert,arrows);
            KQ := PathAlgebra(K,Q); 
            arrows := GeneratorsOfAlgebra(KQ){[num_vert+1..2*num_vert-1]}; 
            rels := [];
            for i in [1..num_vert - 1] do
                if (( admiss_seq[i] < num_vert - i + 1 ) and 
                    ( admiss_seq[i] - 1 <> admiss_seq[i+1] )) then
                    new_rel := One(KQ);
                    for j in [i..i+admiss_seq[i] - 1] do
                        new_rel := new_rel*arrows[j];
                    od;
                    Add(rels,new_rel); 
                fi;
            od;
#
# end of A_n case
#
        else
#
# This is the A_n-tilde case.
            arrows := [];
            for i in [1..num_vert - 1] do
                Add(arrows,[i,i+1]);
            od;
            Add(arrows,[num_vert,1]); 
            Q := Quiver(num_vert,arrows);
            KQ := PathAlgebra(K,Q); 
            arrows := GeneratorsOfAlgebra(KQ){[num_vert+1..2*num_vert]}; 
            rels := [];
#
# Relations starting in the vertices 1..num_vert - 1. 
#
            for i in [1..num_vert - 1] do
                if ( admiss_seq[i] - 1 <> admiss_seq[i+1] ) then
                    new_rel := One(KQ);
                    r := admiss_seq[i] mod num_vert;
                    q := (admiss_seq[i] - r)/num_vert;
                    cycle := One(KQ);
                    for j in [i..num_vert] do
                        cycle := cycle*arrows[j];
                    od;
                    for j in [1..i-1] do
                        cycle := cycle*arrows[j];
                    od;
                    for j in [1..q] do 
                        new_rel := new_rel*cycle;
                    od;
                    if ( i + r - 1 <= num_vert ) then
                        for j in [i..i+r - 1] do
                            new_rel := new_rel*arrows[j];
                        od;
                    else
                        for j in [i..num_vert] do
                            new_rel := new_rel*arrows[j];
                        od;
                        for j in [1..r - (num_vert - i) - 1] do
                            new_rel := new_rel*arrows[j];
                        od;
                    fi;
                    Add(rels,new_rel); 
                fi;
            od;
#
# Relation starting in the vertex num_vert. 
#
            if ( admiss_seq[num_vert] - 1 <> admiss_seq[1] ) then
                new_rel := One(KQ);
                r := admiss_seq[num_vert] mod num_vert;
                q := (admiss_seq[num_vert] - r)/num_vert;
                cycle := One(KQ)*arrows[num_vert];
                for j in [1..num_vert-1] do
                    cycle := cycle*arrows[j];
                od;
                for j in [1..q] do 
                    new_rel := new_rel*cycle;
                od;
                if ( r <> 0 ) then
                    new_rel := new_rel*arrows[num_vert];
                    for j in [1..r - 1] do
                        new_rel := new_rel*arrows[j];
                    od;
                fi;
                Add(rels,new_rel);
            fi;
#
# End of A_n-tilde case.
#
        fi;
        if Length(rels) = 0 then 
            SetFilterObj(KQ, IsNakayamaAlgebra and IsPathAlgebra );
            
            return KQ;
        else
            I := Ideal(KQ,rels);
            gb := GroebnerBasisFunction(KQ)(rels,KQ);
            gbb := GroebnerBasis(I,gb);
            A := KQ/I;
            SetFilterObj(A, IsNakayamaAlgebra and IsAdmissibleQuotientOfPathAlgebra );               
            return A;
        fi; 
    else
#
# A valid admissible sequence was not entered. 
# 
        Print("The admissible sequence is not a valid one.\n");
        return admiss_seq;
    fi; 
    
end
);

InstallOtherMethod ( NakayamaAlgebra,
    "for an admissible sequence and a field",
    [IsList, IsField], 0,
        function( admiss_seq , K)
    
    return NakayamaAlgebra(K,admiss_seq);
end
);

#######################################################################
##
#O  CanonicalAlgebra( <K>, <weights>, <relcoeff> )
##
##  Given a field  K, a sequence of  weights  and a sequence of 
##  coefficients for the relations, this function constructs the 
##  canonical algebra with this data. If the number of weights is
##  two, then the number of coefficients for the relations must be
##  zero, that is, represented by an empty list, [].
##
InstallMethod ( CanonicalAlgebra,
    "for a field, list of lengths of arms, and coefficients for the relations", 
    [IsField, IsList, IsList], 0,
    function( K, weights, relcoeff)
    
    local n, vertices, i, j, arrows, Q, KQ, num_vert, generators, arms,
          start, temp, relations, I, gb, gbb, A; 

        #
        # Checking the input.
        #
    if Length(weights) <> Length(relcoeff) + 2 then
        Error("number of arms in the quiver don't match the number of coefficients for the relations in the quiver,");
    fi;
    if not ForAll( weights, x -> x >= 2) then
        Error("not all weights are greater or equal to 2,");
    fi;
    if not ForAll( relcoeff, x -> x in K ) then
        Error("not all coefficients for the relations are in the entered field,");
    fi;
    if not ForAll( relcoeff, x -> x <> Zero(K) ) then
        Error("not all coefficients for the relations are non-zero,");
    fi;
        #
        # When input is good, start constructing the algebra. 
        # First construct the vertices.
        #
    n := Length(weights);
    vertices := [];
    Add(vertices,"v");
    for i in [1..n] do
        for j in [2..weights[i]] do
            Add(vertices,Concatenation("v",String(i),String(j)));
        od;
    od;
    Add(vertices,"w");
        #
        # Next construct the arrows.
        #
    arrows := [];
    for i in [1..n] do
        Add(arrows,["v",Concatenation("v",String(i),String(2)),Concatenation("a",String(i))]);
        for j in [2..weights[i] - 1] do
            Add(arrows,[Concatenation("v",String(i),String(j)), 
                    Concatenation("v",String(i),String(j+1)),
                    Concatenation("a",String(i),String(j))]);
        od;
        Add(arrows,[Concatenation("v",String(i),String(weights[i])),"w",Concatenation("b",String(i))]);
    od;
        # 
        # Constructing the quiver and path algebra.
        #
    Q := Quiver(vertices,arrows);
    KQ := PathAlgebra(K,Q);
        #
        # If n = 2, it is a path algebra, otherwise construct relations 
        # and quotient of the path algbra just constructed.
        #
    if n = 2 then 
        SetFilterObj(KQ, IsCanonicalAlgebra and IsPathAlgebra );

        return KQ;
    else
        # 
        # First finding all the arrows as elements in the path algebra.
        #
        num_vert := NumberOfVertices(Q);
        generators := GeneratorsOfAlgebra(KQ){[num_vert + 1..num_vert + NumberOfArrows(Q)]};
        #
        # Computing the products of all arrows in an arm of the quiver as 
        # elements in the path algebra. 
        # 
        arms := [];
        start := 1;
        for i in [1..n] do
            temp := One(KQ);
            for j in [start..start + weights[i] - 1] do
                temp := temp*generators[j];
            od;
            start := start + weights[i];
            Add(arms,temp);
        od;
        relations := [];
        # 
        # Constructing the relations and the quotient.
        #
        for i in [3..n] do
            Add(relations,arms[i] - arms[2] + relcoeff[i - 2]*arms[1]);
        od;
        I := Ideal(KQ,relations);
        gb := GroebnerBasisFunction(KQ)(relations,KQ);
        gbb := GroebnerBasis(I,gb);
        A := KQ/I;
        SetFilterObj(A, IsCanonicalAlgebra and IsAdmissibleQuotientOfPathAlgebra );
        return A;
    fi;
end
);

#######################################################################
##
#O  CanonicalAlgebra( <K>, <weights>)
##
##  CanonicalAlgebra with only two arguments, the field and the 
##  two weights.
##
InstallOtherMethod ( CanonicalAlgebra,
    "for a field and a list of lengths of two arms", 
    [IsField, IsList ], 0,
    function( K, weights);

    if Length(weights) <> 2 then 
        Error("the list of weights is different from two, need to enter coefficients for the relations then,");
    fi;
    
    return CanonicalAlgebra(K,weights,[]);
end
);

#######################################################################
##
#O  KroneckerAlgebra( <K>, <n> )
##
##  Given a field  K and a positive integer  n, this function constructs  
##  the Kronecker algebra with  n  arrows. 
##
InstallMethod ( KroneckerAlgebra,
    "for a field and a positive integer", 
    [IsField, IS_INT], 0,
    function( K, n)
    
    local arrows, i, Q, KQ;

    if n < 2 then
        Error("the number of arrows in the Kronecker quiver must be greater or equal to 2,");
    fi;
    arrows := [];
    for i in [1..n] do
        Add(arrows,[1,2,Concatenation("a",String(i))]);
    od;
    Q := Quiver(2,arrows);
    KQ := PathAlgebra(K,Q);
    SetFilterObj(KQ,IsKroneckerAlgebra and IsPathAlgebra);
    return KQ;
end
);


#################################################################
##
#P  IsSpecialBiserialQuiver( <Q> ) 
## 
##  It tests if every vertex in quiver <Q> is a source (resp. target) 
##  of at most 2 arrows.
##
##  NOTE: e.g. path algebra of one loop IS NOT special biserial, but
##        one loop IS special biserial quiver.
##
InstallMethod( IsSpecialBiserialQuiver,
    "for quivers",
    true,
    [ IsQuiver ], 0,
    function ( Q )
    
    local vertex_list, v;
    
    vertex_list := VerticesOfQuiver(Q);
    for v in vertex_list do
        if OutDegreeOfVertex(v) > 2 then
            return false;
        fi;
        if InDegreeOfVertex(v) > 2 then
            return false;
        fi;
    od;
    
    return true;
end
); # IsSpecialBiserialQuiver


#################################################################
##
#P  IsSpecialBiserialAlgebra( <A> ) 
##  <A> = a path algebra
##
##  It tests if a path algebra <A> is an algebra of  a quiver Q which is
##  (IsSpecialBiserialQuiver and IsAcyclicQuiver) and the ideal 0 satisfies
##  the "special biserial" conditions (cf. comment below).
##  
##  NOTE: e.g. path algebra of one loop IS NOT special biserial, but
##        one loop IS special biserial quiver.
##		
InstallMethod( IsSpecialBiserialAlgebra,
    "for path algebras",
    true,
    [ IsPathAlgebra ], 0,
    function ( A )
    
    local Q, alpha; 
    
    Q := QuiverOfPathAlgebra(A);
    if not IsSpecialBiserialQuiver(Q) then
        return false;
    fi;
    
    for alpha in ArrowsOfQuiver(Q) do
        if OutDegreeOfVertex(TargetOfPath(alpha)) > 1 then
            return false;
        fi;          
        if InDegreeOfVertex(SourceOfPath(alpha)) > 1 then  
            return false;
        fi;
    od;
    
    return IsAcyclicQuiver(Q);  # <=> 0 is an admissible ideal
end
); # IsSpecialBiserialAlgebra for IsPathAlgebra


########################################################################
##
#P  IsSpecialBiserialAlgebra( <A> ) 
##  <A> = a quotient of a path algebra
##
##  It tests if an original path algebra is an algebra of a quiver Q which is
##  IsSpecialBiserialQuiver, I is an admissible ideal and I satisfies 
##  the "special biserial" conditions, i.e.:
##  for any arrow a there exists at most one arrow b such that ab\notin I
##  and there exists at most one arrow c such that ca\notin I.
##
  
InstallMethod( IsSpecialBiserialAlgebra,
    "for quotients of path algebras",
    true,
    [ IsQuotientOfPathAlgebra ], 0,
    function ( A )
    
    local Q, PA, I, v, alpha, beta, not_in_ideal;
    
    PA := OriginalPathAlgebra(A);
    Q := QuiverOfPathAlgebra(PA);
    if not IsSpecialBiserialQuiver(Q) then
        return false;
    fi;
    
    I := ElementsFamily(FamilyObj(A))!.ideal;
    
    if not IsAdmissibleIdeal(I) then
        return false;
    fi;
    
    for alpha in ArrowsOfQuiver(Q) do
        not_in_ideal := 0;
        for beta in OutgoingArrowsOfVertex(TargetOfPath(alpha)) do
            if not ElementOfPathAlgebra(PA, alpha*beta) in I then
                not_in_ideal := not_in_ideal + 1;
            fi;
        od;
        if not_in_ideal > 1 then
            return false;
        fi;          
        not_in_ideal := 0;
        for beta in IncomingArrowsOfVertex(SourceOfPath(alpha)) do
            if not ElementOfPathAlgebra(PA, beta*alpha) in I then
                not_in_ideal := not_in_ideal + 1;
            fi;
        od;
        if not_in_ideal > 1 then
            return false;
        fi;
    od;
    
    return true;
end
); # IsSpecialBiserialAlgebra for IsQuotientOfPathAlgebra


#################################################################
##
#P  IsStringAlgebra( <A> )
##  <A> = a path algebra
##
##  Note: A path algebra is a string algebra <=> it is a special biserial algebra
##
InstallMethod( IsStringAlgebra,
    "for quotients of path algebras",
    true,
    [ IsPathAlgebra ], 0,
    function ( A )
    return IsSpecialBiserialAlgebra(A);
end
); # IsStringAlgebra for IsPathAlgebra                      


#################################################################
##
#P  IsStringAlgebra( <A> )
##  <A> = a quotient of a path algebra
##
##  Note: kQ/I is a string algebra <=> kQ/I is a special biserial algebra
##                                     and I is a monomial ideal.                                         
InstallMethod( IsStringAlgebra,
    "for quotients of path algebras",
    true,
    [ IsQuotientOfPathAlgebra ], 0,
    function ( A )
    
    local I;
    
    if not IsSpecialBiserialAlgebra(A) then
        return false;
    fi;
    
    I := ElementsFamily(FamilyObj(A))!.ideal;
    return IsMonomialIdeal(I);
end
  ); # IsStringAlgebra for IsQuotientOfPathAlgebra 

#################################################################
##
#O  PosetAlgebra( <K>, <P> )
##  
##  Takes as arguments a field  <K>  and a poset  <P> and returns
##  the associated poset algebra  <KP>  as a quiver algebra.
##

InstallMethod( PosetAlgebra,
    "for sets",
    [ IsField, IsPoset ],
        
    function( K, P)
    local   minimalrelations,  n,  vertices,  arrows,  i,  j,  Q,  KQ,  
            arrowsofQ,  I,  radKQsquare,  BradKQsquare,  primidemp,  
            V,  lessthan,  temp,  generators,  A,  relations,  r,  
            iVj,  B,  k;
    
    minimalrelations := P!.minimal_relations;
    n := Size(P);     
    # 
    # Setting the names for the vertices of the quiver to have the 
    # same names as the names of the elements of the poset P.
    #
    vertices := List(P!.set, p -> String(p));
    #
    # For each minimal relation  x < y, we define an arrow of the quiver 
    # from the vertex  x  to the vertex  y  with name "x_y".
    #
    arrows := [];
    for i in [ 1..n ] do
        for j in [ 1..n ] do 
            if DirectProductElement( [ P!.set[ i ], P!.set[ j ] ] ) in minimalrelations then 
                Add( arrows, [ String( P!.set[ i ] ), String( P!.set[ j ] ), Concatenation( String( P!.set[ i ] ),"_",String( P!.set[ j ] ) ) ] );
            fi;
        od;
    od;
    # 
    # Creating the quiver and the path algebra thereof.
    #
    Q := Quiver( vertices, arrows );
    KQ := PathAlgebra( K, Q );
    #
    # Finally we find the relations by computing a basis {b_1,..., b_vw} for 
    # vJ^2w for every pair of vertices v and w in Q when  v < w  in P, and if 
    # it has dimension greater or equal to 2, then we add the relations
    # b_1 - b_j for all j >= 2.
    # 
    arrowsofQ := List( ArrowsOfQuiver(Q), a -> One( KQ ) * a ); 
    I := Ideal( KQ, arrowsofQ );
    radKQsquare := ProductSpace( I, I ); 
    BradKQsquare := BasisVectors( Basis( radKQsquare ) ); 
    primidemp := List( VerticesOfQuiver( Q ), v -> v*One( KQ ) );
    V := []; 
    lessthan := PartialOrderOfPoset( P );
    for i in [ 1..n ] do
        for j in [1..n] do
            if ( i <> j ) and ( lessthan( P!.set[ i ], P!.set[ j ] ) ) then
                generators := Filtered( primidemp[ i ] * BradKQsquare * primidemp[ j ], x -> x <> Zero( x ) );
                if Length( generators )  > 1 then 
                    Add( V, [ i, j, Subspace( KQ, generators ) ] );
                fi;
            fi;
        od;
    od;
    if Length( V ) = 0 then
        A := KQ;
    else
        relations := [];
        for r in [ 1..Length( V ) ] do
            i := V[ r ][ 1 ];
            j := V[ r ][ 2 ];
            iVj := V[ r ][ 3 ];
            for j in [1..n] do
                if i <> j and Dimension( iVj ) > 1 then
                    B := BasisVectors( Basis( iVj ) );
                    for k in [ 2..Dimension( iVj ) ] do
                        Add( relations, B[ 1 ] - B[ k ] );
                    od;
                fi;
            od;
        od;
        A := KQ/relations;
    fi;
    SetPosetOfPosetAlgebra( A, P ); 
    SetFilterObj( A, IsPosetAlgebra ); 
    SetFilterObj( A, IsFiniteGlobalDimensionAlgebra );

    return A; 
end
);



#######################################################################
##
#O  BrauerConfigurationAlgebra( K, <brauer_configuration> )
##
##  Given a field  <K>  and a brauer configuration
##  <brauer_configuration> in the form
##  [[vertices], [polygons], [orientations]] this function constructs
##  the Brauer configuration algebra associated to the Brauer
##  configuration. 
##
InstallMethod ( BrauerConfigurationAlgebra,
    "for a field and a valid list",
    [IsField, IsList], 0,
    function( K, brauer_configuration )

    local vertices, polygons, orientations,
        vertex_names,  polygon_names, special_cycles, polygons_containing, num_vertices,
        arrows, quiver_vertices, quiver, path_algebra,
        type1relations, type2relations, type3relations,
        path, vertex1_index, vertex2_index, sc, o, p, i, j, k, a, b, v, A,
        special_cycle;
    #
    #Checking that the input is valid
    #
    if IsField(K) = false then
        Print("The second argument must be a field.\n");
        return fail;
    fi;

    if IsList(brauer_configuration) = false or Length(brauer_configuration) <> 3 or IsList(brauer_configuration[1]) = false or IsList(brauer_configuration[2]) = false or IsList(brauer_configuration[3]) = false then
        Print("Brauer Algebra input must contain 3 arrays: [vertices], [edges or polygons], [orientations].\n");
        return fail;
    fi;

    vertices := brauer_configuration[1];
    vertex_names := [];
    num_vertices := Length(vertices);

    if (num_vertices = 0) then
        Print("There must be at least 1 vertex\n");
        return fail;
    fi;

    for v in vertices do
        if IsList(v) = false or Length(v) <> 2 or (IsString(v[1]) and IsInt(v[2])) = false then
            Print("Vertices should be of the form: [name, multiplicity].\n");
            return fail;
        elif v[2] < 1 then
            Print("Multiplicities must be positive integers.\n");
            return fail;
        else
            Add(vertex_names, v[1]);
        fi;
    od;

    if (IsDuplicateFree(vertex_names)) = false then
        Print("Each vertex must have a unique name.\n");
        return fail;
    fi;

    polygons := brauer_configuration[2];
    polygon_names := [];

    for p in polygons do
        if IsList(p) = false or Length(p) < 2 then
            Print("Edges or Polygons should be lists of length 2 or greater.\n");
            return fail;
        else 
            if (p[1] in vertex_names) then
                Print("An edge or polygon cannot have the same name as a vertex.\n");
                return fail;
            fi;
            Add(polygon_names, p[1]);
            for i in [2.. Length(p)] do
                if (p[i] in vertex_names) = false then
                    Print("An edge or polygon cannot contain a vertex which was not listed.\n");
                    return fail;
                fi;
            od;
        fi;
    od;

    if (IsDuplicateFree(polygon_names)) = false then
        Print("Edges or polygons must have unique names.\n");
        return fail;
    fi;

    orientations := brauer_configuration[3];
    if (Length(orientations) <> num_vertices) then
        Print("There must be an orientation corresponding to each vertex.\n");
        return fail;
    fi;
    for o in orientations do
        if IsList(o) = false  or Length(o) < 1 then
            Print("Orientations must be lists of length 1 or greater.\n");
            return fail;
        fi;
        for p in o do
            if (p in polygon_names) = false then
                Print("Orientations may only contain edges or polygons which were listed.\n");
                return fail;
            fi;
        od;
    od;

    polygons_containing := [];
    for i in [1.. num_vertices] do
        Add(polygons_containing, []);
        for p in polygons do
            if vertex_names[i] in p then
                Add(polygons_containing[i], p[1]);
            fi;
        od;
    od;

    for i in [1.. num_vertices] do
        Sort(orientations[i]);
        Sort(polygons_containing[i]);
        if orientations[i] <> polygons_containing[i] then
            Print("Orientations must contain exactly the edges or polygons which contain the corrosponding vertex.\n");
            return fail;
        fi;
    od;
    #
    #Generating the arrows for the quiver.
    #
    arrows := [];
    special_cycles := [];
    for i in [1.. num_vertices] do
        #
        #We only want to create a self loop if the corresponding vertex has multiplicity
        #
        if Length(orientations[i]) > 1 or vertices[i][2] > 1 then 
            sc := [];
            for j in [1.. Length(orientations[i]) - 1] do
                Add(arrows, [orientations[i][j], orientations[i][j + 1],
                    Concatenation(vertices[i][1], "_", String(j), "_from_",
                    orientations[i][j], "_to_", orientations[i][j + 1])]);
                Add(sc, arrows[Length(arrows)][3]);
            od;
            Add(arrows, [orientations[i][Length(orientations[i])], orientations[i][1],
                Concatenation(vertices[i][1], "_", String(Length(orientations[i])), "_from_",
                orientations[i][Length(orientations[i])], "_to_", orientations[i][1])]);
            Add(sc, arrows[Length(arrows)][3]);
            Add(special_cycles, sc);
        else
            Add(special_cycles, []);
        fi;
    od;

    #
    #Retrieving names of the polygons to be used as vertices in the quiver.
    #
    quiver_vertices := [];
    for p in polygons do
        Add(quiver_vertices, p[1]);
    od;

    quiver := Quiver(quiver_vertices, arrows);
    path_algebra := PathAlgebra(K, quiver);

    #
    #Calculating Type 1 Relations
    #
    type1relations := [];
    for i in [1.. num_vertices - 1] do
        for j in [i + 1.. num_vertices] do
            for a in special_cycles[i] do
                for b in special_cycles[j] do
                    Add(type1relations, path_algebra.(a) * path_algebra.(b));
                    Add(type1relations, path_algebra.(b) * path_algebra.(a));
                od;
            od;
        od;
    od;

    #
    #Helper function, returns the special cycle corresponding to a vertex, beginning at a certain polygon
    #
    special_cycle := function(vertex_index, polygon_index)
        local i, j, path;

        path := path_algebra.(special_cycles[vertex_index][polygon_index]);
        for i in [polygon_index + 1.. Length(special_cycles[vertex_index])] do
            path := path * path_algebra.(special_cycles[vertex_index][i]);
        od;
        for i in [1.. polygon_index - 1] do
            path := path * path_algebra.(special_cycles[vertex_index][i]);
        od;
        return path;
    end;

    #
    #Calculating Type 2 Relations
    #
    type2relations := [];
    for i in [1.. Length(polygons)] do
        if InDegreeOfVertex(VerticesOfQuiver(quiver)[i]) > 1 then
            for j in [2.. Length(polygons[i]) - 1] do
                vertex1_index := Position(vertex_names, polygons[i][j]);
                if IsEmpty(special_cycles[vertex1_index]) = false then
                    for k in [j + 1 .. Length(polygons[i])] do
                        vertex2_index := Position(vertex_names, polygons[i][k]);
                        if IsEmpty(special_cycles[vertex2_index]) = false then
                            Add(type2relations,
                                special_cycle(vertex1_index, Position(orientations[vertex1_index], polygons[i][1]))^vertices[vertex1_index][2]
                                - special_cycle(vertex2_index, Position(orientations[vertex2_index], polygons[i][1]))^vertices[vertex2_index][2]);
                        fi;
                    od;
                fi;
            od;
        fi;
    od;

    #
    #Calculating Type 3 Relations
    #
    type3relations := [];
    for i in [1.. num_vertices] do
        if IsEmpty(special_cycles[i]) = false then
            for j in [1.. Length(special_cycles[i])] do
                Add(type3relations, special_cycle(i, j)^vertices[i][2] * path_algebra.(special_cycles[i][j]));
            od;
        fi;
    od;

    A := path_algebra/Concatenation(type1relations, type2relations, type3relations);
    SetIsSymmetricAlgebra( A, true );

    return A;
end
);

InstallMethod( PreprojectiveAlgebra,
    "for a path algebra module and an integer",
    [ IsPathAlgebraMatModule, IsInt ], 0,
    function( M, n )
    
  local R, P, dimlistP, revdimlistP, top, K, degrees, origin, i, 
        target, V, dimV, products, zero, i_first, j, f, r, first, s, 
        g, gpowers, u, h, temp, i_times_j, ir_first, v, A;

  R := RightActingAlgebra( M );
  if not IsHereditaryAlgebra( R ) then 
    TryNextMethod( );
    # TODO: Implement method for not hereditary algebras.
  fi;
  P := List( [ 0..n ], i -> HomOverAlgebraWithBasisFunction( M, TrD( M, i ) ) );
  if Dimension( TrD( M, n+1 ) ) > 0 then
    Error( "Must increase the second argument.\n" );
  fi;
  dimlistP := List( P, p -> Length( p[ 1 ] ) );
  revdimlistP := Reversed( dimlistP );
  top := n + 2 - PositionNonZero( revdimlistP ); 
  P := P{ [ 1..top ] };
  dimlistP := dimlistP{ [ 1..top ] };
  K := LeftActingDomain( M );
  degrees := [];
  origin := 1;
  for i in [ 1..top ] do
    target := origin + dimlistP[ i ] - 1;
    Add( degrees, [ origin..target ] );
    origin := target + 1;
  od;
  V := FullRowSpace( K, target );
  dimV := Dimension( V );
  products := EmptySCTable( dimV, Zero( K ) );
  zero := List( [ 1..dimV ], u -> [ Zero( K ), u ] );
  zero := Flat( zero );

  for i in [ 1..top ] do 
    if i = 1 then 
      i_first := 0;
    else
      i_first := Maximum( degrees[ i - 1 ] );
    fi;
    for j in degrees[ i ] do
      f := P[ i ][ 1 ][ j - i_first ];
      for r in [ 1..top ] do
        if r = 1 then 
          first := 0;
        else
          first := Maximum( degrees[ r - 1 ] );
        fi;
        for s in degrees[ r ] do
          if i + r - 1 > top then 
            SetEntrySCTable( products, j, s, zero );
          else            
            g := P[ r ][ 1 ][ s - first ];
            gpowers := g;
            for u in [ 1..i - 1 ] do
              gpowers := DualOfModuleHomomorphism( gpowers );
              gpowers := TransposeOfModuleHomomorphism( gpowers );
            od;
            h := f * gpowers;
            temp := P[ i + r - 1 ][ 2 ]( h );
            i_times_j := [ ];
            if i + r = 2 then 
              ir_first := 0;
            else
              ir_first := Maximum( degrees[ i + r - 2 ] );
            fi;
            for v in [ 1..dimV ] do
              if v in degrees[ i + r - 1 ] then
                Append( i_times_j, [ temp[ v - ir_first ], v ] );
              else
                Append( i_times_j, [ Zero( K ), v ] );
              fi;
            od;
            SetEntrySCTable( products, j, s, i_times_j );            
          fi;
        od;
      od;
    od;
  od;
  A := AlgebraByStructureConstants( K, products );
  
  return AlgebraAsQuiverAlgebra( A );
end
);

InstallOtherMethod( PreprojectiveAlgebra,
   "for a path algebra",
   [ IsPathAlgebra ], 0,
   function( A ) 

  local K, Q, Qdouble, preA, arrows, relation, num_arrows, i, 
        vertices, relations;

   K := LeftActingDomain( A );
   Q := QuiverOfPathAlgebra( A );
   if not IsAcyclicQuiver( Q ) then
     Error( "The entered quiver has an oriented cycle.\n" );
   fi;
   Qdouble := DoubleQuiver( Q );
   preA := PathAlgebra( K, Qdouble );
   arrows := ArrowsOfQuiver( Qdouble );
   relation := Zero( preA );
   num_arrows := Length( arrows ) / 2;
   
   for i in [ 1..num_arrows ] do
     relation := relation + One( preA ) * arrows[ 2 * i - 1 ] * arrows[ 2 * i ] - One( preA ) * arrows[ 2 * i ] * arrows[ 2 * i - 1 ];
   od;
   vertices := VerticesOfQuiver( Qdouble );
   relations := List( vertices, v -> ( One( preA ) * v ) * relation * ( One( preA ) * v ) );
   
   return preA / relations;
end
);

#######################################################################
##
#O  AdmissibleSequenceGenerator( <num_vert, height> )
##
##  Constructs admissible sequences of Nakayama algebras with <num_vert>
##  number of vertices, and indecomposable projective modules of length
##  at most <height>. 
## 
InstallMethod( AdmissibleSequenceGenerator, 
    "for two positive integers",
    [ IsPosInt, IsPosInt ], 0,
    function( num_vert, height )
    
    local   set,  sets,  sequences,  admiss_sequences,  m,  test,  i;
    
    set := [ 2..height ];
    sets := List( [ 1..num_vert - 1 ], i -> set );
    Add( sets, [ 1..height ] );
    sequences := Cartesian( sets );
    admiss_sequences := [ ];
    for m in sequences do
        test := true;
        i := 1;
        while ( test and ( i < num_vert ) ) do
            test := ( ( m[ i + 1 ] >= m[ i ] - 1 ) and
                      ( m[ i ] - 1 >= 1 ) );  
            i := i + 1;
        od;
        if test then 
            test := ( m[ 1 ] >= m[ num_vert ] - 1 );
        fi;
        if test then 
            Add( admiss_sequences, m );
        fi;
    od;
    
    return admiss_sequences;
end
  );
