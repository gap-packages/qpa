##########################################################################
##
## DeclareOperations( "CombinatorialMap", [ size, ordering, pairing, markededge ] )
##
## Returns the combinatorial map defined by the set consisting of <size>
## elements, the ordering of the half edges given by <ordering>, the
## paring of the half edges given by the <paring>, and the marked half-edges given
## by <markededge>.  The size of the underlying set given by <size> is given
## a positive integer, the ordering and the pairing are given by two
## permutations <ordering> and <paring> of a set of <size> elements
## ({1, 2,..., <size>}) and the marked half-edges are given as a list of
## of integers in {1, 2,..., <size>}.
##

InstallMethod( CombinatorialMap,

    "for positive integer, two permutations and a list",
    
    [ IsPosInt, IsPerm, IsPerm, IsList ],
               
  function( n, sigma, iota, m )
  
  local G, combmap;
  
  G := SymmetricGroup( n );
  if not ( sigma in G ) then
    Error( "the ordering is not a permutation of proper form, \n" );
  fi;
  if not ( iota in G ) then
    Error( "the paring is not a permutation of proper form, \n" );
  fi;
  if Order( iota ) <> 2 then
    Error( "the paring is not an involution, \n" );
  fi;
  
  combmap := rec( size := n, ordering := sigma, paring := iota, markededge := m );
  ObjectifyWithAttributes( combmap, TheTypeCombinatorialMap );
  SetSize( combmap, n );
  SetOrderingOfCombinatorialMap( combmap, sigma );
  SetPairingOfCombinatorialMap( combmap, iota );
  SetMarkedEdgesOfCombinatorialMap( combmap, m);
  
  return combmap;
end
  );
  
  
######################################################
##
## DeclareAttribute( "FacesOfCombinatorialMap", map )
##
## Returns the list of faces of the combinatorial map <map>. The faces are given
## by a list of half-edges of <map>.
##

InstallMethod( FacesOfCombinatorialMap,

    "for a CombinatorialMap",
    
    [ IsCombinatorialMap ],
               
  function( combmap )

  local faces;

    faces := PairingOfCombinatorialMap( combmap ) * OrderingOfCombinatorialMap( combmap );
    
    return OrbitsPerms( [ faces ], [ 1..Size( combmap ) ] );
end
  );
  
  
########################################################
##
## DeclareAttribute( "GenusOfCombinatorialMap", map)
##
## Returns the genus of the surface represented by the combinatorial map <map>.
##
  
InstallMethod( GenusOfCombinatorialMap,

    "for a CombinatorialMap, returns the genus of the surface it represents",
    
    [ IsCombinatorialMap ],
    
    function( combmap )
    
    local faces, vertnumb, facenumb, edgenumb, euler;
    
    faces := PairingOfCombinatorialMap( combmap ) * OrderingOfCombinatorialMap( combmap );
    vertnumb := Length(OrbitsPerms([OrderingOfCombinatorialMap( combmap )],[ 1..Size( combmap ) ] )) ;
    facenumb := Length(OrbitsPerms( [ faces ], [ 1..Size( combmap ) ] ));
    edgenumb := Size(combmap)/2 ;
    
    euler := vertnumb - edgenumb + facenumb ;
    
    return ((2-euler)/2);
end
);


######################################################################
##
## DeclareAttribute( "DualOfCombinatorialMap", map)
##
## Returns the dual combinatorial map of the combinatoiral map <map>. It has the same size, pairing
## and marked half-edges as <map> and the ordering is given by the faces of <map>.
##

InstallMethod( DualOfCombinatorialMap,

    "for a CombinatorialMap, returns the dual CombinatorialMap",
    
    [ IsCombinatorialMap ],
    
    function( combmap )
    
    local faces, dual, bdots, iota, n;
    
    faces := PairingOfCombinatorialMap( combmap ) * OrderingOfCombinatorialMap( combmap );
    iota := PairingOfCombinatorialMap(combmap);
    n := Size(combmap);
    bdots := MarkedEdgesOfCombinatorialMap(combmap);
    
    dual := CombinatorialMap(n, faces, iota, bdots);
    
    return(dual);
end
);


#######################################################################
##
## DeclareAttribute( "MaximalPathsOfGentleAlgebra", algebra)
##
## Returns the list of maximal paths of the gentle algebra <algebra>.
## Only returns maximal paths of nonzero length.
##

InstallMethod( MaximalPathsOfGentleAlgebra,
    "For a GentleALgebra, returns the list of maximal paths",
    [ IsGentleAlgebra ],
    
  function(A)

  local Q, path, finished, source, Arrows, target, PrevArrow, NextArrow,maxpaths, arrow,
        beta, I,kQ,ListOfArrows,alpha1,alpha2;

    if not IsGentleAlgebra( A ) then
      Error( "The entered algebra is not gentle. \n" );
    fi;

    Q := QuiverOfPathAlgebra(A);
    kQ := OriginalPathAlgebra(A);
    I := RelationsOfAlgebra(A);
    maxpaths := [];
    Arrows := ArrowsOfQuiver(Q);


    for arrow in Arrows do
      path := arrow;
      finished := false;
      while not finished do
        finished := true ;
        source := SourceOfPath(path);
        target := TargetOfPath(path);
        PrevArrow := IncomingArrowsOfVertex(source);
        NextArrow := OutgoingArrowsOfVertex(target);
        ListOfArrows := WalkOfPath(path);
        alpha1 := ListOfArrows[1];
        alpha2 := ListOfArrows[Length(ListOfArrows)];

        for beta in PrevArrow do
          if not (ElementOfPathAlgebra(kQ,beta*alpha1) in I) then
            path := beta*path ;
            finished := false ;
            fi;
        od;

        for beta in NextArrow do
          if not (ElementOfPathAlgebra(kQ,alpha2*beta) in I) then
            path := path*beta;
            finished := false;
          fi;
        od;
      od;
    if not (path in maxpaths) then
      Append(maxpaths, [path]);
    fi;
    od;
    return maxpaths;
end
);


####################################################################
##
## DeclareOperation( "RemoveEdgeOfCombinatorialMap", [map , edge])
##
## Returns the combinatorial map obtained by removing <edge> from the combinatorial
## map <map>. The argument <edge> has to be a list consisting of two half-edges of
## <map> that are paired. The returned combinatorial map does not change size, the two
## removed half-edges are now fixed by the pairing and ordering so do not appear in non-trivial orbits.
##

InstallMethod( RemoveEdgeOfCombinatorialMap,
    "for a CombinatorialMap and a list consisting of two paired half-edges.",
    [IsCombinatorialMap, IsList],
    
    function( combmap, E )
    
    local rho, e1, e2, newrho, newiota;
    
    if Length(E) <> 2 or (not IsSubset([1..Size(combmap)],E)) then
        Error( " The entered list does not correspond to an edge. \n" );
    fi;
    
    rho := OrderingOfCombinatorialMap(combmap);
    e1 := E[1];
    e2 := E[2];
    
    if e2 = e1^rho then
        newrho := rho*(e1,e2)*(e2,e2^rho);
    elif e1 = e2^rho then
        newrho := rho*(e2,e1)*(e1,e1^rho);
    else
        newrho := rho*(e2,e2^rho)*(e1,e1^rho);
    fi;
    
    newiota := (e1,e2)*PairingOfCombinatorialMap(combmap);
    
    return CombinatorialMap(Size(combmap),newrho,newiota,[]);
    end
    );


#######################################################################################
##
## DeclareAttribute( "CombinatorialMapOfGentleAlgebra", algebra)
##
## Returns the combinatorial map corresponding to the Brauer Graph of the trivial extension
## of <algebra> with regard to its dual. This constructions only works if <algebra> is gentle.
##
    
InstallMethod( CombinatorialMapOfGentleAlgebra,
    "for a GentleAlgebra, returns the CombinatorialMap corresponding to the Brauer graph of its trivial extension.",
    
    [ IsGentleAlgebra ],
    
    function(A)
    
    local maxpaths, rho, iota, Q, n, i, V, m, w, startnumb, startvert, startind, arrow, targ,
          targvert, targind , mpoints;
    
    maxpaths := MaximalPathsOfGentleAlgebra(A);
    rho := ();
    iota := ();
    Q := QuiverOfPathAlgebra(A);
    n := Length(VerticesOfQuiver(Q));
    mpoints := [];
    
    
    for i in [1..n] do
        iota := (2*i-1,2*i)*iota;
    od;
    
    
    V := VerticesOfQuiver(Q);
    for m in maxpaths do
        w := WalkOfPath(m);
        startvert := SourceOfPath(m);
        startnumb := Position(V,startvert);
        if (2*startnumb-1)^rho = 2*startnumb-1 then
            startind := 2*startnumb-1 ;
        else
            startind := 2*startnumb;
        fi;
        
        for arrow in w do
            targ := TargetOfPath(arrow);
            targvert := Position(V,targ);
            if (2*targvert-1)^rho = 2*targvert-1 then
                targind := 2*targvert-1 ;
            else
                targind := 2*targvert;
            fi;
            rho := rho*(startind, targind) ;
        od;
        Add(mpoints,startind/rho);
    od;
    
    for i in [1..2*n] do
        if i^rho = i then
            Add(mpoints,i);
        fi;
    od;
    
    return CombinatorialMap(2*n,rho,iota,mpoints);
    
end
);


############################################################################
##
## DeclareAttribute( "MarkedBoundariesOfCombinatorialMap", map)
##
## Returns a list consisting of pairs [Face,numberofmarkedpoints] where Face
## is a face of the combinatorial map <map> and numberofmarkedpoints is the number of marked half-edges
## on the corresponding boundary component.
##

InstallMethod( MarkedBoundariesOfCombinatorialMap,
    "for a CombinatorialMap, returns a list consisting of pairs (F,n) where F is a face of the CombinatorialMap and n the number of marked points on the boundary corresponding to the face F.",
    [ IsCombinatorialMap ],
    
    function(map)
    
    local faces, F, rho, marked, numbmp, halfedge, result;
    
    faces := FacesOfCombinatorialMap(map);
    rho := OrderingOfCombinatorialMap(map);
    marked := MarkedEdgesOfCombinatorialMap(map);
    result := [];
    
    for F in faces do
       numbmp := 0;
       for halfedge in F do
            if halfedge/rho in marked then
                numbmp := numbmp + 1;
            fi;
       od;
       Add(result,[F,numbmp]);
    od;
    Sort(result);
    return result ;
end
);


##############################################################################
##
## DeclareOperation( "WindingNumber", [map,curve])
##
## Returns the combinatorial winding number of the closed curve <curve> in the
## dissected surface represented by the combinatorial map <map>.
##

InstallMethod( WindingNumber,

        "for a CombinatorialMap and a list of half-edges forming a closed curve in minimal position.",
        
        [ IsCombinatorialMap, IsList ],
        
        function(combmap, gamma)
        
        local rho, iota, marked, n, m, w, a, b, localw, i, L;
        
        rho := OrderingOfCombinatorialMap(combmap);
        iota := PairingOfCombinatorialMap(combmap);
        marked := MarkedEdgesOfCombinatorialMap(combmap);
        L := ShallowCopy(gamma);
        n := Length(L)/2;
        
        if n = 0 or (not IsList(L)) then
            Error("The chosen curve in trivial. \n");
        fi;
        
        
        if (L[1] = L[2]^iota) then
            m := Remove(L,1);
            Add(L,m);
        fi;
        
        w := 0;
        
        for i in [1..n] do
            a := L[2*i-1];
            b := L[2*i];
            if a <> b then
                localw := 1 ;
                while a <> b do
                    if a in marked then
                       localw := -1 ;
                    fi;
                    a := a^rho ;
                od;
                w := w + localw ;
            fi;
        od;
        
        return w ;

end
);


##################################################################
##
## DeclareAttribute( "BoundaryCurvesOfCombinatorialMap", map)
##
## Returns a list whose elements are lists of half-edges of the combinatorial map <map>. Each of these
## lists corresponds to a closed curve homotopic to a boundary component of the
## surface represented by <map>. The closed curves are represented by a list of adjacent half-edges.
## The orientation of these curves is chosen so that the corresponding boundary component
## is to the right. The returned curves correspond to the curves in minimal position
## in regards to the dual dissection of the surface.
## In the case of the disk, the boundary curve is trivial and the empty list is returned.
##

InstallMethod( BoundaryCurvesOfCombinatorialMap,

        " for a CombinatorialMap.",
        
        [IsCombinatorialMap],
        
        function(combmap)
        
        local faces, curves, F, halfedge, iota, c, finished, i, j, cur;
        
        faces := MarkedBoundariesOfCombinatorialMap(combmap);
        iota := PairingOfCombinatorialMap(combmap);
        
        curves := [];
        
        for F in faces do
            c := [];
            for halfedge in F[1] do
                Append(c,[halfedge,halfedge^iota]);
            od;
            Add(curves,[c,F[2]]);
        od;
        
        for j in [1..Length(curves)] do
            c := Remove(curves,1);
            cur := ShallowCopy(c[1]);
            finished := false;
            while not finished do
                finished := true;
                if cur[1] = cur[Length(cur)] then
                    Remove(cur,1);
                    Remove(cur,1);
                    Remove(cur);
                    Remove(cur);
                    finished := false;
                else
                    for i in [2..Length(cur)] do
                        if finished and cur[i-1]=cur[i] then
                            Remove(cur,i-2);
                            Remove(cur,i-2);
                            Remove(cur,i-2);
                            Remove(cur,i-2);
                            finished := false;
                        fi;
                    od;
                fi;
                if cur = [] then
                    finished := true;
                fi;
            od;
            Add(curves,[cur,c[2]]);
        od;
        
        return(curves);
        
        end
        );


#############################################################################
##
## DeclareOperation( "DepthSearchCombinatorialMap", [map,half-edge])
##
## Returns a list of pairs of paired half-edges of the combinatorial map <map> corresponding to the cover tree of
## the underlying graph of <map> obtained by a depth first search with root the vertex to
## which is attached <half-edge>. The paired half-edges of the cover tree are oriented towards
## the root of the cover tree.
##

InstallMethod( DepthSearchCombinatorialMap,

            "For a CombinatorialMap and a half-edge of the CombinatorialMap.",
            
            [IsCombinatorialMap, IsPosInt],
            
            function(combmap,x)
            
            local rho, iota, currentvertex, visited, result, pile, current, currentedge, i, vertex;
            
            rho := OrderingOfCombinatorialMap(combmap);
            iota := PairingOfCombinatorialMap(combmap);

            currentvertex := OrbitPerms([rho],x);
            Sort(currentvertex);
            visited := [];
            result := [];
            pile := [[currentvertex,x]];
            
            while pile <> [] do
                current := Remove(pile);
                currentvertex := current[1];
                currentedge := current[2];
                if not (currentvertex in visited) then
                    Add(visited,currentvertex);
                    Add(result,currentedge);
                    for i in currentvertex do
                            vertex :=OrbitPerms([rho],i^iota);
                            Sort(vertex);
                        if not(vertex in visited) then
                            Add(pile,[vertex,i^iota]);
                        fi;
                    od;
                fi;
            od;
            Remove(result,Position(result,x));
            result := OrbitsPerms([iota],result);
            
            return result;
end
);


#############################################################################
##
## DeclareOperation( "WidthSearchCombinatorialMap", [map,half-edge])
##
## Returns a list of pairs of paired half-edges of the combinatorial map <map> corresponding to the cover tree of
## the underlying graph of <map> obtained by a width first search with root the vertex to
## which is attached <half-edge>. The paired half-edges of the cover tree are oriented towards
## the root of the cover tree.
##

InstallMethod( WidthSearchCombinatorialMap,

            "For a CombinatorialMap and a half-edge of the CombinatorialMap.",
            
            [IsCombinatorialMap, IsPosInt],
            
            function(combmap,x)
            
            local rho, iota, currentvertex, visited, result, file, current, currentedge, i, vertex;
            
            rho := OrderingOfCombinatorialMap(combmap);
            iota := PairingOfCombinatorialMap(combmap);

            currentvertex := OrbitPerms([rho],x);
            Sort(currentvertex);
            visited := [];
            result := [];
            file := [[currentvertex,x]];
            
            while file <> [] do
                current := Remove(file,1);
                currentvertex := current[1];
                currentedge := current[2];
                if not (currentvertex in visited) then
                    Add(visited,currentvertex);
                    Add(result,currentedge);
                    for i in currentvertex do
                            vertex :=OrbitPerms([rho],i^iota);
                            Sort(vertex);
                        if not(vertex in visited) then
                            Add(file,[vertex,i^iota]);
                        fi;
                    od;
                fi;
            od;
            Remove(result,Position(result,x));
            result := OrbitsPerms([iota],result);
            
            return result;
end
);


##########################################################################
##
## DeclareAttribute( "NonSeperatingCurve", map)
##
## Returns a list of half-edges of the combinatorial map <map> corresponding to a non-seperating
## closed curve of the surface represented by <map>.
##


InstallMethod( NonSeperatingCurve,
            "for a CombinatorialMap of genus at least one.",
            
            [IsCombinatorialMap],
            
            function(combmap)
            
            local map, dualmap, covertree, edge, iotadual, startedge, dualtree, n, i, e, e1, e2, path1, path2, iota, rho, vertex, alpha, nonsep;
            
            map := combmap;
            dualmap := DualOfCombinatorialMap(map);
            
            covertree := WidthSearchCombinatorialMap(map,1);
            
            for edge in covertree do
                dualmap := RemoveEdgeOfCombinatorialMap(dualmap,edge);
            od;
            
            iotadual := PairingOfCombinatorialMap(dualmap);
            startedge := 1;
            
            while startedge = startedge^iotadual do
                startedge := startedge +1;
            od;
            
            dualtree := DepthSearchCombinatorialMap(dualmap,startedge);
            n := Size(dualmap);
            
            for i in [1..n/2] do
                if not ([2*i-1,2*i] in dualtree or [2*i,2*i-1] in dualtree) and (2*i)^iotadual <> (2*i) then
                    e := [2*i-1,2*i];
                fi;
            od;
            

            e1 := e[1];
            e2 := e[2];
            
            path1 := [];
            path2 := [];
            
            iota := PairingOfCombinatorialMap(combmap);
            rho := OrderingOfCombinatorialMap(combmap);
            
            
            vertex := OrbitPerms([rho],e1);
            while not (1 in vertex) do
                for alpha in vertex do
                    if [alpha,alpha^iota] in covertree then
                        Add(path1,[alpha,alpha^iota]);
                        vertex := OrbitPerms([rho],alpha^iota);
                    fi;
                od;
            od;
            
            vertex := OrbitPerms([rho],e2);
            while not (1 in vertex) do
                for alpha in vertex do
                    if [alpha,alpha^iota] in covertree then
                        Add(path2,[alpha,alpha^iota]);
                        vertex := OrbitPerms([rho],alpha^iota);
                    fi;
                od;
            od;
            
            nonsep := [];
            
            for i in [1..Length(path1)] do
                edge := Remove(path1);
                Add(nonsep,edge[2]);
                Add(nonsep,edge[1]);
            od;
            
            Append(nonsep,e);
            
            for i in [1..Length(path2)] do
                edge := Remove(path2,1);
                Append(nonsep,edge);
            od;
            
            while nonsep[1] = nonsep[Length(nonsep)] do
                Remove(nonsep);
                Remove(nonsep,1);
            od;
            
            return nonsep;
end
);


##############################################################################
##
## DeclareOperation( "CutNonSepCombinatorialMap", [map,curve,index])
##
## Returns a triplet [newmap, [boundary1,boundary2], newindex] where newmap is the combinatorial map
## obtained by cutting the surface represented by <map> along the closed curve <curve>. The lists boundary1
## and boundary2 are the boundaries of newmap that were created by the cut.
## The argument <index> is a list whose i-th element is the half-edge to which i corresponds in
## some original combinatorial map. The result newindex is the updated index for newmap where all
## created half-edges are added.
##

InstallMethod( CutNonSepCombinatorialMap,

            "for a CombinatorialMap, a nonseperating closed curve as a list and an index list for the combinatorial map.",
            
            [IsCombinatorialMap, IsList, IsList],
            
            function(combmap,gamma,index)
            
            local rho, iota, bdots, n, curve, m, l, i, sigma, tau, newmap, newborder, newindex;
            
            rho := OrderingOfCombinatorialMap(combmap);
            iota := PairingOfCombinatorialMap(combmap);
            bdots := MarkedEdgesOfCombinatorialMap(combmap);
            n := Size(combmap)+1;
            
            curve := ShallowCopy(gamma);
            
            newborder := [];
            newindex := ShallowCopy(index);
            
            if (curve[1] <> curve[2]^iota) then
                m := Remove(curve,1);
                Add(curve,m);
            fi;
            
            l := Length(curve)/2;
            
            for i in [1..(l-1)] do
                sigma := (curve[2*i], n+1, n+2);
                tau := (curve[2*i+1], n+2);
                iota := iota*(n,n+1);
                Append(newborder,[n,n+1]);
                n := n+2;
                rho := tau*rho*sigma;
            od;
            
            sigma := (curve[2*l], n+1, Size(combmap)+1);
            tau := (curve[1], Size(combmap)+1);
            iota := iota*(n,n+1);
            Append(newborder,[n,n+1]);
            n := n+1;
            rho := tau*rho*sigma;
            
            for i in curve do
                Add(newindex,index[i]);
            od;
            
            newmap := CombinatorialMap(n,rho,iota,bdots);
            
            return [newmap,[curve,newborder],newindex];
end
);


####################################################################
##
## DeclareOperation( "JoinCurveCombinatorialMap", [map, boundary1, boundary2, index])
##
## Returns a list corresponding to a simple curve of <map> joining the two boundaries <boundary1> and <boundary2>
## in such a way that using <index> it corresponds to a closed curve on the original combinatorial map.
## The list consists of adjacent half-edges of <map> which form the closed curve.
##

InstallMethod( JoinCurveCombinatorialMap,
            "For a CombinatorialMap and two lists detailing the boundaries to join, returns a list corresponding to a simple curve joining the two boundaries, starting and ending at vertices that are identified in the original CombinatorialMap.",
            [IsCombinatorialMap, IsList, IsList, IsList],
    
            function(combmap,bound1,bound2, index)
    
            local rho, iota, b1,b2,b, covertree, vertex, alpha, joincurve, m, i;
    
            rho := OrderingOfCombinatorialMap(combmap);
            iota := PairingOfCombinatorialMap(combmap);
    
            b1 := bound1[1];
    
            covertree := DepthSearchCombinatorialMap(combmap,b1);
    
            joincurve := [];
    
            for b in bound2 do
                if index[b] = index[b1] then
                    b2 := b;
                fi;
            od;
    
            vertex := OrbitPerms([rho],b2);
            while not (b1 in vertex) do
                for alpha in vertex do
                    if [alpha,alpha^iota] in covertree then
                        Add(joincurve,[alpha,alpha^iota]);
                        vertex := OrbitPerms([rho],alpha^iota);
                    fi;
                od;
            od;
    
    
    
            for i in [1..Length(joincurve)] do
                m := Remove(joincurve,1);
                Append(joincurve,m);
            od;
    
            return joincurve;
end
);


###################################################################
##
## DeclareOperation( "CutJoinCurveCombinatorialMap", [map, boundary1, boundary2, curve, index])
##
## Returns the pair [newmap, newindex] where newmap is the combinatorial map obtained
## by cutting the combinatorial map <map> along <curve> which has to be a simple
##ncurve joining <boundary1> and <boundary2>
##


InstallMethod( CutJoinCurveCombinatorialMap,
    " for a CombinatorialMap, two lists detailing boundaries, a list of arcs joining these boundaries and an index for the half-edges, returns the CombinatorialMap obtained by cutting the surface along the curve joining the two boundaries and the updated index corresponding to the cut.",
    [IsCombinatorialMap, IsList, IsList,IsList,IsList],
    
    function(combmap,bound1,bound2,curve,index)
    
    local rho, iota, bdots, n, startedge, b, l, i, sigma, tau, endedge, newindex, newmap, joincurve;
    
    rho := OrderingOfCombinatorialMap(combmap);
    iota := PairingOfCombinatorialMap(combmap);
    bdots := MarkedEdgesOfCombinatorialMap(combmap);
    n := Size(combmap)+1;
    
    joincurve := ShallowCopy(curve);
    
    while joincurve[1] in bound1 or joincurve[1] in bound2 do
        Remove(joincurve,1);
    od;
    
    while joincurve[Length(joincurve)] in bound1 or joincurve[Length(joincurve)] in bound2 do
        Remove(joincurve);
    od;
    
    startedge := joincurve[1];
    b := startedge;
    while not (b in bound1) and not (b in bound2) do
        b := b^rho;
    od;
    
    rho := (startedge,b)*(startedge,n)*rho;
    
    l := Length(joincurve)/2;
    
    for i in [1..(l-1)] do
        sigma := (joincurve[2*i], n, n+1);
        tau := (joincurve[2*i+1], n+1);
        iota := iota*(n,n+1);
        n := n+2;
        rho := tau*rho*sigma;
    od;
    
    newindex := ShallowCopy(index);
    
    for i in joincurve do
        Add(newindex,index[i]);
    od;
    
    endedge := joincurve[Length(joincurve)];
    
    b := endedge;
    while not (b in bound1) and not (b in bound2) do
        b := b^rho;
    od;
    
    rho := (n+1,b)*rho*(endedge,n+1);
    iota := iota*(n,n+1);
    
    newmap := [CombinatorialMap(n+1,rho,iota,bdots),newindex];
    
    return newmap;
    
end
);


#################################################################
##
## DeclareAttribute( "HomologyBasisOfCombinatorialMap", map)
##
## Returns a list of pairs [a_i,b_i] where a_i and b_i are lists of half-edges
## of the combinatorial map <map> forming a closed curve such that the a_i and b_i form a simplectic
## basis of the first homology group of the surface represented by <map> for the
## intersection form.
##

InstallMethod( HomologyBasisOfCombinatorialMap,

        "for a CombinatorialMap",
        
        [IsCombinatorialMap],
        
        function(combmap)
        
        local g,map,basis,index, n, alpha, newalpha, beta, newbeta, newbasis, bound1, bound2, cutresult, i, pair;
        
         g := GenusOfCombinatorialMap(combmap);
         map := ShallowCopy(combmap);
         n := Size(map);
         
         basis := [];
         index := [1..n];
         
         while g > 0 do
            alpha := NonSeperatingCurve(map);
            cutresult := CutNonSepCombinatorialMap(map,alpha,index);
            map := cutresult[1];
            bound1 := cutresult[2,1];
            bound2 := cutresult[2,2];
            index := cutresult[3];
            beta := JoinCurveCombinatorialMap(map,bound1,bound2,index);
            cutresult := CutJoinCurveCombinatorialMap(map,bound1,bound2,beta,index);
            map := cutresult[1];
            index := cutresult[2];
            Add(basis,[alpha,beta]);
            g := GenusOfCombinatorialMap(map);
         od;
         
         newbasis := [];
         
        for pair in basis do
            alpha := pair[1];
            beta := pair[2];
            newalpha := [];
            newbeta := [];
            for i in alpha do
                Add(newalpha,index[i]);
            od;
            for i in beta do
                Add(newbeta,index[i]);
            od;
            Add(newbasis,[newalpha,newbeta]);
         od;
         
        return newbasis;
end
);


#########################################################
##
## DeclareOperation( "AreDerivedEquivalent", [algebra1,algebra2])
##
## Returns true if <algebra1> and <algebra2> are derived equivalent and false otherwise.
## The arguments must be gentle algebras over the same field.
##

InstallMethod( AreDerivedEquivalent ,
    "for two GentleAlgebras, returns true if they are derived equivalent and false if they are not." ,
    [IsQuiverAlgebra , IsQuiverAlgebra],
    
    function(A,B)
    
    local mapA, mapB, gA, gB, bA, bB, boundA, boundB, c, curve, i, wcA, wcB, hombasisA, hombasisB, wnumberBasisA, wnumberBasisB, pair, gcdlistA, gcdlistB, w, oddA, oddB, mod4, arfA, arfB ;
    
    if not (IsGentleAlgebra(A) and IsGentleAlgebra(B)) then
        Error( "At least one of the algebras is not gentle? \n");
    fi;
    
    if LeftActingDomain(A) <> LeftActingDomain(B) then
        return false;
    fi;
    
    mapA := CombinatorialMapOfGentleAlgebra(A);
    mapB := CombinatorialMapOfGentleAlgebra(B);
    
    gA := GenusOfCombinatorialMap(mapA);
    gB := GenusOfCombinatorialMap(mapB);
    
    bA := FacesOfCombinatorialMap(mapA);
    bB := FacesOfCombinatorialMap(mapB);
    
    if gA <> gB then
        return false;
    elif Length(bA) <> Length(bB) then
        return false;
    elif Length(OrbitsPerms([OrderingOfCombinatorialMap(mapA)],[1..Size(mapA)])) <> Length(OrbitsPerms([OrderingOfCombinatorialMap(mapB)],[1..Size(mapB)])) then
        return false;
    fi;
    
    boundA := ShallowCopy(BoundaryCurvesOfCombinatorialMap(mapA));
    boundB := ShallowCopy(BoundaryCurvesOfCombinatorialMap(mapB));
    
    if gA =0 and Length(bA) = 1 then
        return true;
    fi;
    
    for i in [1..Length(boundA)] do
        c := Remove(boundA,1);
        curve := c[1];
        Add(boundA,[-WindingNumber(mapA,curve),c[2]]);
    od;
    
    for i in [1..Length(boundB)] do
        c := Remove(boundB,1);
        curve := c[1];
        Add(boundB,[-WindingNumber(mapB,curve),c[2]]);
    od;
    
    if not(IsSubset(boundA,boundB)) then
        return false;
    fi;
    
    if gA = 0 then
        return true;
    fi;
    
    wcA := [];
    wcB := [];
    
    for i in [1..Length(boundA)] do
        Add(wcA,boundA[i,1]);
        Add(wcB,boundB[i,1]);
    od;
    
    hombasisA := HomologyBasisOfCombinatorialMap(mapA);
    hombasisB := HomologyBasisOfCombinatorialMap(mapB);
    
    wnumberBasisA := [];
    wnumberBasisB := [];
    
    for pair in hombasisA do
        Add(wnumberBasisA,[WindingNumber(mapA,pair[1]),WindingNumber(mapA,pair[2])]);
    od;
    
    for pair in hombasisB do
        Add(wnumberBasisB,[WindingNumber(mapB,pair[1]),WindingNumber(mapB,pair[2])]);
    od;
    
    if gA = 1 then
        gcdlistA := [];
        gcdlistB := [];
        
        for w in wcA do
            Add(gcdlistA,w+2);
        od;
        
        for w in wcB do
            Add(gcdlistB,w+2);
        od;
        
        for pair in wnumberBasisA do
            Append(gcdlistA,pair);
        od;
        
        for pair in wnumberBasisB do
            Append(gcdlistB,pair);
        od;
        
        if AbsInt(Gcd(Integers, gcdlistA)) = AbsInt(Gcd(Integers,gcdlistB)) then
            return true;
        else
            return false;
        fi;
    else
        oddA := false;
        oddB := false;
        
        mod4 := false;
        
        for w in wcA do
            if IsOddInt(w) then
                oddA := true;
            fi;
            if RemInt(w,4) = 0 then
                mod4 := true;
            fi;
        od;
        
        for pair in wnumberBasisA do
            if IsOddInt( pair[1] ) or IsOddInt( pair[2] ) then
                oddA := true;
            fi;
        od;
        
        for w in wcB do
            if IsOddInt(w) then
                oddB := true;
            fi;
        od;
        
        for pair in wnumberBasisA do
            if IsOddInt( pair[1] ) or IsOddInt( pair[2] ) then
                oddA := true;
            fi;
        od;
        
        if oddA and oddB then
            return true;
        elif mod4 then
            return true;
        fi;
        
        arfA := 0;
        arfB := 0;
        
        for pair in wnumberBasisA do
            arfA := arfA + (pair[1]/2 + 1)*(pair[2]/2 + 1);
        od;
        
        for pair in wnumberBasisB do
            arfB := arfB + (pair[1]/2 + 1)*(pair[2]/2 + 1);
        od;
        
        if RemInt(arfA-arfB,2) = 0 then
            return true;
        else
            return false;
        fi;
    fi;
end
);
