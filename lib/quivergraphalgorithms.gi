# GAP Implementation
# This contains various graph algorithms for testing quiver's properties
# (created A. Mroz, 21.02.2013)


#############################################################################
##
#P  IsAcyclicQuiver(<Q>)
##
##  This function returns true if a quiver <Q>  does not 
##  contain an oriented cycle. 
##
InstallMethod( IsAcyclicQuiver,
  "for quivers",
  true,
  [ IsQuiver ], 0,
  function ( Q )
    local Visit, color, GRAY, BLACK, WHITE,
          vert, vertex_list, res, tsorted;
    
    WHITE := 0; GRAY := 1; BLACK := -1;
    
    tsorted := [];

    Visit := function(v)
      local adj, uPos, result; # adjacent vertices
        
      color[v] := GRAY;
        
      adj := List(NeighborsOfVertex(vertex_list[v]),
                  x -> Position(vertex_list, x));
      Info(InfoQuiver, 1, "Adjacent vertices are:", adj);        

      if not IsEmpty(adj) and ForAny(adj, x -> color[x] = GRAY) then
        Info(InfoQuiver, 1, "Found a GRAY vertex.");
        return false;
      fi;

      for uPos in adj do
        if color[uPos] = WHITE then
          result := Visit( uPos );
          if not result then
            return false;
          fi;
        fi;
      od;

      Add(tsorted, vertex_list[v]);
      color[v] := BLACK;
      return true;
    end;
    
    color := [];
    vertex_list := VerticesOfQuiver(Q);

    for vert in [1 .. Length(vertex_list) ] do
      color[vert] := WHITE;
    od;
    
    for vert in [1 .. Length(vertex_list) ] do
      if color[vert] = WHITE then
        res := Visit(vert);
        if not res then
          return false;
        fi;
      fi;
    od;
    
    tsorted := Reversed(tsorted);
    Q!.tsorted := tsorted;

    return true;

  end
); # IsAcyclicQuiver



#############################################################################
##
#P  IsConnectedQuiver(<Q>)
##
##  This function returns true if a quiver <Q> is a connected graph
##  (i.e. each pair of vertices is connected by an unoriented path).
##
InstallMethod( IsConnectedQuiver,
  "for quivers",
  true,
  [ IsQuiver ], 0,
  function ( Q )
    local  Visit, color, GRAY, BLACK, WHITE,
          vertex_list, res, vert;
    
    WHITE := 0; GRAY := 1; BLACK := -1;
    
    Visit := function(v)
      local adj, uPos, result, vertex, arrow; 
        
      color[v] := GRAY;
        
      adj := List(NeighborsOfVertex(vertex_list[v]),
                  x -> Position(vertex_list, x));
      for arrow in IncomingArrowsOfVertex(vertex_list[v]) do
        vertex := Position(vertex_list, SourceOfPath(arrow));
        if not (vertex in adj) then
          Add(adj, vertex);
        fi;
      od; # Now adj contains all vertices adjacent to v (by arrows outgoing from and incoming to v)
      
      for vertex in adj do
        if color[vertex] = WHITE then
          Visit(vertex);
        fi;
      od;
      color[v] := BLACK;
      
      return true;
    end;
    
    color := [];
    vertex_list := VerticesOfQuiver(Q);
    
    if Length(vertex_list) = 0 then
      return true;
    fi;
    
    for vert in [1 .. Length(vertex_list) ] do
      color[vert] := WHITE;
    od;
    
    Visit(1);
    
    return not (WHITE in color);

  end
); #IsConnectedQuiver




#############################################################################
##
#P  IsUAcyclicQuiver(<Q>)
##
##  This function returns true if a quiver <Q>  does not 
##  contain an unoriented cycle. 
##  Note: an oriented cycle is also an unoriented cycle
##
InstallMethod( IsUAcyclicQuiver,
  "for quivers",
  true,
  [ IsQuiver ], 0,
  function ( Q )
    local Visit, color, GRAY, BLACK, WHITE,
          vertex_list, res, vert, i;
    
    WHITE := 0; GRAY := 1; BLACK := -1;
    
    
    Visit := function(v, predecessor)
      local adj, result, vertex, arrow; 
        
      color[v] := GRAY;
        
      adj := [];
      for arrow in IncomingArrowsOfVertex(vertex_list[v]) do
        vertex := Position(vertex_list, SourceOfPath(arrow));
        if vertex in adj then
          return false; # (un)oriented cycle of length 1 or 2 detected!
        fi;
        Add(adj, vertex);
      od;
      for arrow in OutgoingArrowsOfVertex(vertex_list[v]) do
        vertex := Position(vertex_list, TargetOfPath(arrow));
        if vertex in adj then
          return false; # (un)oriented cycle of length 1 or 2 detected!
        fi;
        Add(adj, vertex);
      od; 
      # Now adj contains all vertices adjacent to v by arrows incoming to v and outgoing from v 
      
      for vertex in adj do
        if (color[vertex] = GRAY) and (predecessor <> vertex) then
          #cycle detected!
          return false;
        fi;
        if color[vertex] = WHITE then
          if not Visit(vertex, v) then
            return false; 
          fi;
        fi;
      od;
      
      color[v] := BLACK;
      
      return true;
    end;
    
    color := [];
    vertex_list := VerticesOfQuiver(Q);
    
    for i in [1 .. Length(vertex_list) ] do
      color[i] := WHITE;
    od;
     
    for vert in [1 .. Length(vertex_list) ] do
      if color[vert] = WHITE then
        if not Visit(vert, -1) then
          return false;
        fi;
      fi;
    od;
    
    return true;

  end
); #IsUAcyclicQuiver










#############################################################################
##
#P  IsTreeQuiver(<Q>)
##
##  This function returns true if a quiver <Q> is a tree as a graph
##  (i.e. it is connected and contains no unoriented cycles).
##
InstallMethod( IsTreeQuiver,
  "for quivers",
  true,
  [ IsQuiver ], 0,
  function ( Q )
    
    return (NumberOfArrows(Q) = NumberOfVertices(Q) - 1) and IsConnectedQuiver(Q);

  end
); #IsTreeQuiver



#############################################################################
##
#P  IsDynkinQuiver(<Q>)
##
##  This function returns true if a quiver <Q>  is a Dynkin quiver,
##  i.e. an underlying graph is a Dynkin diagram.
##  Note: function works only for connected quivers.
##
InstallMethod( IsDynkinQuiver,
  "for quivers",
  true,
  [ IsQuiver ], 0,
  function ( Q )
  
    local arr, v, vertices, degrees, deg, fork_vertex, 
          fork_arrows, GoThrough, star_type;
  
    # uaxiliary function computing the max. length of an unoriented path
    # starting from vertex s, going through the arrow
    # It makes sense only in this particular local situation of a "star" quiver 
    GoThrough := function(s, arrow)
      local temp_s, temp_arr, can_go, counter, arrows;
      
      counter := 0;
      can_go := true;
      temp_s := s;
      temp_arr := arrow;
      while can_go do
        if TargetOfPath(temp_arr) <> temp_s then
          temp_s := TargetOfPath(temp_arr);
        else
          temp_s := SourceOfPath(temp_arr);
        fi;
        counter := counter + 1;
        arrows := IncomingArrowsOfVertex(temp_s);
        Append(arrows, OutgoingArrowsOfVertex(temp_s));
        arrows := Difference(arrows, [temp_arr]);
        if Length(arrows) = 1 then
          temp_arr := arrows[1];
        else
          can_go := false;
        fi;
      od;
      
      return counter;
    end; # GoThrough
    
    if not IsConnectedQuiver(Q) then
      Print("Quiver is not connected.\n");
      return false;
    fi;
    if not (NumberOfArrows(Q) = NumberOfVertices(Q) - 1) then
      Print("Quiver contains an (un)oriented cycle.\n");
      return false;
    fi;
    
    vertices := VerticesOfQuiver(Q);
    degrees := [];
    
    # collecting info on degrees of vertices
    for v in vertices do
      deg := InDegreeOfVertex(v) + OutDegreeOfVertex(v);
      if deg > 2 then
        fork_vertex := v;
        Add(degrees, deg);
      fi;
    od;
    
    if Length(degrees) = 0 then
      # No vertex of degree > 2 and Q is tree => Q is  A_n 
      Print("A_",Length(vertices),"\n");
      return true;
    else # TESTING FOR REMAINING DYNKINS
    
      if degrees = [3] then
      # necessary condition for being D_n, E_6,7,8 satisfied (
      # (i.e. we have a "star" quiver with 3 arms)  
        fork_arrows := IncomingArrowsOfVertex(fork_vertex);
        Append(fork_arrows, OutgoingArrowsOfVertex(fork_vertex));
        star_type := [];
        for arr in fork_arrows do
          Add(star_type, GoThrough(fork_vertex, arr));
        od;
        Sort(star_type);
        #Print(star_type,"\n");
        
        if star_type{[1,2]} = [1,1] then
          Print("D_",Length(vertices),"\n");
          return true;
        fi;
        if star_type in [ [1,2,2], [1,2,3], [1,2,4] ] then
          Print("E_",Length(vertices),"\n");  
          return true;
        fi;
        
      fi;
      
    fi;
    
    return false;
    
    
    
    end
); #IsDynkinQuiver




######################################################
##
#O FullSubquiver( <Q> , <L>)
##
## This function returns a quiver which is a fullsubquiver
## of a quiver <Q> induced by the list of vertices <L>.
## The names of vertices and arrows in resulting (sub)quiver 
## remain the same as in original one.
##
InstallMethod( FullSubquiver,
               
  [ IsQuiver, IsList ],
  function( Q, L )
    local arr, vertices, arrows, result, src, trg;
    
    if not ForAll(L, x -> (IsQuiverVertex(x) and (x in Q)) ) then
      Error("L should be a list of vertices of Q!");
    fi;
    
    vertices := List(L, String);
    
    arrows := [];
    for arr in ArrowsOfQuiver(Q) do
      src := SourceOfPath(arr);
      trg := TargetOfPath(arr);
      if (src in L) and (trg in L) then
        Add(arrows, [ String(src),
                      String(trg),
                      String(arr) ] );
      fi;
    od;
    
    result := Quiver(vertices, arrows);
    
    return result;
  
  end
); # FullSubquiver




######################################################
##
#O ConnectedComponentsOfQuiver( <Q> )
##
## This function returns a list of quivers which are 
##  all connected components a quiver <Q>.
##
InstallMethod( ConnectedComponentsOfQuiver,
               
  [ IsQuiver ],
  function( Q )
    local  Visit, color, GRAY, BLACK, WHITE,
          vertex_list, components, comp_q, vert, comp_verts;
    
    WHITE := 0; GRAY := 1; BLACK := -1;
    
    Visit := function(v)
      local adj, result, vertex, arrow; 
        
      color[v] := GRAY;
        
      adj := List(NeighborsOfVertex(vertex_list[v]),
                  x -> Position(vertex_list, x));
      for arrow in IncomingArrowsOfVertex(vertex_list[v]) do
        vertex := Position(vertex_list, SourceOfPath(arrow));
        if not (vertex in adj) then
          Add(adj, vertex);
        fi;
      od; # Now adj contains all vertices adjacent to v (by arrows outgoing from and incoming to v)
      
      for vertex in adj do
        if color[vertex] = WHITE then
          Visit(vertex);
        fi;
      od;
      color[v] := BLACK;
      Add(comp_verts, vertex_list[v]);
      
      return true;
    end;
    
    color := [];
    vertex_list := VerticesOfQuiver(Q);
    components := [];
   
    for vert in [1 .. Length(vertex_list) ] do
      color[vert] := WHITE;
    od;
    
    for vert in [1 .. Length(vertex_list) ] do
      if color[vert] = WHITE then
        comp_verts := [];
        Visit(vert);
        comp_q := FullSubquiver(Q, comp_verts);
        SetIsConnectedQuiver(comp_q, true);
        Add(components, comp_q );
      fi;
    od;
    
    return components;
  end
); # ConnectedComponentsOfQuiver