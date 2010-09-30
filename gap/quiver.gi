# GAP Implementation
# This file was generated from
# $Id: quiver.gi,v 1.2 2010/09/30 14:02:22 oysteini Exp $

InstallGlobalFunction(
  Path,
  function( arg )
    local path, vertex_name, u, v, arrow_name, arrow_list;

    # single vertex name
    if Length( arg ) = 2 and IsFamily( arg[1] ) and IsString( arg[2] ) then
      vertex_name := arg[2];
      path:= Objectify( NewType( arg[1],
                        IsVertex and IsVertexRep and IsAttributeStoringRep ),
                        rec( vertex_name := vertex_name ) );
      SetSourceOfPath( path, path );
      SetTargetOfPath( path, path );
      SetLengthOfPath( path, 0 );
      SetWalkOfPath( path, [] );
      SetIsZeroPath( path, false );
      SetIncomingArrowsOfVertex( path, [] );
      SetOutgoingArrowsOfVertex( path, [] );
      SetOutDegreeOfVertex( path, 0 );
      SetInDegreeOfVertex( path, 0 );
      SetNeighborsOfVertex( path, [] );

      return path;

    # arrow from two vertices
    elif Length( arg ) = 4 and IsFamily( arg[1]) and IsVertex( arg[2] )
                         and IsVertex( arg[3] ) and IsString( arg[4] ) then
      u := arg[2];
      v := arg[3];
      arrow_name := arg[4];
      path:= Objectify( NewType( arg[1], 
                        IsArrow and IsArrowRep and IsAttributeStoringRep ),
                        rec( arrow_name := arrow_name ) );
      SetSourceOfPath( path, u );
      SetTargetOfPath( path, v );
      SetLengthOfPath( path, 1 );
      SetWalkOfPath( path, [ path ] );
      SetIsZeroPath( path, false );
      SetOutgoingArrowsOfVertex( u,
              Concatenation( OutgoingArrowsOfVertex(u), [ path ] ) );
      SetIncomingArrowsOfVertex( v,
              Concatenation( IncomingArrowsOfVertex(v), [ path ] ) );
      SetOutDegreeOfVertex( u, Length( OutgoingArrowsOfVertex(u) ) );
      SetInDegreeOfVertex(  v, Length( IncomingArrowsOfVertex(v) ) );

      # Don't add neighbor more than once.
      if not (v in NeighborsOfVertex(u)) then
          SetNeighborsOfVertex( u, Concatenation(NeighborsOfVertex(u), [ v ]) );
      fi;

      return path;

    # list of arrows
    elif Length( arg ) = 2 and IsFamily( arg[1] ) and IsList( arg[2] ) then
      arrow_list := arg[2];

      return Product( arrow_list );

    fi;

    # no argument given, error
    Error("usage: Path(<Fam>, <vertex_name>), \
       (<Fam>, <U>, <V>, <arrow_name>), (<Fam>, < [ A1, A2, ... ] >)");

  end
);


InstallMethod(ExtRepOfObj, 
  "for zero paths",
  true, [IsZeroPath], 0,
  P -> [0]
);


InstallMethod(ExtRepOfObj,
  "for vertices",
  true, [IsVertex], 0,
  V -> [V!.gen_pos]
);


InstallMethod(ExtRepOfObj,
  "for paths",
  true, [IsPath], 0,
  function(P)
    local g, rep;
    rep := [];
    for g in WalkOfPath(P) do
      Add(rep, g!.gen_pos);
    od;
    return rep;
  end
);


InstallMethod( ObjByExtRep,
  "for quiver elements family and lists",
  true,
  [IsPathFamily, IsList], 0,
  function( Fam, descr )
    local Q, gens, walk;

    Q := Fam!.quiver;
    gens := GeneratorsOfQuiver(Q);

    walk := List( descr, function(i)
      if i <> 0 then
        return gens[i];
      else
        return Zero(Q);
      fi;
    end );

    return Product(walk);
  end
);


InstallMethod( \*,
  "for vertices",
  IsIdenticalObj,
  [ IsVertex, IsVertex ], 0,
  function( u, v )
    if( IsIdenticalObj(u, v) ) then
      return u;
    else
      return Zero(FamilyObj(u));
    fi;
  end
);


InstallMethod( \*,
  "for a vertex and a path",
  IsIdenticalObj,
  [ IsVertex, IsPath ], 0,
  function( vert, path )
    if( IsIdenticalObj(vert, SourceOfPath(path)) ) then
      return path;
    else
      return Zero(FamilyObj(vert));
    fi;
  end
);


InstallMethod( \*,
  "for a path and a vertex",
  IsIdenticalObj,
  [ IsPath, IsVertex ], 0,
  function( path, vert )
    if( IsIdenticalObj(TargetOfPath(path), vert) ) then
      return path;
    else
      return Zero(FamilyObj(path));
    fi;
  end
);


InstallMethod( \*,
  "for the ZeroPath and a path",
  IsIdenticalObj,
  [ IsZeroPath, IsPath ], 0,
  function( zero, path )
    return zero;
  end
);


InstallMethod( \*,
  "for a path and the ZeroPath",
  IsIdenticalObj,
  [ IsPath, IsZeroPath ], 0,
  function( path, zero )
    return zero;
  end
);


InstallMethod( \*,
  "for a path and a path",
  IsIdenticalObj,
  [ IsPath, IsPath ], 0,
  function( path1, path2 )

    local result;

    if( IsIdenticalObj(TargetOfPath(path1), SourceOfPath(path2) ) ) then

      result:= Objectify( NewType( FamilyObj(path1),
                                   IsPath and IsAttributeStoringRep ),
                          rec() );
      SetSourceOfPath( result, SourceOfPath(path1) );
      SetTargetOfPath( result, TargetOfPath(path2) );
      SetLengthOfPath( result, LengthOfPath(path1) + LengthOfPath(path2) );
      SetWalkOfPath( result,
                     Concatenation( WalkOfPath(path1), WalkOfPath(path2) ) );
      SetIsZeroPath( result, false );

      return result;

    else

      return Zero(FamilyObj(path1));

    fi;
  end
);


InstallMethod( \=,
  "for vertices",
  IsIdenticalObj,
  [ IsVertex, IsVertex ], 0,
  function( u, v )
    return IsIdenticalObj(u, v);
  end
);


InstallMethod( \=,
  "for arrows",
  IsIdenticalObj,
  [ IsArrow, IsArrow ], 0,
  function( a, b )
    return IsIdenticalObj(a, b);
  end
);


InstallMethod( \=,
  "for the ZeroPath and the ZeroPath",
  IsIdenticalObj,
  [ IsZeroPath, IsZeroPath ], 0,
  function( a, b )
    return true;
  end
);


InstallMethod( \=,
  "for the ZeroPath and a path",
  IsIdenticalObj,
  [ IsZeroPath, IsPath ], 0,
  function( a, b )
    return false;
  end
);


InstallMethod( \=,
  "for a path and the ZeroPath",
  IsIdenticalObj,
  [ IsPath, IsZeroPath ], 0,
  function( a, b )
    return false;
  end
);


InstallMethod( \=,
  "for path",
  IsIdenticalObj,
  [ IsPath, IsPath ], 0,
  function( p1, p2 )

    if( not IsIdenticalObj(SourceOfPath(p1) , SourceOfPath(p2)) ) then
      return false;

    elif( not IsIdenticalObj(SourceOfPath(p1), SourceOfPath(p2)) ) then
      return false;

    elif( LengthOfPath(p1) <> LengthOfPath(p2)) then
      return false;

    else
      return WalkOfPath(p1) = WalkOfPath(p2);

    fi;

  end
);


InstallMethod( \<, "for paths",
  IsIdenticalObj,
  [IsPath, IsPath], 0,
  function(p, q)

    local O;

    if IsZeroPath(p) and not IsZeroPath(q) then
      return true;

    elif (IsZeroPath(q) and not IsZeroPath(p)) or (p = q) then
      return false;

    else
      O := OrderingOfQuiver(FamilyObj(p)!.quiver);
      return LessThanByOrdering(O, p, q);
    fi;

  end
);


InstallMethod( PrintObj,
  "for ZeroPath",
  true,
  [ IsZeroPath ], 0,
  function( obj )
    Print( "0" );;
  end
);


InstallMethod( String,
  "for vertex",
  true,
  [ IsVertex and IsVertexRep ], 0,
  function( obj )
    return obj!.vertex_name;
  end
);


InstallMethod( PrintObj,
  "for vertex",
  true,
  [ IsVertex and IsVertexRep ], 0,
  function( obj )
    Print( obj!.vertex_name );
  end
);


InstallMethod( String,
  "for arrow",
  true,
  [ IsArrow and IsArrowRep ], 0,
  function( obj )
    return obj!.arrow_name;
  end
);


InstallMethod( PrintObj,
  "for arrow",
  true,
  [ IsArrow and IsArrowRep ], 0,
  function( obj )
    Print( obj!.arrow_name );
  end
);


InstallMethod( PrintObj,
  "for paths of length > 1",
  true,
  [ IsPath ], 0,
  function( obj )

    local walk, i, l, exponent;

    walk := WalkOfPath(obj);
    l := Length(walk);

    Print(walk[1]);
    exponent := 1;

    for i in [2..l] do
      if walk[i] <> walk[i-1] then
        if exponent > 1 then
          Print("^", exponent);
        fi;
        Print("*", walk[i]);
        exponent := 1;
      else
        exponent := exponent + 1;
      fi;
    od;

    if exponent > 1 then
      Print("^", exponent);
    fi;

  end
);


InstallGlobalFunction( Quiver,
  function( arg )

    local vertices, arrows, vertices_by_name, arrow_spec_size,
          name, u, v, i, j, k, msg, arrow_count, Q,
          matrix, record, Fam, zero, frompos, topos;

    vertices := [];
    arrows := [];
    
    # Create a new family for paths in this quiver.
    Fam := NewFamily( Concatenation( "FamilyOfPathsWithin", UniqueQuiverName() ), IsPath ); 
    SetFilterObj( Fam, IsPathFamily );
    zero := Objectify( NewType( Fam, IsPath and IsAttributeStoringRep ), rec() );
    SetIsZeroPath( zero, true );
    SetSourceOfPath( zero, zero );
    SetTargetOfPath( zero, zero );
    SetZero(Fam, zero);

    # Quiver(N, [arrow_spec])
    if Length( arg ) = 2 and IsPosInt( arg[1] ) and IsList( arg[2] ) then
      matrix := [];

      for i in [ 1 .. arg[1] ] do
        name := Concatenation( "v", String(i) );
        vertices[i] := Path( Fam, name );
        matrix[i] := [];
        for j in [ 1 .. arg[1] ] do
          matrix[i][j] := 0;
        od;
      od;
            
      if( Length( arg[2] ) > 0 ) then
        arrow_spec_size := Length( arg[2][1] );
      else
        arrow_spec_size := -1;
      fi;

      for i in [ 1 .. Length(arg[2]) ] do

        if( not Length( arg[2][i] ) = arrow_spec_size ) then
          Error("All of the entries in the [arrow_spec] list must be of the same size.");
          return 0;
        fi;

        if( not IsPosInt( arg[2][i][1] ) or not IsPosInt( arg[2][i][2]) ) then
          Error("You must use vertex numbers in arrow_spec \
                 when using Quiver( <N>, <Arrow_spec> ).");
          return 0;
        fi;

        frompos := arg[2][i][1];
        topos := arg[2][i][2];
        matrix[frompos][topos] := matrix[frompos][topos] + 1;
        u := vertices[ frompos ];
        v := vertices[ topos ];

        if( arrow_spec_size = 3 ) then
          name := arg[2][i][3];
        else
          name := Concatenation( "a", String(i) );
        fi;

        arrows[i] := Path( Fam, u, v, name );
      od;

    # Quiver([vertex_name, ...], [arrow_spec, ...])
    elif Length( arg ) = 2 and IsList( arg[1] ) and IsList( arg[2] ) then
      matrix := [];

      for i in [ 1 .. Length( arg[1] ) ] do

        vertices[i] := Path( Fam, arg[1][i] );
        matrix[i] := [];

        for j in [ 1 .. Length( arg[1] ) ] do
          matrix[i][j] := 0;
        od;

      od;
            
      if( Length( arg[2] ) > 0 ) then
        arrow_spec_size := Length( arg[2][1] );
        vertices_by_name := not IsPosInt( arg[2][1][1] );
      fi;

      for i in [ 1 .. Length(arg[2]) ] do
        if( not Length( arg[2][i] ) = arrow_spec_size ) then
            Error("All of the entries in the [arrow_spec] list must be of the same size.");
            return 0;
        fi;

        if( vertices_by_name ) then
          j := Position( arg[1], arg[2][i][1] );

          if(j = fail) then
            msg := "Cannot find vertex: ";
            Append( msg, arg[2][i][1]);
            Append( msg, " in the list of vertices:\n");
            Append( msg, String( arg[1] ) );
            Error( msg );
            return 0;
          fi;

          frompos := j;
          u := vertices[ j ];
          j := Position( arg[1], arg[2][i][2] );

          if(j = fail) then
            msg := "Cannot find vertex: ";
            Append( msg, arg[2][i][2]);
            Append( msg, " in the list of vertices:\n");
            Append( msg, String( arg[1] ) );
            Error( msg );
            return 0;
          fi;

          topos := j;
          v := vertices[ j ];
          matrix[frompos][topos] := matrix[frompos][topos] + 1;

        else

          frompos := arg[2][i][1];
          topos := arg[2][i][2];
          u := vertices[ frompos ];
          v := vertices[ topos ];
          matrix[frompos][topos] := matrix[frompos][topos] + 1;

        fi;
         
        if( arrow_spec_size = 3 ) then
          name := arg[2][i][3];
        else
          name := Concatenation( "a", String(i) );
        fi;

        arrows[i] := Path( Fam, u, v, name );
      od;

    # Quiver(adj_matrix)
    elif Length( arg ) = 1 and IsMatrix( arg[1] ) then
      matrix := arg[1];

      for i in [ 1 .. Length( matrix ) ] do
        name := Concatenation( "v", String(i) );
        vertices[i] := Path( Fam, name );
      od;

      arrow_count := 0;

      for i in [ 1 .. Length( matrix ) ] do
        if( not Length( matrix ) = Length( matrix[i] ) ) then
            Error("The adjacency matrix must be square.");
            return 0;
        fi;

        for j in [ 1 .. Length( matrix[i] ) ] do
          for k in [ 1 .. matrix[i][j] ] do
            arrow_count := arrow_count + 1;
            u := vertices[i];
            v := vertices[j];
            name := Concatenation( "a", String(arrow_count) );
            arrows[arrow_count] := Path( Fam, u, v, name );
          od;
        od;

      od;
    
    else
      # no argument given, error
      Error("usage: Quiver(<N>, <[arrow_spec, ...]>), \
            (<[vertex_name, ...]>, [arrow_spec, ...]>), (<adj_matrix>)");
    fi;

    record := rec();

    for i in vertices do
      record.( String(i) ) := i;
    od;

    for i in arrows do
      record.( String(i) ) := i;
    od;

    Q:= Objectify( NewType( CollectionsFamily(Fam), 
                            IsQuiverRep and IsAttributeStoringRep ),
                   rec( pieces := record ) );

    Fam!.quiver := Q;
    SetVerticesOfQuiver( Q, vertices );
    SetArrowsOfQuiver( Q, arrows );
    SetOrderOfQuiver( Q, Length( vertices ) );
    SetSizeOfQuiver( Q, Length( arrows ) );
    SetGeneratorsOfMagma( Q, Concatenation( vertices, arrows ) );
    SetAdjacencyMatrixOfQuiver( Q, matrix );
    SetIsWholeFamily(Q, true);
    SetIsAssociative(Q, true);

    # This is used to construct the external representations of elements
    for i in [1..Length(GeneratorsOfMagma(Q))] do
      GeneratorsOfMagma(Q)[i]!.gen_pos := i;
    od;

    SetOrderingOfQuiver( Q, 
      LengthOrdering(Q,LeftLexicographicOrdering(Q,GeneratorsOfMagma( Q )))
    );

    SetZero(Q, zero);

    return Q;

  end
);


InstallMethod( \.,
  "for a quiver",
  true,
  [ IsQuiver, IsPosInt ], 0,
  function( Q, nam )
    return Q!.pieces.( NameRNam(nam)  );
  end
);


GlobalQuiverCount := 0;


InstallGlobalFunction( UniqueQuiverName,
  function( arg )
    GlobalQuiverCount := GlobalQuiverCount + 1;
    return Concatenation( "Quiver_", String(GlobalQuiverCount) );
  end
);


InstallMethod( OrderedBy, 
  "for a quiver and an ordering",
  true,
  [IsQuiver, IsQuiverOrdering], 0,
  function(Q, O)
    local Fam, new_vertices, new_arrows, new_Q, v, a, new_v, new_a,
          record, zero;
    
    # Create a new family for paths in this quiver.
    Fam := NewFamily( Concatenation( "FamilyOfPathsWithin", UniqueQuiverName() ), IsPath ); 
    SetFilterObj( Fam, IsPathFamily );
    zero := Objectify( NewType( Fam, IsPath and IsAttributeStoringRep ), rec() );
    SetIsZeroPath( zero, true );
    SetSourceOfPath( zero, zero );
    SetTargetOfPath( zero, zero );
    SetZero(Fam, zero);

    new_vertices := [];

    for v in VerticesOfQuiver(Q) do
      new_v := Path(Fam, v!.vertex_name);
      new_v!.gen_pos := v!.gen_pos;
      Add(new_vertices, new_v);
    od;

    new_arrows := [];

    for a in ArrowsOfQuiver(Q) do
      new_a := Path(Fam, new_vertices[SourceOfPath(a)!.gen_pos],
                         new_vertices[TargetOfPath(a)!.gen_pos],
                         a!.arrow_name);
      new_a!.gen_pos := a!.gen_pos;
      Add(new_arrows, new_a);
    od;

    record := rec();

    for v in new_vertices do
      record.( String(v) ) := v;
    od;

    for a in new_arrows do
      record.( String(a) ) := a;
    od;

    new_Q:= Objectify( NewType( CollectionsFamily(Fam), 
                                IsQuiverRep and IsAttributeStoringRep ),
                       rec( pieces := record ) );

    Fam!.quiver := new_Q;

    SetVerticesOfQuiver( new_Q, new_vertices );
    SetArrowsOfQuiver( new_Q, new_arrows );
    SetOrderOfQuiver( new_Q, Length( new_vertices ) );
    SetSizeOfQuiver( new_Q, Length( new_arrows ) );
    SetGeneratorsOfMagma( new_Q, Concatenation( new_vertices, new_arrows ) );
    SetAdjacencyMatrixOfQuiver( new_Q, AdjacencyMatrixOfQuiver(Q) );
    SetOrderingOfQuiver( new_Q, O ); 
    SetZero(new_Q, zero);
    SetIsWholeFamily(new_Q, true);

    return new_Q;

  end
);


InstallMethod( IsAcyclic,
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
);


InstallMethod( IsFinite,
  "for quivers",
  true,
  [ IsQuiver ], 0,
  IsAcyclic
);


InstallMethod( Size,
  "for quivers",
  true,
  [IsQuiver], 0,
  function(Q)
    local tsorted, node, arrow, count, nodeNum, child, childNum, total;

    if IsFinite(Q) then
      count := [];
      total := 1; # for 0 element
      tsorted := Q!.tsorted;

      for node in tsorted do
        count[ExtRepOfObj(node)[1]] := 0;
      od;

      for node in tsorted do
        nodeNum := ExtRepOfObj(node)[1];

        for arrow in OutgoingArrowsOfVertex(node) do
          child := TargetOfPath(arrow);
          childNum := ExtRepOfObj(child)[1];
          # Add one for the arrow
          count[childNum] := count[childNum] + count[nodeNum] + 1;
        od;

        total := total + count[nodeNum] + 1; # One for the vertex
      od;

      return total;
    else
      return infinity;
    fi;

  end
);


InstallMethod( Iterator,
  "method for quivers",
  true,
  [ IsQuiver ], 0,
  function( Q )

    local iter;
    
    iter := Objectify( NewType( IteratorsFamily, 
                       IsIterator and IsMutable and IsQuiverIteratorRep ),
                       rec( quiver := Q, position := 0 ) );
    return iter;

  end
);


InstallMethod( Enumerator,
  "method for quivers",
  true,
  [ IsQuiver ], 0,
  function( Q )

    local enum;
    
    enum := Objectify( NewType( FamilyObj( Q ), IsQuiverEnumerator ),
                       rec( quiver := Q, cache := WeakPointerObj( [] ) ) );
    return enum;

  end
);


InstallMethod( IsDoneIterator,
  "method for iterator of quivers",
  true,
  [ IsIterator and IsMutable and IsQuiverIteratorRep ], 0,
  function( iter )
    return NextPath( iter!.quiver, iter!.position ) = fail;
  end
);


InstallMethod( NextIterator,
  "method for iterator of quivers",
  true,
  [ IsIterator and IsMutable and IsQuiverIteratorRep ], 0,
  function( iter )
    local next;

    next := NextPath( iter!.quiver, iter!.position );
    if( next = fail ) then
      return fail;
    else
      iter!.position := next[1];
      return next[2];
    fi;

  end
);


InstallMethod( \[\],
  true, 
  [ IsQuiverEnumerator, IsPosInt ],
  0,
  function( enum, number )
    local path_length, paths_of_length, matrix, pos, i, j, result, orig_num;

    result := ElmWPObj( enum!.cache, number);

    if result <> fail then
        return result;
    fi;

    orig_num := number;

    if( number <= Length(VerticesOfQuiver(enum!.quiver) ) ) then
      return VerticesOfQuiver(enum!.quiver)[number];
    else
      number := number - Length(VerticesOfQuiver(enum!.quiver));
    fi;

    path_length := 1;

    while number >= 1 do
      paths_of_length := 0;
      matrix := AdjacencyMatrixOfQuiver(enum!.quiver) ^ path_length;

      if IsZero(matrix) then
        return fail;
      fi;

      for i in [ 1 .. Length(matrix) ] do
        for j in [ 1 .. Length( matrix[i] ) ] do
          paths_of_length := paths_of_length + matrix[i][j];
        od;
      od;

      if number <= paths_of_length then
        if path_length = 1 then
          pos := Length( VerticesOfQuiver(enum!.quiver) );
        else
          pos := [];
          for i in [ 1 .. path_length - 1 ] do
            pos[i] := Length( ArrowsOfQuiver(enum!.quiver) );
          od;
        fi;

        for i in [ 1 .. number ] do
          result := NextPath( enum!.quiver, pos );
          if result = fail then
            Error("This should never happen!");
          fi;
          pos := result[1];
        od;

        SetElmWPObj( enum!.cache, orig_num, result[2] );
        return result[2];

      else
        path_length := path_length + 1;
        number := number - paths_of_length;
      fi;

    od;

  end
);


InstallMethod( NextPath,
  "helper method for iterator and enumerator of quivers",
  true,
  [ IsQuiver, IsObject ], 0,
  function( Q, list )
    local arrows, word_size, current_position,
          BuildPath, IncrementPosition;
    
    current_position := ShallowCopy( list );

    if( current_position = fail ) then
      return fail;
    fi;

    arrows := ArrowsOfQuiver(Q);
    IncrementPosition := function( pos_list )
      local last_letter;

      last_letter := Length( pos_list );
      # if you have run out of arrows, increment the at last_letter-1
      if( pos_list[last_letter] = Length(arrows) ) then
        pos_list[last_letter] := 1; 
        if( last_letter > 1 ) then
          pos_list{[1 .. last_letter - 1]} :=
            IncrementPosition( pos_list{[1 .. last_letter - 1]} ){[1 .. last_letter - 1]};
        fi;
        # If you have exhausted all paths of this length, add an arrow
        if( ForAll( pos_list, p -> p = 1 ) ) then
          Add(pos_list, 1);
        fi; 
      else
        pos_list[last_letter] := pos_list[last_letter] + 1;
      fi;
      return pos_list;
    end;

    BuildPath := function( pos_list )
      local i,  path;  # Used to iterate through position list.
          
      path := arrows[ pos_list[1] ];
      for i in [2 .. Length(pos_list) ] do
        path := path * arrows[ pos_list[i] ];
      od;
      return path;
    end;
    
    # first loop through vertices of quiver
   if( not IsList( current_position ) ) then
      current_position := current_position + 1;
      if( current_position <= Length( VerticesOfQuiver(Q) ) ) then
        return [current_position, VerticesOfQuiver(Q)[current_position]];
      else
        # if you are out of vertices, start building paths
        current_position := [ 1 ];
      fi;
    else
      current_position := IncrementPosition(current_position);
    fi;
    
    while( BuildPath(current_position) = Zero(Q) ) do
      if( ForAll( current_position, p -> p = 1 ) ) then
        if( IsZero(AdjacencyMatrixOfQuiver(Q) ^ Length( current_position)) ) then
          return fail;
        fi;
      fi;
      current_position := IncrementPosition(current_position);
    od;

    return [current_position, BuildPath(current_position)];

  end
);


InstallMethod( ViewObj,
  "for quiver",
  true,
  [ IsQuiver ], 0,
  function( Q )
    Print( "<quiver with " );
    Print( OrderOfQuiver(Q) );
    Print( " vertices and " );
    Print( SizeOfQuiver(Q) );
    Print( " arrows>" );
  end
);


InstallMethod( PrintObj,
  "for quiver",
  true,
  [ IsQuiverRep ], 0,
  function( Q )
    local i, first;

    Print( "Quiver( [" );
    first := true;
    for i in VerticesOfQuiver(Q) do
      if( not first ) then
        Print(",");
      fi;
      Print("\"");
      Print( i );
      Print( "\"");
      first := false;
    od;
    
    Print( "], [" );
    
    first := true;
    for i in ArrowsOfQuiver(Q) do
      if( not first ) then
        Print(",");
      fi;
      Print( "[\"" );
      Print( SourceOfPath(i) );
      Print( "\",\"" );
      Print( TargetOfPath(i) );
      Print( "\",\"" );
      Print( i );
      Print( "\"]" );
      first := false;
    od;

    Print( "] )" );

  end
);
  
  
InstallMethod( QuiverContainingPath,
        "for a path",
        [ IsPath ],
        function( p )
    return FamilyObj( p )!.quiver;
end );
