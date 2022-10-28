# GAP Implementation
# This file was generated from
# $Id: quiver.gi,v 1.5 2012/06/09 07:51:54 sunnyquiver Exp $

InstallGlobalFunction(
  Path,
  function( arg )
    local path, vertex_name, u, v, arrow_name, arrow_list;

    # single vertex name
    if Length( arg ) = 2 and IsFamily( arg[1] ) and IsString( arg[2] ) then
      vertex_name := arg[2];
      path:= Objectify( NewType( arg[1],
                        IsQuiverVertex and IsQuiverVertexRep and IsAttributeStoringRep ),
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
    elif Length( arg ) = 4 and IsFamily( arg[1]) and IsQuiverVertex( arg[2] )
                         and IsQuiverVertex( arg[3] ) and IsString( arg[4] ) then
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
  true, [IsQuiverVertex], 0,
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
  [ IsQuiverVertex, IsQuiverVertex ], 0,
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
  [ IsQuiverVertex, IsPath ], 0,
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
  [ IsPath, IsQuiverVertex ], 0,
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
  [ IsQuiverVertex, IsQuiverVertex ], 0,
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
  [ IsQuiverVertex and IsQuiverVertexRep ], 0,
  function( obj )
    return obj!.vertex_name;
  end
);


InstallMethod( PrintObj,
  "for vertex",
  true,
  [ IsQuiverVertex and IsQuiverVertexRep ], 0,
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

    local temp, vertices, arrows, vertices_by_name, arrow_spec_size,
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

    temp := 0;

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

      temp := 1;

      for i in [1..Length(arg[2])] do

        if String(arrows[i]) <> [CharInt(i+96)] then
          temp := 0;
          break;
        fi;
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

    if temp = 1 then
      Q:= Objectify( NewType( CollectionsFamily(Fam),
                              IsQuiverSA and IsQuiverRep and IsAttributeStoringRep ),
                     rec( pieces := record ) );
    else
      Q:= Objectify( NewType( CollectionsFamily(Fam),
                              IsQuiverRep and IsAttributeStoringRep ),
                     rec( pieces := record ) );
    fi;

    Fam!.quiver := Q;
    SetVerticesOfQuiver( Q, vertices );
    SetArrowsOfQuiver( Q, arrows );
    SetNumberOfVertices( Q, Length( vertices ) );
    SetNumberOfArrows( Q, Length( arrows ) );
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
    SetNumberOfVertices( new_Q, Length( new_vertices ) );
    SetNumberOfArrows( new_Q, Length( new_arrows ) );
    SetGeneratorsOfMagma( new_Q, Concatenation( new_vertices, new_arrows ) );
    SetAdjacencyMatrixOfQuiver( new_Q, AdjacencyMatrixOfQuiver(Q) );
    SetOrderingOfQuiver( new_Q, O );
    SetZero(new_Q, zero);
    SetIsWholeFamily(new_Q, true);

    return new_Q;

  end
);



InstallMethod( IsFinite,
  "for quivers",
  true,
  [ IsQuiver ], 0,
  IsAcyclicQuiver
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
                       rec( quiver := Q ) );
    return enum;

  end
);

InstallMethod( IsDoneIterator,
  "method for iterator of quivers",
  true,
  [ IsIterator and IsMutable and IsQuiverIteratorRep ], 0,
  function( iter )
    return iter!.position = [];
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
      iter!.position := next[2];
      return next[1];
    fi;

  end
);


InstallMethod( \[\],
  true,
  [ IsQuiverEnumerator, IsPosInt ],
  0,
  function( enum, number )
    local i, pos, res;

    i := 0;
    pos := 0;
    for i in [1..number] do
      res := NextPath(enum!.quiver, pos);

      if res = fail then
        return fail;
      fi;

      pos := res[2];
    od;

    return res[1];
  end
);


InstallMethod( NextPath,
  "helper method for iterator and enumerator of quivers",
  true,
  [ IsQuiver, IsObject ], 0,
  function( Q, list )
    local current_position, res, next_arrow;

    current_position := ShallowCopy( list );

    if current_position = 0 then
      current_position := ShallowCopy( VerticesOfQuiver(Q) );
    fi;
    if current_position = [] then
      return fail;
    fi;

    for next_arrow in OutgoingArrowsOfVertex( TargetOfPath(current_position[1]) ) do
      Add(current_position, current_position[1]*next_arrow);
    od;
    res := current_position[1];
    Remove(current_position, 1);
    return [res, current_position];

  end
);


InstallMethod( ViewObj,
  "for quiver",
  true,
  [ IsQuiver ], NICE_FLAGS + 1,
  function( Q )

    Print( "<quiver with " );
    Print( NumberOfVertices(Q) );
    Print( " vertices and " );
    Print( NumberOfArrows(Q) );
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


InstallMethod( VertexIndex,
        "for a vertex",
        true,
        [ IsQuiverVertex ],
        function( v )
    return v!.gen_pos;
end );


InstallMethod( ArrowIndex,
        "for an arrow",
        true,
        [ IsArrow ],
        function( a )
    return a!.gen_pos - Length( VerticesOfQuiver( QuiverContainingPath( a ) ) );
end );


InstallMethod( GeneratorIndex,
        "for a vertex",
        true,
        [ IsQuiverVertex ],
        function( v )
    return v!.gen_pos;
end );


InstallMethod( GeneratorIndex,
        "for an arrow",
        true,
        [ IsArrow ],
        function( a )
    return a!.gen_pos;
end );

InstallMethod( OppositeQuiver,
        "for a quiver",
        [ IsQuiver ],
        function( Q )
    local op_name, deop_name, op_vertex, op_arrow, Q_op;

    op_name :=
      function( name )
        return Concatenation( name, "_op" );
    end;

    deop_name :=
      function( name )
        return name{ [ 1 .. Length(name)-3 ] };
    end;

    op_vertex :=
      function( v )
        return op_name( String( v ) );
    end;

    op_arrow :=
      function( a )
        return [ op_vertex( TargetOfPath( a ) ),
                 op_vertex( SourceOfPath( a ) ),
                 op_name( String( a ) ) ];
    end;

    Q_op := Quiver( List( VerticesOfQuiver( Q ), op_vertex ),
                    List( ArrowsOfQuiver( Q ), op_arrow ) );

    SetOppositeQuiver( Q_op, Q );
    SetOppositeQuiverNameMap( Q, op_name );
    SetOppositeQuiverNameMap( Q_op, deop_name );

    return Q_op;
end );


InstallMethod( OppositePath,
        "for a path",
        [ IsPath ],
        function( p )
    local Q, Q_op, walk_names, walk_op, walk_op_names;

    Q := QuiverContainingPath( p );
    Q_op := OppositeQuiver( Q );

    walk_names := List( WalkOfPathOrVertex( p ), String );
    walk_op_names := List( Reversed( walk_names ), OppositeQuiverNameMap( Q ) );
    walk_op := List( walk_op_names,
                     function( name ) return Q_op.( name ); end );
    return Product( walk_op );
end );

InstallMethod( SeparatedQuiver,
    "for a quiver",
    true,
    [ IsQuiver ],
    0,
    function( Q )

    local oldarrows, vertices, arrows;
    #
    # Calling the vertices in the separated quiver for <label> and <label>', when
    # <label> is the name of a vertex in the old quiver. For an arrow  a : v ---> w
    # in the quiver, we call the corresponding arrow  v ---> w'  also  a  in the
    # separated quiver.
    #
    vertices := List(VerticesOfQuiver(Q), v -> String(v));
    Append(vertices, List(vertices, v -> Concatenation(v,"'")));
    oldarrows := ArrowsOfQuiver(Q);
    arrows := List(oldarrows, a -> [String(SourceOfPath(a)),Concatenation(String(TargetOfPath(a)),"'"),String(a)]);

    return Quiver(vertices,arrows);
end
  );

##########################################################################
##
#O DeclearOperations( "DynkinQuiverAn", [ n, orientation ] )
##
## Returns the Dynkin quiver A_n  with  n  vertices names 1, 2, 3,..., n,
## and arrows with names 1 -- a_1 -- 2 -- a_2 -- .... -- a_{n-1} -- n.
## The orientation is given by a list of  l's  and r's, ["r","l","l",...],
## where an  l  in coordinate  i  means that the  i-th arrow is oriented
## to the left and an  r  means oriented to the right.
##
InstallMethod( DynkinQuiverAn,
    "for a positive integer and a list",
    [ IS_INT, IsList ],

    function( n, orientation )
    local vertices, arrows, i, Q;

    if n <= 0 then
        Error("the quiver must have at least one vertex,\n ");
    fi;
    if Length(orientation) > n - 1 then
        Error("too many orientation parameters compared to the number of arrows,\n ");
    fi;
    if Length(orientation) < n - 1 then
        Error("too few orientation parameters compared to the number of arrows,\n ");
    fi;
    if not IsSubset(Set(["l","r"]), Set(orientation)) then
        Error("the orientation parameters are not in the set [\"l\",\"r\"],\n");
    fi;
    #
    # Naming the verties and the arrows in A_n.  The vertices are named
    # "1", "2", "3", and so on.  The arrows are named as follows
    #         1 -- a_1 -- 2 -- a_2 -- 3 .... n - 1 -- a_{n-1} -- n
    #
    vertices := List([1..n], i -> String(i));
    arrows := [];
    for i in [1..n-1] do
        if orientation[i] = "l" then
            Add(arrows, [vertices[i+1], vertices[i], Concatenation("a_",String(i))]);
        else
            Add(arrows, [vertices[i], vertices[i+1], Concatenation("a_",String(i))]);
        fi;
    od;
    Q := Quiver(vertices, arrows);
    SetIsDynkinQuiver(Q, true);
    return Q;
end
  );

##########################################################################
##
#O DeclearOperations( "DynkinQuiverEn", [ n, orientation ] )
##
## Returns the Dynkin quiver E_n  with  n + 1  vertices, for n = 6, 7, 8,
## names 1, 2, 3,..., n + 1, and arrows with names
##                         n
##                         |
##                       a_{n-1}
##                         |
## 1 -- a_1 -- 2 -- a_2 -- 3 -- a_3 -- 4 -- a_4 -- .... -- a_{n-2} -- n - 1.
##
## The orientation is given by a list of  n - 2  characters l's  and r's
## for the orientation of the arrows {a_1, a_2,..., a_{n-2}} and one d
## or  u  for the orientation down or up for the arrow  a_{n-1}, for instance,
##                  ["r","l","l",...,"r","d"],
## where an  l  in coordinate  i  means that the  i-th arrow is oriented
## to the left and an  r  means oriented to the right.
##
InstallMethod( DynkinQuiverEn,
    "for a positive integer and a list",
    [ IS_INT, IsList ],

    function( n, orientation )
    local vertices, arrows, i, Q;

    if not n in [6, 7, 8] then
        Error("the entered value  n  is not 6, 7 or 8,\n");
    fi;
    if Length(orientation) > n - 1 then
        Error("too many orientation parameters compared to the number of arrows,\n ");
    fi;
    if Length(orientation) < n - 1 then
        Error("too few orientation parameters compared to the number of arrows,\n ");
    fi;
    if not IsSubset(Set(["l","r"]), Set(orientation{[1..n - 2]})) then
        Error("the ",n - 1," first orientation parameters are not in the set [\"l\",\"r\"],\n");
    fi;
    if not orientation[n - 1] in Set(["d","u"]) then
        Error("the last orientation parameter is not in the set [\"d\",\"u\"],\n");
    fi;
    #
    # Naming the verties and the arrows in E_n.  The vertices are named
    # "1", "2", "3", and so on in the following way.  The arrows are named as follows
    #                             n
    #                             |
    #                            a_{n-1}
    #                             |
    #     1 -- a_1 -- 2 -- a_2 -- 3 -- a_3 -- 4 -- a_4 -- 5 .... -- a_{n-2} -- n - 1.
    #
    vertices := List([1..n], i -> String(i));
    arrows := [];
    for i in [1..n - 2] do
        if orientation[i] = "l" then
            Add(arrows, [vertices[i+1], vertices[i], Concatenation("a_",String(i))]);
        else
            Add(arrows, [vertices[i], vertices[i+1], Concatenation("a_",String(i))]);
        fi;
    od;
    if orientation[n - 1] = "d" then
        Add(arrows, [vertices[n], vertices[3], Concatenation("a_",String(n - 1))]);
    else
        Add(arrows, [vertices[3], vertices[n], Concatenation("a_",String(n - 1))]);
    fi;
    Q := Quiver(vertices, arrows);
    SetIsDynkinQuiver(Q, true);
    return Q;
end
  );

##########################################################################
##
#O DeclearOperations( "DynkinQuiverDn", [ n, orientation ] )
##
## Returns the Dynkin quiver D_n  with  n  vertices, for n >= 4, with
## names 4, 5, 6,..., n, and arrows with names
##
##     1
##       \ 
##        a_1
##         \
##          3 -- a_3 -- 4 -- a_2 -- .... n - 1 -- a_{n-1} -- n.
##         /
##        a_2
##       /
##      2
##
## The orientation is given by a list of  n - 1 characters l's  and r's for the
## orientation of the arrows {a_1, a_2,..., a_{n-1}}, for instance,
## ["r","l","l",...,"r"], where an  l  in coordinate  i  means that the
## i-th arrow is oriented to the left and an  r  means oriented to the right.
##
InstallMethod( DynkinQuiverDn,
    "for a positive integer and a list",
    [ IS_INT, IsList ],

    function( n, orientation )
    local vertices, arrows, i, Q;

    if not n >= 4 then
        Error("the entered value  n  is not greater or equal to 4,\n");
    fi;
    if Length(orientation) > n - 1 then
        Error("too many orientation parameters compared to the number of arrows,\n ");
    fi;
    if Length(orientation) < n - 1 then
        Error("too few orientation parameters compared to the number of arrows,\n ");
    fi;
    if not IsSubset(Set(["l","r"]), Set(orientation{[1..n-1]})) then
        Error("the orientation parameters are not in the set [\"l\",\"r\"],\n");
    fi;
    #
    # Naming the verties and the arrows in D_n.  The vertices are named
    # "1", "2", "3", and so on in the following way.  The arrows are named as follows
    #     1
    #       \ 
    #        a_1
    #         \
    #          3 -- a_3 -- 4 -- a_2 -- .... n - 1 -- a_{n-1} -- n.
    #         /
    #        a_2
    #       /
    #      2
    #
    vertices := List([1..n], i -> String(i));
    arrows := [];
    for i in [1..2] do
        if orientation[i] = "l" then
            Add(arrows, [vertices[3], vertices[i], Concatenation("a_",String(i))]);
        else
            Add(arrows, [vertices[i], vertices[3], Concatenation("a_",String(i))]);
        fi;
    od;
    for i in [3..n-1] do
        if orientation[i] = "l" then
            Add(arrows, [vertices[i+1], vertices[i], Concatenation("a_",String(i))]);
        else
            Add(arrows, [vertices[i], vertices[i+1], Concatenation("a_",String(i))]);
        fi;
    od;
    Q := Quiver(vertices, arrows);
    SetIsDynkinQuiver(Q, true);
    return Q;
end
  );

##########################################################################
##
#O DeclearOperation( "DynkinQuiver", [ Delta, n, orientation ] )
##
## Returns the Dynkin quiver of type  Delta, where  Delta  is A, D or E.
##
InstallMethod( DynkinQuiver,
    "for a positive integer and a list",
    [ IsString, IS_INT, IsList ],

    function( Delta, n, orientation );

    if String(Delta) = "A" then
        return DynkinQuiverAn(n, orientation);
    fi;
    if String(Delta) = "D" then
        return DynkinQuiverDn(n, orientation);
    fi;
    if String(Delta) = "E" then
        return DynkinQuiverEn(n, orientation);
    fi;
end
  );

##########################################################################
##
#A DeclearAttribute( "DoubleQuiver", [ Q ] )
##
## Returns the double of the quiver  Q  where each arrow  a  in  Q  gets a
## companion  abar  going in the opposite direction.
##
InstallMethod( DoubleQuiver,
   "for a quiver",
   [ IsQuiver ], 0,
   function( Q )

  local vertices, arrows, newvertices, newarrows, a, origin, terminus;

   vertices := VerticesOfQuiver( Q );
   arrows := ArrowsOfQuiver( Q );
   newvertices := List( vertices, String );
   newarrows := [];
   for a in arrows do
       origin := String( SourceOfPath( a ) );
       terminus := String( TargetOfPath( a ) );
       Add( newarrows, [ origin, terminus, String( a ) ] );
       Add( newarrows, [ terminus, origin, Concatenation( String( a ), "bar" ) ] );
   od;

   return Quiver( newvertices, newarrows );
end
);
