# GAP Implementation
# $Id: opposite.gi,v 1.4 2012/02/27 12:26:34 sunnyquiver Exp $

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


InstallMethod( OppositePathAlgebra,
        "for a path algebra",
        [ IsPathAlgebra ],
        function( PA )
    local PA_op, field, quiver_op;

    field := LeftActingDomain( PA );
    quiver_op := OppositeQuiver( QuiverOfPathAlgebra( PA ) );
    PA_op := PathAlgebra( field, quiver_op );

    SetOppositePathAlgebra( PA_op, PA );

    return PA_op;
end );


InstallMethod( OppositeAlgebra,
        "for a path algebra",
        [ IsPathAlgebra ],
        OppositePathAlgebra );


InstallGlobalFunction( OppositePathAlgebraElement,
        function( elem )
    local PA, PA_op, terms, op_term;

    PA := PathAlgebraContainingElement( elem );
    PA_op := OppositePathAlgebra( PA );

    if elem = Zero(PA) then 
       return Zero(PA_op);
    else 
       op_term :=
       function( t )
          return t.coeff * One( PA_op ) * OppositePath( t.monom );
       end;

       terms := PathAlgebraElementTerms( elem );
       return Sum( terms, op_term );
    fi;
end );


InstallMethod( OppositeRelations,
        "for a list of relations",
        [ IsDenseList ],
        function( rels )
    return List( rels, OppositePathAlgebraElement );
end );
        
        
InstallMethod( OppositePathAlgebra,
        "for a quotient of a path algebra",
        [ IsQuotientOfPathAlgebra ],
        function( quot )
    local PA, PA_op, rels, rels_op, I_op, gb, gbb, quot_op;

    PA := OriginalPathAlgebra( quot );
    PA_op := OppositePathAlgebra( PA );
    rels  := RelatorsOfFpAlgebra( quot );
    rels_op := OppositeRelations( rels );
    I_op := Ideal( PA_op, rels_op );
    gb   := GBNPGroebnerBasis(rels_op,PA_op);
    gbb  := GroebnerBasis(I_op,gb);
    quot_op := PA_op / I_op;

    SetOppositePathAlgebra( quot_op, quot );

    return quot_op;
end );


InstallMethod( OppositeAlgebra,
        "for a quotient of a path algebra",
        [ IsQuotientOfPathAlgebra ],
        OppositePathAlgebra );
