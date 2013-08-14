# GAP Implementation
# $Id: util.gi,v 1.1 2011/04/28 09:36:47 oysteini Exp $

InstallMethod( PositionsProperty,
        "for a list and a function",
        [ IsList, IsFunction ],
        function( list, fun )
    local positions, i;

    positions := [ ];
    for i in [ 1 .. Length( list ) ] do
        if fun( list[ i ] ) then
            Append( positions, [ i ] );
        fi;
    od;

    return positions;
end );


InstallMethod( PositionsNonzero,
        "for a list",
        [ IsList ],
        function( list )
    return PositionsProperty(
                   list,
                   function( x ) return not IsZero( x ); end );
end );


#InstallMethod( Square,
#        "for a rational number",
#        [ IsRat ],
#        function( x ) return x * x; end);
        
        
InstallMethod( NullList,
        "for a positive integer and a field",
        [ IsPosInt, IsField ],
        function( length, field )
    return NullMat( 1, length, field )[ 1 ];
end );

