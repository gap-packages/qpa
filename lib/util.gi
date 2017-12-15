# GAP Implementation
# $Id: util.gi,v 1.1 2011/04/28 09:36:47 oysteini Exp $


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

