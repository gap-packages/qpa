# GAP Declarations
# $Id: util.gd,v 1.1 2011/04/28 09:36:47 oysteini Exp $

# util.gd: General utilities.

# Declares functions which are not specific to QPA.  That is, they
# could have been a part of GAP itself, but they are not.

# PositionsNonzero( list ):

# Gives the indices of all nonzero elements in list.

DeclareOperation( "PositionsNonzero", [ IsList ] );


# Square( x ):

# Returns the square of the number x.

# DeclareOperation( "Square", [ IsRat ] );


# NullList( length, field )

# Returns a list with the given length where each element is the zero
# element of the field.

DeclareOperation( "NullList", [ IsPosInt, IsField ] );

