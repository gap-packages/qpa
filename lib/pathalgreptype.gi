# GAP Implementation
# This contains various tools for determining a representation type of a (quotient of) path algebra
# (created A. Mroz, 08.01.2013)



#################################################################
##
#P  IsFiniteTypeAlgebra( <A> ) 
##  <A> = a path algebra
##
##	Returns true if <A> is of finite representation type.
##  Returns false if <A> is of infinite representation type.
##  Returns fail if we can not determine the representation type
##  (i.e. it impossible from theoretical/algorithmic point of view
##  or a suitable criterion has not been implemented yet).
##
##  Note for developers: this function can be a collection
##  of known facts from representation theory for determining
##  a representation type of an algebra (e.g. connected with quadratic forms etc.). 
##  Now it is only a beginning of a large (and hard) project.
##	
InstallMethod( IsFiniteTypeAlgebra,
    "for path algebras",
    true,
    [ IsPathAlgebra ], 0,
    function ( A )
    
    local Q, comps; 
    
    Q := QuiverOfPathAlgebra(A);
    if NumberOfArrows(Q) = 0 then
      Print("Finite type!\nSemisimple algebra (Q has no arrows).\n");
      return true;
    fi;
    
    comps := ConnectedComponentsOfQuiver(Q);
    if ForAll(comps, x -> IsDynkinQuiver(x)) then
      Print("Finite type!\nQuiver is a (union of) Dynkin quiver(s).\n");
      return true;
    else
      Print("Infinite type!\nQuiver is not a (union of) Dynkin quiver(s).\n");
      return false;
    fi;
    
end
); # IsFiniteTypeAlgebra for IsPathAlgebra

InstallMethod ( BongartzTest, 
"for a finite dimensional algebra",
[ IsQuiverAlgebra, IS_INT ],
function( A, bound )
        
    local   S,  P,  I,  n,  maxdim,  num,  DTrS,  TrDS,  TrDP,  DTrI,  
            i;

    if not IsAdmissibleQuotientOfPathAlgebra( A ) then
        Error( "The entered algebra is not an admissible quotient of a path algebra.\n" );
    fi;

    Print( "Entering Bongartz' test for infinite type.\n" );
    S := SimpleModules( A );
    P := IndecProjectiveModules( A );
    I := IndecInjectiveModules( A );
    n := Length( S );
    maxdim := Maximum( 2 * Dimension( A ), 30 );
    num := 0;
    DTrS := ShallowCopy( S );
    TrDS := ShallowCopy( S );
    TrDP := ShallowCopy( P );
    DTrI := ShallowCopy( I );
    while ( num < bound ) do
        for i in [ 1..n ] do
            DTrS[ i ] := DTr( DTrS[ i ] );
            if Dimension( DTrS[ i ] ) > maxdim then
                return false;
            fi;
            TrDS[ i ] := TrD( TrDS[ i ] );
            if Dimension( TrDS[ i ] ) > maxdim then
                return false;
            fi;
            TrDP[ i ] := TrD( TrDP[ i ] );
            if Dimension( TrDP[ i ] ) > maxdim then
                return false;
            fi;
            DTrI[ i ] := DTr( DTrI[ i ] );
            if Dimension( DTrI[ i ] ) > maxdim then
                return false;
            fi;
        od;
        num := num + 1;
    od;
    
    return fail;
end
  );

########################################################################
##
#P  IsFiniteTypeAlgebra( <A> ) 
##  <A> = a quotient of a path algebra
##
##	Returns true if <A> is of finite representation type.
##  Returns false if <A> is of infinite representation type.
##  Returns fail if we can not determine the representation type
##  (i.e. it impossible from theoretical/algorithmic point of view
##  or a suitable criterion has not been implemented yet).
##
##  Note for developers: this function can be a collection
##  of known facts from representation theory for determining
##  a representation type of an algebra (e.g. connected with quadratic forms etc.). 
##  Now it is only a beginning of a large (and hard) project.
##	
InstallMethod( IsFiniteTypeAlgebra,
    "for quotients of path algebras",
    true,
    [ IsQuotientOfPathAlgebra ], 0,
    function ( A )
    
    local Q, Qbar, comps;
    
    if not IsAdmissibleQuotientOfPathAlgebra(A) then 
        Error(" the entered algebra is not an admissible quotient of a path algebra,\n");
    fi; 
    Q := QuiverOfPathAlgebra(A);
    #
    # Checking if  A  is a semisimple algebra, that is, the quiver of  A
    # has no arrows.
    #
    if NumberOfArrows(Q) = 0 then
        Print("Finite type!\nSemisimple algebra (Q has no arrows).\n");
        return true;
    fi;
    #
    # Checking if  A  is a quotient of a path algebra of finite type, 
    # that is, the quiver of  A  is a union of Dynkin quivers.
    #
    comps := ConnectedComponentsOfQuiver(Q);
    if ForAll(comps, x -> IsDynkinQuiver(x)) then
        Print("Finite type!\nQuiver is a (union of) Dynkin quiver(s).\n");
        return true;
    fi;
    #
    # Checking if  A/rad^2  is of infinite type. If so, then  A  is of 
    # infinite type. The algebra  A/rad^2  is stably equivalent to  kQbar,
    # where  Qbar  is the separated quiver of  Q.  The algebra  kQbar  is 
    # of infinite type if and only if there is a connected component of  
    # Qbar which is not Dynkin. If  A  is a radical square zero algebra, 
    # then  A = A/rad^2 and  A  is of finite type if and only if  kQbar is. 
    #
    Qbar := SeparatedQuiver(Q);
    if not ForAll( ConnectedComponentsOfQuiver(Qbar), IsDynkinQuiver ) then
        Print("Have checked representation type of the algebra modulo the square of the radical.\n");
        return false;
    elif IsRadicalSquareZeroAlgebra(A) then
        return true;
    fi;
    #
    # By a result of Jans it is known that if the algebra is of finite representation type, 
    # then the algebra is distributive. Checking this next.
    #
    if not IsDistributiveAlgebra(A) then
        return false;
    fi;
    #
    # If  A  is a special biserial algebra, we do the following.
    #
    if IsSpecialBiserialAlgebra(A) then
        if IsUAcyclicQuiver(Q) then
            # Follows from [A. Skowronski, J.Waschb¨usch. Representation-finite biserial algebras. 
            # Journal f¨ur Mathematik. Band 345, 1983.]: Special biserial algebra is of fin. rep. type
            # <=> it contains no "primitive strings".
            # TODO: full application of above thm, i.e. not only for quivers with no unoriented cycles.
            # But even now it covers quite large class of "not very complicated" algebras of finite type.
            Print("Finite type!\nSpecial biserial algebra with no unoriented cycles in Q.\n");
            return true;
        fi;
    fi;
    #
    # Finally applying the test for infinite type given by Bongartz with bound 100.
    #
    if BongartzTest( A, 100 ) = false then
       return false;
    fi;
    
    Print( "Can not determine the representation type.\n" );
    return fail;
end
); # IsFiniteTypeAlgebra for IsQuotientOfPathAlgebra