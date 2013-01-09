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
##  (i.e. it imposible from theoretical/algorithmic point of view
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
    
    local Q; 
    
    Q := QuiverOfPathAlgebra(A);
    
    if IsSpecialBiserialAlgebra(A) then
      if IsTreeQuiver(Q) then 
        # Follows from [A. Skowronski, J.Waschb¨usch. Representation-finite biserial algebras. 
        # Journal f¨ur Mathematik. Band 345, 1983.]: Special biserial algebra is of fin. rep. type
        # <=> it contains no "primitive strings".
        # TODO: full application of above thm, i.e. not only for trees.
        # But even now it covers quite large class of "simple" algebras of finite type.
        return true;
      fi;
    fi;
    
    Print("Can not determine the representation type.\n");
    return fail;  
end
); # IsFiniteTypeAlgebra for IsPathAlgebra


########################################################################
##
#P  IsFiniteTypeAlgebra( <A> ) 
##  <A> = a quotient of a path algebra
##
##	Returns true if <A> is of finite representation type.
##  Returns false if <A> is of infinite representation type.
##  Returns fail if we can not determine the representation type
##  (i.e. it imposible from theoretical/algorithmic point of view
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
    
    local Q, PA, I;
    
    PA := OriginalPathAlgebra(A);
    Q := QuiverOfPathAlgebra(PA);
    I := ElementsFamily(FamilyObj(A))!.ideal;
    
    if IsSpecialBiserialAlgebra(A) then
      if IsTreeQuiver(Q) then
        # Follows from [A. Skowronski, J.Waschb¨usch. Representation-finite biserial algebras. 
        # Journal f¨ur Mathematik. Band 345, 1983.]: Special biserial algebra is of fin. rep. type
        # <=> it contains no "primitive strings".
        # TODO: full application of above thm, i.e. not only for trees.
        # But even now it covers quite large class of "simple" algebras of finite type.
        return true;
      fi;
    fi;

    Print("Can not determine the representation type.\n");
    return fail;
end
); # IsFiniteTypeAlgebra for IsQuotientOfPathAlgebra
