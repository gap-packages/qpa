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
    
    local Q, comps; 
    
    Q := QuiverOfPathAlgebra(A);
    if NumberOfArrows(Q) = 0 then
      Print("Finite type!\nSemisimple algebra (Q has no arrows).\n");
      return true;
    fi;
    
    comps := ConnectedComponents(Q);
    if ForAll(comps, x -> IsDynkinQuiver(x)) then
      Print("Finite type!\nQuiver is a (union of) Dynkin quiver(s).\n");
      return true;
    else
      Print("Infinite type!\nQuiver is not a (union of) Dynkin quiver(s).\n");
      return false;
    fi;
    
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
    
    local Q, PA, I, comps;
    
    PA := OriginalPathAlgebra(A);
    Q := QuiverOfPathAlgebra(PA);
    
    if NumberOfArrows(Q) = 0 then
      Print("Finite type!\nSemisimple algebra (Q has no arrows).\n");
      return true;
    fi;
    
    I := ElementsFamily(FamilyObj(A))!.ideal;
    
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

    comps := ConnectedComponents(Q);
    if ForAll(comps, x -> IsDynkinQuiver(x)) then
      Print("Finite type!\nQuiver is a (union of) Dynkin quiver(s).\n");
      return true;
    fi;
    
    Print("Can not determine the representation type.\n");
    return fail;
end
); # IsFiniteTypeAlgebra for IsQuotientOfPathAlgebra

InstallMethod( IsFiniteTypeAlgebra, 
    "for a radical square zero quotient of a path algebra",
    [ IsRadicalSquareZeroAlgebra and IsAdmissibleQuotientOfPathAlgebra ],
    0,
    function( A )
    
    local Q;
    
    Q := SeparatedQuiver(QuiverOfPathAlgebra(A));
    
    if ForAll( ConnectedComponents(Q), IsDynkinQuiver ) then
        return true;
    else
        return false;
    fi;
end
  );
