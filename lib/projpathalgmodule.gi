# GAP Implementation
# $Id: present.gi,v 1.8 2012/09/27 08:55:07 sunnyquiver Exp $

InstallMethod( IsNormalForm,
  "for path algebra vectors",
  true,
  [IsPathAlgebraVector], 0,
  ReturnFalse
);

#######################################################################
##
#O  PathAlgebraVector( <fam>, <components> ) 
##
##  This function constructs a PathAlgebraVector, where the tip is 
##  computed using the following order: the tip is computed for each 
##  coordinate, if the largest of these occur as a tip of several 
##  coordinates, then the coordinate with the smallest index from 1 to 
##  the length of vector is chosen.
##
##  This function is typically used when constructing elements of a module 
##  constructed by the command <C>RightProjectiveModule</C>. If <C>P</C> is 
##  constructed as say, <C>P := RightProjectiveModule(KQ, [KQ.v1, KQ.v1, KQ.v2])</C>, 
##  then <C>ExtRepOfObj(p)</C>, where <C>p</C> is an element if <C>P</C> is
##  a <C>PathAlgebraVector</C>. 
##
InstallMethod( PathAlgebraVector,
  "for families of path algebra vectors and hom. list",
  true,
  [IsPathAlgebraVectorFamily, IsHomogeneousList], 0,
  function(fam, components)
    local i, j, c, pos, tipmonomials, sortedtipmonomials, tippath;

    if not ( IsIdenticalObj( fam!.componentFam, FamilyObj(components) )
             and fam!.vectorLen = Length(components) )
    then
        Error("components are not in the proper family");
    fi;
#
#  Checking if the entered vector, components, is a vector satisfying the 
#  vertex conditions of elements in the family fam.
# 
    for i in [fam!.vectorLen,(fam!.vectorLen-1)..1] do
        c := components[i];

        if not IsLeftUniform(c) then
            Error("components must be left uniform");
        fi;

        if not IsZero(c) then
            if SourceOfPath(TipMonomial(c)) <> fam!.vertices[i] then
              Error("component ", i, " starts with an improper vertex.");
            fi;
        fi;
    od;
#
#   Finding the position of the tip of the vector.
# 
    pos := 1;
    tipmonomials := List(components,TipMonomial);
    sortedtipmonomials := ShallowCopy(tipmonomials);
    Sort(sortedtipmonomials,\<);

    if Length(tipmonomials) > 0 then
       tippath := sortedtipmonomials[Length(sortedtipmonomials)];
       pos := Minimum(Positions(tipmonomials,tippath));
    fi;

    return Objectify(fam!.defaultType, [Immutable(components), pos]);
end
);

InstallOtherMethod( PathAlgebraVector,
  "for families of path algebra vectors and hom. list (normal form)",
  true,
  [IsPathAlgebraVectorFamily and HasNormalFormFunction, 
   IsHomogeneousList], 0,
  function(fam, components)
    local i, j, c, pos, tipmonomials, sortedtipmonomials, tippath;

    if not ( IsIdenticalObj( fam!.componentFam, FamilyObj(components) )
             and fam!.vectorLen = Length(components) )
    then
        Error("components are not in the proper family");
    fi;
#
#  Checking if the entered vector, components, is a vector satisfying the 
#  vertex conditions of elements in the family fam.
# 
    for i in [fam!.vectorLen,(fam!.vectorLen-1)..1] do
        c := components[i];

        if not IsLeftUniform(c) then
            Error("components must be left uniform");
        fi;

        if not IsZero(c) then
            if SourceOfPath(TipMonomial(c)) <> fam!.vertices[i] then
              Error("component ", i, " starts with an improper vertex.");
            fi;
        fi;
    od;
#
#   Finding the position of the tip of the vector.
# 
    pos := 1;
    tipmonomials := List(components,TipMonomial);
    sortedtipmonomials := ShallowCopy(tipmonomials);
    Sort(sortedtipmonomials,\<);

    if Length(tipmonomials) > 0 then
       tippath := sortedtipmonomials[Length(sortedtipmonomials)];
       pos := Minimum(Positions(tipmonomials,tippath));
    fi;

    return NormalFormFunction(fam)(fam, components, pos);
end
);


InstallMethod( PathAlgebraVectorNC,
  "for path algebra vector family and hom. list (nonchecking)",
  true,
  [ IsPathAlgebraVectorFamily, IsHomogeneousList, IsInt, IsBool ], 0,
  function( fam, components, pos, normal )
    return Objectify( fam!.defaultType, [Immutable(components), pos] );
  end
);


InstallMethod( PathAlgebraVectorNC,
  "for path algebra vector family and hom. list (nonchecking)",
  true,
  [ IsPathAlgebraVectorFamily and HasNormalFormFunction,
    IsHomogeneousList, IsInt, IsBool ], 0,
  function( fam, components, pos, normal )
    if normal then
      return Objectify( fam!.defaultType, [Immutable(components), pos]);
    else
      return NormalFormFunction(fam)(fam, components, pos);
    fi;
  end
);

#######################################################################
##
#O  Vectorize( <M>, <vector> ) 
##
##  This function constructs an element in the module  <M>, which is a  
##  PathAlgebraModule, from the entered vector <vector>.
##  
InstallOtherMethod( Vectorize,
  "for path algebra modules and a homogeneous list",
  true,
  [IsPathAlgebraModule, IsHomogeneousList], 0,
  function( M, vector )
    local vecElmFam, modElmFam;

    vecElmFam := FamilyObj(ExtRepOfObj(Zero(M)));
    modElmFam := ElementsFamily(FamilyObj(M));

    return ObjByExtRep(modElmFam, PathAlgebraVector(vecElmFam, vector));
  end
);

#######################################################################
##
#O  \+ ( <u>, <v> ) 
##
##  This function implements the addition of two PathAlgebraVectors.
## 
InstallMethod( \+,
  "for two path algebra vectors",
  IsIdenticalObj,
  [IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep,
   IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep], 0,
  function( u, v )
    local fam, pos, c, i, normal, tipmonomials, 
          sortedtipmonomials, tippath;

    c := u![1] + v![1];
    fam := FamilyObj(u);

    pos := 1;
    tipmonomials := List(c,TipMonomial);
    sortedtipmonomials := ShallowCopy(tipmonomials);
    Sort(sortedtipmonomials,\<);

    if Length(tipmonomials) > 0 then
       tippath := sortedtipmonomials[Length(sortedtipmonomials)];
       pos := Minimum(Positions(tipmonomials,tippath));
    fi;

    normal := IsNormalForm(u) and IsNormalForm(v);
       
    return PathAlgebraVectorNC(fam, c, pos, normal);
end
);

#######################################################################
##
#O  \- ( <u>, <v> ) 
##
##  This function implements the difference of two PathAlgebraVectors.
## 
InstallMethod( \-,
  "for two path algebra vectors",
  IsIdenticalObj,
  [IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep,
   IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep], 0,
  function( u, v )
    local fam, pos, c, i, start, normal, tipmonomials, 
          sortedtipmonomials, tippath;

    c := u![1] - v![1];
    fam := FamilyObj(u);

    pos := 1;
    tipmonomials := List(c,TipMonomial);
    sortedtipmonomials := ShallowCopy(tipmonomials);
    Sort(sortedtipmonomials,\<);

    if Length(tipmonomials) > 0 then
       tippath := sortedtipmonomials[Length(sortedtipmonomials)];
       pos := Minimum(Positions(tipmonomials,tippath));
    fi;

    normal := IsNormalForm(u) and IsNormalForm(v);

    return PathAlgebraVectorNC(fam, c, pos, normal);
end
);

#######################################################################
##
#O  AdditiveInverseOp ( <x> ) 
##
##  This functionimplements the addition inverse of a PathAlgebraVector.
## 
InstallMethod( AdditiveInverseOp,
  "for path algebra vectors",
  true,
  [IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep], 0,
  function( x )
    return PathAlgebraVectorNC(FamilyObj(x), List(x![1], AdditiveInverse), 
                               x![2], IsNormalForm(x));
end
);

#######################################################################
##
#O  ZeroOp ( <x> ) 
##
##  This function returns zero element of the module where a 
##  PathAlgebraVector lives.
##
InstallMethod( ZeroOp,
  "for path algebra vectors",
  true,
  [IsPathAlgebraVector], 0,
  function( x )
    return FamilyObj(x)!.zeroVector;
  end
);

#######################################################################
##
#O  \* ( <n>, <x> ) 
##
##  This function implements scalar multiplication from the left.
##
InstallOtherMethod( \*,
  "for path algebra vectors and field element",
  true,
  [IsScalar, IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep], 0,
  function( n, x )
    local pos;

    if IsZero(n) then
      return FamilyObj(x)!.zeroVector;
    else
      pos := x![2];
    fi;
      return PathAlgebraVectorNC(FamilyObj(x), n*x![1], pos,
                                 IsNormalForm(x));
  end
);

#######################################################################
##
#O  \* ( <x>, <n> ) 
##
##  This function implements scalar multiplication from the right.
##
InstallMethod( \*,
  "for path algebra vectors and field element",
  true,
  [IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep, IsScalar], 0,
  function( x, n )
    local pos;        
    if IsZero(n) then
      return FamilyObj(x)!.zeroVector;
    else
      pos := x![2];
    fi;
    return PathAlgebraVectorNC(FamilyObj(x), n*x![1], pos,
                               IsNormalForm(x));
  end
);


InstallMethod( PrintObj,
  "for path algebra vectors",
  true,
  [IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep], 0,
  function( x )
    Print( x![1] );
  end
);

#######################################################################
##
#O  LeadingTerm ( <x> ) 
##
##  This function finds the leading term/tip of a PathAlgebraVector, where
##  the tip is computed using the following order: the tip is computed
##  for each coordinate, if the largest of these occur as a tip of 
##  several coordinates, then the coordinate with the smallest index 
##  from 1 to the length of vector is chosen. The position of the tip
##  was computed when the PathAlgebraVector was created. 
##  
InstallOtherMethod( LeadingTerm,
  "for path algebra vectors",
  true,
  [IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep], 0,
  function( x )
    local pos, lt, elem, zero, fam;

    pos := x![2];
    fam := FamilyObj(x);

    if pos <= fam!.vectorLen then
      lt := LeadingTerm(x![1][pos]);
      zero := Zero(x![1][1]); # of path algebra
      elem := ListWithIdenticalEntries(Length(x![1]), zero);
      elem[pos] := lt;

      return PathAlgebraVectorNC(FamilyObj(x), elem, pos,
                                 IsNormalForm(x));
    else
        Print("YEEEOUCH!\n");
        return x;
    fi;
end
);

#######################################################################
##
#O  LeadingCoefficient ( <x> ) 
##
##  The LeadingCoefficient of a PathAlgebraVector  <x>  is the leading
##  coefficient of the tip of the PathAlgebraVector  <x>.
##  
InstallOtherMethod( LeadingCoefficient,
  "for path algebra vectors",
  true,
  [IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep], 0,
  function( x )
    local pos, fam;

    pos := x![2];
    fam := FamilyObj(x);
    if pos <= fam!.vectorLen then
      return LeadingCoefficient(x![1][pos]);
    else
      return LeadingCoefficient(x![1][1]);
    fi;
end
);

#######################################################################
##
#A  TargetVertex( <v> ) 
##
##  Given a PathAlgebraVector  <v>, if  <v>  is right uniform, this 
##  function finds the vertex  w  such that  <v>*w = <v>  whenever <v>  
##  is non-zero, and returns the zero path otherwise. If  <v>  is not 
##  right uniform it returns fail. 
##  
InstallOtherMethod( TargetVertex,
  "for path algebra vectors",
  true,
  [ IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep ], 0,
  function( v )
    if IsRightUniform(v) then
      if IsZero(v) then
        return FamilyObj(v)!.zeroPath;
      else
        return TargetOfPath(LeadingMonomial(v![1][v![2]]));
      fi;
    fi;
    return fail;
  end
);

#######################################################################
##  Don't understand what this function really is supposed to do.
##  Cannot see that it is used.
##  
#A  SourceVertex( <v> ) 
##
##  Given a PathAlgebraVector  <v>, if  <v>  is right uniform, this 
##  function finds the vertex  w  such that  w*LeadingTerm(<v> = 
##  LeadingTerm(<v>)  whenever <v> is non-zero, and returns the zero 
##  path otherwise. If  <v>  is not right uniform it returns fail.
## 
InstallOtherMethod( SourceVertex,
  "for path algebra vectors",
  true,
  [ IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep ], 0,
  function( v )
    if IsRightUniform(v) then
      if IsZero(v) then
        return FamilyObj(v)!.zeroPath;
      else
        return SourceOfPath(LeadingMonomial(v![1][v![2]]));
      fi;
    fi;

    return fail;
end
);

#######################################################################
##
#A  LeadingComponent( <v> ) 
##
##  Given a PathAlgebraVector  <v>, this function finds the vertex 
##  coordinate  pos  where the tip of the vector occurs, and returns
##  v[pos], whenever <v> is non-zero. It returns the zero otherwise. 
## 
InstallMethod( LeadingComponent,
  "for path algebra vectors",
  true,
  [ IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep ], 0,
  function( v )
    if IsZero(v) then
      return v![1][ 1 ];
    fi;
    return v![1][ v![2] ];
  end
);

#######################################################################
##
#A  LeadingPosition( <v> ) 
##
##  Given a PathAlgebraVector  <v>, the LeadingPosition is in which 
##  coordinate the tip of the vector occurs.
## 
InstallMethod( LeadingPosition,
  "for path algebra vectors",
  true,
  [ IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep ], 0,
  function( v )
    return v![2];
  end
);

#######################################################################
##
#O  \< ( <x>, <y> ) 
##
##  This operation defined the ordering on the elements of a 
##  PathAlgebraModule: (1) If the tipmonomial of the leading term of  
##  <x>  is strictly less than that of  <y>, then  <x> < <y>, (2) if 
##  the tipmonomial of the leading term of  <x>  and of  <y>  are 
##  identical, then if the position in which coordinate they occor  
##  is strictly larger for  <x>  than  for  <y>, then  <x> < <y>, (3)
##  otherwise  <x>  is not less than  <y>.  
## 
InstallMethod( \<,
  "for path algebra vectors",
  IsIdenticalObj,
  [IsPathAlgebraVector 
     and IsPathAlgebraVectorDefaultRep
     and IsNormalForm,
   IsPathAlgebraVector 
     and IsPathAlgebraVectorDefaultRep
     and IsNormalForm], 0,
  function(x, y)
    local fam;

    fam := FamilyObj(x);

    if TipMonomial(x![1][x![2]]) < TipMonomial(y![1][y![2]]) then
	return true;
    elif TipMonomial(x![1][x![2]]) = TipMonomial(y![1][y![2]]) then
	return x![2] > y![2];
    else
	return false;
    fi; 
end
);

#######################################################################
##
#O  \=( <x>, <y> ) 
##
##  Given two PathAlgebraVectors  <x>  and  <y>, then they are considered
##  equal if all the coordinates are the same.
##  
InstallMethod( \=,
  "for path algebra vectors",
  IsIdenticalObj,
  [IsPathAlgebraVector
     and IsPathAlgebraVectorDefaultRep 
     and IsNormalForm,
   IsPathAlgebraVector 
     and IsPathAlgebraVectorDefaultRep 
     and IsNormalForm], 0,
  function(x, y)
    return x![1] = y![1];
  end
);
  
#######################################################################
##
#O  IsLeftDivisible( <x>, <y> ) 
##
##  Given two PathAlgebraVectors  <x>  and  <y>,  then  <y>  is said to 
##  left divide  <x>, if the tip of  <x>  and the tip of  <y>  occur in
##  the same coordinate, and the tipmonomial of the tip of  <y>  left-
##  divides the tipmonomial of the tip of  <x>. 
##
InstallMethod( IsLeftDivisible,
  "for path algebra vectors",
  IsIdenticalObj,
  [IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep,
   IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep], 0,
  function( x, y )
    local xTipPos, xLeadingTerm, xCoeff, xMon, xWalk,
          yTipPos, yLeadingTerm, yCoeff, yMon, yWalk;

    xTipPos := x![2];
    xLeadingTerm := LeadingTerm(x)![1][xTipPos];
    xMon := TipMonomial(xLeadingTerm);
    xWalk := WalkOfPath(xMon);
    yTipPos := y![2];
    yLeadingTerm := LeadingTerm(y)![1][yTipPos];
    yMon := TipMonomial(yLeadingTerm);
    yWalk := WalkOfPath(yMon);

    return xTipPos = yTipPos and not IsEmpty(xWalk)
                             and PositionSublist(xWalk, yWalk) = 1;
  end
);

#######################################################################
##
#O  IsLeftUniform( <v> ) 
##
##  Given a PathAlgebraVectors  <v>, this function tests if there is a
##  vertex  <e>  such that  <e>*<v> = <v>.
##
InstallOtherMethod( IsLeftUniform,
  "for path algebra vectors",
  true,
  [ IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep ], 0,
  function( v )
    local s, verts;

    if IsZero(v) then
      return true;
    fi;

    if ForAll(v![1], IsLeftUniform) then
      verts := Filtered( v![1], x -> not IsZero(x) );
      verts := List( verts, function(x)
                     if IsPackedElementDefaultRep(x) then
                        return CoefficientsAndMagmaElements(x![1])[1];
                     else
                        return CoefficientsAndMagmaElements(x)[1];
                     fi;
               end );
      verts := List( verts, x -> SourceOfPath(x) );
      verts := AsSet(verts);
      return Length(verts) = 1;
    fi;

    return false;        
end
);

#######################################################################
##
#O  IsRightUniform( <v> ) 
##
##  Given a PathAlgebraVectors  <v>, this function tests if there is a
##  vertex  <e>  such that  <v>*<e> = <v>.
##
InstallOtherMethod( IsRightUniform,
  "for path algebra vectors",
  true,
  [ IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep ], 0,
  function( v )
    local s, verts;

    if IsZero(v) then
      return true;
    fi;

    if ForAll(v![1], IsRightUniform) then
      verts := Filtered( v![1], x -> not IsZero(x) );
      verts := List( verts, function(x)
                if IsPackedElementDefaultRep(x) then
                    return CoefficientsAndMagmaElements(x![1])[1];
                else
                    return CoefficientsAndMagmaElements(x)[1];
                fi;
            end );
      verts := List( verts, x -> TargetOfPath(x) );
      verts := AsSet(verts);
      return Length(verts) = 1;
    fi;

    return false;        

  end
);

#######################################################################
##
#O  IsUniform( <v> ) 
##
##  Given a PathAlgebraVectors  <v>, this function tests if there is a
##  vertex  <e>  such that  <e>*<v> = <v> = <v>*<e>, that is, if  <v> 
##  is both left and right uniform.
##
InstallOtherMethod( IsUniform,
  "for path algebra vectors",
  true,
  [ IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep ], 0,
  function( v )
    local s, verts;

    if IsZero(v) then
      return true;
    fi;

    return IsLeftUniform(v) and IsRightUniform(v);
end
);


InstallMethod( NewBasis,
  "for a space of path algebra vectors and a homogeneous list",
  IsIdenticalObj,
  [ IsFreeLeftModule and IsPathAlgebraVectorCollection, IsList ], 0,
  function( V, vectors )
    local B;  # The basis

    B := Objectify( NewType( FamilyObj( V ),
                             IsBasis and
                             IsBasisOfPathAlgebraVectorSpace and
                             IsAttributeStoringRep ),
                    rec() );
    SetUnderlyingLeftModule( B, V );
    SetBasisVectors( B, vectors );

    return B;
end
  );

InstallMethod( NewBasis,
  "for a space of path algebra vectors and a homogeneous list",
  IsIdenticalObj,
  [ IsFreeLeftModule and IsPathAlgebraModule, IsList ], 0,
  function( V, vectors )
    local B;  # The basis

    B := Objectify( NewType( FamilyObj( V ),
                             IsBasis and
                             IsBasisOfPathAlgebraVectorSpace and
                             IsAttributeStoringRep ),
                    rec() );
    SetUnderlyingLeftModule( B, V );
    SetBasisVectors( B, vectors );

    return B;
end
);


HOPF_TriangularizePathAlgebraVectorList := function( vectors )
  local basechange, heads, k, head, cf, pos, i, one, testHead,
      b, b1, sort, hpos;

  vectors := Filtered( vectors, x -> not IsZero(x) );

  if Length(vectors) = 0 then
    return rec( echelonbas := vectors, heads := [], basechange := [] );
  fi;

  one := One(LeadingCoefficient(vectors[1]));
  basechange := List( [1..Length(vectors)], x -> [ [x, one] ] );
  SortParallel( vectors, basechange );
  heads := [];
  k := Length(vectors);

  while k > 0 do

    sort := false;
    cf := LeadingCoefficient(vectors[k]);
    vectors[k] := vectors[k]/cf;

    for i in [1..Length(basechange[k])] do
      basechange[k][i][2] := basechange[k][i][2]/cf;
    od;

    head := LeadingTerm(vectors[k]);
    hpos := LeadingPosition(head);
    Add(heads, head);

    for i in Reversed([1..k-1]) do
        pos := LeadingPosition(vectors[i]);
        testHead := LeadingTerm(vectors[i]);

        if pos = hpos and head![1][pos] = testHead![1][pos] then
          sort := true;
          cf := LeadingCoefficient(vectors[i]);
          vectors[i] := vectors[i] - cf*vectors[k];

          for b in basechange[k] do
            b1 := [ b[1], -cf*b[2] ];
            pos := PositionSorted(basechange[i], b1,
                           function( x, y ) return x[1] < y[1];
                   end );

            if Length(basechange[i]) < pos or basechange[i][pos][1] <> b1[1] then
              InsertElmList( basechange[i], pos, b1 );
            else
              basechange[i][pos][2] := basechange[i][pos][2] + b1[2];
            fi;

          od;
        fi;
      od;

      if sort then
        # Remove any zeroes and sort the list again
        for i in [1..k-1] do
          if IsZero(vectors[i]) then
            Unbind(vectors[i]);
            Unbind(basechange[i]);
            k := k-1;
          fi;
        od;
        vectors := Filtered(vectors, x -> IsBound(x));
        basechange := Filtered(basechange, x -> IsBound(x));
        SortParallel(vectors, basechange);
      fi; 

      k := k-1;
   od;

   return rec( echelonbas := vectors, 
               heads := Reversed(heads), 
               basechange := basechange );
end;


InstallMethod( BasisNC,
  "for a space spanned by path alg. vec., and a list of path alg. vec.",
  IsIdenticalObj,
  [ IsFreeLeftModule and IsPathAlgebraVectorCollection,
    IsPathAlgebraVectorCollection and IsList ], 0,
  function( V, vectors )
    local B, info;

    vectors := Filtered(vectors, x -> not IsZero(x));
    info := HOPF_TriangularizePathAlgebraVectorList( vectors );

    if Length(info.echelonbas) <> Length(vectors) then 
      return fail; 
    fi;

    B := NewBasis( V, vectors );
    B!.echelonBasis := info.echelonbas;
    B!.heads := info.heads;
    B!.baseChange := info.basechange;

    return B;

  end
);


InstallMethod( Basis,
  "for a space spanned by path alg. vec., and a list of path alg. vec.",
  IsIdenticalObj,
  [ IsFreeLeftModule and IsPathAlgebraVectorCollection,
    IsPathAlgebraVectorCollection and IsList ], NICE_FLAGS+1,
  function( V, vectors )
    if not ForAll( vectors, x -> x in V ) then
      return fail;
    fi;

    return BasisNC(V, vectors);

  end
);


InstallMethod( BasisOfDomain,
  "for a space spanned by path alg. vec.",
  true,
  [IsFreeLeftModule and IsPathAlgebraVectorCollection], NICE_FLAGS+1,
  function( V )
    local B, info, vectors;

    vectors := ShallowCopy(GeneratorsOfLeftModule(V));
    info := HOPF_TriangularizePathAlgebraVectorList( vectors );
    B := NewBasis( V, info.echelonbas );
    B!.echelonBasis := info.echelonbas;
    B!.heads := info.heads;
    B!.baseChange := List( [1..Length(info.echelonbas)],
                           x -> [[ x, One(LeftActingDomain(V)) ]] );
    return B;

  end
);


InstallMethod( Basis,
  "for a space spanned by path alg. vec.",
  true,
  [IsFreeLeftModule and IsPathAlgebraVectorCollection], NICE_FLAGS+1,
  function( V )
    local B, info, vectors;
    vectors := ShallowCopy(GeneratorsOfLeftModule(V));
    info := HOPF_TriangularizePathAlgebraVectorList( vectors );
    B := NewBasis( V, info.echelonbas );
    B!.echelonBasis := info.echelonbas;
    B!.heads := info.heads;
    B!.baseChange := List( [1..Length(info.echelonbas)],
                           x -> [[ x, One(LeftActingDomain(V)) ]] );
    return B;

  end
);

#O
InstallMethod( Coefficients,
  "for basis of path alg. vec, and path alg. vec.",
  true,
  [IsBasisOfPathAlgebraVectorSpace, 
   IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep], 0,
  function( B, v )
    local w, cf, i, b, c, zero, pos, lm, hlm, hpos;

    w := v;
    zero := Zero(LeadingCoefficient(v));
    cf := List( BasisVectors(B), x -> zero );
    i := Length(B!.heads);

    while i > 0 and not IsZero(w) do
      lm := LeadingMonomial(LeadingComponent(w));
      pos := LeadingPosition(w);

      # find a matching basis vector
      repeat
        hlm := LeadingMonomial(LeadingComponent(B!.heads[i]));
        hpos := LeadingPosition(B!.heads[i]);
        i := i - 1;
      until (hpos = pos and hlm <= lm) or pos < hpos or i = 0;

      # Either no match was found, so we're done or
      # match was found, so we update the coefficients.
      if lm <> hlm then
        return fail;
      else
        # the leading monomials match
        c := LeadingCoefficient(w);
        w := w - c*B!.echelonBasis[i + 1];
        for b in B!.baseChange[i + 1] do
          cf[b[1]] := cf[b[1]] + b[2]*c;
        od;
      fi;

    od;

    if IsZero(w) then
      return cf;
    fi;

    return fail;
  end
);


InstallMethod( \in,
  "for path algebra vector and space generated by path alg. vec.",
  IsElmsColls,
  [IsPathAlgebraVector, 
   IsFreeLeftModule and IsPathAlgebraVectorCollection], NICE_FLAGS+1,
  function( v, V )
    local B, cf;

    B := BasisOfDomain(V);
    cf := Coefficients(B, v);

    return cf <> fail;

  end
);



#######################################################################
##
#O  RightProjectiveModule( <A>, <els> )
##
##  This function takes as input a  PathAlgebra  <A>  and a list of 
##  uniform elements  <els>  in  <A>, and constructs the projective 
##  <A>-module  <Math>\oplus_{e \in els} eA</Math>. The filters 
##  IsPathAlgebraModule and IsVertexProjectiveModule are set true for 
##  the output. Given an element  x  in module constructed by this 
##  command, then  ExtRepOfObj(x)  is a PathAlgebraVector.  
##
InstallMethod( RightProjectiveModule,
  "for path algebra A and list of uniform elements of A",
  IsIdenticalObj,
  [ IsPathAlgebra, IsHomogeneousList ], 0,
  function( A, els )
    local i, fam, zero, gen, gens, M;

    if not ( ForAll(els, IsUniform) and ForAll(els, x-> \in(x, A)) ) then
      TryNextMethod();
      Error("<els> should be a list of right uniform elements in ",A,"\n");
    fi;

    fam := NewFamily( "IsPathAlgebraVectorFamily", IsPathAlgebraVector );
    fam!.defaultType := NewType( fam, IsPathAlgebraVectorDefaultRep
                                      and IsNormalForm );
    fam!.vectorLen := Length(els);
    fam!.elements := Immutable(els);
    fam!.componentFam := FamilyObj(A);
    fam!.zeroPath := Zero(QuiverOfPathAlgebra(A));
    zero := Zero(A);
    fam!.zeroVector := PathAlgebraVectorNC(fam,
        ListWithIdenticalEntries(fam!.vectorLen, zero),
        fam!.vectorLen + 1, true);

    gens := [];
    for i in [1..fam!.vectorLen] do
      gen := ListWithIdenticalEntries(fam!.vectorLen, zero);
      gen[i] := els[i];
      Add(gens, PathAlgebraVectorNC( fam, gen, i, true));
    od;

    M := RightAlgebraModuleByGenerators( A, \^, gens );
    SetIsPathAlgebraModule( M, true );
    SetIsVertexProjectiveModule( M, true );
    SetIsWholeFamily( M, true );
    fam!.wholeModule := M;

    return M;

  end
);

#######################################################################
##
#O  RightProjectiveModule( <A>, <els> )
##
##  This function takes as input a  PathAlgebra  <A>  and a list of 
##  trivial paths in  <A>, and constructs the projective <A>-module  
##  <Math>\oplus_{e \in els} eA</Math>. The filters IsPathAlgebraModule
##  and IsVertexProjectiveModule are set true for the output. Given an 
##  element  x  in module constructed by this command, then  
##  ExtRepOfObj(x)  is a PathAlgebraVector.  
##
InstallMethod( RightProjectiveModule,
    "for path algebra and list of vertices",
    IsIdenticalObj,
    [ IsPathAlgebra, IsHomogeneousList ], 0,
    function( A, verts )
    local i, fam, zero, gen, gens, M;

    if not ForAll(verts, function(x)
        return x = Tip(x) and IsQuiverVertex(TipMonomial(x))
               and One(TipCoefficient(x)) = TipCoefficient(x);
         end ) then
         TryNextMethod();
         Error("<verts> should be a list of embedded vertices");
     fi;

     fam := NewFamily( "IsPathAlgebraVectorFamily", IsPathAlgebraVector );
     fam!.defaultType := NewType( fam, IsPathAlgebraVectorDefaultRep
                                 and IsNormalForm );
     fam!.vectorLen := Length(verts);
     fam!.vertices := Immutable(List(verts, TipMonomial));
     fam!.componentFam := FamilyObj(A);
     fam!.zeroPath := Zero(QuiverOfPathAlgebra(A));
     zero := Zero(A);
     fam!.zeroVector := PathAlgebraVectorNC(fam,
                                ListWithIdenticalEntries(fam!.vectorLen, zero),
                                fam!.vectorLen + 1, true);
     
     gens := [];
     for i in [1..fam!.vectorLen] do
         gen := ListWithIdenticalEntries(fam!.vectorLen, zero);
         gen[i] := verts[i];
         Add(gens, PathAlgebraVectorNC( fam, gen, i, true));
     od;
     
     M := RightAlgebraModuleByGenerators( A, \^, gens );
     SetIsPathAlgebraModule( M, true );
     SetIsVertexProjectiveModule( M, true );
     SetIsWholeFamily( M, true );
     fam!.wholeModule := M;
     
     return M;
end
);

#
#  Not sure what this version is doing.
#
InstallMethod( RightProjectiveModule,
    "for f.p. path algebras and vertices",
    IsIdenticalObj,
    [ IsQuotientOfPathAlgebra, IsHomogeneousList ], 0,
    function( A, verts )
    local i, fam, zero, gen, gens, M;

    if not ForAll(verts, function(x)
        return x = Tip(x)
               and IsQuiverVertex(TipMonomial(x))
               and One(TipCoefficient(x)) = TipCoefficient(x);
        end ) then
        TryNextMethod();
        Error("<verts> should be a list of embedded vertices");
    fi;

    fam := NewFamily( "IsPathAlgebraVectorFamily", IsPathAlgebraVector );
    
    if HasNormalFormFunction(ElementsFamily(FamilyObj(A))) then
        fam!.defaultType := NewType( fam, IsPathAlgebraVectorDefaultRep
                                        and IsNormalForm );
    else
        fam!.defaultType := NewType( fam, IsPathAlgebraVectorDefaultRep );
    fi;

    fam!.vectorLen := Length(verts);
    fam!.vertices := Immutable(List(verts, TipMonomial));
    fam!.componentFam := FamilyObj(A);
    fam!.zeroPath := Zero(QuiverOfPathAlgebra(A));
    zero := Zero(A);

    fam!.zeroVector := PathAlgebraVectorNC(fam,
                               ListWithIdenticalEntries(fam!.vectorLen, zero),
                               fam!.vectorLen + 1, true);
    
    gens := [];
    for i in [1..fam!.vectorLen] do
        gen := ListWithIdenticalEntries(fam!.vectorLen, zero);
        gen[i] := verts[i];
        Add(gens, PathAlgebraVectorNC( fam, gen, i, true));
    od;
    
    M := RightAlgebraModuleByGenerators( A, \^, gens );
    SetIsPathAlgebraModule( M, true );
    SetIsWholeFamily( M, true );
    SetIsVertexProjectiveModule( M, true );
    fam!.wholeModule := M;
    
    return M;
end
);


InstallMethod( SubAlgebraModule,
  "for path algebra module and a list of submodule generators",
  IsIdenticalObj,
  [ IsFreeLeftModule and IsPathAlgebraModule,
    IsAlgebraModuleElementCollection and IsList], 0,
  function( V, gens )
    local sub;

    sub := Objectify( NewType( FamilyObj( V ),
                      IsLeftModule and IsAttributeStoringRep),
                      rec() );
    SetParent( sub, V );
    SetIsAlgebraModule( sub, true );
    SetIsPathAlgebraModule( sub, true );
    SetLeftActingDomain( sub, LeftActingDomain(V) );
    SetGeneratorsOfAlgebraModule( sub, gens );

    if HasIsFiniteDimensional( V ) and IsFiniteDimensional( V ) then
      SetIsFiniteDimensional( sub, true );
    fi;

    if IsLeftAlgebraModuleElementCollection( V ) then
#Print("left\n");
      SetLeftActingAlgebra( sub, LeftActingAlgebra( V ) );
    fi;

    if IsRightAlgebraModuleElementCollection( V ) then
#Print("right\n");
      SetRightActingAlgebra( sub, RightActingAlgebra( V ) );
    fi;

    return sub;

  end
);


InstallOtherMethod( SubAlgebraModule,
  "for path algebra module and a list of submodule generators",
  function(F1,F2,F3) return IsIdenticalObj( F1, F2 ); end,
  [ IsFreeLeftModule and IsPathAlgebraModule,
    IsAlgebraModuleElementCollection and IsList,
    IsString ], 0,
  function( V, gens, str )
    local sub;

    if str <> "basis" then
      Error("Usage: SubAlgebraModule( <V>, <gens>, <str> ) where the last argument is string \"basis\"" );
    fi;

    sub := SubAlgebraModule( V, gens );
    SetGeneratorsOfLeftModule( sub, gens );
#                          SetBasisOfDomain( sub, NewBasis( sub, gens ) );
    SetBasis( sub, NewBasis( sub, gens ) );
    SetDimension( sub, Length( gens ) );

    return sub;

  end
);

  
#P
InstallMethod( IsVertexProjectiveModule,
    "for a path algebra module",
    true,
    [ IsPathAlgebraModule ], 0,
    function( M )
    local vecFam, A, verts, gen, zero, i, vs;

    A := ActingAlgebra( M );
    zero := Zero(A);
    vecFam := ElementsFamily(FamilyObj(M))!.underlyingModuleEltsFam;
    verts := vecFam!.vertices;

    for i in [1..vecFam!.vectorLen] do
        gen := ListWithIdenticalEntries(vecFam!.vectorLen, zero);
        gen[i] := verts[i]*One(A);
        gen := PathAlgebraVectorNC(vecFam, gen, i, true);
        vs := UnderlyingLeftModule(Basis(M)!.delegateBasis);

        if not gen in vs then
            return false;
        fi;
    od;

    return true;
end
);


InstallOtherMethod( \^,
  "for path algebra vectors and path algebra",
  function(x, a)
    return IsIdenticalObj( ElementsFamily(x!.componentFam), a );
  end,
  [IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep, IsRingElement], 0,
  function( x, a )
    local i, fam, components, pos;

    fam := FamilyObj(x);
    components := x![1]*a;
    pos := fam!.vectorLen+1;

#    for i in [x![2]..fam!.vectorLen] do
#      if not IsZero(components[i]) then
#        pos := i;
#        break;
#      fi;
#    od;

#    return PathAlgebraVectorNC(FamilyObj(x), components, pos, false);
    return PathAlgebraVector(FamilyObj(x), components);
end
);


InstallMethod( CanonicalBasis,
  "for vertex projective path algebra modules",
  true,
  [ IsFreeLeftModule 
    and IsPathAlgebraModule 
    and IsVertexProjectiveModule 
    and IsRightAlgebraModuleElementCollection ], NICE_FLAGS+1,
  function(M)
    local A, abv, vecFam, mbv, MB, elem, bv, zero, i;

    A := ActingAlgebra(M);
    zero := Zero(A);
    abv := BasisVectors(CanonicalBasis(A));

    vecFam := ElementsFamily(FamilyObj(M))!.underlyingModuleEltsFam;
    mbv := []; 

    for i in [1..vecFam!.vectorLen] do
      for elem in abv do
        if not IsZero(vecFam!.vertices[i]*elem) then
          bv := ListWithIdenticalEntries(vecFam!.vectorLen, zero);
          bv[i] := elem;
          Add(mbv, Vectorize(M, bv));
        fi;
      od;
    od;

    Sort(mbv);
    MB := BasisNC(M, mbv);
    SetIsCanonicalBasis( MB, true );

    return MB;
  end
);


InstallMethod( Basis,
    "for vertex projective path algebra modules",
    true,
    [ IsFreeLeftModule 
      and IsPathAlgebraModule 
      and IsVertexProjectiveModule
      and IsRightAlgebraModuleElementCollection ], NICE_FLAGS+1,
      CanonicalBasis
);


InstallMethod( Basis,
  "for a path algebra module",
  true,
  [ IsFreeLeftModule 
    and IsPathAlgebraModule
    and IsRightAlgebraModuleElementCollection ], NICE_FLAGS+1,
  function( M )
    local fam, gens, gb, W, vecs, newbVecs;

    if not IsPathAlgebra(ActingAlgebra(M)) then
      TryNextMethod();
    fi;

    fam := ElementsFamily(FamilyObj(M))!.underlyingModuleEltsFam;
    gens := GeneratorsOfAlgebraModule( M );
    gb := RightGroebnerBasisOfModule( M );
    W := fam!.wholeModule;
    vecs := BasisVectors( CanonicalBasis(W) );

    newbVecs := List(vecs, x -> x - CompletelyReduce( gb, x ));
    newbVecs := AsSet(Filtered(newbVecs, x -> not IsZero(x)));

    SetDimension( M, Length(newbVecs) );

    return NewBasis( M, newbVecs );

end
);

#A
#######################################################################
##
#A  UniformGeneratorsOfModule( <M> )
##
##  This function takes as input a  PathAlgebraModule  <M>  and 
##  constructs a set of right uniform generators of the module  <M>. 
##  If  <M>  is the zero module, then it returns an empty list.
##
InstallMethod( UniformGeneratorsOfModule,
  "for path algebra modules",
  true,
  [ IsPathAlgebraModule and 
    IsRightAlgebraModuleElementCollection ], 0,
  function( M )
    local A, Q, uniformGens, v, uniformGen, gen;

    A := RightActingAlgebra(M);
    Q := QuiverOfPathAlgebra(A);
    uniformGens := [];
    for v in VerticesOfQuiver(Q) do
      for gen in GeneratorsOfAlgebraModule(M) do
        uniformGen := gen^(v*One(A));
        if not IsZero(uniformGen) then
          AddSet(uniformGens, uniformGen);
        fi;
      od;
    od;

    return uniformGens;
end
);


# reduce x by y
ReduceRightModuleElement := function(x, y)
  local xTipPos, xLeadingTerm, xCoeff, xMon, xWalk,
        yTipPos, yLeadingTerm, yCoeff, yMon, yWalk,
        first, last, path;

  xTipPos := x![2];
  xLeadingTerm := LeadingTerm(x)![1][xTipPos];
  xCoeff := TipCoefficient(xLeadingTerm);
  xMon := TipMonomial(xLeadingTerm);
  xWalk := WalkOfPath(xMon);
  yTipPos := y![2];
  yLeadingTerm := LeadingTerm(y)![1][yTipPos];
  yCoeff := TipCoefficient(yLeadingTerm);
  yMon := TipMonomial(yLeadingTerm);
  yWalk := WalkOfPath(yMon);

  path := One(xLeadingTerm);
  first := Length(yWalk) + 1;
  last := Length(xWalk);

  if first <= last then
    path := path*Product(xWalk{[first..last]});
  fi;

  path := xCoeff/yCoeff * path;
  x := x - y^path;

  return x;
end;

########################################################
#
#  Doesn't seem to be used anywhere.
#
NormalizeRightModuleElement := function(G, x)
  local g, newX, reduced;

  newX := Zero(x);

  while not IsZero(x) do
    reduced := false;

    for g in G do
      if IsLeftDivisible(x, g) then
        x := ReduceRightModuleElement( x, g );
        reduced := true;
        break;
      fi;
    od;

    if not reduced then
      newX := newX + LeadingTerm(x);
      x := x - LeadingTerm(x);
    fi;

  od;

  return newX;
end;

#######################################################################
##
#A  RightGroebnerBasisOfModule( <M> ) 
##
##  This function takes as input a  PathAlgebraModule  and constructs
##  a right Groebner basis. 
##  
InstallMethod( RightGroebnerBasisOfModule,
  "for a path algebra module",
  true,
  [ IsPathAlgebraModule and 
    IsRightAlgebraModuleElementCollection ], 0,
  function( M )
    local A, Q, uniformGens,
          H, Hlen, redset, i, j, r, s, reducible, tipPos, toReduce,
          staticDictionaries, gbasisElems, gb, fam, flag;

    A := RightActingAlgebra(M);

    if not IsPathAlgebra(A) then
      Error("The acting domain must be a path algebra");
    fi;

    Q := QuiverOfPathAlgebra(A);
    uniformGens := UniformGeneratorsOfModule( M );
    H := List(uniformGens, ExtRepOfObj);

# Print("H: ",H,"\n");

    # Tip reduce H to create a Right Groebner Basis:
    reducible := true;
    while reducible do

      # remove zeros:
      H := Filtered(H, x -> not IsZero(x));

      reducible := false;

      Hlen := Length(H);

      for i in [1..Hlen] do
        redset := Difference([1..Hlen],[i]);
        for j in redset do
          # if H[j] divides H[i]:
          if IsLeftDivisible(H[i], H[j]) then
            H[i] := ReduceRightModuleElement( H[i], H[j] );
            reducible := true;
            break;
          fi;        
        od;
        if reducible then
          break;
        fi;
      od;

    od;

    fam := ElementsFamily(FamilyObj(H));
    staticDictionaries := [];
    gbasisElems := List([1..fam!.vectorLen], x -> []);

    for i in H do
      Add(gbasisElems[i![2]], i);
    od;

    for i in [1..fam!.vectorLen] do
      if not IsEmpty(gbasisElems[i]) then
        staticDictionaries[i] := QuiverStaticDictionary(Q, List(gbasisElems[i], 
                         x -> LeadingMonomial(x![1][x![2]])) );
      fi;
    od;

    gb := Objectify( NewType(FamilyObj(M),
                     IsPathAlgebraModuleGroebnerBasis and
                     IsPathAlgebraModuleGroebnerBasisDefaultRep),
                     rec( gbasisElems := Immutable(gbasisElems),
                          staticDictionaries := Immutable(staticDictionaries) )
                   );

    gbasisElems := List(Flat(gbasisElems), 
                        x -> ObjByExtRep(ElementsFamily(FamilyObj(M)), x));
    SetBasisVectors(gb, gbasisElems);
    SetUnderlyingModule(gb, M);
    SetIsRightPathAlgebraModuleGroebnerBasis(gb, true);

    return gb;
end
);


InstallMethod( ViewObj,
  "for a Groebner basis of a path algebra module",
  true,
  [ IsPathAlgebraModuleGroebnerBasis and HasBasisVectors ],
  0,
  function( B )
    Print( "<Groebner basis of " );
    View( UnderlyingModule( B ) );
    Print(", ");
    View( BasisVectors( B ) );
    Print( " >" );
  end
);


InstallMethod( ViewObj,
  "for a basis of a right module",
  true,
  [ IsPathAlgebraModuleGroebnerBasis ],
  0,
  function( B )
    Print( "<Groebner basis of " );
    View( UnderlyingModule( B ) );
    Print(", ...>");
  end
);

#######################################################################
##
#O  CompletelyReduce( <gb>, <v> ) 
##
##  This function takes as input a right Groebner basis for a module 
##  constructed by  RightProjectiveModule  and an element thereof, and 
##  completely reduces the element  <v>  modulo this Groebner basis. 
##
InstallOtherMethod( CompletelyReduce,
  "for a path algebra module Groebner basis",
  IsCollsElms,
  [ IsRightPathAlgebraModuleGroebnerBasis 
    and IsPathAlgebraModuleGroebnerBasisDefaultRep,
    IsRightAlgebraModuleElement and IsPackedElementDefaultRep ], 0,
  function( gb, v )
    local fam;

    v := ExtRepOfObj(v);
    fam := ElementsFamily(FamilyObj(gb));

    return ObjByExtRep(fam, CompletelyReduce(gb, v));
  end
);

#######################################################################
##
#O  CompletelyReduce( <gb>, <v> ) 
##
##  This function takes as input a right Groebner basis  <gb>  for a 
##  module constructed by  RightProjectiveModule  and an element  <v>  
##  in the underlying  PathAlgebraModule  and completely reduces the 
##  element  <v> modulo this Groebner basis. This is the companion of 
##  operation defined above.   
##
InstallOtherMethod( CompletelyReduce,
  "for a path algebra module Groebner basis",
  true,
  [ IsRightPathAlgebraModuleGroebnerBasis 
    and IsPathAlgebraModuleGroebnerBasisDefaultRep,
    IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep ], 0,
  function( gb, v )
    local matches, c, newVec, tipMon, match, redVector, fam;

    if IsZero(v) then
      return v;
    fi;

    newVec := Zero(v);

    repeat
      c := v![2];
      tipMon := LeadingMonomial( v![1][c] );

      if IsBound(gb!.staticDictionaries[c]) then
        matches := PrefixSearch( gb!.staticDictionaries[c], tipMon );
      else
        matches := [];
      fi;

      if not IsEmpty(matches) then
        match := matches[1][2][1];
        redVector := gb!.gbasisElems[c][match];
        v := ReduceRightModuleElement( v, redVector );
      else
        newVec := newVec + LeadingTerm(v);
        v := v - LeadingTerm(v);
      fi;

    until IsZero(v);

    return newVec;

end
);


InstallMethod( PrintObj,
  "for f.p. path algebra vectors", 
  true,
  [IsFpPathAlgebraVector 
   and IsPathAlgebraVectorDefaultRep ], 0,
  function( x )
    Print( "[ ", x![1], " ]" );
  end
);

#O
InstallMethod( LiftPathAlgebraModule,
  "for a right path algebra module over a f.p. path algebra",
  true,
  [ IsPathAlgebraModule and IsRightAlgebraModuleElementCollection ], 0,
  function( M )
    local A, fpFam, I, gb, vecFam, liftedParent, parentFam, zero,
          rgbElems, gens, gen, i, Lift, KGamma, liftedVerts,
          liftedFam, elem, component, liftedM;

    A := ActingAlgebra( M );

    if not IsQuotientOfPathAlgebra( A ) then
      return M;
    fi;

    # Verify that the ideal has a Groebner basis
    fpFam := ElementsFamily(FamilyObj(A));
    I := fpFam!.ideal;
    if not HasGroebnerBasisOfIdeal( I ) then
      TryNextMethod();
    fi;

    # Make sure the Groebner basis is complete
    gb := GroebnerBasisOfIdeal( I );
    if not IsCompleteGroebnerBasis( gb ) then
      TryNextMethod();
    fi;

    Lift := function(type, x)
      return Objectify(type, [List(x![1], y -> y![1]), x![2]]);
    end;

    KGamma := fpFam!.pathAlgebra;
    vecFam := ElementsFamily(FamilyObj(M))!.underlyingModuleEltsFam;
    liftedVerts := List(vecFam!.vertices, x -> x*One(KGamma));
    liftedParent := RightProjectiveModule(KGamma, liftedVerts);
    parentFam := ElementsFamily(FamilyObj(liftedParent));
    liftedFam := parentFam!.underlyingModuleEltsFam;

    if not IsVertexProjectiveModule(M) then
      gens := List( GeneratorsOfAlgebraModule(M), 
                    x -> ObjByExtRep(parentFam, 
                    Lift(liftedFam!.defaultType, x![1]) ) );
    else
      gens := [];
    fi;

    rgbElems := Enumerator(RightGroebnerBasis(I));
    zero := Zero(KGamma);

    for i in [1..liftedFam!.vectorLen] do

      for elem in rgbElems do

        component := liftedFam!.vertices[i]*elem;

        if not IsZero(component) then
          gen := ListWithIdenticalEntries(liftedFam!.vectorLen, zero);
          gen[i] := component;
          Add(gens, Vectorize(liftedParent, gen));
        fi;

      od;

    od;

    if not IsEmpty(gens) then
      liftedM := SubAlgebraModule(liftedParent, gens);
      return liftedM;
    else
      return liftedParent;
    fi;

  end
);

#######################################
#  
#  Doesn't seem to be used.   
##  
InstallMethod( NaturalHomomorphismBySubAlgebraModule,
  "for vertex projective path algebra module and path algebra module",
  IsIdenticalObj,
  [ IsPathAlgebraModule and IsVertexProjectiveModule,
    IsPathAlgebraModule ], 0,
  function( M, N )
    local fam, hom, gb, U, A, zero, nontips, B, i, bv, gens, gen, baslist;

    if IsLeftAlgebraModuleElementCollection( M ) then
      TryNextMethod();
    fi;

    fam := NewFamily( "IsFpPathAlgebraVectorFamily", IsFpPathAlgebraVector );
    fam!.defaultType := NewType( fam, IsPathAlgebraVectorDefaultRep );
    fam!.preimageFam := ElementsFamily(FamilyObj(M))!.underlyingModuleEltsFam;
    fam!.vectorLen := fam!.preimageFam!.vectorLen;
    fam!.vertices := fam!.preimageFam!.vertices;
    fam!.componentFam := fam!.preimageFam!.componentFam;
    fam!.bottomModule := N;
    fam!.needsLift := IsQuotientOfPathAlgebra(ActingAlgebra(N));
    fam!.zeroPath := Zero(QuiverOfPathAlgebra(ActingAlgebra(N)));

    N := LiftPathAlgebraModule(N);

    # Get right groebner basis of module:
    gb := RightGroebnerBasisOfModule(N);

    if gb <> fail then

      if fam!.needsLift then
        fam!.liftedFam := ElementsFamily(FamilyObj(N))!.underlyingModuleEltsFam;
      fi;

      fam!.gb := gb;
      fam!.defaultType := NewType(fam, IsPathAlgebraVectorDefaultRep
                                       and IsNormalForm );

      SetNormalFormFunction(fam, function(fam, components, pos)
            local nf, preVec, fpFam, vecFam;

            if fam!.needsLift then
              fpFam := ElementsFamily(FamilyObj(components));
              vecFam := fam!.liftedFam;
              components := List(components, x -> x![1]);
            else
              vecFam := fam!.preimageFam;
            fi;

            preVec := PathAlgebraVectorNC(vecFam, components, pos, false);

            nf := CompletelyReduce(fam!.gb, preVec);
              
            if fam!.needsLift then
              nf![1] := List(nf![1], x->ElementOfQuotientOfPathAlgebra(fpFam, x, true));
            fi;

#Print("BANG ",Objectify( fam!.defaultType, [nf![1], nf![2]])," \n");
#            return Objectify( vecFam!.defaultType, [nf![1], nf![2]] );
            return Objectify( fam!.defaultType, [nf![1], nf![2]] );
        end );

      nontips := Filtered(BasisVectors(CanonicalBasis(M)), function(x)
            local er;
            er := ExtRepOfObj(x);
            return er![1] = NormalFormFunction(fam)(fam,er![1],er![2])![1];
        end );

      nontips := List( nontips, function(x)
            local er;
            er := ExtRepOfObj(x);
            return PathAlgebraVectorNC(fam,er![1],er![2],true);
        end );

    fi;

    A := ActingAlgebra(M);
    zero := Zero(A);
    fam!.zeroVector := PathAlgebraVectorNC(fam, 
                        ListWithIdenticalEntries(fam!.vectorLen, zero),
                        fam!.vectorLen + 1, true );

    gens := [];
    for i in [1..fam!.vectorLen] do
      gen := ListWithIdenticalEntries(fam!.vectorLen, zero);
      gen[i] := fam!.vertices[i]*One(A);
      Add(gens, PathAlgebraVectorNC( fam, gen, i, true ));
    od;

    gens := Filtered(gens, x -> not IsZero(x));
    U := RightAlgebraModuleByGenerators( A, \^, gens );
    SetIsFpPathAlgebraModule( U, true );
    SetIsPathAlgebraModule( U, true );
    SetIsWholeFamily(U, true);
    fam!.wholeModule := U;
    if IsBound(nontips) then
      Sort(nontips);
      baslist := List(nontips, x->ObjByExtRep(ElementsFamily(FamilyObj(U)), x) );

#Print(nontips,"\n",IsFreeLeftModule(U),"\n",IsPathAlgebraVectorCollection(baslist));

      B := NewBasis( U, baslist ); 
      SetBasis( U, B );
      SetDimension( U, Length(nontips) );
    fi;

    bv := BasisVectors(Basis(M));
    hom := LeftModuleHomomorphismByImagesNC( M, U, bv,
               List(bv, x -> Vectorize(U, ExtRepOfObj(x)![1]) ) );
    SetIsSurjective(hom, true);

    return hom;   

  end
);


#################################################################################
#  Does not seemed to be used!
#  
# VertexProjectivePresentation
#
# dependencies:
#   present.gi
#     1) SubAlgebraModule      [ IsFreeLeftModule and IsPathAlgebraModule,
#                                IsAlgebraModuleElementCollection and IsList]
#     2) SourceVertex          [ IsPathAlgebraVector and IsPathAlgebraVectorDefaultRep ]
#     3) RightProjectiveModule [ IsPathAlgebra, IsHomogeneousList ]
#     4) Vectorize             [ IsPathAlgebraModule, IsHomogeneousList]
#     5) NaturalHomomorphismBySubAlgebraModule, (used by the quotient P/N)
#                              [ IsPathAlgebraModule and IsVertexProjectiveModule,
#                                IsPathAlgebraModule ]
#
#   algpath.gi
#     1) LeadingMonomial [ IsElementOfMagmaRingModuloRelations ]
#     2) IsLeftUniform   [ IsElementOfMagmaRingModuloRelations ]
#
#   quiver.gi
#     1) IsZeroPath (property)
#     
InstallMethod( VertexProjectivePresentation,
  "For a path algebra and a list of lists of path algebra elements",
  IsElmsColls,
  [ IsAlgebra, IsRingElementTable ], 0,
  function( A, map )
    local gen, len, verts, gVerts, i, P, N;

    if not (IsPathAlgebra(A) or IsQuotientOfPathAlgebra(A)) then
        TryNextMethod();
    fi;

    # Make sure generators are not empty
    if IsEmpty(map) then
        Error("Usage: VertexProjectivePresentation( <A>, <map> )",
              " <map> must be nonempty.");
    fi;

    # Make sure all entries are left uniform
    for gen in map do
      if not ForAll(gen, IsLeftUniform) then
        Error("Usage: VertexProjectivePresentation( <A>, <map> )",
              " entries in <map> must be left uniform.");
      fi;
    od;

    # Verify we really have proper presentation:
    verts := List(map[1], x -> SourceVertex(LeadingMonomial(x)));

    for gen in map do
      gVerts := List(gen, x -> SourceVertex(LeadingMonomial(x)));

      for i in [1..Length(verts)] do
        if IsZeroPath(verts[i]) and not IsZeroPath(gVerts[i]) then
          verts[i] := gVerts[i];
        elif verts[i] <> gVerts[i] and not IsZeroPath(gVerts[i]) then
          Error("Usage: VertexProjectivePresentation( <A>, <map> )",
                " <map> contains mismatched starting vertices.");
        fi;
      od;

    od;

    # Check columns are non-zero:
    if not ForAll(verts, x -> not IsZeroPath(x)) then
      Error("Usage: VertexProjectivePresentation( <A>, <map> )",
            " <map> contains all zeroes in some column.");
    fi;

    # Ok, everything is good, construct the module
    P := RightProjectiveModule( A, List(verts, x -> x*One(A) ) );

    # Convert map:
    map := List( map, x -> Vectorize(P, x) );

    # Determine...
    N := SubAlgebraModule(P, map);

    # Return quotient module:
    return P/N;

  end
);

###############################################################  
#  Does not seem to be used!
#  
InstallOtherMethod( VertexProjectivePresentation,
  "For a path algebra and a list of lists of path algebra elements",
  function(F1, F2, F3) 
    return IsElmsColls(F1, F2) and IsIdenticalObj(F1,F3); 
  end,
  [ IsAlgebra, IsRingElementTable, IsList ], 0,
  function( A, map, verts )
    local gen, len, gVerts, i, P, N;

    if not (IsPathAlgebra(A) or IsQuotientOfPathAlgebra(A)) then
      TryNextMethod();
    fi;

    # Make sure generators are not empty
    if IsEmpty(map) then
      Error("Usage: VertexProjectivePresentation( <A>, <map>, <verts> ) <map> must be nonempty.");
    fi;

    # Make sure generators have the same length
    len := Length(verts);
    if not ForAll(map, x -> Length(x) = len ) then
      Error("Usage: VertexProjectivePresentation( <A>, <map>, <verts> ) rows",
            " in <map> and <verts> must have the same length.");
    fi;

    # Make sure all entries are left uniform
    for gen in map do
      if not ForAll(gen, IsLeftUniform) then
        Error("Usage: VertexProjectivePresentation( <A>, <map>, <verts> )",
              " entries in <map> must be left uniform.");
      fi;
    od;

    if not ForAll(verts, function(x)
                   return x = Tip(x) and IsQuiverVertex(TipMonomial(x))
                   and One(TipCoefficient(x)) = TipCoefficient(x);
                 end )
      then
        Error("<verts> should be a list of embedded vertices");
    fi;

    verts := List(verts, TipMonomial);

    # Verify vertices
    for gen in map do
      gVerts := List(gen, x -> SourceVertex(LeadingMonomial(x)));

      for i in [1..Length(verts)] do

        if verts[i] <> gVerts[i] and not IsZeroPath(gVerts[i]) then
          Error("Usage: VertexProjectivePresentation( <A>, <map>, <verts> )",
                " <map> contain mismatched starting vertices.");
        fi;

      od;
    od;

    # Ok, everything is good, construct the module
    P := RightProjectiveModule( A, List(verts, x -> x*One(A) ) );
    map := List( map, x -> Vectorize(P, x) );
    N := SubAlgebraModule(P, map);

    return P/N;

  end
);

###############################################################  
#  Maybe not used?
# 
InstallOtherMethod( Vectorize,
  "for path algebra modules and a homogeneous list",
  true,
  [IsFpPathAlgebraModule, IsHomogeneousList], 0,
  function( M, vector )
    local vecElmFam, modElmFam;

    vecElmFam := FamilyObj(ExtRepOfObj(Zero(M)));
    modElmFam := ElementsFamily(FamilyObj(M));

    return ObjByExtRep(modElmFam, PathAlgebraVector(vecElmFam, vector));
  end
);


###############################################################  
#  Not gotten it to work. Unsure about what to input.
# 
HOPF_FpPathAlgebraModuleHom := function(A, M, N)
  local B, bv, F, Q, verts, partVec, vec, mon, i, j, k, l, rt,
      arrows, src, tgt, map, maps, Mfam, mapArrows, eqns, uGens,
      rows, cols, currentRow, currentCol, colLabels, cokernel,
      tensElems, Ndim, endo, componentBasis, nf, lc, gens,
      r, s, coeffs, term, ebv, Nfam, mbv, embv, Mdim, tmpfunc,bvpart;

  rt := Runtime();
  Info( InfoPathAlgebraModule, 1, "Starting the computation" );

  F := LeftActingDomain(A);
  Q := QuiverOfPathAlgebra(A);
  B := Basis(N);

  Info(InfoPathAlgebraModule, 1, "START: Partitioning vertices: ", Runtime()-rt);

  bv := BasisVectors(B);
  ebv := List(bv, ExtRepOfObj);
  verts := VerticesOfQuiver(Q);
  partVec := List(verts, x -> []);

  for i in [1..Length(bv)] do
    mon := ExtRepOfObj( TargetVertex( ebv[i] ) )[1];
    Add(partVec[mon], i);
  od;

  Info(InfoPathAlgebraModule, 1, "FINISH: Partitioning vertices: ", Runtime()-rt);
  Info(InfoPathAlgebraModule, 1, "START: Linear maps: ", Runtime()-rt);

  arrows := GeneratorsOfAlgebraWithOne(A);
#Print(arrows,"\n");
  maps := [];

  for i in [1..Length(arrows)] do
    src := ExtRepOfObj(SourceOfPath( LeadingMonomial(arrows[i]) ))[1];
    tgt := ExtRepOfObj(TargetOfPath( LeadingMonomial(arrows[i]) ))[1];
    bvpart := bv{partVec[src]};
#Print(src," ",tgt,": ",IsBasisOfPathAlgebraVectorSpace(B)," ",IsPathAlgebraVector(bvpart[1]^(arrows[i]))," ",partVec[tgt],"\n");
    tmpfunc := function (x)
                 local xx;
                 xx := Coefficients(B, x^(arrows[i]));
                 return xx{partVec[tgt]};
               end;
#    maps[i] := List( bvpart, x ->
#                     Coefficients(B, x^(arrows[i])){partVec[tgt]} );
    maps[i] := List( bvpart, tmpfunc );
  od;

  Info(InfoPathAlgebraModule, 1, "FINISH: Linear maps: ", Runtime()-rt);
  Info(InfoPathAlgebraModule, 1, "START: Linear equations: ", Runtime()-rt);

  Mfam := ElementsFamily(FamilyObj(M))!.underlyingModuleEltsFam;
  Nfam := ElementsFamily(FamilyObj(N))!.underlyingModuleEltsFam;
  uGens := UniformGeneratorsOfModule(Mfam!.bottomModule);
  cols := List( Mfam!.vertices, x -> Length(partVec[ ExtRepOfObj(x)[1] ]) );
  rows := List( uGens, 
                function(x)
                  local v;
                  v := TargetVertex( ExtRepOfObj(x) );
                  return Length(partVec[ ExtRepOfObj(v)[1] ]);
                end );

  # rows and cols are reversed because we directly
  # construct the transposed eqns matrix to avoid
  # unnecessary copying.    
  eqns := NullMat( Sum(cols), Sum(rows), F );

  currentRow := 1;
  for i in [1..Length(uGens)] do
    currentCol := 1;
    for j in [1..Mfam!.vectorLen] do
      term := ExtRepOfObj(uGens[i])![1][j];

      if not IsZero(term) then
        eqns{[currentCol..currentCol+cols[j]-1]}
            {[currentRow..currentRow+rows[i]-1]}:=
                    MappedExpression(term, arrows, maps);
      fi;

      currentCol := currentCol + cols[j];
    od;
    currentRow := currentRow + rows[i];
  od;

  Info(InfoPathAlgebraModule, 1, "FINISH: Linear equations: ", Runtime()-rt);
  Info(InfoPathAlgebraModule, 1, "START: CoKernelnel: ", Runtime()-rt);
  Info(InfoPathAlgebraModule, 2, "Dimensions: ", DimensionsMat(eqns));

  cokernel := NullspaceMatDestructive( eqns );

  Info(InfoPathAlgebraModule, 1, "FINISH: CoKernelnel: ", Runtime()-rt);
  Info(InfoPathAlgebraModule, 1, "START: Interpret basis: ", Runtime()-rt);

  mbv := BasisVectors(Basis(M));
  componentBasis := List(GeneratorsOfAlgebraModule(M), 
                             x -> PositionSet(mbv, x) );

  colLabels := Flat( List( [1..Mfam!.vectorLen], 
            x -> ListWithIdenticalEntries( cols[x], componentBasis[x] ) ) );
  colLabels := [colLabels, Flat( List( Mfam!.vertices, 
            x -> partVec[ ExtRepOfObj(x)[1] ]))];

  Ndim := Dimension(N);
  Mdim := Dimension(M);
  tensElems := List( [1..Length(cokernel)], x -> [] );

  endo := [];
  for i in [1..Length(cokernel)] do
    endo[i] := NullMat(Mdim, Ndim, F);

    for j in [1..Sum(cols)] do
      k := colLabels[1][j];
      l := colLabels[2][j];
      if not IsZero(cokernel[i][j]) then
        Add(tensElems[i], [cokernel[i][j], k, l]);
      fi;
    od;
  od;

  Info(InfoPathAlgebraModule, 1, "FINISH: Interpret basis: ", Runtime()-rt);
  Info(InfoPathAlgebraModule, 1, "START: Homomorphisms: ", Runtime()-rt);

  embv := List(mbv, ExtRepOfObj);
  for j in [1..Mdim] do
    map := MappedExpression( LeadingComponent(embv[j]), arrows, maps );
    src := ExtRepOfObj( SourceVertex(LeadingMonomial(LeadingComponent(embv[j]))))[1];
    tgt := ExtRepOfObj( TargetVertex(LeadingMonomial(LeadingComponent(embv[j]))))[1];
    for i in [1..Length(cokernel)] do
      for k in tensElems[i] do
        if LeadingPosition( embv[j] ) = LeadingPosition( embv[ k[2] ] ) 
                     and ExtRepOfObj(TargetVertex( ebv[k[3]] ))[1] = src
           then
             endo[i][j]{partVec[tgt]} := endo[i][j]{partVec[tgt]} + 
                                     k[1]*map[PositionSet(partVec[src], k[3])];
        fi;
      od;    
    od;
  od;

  for map in endo do
    MakeImmutable(map);
    ConvertToMatrixRep(map, F);
  od;

  Info(InfoPathAlgebraModule, 1, "FINISH: Homomorphisms: ", Runtime()-rt);
  Info(InfoPathAlgebraModule, 1, "START: Hom gens: ", Runtime()-rt);

  gens := List(endo, img -> LeftModuleHomomorphismByMatrix(Basis(M), img, Basis(N)) );

  for i in [1..Length(gens)] do
    SetFilterObj(gens[i], IsAlgebraModuleHomomorphism);
  od;

  Info(InfoPathAlgebraModule, 1, "FINISH: Hom gens: ", Runtime()-rt);

  return gens;
end;

###############################################################  
#  Not gotten it to work. Unsure about what to input.
# 
InstallMethod(Hom,
  "for a path algebra and two f.p. path algebra modules",
  true,
  [ IsRing, IsFpPathAlgebraModule, IsFpPathAlgebraModule ], 0,
  function( A, M, N )
    local F, gens;

    if not IsIdenticalObj(A, ActingAlgebra(M)) then
      Error("The path algebra is not compatible with the module.");
    fi;
 
    F := LeftActingDomain(A);
    gens := HOPF_FpPathAlgebraModuleHom(A, M, N);

    return VectorSpace( F, gens, "basis" );
  end
);

###############################################################  
#  Not gotten it to work. Unsure about what to input.
# 
InstallMethod( End,
  "for a path algebra and a path algebra module",
  true,
  [IsRing, IsFpPathAlgebraModule], 0,
  function( A, M )
    local F, gens;

    if not IsIdenticalObj(A, ActingAlgebra(M)) then
      Error("The path algebra is not compatible with the module.");
    fi;
 
    F := LeftActingDomain(A);
    gens := HOPF_FpPathAlgebraModuleHom(A, M, M);

    return AlgebraWithOne( F, gens, "basis" );
  end
);

###############################################################  
#  Doesn't seem to be used!
# 
InstallMethod( ImagesSource,
  "for a linear g.m.b.i that is an algebra module homomorphism",
  true,
  [ IsGeneralMapping and IsLinearGeneralMappingByImagesDefaultRep
    and IsAlgebraModuleHomomorphism ],
  function(map)
    if IsBound( map!.basisimage ) then
      return UnderlyingLeftModule( map!.basisimage );
    else
      return SubAlgebraModule(Range(map), map!.genimages);
    fi;
  end
);

###############################################################  
#  Doesn't seem to be used!
# 
InstallMethod( ImagesSource,
  "for a linear g.m.b.i that is an algebra module homomorphism",
  true,
  [ IsGeneralMapping and IsLinearMappingByMatrixDefaultRep
    and IsAlgebraModuleHomomorphism ],
  function(map)
    local images;

    if IsBound( map!.basisimage ) then
      return UnderlyingLeftModule( map!.basisimage );
    else
      images := List( map!.matrix, 
                      row -> LinearCombination( map!.basisrange, row ));
                           
      return SubAlgebraModule(Range(map), images);
    fi;
  end
);
  
#######################################################################
##
#O  CompletelyReduceGroebnerBasisForModule( <GB> ) 
##
##  This function takes as input an right Groebner basis for a vertex 
##  projective module or a submodule thereof, an constructs completely
##  reduced right Groebner from it.
##  
InstallMethod( CompletelyReduceGroebnerBasisForModule,
        "for complete groebner bases",
        true,
        [IsPathAlgebraModuleGroebnerBasis], 0,
        function( GB )
    local grb, i, n, tip, rel;

    grb := TipReduce(BasisVectors(GB));
    # now completely reduce each relation
    n := Length(grb);
    for i in [1..n] do
        rel := grb[i];
        tip := LeadingTerm(rel);
        grb[i] := tip + CompletelyReduce(GB, rel - tip);
    od;
    SetIsCompletelyReducedGroebnerBasis(grb, true);

    return grb;
end
);

#######################################################################
##
#O  ProjectivePathAlgebraPresentation( <M> ) 
##
##  This function takes as input a  PathAlgebraMatModule  and constructs
##  a projective presentation of this module over the path algebra over
##  which it is defined, ie. a projective resolution of length 1. It 
##  returns a list of five elements: (1) a projective module  P  over 
##  the path algebra, which modulo the relations induced the projective 
##  cover of  <M>, (2) a submodule  U  of  P  such that P/U is isomorphic
##  to  <M>, (3) module generators of  P, (4) module generators for  U
##  which forms a completely reduced right Groebner basis for  U, and 
##  (5) a matrix with entries in the path algebra which gives the map 
##  from  U  to  P, if  U  were considered a direct sum of vertex projective
##  modules over the path algebra. 
##
InstallMethod( ProjectivePathAlgebraPresentation,
        "for a finite dimensional module over a path algebra",
        true, 
        [ IsPathAlgebraMatModule ], 0,
        function( M ) 

    local   A,  B,  Q,  fam,  KQ,  gens,  K,  I,  gb,  gbb,  rtgbofI,  
            num_vert,  verticesinpathalg,  vertices,  verticesinalg,  
            arrows,  BB,  B_M,  G,  MM,  projective_cover,  i,  r,  j,  
            matrix,  solutions,  Solutions,  projective_vector,  s,  
            vertex_list,  P,  f0,  PI_list,  zero_in_P,  g,  temp,  
            syzygyoverKQ,  U,  RG,  redRG,  f1,  mat;

    A := RightActingAlgebra(M);
    B := CanonicalBasis(A);
    Q := QuiverOfPathAlgebra(A); 
    fam := ElementsFamily(FamilyObj( UnderlyingLeftModule( B ) ));
    KQ := OriginalPathAlgebra(A);
    gens := GeneratorsOfTwoSidedIdeal(fam!.ideal);
    K := LeftActingDomain(A);
    I := Ideal(KQ,gens);
    gb := GroebnerBasisFunction(KQ)(gens,KQ);
    gbb := GroebnerBasis(I,gb);
    rtgbofI := RightGroebnerBasis(I);

    num_vert := NumberOfVertices(Q);
    verticesinpathalg := GeneratorsOfAlgebra(KQ){[1..num_vert]};
    vertices := VerticesOfQuiver(Q);
    if IsPathAlgebra(A) then 
        verticesinalg := GeneratorsOfAlgebra(A){[1..num_vert]};
        arrows := GeneratorsOfAlgebra(A){[1+num_vert..num_vert+NumberOfArrows(Q)]};  
    else 
        verticesinalg := GeneratorsOfAlgebra(A){[2..1+num_vert]};
        arrows := GeneratorsOfAlgebra(A){[2+num_vert..1+num_vert+NumberOfArrows(Q)]};
    fi;

    BB := BasisOfProjectives(A);
    BB := List(BB, x -> List(x, y -> Flat(y)));
    BB := List(BB, x -> Flat(x));
    
    B_M := CanonicalBasis(M);
    G := MinimalGeneratingSetOfModule(M);
#
# Finding the products of the generators of M and the basis of the projective covering
# this generator.
#
    MM := [];
    projective_cover := List(G, x -> VertexPosition(SupportModuleElement(x)![1]));
    for i in [1..Length(G)] do
        r := projective_cover[i];
        for j in [1..Length(BB[r])] do 
            Add(MM,G[i]^BB[r][j]); 
        od;
    od;
#
#  Viewing the module as a vector space, flattening the vectors into K^dim(M).
#
    matrix := NullMat(Length(MM),Dimension(M),K);
    for i in [1..Length(MM)] do
        matrix[i] := Flat(ExtRepOfObj(MM[i])![1]);
    od;
#
#  Finding the kernel of the projective cover as a vectorspace 
#  
    solutions := NullspaceMat(matrix);
#
#  Finding the kernel of the projective cover as a submodule of the 
#  projective cover. 
#
    Solutions := [];
    for i in [1..Length(solutions)] do
        projective_vector := [];
        s := 0;
        for j in [1..Length(projective_cover)] do
            r := Length(BB[projective_cover[j]]);
            Add(projective_vector,LinearCombination(BB[projective_cover[j]],solutions[i]{[s+1..s+r]}));
            s := s + Length(BB[projective_cover[j]]);
        od;
        Add(Solutions,projective_vector);
    od;
#  
#  Finding the projective  P  over KQ inducing the projective cover.
#
    vertex_list := List(G, x -> SupportModuleElement(x)![1]![1]);
    P := RightProjectiveModule(KQ,vertex_list);
    f0 := GeneratorsOfAlgebraModule(P);
#
#  Finding  PI
#
    PI_list := [];
    zero_in_P := List(vertex_list, x -> Zero(KQ));
    for i in [1..Length(vertex_list)] do
        for g in rtgbofI do
            if vertex_list[i]*g <> Zero(KQ) then
                temp := ShallowCopy(zero_in_P);
                temp[i] := vertex_list[i];
                Add(PI_list,Vectorize(P,temp*g));
            fi;
        od;
    od;
#
#  Constructing the syzygy over KQ.
#
    syzygyoverKQ := List(Solutions, x -> List(x, y -> y![1]));
    syzygyoverKQ := List(syzygyoverKQ, x -> Vectorize(P,x));
    Append(syzygyoverKQ,PI_list);
    
    if Length(syzygyoverKQ) = 0 then
        return [P,[],[],[],[]];
    else 
        fam := ElementsFamily(FamilyObj(P));
        U := SubAlgebraModule(P,syzygyoverKQ);
        RG := RightGroebnerBasisOfModule(U);
        redRG := CompletelyReduceGroebnerBasisForModule(RG);
        f1 := List( redRG, x -> ObjByExtRep(fam, x));
        mat := TransposedMat(List(f1, x -> x![1]![1]));
        return [ P, U, f0, f1, mat, redRG ];
    fi;
end
);
