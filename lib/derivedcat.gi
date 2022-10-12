# projective and injective complexes -- property and printing

###################################################
##
#P IsProjectiveComplex( <C> )
##
## Returns if <C> is a finite complex
## with only projective objects. If <C> is infinite,
## it is not checked whether the complex actually 
## consists of projectives.
##
InstallMethod( IsProjectiveComplex,
               [ IsQPAComplex ],
               function( C )

    local i, obj;

    if( IsInt(UpperBound(C)) and IsInt(LowerBound(C)) ) then
        for i in [LowerBound(C)..UpperBound(C)] do
            obj := ObjectOfComplex(C,i);
            if( not IsProjectiveModule(obj) ) then
                return false;
            fi;
        od;
        return true;
    fi;
    return false;
end);

###################################################
##
#P IsInjectiveComplex( <C> )
##
## Returns if <C> is a finite complex
## with only injective objects. If <C> is infinite,
## it is not checked whether the complex actually 
## consists of injectives.
##
InstallMethod( IsInjectiveComplex,
               [ IsQPAComplex ],
               function( C )
    local i, obj;

    if( IsInt(UpperBound(C)) and IsInt(LowerBound(C)) ) then
        for i in [LowerBound(C)..UpperBound(C)] do
            obj := ObjectOfComplex(C,i);
            if( not IsInjectiveModule(obj) ) then
                return false;
            fi;
        od;
        return true;
    fi;
    return false;
end);

# print methods

###################################################
##
#O PrintObj( <C> )
##
## For a projective complex <C>
## Prints a projective complex such that the objects
## are displayed as "Pi" for the correct vertex i
##
## value 1: ensures that a complex which is both 
## projective and injective should be printed as projective.
InstallMethod( PrintObj,
               "for a projective complex",
               [ IsProjectiveComplex ],
               1,
               function( C )

    local list, start, stop, l, t, symbol, finitetest, upperlimit, x;
    
    # check if C is finite or not
    if (IsInt(UpperBound(C))) then
        finitetest := true;
        start := UpperBound(C);
        stop := LowerBound(C);
    else
        finitetest := false;
        start := HighestKnownDegree(C);
        stop := LowerBound(C);
    fi;
        
    symbol := "P";

    if finitetest then
        list := Reversed(DescriptionOfFiniteProjComplex(C));
    else
        list := [stop..start];
        list := Reversed(List( list, x ->
                               DescriptionOfProjComplexInDegree(C,x) ));
    fi;

    # do the printing
    if (not finitetest) then
        Print("--- ->");
    else
        Print("0 ->");
    fi;
    for l in list do
        if (not IsEmpty(l)) then
            Print(" ",start,": ");
            for t in [1..Length(l)] do
                if (t = 1) then
                    Print(symbol,l[t]);
                else
                    Print(" + ",symbol,l[t]);
                fi;
            od;
            Print(" ->");
            start := start - 1;
        fi;
    od;
    Print(" 0");

end);

###################################################
##
#O PrintObj( <C> )
##
## For an injective complex <C>
## Prints an injective complex such that the objects
## are displayed as "Ii" for the correct vertex i
##
## value 0: ensures that a complex which is both 
## projective and injective should be printed as projective.
InstallMethod( PrintObj,
               "for an injective complex",
               [ IsInjectiveComplex ],
               0,
               function( C )

    local list, start, stop, l, t, symbol, finitetest, upperlimit, x;
    # check if C is finite or not
    if (IsInt(UpperBound(C))) then
        finitetest := true;
        start := UpperBound(C);
        stop := LowerBound(C);
    else
        finitetest := false;
        start := HighestKnownDegree(C);
        stop := LowerBound(C);
    fi;
        
    symbol := "I";

    if finitetest then
        list := Reversed(DescriptionOfFiniteInjComplex(C));
    else
        list := [stop..start];
        list := Reversed(List( list, x ->
                               DescriptionOfInjComplexInDegree(C,x) ));
    fi;

    # do the printing
    if (not finitetest) then
        Print("--- ->");
    else
        Print("0 ->");
    fi;
    for l in list do
        if (not IsEmpty(l)) then
            Print(" ",start,": ");
            for t in [1..Length(l)] do
                if (t = 1) then
                    Print(symbol,l[t]);
                else
                    Print(" + ",symbol,l[t]);
                fi;
            od;
            Print(" ->");
            start := start - 1;
        fi;
    od;
    Print(" 0");

end);


# projective resolution of a complex

###################################################
##
#O ProjectiveResolutionOfComplex( <C> )
##
## for a complex C in D^b(mod A) where mod A
## may have infinite global dimension.
##
## The result is a new complex C' (possibly of the
## infinite type), which is quasi-isomorphic to C. 
## C' is an element of the homotopty category K^(b)(P).
##
InstallMethod( ProjectiveResolutionOfComplex,
               [ IsQPAComplex ],
               function( C )

    local cat, dList, tList, t1, PB1, PC1, pos, ker, h, PB, PC, lastKernel, i,
          PCompl, j, tempCompl, projres, isos, reducedModules, v, d, u, diffs,
          pospart, newpart, negpart, newdiffs;

    # if C is projective, return C itself:
    if IsProjectiveComplex(C) then
        return C;
    fi;


    cat := CatOfComplex(C);
    i := LowerBound(C);
    # where the loop should stop
    if(IsInt(UpperBound(C))) then
        j := UpperBound(C) - LowerBound(C) + 1;
    else
        j := HighestKnownDegree(C)-LowerBound(C)+1;
    fi;
    dList := [];
    tList := [];

    # special case for stalk complexes
    if (LengthOfComplex(C) = 1) then
        PCompl := ProjectiveResolution(ObjectOfComplex(C,i));
        PCompl := BrutalTruncationBelow(PCompl, 0);
        PCompl := ShiftUnsigned(PCompl, -i);
        SetIsProjectiveComplex(PCompl,true);
        return PCompl;
    fi;

    # do the first position
    t1 := ProjectiveCover(ObjectOfComplex(C,i));
    Append(tList, [t1]);
    PB1 := PullBack(DifferentialOfComplex(C,i+1) , t1);
    PC1 := ProjectiveCover(Source(PB1[1]));
    Append(dList, [PC1*PB1[1]]);
    Append(tList, [PC1*PB1[2]]);

    # check if the first differential is an isomorphism
    isos := CheckForIsomorphism(dList[1]);
    while not(IsEmpty(isos)) do
        reducedModules := [DirectSumMinusSummands(Source(dList[1]),[isos[3]]),
                           DirectSumMinusSummands(Range(dList[1]),[isos[4]])];
        v := ReduceDifferential(dList[1],isos,[Source(reducedModules[1][1]),
                                                  Source(reducedModules[2][1])]);
        tList[2] := ReduceTMap(tList[2],dList[1],isos, reducedModules[1][1]);
        tList[1] := ReducePreviousTMap(tList[1],reducedModules[2][1]);
        dList[1] := v;
        isos := CheckForIsomorphism(dList[1]);
    od;

    # build the rest of the complex, until we get zero
    for pos in [1..j] do
        ker := KernelInclusion(dList[pos]);
        h := ker*tList[pos+1];

        # note: the following check shouldn't really be in this
        # method, but is needed due to the fact that DirectSumOfQPAModules
        # can't handle zero objects.
        if(IsZero(DimensionVector(Source(ker))) and 
           IsZero(DimensionVector(ObjectOfComplex(C,pos+1+i)))) then
            PCompl := FiniteComplex(cat, i+1, dList);
            SetIsProjectiveComplex(PCompl, true);
            return PCompl;
        fi;

        PB := PullBack(DifferentialOfComplex(C,pos+1+i), h);
        PC := ProjectiveCover(Source(PB[1]));

        # if we did get zero: return a finite complex
        if(IsZero(DimensionVector(Source(PC)))) then
            PCompl := FiniteComplex(cat, i+1, dList);
            SetIsProjectiveComplex(PCompl, true);
            return PCompl;
        fi;

        Append(dList, [PC*PB[1]*ker]);
        Append(tList, [PC*PB[2]]);
        
        # check if the newly created differential includes isomorphisms
        isos := CheckForIsomorphism(dList[pos + 1]);
        while not(IsEmpty(isos)) do
            reducedModules := [DirectSumMinusSummands(Source(dList[pos+1]),[isos[3]]),
                               DirectSumMinusSummands(Range(dList[pos+1]),[isos[4]])];
            v := ReduceDifferential(dList[pos+1],isos,[Source(reducedModules[1][1]),
                                                          Source(reducedModules[2][1])]);
            tList[pos+2] := ReduceTMap(tList[pos+2],dList[pos+1],isos, reducedModules[1][1]);
            tList[pos+1] := ReducePreviousTMap(tList[pos+1],reducedModules[2][1]);
            dList[pos+1] := v;
            dList[pos] := ReducePreviousDifferential(dList[pos],reducedModules[2][1]);
            isos := CheckForIsomorphism(dList[pos + 1]);
        od;

    od;
    # return an infinite complex, where the part following the above
    # is the projective resolution of the kernel of the last differential.

    tempCompl := FiniteComplex(cat, i+1, dList);
    tempCompl := SyzygyTruncation(tempCompl, j+1+i);
    projres := ProjectiveResolution(ObjectOfComplex(tempCompl, j+2+i));
    PCompl := YonedaProduct(projres, tempCompl);

    # take care of isos in the "overlap differential" (note: checks for only 
    # one iso here … possible bug)
    isos := CheckForIsomorphism(DifferentialOfComplex(PCompl, pos+2+i));
    if( not IsEmpty(isos)) then
        d := DifferentialOfComplex(PCompl, pos+2+i);
        reducedModules := [DirectSumMinusSummands(Source(d),[isos[3]]),
                           DirectSumMinusSummands(Range(d),[isos[4]])];
        v := ReduceDifferential(d,isos,[Source(reducedModules[1][1]),
                                           Source(reducedModules[2][1])]);
        u := ReducePreviousDifferential(DifferentialOfComplex(PCompl, pos+1+i),
                                        reducedModules[2][1]);
        # edit the "overlap differential" and the previous, make new complex
        diffs := DifferentialsOfComplex(PCompl);
        pospart := PositivePartFrom( diffs, pos+3+i);
        newpart := FiniteInfList( pos+1+i, [u,v]);;
        negpart := NegativePartFrom( diffs, pos + i);
        newdiffs := InfConcatenation(pospart, newpart, negpart);
        PCompl := ComplexByDifferentialList(cat, newdiffs);
    fi;
    SetIsProjectiveComplex(PCompl, true);
    return PCompl;
end);

###################################################
##
#O CheckForIsomorphism( <d> )
##
## for a map <d>, typically a differential, between
## two modules (+Pi) and (+Qj), each Pi and Qj indec 
## (but not necessarily non-iso) projective.
##
## Checks if the _same_ indec projective, say P, 
## appears as a direct summand in both (+Pi) and (+Qj).
## This can possibly happen for more than one summand,
## but as soon as one is found, it is returned.
##
## Returns a list consisting of
##  - the (described above) map from P --> P
##  - the vertex number associated to P (an integer)
##  - the position of P in (+Pi) (an integer)
##  - the position of Pi in (+Qj) (an integer)
##
InstallMethod( CheckForIsomorphism,
               [ IsPathAlgebraMatModuleHomomorphism ],
               function( d )

    local descr1, descr2, mutualIndecs, posList1, posList2, i,
          allMaps, inclusions, projections, j, k, u, u2,
          identities;

    descr1 := DescriptionOfProjectiveModule(Source(d));
    descr2 := DescriptionOfProjectiveModule(Range(d));
    
    # check if there are common indec projectives (+Pi) and (+Qj)
    mutualIndecs := Intersection(descr1,descr2);
    if(IsEmpty(mutualIndecs)) then
        return [];
    fi;
    # find the positions of those in the modules:
    posList1 := List(mutualIndecs, i -> Positions(descr1, i));
    posList2 := List(mutualIndecs, i -> Positions(descr2, i));

    # find all maps between isomorphic direct summands of (+Pi) and (+Qj)
    allMaps := FindAllMapComponents(d);

    identities := [];
    for i in [1..Length(mutualIndecs)] do
        for j in [1..Length(posList1[i])] do
            for k in [1..Length(posList2[i])] do
                u2 := allMaps[posList1[i][j]][posList2[i][k]];
                if IsIsomorphism(u2) then
                    return [u2,mutualIndecs[i],posList1[i][j],posList2[i][k]];
                fi;
            od;
        od;
    od;
    identities := Filtered(identities, i -> IsIsomorphism(i[1]));

    return identities;
end);

###################################################
##
#O ReduceDifferential( <d>, <list>, <reducedModules> )
##
## for a differential <d> between two modules (+Pi) 
## and (+Qj), each Pi and Qj indec (but not neces-
## sarily non-iso) projective, a list <list> which is 
## output from CheckForIsomorphism(d) and a list 
## <reducedModules> consisting of (+Pi~) and (+Qj~).
##
## Returns: a differential (+Pi~) --> (+Qj~) which 
## is reduced with respect to the map in <list>.
##
InstallMethod( ReduceDifferential,
               [ IsPathAlgebraMatModuleHomomorphism, IsList, IsList ],
               function( d, list, reducedModules )
    local n, m, allMaps, newMaps, map, leftMap, rightMap, inverse, A,
          newAllMaps;
    allMaps := FindAllMapComponents(d);
    n := Length(allMaps);
    m := Length(allMaps[1]);
    A := RightActingAlgebra(Source(d));

    # check for zero: if the new modules are zero, return zero map.
    if (IsZero(DimensionVector(reducedModules[2]))) then
        return ZeroMapping(reducedModules[1],reducedModules[2]);
    elif(IsZero(DimensionVector(reducedModules[1]))) then
        return ZeroMapping(reducedModules[1],reducedModules[2]);
    fi;

    # create the new maps: first, remove row j
    newMaps := Concatenation(allMaps{[1..(list[3]-1)]},
                             allMaps{[(list[3]+1)..n]});
    # second, remove column i
    newMaps := List(newMaps, map -> Concatenation(map{[1..(list[4]-1)]},
                                                  map{[(list[4]+1)..m]}));
    

    # compute the maps that is to be subtracted from newMaps
    leftMap := List(allMaps, map -> map{[list[4]]});
    Remove(leftMap,list[3]);
    rightMap := [Concatenation(allMaps[list[3]]{[1..(list[4]-1)]},
                              allMaps[list[3]]{[(list[4]+1)..m]})];
    inverse := InverseOfIsomorphism(allMaps[list[3]][list[4]]);

    # compute all map components for the new differential
    newAllMaps := newMaps - leftMap*inverse*rightMap;

    return DirectSumProjections(reducedModules[1])*newAllMaps*
           DirectSumInclusions(reducedModules[2]);
end);

###################################################
##
#O ReducePreviousDifferential( <d>, <inclusion> )
##
## for a differential <d> between two modules (+Qj) 
## and (+Rk), each Qj and Rk indec (but not neces-
## sarily non-iso) projective, and a map <inclusion> which 
## is the inclusion of the new source (+Qj~) into
## (+Qj).
## 
## Returns the reduced differential of <d>.
##
InstallMethod( ReducePreviousDifferential,
               [ IsPathAlgebraMatModuleHomomorphism,
                 IsPathAlgebraMatModuleHomomorphism ],
               function( d, inclusion )
    return inclusion*d;

end);

###################################################
##
#O ReduceTMap( <t>, <d>, <list>, <incl> )
##
## Similar for <t> as ReduceDifferential is for <d>.
## Input: <t> is the degree i component of the qis
## between C and C', where C' is projective resolution
## of C. <d> is the differential of C' in the same
## degree. <list> is output from CheckForIsomorphism(<d>),
## and <incl> is the inclusion of the reduced version of
## Source(d) into Source(d). 
##
InstallMethod( ReduceTMap,
               [ IsPathAlgebraMatModuleHomomorphism,
                 IsPathAlgebraMatModuleHomomorphism, IsList,
                 IsPathAlgebraMatModuleHomomorphism ],
               function( t, d, list, incl )
    local tQ, tP, newd, inclusions, projections, allMaps, leftMap, inverse, newAllMaps, subtract;
    tQ := incl*t;
    tP := DirectSumInclusions(Source(d))[list[3]]*t;
    newd := incl*d;
    allMaps := FindAllMapComponents(newd);
    leftMap := List(allMaps, map -> map{[list[4]]});
    inverse := InverseOfIsomorphism(FindAllMapComponents(d)[list[3]][list[4]]);
    projections := ShallowCopy(DirectSumProjections(Source(d)));
    Remove(projections,list[3]);
    subtract := incl*projections*leftMap*inverse*tP;
    if IsEmpty(subtract) then
        newAllMaps := tQ;
    else
        newAllMaps := tQ - subtract;
    fi;
    return newAllMaps[1];
end);

###################################################
##
#O ReducePreviousTMap( <t>, <incl> )
##
## Similar for <t> as ReducePreviousDifferential is 
## for <d>. 
##
InstallMethod( ReducePreviousTMap,
               [IsPathAlgebraMatModuleHomomorphism,
                IsPathAlgebraMatModuleHomomorphism],
               function( t, incl)

    return incl*t;
end);


# tau of complex

######################################################
##
#O TauOfComplex( <C> )
##
## <C> is a bounded complex over.
##
## Computes tau of a complex in D^b(mod A), where A has
## finite global dimension.
##
InstallMethod( TauOfComplex,
               [ IsQPAComplex ],
               function( C )
    
    local cat, projVersion, tau, injVersion;

    projVersion := Shift(ProjectiveResolutionOfComplex(C),1);
    injVersion := ProjectiveToInjectiveComplex(projVersion);
    injVersion := CutComplexAbove(injVersion);
    tau := CutComplexAbove(ProjectiveResolutionOfComplex(injVersion));
    SetIsProjectiveComplex(tau,true);
    
    return tau;

end);

######################################################
##
#O ProjectiveToInjectiveComplex( <C> )
##
## <C> is a (possibly finite) complex of projectives.
##
## Applies DHom_{A}(-,A) to <C>. Output is a complex
## of injectives. This is infinite if and only if <C>
## is infinite. If <C> is discovered to be finite, 
## the next method is used.
##
InstallMethod( ProjectiveToInjectiveComplex,
               [ IsQPAComplex ],
               function( PCompl )
    
    local A, cat, start, stop, injectives, inj1, inj2,
          descr1, descr2, i, computeDifferential,
          middle, f, u, mats, injmap, ICompl;

    # check whether PCompl is finite, in which case use below method
    if (IsInt(UpperBound(PCompl))) then
        return ProjectiveToInjectiveFiniteComplex(PCompl);
    fi;

    # find information about the complex
    start := LowerBound(PCompl);
    stop := HighestKnownDegree(PCompl);
    if start = stop then
        stop := stop + 1;
    fi;
    A := RightActingAlgebra(ObjectOfComplex(PCompl, start));
    cat := CatOfRightAlgebraModules(A);
    injectives := IndecInjectiveModules(A);

    # construct the known part of the complex
    middle := [];
    descr1 := DescriptionOfProjComplexInDegree(PCompl, start);
    inj1 := List(descr1, x -> injectives[x]);
    inj1 := DirectSumOfQPAModules(inj1);

    for i in [(start+1)..stop] do
        descr2 := DescriptionOfProjComplexInDegree(PCompl, i);
        if IsEmpty(descr2) then
            break;
        fi;
        # apply DHom_{A}(-,A) to the differentials
        f := DifferentialOfComplex(PCompl, i);
        u := StarOfMapBetweenProjectives(f,descr2,
                                         descr1);
        mats := MatricesOfDualMap(u);
        
        # constructing the correct injectives
        inj2 := List(descr2, x -> injectives[x]);
        inj2 := DirectSumOfQPAModules(inj2);

        # constructing the map between the injectives
        injmap := RightModuleHomOverAlgebra(inj2,inj1,mats);
        Append(middle, [injmap]);

        # updating injectives
        inj1 := inj2;
        descr1 := descr2;
    od;

    # the function finding later differentials ("unknown" part of the complex)
    computeDifferential := function( C,i )
        local rangeobj,descr1,descr2, f, u, mats, inj1, injmap, x;
        rangeobj := ObjectOfComplex(C,i-1);
        descr1 := DescriptionOfProjComplexInDegree(PCompl, i);
        descr2 := DescriptionOfProjComplexInDegree(PCompl, i-1);
        f := DifferentialOfComplex(PCompl, i);
        if IsZero(f) then
            return cat.zeroMap(cat.zeroObj, rangeobj);
        fi;
        u := StarOfMapBetweenProjectives(f, descr1, descr2);
        mats := MatricesOfDualMap(u);
        inj1 := List(descr1, x -> injectives[x]);
        inj1 := DirectSumOfQPAModules(inj1);
        return RightModuleHomOverAlgebra(inj1, rangeobj, mats);
    end;
    ICompl := Complex( cat,
                       start+1,
                       middle,
                       [ "pos", computeDifferential, true ],
                       "zero" );

    SetIsInjectiveComplex(ICompl, true);
    return ICompl;
end);

######################################################
##
#O ProjectiveToInjectiveFiniteComplex( <C> )
##
## <C> is a finite complex of projectives.
##
## Applies DHom_{A}(-,A) to <C>. Output is a finite complex
## of injectives. Note that the previous (and more general)
## method checks for infinity, so there is no need for
## explicitly calling this method.
##
InstallMethod( ProjectiveToInjectiveFiniteComplex,
               [ IsQPAComplex ],
               function( PCompl )

    local A, cat, start, stop, descr, maplist, i, f, u, mats, x, inj1,
          inj2, injmap, injectives, ICompl;
    start := LowerBound(PCompl);
    stop := UpperBound(PCompl);

    A := RightActingAlgebra(ObjectOfComplex(PCompl,start));
    cat := CatOfRightAlgebraModules(A);
    injectives := IndecInjectiveModules(A);

    # get a description of which projectives are in the complex
    descr := DescriptionOfFiniteProjComplex(PCompl);
    maplist := [];

    # note: need special case for when PCompl is stalk complex!
    if (LengthOfComplex(PCompl) = 1) then
        ICompl := StalkComplex(cat, injectives[descr[1][1]], start);
        SetIsInjectiveComplex(ICompl, true);
        return ICompl;
    fi;
    
    # the first injective
    inj1 := List(descr[1], x -> injectives[x]);
    inj1 := DirectSumOfQPAModules(inj1);


    for i in [(start+1)..stop] do
#        Print("ser nå på avbildningen i grad ", i , ".\n");
        f := DifferentialOfComplex(PCompl, i);
        u := StarOfMapBetweenProjectives(f,descr[i+1-start],descr[i-start]);
        mats := MatricesOfDualMap(u);
        
        # constructing the correct injectives
        inj2 := List(descr[i+1-start], x -> injectives[x]);
        inj2 := DirectSumOfQPAModules(inj2);

        # constructing the map between the injectives
        injmap := RightModuleHomOverAlgebra(inj2,inj1,mats);
        Append(maplist, [injmap]);

        # updating injectives
        inj1 := inj2;
    od;
    ICompl := FiniteComplex(cat, start+1, maplist);
    SetIsInjectiveComplex(ICompl, true);
    return ICompl;

end);

######################################################
##
#O StarOfMapBetweenProjectives( <f>, <list_i>, <list_j> )
##
## <f> is a map between two projective modules (+Pi) and (+Qj),
## where <list_i> and <list_j> give the indecomposable 
## summands of (+Pi) and (+Qj) (in terms of vertices of the quiver).
##
## Output is a map between the corresponding projectives
## in the opposite algebra, that is, Hom_{A}(-,A) is applied.
##
## This function is supposed to be used only when called
## from the "ProjectiveToInjectiveComplex"-methods. It cannot,
## for instance, handle zero maps.
##
InstallMethod( StarOfMapBetweenProjectives,
               [ IsPathAlgebraMatModuleHomomorphism, IsList, IsList ],
               function( f, list_i, list_j )
    local maplist, u, j, i, map, source, projections, range,
          inclusions, outermaplist;
      
    maplist := [];
    map := 0;
    if(not Length(list_i) > 1 ) then
        # P module indec
        map := StarOfMapBetweenIndecProjectives(f, list_i[1], list_j);
    else
        # P module not indec
        map := StarOfMapBetweenDecompProjectives(f, list_i, list_j);
    fi;
    return map; 
end);

######################################################
##
#O StarOfMapBetweenIndecProjectives( <f>, <i>, <j_list> )
##
## <f> is a map between two projective modules Pi and +(Qj),
## where Pi is indecomposable and <j_list> gives the indecom-
## posable summands of +(Qj) (in terms of vertices of the 
## quiver).
##
## Output is a map between the corresponding projectives
## in the opposite algebra, that is, Hom_{A}(-,A) is applied.
##
## This function is supposed to be used only when called
## from the "StarOfMapBetweenProjectives"-method. It cannot,
## for instance, handle zero maps.
##
InstallMethod( StarOfMapBetweenIndecProjectives,
                    [ IsPathAlgebraMatModuleHomomorphism, IsInt, IsList ],
                    function( f, i, j_list )

    local A, BP, e_i, elem, BasisOfVertex, eleminalg,
          A_op, P_op, a_op, e_i_op, e_i_op_a_op, s, u,
          vertices, jvertice, support,
          u_list, source, projections,
	  summands, summand, k;
    
    # finding the algebra and the indec projectives
    A := RightActingAlgebra(Source(f));
    BP := BasisOfProjectives(A);
 
    # the opposite algebra
    A_op := OppositeAlgebra(A);
    P_op := IndecProjectiveModules(A_op);
	  
    u_list := [];
    
    # in case Range(f) is not a direct sum:
    if not IsDirectSumOfModules(Range(f)) then
        summands := [f];
    else
        summands := f*DirectSumProjections(Range(f));
    fi;

    for k in [1..Length(j_list)] do
        if(IsZero(summands[k])) then
            Append(u_list,[ZeroMapping(P_op[j_list[k]],P_op[i])]);
        else

        # constructing f(e_i) as an element in the path algebra
            e_i := MinimalGeneratingSetOfModule(Source(summands[k]))[1];
            elem := ImageElm(summands[k], e_i);
            BasisOfVertex := BP[j_list[k]][i];
            eleminalg := LinearCombination(BasisOfVertex, elem![1]![1][i]);
            a_op := OppositePathAlgebraElement(eleminalg);
    
        # aop as element in the module
            e_i_op :=  Filtered(MinimalGeneratingSetOfModule(P_op[i]),
                                m -> not(IsZero(m![1]![1][i])))[1];
            e_i_op_a_op := e_i_op^a_op;
         # the map Qj* --> Pi*
            vertices := VerticesOfPathAlgebra(A_op);
            jvertice := vertices[j_list[k]]; #?
            support := e_i_op_a_op^jvertice;
            
        # zero maps must also be included
            if (IsZero(support)) then
                u := ZeroMapping(P_op[j_list[k]],P_op[i]);
            else
                u := HomFromProjective(support, P_op[i]);
            fi;
            Append(u_list, [u]);
        fi;
    od;
    source := List(u_list, u -> Source(u));
    source := DirectSumOfQPAModules(source);
    projections := DirectSumProjections(source);
    return projections*u_list;

end);

######################################################
##
#O StarOfMapBetweenDecompProjectives( <f>, <i_list>, <j_list> )
##
## <f> is a map between two projective modules (+Pi) and (+Qj),
## where <i_list> and <j_list> give the indecomposable summands
## of (+Pi) and (+Qj) (in terms of vertices of the quiver).
##
## Output is a map between the corresponding projectives
## in the opposite algebra, that is, Hom_{A}(-,A) is applied.
##
## This function is supposed to be used only when called
## from the "StarOfMapBetweenProjectives"-method. It cannot,
## for instance, handle zero maps.
##
InstallMethod( StarOfMapBetweenDecompProjectives,
                    [ IsPathAlgebraMatModuleHomomorphism, IsList, IsList ],
                    function( f, i_list, j_list )

    local maplist, i, range, m, inclusions, i_inclusions, A, BP, A_op, P_op, g, summands,
          k, u_list, e_i, elem, BasisOfVertex, eleminalg, a_op, e_i_op, e_i_op_a_op,
          vertices, jvertices, j, supportList, s, u, source, projections;

    maplist := [];
    i_inclusions := DirectSumInclusions(Source(f));

    for i in [1..Length(i_list)] do
        g := i_inclusions[i]*f;
        Append(maplist, [StarOfMapBetweenIndecProjectives(g, i_list[i], j_list)]);

    od;

    # cheating to get same source of the maps:
    if (Length(maplist) > 1) then
        for j in [2..Length(maplist)] do
            maplist[j] := DirectSumProjections(Source(maplist[1]))*
                          DirectSumInclusions(Source(maplist[j]))*maplist[j];
       od;
    fi;

    range := List(maplist, m-> Range(m));
    range := DirectSumOfQPAModules(range);
    inclusions := DirectSumInclusions(range);
    return maplist*inclusions;

end);

# various methods for describing proj/inj complexes
# and printing them

###################################################
##
#O DescriptionOfProjectiveModule( <M> )
##
## <M> is a projective module, either indecomposable
## (M = Pi for some vertex i) or a direct sum (+Pi) 
## of indecomposable projectives.
##
## Returns a list with the vertex numbers of the
## indecomposable direct summands of <M>.
##
InstallMethod( DescriptionOfProjectiveModule,
               [ IsPathAlgebraMatModule ],
               function( M )

    local list, comp, incls, incl;

    list := [];
    
    if IsZero(DimensionVector(M)) then
        return list;
    fi;

    if(not IsDirectSumOfModules(M)) then
        comp := CompareWithIndecProjective(M);
        Append(list, [comp]);
    else
        incls := DirectSumInclusions(M);
        for incl in incls do
            comp := CompareWithIndecProjective(Source(incl));
            Append(list, [comp]);
        od;
    fi;
    return list;
end);



###################################################
##
#O DescriptionOfProjComplexInDegree( <C>,<i> )
##
## <C> is a complex consisting only of projectives,
## and i is an integer.
##
## Returns a list with the vertex numbers of the
## projectives in degree <i> of <C>.
##
InstallMethod( DescriptionOfProjComplexInDegree,
               [ IsQPAComplex, IsInt ],
               function( C, i )
    return DescriptionOfProjOrInjComplexInDegree(C, i, true);
end);

###################################################
##
#O DescriptionOfInjComplexInDegree( <C>,<i> )
##
## <C> is a complex consisting only of injectives,
## and i is an integer.
##
## Returns a list with the vertex numbers of the
## injectives in degree <i> of <C>.
##
InstallMethod( DescriptionOfInjComplexInDegree,
               [ IsQPAComplex, IsInt ],
               function( C, i )
    return DescriptionOfProjOrInjComplexInDegree(C, i, false);
end);

###################################################
##
#O DescriptionOfProjOrInjComplexInDegree( <C>,<i>,<test> )
##
## <C> is a complex consisting only of projectives OR
## only of injectives, and <i> is an integer.
##
## <test> is a boolean, with values interpreted as
##     <test> = true  <->  C consists of projectives
##     <test> = false <->  C consists of injectives
##
## Returns a list with the vertex numbers of the
## projectives or injectives in degree <i> of <C>.
##
InstallMethod( DescriptionOfProjOrInjComplexInDegree,
               [ IsQPAComplex, IsInt, IsBool ],
               function( C, i, test )
    local obj,comp, list, incls, incl;
    
    obj := ObjectOfComplex(C, i);

    if test then
        return DescriptionOfProjectiveModule(obj);
    else    
        list := [];
        
        if IsZero(DimensionVector(obj)) then
            return list;
        fi;

        if(not IsDirectSumOfModules(obj)) then
            comp := CompareWithIndecInjective(obj);
            Append(list, [comp]);
        else
            incls := DirectSumInclusions(obj);
            for incl in incls do
                comp := CompareWithIndecInjective(Source(incl));
                Append(list, [comp]);
            od;
        fi;
        return list;
    fi;
end);

######################################################
##
#O DescriptionOfFiniteProjComplex( <C> )
##
## <C> is a finite complex consisting only of 
## projectives.
##
## Returns a list of lists, one list for each degree
## of <C>, containing the vertex numbers of the
## projectives in this degree.
##
InstallMethod( DescriptionOfFiniteProjComplex,
               [ IsQPAComplex ],
               function( C )
    return DescriptionOfFiniteProjOrInjComplex(C,true);
end);

######################################################
##
#O DescriptionOfFiniteInjComplex( <C> )
##
## <C> is a finite complex consisting only of 
## injectives.
##
## Returns a list of lists, one list for each degree
## of <C>, containing the vertex numbers of the
## injectives in this degree.
##
InstallMethod( DescriptionOfFiniteInjComplex,
               [ IsQPAComplex ],
               function( C )
    return DescriptionOfFiniteProjOrInjComplex(C,false);
end);

###################################################
##
#O DescriptionOfFiniteProjOrInjComplex( <C>,<test> )
##
## <C> is a complex consisting only of projectives OR
## only of injectives.
##
## <test> is a boolean, with values interpreted as
##     <test> = true  <->  C consists of projectives
##     <test> = false <->  C consists of injectives
##
## Returns a list of lists, one list for each degree
## of <C>, containing the vertex numbers of the
## projectives/injectives in this degree.
##
InstallMethod( DescriptionOfFiniteProjOrInjComplex,
               [ IsQPAComplex, IsBool ],
               function( C, test )

    local i,obj,comp, incls, incl, list, templist, start, stop;
    
    start := LowerBound(C);
    stop := UpperBound(C);
    list := [];
    
    if(not(IsInt(UpperBound(C)))) then
        Error("complex entered is not finite!");
    fi;

    for i in [start..stop] do
        templist := DescriptionOfProjOrInjComplexInDegree(C, i, test);
        Append(list, [templist]);
    od;

    return list;
end);

######################################################
##
#O CompareWithIndecProjective( <M> )
##
## <M> is a indecomposable projective module.
##
## Returns the number of the vertex <M> belongs to.
##
InstallMethod( CompareWithIndecProjective,
               [ IsPathAlgebraMatModule ],
               function( M )

    local A, projlist, i;

    #first find the algebra of which M is a module
    A := RightActingAlgebra(M);
    
    projlist := IndecProjectiveModules(A);

    for i in [1..Length(projlist)] do
        if(DimensionVector(M) = DimensionVector(projlist[i])) then
            return i;
        fi;
    od;

    return false;
end);

######################################################
##
#O CompareWithIndecInjective( <M> )
##
## <M> is a indecomposable injective module.
##
## Returns the number of the vertex <M> belongs to.
##
InstallMethod( CompareWithIndecInjective,
               [ IsPathAlgebraMatModule ],
               function( M )

    local A, injlist, i;

    #first find the algebra of which M is a module
    A := RightActingAlgebra(M);
    
    injlist := IndecInjectiveModules(A);

    for i in [1..Length(injlist)] do
        if(DimensionVector(M) = DimensionVector(injlist[i])) then
            return i;
        fi;
    od;

    return false;

end);

# other methods used, not directly connected to the topic
# but not good for general use

# moved to modulehom.g?
######################################################
##
#O MultiplyListsOfMaps( <projections>, <matrix>, <inclusions> )
##
## <projections> is a list of m maps, <matrix> is a list of
## m lists of n maps, <inclusions> is a list of n maps. 
## Considering <projections> as a 1xm-matrix, <matrix> as an
## mxn-matrix and and <inclusions> as an nx1-matrix, the 
## matrix product is computed. Naturally, the maps in the 
## matrices must be composable.
##
## Output is a map (not a 1x1-matrix).
##
## Utility method not supposed to be here, but there seemed
## to be no existing GAP method for this.
##
#InstallMethod( MultiplyListsOfMaps,
#                    [ IsList, IsList, IsList ],
#                    function( projections, matrix, inclusions )
#    local sum, list, n, m, i, j;
#
#    n := Length(projections);
#    m := Length(inclusions);
#    list := [1..m];
#    sum := 0;
#
#    for i in [1..m] do
#       list[i] := projections*matrix[i];
#    od;
#
#    sum := list*inclusions;
#    return sum;
#end);

######################################################
##
#O FindAllMapComponents( <inclusions>, <d>, <projections> )
##
## <inclusions> is the direct sum inclusions of a module
## (+Pi), <projections> is the direct sum projections of a
## module (+Qj). <d> is a map (+Pi) ---> (+Qj).
##
## Output is all components of d (all elements in the matrix
## d consists of). Order on resulting set of maps:
##
## [ [P1 --> Q1, P1 --> Q2, ..., P1 --> Qm],
## [P2 --> Q1, ...,P2 --> Qm], ..., [Pn --> Q1, ...,Pn --> Qm] ]
##
InstallMethod( FindAllMapComponents,
               [ IsPathAlgebraMatModuleHomomorphism ],
               function( d )

    local firstPart, n, fullPart, i, inclusions, projections;
    inclusions := DirectSumInclusions(Source(d));
    projections := DirectSumProjections(Range(d));
    firstPart := inclusions*d;
    n := Length(inclusions);
    fullPart := [];
    for i in [1..n] do
        Append(fullPart, [List(projections, f -> firstPart[i]*f)]);
    od;
    return fullPart;

end);

######################################################
##
#O MatricesOfDualMap( <f> )
##
## <f> is a map Pj* --> Pi*, where Pj* and Pi* are
## projective A^{op}-modules.
##
## Returns the matrices for the dual map Df: Ii --> Ij
## (not the map itself!) in A.
##
InstallMethod( MatricesOfDualMap,
               [ IsPathAlgebraMatModuleHomomorphism ],
               function( f )
    local maps, origmaps, x;
    origmaps := f!.maps;
    maps := List(origmaps, x -> TransposedMat(x));
    return maps;
end);



######################################################
##
#O DirectSumMinusSummands( <M>, <list> )
##
## <M> = (+Mi) a direct sum M_1 + M_2 + ... + M_n, not nec. non-iso
## <list> is a list of integers [a1, a2, ... , ar], r <= n
##
## Output:
## let M' = M/(+Mj) where Mj is M_a1 + M_a2 + ... + M_ar
## the output is [inclusion, projection] where inclusion is the canonical
## inclusion M' --> M and projection is the canonical projection M --> M'
## if M' = 0, then a pair of appropriate zero maps is returned.
##
InstallMethod( DirectSumMinusSummands,
               [ IsPathAlgebraMatModule, IsList ],
               function( M, list )
    local inclusions, tempSource, redundantSummand, test, i, complementList,
          newModule, sourceProjections, j, inclusion, projections, rangeInclusions, projection;
    
    inclusions := DirectSumInclusions(M);
    tempSource := List(inclusions, map -> Source(map));

    # keep the redundant summands
    redundantSummand := tempSource{list};

    # remove them from tempSource
    test := function(i)
    	 if i in list then
	    return 0;
	 fi;
	 return i;
    end;
    complementList := PositionsNonzero(List([1..Length(tempSource)],
    		      l -> test(l)));
    tempSource := tempSource{complementList};

    if IsEmpty(tempSource) then
        return[ZeroMapping(ZeroModule(RightActingAlgebra(M)),M),
               ZeroMapping(M, ZeroModule(RightActingAlgebra(M)))];
    fi;
    newModule := DirectSumOfQPAModules(tempSource);

    # projections from new module to all dir. summands of the original module
    sourceProjections := DirectSumProjections(newModule);
    j := 1;
    for i in list do
        sourceProjections := Concatenation(sourceProjections{[1..(i-1)]},
                         [ZeroMapping(newModule,redundantSummand[j])],
                         sourceProjections{[(i)..Length(sourceProjections)]});
        j := j + 1;
    od;
    inclusion := sourceProjections*inclusions;

    # doing the same thing for projections
    projections := DirectSumProjections(M);

    rangeInclusions := DirectSumInclusions(newModule);
    j := 1;
    for i in list do
        rangeInclusions := Concatenation(rangeInclusions{[1..(i-1)]},
                         [ZeroMapping(redundantSummand[j],newModule)],
                         rangeInclusions{[(i)..Length(rangeInclusions)]});
        j := j + 1;
    od;
    projection := projections*rangeInclusions;
    return [inclusion,projection];
end);