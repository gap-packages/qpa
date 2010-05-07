# GAP Implementation
# This file was generated from 
# $Id: groebner.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallMethod( GroebnerBasis, 
  "for an ideal of a path algebra and a collection of elements",
  IsIdenticalObj,
  [ IsFLMLOR, IsElementOfMagmaRingModuloRelationsCollection ], 0,
  function(I, rels)
    local GB, R;

    if not IsElementOfPathRing(ElementsFamily(FamilyObj(I))) then
        TryNextMethod();
    fi;

    R := LeftActingRingOfIdeal(I);

    GB := Objectify( NewType( FamilyObj(rels),
                       IsGroebnerBasis and IsGroebnerBasisDefaultRep ),
                       rec( ideal := I,
                            relations := AsSortedList(rels),
                            algebra := R
                          )
	            );

    SetIsCompleteGroebnerBasis(GB, true);
    SetGroebnerBasisOfIdeal(I, GB);

    return GB;

  end
);


InstallMethod( CompletelyReduce,
  "for complete tip reduced groebner bases",
  true, # FIXME
  [ IsGroebnerBasis and IsCompleteGroebnerBasis 
    and IsTipReducedGroebnerBasis, IsRingElement], 0,
  function( GB, rel )
    local pats, st, tip, redRel, zero;

    zero := Zero(rel);
    redRel := zero;

    while rel <> zero do
      rel := TipReduce(GB, rel);
      tip := Tip(rel);
      rel := rel - tip;
      redRel := redRel + tip;
    od;

    return redRel;

  end
);


InstallMethod( CompletelyReduceGroebnerBasis,
  "for complete groebner bases",
  true,
  [IsGroebnerBasis and IsCompleteGroebnerBasis], 0,
  function( GB )
    local i, n, tip, rel;

    # first tip reduce it
    GB := TipReduceGroebnerBasis(GB);

    # now completely reduce each relation
    n := Length(GB!.relations);

    for i in [1..n] do
      rel := GB!.relations[i];
      tip := Tip(rel);
      GB!.relations[i] := tip + CompletelyReduce(GB, rel - tip);
    od;

    SetIsCompletelyReducedGroebnerBasis(GB, true);

    return GB;

  end
);


ReduceRelation := function( rel, redRel, extTip, start, rend )
  local redCoeff, fam, n, left, right, q, qFam, tipCoeff;
 
  fam := FamilyObj(rel);
  q := QuiverOfPathAlgebra(fam!.pathRing);
  qFam := ElementsFamily(FamilyObj(q));

  tipCoeff := TipCoefficient(rel);
  redCoeff := tipCoeff / TipCoefficient(redRel) ;
  if start <> 1 then
    left := extTip{[1..start-1]};
    left := ElementOfMagmaRing(fam, fam!.zeroRing, 
                               [ One(fam!.zeroRing)],
                               [ ObjByExtRep(qFam, left) ]);
  else
    left := One(rel);
  fi;

  if  rend < Length(extTip) then
    right := extTip{[rend+1..Length(extTip)]};
    right := ElementOfMagmaRing( fam, fam!.zeroRing,
                                [ One(fam!.zeroRing)],
                                [ ObjByExtRep(qFam, right) ]);
  else
    right := One(rel);
  fi;

  Info(InfoGroebnerBasis, 1, "redCoeff: ", redCoeff);
  Info(InfoGroebnerBasis, 1, "left: ", left);
  Info(InfoGroebnerBasis, 1, "right: ", right);
  Info(InfoGroebnerBasis, 1, "rel: ", rel);
  Info(InfoGroebnerBasis, 1, "redRel: ", redRel);

  rel := rel - redCoeff * left * redRel * right;

  return rel;

end;


InstallMethod( TipReduce,
  "for complete tip reduced groebner bases",
  true, # FIXME
  [IsGroebnerBasis and IsCompleteGroebnerBasis 
   and IsTipReducedGroebnerBasis, IsRingElement], 0,
  function( GB, rel )
    local matches, st, tip, tipMon, extTip, redMatch, redPat, 
          zero, start, rend, patLen;

    st := GB!.staticDict;
    zero := Zero(rel);

    repeat

      tip := Tip(rel);
      tipMon := TipMonomial(rel);
      extTip := ExtRepOfObj(tipMon);

      matches := Search(st, tipMon);
      if not IsEmpty(matches) then
        redMatch := matches[1];
        redPat := Patterns(st)[redMatch[2][1]];
        patLen := Length(WalkOfPath(redPat));
        start := redMatch[1] - patLen + 1;
        rend := redMatch[1];
                
        rel := ReduceRelation(rel, GB!.relations[redMatch[2][1]],
                                  extTip, start, rend);
      fi;

    until IsEmpty(matches) or rel = zero;

    return rel;

  end
);


InstallMethod( TipReduce,
  "for right groebner bases and an element",
  true, # FIXME
  [ IsRightGroebnerBasis and IsGroebnerBasisDefaultRep
    and IsTipReducedGroebnerBasis,
    IsRingElement ], 0,
  function( gb, rel )
    local zero, dict, relations, tip, extTip, matches,
          redMatch, start, rend;

    zero := Zero(rel);
    dict := gb!.staticDict;
    relations := gb!.relations;

    repeat

      tip := TipMonomial(rel);
      extTip := ExtRepOfObj(tip);
      matches := PrefixSearch(dict, tip);
      if not IsEmpty(matches) then
        redMatch := matches[1];
        start := 1;
        rend := redMatch[1];
        rel := ReduceRelation(rel, relations[redMatch[2][1]],
                              extTip, start, rend);
      fi;

    until IsEmpty(matches) or rel = zero;

    return rel;

  end
);


InstallMethod( TipReduceGroebnerBasis, 
  "for groebner bases in default representation",
  true,
  [IsGroebnerBasis 
   and IsCompleteGroebnerBasis 
   and IsGroebnerBasisDefaultRep], 0,
  function( GB )
    local gbRels, trgbRels, rel, tip, extTip, pats,
          st, redPat, n, redCoeff, extras, extra, t, quiver;

    if not (HasIsTipReducedGroebnerBasis(GB)
            and IsTipReducedGroebnerBasis(GB))
    then
      gbRels := Set(GB!.relations);
      trgbRels := [];
      st := CreateSuffixTree();
 
      while not IsEmpty(gbRels) do

        rel := gbRels[1];
        RemoveSet(gbRels, rel);
        Info(InfoGroebnerBasis, 1, "Removed relation: ", rel);

        if rel <> Zero(rel) then

          tip := TipMonomial(rel);
          extTip := ExtRepOfObj(tip);
          pats := AllPatternsInSequence(st, extTip);

          if not IsEmpty(pats) then
            redPat := pats[1];
            rel := ReduceRelation(rel, trgbRels[redPat[1]], extTip, redPat[2],
                                  redPat[2] + Length(st[redPat[1]])-1);
 
            if rel <> Zero(rel) then
              Info(InfoGroebnerBasis,1, "Adding relation: ", rel);
              AddSet(gbRels, rel);
            fi;

          else
            extras := SequenceInPatterns(st, extTip);

            if not IsEmpty(extras) then
              for extra in extras do
                if IsBound(trgbRels[extra[1]]) then
                  AddSet(gbRels, trgbRels[extra[1]]);
                  Unbind(trgbRels[extra[1]]);
                  DeletePatternFromSuffixTree(st, extra[1]);
                fi;
              od;
            fi;

            t := InsertPatternIntoSuffixTree(st, extTip);
            trgbRels[t] := rel;

          fi;
        fi;
      od;

      Sort(trgbRels);
      GB!.relations := trgbRels;

    fi;

    if not IsBound(GB!.staticDict) then
        quiver := 
          QuiverOfPathAlgebra(FamilyObj(GB!.relations[1])!.pathRing);

        # Construct a static dictionary for later reductions
        GB!.staticDict := QuiverStaticDictionary(quiver,
                              List(GB!.relations, TipMonomial));
    fi;

    SetIsTipReducedGroebnerBasis(GB, true);

    return GB;

  end
);


InstallMethod( TipReduceGroebnerBasis,
  "for right groebner bases",
  true,
  [IsRightGroebnerBasis and IsGroebnerBasisDefaultRep], 0,
  function( gb )
    local gbRels, trgbRels, quiver, i, relTip, relWalk, iWalk,
            reduced, trel, rel;


    if not (HasIsTipReducedGroebnerBasis(gb)
            and IsTipReducedGroebnerBasis(gb))
    then
      gbRels := Set(gb!.relations);
      trgbRels := [];

      while not IsEmpty(gbRels) do
        rel := gbRels[1];
        RemoveSet(gbRels, rel);

        if rel <> Zero(rel) then
          reduced := false;
          relWalk := WalkOfPath(TipMonomial(rel));

          for i in [1..Length(trgbRels)] do
            iWalk := WalkOfPath(TipMonomial(trgbRels[i]));

            if PositionSublist(iWalk, relWalk) = 1 then
              rel := ReduceRelation(rel, trgbRels[i],
                         ExtRepOfObj(TipMonomial(rel)),
                         1, Length(trgbRels[i]));
              reduced := true;
              break;
            fi;
          od;

          if not reduced then
            for i in [1..Length(trgbRels)] do
              trel := trgbRels[i];
              iWalk := WalkOfPath(TipMonomial(trel));
              if PositionSublist(relWalk, iWalk) = 1 then
                trel := ReduceRelation(trel, rel,
                            ExtRepOfObj(TipMonomial(trel)),
                            1, Length(rel));
                AddSet(gbRels, trel);
                Unbind(trgbRels[i]);
              fi;
            od;

            trgbRels := Set(trgbRels);
            AddSet(trgbRels, rel);

          fi;
        fi;
      od;

      Sort(trgbRels);
      gb!.relations := trgbRels;
      SetIsTipReducedGroebnerBasis(gb, true);        

    fi;

    if not IsBound(gb!.staticDict) then
      quiver := 
        QuiverOfPathAlgebra(FamilyObj(gb!.relations[1])!.pathRing);

        # Construct a static dictionary for later reductions
      gb!.staticDict := QuiverStaticDictionary(quiver,
                          List(gb!.relations, TipMonomial));
    fi;

    return gb;

  end
);


InstallTrueMethod(IsTipReducedGroebnerBasis, 
  IsGroebnerBasis and IsCompletelyReducedGroebnerBasis);


InstallMethod( Iterator,
  "for groebner bases",
  true,
  [IsGroebnerBasis and IsGroebnerBasisDefaultRep], 0,
  function(GB)
    return Objectify( NewType( IteratorsFamily,
        IsIterator and IsMutable and IsGroebnerBasisIteratorRep ),
        rec( relations := GB!.relations, position := 1 ) );
  end
);


InstallMethod( IsDoneIterator,
  "for iterators of groebner bases",
  true,
  [IsIterator and IsGroebnerBasisIteratorRep], 0,
  function(iterator)
    return iterator!.position > Length(iterator!.relations);
  end
);


InstallMethod( NextIterator,
  "for iterators of groebner bases",
  true,
  [IsIterator and IsMutable and IsGroebnerBasisIteratorRep], 0,
  function(iterator)
    local elem;
    elem := iterator!.relations[iterator!.position];
    iterator!.position := iterator!.position + 1;
    return elem;
  end
);


InstallMethod( Enumerator,
  "for groebner bases",
  true,
  [IsGroebnerBasis and IsGroebnerBasisDefaultRep], 0,
  function(GB)
    return Immutable(GB!.relations);
  end
);


InstallOtherMethod( SetEnumerator,
  "for groebner bases",
  true,
  [ IsGroebnerBasisDefaultRep and IsAttributeStoringRep, IsList ],
  SUM_FLAGS+1,
  function( GB, relations )
    if HasIsCompletelyReducedGroebnerBasis( GB ) and
       IsCompletelyReducedGroebnerBasis( GB )
    then
        GB!.Enumerator := relations;
        SetFilterObj( GB, HasEnumerator );
    fi;
  end
);


InstallOtherMethod( SetAsList,
  "for groebner bases",
  true,
  [ IsGroebnerBasisDefaultRep, IsList ],
  SUM_FLAGS+1,
  function( GB, list )
    if HasIsCompletelyReducedGroebnerBasis( GB ) and
       IsCompletelyReducedGroebnerBasis( GB )
    then
        GB!.AsList := list;
        SetFilterObj( GB, HasAsList );
    fi;
  end
);


InstallMethod( Nontips,
  "for complete Groebner bases",
  true,
  [IsCompleteGroebnerBasis and IsGroebnerBasisDefaultRep], 0,
  function( GB )
    TipReduceGroebnerBasis( GB ); # only does it if not reduced
    return DifferenceWords(GB!.staticDict);
  end
);


InstallMethod( AdmitsFinitelyManyNontips,
  "for complete Groebner bases",
  true,
  [IsCompleteGroebnerBasis and IsGroebnerBasisDefaultRep], 0,
  function( GB )
    TipReduceGroebnerBasis(GB);
    return IsFiniteDifference(GB!.staticDict);
  end
);


InstallMethod( NontipSize,
  "for complete Groebner bases",
  true,
  [ IsCompleteGroebnerBasis and IsGroebnerBasisDefaultRep ], 0,
  function(GB)
    TipReduceGroebnerBasis(GB);
    return DifferenceSize(GB!.staticDict);
  end
);


InstallMethod( ViewObj,
  "for groebner bases",
  true,
  [IsGroebnerBasis and IsGroebnerBasisDefaultRep], 0,
  function(GB)
    Print("<");

    if not HasIsCompleteGroebnerBasis(GB) 
       or not IsCompleteGroebnerBasis(GB) then
      Print("partial ");
    else
      Print("complete ");
    fi;

    if IsRightGroebnerBasis(GB) then
      Print("right ");
    else
      Print("two-sided ");
    fi;

    Print("Groebner basis containing ", 
          Length(GB!.relations), " elements>");

  end
);


InstallMethod( PrintObj,
  "for groebner bases",
  true,
  [IsGroebnerBasis and IsGroebnerBasisDefaultRep], 0,
  function(GB)
    if IsRightGroebnerBasis(GB) then
      Print("Right");
    fi;

    Print("GroebnerBasis( ", GB!.algebra, ", ", GB!.relations, " )");

  end
);


InstallMethod( IsPrefixOfTipInTipIdeal,
  "for a groebner basis and path algebra element",
  IsCollsElms,
  [ IsCompleteGroebnerBasis, IsElementOfMagmaRingModuloRelations ], 0,
  function( GB, R )
    local dict, matches, tip, extTip;

    tip := TipMonomial(R);

    if IsZeroPath(tip) then
      return true;
    else
      GB := TipReduceGroebnerBasis(GB); # create dictionary
      dict := GB!.staticDict;
      matches := Search(dict, tip);

      if not IsEmpty(matches) 
         and matches[1][1] < Length(WalkOfPath(tip))
      then
          return true;
      else
          return false;
      fi;

    fi;

  end
);


InstallOtherMethod( RightGroebnerBasis,
  "for two-sided ideals with a complete Groebner basis",
  true,
  [IsRing and HasGroebnerBasisOfIdeal], 0,
  function( I )
    local gb, relations, p, g, x, rightgb, ntips, R;

    gb := GroebnerBasisOfIdeal( I );
    R := LeftActingRingOfIdeal(I);

    if not IsCompleteGroebnerBasis( gb ) then
      TryNextMethod();
    fi;

    CompletelyReduceGroebnerBasis(gb);

    relations := [];
    ntips := Nontips(gb);

    # we have finite number of nontips in first case:
    if (ntips <> fail) then
      for p in ntips do
        for g in AsList(gb) do
          x := p * g;
          if not IsPrefixOfTipInTipIdeal( gb, x ) then
            Add(relations, x);
          fi;
        od;
      od;
    else
      Error("The quotient of R by the given ideal I must be finite dimensional.\n");
    fi;

    Sort(relations);

    rightgb := Objectify( NewType(FamilyObj(relations), 
                                  IsRightGroebnerBasis
                                  and IsGroebnerBasisDefaultRep),
                          rec( relations := relations ) );

    SetIsCompleteGroebnerBasis(rightgb, true);
    SetIsTipReducedGroebnerBasis(rightgb, true);

    return rightgb;

  end
);
