# GAP Implementation
# This file was generated from 
# $Id: ordering.gi,v 1.2 2012/02/27 12:26:34 sunnyquiver Exp $

InstallMethod(LessThanByOrdering,
 "for orderings and two paths",
 true,
 [IsQuiverOrdering, IsPath, IsPath], 0,
 function(O,a,b)
     return ComparisonFunction(O)(a,b) < 0;
 end
);

InstallMethod(IsWellOrdering,
 "for an ordering",
 true,
 [IsQuiverOrdering],0,
 function(O)
     return IsWellSubset(O,GeneratorsOfQuiver(O!.quiver));
  end
);


InstallMethod(IsWellReversedOrdering,
 "for an ordering",
 true,
 [IsQuiverOrdering],0,
 function(O)
     return IsWellReversedSubset(O,GeneratorsOfQuiver(O!.quiver));
  end
);


InstallMethod(IsTotalOrdering,
 "for an ordering",
 true,
 [IsQuiverOrdering],0,
 function(O)
     return IsTotalSubset(O,GeneratorsOfQuiver(O!.quiver));
  end
);


InstallMethod(IsAdmissibleOrdering,
 "for an ordering",
 true,
 [IsQuiverOrdering],0,
 function(O)
     #Check IsTotalOrdering() first, because it is often faster
     return IsTotalOrdering(O) and IsWellOrdering(O);
  end
);


DeclareGlobalFunction("CheckQuiverSubsetForCycles");

InstallGlobalFunction(CheckQuiverSubsetForCycles,
 function(Q,L)
     local h,n,a,b,i,r,visit;

     #Possibly a big shortcut
     if HasIsAcyclicQuiver(Q) and IsAcyclicQuiver(Q) then
         return true;
     fi;

     #Identify every arrow in the list, as well as every vertex we
     #can start from
     h:=[];
     for a in L do
         if IsArrow(a) then
             i:=ExtRepOfObj(a)[1];
             h[i]:=a;
             b:=SourceOfPath(a);
             h[ExtRepOfObj(b)[1]]:=b;
         fi;
     od;

     #If we have some arrows
     if Length(h)>0 then
         n:=[];
         i:=1;
         #Visit a vertex
         #v is the vertex
         #j is 3 * the generator position of the vertex
         #c is the label of the source vertex for the last arrow which was in the subset
         #n[j]   = the depth-first search label of the vertex
         #n[j-1] = true: the vertex was visited successfully; false: vertex is being visited
         #n[j-2] = index of the highest vertex reachable from cycles containing this vertex
         visit:=function(v,j,c)
             local adj,k,t,p,m;
             n[j]:=i;      #Label the vertex
             n[j-1]:=true; #Mark us in the process of being visited
             i:=i+1;
             #Check all the arrows 
             adj:=OutgoingArrowsOfVertex(v);
             for p in adj do
                 t:=TargetOfPath(p);
                 k:=ExtRepOfObj(t)[1]*3;
                 #Arrow is inside our subset
                 if IsBound(h[ExtRepOfObj(p)[1]]) then
                     #Arrow points to an unvisited vertex
                     if not IsBound(n[k]) then
                         if not visit(t,k,n[j]) then
                             return false;
                         fi;
                     #The vertex can reach us
                     #Or the vertex is part of a cycle that reaches us
                     elif n[k-1] or (IsBound(n[k-2]) and n[n[k-2]-1]) then
                         return false;
                     fi;
                 #Arrow is outside our subset
                 else
                     #Arrow points to an unvisited vertex
                     if not IsBound(n[k]) then
                         if not visit(t,k,c) then
                             return false;
                         #If it's part of a cycle that reaches us, mark us
                         elif IsBound(n[k-2]) and n[n[k-2]-1] and
                          ((not IsBound(n[j-2])) or n[n[j-2]]>n[n[k-2]]) then
                             n[j-2]:=n[k-2];
                         fi;
                     #The vertex can reach us
                     elif n[k-1] or (IsBound(n[k-2]) and n[n[k-2]-1]) then
                         if IsBound(n[k-2]) then
                             m:=n[k-2];
                         else
                             m:=k;
                         fi;
                         #The cycle contains something from our set
                         if n[m]<=c then
                             return false;
                         #Mark us as belonging to the cycle
                         elif (not IsBound(n[j-2])) or n[n[j-2]]>n[m] then
                             n[j-2]:=m;
                         fi;
                     fi;
                 fi;
             od;
             n[j-1]:=false;
             return true;
         end;
         for b in h do
             if not IsQuiverVertex(b) then
                 break;
             fi;
             r:=ExtRepOfObj(b)[1]*3;
             if not IsBound(n[r]) then
                 if not visit(b,r,0) then
                     return false;
                 fi;
             fi;
         od;
     fi;
     return true;
  end
);


InstallGlobalFunction(LengthOrdering,function(arg)
    local O,fam;

    if not ((Length(arg) = 1 or (Length(arg) = 2 and IsQuiverOrdering(arg[2])))
        and IsQuiver(arg[1])) then
        Error("usage: LengthOrdering(Q) or LengthOrdering(Q,O)",
	      " where Q is a quiver and O is an ordering");
    fi;

    fam := NewFamily("FamilyLengthOrdering", IsLengthOrdering);
    O := Objectify(NewType(fam,IsLengthOrdering and IsAttributeStoringRep),
                   rec() );

    SetOrderingName(O, "length");

    # Set quiver in record:
    O!.quiver:=arg[1];

    # We already know that the argument is an ordering
    if Length(arg)=2 then
        SetNextOrdering(O, arg[2]);
        if HasLexicographicTable(arg[2]) then
            SetLexicographicTable(O,LexicographicTable(arg[2]));
            O!.lex:=arg[2]!.lex;
        fi;
    fi;

    SetComparisonFunction(O,function(a,b)
      local lena, lenb;

      if HasWalkOfPath(a) then
          lena := Length(WalkOfPath(a));
      else
          lena := 0;
      fi;

      if HasWalkOfPath(b) then
          lenb := Length(WalkOfPath(b));
      else
          lenb := 0;
      fi;
#Print("lena: ",lena,", lenb: ",lenb,"\n");
      if lena < lenb then
          return -1;
      elif lena = lenb then
          if HasNextOrdering(O) then
              return ComparisonFunction(NextOrdering(O))(a,b);
          else
              return 0;
          fi;
      else
          return 1;
      fi;

      end
    );

    return O;
  end
);


InstallMethod(IsWellReversedOrdering,
 "for a length ordering",
 true,
 [IsLengthOrdering],0,
 function(O)
     return IsFinite(O!.quiver);
  end
);


InstallMethod(IsTotalOrdering,
 "for a length ordering",
 true,
 [IsLengthOrdering],0,
 function(O)
     if Length(VerticesOfQuiver(O!.quiver))<1 then
         return true;
     elif HasNextOrdering(O) then
         return IsTotalOrdering(NextOrdering(O));
     else
         return false;
     fi;
  end
);


InstallMethod(IsWellSubset,
 "for a length ordering and a list of generators",
 true,
 [IsLengthOrdering,IsHomogeneousList],0,
 function(O,L)
     return true;
  end
);


InstallMethod(IsWellReversedSubset,
 "for a length ordering",
 true,
 [IsLengthOrdering,IsHomogeneousList],0,
 function(O,L)
     return CheckQuiverSubsetForCycles(O!.quiver,L);
  end
);


InstallMethod(IsTotalSubset,
 "for a length ordering and a list of generators",
 true,
 [IsLengthOrdering,IsHomogeneousList],0,
 function(O,L)
     if Length(L)<1 then
         return true;
     elif HasNextOrdering(O) then
         return IsTotalSubset(NextOrdering(O),L);
     else
         return false;
     fi;
  end
);


InstallMethod(IsWellOrdering,
 "for a lexicographic ordering and a quiver",
 true,
 [IsLexicographicOrdering],0,
 function(O)
     local i,l,a,v,m,ret;
     if IsFinite(O!.quiver) then
         return true;
     else
         l:=LexicographicTable(O);
         #This must be an arrow, or we wouldn't be infinite
         a:=l[Length(l)];
         v:=SourceOfPath(a);
         i:=ExtRepOfObj(v)[1];
         if i<>ExtRepOfObj(TargetOfPath(a))[1] or InDegreeOfVertex(v)>1 then
             return false;
         fi;
         m:=AdjacencyMatrixOfQuiver(O!.quiver);
         m[i][i]:=0;
         ret:=IsFinite(Quiver(m));
         m[i][i]:=1;
         return ret;
     fi;
  end
);


InstallMethod(IsWellOrdering,
 "for a lexicographic ordering and a quiver",
 true,
 [IsLexicographicOrdering],0,
 function(O)
     local i,l,a,v,m,ret;
     if IsFinite(O!.quiver) then
         return true;
     else
         l:=LexicographicTable(O);
         #This must be an arrow, or we wouldn't be infinite
         a:=l[Length(VerticesOfQuiver(O!.quiver))+1];
         v:=SourceOfPath(a);
         i:=ExtRepOfObj(v)[1];
         if i<>ExtRepOfObj(TargetOfPath(a))[1] or InDegreeOfVertex(v)>1 then
             return false;
         fi;
         m:=AdjacencyMatrixOfQuiver(O!.quiver);
         m[i][i]:=0;
         ret:=IsFinite(Quiver(m));
         m[i][i]:=1;
         return ret;
     fi;
  end
);


InstallMethod(IsWellSubset,
 "for a lexicographic ordering and a list of generators",
 true,
 [IsLexicographicOrdering,IsHomogeneousList],0,
 function(O,L)
     local h,n,a,b,i,r,lgst,lgst_lex,visit;
     #Possibly a big shortcut
     if HasIsAcyclicQuiver(O!.quiver) and IsAcyclicQuiver(O!.quiver) then
         return true;
     fi;
     #Identify every arrow in the list, as well as every vertex we
     #can start from, and the largest arrow
     h:=[];
     for a in L do
         if IsArrow(a) then
             i:=ExtRepOfObj(a)[1];
             h[i]:=a;
             if not IsBound(lgst) or
              lgst_lex<O!.lex[i] then
                 lgst:=a;
                 lgst_lex:=O!.lex[i];
             fi;
             b:=SourceOfPath(a);
             h[ExtRepOfObj(b)[1]]:=b;
         fi;
     od;
     #If we have some arrows
     if IsBound(lgst) then
         #The largest arrow may appear in cycles
         Unbind(h[ExtRepOfObj(lgst)[1]]);
         n:=[];
         i:=1;
         visit:=#Visit a vertex
                #v is the vertex
                #j is 3 * the generator position of the vertex
                #p is the 3 * the generator position of the vertex's parent, or -1
                #c is the label of the source vertex for the last arrow which was a parent
                #n[j]   = the depth-first search label of the vertex
                #n[j-1] = -1: being visited, -2: and has cycle, 1: done being visited, 2: and has cycle
                #n[j-2] = index of the highest vertex reachable from cycles containing this vertex
                function(v,j,c)
                    local adj,k,t,p,m;
                    n[j]:=i;    #Label the vertex
                    n[j-1]:=-1; #Mark us as being visited
                    i:=i+1;
                    #Check all the arrows 
                    adj:=OutgoingArrowsOfVertex(v);
                    for p in adj do
                        t:=TargetOfPath(p);
                        k:=ExtRepOfObj(t)[1]*3;
                        #Arrow is inside our subset, and not the largest
                        if IsBound(h[ExtRepOfObj(p)[1]]) then
                            #Arrow points to an unvisited vertex
                            if not IsBound(n[k]) then
                                if not visit(t,k,n[j]) then
                                    return false;
                                fi;
                            #Arrow points to the largest arrow and largest arrow belongs to a cycle
                            #Or the vertex can reach us
                            #Or the vertex is part of a cycle that reaches us
                            elif (n[k]=1 and n[k-1]=2) or n[k-1]<0 or
                             (IsBound(n[k-2]) and n[n[k-2]-1]<0) then
                                return false;
                            fi;
                        #Arrow is outside our subset, or the largest
                        else
                            #Arrow points to an unvisited vertex
                            if not IsBound(n[k]) then
                                if p=lgst then
                                    if not visit(t,k,-1) then
                                        return false;
                                    fi;
                                else
                                    if not visit(t,k,c) then
                                        return false;
                                    fi;
                                fi;
                                if IsBound(n[k-2]) and (n[n[k-2]-1]<0) and
                                 ((not IsBound(n[j-2])) or n[n[j-2]]>n[n[k-2]]) then
                                    n[j-2]:=n[k-2];
                                fi;
                            #The vertex can possibly reach us
                            elif n[k-1]<0 or
                             (IsBound(n[k-2]) and n[n[k-2]-1]<0) then
                                if IsBound(n[k-2]) then
                                    m:=n[k-2];
                                else
                                    m:=k;
                                fi;
                                #The cycle contains something from our set
                                if n[m]<=c then
                                    return false;
                                #The cycle is outside our set
                                elif (not IsBound(n[j-2])) or n[n[j-2]]>n[m] then
                                    n[j-2]:=m;
                                fi;
                                #The cycle descends directly from the largest arrow
                                if c<0 or p=lgst then
                                    n[m-1]:=-2;
                                fi;
                            fi;
                        fi;
                    od;
                    n[j-1]:=-n[j-1]; #Mark us as done being visited
                    return true;
                end;
         #Visit source of largest arrow first
         b:=SourceOfPath(lgst);
         r:=ExtRepOfObj(b)[1]*3;
         if not visit(b,r,0) then
             return false;
         fi;
         #Visit any unconnected components reachable from our subset
         for b in h do
             if not IsQuiverVertex(b) then
                 break;
             fi;
             r:=ExtRepOfObj(b)[1]*3;
             if not IsBound(n[r]) then
                 if not visit(b,r,0) then
                     return false;
                 fi;
             fi;
         od;
     fi;
     return true;
  end
);


InstallMethod(IsWellReversedSubset,
 "for a lexicographic ordering and a list of generators",
 true,
 [IsLexicographicOrdering,IsHomogeneousList],0,
 function(O,L)
     local h,n,a,b,i,r,lgst,lgst_lex,visit;
     #Possibly a big shortcut
     if HasIsAcyclicQuiver(O!.quiver) and IsAcyclicQuiver(O!.quiver) then
         return true;
     fi;
     #Identify every arrow in the list, as well as every vertex we
     #can start from, and the largest arrow
     h:=[];
     for a in L do
         if IsArrow(a) then
             i:=ExtRepOfObj(a)[1];
             h[i]:=a;
             if not IsBound(lgst) or
              lgst_lex>O!.lex[i] then
                 lgst:=a;
                 lgst_lex:=O!.lex[i];
             fi;
             b:=SourceOfPath(a);
             h[ExtRepOfObj(b)[1]]:=b;
         fi;
     od;
     #If we have some arrows
     if IsBound(lgst) then
         #The largest arrow may appear in cycles
         Unbind(h[ExtRepOfObj(lgst)[1]]);
         n:=[];
         i:=1;
         visit:=#Visit a vertex
                #v is the vertex
                #j is 3 * the generator position of the vertex
                #p is the 3 * the generator position of the vertex's parent, or -1
                #c is the label of the source vertex for the last arrow which was a parent
                #n[j]   = the depth-first search label of the vertex
                #n[j-1] = -1: being visited, -2: and has cycle, 1: done being visited, 2: and has cycle
                #n[j-2] = index of the highest vertex reachable from cycles containing this vertex
                function(v,j,c)
                    local adj,k,t,p,m;
                    n[j]:=i;    #Label the vertex
                    n[j-1]:=-1; #Mark us as being visited
                    i:=i+1;
                    #Check all the arrows 
                    adj:=OutgoingArrowsOfVertex(v);
                    for p in adj do
                        t:=TargetOfPath(p);
                        k:=ExtRepOfObj(t)[1]*3;
                        #Arrow is inside our subset, and not the largest
                        if IsBound(h[ExtRepOfObj(p)[1]]) then
                            #Arrow points to an unvisited vertex
                            if not IsBound(n[k]) then
                                if not visit(t,k,n[j]) then
                                    return false;
                                fi;
                            #Arrow points to the largest arrow and largest arrow belongs to a cycle
                            #Or the vertex can reach us
                            #Or the vertex is part of a cycle that reaches us
                            elif (n[k]=1 and n[k-1]=2) or n[k-1]<0 or
                             (IsBound(n[k-2]) and n[n[k-2]-1]<0) then
                                return false;
                            fi;
                        #Arrow is outside our subset, or the largest
                        else
                            #Arrow points to an unvisited vertex
                            if not IsBound(n[k]) then
                                if p=lgst then
                                    if not visit(t,k,-1) then
                                        return false;
                                    fi;
                                else
                                    if not visit(t,k,c) then
                                        return false;
                                    fi;
                                fi;
                                if IsBound(n[k-2]) and (n[n[k-2]-1]<0) and
                                 ((not IsBound(n[j-2])) or n[n[j-2]]>n[n[k-2]]) then
                                    n[j-2]:=n[k-2];
                                fi;
                            #The vertex can possibly reach us
                            elif n[k-1]<0 or
                             (IsBound(n[k-2]) and n[n[k-2]-1]<0) then
                                if IsBound(n[k-2]) then
                                    m:=n[k-2];
                                else
                                    m:=k;
                                fi;
                                #The cycle contains something from our set
                                if n[m]<=c then
                                    return false;
                                #The cycle is outside our set
                                elif (not IsBound(n[j-2])) or n[n[j-2]]>n[m] then
                                    n[j-2]:=m;
                                fi;
                                #The cycle descends directly from the largest arrow
                                if c<0 or p=lgst then
                                    n[m-1]:=-2;
                                fi;
                            fi;
                        fi;
                    od;
                    n[j-1]:=-n[j-1]; #Mark us as done being visited
                    return true;
                end;
         #Visit source of largest arrow first
         b:=SourceOfPath(lgst);
         r:=ExtRepOfObj(b)[1]*3;
         if not visit(b,r,0) then
             return false;
         fi;
         #Visit any unconnected components reachable from our subset
         for b in h do
             if not IsQuiverVertex(b) then
                 break;
             fi;
             r:=ExtRepOfObj(b)[1]*3;
             if not IsBound(n[r]) then
                 if not visit(b,r,0) then
                     return false;
                 fi;
             fi;
         od;
     fi;
     return true;
  end
);


InstallMethod(IsTotalSubset,
 "for a lexicographic ordering and a list of generators",
 true,
 [IsLexicographicOrdering,IsHomogeneousList],0,
 function(O,L)
     return true;
 end);
InstallGlobalFunction(LeftLexicographicOrdering,function(arg)
    local fam,O,Q,lex,i,j,rep,gens,oldgens;

    if not (Length(arg)=1 or (Length(arg)=2 and IsHomogeneousList(arg[2])) and IsQuiver(arg[1])) then
        Error("usage: LeftLexicographicOrdering(Q,L), Q must be q quiver, L is an optional list of gens");
    fi;

    fam:=NewFamily("FamilyLeftLexicographicOrderings",
                   IsLeftLexicographicOrdering);
    O:=Objectify(NewType(fam,IsLeftLexicographicOrdering and
                             IsAttributeStoringRep),rec());
    Q:=arg[1];
    O!.quiver:=Q;

    oldgens:=GeneratorsOfQuiver(Q);
     lex:=[];
     gens:=[];
     j:=1;
     if Length(arg)>=2 then
         for i in arg[2] do
             if IsQuiverVertex(i) then
                 rep:=ExtRepOfObj(i)[1];
                 if not IsBound(lex[rep]) then
                     lex[rep]:=j;
                     j:=j+1;
                     Add(gens,i);
                 fi;
             fi;
         od;
     fi;
     if Length(VerticesOfQuiver(Q))>=j then
         for i in oldgens do
             if not IsQuiverVertex(i) then
                 break;
             fi;
             rep:=ExtRepOfObj(i)[1];
             if not IsBound(lex[rep]) then
                 lex[rep]:=j;
                 j:=j+1;
                 Add(gens,i);
             fi;
         od;
     fi;
     if Length(arg)>=2 then
         for i in arg[2] do
             if IsArrow(i) then
                 rep:=ExtRepOfObj(i)[1];
                 if not IsBound(lex[rep]) then
                     lex[rep]:=j;
                     j:=j+1;
                     Add(gens,i);
                fi;
             fi;
         od;
     fi;
     if Length(oldgens)>=j then
         for i in [Length(VerticesOfQuiver(Q))+1..Length(oldgens)] do
             rep:=ExtRepOfObj(oldgens[i])[1];
             if not IsBound(lex[rep]) then
                 lex[rep]:=j;
                 j:=j+1;
                 Add(gens,oldgens[i]);
             fi;
         od;
     fi;

    SetOrderingName(O,"left lexicographic");
    SetLexicographicTable(O,gens);
    O!.lex:=lex;
    SetComparisonFunction(O,function(a,b)
        local m,n,l,r,s,j;
        if IsZeroPath(a) then
            if IsZeroPath(b) then
                return 0;
            else
                return -1;
            fi;
        elif IsZeroPath(b) then
            return 1;
        fi;
        r:=ExtRepOfObj(a); # remember, these are paths
        s:=ExtRepOfObj(b);
        m:=Length(r);
        n:=Length(s);
        l:=Minimum(m, n);
        for j in [1..l] do
            if lex[r[j]]<>lex[s[j]] then
                if lex[r[j]]<lex[s[j]] then
                    return -1;
                else
                    return 1;
                fi;
            fi;
        od;
        if m<n then
            return -1;
        elif m=n then
            if HasNextOrdering(O) then
                return ComparisonFunction(NextOrdering(O))(a,b);
            else
                return 0;
            fi;
        else
            return 1;
        fi;
    end);

    return O;

  end
);


InstallGlobalFunction(RightLexicographicOrdering,function(arg)
    local fam,O,Q,lex,i,j,rep,gens,oldgens;

    if not (Length(arg)=1 or (Length(arg)=2 and IsHomogeneousList(arg[2])) and IsQuiver(arg[1])) then
        Error("usage: RightLexicographicOrdering(Q,L), Q must be q quiver, L is an optional list of gens");
    fi;

    fam:=NewFamily("FamilyRightLexicographicOrderings",
                   IsRightLexicographicOrdering);
    O:=Objectify(NewType(fam,IsRightLexicographicOrdering and
                             IsAttributeStoringRep),rec());
    Q:=arg[1];
    O!.quiver:=Q;

    oldgens:=GeneratorsOfQuiver(Q);
     lex:=[];
     gens:=[];
     j:=1;
     if Length(arg)>=2 then
         for i in arg[2] do
             if IsQuiverVertex(i) then
                 rep:=ExtRepOfObj(i)[1];
                 if not IsBound(lex[rep]) then
                     lex[rep]:=j;
                     j:=j+1;
                     Add(gens,i);
                 fi;
             fi;
         od;
     fi;
     if Length(VerticesOfQuiver(Q))>=j then
         for i in oldgens do
             if not IsQuiverVertex(i) then
                 break;
             fi;
             rep:=ExtRepOfObj(i)[1];
             if not IsBound(lex[rep]) then
                 lex[rep]:=j;
                 j:=j+1;
                 Add(gens,i);
             fi;
         od;
     fi;
     if Length(arg)>=2 then
         for i in arg[2] do
             if IsArrow(i) then
                 rep:=ExtRepOfObj(i)[1];
                 if not IsBound(lex[rep]) then
                     lex[rep]:=j;
                     j:=j+1;
                     Add(gens,i);
                fi;
             fi;
         od;
     fi;
     if Length(oldgens)>=j then
         for i in [Length(VerticesOfQuiver(Q))+1..Length(oldgens)] do
             rep:=ExtRepOfObj(oldgens[i])[1];
             if not IsBound(lex[rep]) then
                 lex[rep]:=j;
                 j:=j+1;
                 Add(gens,oldgens[i]);
             fi;
         od;
     fi;

    SetOrderingName(O,"right lexicographic");
    SetLexicographicTable(O,gens);
    O!.lex:=lex;
    SetComparisonFunction(O,function(a,b)
        local m,n,l,r,s,j;

        if IsZeroPath(a) then
            if IsZeroPath(b) then
                return 0;
            else
                return -1;
            fi;
        elif IsZeroPath(b) then
            return 1;
        fi;
        r:=ExtRepOfObj(a); # remember, these are paths
        s:=ExtRepOfObj(b); 
        m:=Length(r);
        n:=Length(s);
        l:=Minimum(m, n);
        for j in [0..l-1] do
            if lex[r[m-j]]<>lex[s[n-j]] then
                if lex[r[m-j]]<lex[s[n-j]] then
                    return -1;
                else
                    return 1;
                fi;
            fi;
        od;
        if m<n then
            return -1;
        elif m=n then
            if HasNextOrdering(O) then
                return ComparisonFunction(NextOrdering(O))(a,b);
            else
                return 0;
            fi;
        else
            return 1;
        fi;
    end);
    SetIsTotal(O,true);

    return O;

  end
);


InstallGlobalFunction(ReverseOrdering,function(arg)
    local fam,O;
    if not (Length(arg)=2 and IsQuiverOrdering(arg[2]) and IsQuiver(arg[1])) then
        Error("Usage: ReverseOrdering(Q,O) where Q is a quiver and O is an ordering");
    fi;

    fam:=NewFamily("FamilyReverseOrderings",IsReverseOrdering);
    O:=Objectify(NewType(fam,IsReverseOrdering and IsAttributeStoringRep),
                 rec());

    O!.quiver:=arg[1];
    SetOrderingName(O,"reverse");
    SetNextOrdering(O,arg[2]);
    if HasLexicographicTable(arg[2]) then
        SetLexicographicTable(O,LexicographicTable(arg[2]));
        O!.lex:=arg[2]!.lex;
    fi;

    SetComparisonFunction(O,function(a,b)
        local ret;
        if IsZeroPath(a) then
            if IsZeroPath(b) then
                return 0;
            else
                return ComparisonFunction(NextOrdering(O))(a,b);
            fi;
        elif IsZeroPath(b) then
            return ComparisonFunction(NextOrdering(O))(a,b);
        fi;
        if IsQuiverVertex(a) then
            if not IsQuiverVertex(b) then
                return ComparisonFunction(NextOrdering(O))(a,b);
            fi;
        elif IsQuiverVertex(b) then
            return ComparisonFunction(NextOrdering(O))(a,b);
        fi;
        return -ComparisonFunction(NextOrdering(O))(a,b);
    end);

    return O;
  end
);


InstallMethod(IsWellOrdering,
 "for a reverse ordering and a quiver",
 true,
 [IsReverseOrdering],0,
 function(O)
     return IsWellReversedOrdering(NextOrdering(O));
  end
);


InstallMethod(IsWellReversedOrdering,
 "for a reverse ordering and a quiver",
 true,
 [IsReverseOrdering],0,
 function(O)
     return IsWellOrdering(NextOrdering(O));
  end
);


InstallMethod(IsTotalOrdering,
 "for a reverse ordering",
 true,
 [IsReverseOrdering],0,
 function(O)
     return IsTotalOrdering(NextOrdering(O));
  end
);


InstallMethod(IsWellSubset,
 "for a reverse ordering and a list of generators",
 true,
 [IsReverseOrdering,IsHomogeneousList],0,
 function(O,L)
     return IsWellReversedSubset(NextOrdering(O),L);
  end
);


InstallMethod(IsWellReversedSubset,
 "for a reverse ordering and a list of generators",
 true,
 [IsReverseOrdering,IsHomogeneousList],0,
 function(O,L)
     return IsWellSubset(NextOrdering(O),L);
  end
);


InstallMethod(IsTotalSubset,
 "for a reverse ordering and a list of generators",
 true,
 [IsReverseOrdering,IsHomogeneousList],0,
 function(O,L)
     return IsTotalSubset(NextOrdering(O),L);
  end
);


InstallMethod(IsWellReversedOrdering,
 "for a vector ordering",
 true,
 [IsVectorOrdering],0,
 function(O)
     return IsFinite(O!.quiver);
  end
);


InstallMethod(IsTotalOrdering,
 "for a vector ordering",
 true,
 [IsVectorOrdering],0,
 function(O)
     if Length(VerticesOfQuiver(O!.quiver))<1 then
         return true;
     elif HasNextOrdering(O) then
         return IsTotalOrdering(NextOrdering(O));
     else
         return false;
     fi;
  end
);


InstallMethod(IsWellSubset,
 "for a vector ordering and a list of generators",
 true,
 [IsVectorOrdering,IsHomogeneousList],0,
 function(O,L)
     return true;
  end
);


InstallMethod(IsWellReversedSubset,
 "for a vector ordering and a list of generators",
 true,
 [IsVectorOrdering,IsHomogeneousList],0,
 function(O,L)
     return CheckQuiverSubsetForCycles(O!.quiver,L);
  end
);


InstallMethod(IsTotalSubset,
 "for a vector ordering and a list of generators",
 true,
 [IsVectorOrdering,IsHomogeneousList],0,
 function(O,L)
     if Length(L)<1 then
         return true;
     elif HasNextOrdering(O) then
         return IsTotalSubset(NextOrdering(O),L);
     else
         return false;
     fi;
  end
);


InstallGlobalFunction(LeftVectorOrdering,function(arg)
    local fam,O,Q,lex,i,j,rep,gens,oldgens;
    if not (Length(arg)=1 or (Length(arg)=2 or (Length(arg)=3 and
     IsQuiverOrdering(arg[3])) and IsHomogeneousList(arg[2])) and IsQuiver(arg[1]))
     then
        Error("usage: LeftVectorOrdering(Q,L,O), the first arg must be a quiver, the second a list of generators, and the optional third argument an ordering");
    fi;
    fam:=NewFamily("FamilyLeftVectorOrderings",IsLeftVectorOrdering);
    O:=Objectify(NewType(fam,IsLeftVectorOrdering and IsAttributeStoringRep),
                 rec());
    Q:=arg[1];
    O!.quiver:=Q;

    if Length(arg)>=3 and HasLexicographicTable(arg[3]) then
        oldgens:=LexicographicTable(arg[3]);
    else
        oldgens:=GeneratorsOfQuiver(Q);
    fi;
     lex:=[];
     gens:=[];
     j:=1;
     if Length(arg)>=2 then
         for i in arg[2] do
             if IsQuiverVertex(i) then
                 rep:=ExtRepOfObj(i)[1];
                 if not IsBound(lex[rep]) then
                     lex[rep]:=j;
                     j:=j+1;
                     Add(gens,i);
                 fi;
             fi;
         od;
     fi;
     if Length(VerticesOfQuiver(Q))>=j then
         for i in oldgens do
             if not IsQuiverVertex(i) then
                 break;
             fi;
             rep:=ExtRepOfObj(i)[1];
             if not IsBound(lex[rep]) then
                 lex[rep]:=j;
                 j:=j+1;
                 Add(gens,i);
             fi;
         od;
     fi;
     if Length(arg)>=2 then
         for i in arg[2] do
             if IsArrow(i) then
                 rep:=ExtRepOfObj(i)[1];
                 if not IsBound(lex[rep]) then
                     lex[rep]:=j;
                     j:=j+1;
                     Add(gens,i);
                fi;
             fi;
         od;
     fi;
     if Length(oldgens)>=j then
         for i in [Length(VerticesOfQuiver(Q))+1..Length(oldgens)] do
             rep:=ExtRepOfObj(oldgens[i])[1];
             if not IsBound(lex[rep]) then
                 lex[rep]:=j;
                 j:=j+1;
                 Add(gens,oldgens[i]);
             fi;
         od;
     fi;

    if Length(arg)>=3 then
        SetNextOrdering(O,arg[3]);
        if HasLexicographicTable(arg[3]) then
            if LexicographicTable(arg[3])<>gens then
                Error("Must use same lexicographic table");
            fi;
            O!.lex:=arg[3]!.lex;
            SetLexicographicTable(O,LexicographicTable(arg[3]));
        else
            O!.lex:=lex;
            SetLexicographicTable(O,gens);
        fi;
    else
        O!.lex:=lex;
        SetLexicographicTable(O,gens);
    fi;
    SetOrderingName(O,"left vector");
    SetComparisonFunction(O,function(a,b)
        local v,w,r,s,i,j,o;
        if IsQuiverVertex(a) or IsZeroPath(a) then
            if IsQuiverVertex(b) or IsZeroPath(b) then
                if HasNextOrdering(O) then
                    return ComparisonFunction(NextOrdering(O))(a,b);
                else
                    return 0;
                fi;
            else
                return -1;
            fi;
        elif IsQuiverVertex(b) or IsZeroPath(b) then
            return 1;
        fi;
        # At this point, both items are paths
        o:=Length(VerticesOfQuiver(O!.quiver));
        v:=ListWithIdenticalEntries(Length(gens)-o,0);
        r:=ExtRepOfObj(a);
        for i in [1..Length(r)] do
            j:=lex[r[i]]-o;
            v[j]:=v[j]+1;
        od;
        w:=ListWithIdenticalEntries(Length(gens)-o,0);
        s:=ExtRepOfObj(b);
        for i in [1..Length(s)] do
            j:=lex[s[i]]-o;
            w[j]:=w[j]+1;
        od;
        for i in [1..Length(v)] do
            if v[i]<w[i] then
                return -1;
            elif v[i]>w[i] then
                return 1;
            fi;
        od;
        if HasNextOrdering(O) then
            return ComparisonFunction(NextOrdering(O))(a,b);
        else
            return 0;
        fi;
    end);

    return O;

  end
);


InstallGlobalFunction(RightVectorOrdering,function(arg)
    local fam,O,Q,lex,i,j,rep,gens,oldgens;
    if not (Length(arg)=1 or (Length(arg)=2 or (Length(arg)=3 and
     IsQuiverOrdering(arg[3])) and IsHomogeneousList(arg[2])) and IsQuiver(arg[1]))
     then
        Error("usage: RightVectorOrdering(Q,L,O), the first arg must be a quiver, the second a list of generators, and the optional third argument an ordering");
    fi;
    fam:=NewFamily("FamilyRightVectorOrderings",IsRightVectorOrdering);
    O:=Objectify(NewType(fam,IsRightVectorOrdering and IsAttributeStoringRep),
                 rec());
    Q:=arg[1];
    O!.quiver:=Q;

    if Length(arg)>=3 and HasLexicographicTable(arg[3]) then
        oldgens:=LexicographicTable(arg[3]);
    else
        oldgens:=GeneratorsOfQuiver(Q);
    fi;
     lex:=[];
     gens:=[];
     j:=1;
     if Length(arg)>=2 then
         for i in arg[2] do
             if IsQuiverVertex(i) then
                 rep:=ExtRepOfObj(i)[1];
                 if not IsBound(lex[rep]) then
                     lex[rep]:=j;
                     j:=j+1;
                     Add(gens,i);
                 fi;
             fi;
         od;
     fi;
     if Length(VerticesOfQuiver(Q))>=j then
         for i in oldgens do
             if not IsQuiverVertex(i) then
                 break;
             fi;
             rep:=ExtRepOfObj(i)[1];
             if not IsBound(lex[rep]) then
                 lex[rep]:=j;
                 j:=j+1;
                 Add(gens,i);
             fi;
         od;
     fi;
     if Length(arg)>=2 then
         for i in arg[2] do
             if IsArrow(i) then
                 rep:=ExtRepOfObj(i)[1];
                 if not IsBound(lex[rep]) then
                     lex[rep]:=j;
                     j:=j+1;
                     Add(gens,i);
                fi;
             fi;
         od;
     fi;
     if Length(oldgens)>=j then
         for i in [Length(VerticesOfQuiver(Q))+1..Length(oldgens)] do
             rep:=ExtRepOfObj(oldgens[i])[1];
             if not IsBound(lex[rep]) then
                 lex[rep]:=j;
                 j:=j+1;
                 Add(gens,oldgens[i]);
             fi;
         od;
     fi;

    if Length(arg)>=3 then
        SetNextOrdering(O,arg[3]);
        if HasLexicographicTable(arg[3]) then
            if LexicographicTable(arg[3])<>gens then
                Error("Must use same lexicographic table");
            fi;
            O!.lex:=arg[3]!.lex;
            SetLexicographicTable(O,LexicographicTable(arg[3]));
        else
            O!.lex:=lex;
            SetLexicographicTable(O,gens);
        fi;
    else
        O!.lex:=lex;
        SetLexicographicTable(O,gens);
    fi;
    SetOrderingName(O,"right vector");
    SetComparisonFunction(O,function(a,b)
        local v,w,r,s,i,j,l,o;
        if IsQuiverVertex(a) or IsZeroPath(a) then
            if IsQuiverVertex(b) or IsZeroPath(b) then
                if HasNextOrdering(O) then
                    return ComparisonFunction(NextOrdering(O))(a,b);
                else
                    return 0;
                fi;
            else
                return -1;
            fi;
        elif IsQuiverVertex(b) or IsZeroPath(b) then
            return 1;
        fi;
        # At this point, both items are paths
        o:=Length(VerticesOfQuiver(O!.quiver));
        v:=ListWithIdenticalEntries(Length(gens)-o,0);
        r:=ExtRepOfObj(a);
        for i in [1..Length(r)] do
            j:=lex[r[i]]-o;
            v[j]:=v[j]+1;
        od;
        w:=ListWithIdenticalEntries(Length(gens)-o,0);
        s:=ExtRepOfObj(b);
        for i in [1..Length(s)] do
            j:=lex[s[i]]-o;
            w[j]:=w[j]+1;
        od;
        l:=Length(v);
        for i in [0..l-1] do
            if v[l-i]<w[l-i] then
                return -1;
            elif v[l-i]>w[l-i] then
                return 1;
            fi;
        od;
        if HasNextOrdering(O) then
            return ComparisonFunction(NextOrdering(O))(a,b);
        else
            return 0;
        fi;
    end);
    return O;
  end
);


InstallGlobalFunction(WeightOrdering,function(arg)
    local fam,O,Q,weight,i,j,rep,gens,p;
    if not (Length(arg)=1 or (Length(arg)=2 or (Length(arg)=3 and
     IsQuiverOrdering(arg[3])) and IsHomogeneousList(arg[2])) and IsQuiver(arg[1]))
     then
        Error("usage: WeightOrdering(Q,L,O), the first arg must be a quiver, the second a list of weights, one for each arrow in the quiver, and the optional third argument an ordering");
    fi;
    fam:=NewFamily("FamilyWeightOrderings",IsWeightOrdering);
    O:=Objectify(NewType(fam,IsWeightOrdering and IsAttributeStoringRep),
                 rec());
    Q:=arg[1];
    O!.quiver:=Q;
    weight:=ListWithIdenticalEntries(Length(ArrowsOfQuiver(Q)),0);
    if Length(arg)>=2 and Length(arg[2])>0 then
        for p in [1..Length(arg[2])] do
            weight[p]:=arg[2][p];
        od;
    fi;
    O!.weight:=weight;

    if Length(arg)>=3 then
        SetNextOrdering(O,arg[3]);
        if HasLexicographicTable(arg[3]) then
            SetLexicographicTable(O,LexicographicTable(arg[3]));
            O!.lex:=arg[3]!.lex;
        fi;
    fi;
    SetOrderingName(O,"weight");
    SetComparisonFunction(O,function(a,b)
        local v,w,r,s,i,o;
        if IsQuiverVertex(a) or IsZeroPath(a) then
            if IsQuiverVertex(b) or IsZeroPath(b) then
                if HasNextOrdering(O) then
                    return ComparisonFunction(NextOrdering(O))(a,b);
                else
                    return 0;
                fi;
            else
                return -1;
            fi;
        elif IsQuiverVertex(b) or IsZeroPath(b) then
            return 1;
        fi;
        # At this point, both items are paths
        o:=Length(VerticesOfQuiver(O!.quiver));
        r:=ExtRepOfObj(a);
        s:=ExtRepOfObj(b);
        w:=0;
        for i in r do
            w:=w+weight[i-o];
        od;
        v:=0;
        for i in s do
            v:=v+weight[i-o];
        od;
        if w=v then
            if HasNextOrdering(O) then
                return ComparisonFunction(NextOrdering(O))(a,b);
            else
                return 0;
            fi;
        else
            return w-v;
        fi;
    end);
    return O;
  end
);


InstallMethod(IsWellReversedOrdering,
 "for a weight ordering",
 true,
 [IsWeightOrdering],0,
 function(O)
     local z,o;
     if IsFinite(O!.quiver) then
         return true;
     fi;
     o:=Length(VerticesOfQuiver(O!.quiver));
     z:=Filtered(ArrowsOfQuiver(O!.quiver),function(a)
         return O!.weight[ExtRepOfObj(a)[1]-o]<>0;
     end);
     if Length(z)>0 and not CheckQuiverSubsetForCycles(O!.quiver,z) then
         return false;
     fi;
     return not HasNextOrdering(O) or
      IsWellReversedSubset(NextOrdering(O),z);
  end
);


InstallMethod(IsTotalOrdering,
 "for a weight ordering",
 true,
 [IsWeightOrdering],0,
 function(O)
     if Length(VerticesOfQuiver(O!.quiver))<1 then
         return true;
     elif HasNextOrdering(O) then
         return IsTotalOrdering(NextOrdering(O));
     else
         return true;
     fi;
  end
);


InstallMethod(IsWellSubset,
 "for a weight ordering and a list of generators",
 true,
 [IsWeightOrdering,IsHomogeneousList],0,
 function(O,L)
     return true;
  end
);


InstallMethod(IsWellReversedSubset,
 "for a weight ordering and a list of generators",
 true,
 [IsWeightOrdering,IsHomogeneousList],0,
 function(O,L)
     local z,o;
     if HasIsAcyclicQuiver(O!.quiver) and IsAcyclicQuiver(O!.quiver) then
         return true;
     fi;
     o:=Length(VerticesOfQuiver(O!.quiver));
     z:=Filtered(L,function(a)
         return IsArrow(a) and O!.weight[ExtRepOfObj(a)[1]-o]<>0;
     end);
     if Length(z)>0 and not CheckQuiverSubsetForCycles(O!.quiver,z) then
         return false;
     fi;
     z:=Filtered(L,function(a)
         return IsArrow(a) and O!.weight[ExtRepOfObj(a)[1]-o]=0;
     end);
     if Length(z)>0 and not CheckQuiverSubsetForCycles(O!.quiver,z) then
         return not HasNextOrdering(O) or
          IsWellReversedSubset(NextOrdering(O),z);
     fi;
     return true;
  end
);


InstallMethod(IsTotalSubset,
 "for a weight ordering and a list of generators",
 true,
 [IsWeightOrdering,IsHomogeneousList],0,
 function(O,L)
     if Length(L)<1 or
      (Length(L)=1 and IsArrow(L[1]) and
       O!.weight[ExtRepOfObj(L[1])[1]-
                 Length(VerticesOfQuiver(O!.quiver))]<>0) then
         return true;
     elif HasNextOrdering(O) then
         return IsTotalSubset(NextOrdering(O,L));
     else
         return true;
     fi;
  end
);


DeclareOperation("PathFromArrowList",[IsQuiver,IsHomogeneousList]);

InstallMethod(PathFromArrowList,
 "for a list of arrows",
 true,
 [IsQuiver,IsHomogeneousList],0,
 function(Q,L)
     local ret;
     if Length(L)<1 then
         return Zero(Q);
     elif Length(L)=1 then
         return L[1]; #This handles vertices, too
     else
         ret:=Objectify(NewType(FamilyObj(Q),
                                IsPath and IsAttributeStoringRep),
                        rec());
         SetSourceOfPath(ret,SourceOfPath(L[1]));
         SetTargetOfPath(ret,TargetOfPath(L[Length(L)]));
         SetLengthOfPath(ret,Length(L));
         SetWalkOfPath(ret,L);
         SetIsZeroPath(ret,false);
         return ret;
     fi;
  end
);


DeclareOperation("CompareArrowLists",[IsQuiverOrdering,IsHomogeneousList,IsHomogeneousList]);

InstallMethod(CompareArrowLists,
 "for an ordering and two lists of arrows",
 true,
 [IsQuiverOrdering,IsHomogeneousList,IsHomogeneousList],0,
 function(O,a,b)
     return ComparisonFunction(O)(PathFromArrowList(O!.quiver,a),
                                  PathFromArrowList(O!.quiver,b));
  end
);


InstallGlobalFunction(BlockOrdering,function(arg)
    local fam,O,Q,L,block,orders,i,j,rep;
    if not (Length(arg)=1 or (Length(arg)=2 or (Length(arg)=3 and
     IsQuiverOrdering(arg[3])) and IsHomogeneousList(arg[2])) and IsQuiver(arg[1]))
     then
        Error("usage: BlockOrdering(Q,L,O), the first arg must be a quiver, the second a list of [ordering,[generator list]] pairs, and the optional third argument an ordering");
    fi;
    L:=arg[2];
    if Length(L)<1 then
        Error("BlockOrdering: Subset list cannot be empty");
    fi;

    fam:=NewFamily("FamilyBlockOrderings",IsBlockOrdering);
    O:=Objectify(NewType(fam,IsBlockOrdering and IsAttributeStoringRep),
                 rec());
    Q:=arg[1];
    O!.quiver:=Q;

    block:=[];
    orders:=[];
    L:=Flat(L);
    i:=1;
    while i<Length(L) do
        if not IsQuiverOrdering(L[i]) then
            Error("Second argument must be a list of [ordering,[generator list]] pairs");
        fi;
        Add(orders,L[i]);
        if HasLexicographicTable(L[i]) then
            if HasLexicographicTable(O) then
                if LexicographicTable(O)<>LexicographicTable(L[i]) then
                    Error("Conflicting lexicographic tables");
                fi;
            else
                SetLexicographicTable(O,LexicographicTable(L[i]));
                O!.lex:=L[i]!.lex;
            fi;
        fi;
        i:=i+1;
        while i<Length(L) do
            if (not IsArrow(L[i])) and (not IsQuiverVertex(L[i])) then
                break;
            fi;
            block[ExtRepOfObj(L[i])[1]]:=Length(orders);
            i:=i+1;
        od;
    od;
    L:=Length(GeneratorsOfQuiver(Q));
    if L>0 then
        for i in [1..L] do
            if not IsBound(block[i]) then
                block[i]:=Length(orders);
            fi;
        od;
    fi;
    O!.block:=block;
    O!.orders:=orders;

    if Length(arg)>=3 then
        SetNextOrdering(O,arg[3]);
        if HasLexicographicTable(arg[3]) then
            if HasLexicographicTable(O) then
                if LexicographicTable(O)<>LexicographicTable(arg[3]) then
                    Error("Conflicting lexicographic tables");
                fi;
            else
                SetLexicographicTable(O,LexicographicTable(arg[3]));
                O!.lex:=arg[3]!.lex;
            fi;
        fi;
    fi;

    SetOrderingName(O,"block");
    SetComparisonFunction(O,function(a,b)
        local ret,e,f,p,q,o,i,g;
        o:=Length(O!.orders);
        p:=[];
        for i in [1..o] do
            Add(p,[]);
        od;
        q:=StructuralCopy(p);
        g:=GeneratorsOfQuiver(O!.quiver);
        if not IsZeroPath(a) then
            e:=ExtRepOfObj(a);
            for i in e do
                Add(p[block[i]],g[i]);
            od;
        fi;
        if not IsZeroPath(b) then
            f:=ExtRepOfObj(b);
            for i in f do
                Add(q[block[i]],g[i]);
            od;
        fi;
        while o>0 do
            ret:=CompareArrowLists(O!.orders[o],p[o],q[o]);
            if ret<>0 then
                return ret;
            fi;
            o:=o-1;
        od;
        if HasNextOrdering(O) then
            return ComparisonFunction(NextOrdering(O))(a,b);
        fi;
        return 0;
    end);
    return O;
  end
);


DeclareGlobalFunction("CheckBlocksForCycles");

InstallGlobalFunction(CheckBlocksForCycles,
 function(O,block)
     local i,r,n,visit;
     #Possibly a big shortcut
     if HasIsAcyclicQuiver(O!.quiver) and IsAcyclicQuiver(O!.quiver) then
         return true;
     fi;
     n:=[];
     i:=1;
     #Visit a vertex
     #v is the vertex
     #j is 4 * the generator position of the vertex
     #c is the label of the source vertex for the last arrow which transitioned blocks
     #b is the current block number, or 0
     #a is the label of the source vertex for highest previous 0 arrow in a row
     #n[j]   = the depth-first search label of the vertex
     #n[j-1] = true: the vertex was visited successfully; false: vertex is being visited
     #n[j-2] = index of the highest vertex reachable from cycles containing this vertex
     #n[j-3] = the block number of any cycle this vertex belongs to, or 0
     visit:=function(v,j,c,b,a)
         local adj,k,t,p,d,e,f,m;
         n[j]:=i;
         n[j-1]:=true;
         n[j-3]:=0;
         i:=i+1;
         #Check all the arrows
         adj:=OutgoingArrowsOfVertex(v);
         for p in adj do
             t:=TargetOfPath(p);
             k:=ExtRepOfObj(t)[1]*4;
             #Look for block transition
             e:=block[ExtRepOfObj(p)[1]];
             if e>0 then
                 if b>0 and e<>b then
                     f:=a;
                 else
                     f:=c;
                 fi;
                 d:=i;
             else
                 e:=b;
                 f:=c;
                 d:=a;
             fi;
             #Arrow points to an unvisited vertex
             if not IsBound(n[k]) then
                 if not visit(t,k,f,e,d) then
                     return false;
                 #If it's part of a cycle that reaches us
                 elif IsBound(n[k-2]) and n[n[k-2]-1] then
                     #Mark how high we can reach
                     if (not IsBound(n[j-2])) or n[n[j-2]]>n[n[k-2]] then
                         n[j-2]:=n[k-2];
                     fi;
                     #Mark what block the cycle is in
                     if n[k-3]<>0 then
                         n[j-3]:=n[k-3];
                     fi;
                 fi;
             #The vertex can reach us
             elif n[k-1] or (IsBound(n[k-2]) and n[n[k-2]-1]) then
                 if IsBound(n[k-2]) then
                     m:=[k-2];
                 else
                     m:=k;
                 fi;
                 #The cycle contains something outside our block
                 if n[m]<f or (n[k-3]>0 and e>0 and n[k-3]<>e) then
                     return false;
                 #Mark us a belonging to the cycle
                 elif (not IsBound(n[j-2])) or n[n[j-2]]>n[m] then
                     n[j-2]:=m;
                 fi;
             fi;
         od;
         n[j-1]:=false;
         return true;
     end;
     for i in [1..Length(VerticesOfQuiver(O!.quiver))] do
         if block[i]>0 then
             r:=i*4;
             if not IsBound(n[r]) then
                 if not visit(VerticesOfQuiver(O!.quiver)[i],r,0,0,0) then
                     return false;
                 fi;
             fi;
         fi;
     od;
     return true;
  end
);


InstallMethod(IsTotalOrdering,
 "for a block ordering",
 true,
 [IsBlockOrdering],0,
 function(O)
     local l,s,i,block;
     l:=GeneratorsOfQuiver(O!.quiver);
     s:=[];
     block:=O!.block;
     for i in [1..Length(O!.orders)] do
         Add(s,[]);
     od;
     for i in l do
         Add(s[block[ExtRepOfObj(i)[1]]],i);
     od;
     for i in [1..Length(O!.orders)] do
         if not IsTotalSubset(O!.orders[i],s[i]) then
             if HasNextOrdering(O) then
                 return IsTotalOrdering(NextOrdering(O));
             else
                 return false;
             fi;
         fi;
     od;
     if Length(block)>0 then
         if not CheckBlocksForCycles(O,block) then
            if HasNextOrdering(O) then
                return IsTotalOrdering(NextOrdering(O));
            else
                return false;
            fi;
         fi;
     fi;
     return true;
  end
);


InstallMethod(IsWellSubset,
 "for a block ordering and a list of generators",
 true,
 [IsBlockOrdering,IsHomogeneousList],0,
 function(O,L)
     local s,i,block;
     s:=[];
     block:=O!.block;
     for i in [1..Length(O!.orders)] do
         Add(s,[]);
     od;
     for i in L do
         Add(s[block[ExtRepOfObj(i)[1]]],i);
     od;
     for i in [1..Length(O!.orders)] do
         if not IsWellSubset(O!.orders[i],s[i]) then
             return false;
         fi;
     od;
     return true;
  end
);


InstallMethod(IsWellReversedSubset,
 "for a block ordering and a list of generators",
 true,
 [IsBlockOrdering,IsHomogeneousList],0,
 function(O,L)
     local s,i,block;
     s:=[];
     block:=O!.block;
     for i in [1..Length(O!.orders)] do
         Add(s,[]);
     od;
     for i in L do
         Add(s[block[ExtRepOfObj(i)[1]]],i);
     od;
     for i in [1..Length(O!.orders)] do
         if not IsWellReversedSubset(O!.orders[i],s[i]) then
             return false;
         fi;
     od;
     return true;
  end
);


InstallMethod(IsTotalSubset,
 "for a block ordering and a list of generators",
 true,
 [IsBlockOrdering,IsHomogeneousList],0,
 function(O,L)
     local l,s,i,j,block;
     s:=[];
     block:=ListWithIdenticalEntries(Length(O!.block),0);
     for i in [1..Length(O!.orders)] do
         Add(s,[]);
     od;
     for i in L do
         if IsArrow(i) then
             j:=ExtRepOfObj(SourceOfPath(i))[1];
             block[j]:=O!.block[j];
         fi;
         j:=ExtRepOfObj(i)[1];
         block[j]:=O!.block[j];
         Add(s[block[j]],i);
     od;
     for i in [1..Length(O!.orders)] do
         if not IsTotalSubset(O!.orders[i],s[i]) then
             if HasNextOrdering(O) then
                 return IsTotalSubset(NextOrdering(O),L);
             else
                 return false;
             fi;
         fi;
     od;
     if Length(block)>0 then
         if not CheckBlocksForCycles(O,block) then
             if HasNextOrdering(O) then
                 return IsTotalOrdering(NextOrdering(O));
             else
                 return false;
             fi;
         fi;
     fi;
     return true;
  end
);


InstallGlobalFunction(WreathOrdering,function(arg)
    local fam,O,Q,L,block,orders,i,j,rep;
    if not (Length(arg)=1 or (Length(arg)=2 or (Length(arg)=3 and
     IsQuiverOrdering(arg[3])) and IsHomogeneousList(arg[2])) and IsQuiver(arg[1]))
     then
        Error("usage: WreathOrdering(Q,L,O), the first arg must be a quiver, the second a list of [ordering,[generator list]] pairs, and the optional third argument an ordering");
    fi;
    L:=arg[2];
    if Length(L)<1 then
        Error("WreathOrdering: Subset list cannot be empty");
    fi;

    fam:=NewFamily("FamilyWreathOrderings",IsWreathOrdering);
    O:=Objectify(NewType(fam,IsWreathOrdering and IsAttributeStoringRep),
                 rec());
    Q:=arg[1];
    O!.quiver:=Q;

    block:=[];
    orders:=[];
    L:=Flat(L);
    i:=1;
    while i<Length(L) do
        if not IsQuiverOrdering(L[i]) then
            Error("Second argument must be a list of [ordering,[generator list]] pairs");
        fi;
        Add(orders,L[i]);
        if HasLexicographicTable(L[i]) then
            if HasLexicographicTable(O) then
                if LexicographicTable(O)<>LexicographicTable(L[i]) then
                    Error("Conflicting lexicographic tables");
                fi;
            else
                SetLexicographicTable(O,LexicographicTable(L[i]));
                O!.lex:=L[i]!.lex;
            fi;
        fi;
        i:=i+1;
        while i<Length(L) do
            if (not IsArrow(L[i])) and (not IsQuiverVertex(L[i])) then
                break;
            fi;
            block[ExtRepOfObj(L[i])[1]]:=Length(orders);
            i:=i+1;
        od;
    od;
    L:=Length(GeneratorsOfQuiver(Q));
    if L>0 then
        for i in [1..L] do
            if not IsBound(block[i]) then
                block[i]:=Length(orders);
            fi;
        od;
    fi;
    O!.block:=block;
    O!.orders:=orders;

    if Length(arg)>=3 then
        SetNextOrdering(O,arg[3]);
        if HasLexicographicTable(arg[3]) then
            if HasLexicographicTable(O) then
                if LexicographicTable(O)<>LexicographicTable(arg[3]) then
                    Error("Conflicting lexicographic tables");
                fi;
            else
                SetLexicographicTable(O,LexicographicTable(arg[3]));
                O!.lex:=arg[3]!.lex;
            fi;
        fi;
    fi;

    SetOrderingName(O,"wreath");
    SetComparisonFunction(O,function(a,b)
        local ret,e,f,o,s,t,i,j,c,d,p,q,m,n,x,y,g;
        if IsZeroPath(a) then
            e:=[];
        else
            e:=ExtRepOfObj(a);
        fi;
        if IsZeroPath(b) then
            f:=[];
        else
            f:=ExtRepOfObj(b);
        fi;
        o:=1;
        i:=Minimum(Length(e),Length(f));
        # Skip equal elements
        while o<=i do
            if e[o]<>f[o] then
                break;
            fi;
            o:=o+1;
        od;
        # If we skipped everything, the paths are identical
        if o>Length(e) and o>Length(f) then
            return 0;
        fi;
        g:=GeneratorsOfQuiver(O!.quiver);
        s:=0;
        m:=0;
        p:=[];
        c:=[];
        x:=[];
        if o<=Length(e) then
            for i in [o..Length(e)] do
                j:=block[e[i]];
                if j>m then
                    s:=s+1;
                    m:=j;
                    Add(c,m);
                    Add(p,[g[e[i]]]);
                elif j=m then
                    Add(p[s],g[e[i]]);
                else
                    Add(x,g[e[i]]);
                fi;
            od;
        fi;
        t:=0;
        n:=0;
        q:=[];
        d:=[];
        y:=[];
        if o<=Length(f) then
            for i in [o..Length(f)] do
                j:=block[f[i]];
                if j>n then
                    t:=t+1;
                    n:=j;
                    Add(d,n);
                    Add(q,[g[f[i]]]);
                elif j=n then
                    Add(q[t],g[f[i]]);
                else
                    Add(y,g[f[i]]);
                fi;
            od;
        fi;
        while s>0 and t>0 do
            if c[s]=d[t] then
                ret:=CompareArrowLists(orders[c[s]],p[s],q[t]);
                s:=s-1;
                t:=t-1;
            elif c[s]>d[t] then
                ret:=CompareArrowLists(orders[c[s]],p[s],[]);
                s:=s-1;
            else
                ret:=CompareArrowLists(orders[d[t]],[],q[t]);
                t:=t-1;
            fi;
            if ret<>0 then
                return ret;
            fi;
        od;
        while s>0 do
            ret:=CompareArrowLists(orders[c[s]],p[s],[]);
            if ret<>0 then
                return ret;
            fi;
            s:=s-1;
        od;
        while t>0 do
            ret:=CompareArrowLists(orders[d[t]],[],q[t]);
            if ret<>0 then
                return ret;
            fi;
            t:=t-1;
        od;
        ret:=CompareArrowLists(O,x,y);
        if ret=0 and HasNextOrdering(O) then
            ret:=ComparisonFunction(NextOrdering(O))(a,b);
        fi;
        return ret;
    end);

    return O;

  end
);


PrintOrdering :=  function(O)
    local done;
    
    done := false;
    Print("<");
    while not done do
        if HasOrderingName(O) then
            Print(OrderingName(O));
        else
            Print("unnamed");
        fi;
        Print(" ");
        if HasNextOrdering(O) then
            O := NextOrdering(O);
        else
            done := true;
        fi;
    od;
    Print("ordering>");
end;


InstallMethod( ViewObj,
    "for orderings",
    true,
    [IsQuiverOrdering], 0,
    PrintOrdering
);


InstallMethod( PrintObj,
    "for orderings",
    true,
    [IsQuiverOrdering], 0,
    PrintOrdering
);
