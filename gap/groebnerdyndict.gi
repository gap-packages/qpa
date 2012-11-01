# GAP Implementation
# This file was generated from 
# $Id: dyndict.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallGlobalFunction( CreateSuffixTree, function()
    local tree, family, elementFamily, root;

    family := NewFamily("SuffixTreeFamily", IsSuffixTree);
    elementFamily := NewFamily("SuffixTreeElementFamily",
                               IsSuffixTreeElement);
    root := CreateInternalSuffixTreeNode(elementFamily, [], 0, 0);
    root!.dynTree := CreateDynamicTree(root);
    tree := Objectify( NewType(family,
                               IsSuffixTree and IsDefaultSuffixTreeRep),
                       rec( family := elementFamily,
                            patterns := [],
                            root := root,
                            deadLength := 0 ) );
     return tree;
end );
InstallGlobalFunction( CreateInternalSuffixTreeNode, 
function(family, pattern, start, lend)
    local node;
    node := Objectify( NewType( family,
        IsInternalSuffixTreeNode and IsDefaultInternalSuffixTreeNodeRep ),
        rec( marked := false,
             parentEdge := 0,
             label := [pattern, start, lend],
             suffixLink := 0,
             edges := [] ) );
    return node;
end );
InstallGlobalFunction( CreateExternalSuffixTreeNode,
function(family, pattern, start)
    local node;
    node := Objectify( NewType( family,
        IsExternalSuffixTreeNode and IsDefaultSuffixTreeNodeRep ),
        rec( parentEdge := 0,
             label := [ [pattern,start] ] ) );
    return node;
end );
InstallGlobalFunction( CreateSuffixTreeEdge,
function( family, parent, child, pattern, start, lend )
    local edge;

    edge := Objectify( NewType( family,
        IsSuffixTreeEdge and IsDefaultSuffixTreeEdgeRep),
        rec( parent := parent,
             child := child,
             label := [pattern, start, lend] ) );
   
   # update the parent and child
   Add(parent!.edges, edge);
   child!.parentEdge := edge;
   return edge;
end );
InstallMethod( ParentNode,
    "for suffix tree nodes",
    true,
    [ IsSuffixTreeNode and IsDefaultSuffixTreeNodeRep ], 0,
    function( v )
        if v!.parentEdge = 0 then
            return 0;
        else
            return v!.parentEdge!.parent;
        fi;
    end );
InstallMethod( \=,
    "for suffix tree elements",
    IsIdenticalObj,
    [IsSuffixTreeElement, IsSuffixTreeElement], 0,
    IsIdenticalObj );
InstallOtherMethod( \[\],
    "for suffix trees and positive integers",
    true,
    [IsSuffixTree and IsDefaultSuffixTreeRep,
     IsPosInt], 0,
    function( st, i )
        local n;
        n := Length(st!.patterns[i][2]);
        return st!.patterns[i][2]{[1..n-1]};
    end );
InstallMethod( AddInternalSuffixTreeNode,
    "for suffix tree edges, integers, and internal suffix tree nodes",
    true,
    [IsSuffixTreeEdge and IsDefaultSuffixTreeEdgeRep,
     IsPosInt,
     IsInternalSuffixTreeNode and IsDefaultInternalSuffixTreeNodeRep], 0,
    function( edge, location, node )
        local oldChild, pattern, start, lend;

        # get the information needed to construct the label on the
        # new edge. This must happen before updating the old edge
        pattern := edge!.label[1];
        start := edge!.label[2] + location;
        lend := edge!.label[3];

        # update the old edge
        oldChild := edge!.child;
        edge!.child := node;
        node!.parentEdge := edge;
        edge!.label[3] := edge!.label[2] + location - 1;

        # create a new edge
        CreateSuffixTreeEdge(FamilyObj(edge), node, 
                             oldChild, pattern, start, lend);
   end);
InstallMethod( AddExternalSuffixTreeNode,
    "for internal and external suffix tree nodes",
    IsIdenticalObj,
    [IsInternalSuffixTreeNode and IsDefaultInternalSuffixTreeNodeRep,
     IsExternalSuffixTreeNode and IsDefaultSuffixTreeNodeRep], 0,
    function( parent, child )
        local pattern, start, lend;

        pattern := child!.label[1][1];
        start := Length(LabelOfSuffixTreeNode(parent)) + child!.label[1][2];
        lend := -1;

        CreateSuffixTreeEdge(FamilyObj(parent), parent, child,
                             pattern, start, lend);
    end );
InstallMethod( LabelOfSuffixTreeNode,
    "for external suffix tree nodes",
    true,
    [ IsExternalSuffixTreeNode and IsDefaultSuffixTreeNodeRep ], 0,
    function( node )
        local pattern, start, lend;

        pattern := node!.label[1][1][2];
        start := node!.label[1][2];
        lend := Length(pattern);

        return pattern{[start..lend]};
    end );
InstallMethod( LabelOfSuffixTreeNode,
    "for internal suffix tree nodes",
    true,
    [ IsInternalSuffixTreeNode and IsDefaultInternalSuffixTreeNodeRep ], 0,
    function( node )
        local pattern;

        if Length(node!.label[1]) = 0 then
            # special case for root node
            return [];
        fi;

        pattern := node!.label[1][2];

        return pattern{[node!.label[2]..node!.label[3]]};
    end );
InstallMethod( LabelOfSuffixTreeEdge,
    "for suffix tree edges",
    true,
    [ IsSuffixTreeEdge and IsDefaultSuffixTreeEdgeRep ], 0,
    function(edge)
        local pattern, start, lend;

        pattern := edge!.label[1][2];
        start := edge!.label[2];
        lend := edge!.label[3];
        if lend < start then
           # add on the termination symbol
           lend := Length(pattern);
        fi;
        Assert( 1, start <= Length(pattern), 
                "Start of pattern is longer than pattern" );
        return pattern{[start..lend]};
    end );
InstallMethod( SuffixTreeEdgeStartingAt,
    "for an internal suffix tree node and an object",
    true,
    [IsInternalSuffixTreeNode and IsDefaultInternalSuffixTreeNodeRep,
     IsObject], 0,
    function( node, object )
        local edge, label;

        for edge in node!.edges do
            label := LabelOfSuffixTreeEdge(edge);
            if Length(label) > 0 and object = label[1] then
                return edge;
            fi;
        od;
        return fail;
    end );
InsertIntoDynamicForest := function(u)
    local v, w;

    v := ParentNode(u);

    if IsExternalSuffixTreeNode(u) then
        u!.dynTree := CreateDynamicTree(u);
        LinkDynamicTrees(u!.dynTree, v!.dynTree);
    elif u!.marked and IsBound(u!.dynTree) then
        Assert(1, IsInternalSuffixTreeNode(u), 
               "$u$ is marked and old, should be internal");
        CutDynamicTree(u!.dynTree);
    else
        Assert(1, IsInternalSuffixTreeNode(u) and not IsBound(u!.dynTree),
               "$u$ should be a new internal node here!");

        u!.dynTree := CreateDynamicTree(u);
        # NOTE: Major assumption here! The edge containing the child
        # that was not newly created is the first one. If the order in
        # which edges are created when a new internal node is put into
        # the suffix tree, then this will probably fail to be true
        w := u!.edges[1]!.child;

        if not u!.marked then
            if RootOfDynamicTree(v!.dynTree) = RootOfDynamicTree(w!.dynTree)
            then
                CutDynamicTree(w!.dynTree);
                LinkDynamicTrees(u!.dynTree, v!.dynTree);
                LinkDynamicTrees(w!.dynTree, u!.dynTree);
            else
                LinkDynamicTrees(u!.dynTree, v!.dynTree);
            fi;
        else
            if RootOfDynamicTree(v!.dynTree) = RootOfDynamicTree(w!.dynTree)
            then
                CutDynamicTree(w!.dynTree);
                LinkDynamicTrees(w!.dynTree, u!.dynTree);
            fi;
        fi;
    fi;
end;
ApplySuffixTreeExtension := function(where, pattern, i, j)
    local currentNode, currentEdge, location, label, newNode,
          nextChar, family, retval, termEdge, justMarked;

    justMarked := false;
    currentNode := 0; currentEdge := 0; location := 0; newNode := 0;
    nextChar := pattern[2][i+1];
    retval := rec( internal := false, currentEdge := fail, location := 0 );

    if IsList(where) then
        currentEdge := where[1];
        location := where[2];
        label := LabelOfSuffixTreeEdge(currentEdge);
        if label[location+1] <> nextChar then
            family := FamilyObj(currentEdge);
            # Have to add a new internal node and external node
            currentNode := CreateInternalSuffixTreeNode(family, pattern, j, i);
            AddInternalSuffixTreeNode(currentEdge, location, currentNode);
            Info( InfoDynamicDictionary, 1, 
                  "Added new internal node with edge label ",
                  LabelOfSuffixTreeEdge(currentEdge) );
            newNode := CreateExternalSuffixTreeNode(family, pattern, j);
            AddExternalSuffixTreeNode(currentNode, newNode);
            Info( InfoDynamicDictionary, 1,
                  "Added new external node with label ",
                  LabelOfSuffixTreeNode(newNode) );
            retval.internal := true;
            if nextChar = [] then
                termEdge := SuffixTreeEdgeStartingAt(currentNode, []);
                if ForAny(termEdge!.child!.label, x -> x[2] = 1) then
                    currentNode!.marked := true;
                    Info(InfoDynamicDictionary, 1, "Marked node: ",
                         LabelOfSuffixTreeNode(currentNode));
                    justMarked := true;
                fi;
            fi;
            InsertIntoDynamicForest(currentNode);
            InsertIntoDynamicForest(newNode);
        else
            currentNode := currentEdge!.parent;
            if location+1 < Length(label) then
                retval.currentEdge := currentEdge;
                retval.location := location+1;
                retval.nextNode := currentNode;
            else
                retval.nextNode := currentEdge!.child;
            fi;
        fi;
    else
        currentNode := where;
        currentEdge := SuffixTreeEdgeStartingAt(currentNode, nextChar);
        if currentEdge = fail then
            # Have to add a new external node
            Assert(1, not IsExternalSuffixTreeNode(currentNode),
                   "Trying to extend an external node!");
            family := FamilyObj(currentNode);
            newNode := CreateExternalSuffixTreeNode(family, pattern, j);
            AddExternalSuffixTreeNode(currentNode, newNode);
            Info( InfoDynamicDictionary, 1,
                  "Added new external node with label ",
                  LabelOfSuffixTreeNode(newNode) );
            retval.nextNode := currentNode;
            InsertIntoDynamicForest(newNode);
        else
            if Length(LabelOfSuffixTreeEdge(currentEdge)) > 1 then
                retval.currentEdge := currentEdge;
                retval.location := 1;
                retval.nextNode := currentNode;
            else
                retval.nextNode := currentEdge!.child;
            fi;
        fi;
    fi;

    if nextChar = [] and newNode = 0 then
        Assert(1, IsExternalSuffixTreeNode(currentEdge!.child),
               "Match of terminating character not at leaf!");
        Add(currentEdge!.child!.label, [pattern, j]);
        # Mark an old node
        if j = 1 and
           Length(LabelOfSuffixTreeNode(currentEdge!.parent)) = i then
            currentEdge!.parent!.marked := true;
            Info(InfoDynamicDictionary, 1, "Marked node: ",
                 LabelOfSuffixTreeNode(currentEdge!.parent));
            InsertIntoDynamicForest(currentEdge!.parent);
        fi;
        # also, fake as if we added a new leaf so the phase doesn't
        # prematurely end.
        newNode := currentEdge!.child;
    fi;

    retval.currentNode := currentNode;
    retval.newNode := newNode;

    return retval;
end;    
InstallMethod( InsertPatternIntoSuffixTree,
    "for suffix trees in default representation and lists",
    true,
    [IsSuffixTree and IsDefaultSuffixTreeRep, IsList], 0,
    function(suffixTree, pattern)
        local i, j, j_i, m, label, location,
              phaseDone, extInfo, newPattern,
              currentEdge, currentNode, nextNode,
              v, w, x, g, g1, h;

        pattern := Concatenation(pattern, [[]]);
        m := Length(pattern);
        newPattern := [Length(suffixTree!.patterns)+1, pattern];

        i := 0;
        currentNode := suffixTree!.root;
        nextNode := currentNode;
        while i < m-1 do
            location := 0;
            currentEdge := SuffixTreeEdgeStartingAt(currentNode, pattern[i+1]);
            if currentEdge <> fail then
                label := LabelOfSuffixTreeEdge(currentEdge);
                nextNode := currentEdge!.child;
                if i < m then
                    i := i + 1;
                    location := 1;
                fi;
            else
                # failed to match on any of the children, so get out
                break; 
            fi;

            # move as far as we can down the edge
            while location < Length(label) and i < m-1 
                  and pattern[i+1] = label[location+1] do
                i := i + 1;
                location := location + 1;
            od;

            # If we finished matching the current label, move on to the next
            # node. Note that we do not need to check if this is an internal
            # or external node since external edges have the terminating
            # symbol.
            if location = Length(label) then
                currentNode := nextNode;
                location := 0;
            else
                # we ended in the middle of an edge so break out.
                break;
            fi;
        od;

        Info( InfoDynamicDictionary, 1, "Matched ", i, " items.");

        j_i := 0;
        while i < m do
            phaseDone := false;
            j := j_i+1;
            Info( InfoDynamicDictionary, 1,
                  "Begin phase ", i, " extension ", j );
            Info( InfoDynamicDictionary, 1,
                  "Current node: ", LabelOfSuffixTreeNode(currentNode) );
            if currentEdge <> fail then
                Info( InfoDynamicDictionary, 1,
                      "Current edge: ", LabelOfSuffixTreeEdge(currentEdge) );
                Info( InfoDynamicDictionary, 1,
                      "Location: ", location );
            fi;
            if location <> 0 then
                extInfo := ApplySuffixTreeExtension([currentEdge, location], 
                                                    newPattern, i, j);
            else
                extInfo := ApplySuffixTreeExtension(currentNode,
                                                    newPattern, i, j);
            fi;
            if extInfo.internal then
                w := extInfo.currentNode;
            else
                w := 0;
            fi;
            v := extInfo.currentNode;
            if extInfo.newNode = 0 then
                j_i := j - 1;
                currentNode := extInfo.nextNode;
                currentEdge := extInfo.currentEdge;
                location := extInfo.location;
                phaseDone := true;
            else
                # increment j_i since we added a new leaf
                j_i := j_i + 1;
                currentNode := suffixTree!.root;
                location := 0;
            fi;
            j := j + 1;
            while j <= i+1 and not phaseDone do
                Info( InfoDynamicDictionary, 1,
                      "Begin phase ", i, " extension ", j );
                Info( InfoDynamicDictionary, 1,
                     "Current node: ", LabelOfSuffixTreeNode(currentNode) );
                if currentEdge <> fail then
                    Info( InfoDynamicDictionary, 1,
                          "Current edge: ", 
                          LabelOfSuffixTreeEdge(currentEdge) );
                    Info( InfoDynamicDictionary, 1,
                          "Location: ", location );
                fi;
                # assumption, v is not the root right now
                if v!.suffixLink = 0 and v <> suffixTree!.root then
                    label := LabelOfSuffixTreeEdge(v!.parentEdge);
                    v := v!.parentEdge!.parent;
                else
                    label := [];
                fi;
                if v <> suffixTree!.root then
                    currentNode := v!.suffixLink;
                else
                    label := pattern{[j..i]};
                    currentNode := v;
                fi;
                g := Length(label);
                h := 1;
                location := 0;

                while g > 0 do
                    currentEdge := SuffixTreeEdgeStartingAt(currentNode, label[h]);
                    g1 := Length(LabelOfSuffixTreeEdge(currentEdge));
                    if g1 < g then
                        h := h + g1;
                        currentNode := currentEdge!.child;
                    elif g = g1 then
                        currentNode := currentEdge!.child;
                        nextNode := currentNode;
                        location := 0;
                    else
                        nextNode := currentEdge!.child;
                        location := g;
                    fi;
                    g := g - g1;
                od;
                if location <> 0 then
                    extInfo := ApplySuffixTreeExtension([currentEdge, location], 
                                                        newPattern, i, j);
                else
                    extInfo := ApplySuffixTreeExtension(currentNode,
                                                        newPattern, i, j);
                fi;
                if extInfo.internal then
                    x := extInfo.currentNode;
                else
                    x := 0;
                fi;
                v := extInfo.currentNode;
                if extInfo.newNode = 0 then
                    j_i := j - 1;
                    currentNode := extInfo.nextNode;
                    currentEdge := extInfo.currentEdge;
                    location := extInfo.location;
                    phaseDone := true;
                else
                    # increment j_i since we added a new leaf
                    j_i := j_i + 1;
                    currentNode := suffixTree!.root;
                    location := 0;
                fi;
                if w <> 0 then
                    Assert(1, LabelOfSuffixTreeNode(v) = pattern{[j..i]}, 
                           "Incorrectly setting suffix link!");
                    w!.suffixLink := v;
                fi;
                w := x;
                j := j + 1;
            od;
            i := i+1;
        od;

        Add(suffixTree!.patterns, newPattern);
        return Length(suffixTree!.patterns);
    end );
DeleteFromDynamicForest := function(u, deleted, q)
    local v;

    v := ParentNode(u);

    if IsExternalSuffixTreeNode(u) then
        CutDynamicTree(u!.dynTree);
    elif deleted then
        CutDynamicTree(u!.dynTree);
        CutDynamicTree(q!.dynTree);
        LinkDynamicTrees(q!.dynTree, v!.dynTree); 
    else
        LinkDynamicTrees(u!.dynTree, v!.dynTree);
    fi;
end;
InstallMethod( DeletePatternFromSuffixTree,
    "for suffix trees and positive integers",
    true,
    [IsSuffixTree and IsDefaultSuffixTreeRep, IsPosInt], 0,
    function( st, patternNum )
        local currentNode, currentEdge, nextNode, location, pos,
              g, g1, h, pattern, label, i, m, u, v, w, children,
              parent, n;

        pattern := st!.patterns[patternNum];

        currentNode := st!.root;
        label := pattern[2];
        m := Length(label);

        for i in [1..m] do
            g := Length(label);
            h := 1;
            location := 0;

            while g > 0 do
                currentEdge := SuffixTreeEdgeStartingAt(currentNode, label[h]);
                g1 := Length(LabelOfSuffixTreeEdge(currentEdge));
                if g1 < g then
                    h := h + g1;
                    currentNode := currentEdge!.child;
                elif g = g1 then
                    currentNode := currentEdge!.child;
                    nextNode := currentNode;
                    location := 0;
                else
                    nextNode := currentEdge!.child;
                    location := g;
                fi;
                g := g - g1;
            od;
            u := currentNode;
            v := ParentNode(u);
            children := Sum(v!.edges, function(x)
                if IsExternalSuffixTreeNode(x!.child) then
                    return Length(x!.child!.label);
                else
                    return 1;
                fi;
            end );
            if v = st!.root then
                Assert(1, IsExternalSuffixTreeNode(u), "$u$ should be an external node!");
                Info(InfoDynamicDictionary, 1, "Deleting node $u$: ", 
                     LabelOfSuffixTreeNode(u));
                Info(InfoDynamicDictionary, 1, "Length of $u$ label: ", Length(u!.label));
                if Length(u!.label) > 1 then
                    u!.label := Filtered(u!.label, x -> x[1][1] <> patternNum);
                else
                    DeleteFromDynamicForest(u, true, 0);
                    v!.edges := Filtered(v!.edges, x -> x <> u!.parentEdge);
                    u!.parentEdge!.parent := 0;
                    u!.parentEdge!.child := 0;
                fi;
                if i = 1 
                   and LabelOfSuffixTreeEdge(u!.parentEdge) = [[]] 
                   and v!.marked
                then
                    Info(InfoDynamicDictionary, 1, "Unmarked: ", LabelOfSuffixTreeNode(v));
                    v!.marked := false;
                    DeleteFromDynamicForest(v, false, 0);
                fi;
            elif children > 2 then
                Assert(1, IsExternalSuffixTreeNode(u), "$u$ should be an external node!");
                Info(InfoDynamicDictionary, 1, "Deleting node $u$: ", 
                     LabelOfSuffixTreeNode(u));
                Info(InfoDynamicDictionary, 1, "Length of $u$ label: ", Length(u!.label));
                if Length(u!.label) > 1 then
                    u!.label := Filtered(u!.label, x -> x[1][1] <> patternNum);
                else
                    DeleteFromDynamicForest(u, true, 0);
                    v!.edges := Filtered(v!.edges, x -> x <> u!.parentEdge);
                    u!.parentEdge!.parent := 0;
                    u!.parentEdge!.child := 0;
                fi;
                if i = 1 
                   and LabelOfSuffixTreeEdge(u!.parentEdge) = [[]] 
                   and v!.marked
                then
                    Info(InfoDynamicDictionary, 1, "Unmarked: ", LabelOfSuffixTreeNode(v));
                    v!.marked := false;
                    DeleteFromDynamicForest(v, false, 0);
                fi;
            else
                Assert(1, IsExternalSuffixTreeNode(u), "$u$ should be an external node!");
                Info(InfoDynamicDictionary, 1, "Deleting node $u$: ", 
                     LabelOfSuffixTreeNode(u));
                Info(InfoDynamicDictionary, 1, "Length of $u$ label: ", Length(u!.label));
                if Length(u!.label) > 1 then
                    u!.label := Filtered(u!.label, x -> x[1][1] <> patternNum);
                else
                    DeleteFromDynamicForest(u, true, 0);
                    v!.edges := Filtered(v!.edges, x -> x <> u!.parentEdge);
                    u!.parentEdge!.parent := 0;
                    u!.parentEdge!.child := 0;
                fi;
                if i = 1 
                   and LabelOfSuffixTreeEdge(u!.parentEdge) = [[]] 
                   and v!.marked
                then
                    Info(InfoDynamicDictionary, 1, "Unmarked: ", LabelOfSuffixTreeNode(v));
                    v!.marked := false;
                    DeleteFromDynamicForest(v, false, 0);
                fi;
                Assert(1, Length(v!.edges) = 1, "$v$ is only supposed to have one child!");
                Info(InfoDynamicDictionary, 1, "Deleting node $v$: ", 
                     LabelOfSuffixTreeNode(v));
                w := v!.edges[1]!.child;
                DeleteFromDynamicForest(v, true, w);
                parent := ParentNode(v);
                n := Length(LabelOfSuffixTreeEdge(v!.parentEdge));
                v!.parentEdge!.label := w!.parentEdge!.label;
                v!.parentEdge!.label[2] := v!.parentEdge!.label[2] - n;
                v!.parentEdge!.child := w;
                w!.parentEdge := v!.parentEdge;
                v!.parentEdge := 0;
                Info(InfoDynamicDictionary, 1, "Deleted edge, new label is",
                     LabelOfSuffixTreeEdge(w!.parentEdge) );
            fi;
            if i < m then
                if v = st!.root then
                    currentNode := v;
                    label := label{[2..Length(label)]};
                else
                    currentNode := v!.suffixLink;
                    label := LabelOfSuffixTreeEdge(u!.parentEdge);
                fi;
            fi;
       od;
    end );
InstallMethod( AllPatternsInSequence,
    "for suffix trees and lists",
    true,
    [IsSuffixTree and IsDefaultSuffixTreeRep, IsList], 0,
    function( st, s )
        local clocus, elocus, alpha, beta, beta_h, j, k, n, x, y, u, v,
              occ, currentEdge, lastNode, currentLabel, location;

        n := Length(s);
        occ := [];
        clocus := st!.root;
        beta := [];
        k := 1;

        for j in [1..n] do
            if clocus <> st!.root then
                x := clocus!.suffixLink;
                beta_h := beta;
                alpha := [];
                while Length(alpha) < Length(beta_h) do
                    beta_h := beta_h{[Length(alpha)+1..Length(beta_h)]};
                    currentEdge := SuffixTreeEdgeStartingAt(x, beta_h[1]);
                    Assert(1, currentEdge <> fail, "Failed to find a guaranteed path!");
                    alpha := LabelOfSuffixTreeEdge(currentEdge);
                    if Length(alpha) < Length(beta_h) then
                        x := currentEdge!.child;
                    fi;
                od;
                if Length(alpha) > Length(beta_h) then
                    clocus := x;
                    beta := beta_h;
                else
                    if alpha = [] then
                        y := x;
                    else
                        y := currentEdge!.child;
                    fi;
                    lastNode := y;
                    k := j + Length(LabelOfSuffixTreeNode(y));
                    if k <= Length(s) then
                        currentEdge := SuffixTreeEdgeStartingAt(y, s[k]);
                    else
                        currentEdge := fail;
                    fi;
                    while currentEdge <> fail do
                        y := currentEdge!.child;
                        location := 2;
                        k := k + 1;
                        currentLabel := LabelOfSuffixTreeEdge(currentEdge);
                        while location <= Length(currentLabel) and k <= Length(s) do
                            if currentLabel[location] = s[k] then
                                k := k + 1;
                                location := location + 1;
                            else
                                break;
                            fi;
                        od;
                        if location > Length(currentLabel) then
                            lastNode := y;
                            if k <= Length(s) then
                                currentEdge := SuffixTreeEdgeStartingAt(y, s[k]);
                            else
                                currentEdge := fail;
                            fi;
                        else
                            currentEdge := fail;
                        fi;
                    od;
                    clocus := lastNode;
                    if j+Length(LabelOfSuffixTreeNode(clocus)) <= Length(s) then
                        beta := s{[j+Length(LabelOfSuffixTreeNode(clocus))..k-1]};
                    else
                        beta := [];
                    fi;
                fi;
            else
                y := st!.root;
                beta_h := beta{[2..Length(beta)]};
                lastNode := y;
                k := j + Length(LabelOfSuffixTreeNode(y));
                if k <= Length(s) then
                    currentEdge := SuffixTreeEdgeStartingAt(y, s[k]);
                else
                    currentEdge := fail;
                fi;
                while currentEdge <> fail do
                    y := currentEdge!.child;
                    location := 2;
                    k := k + 1;
                    currentLabel := LabelOfSuffixTreeEdge(currentEdge);
                    while location <= Length(currentLabel) and k <= Length(s) do
                        if currentLabel[location] = s[k] then
                            k := k + 1;
                            location := location + 1;
                        else
                            break;
                        fi;
                    od;
                    if location > Length(currentLabel) then
                        lastNode := y;
                        if k <= Length(s) then
                            currentEdge := SuffixTreeEdgeStartingAt(y, s[k]);
                        else
                            currentEdge := fail;
                        fi;
                    else
                        currentEdge := fail;
                    fi;
                od;
                clocus := lastNode;
                if j+Length(LabelOfSuffixTreeNode(clocus)) <= Length(s) then
                    beta := s{[j+Length(LabelOfSuffixTreeNode(clocus))..k-1]};
                else
                    beta := [];
                fi;
            fi;
            Info(InfoDynamicDictionary, 1, "Label of clocus: ", 
                 LabelOfSuffixTreeNode(clocus));
            Info(InfoDynamicDictionary, 1, "Beta: ", beta);
            if not IsEmpty(beta) then
                currentEdge := SuffixTreeEdgeStartingAt(clocus, beta[1]);
                if currentEdge <> fail then
                    elocus := currentEdge!.child;
                    if IsExternalSuffixTreeNode(elocus) 
                       and Length(LabelOfSuffixTreeEdge(currentEdge)) = Length(beta)+1
                    then
                        Append(occ, List(Filtered(elocus!.label, z -> z[2] = 1),
                                         t -> [t[1][1], j]));
                    fi;
                fi;
            fi;
            u := clocus;
            while u <> st!.root do
                v := RootOfDynamicTree(u!.dynTree)!.item;
                if v <> st!.root then
                    currentEdge := SuffixTreeEdgeStartingAt(v, []);
                    elocus := currentEdge!.child;
                    Append(occ, List(Filtered(elocus!.label, z -> z[2] = 1),
                                     t -> [t[1][1], j]));
                    u := ParentNode(v);
                else
                    u := st!.root;
                fi;
            od;
        od;

        return occ;
    end );
InstallMethod( SequenceInPatterns,
    "for suffix trees and lists",
    true,
    [IsSuffixTree and IsDefaultSuffixTreeRep, IsList], 0,
    function( st, seq )
        local currentNode, currentEdge, i, m, n, location,
              currentLabel, occ, InOrder;

        occ := [];
        currentNode := st!.root;
        m := Length(seq);
        i := 0;
        location := 0;
        while i < m do
            if location = 0 then
                currentEdge := SuffixTreeEdgeStartingAt(currentNode, seq[i+1]);
                if currentEdge = fail then
                    break;
                fi;
                location := 1;
                currentLabel := LabelOfSuffixTreeEdge(currentEdge);
                n := Length(currentLabel);
            elif seq[i+1] = currentLabel[location+1] then
                location := location + 1;
            else
                break;
            fi;
            if location >= n then
                currentNode := currentEdge!.child;
                location := 0;
            fi;
            i := i + 1;
        od;

        if i = m then
            if location <> 0 then
                currentNode := currentEdge!.child;
            fi;

            InOrder := function(node)
                local edge;

                if IsExternalSuffixTreeNode(node) then
                    Append(occ, List(node!.label, x -> [x[1][1], x[2]]));
                else
                    for edge in node!.edges do
                        InOrder(edge!.child);
                    od;
                fi;
            end;

            InOrder(currentNode);
        fi;

        return occ;
    end );
