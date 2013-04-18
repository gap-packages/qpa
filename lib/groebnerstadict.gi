# GAP Implementation
# This file was generated from 
# $Id: stadict.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallGlobalFunction( StaticDictionary, function(patterns)
    local pattern, tree, dict, i, l, newNode, char, v, w, q, fam,
          parent, currentNode;

    if not IsDenseList(patterns) then
        Error("<patterns> must be a dense list");
    fi;

    tree := rec( labels := [],
                children := [],
                patterns := [] );

    tree.failureNode := tree;

    for i in [1..Length(patterns)] do
        pattern := patterns[i];
        currentNode := tree;
        for char in pattern do
            l := Position(currentNode.labels, char);
            if l = fail then
                Add(currentNode.labels, char);
                newNode := rec( label := char,
                                labels := [],
                                children := [],
                                parent := currentNode,
                                patterns := [] );
                Add(currentNode.children, newNode);
                currentNode := newNode;
            else
                currentNode := currentNode.children[l];
            fi;
        od;
        Add(currentNode.patterns, i);
    od;
    q := [];
    List(tree.children, function(x)
        x.failureNode := tree;
        Append(q, x.children);
        return x;
    end );

    while not IsEmpty(q) do
        v := q[1];
        q := q{[2..Length(q)]};
        parent := v.parent;
        char := v.label;
        w := parent.failureNode;
        while w <> tree 
              and Position(w.labels, char) = fail
        do
            w := w.failureNode;
        od;
        l := Position(w.labels, char);
        if l = fail then
            v.failureNode := tree;
        else
            v.failureNode := w.children[l];
        fi;
        if not IsEmpty(v.failureNode.patterns) then
            v.outputNode := v.failureNode;
        elif IsBound(v.failureNode.outputNode) then
            v.outputNode := v.failureNode.outputNode;
        fi;
        Append(q, v.children);
    od;
        fam := NewFamily( "StaticDictionaryFamily", IsStaticDictionary );
        dict := Objectify( 
                    NewType( fam, IsStaticDictionary
                                  and IsStaticDictionaryDefaultRep ),
                    rec( tree := tree ) ); 
        SetPatterns( dict, Immutable(patterns) );

    return dict;
end );
InstallMethod( Search,
    "for static dictionaries in their default rep.",
    true,
    [IsStaticDictionary and IsStaticDictionaryDefaultRep,
     IsDenseList], 10,
    function(dict, text)
        local c, w, m, j, k, currentNode, matches, retval;

        c := 1;
        w := dict!.tree;
        m := Length(text);
        retval := [];

        repeat
            j := Position(w.labels, text[c]);
            while j <> fail do
                currentNode := w.children[j];
                matches := [];
                Append(matches,currentNode.patterns);
                while IsBound(currentNode.outputNode) do
                    currentNode := currentNode.outputNode;
                    Append(matches,currentNode.patterns);
                od;
                if not IsEmpty(matches) then
                    Add(retval, [c, matches]);
                fi;
                w := w.children[j];
                c := c + 1;
                if c <= m then
                    j := Position(w.labels, text[c]);
                else
                    j := fail;
                fi;
            od;
            if w = dict!.tree then
                c := c + 1;
            fi;
            w := w.failureNode;
         until c > m;

         return retval;
    end );
InstallMethod( ViewObj,
    "for static dictionaries",
    true,
    [IsStaticDictionary], 0,
    function(dict)
        local l;
        Print("<static dictionary containing ");
        l := Length(Patterns(dict));
        Print(l);
        Print(" pattern");
        if l > 1 then
            Print("s");
        fi;
        Print(">");
    end );
InstallMethod( ViewObj,
    "for static dictionaries",
    true,
    [IsQuiverStaticDictionary], 0,
    function(dict)
        local l;
        Print("<static dictionary for quivers containing ");
        l := Length(Patterns(dict));
        Print(l);
        Print(" path");
        if l > 1 then
            Print("s");
        fi;
        Print(">");
  end
);


InstallMethod( IsFiniteDifference,
  "for static dictionaries for quivers",
  true,
  [IsQuiverStaticDictionary and IsQuiverStaticDictionaryDefaultRep], 0,
  function( dict )
    local tree, q, cycleExists, target, arrow, parent, outgoing,
          l, fNode, vertex, gray, black, visit, tsorted, i, currentNode;

    tree := dict!.tree;
    Info( InfoStaticDictionary, 1, "Creating finite automata...");
    q := [];

    for i in [1..Length(tree.children)] do
      currentNode := tree.children[i];
      vertex := tree.labels[i];
      Append(q, Filtered(currentNode.children, x -> not IsEmpty(x.children)));

      for outgoing in OutgoingArrowsOfVertex(vertex) do
        l := Position(currentNode.labels, outgoing);
        if l = fail then
          Add(currentNode.labels, outgoing);
          l := Position(tree.labels, TargetOfPath(outgoing));
          Add(currentNode.children, tree.children[l]);
        fi;
      od;
    od;

    while not IsEmpty(q) do
      currentNode := q[1];
      q := q{[2..Length(q)]};
      Append(q, Filtered(currentNode.children, x -> IsEmpty(x.patterns)));
      parent := currentNode.parent;

      arrow := currentNode.label;
      target := TargetOfPath(arrow);
      for outgoing in OutgoingArrowsOfVertex(target) do
        l := Position(currentNode.labels, outgoing);
        if l = fail then
          fNode := currentNode.failureNode;
          l := Position(fNode.labels, outgoing);
          Assert(1, l <> fail, "Outgoing arrow wasn't found!");
          Add(currentNode.labels, outgoing);
          Add(currentNode.children, fNode.children[l]);
        fi;
      od;
    od;

    Info( InfoStaticDictionary, 1, "Checking for cycle...");
    gray := 1;
    black := 2;
    cycleExists := false;
    tsorted := [];

    visit := function(u)
      local v;

      u.color := gray;
      for v in u.children do
        if not IsBound(v.color) then
          visit(v);
        elif v.color = gray then
          cycleExists := true;
        fi;
      od;

      Add(tsorted, u);
      u.color := black;
    end;

    visit(tree);

    if not cycleExists then
      tsorted := Reversed(tsorted);
      dict!.tsorted := tsorted;
    fi;

    return not cycleExists;

  end
);


InstallMethod( DifferenceSize,
    "for static dictionaries for quivers",
    true,
    [IsQuiverStaticDictionary and IsQuiverStaticDictionaryDefaultRep], 0,
    function(dict)
        local tree, tsorted, total, node, child;

        if IsFiniteDifference(dict) then
            tree := dict!.tree;
            tsorted := dict!.tsorted;
            total := -1;                # don't count the root node

            for node in tsorted do
                 node.pathCount := 0;
            od;

            tree.pathCount := 1;

            for node in tsorted do
                for child in node.children do
                    child.pathCount := child.pathCount + node.pathCount;
                od;
                # Leaves accepting patterns are tips, so don't count them.
                if IsEmpty(node.patterns) 
                then
                    total := total + node.pathCount;
                fi;
                # Free memory (later)
                #Unbind(node.pathCount);
            od;
            return total;
        else
            return infinity;
        fi;
    end );
InstallMethod( DifferenceWords,
    "for static dictionaries for quivers",
    true,
    [IsQuiverStaticDictionary and IsQuiverStaticDictionaryDefaultRep], 0,
    function( dict )
        local tree, tsorted, paths, node, arrow, path, i, newPaths;

        if IsFiniteDifference(dict) then
            tree := dict!.tree;
            tsorted := dict!.tsorted;
            paths := [];

            for node in tsorted do
                 node.paths := [];
            od;

            for i in [1..Length(tree.children)] do
                tree.children[i].paths := [tree.labels[i]];
            od;

            for node in tsorted do
                for i in [1..Length(node.children)] do
                    newPaths := List(node.paths, x -> x*node.labels[i]);
                    Append(node.children[i].paths, newPaths);
                od;

                # Leaves accepting patterns are tips, so don't count them.
                if IsEmpty(node.patterns)
                then
                    Append(paths, node.paths);
                fi;
                # Free memory
                # Unbind(node.paths);
            od;
            Sort(paths);
            return paths;
        else
            return fail;
        fi;

    end
);


InstallGlobalFunction( QuiverStaticDictionary, function(quiver, patterns)
    local pattern, tree, dict, i, l, newNode, char, v, w, q, fam,
          parent, currentNode, elemFam;

    if not IsQuiver(quiver) then
        Error("<quiver> must be a quiver");
    fi;

    elemFam := ElementsFamily(FamilyObj(quiver));

    if not IsDenseList(patterns) 
       or not ForAll(patterns, x -> IsIdenticalObj(FamilyObj(x), elemFam))
    then
        Error("<patterns> must be a dense list of elements in <quiver>");
    fi;

    tree := rec( labels := [],
                children := [],
                patterns := [] );

    tree.failureNode := tree;

    # Initialize for the vertices
    for v in VerticesOfQuiver(quiver) do
        Add(tree.labels, v);
        newNode := rec( label := v,
                        labels := [],
                        children := [],
                        parent := tree,
                        patterns := [] );
        Add(tree.children, newNode);
    od;

    for i in [1..Length(patterns)] do
        pattern := WalkOfPath(patterns[i]);

        # starting node is a transition on the starting vertex.
        v := SourceOfPath(patterns[i]);
        currentNode := tree.children[Position(tree.labels, v)];

        for char in pattern do
            l := Position(currentNode.labels, char);
            if l = fail then
                Add(currentNode.labels, char);
                newNode := rec( label := char,
                                labels := [],
                                children := [],
                                parent := currentNode,
                                patterns := [] );
                Add(currentNode.children, newNode);
                currentNode := newNode;
            else
                currentNode := currentNode.children[l];
            fi;
        od;
        Add(currentNode.patterns, i);
    od;
    q := [];
    List(tree.children, function(x)
        x.failureNode := tree;
        Append(q, x.children);
        return x;
    end );

    while not IsEmpty(q) do
        v := q[1];
        q := q{[2..Length(q)]};
        parent := v.parent;
        char := parent.labels[Position(parent.children, v)];
        w := parent.failureNode;
        while w <> tree 
              and Position(w.labels, char) = fail
        do
            w := w.failureNode;
        od;
        l := Position(w.labels, char);
        if l = fail then
            v.failureNode := tree.children[Position(tree.labels, TargetOfPath(char))];
        else
            v.failureNode := w.children[l];
        fi;
        if not IsEmpty(v.failureNode.patterns) then
            v.outputNode := v.failureNode;
        elif IsBound(v.failureNode.outputNode) then
            v.outputNode := v.failureNode.outputNode;
        fi;
        Append(q, v.children);
    od;
        fam := NewFamily( "QuiverStaticDictionaryFamily",
                          IsQuiverStaticDictionary );
        dict := Objectify( 
                    NewType( fam, IsQuiverStaticDictionary
                                  and IsQuiverStaticDictionaryDefaultRep ),
                    rec( tree := tree ) ); 
        SetPatterns( dict, ShallowCopy(patterns) );

    return dict;
end );
InstallOtherMethod( Search,
    "for static dictionaries for quivers in their default rep.",
    true,
    [IsQuiverStaticDictionary and IsStaticDictionaryDefaultRep,
     IsPath], 10,
    function(dict, text)
        local c, w, m, j, k, currentNode, matches, retval, tree, l;

        tree := dict!.tree;
        w := tree.children[Position(tree.labels, SourceOfPath(text))];
        text := WalkOfPath(text);
        c := 1;
        m := Length(text);
        retval := [];

        while c <= m do
            j := Position(w.labels, text[c]);
            while j <> fail do
                currentNode := w.children[j];
                matches := [];
                Append(matches,currentNode.patterns);
                while IsBound(currentNode.outputNode) do
                    currentNode := currentNode.outputNode;
                    Append(matches,currentNode.patterns);
                od;
                if not IsEmpty(matches) then
                    Add(retval, [c, matches]);
                fi;
                w := w.children[j];
                c := c + 1;
                if c <= m then
                    j := Position(w.labels, text[c]);
                else
                    j := fail;
                fi;
            od;
            if w = tree then
                c := c + 1;
                if c <= m then
                    l := Position(tree.labels, SourceOfPath(text[c]));
                    w := tree.children[l];
                fi;
            else
                w := w.failureNode;
            fi;
         od;

         return retval;
    end );
InstallMethod( PrefixSearch,
    "for quiver static dictionaries and paths",
    true,
    [ IsQuiverStaticDictionary and IsStaticDictionaryDefaultRep, IsPath ], 0,
    function( dict, text )
        local w, c, m, retval, j, matches, currentNode;

        w := dict!.tree;
        j := Position(w.labels, SourceOfPath(text));
        text := WalkOfPath(text);
        c := 0;
        m := Length(text);
        retval := [];

        while j <> fail do
            currentNode := w.children[j];
            matches := [];
            Append(matches,currentNode.patterns);
            if not IsEmpty(matches) then
                Add(retval, [c, matches]);
            fi;
            w := w.children[j];
            c := c + 1;
            if c <= m then
                j := Position(w.labels, text[c]);
            else
                j := fail;
            fi;
        od;
        return retval;
     end );
InstallMethod( ViewObj,
    "for static dictionaries",
    true,
    [IsStaticDictionary], 0,
    function(dict)
        local l;
        Print("<static dictionary containing ");
        l := Length(Patterns(dict));
        Print(l);
        Print(" pattern");
        if l > 1 then
            Print("s");
        fi;
        Print(">");
    end );
InstallMethod( ViewObj,
    "for static dictionaries",
    true,
    [IsQuiverStaticDictionary], 0,
    function(dict)
        local l;
        Print("<static dictionary for quivers containing ");
        l := Length(Patterns(dict));
        Print(l);
        Print(" path");
        if l > 1 then
            Print("s");
        fi;
        Print(">");
    end );
