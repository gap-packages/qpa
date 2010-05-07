# GAP Implementation
# This file was generated from 
# $Id: dyntree.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallMethod( \=,
    "for splay trees",
    IsIdenticalObj,
    [IsDynamicTree, IsDynamicTree], 0,
    IS_IDENTICAL_OBJ );
InstallGlobalFunction( CreateDynamicTree, function(arg)
    local dynTree;

    if Length(arg) = 0 or Length(arg) = 1 then
        dynTree := Objectify(NewType(DynamicTreeFamily, 
                IsDynamicTree and IsDefaultDynamicTreeRep),
            rec( parent := 0, left := 0, right := 0 ) );
        if Length(arg) = 1 then
            dynTree!.item := arg[1];
        fi;
        return dynTree;
    else
        Error("must have 0 or 1 arguments");
    fi;
    end );
InstallMethod( RightRotateDynamicTree,
    "for dynamic trees in default representation",
    true,
    [IsDynamicTree and IsDefaultDynamicTreeRep], 0,
    function(y)
        local x, z;

        if y <> 0 then        
            x := y!.left;
            z := y!.parent;
            if z <> 0 then
                if z!.left = y then
                    z!.left := x;
                elif z!.right = y then
                    z!.right := x;
                fi;
            fi;
            y!.left := x!.right;
            x!.right := y;
            x!.parent := z;
            y!.parent := x;
            if y!.left <> 0 then
                y!.left!.parent := y;
            fi;
        fi;
    end );
InstallMethod( LeftRotateDynamicTree,
    "for dynamic trees in default representation",
    true,
    [IsDynamicTree and IsDefaultDynamicTreeRep], 0,
    function(y)
        local x, z;

        if y <> 0 then        
            x := y!.right;
            z := y!.parent;
            if z <> 0 then
                if z!.left = y then
                    z!.left := x;
                elif z!.right = y then
                    z!.right := x;
                fi;
            fi;
            y!.right := x!.left;
            x!.left := y;
            x!.parent := z;
            y!.parent := x;
            if y!.right <> 0 then
                y!.right!.parent := y;
            fi;
        fi;
    end );
InstallMethod( SpliceDynamicTree,
    "for dynamic trees in default representation",
    true,
    [IsDynamicTree and IsDefaultDynamicTreeRep], 0,
    function(v)
        local u, w;

        w := v!.parent;
        if w <> 0 then
            Assert(1, v <> w!.left and v <> w!.right, 
                   "The node is not a middle child");
            w!.left := v;
        fi;
    end );
InstallMethod( SplayDynamicTree,
    "for dynamic trees in default representation",
    true,
    [ IsDynamicTree and IsDefaultDynamicTreeRep ], 0,
    function( x )
        local Left, Right, Splay, Parent, Grandparent, 
              p, w;

        Left := function(t)
            if t = 0 then
                return 0;
            else
                return t!.left;
            fi;
        end;
        Right := function(t)
            if t = 0 then
                return 0;
            else
                return t!.right;
            fi;
        end;
        Parent := function(t)
            if t = 0 then
                return 0;
            else
                return t!.parent;
            fi;
        end;
        Grandparent := t -> Parent(Parent(t));
        Splay := function(t)
            local p, g;

            p := Parent(t);
            while t = Left(p) or t = Right(p) do
                g := Grandparent(t);
                if t = Left(p) then
                    if p = Left(g) then 
                        RightRotateDynamicTree(g);
                    fi;
                    RightRotateDynamicTree(Parent(t));
                elif t = Right(p) then
                    if p = Right(g) then
                        LeftRotateDynamicTree(g);
                    fi;
                    LeftRotateDynamicTree(Parent(t));
                fi;
                p := Parent(t);
            od;
        end;

        w := x;
        while Parent(w) <> 0 do
            Splay(w);
            w := Parent(w);
        od;
        w := x;
        while Parent(w) <> 0 do
            SpliceDynamicTree(w);
            w := Parent(w);
        od;
        Splay(x);
    end );
InstallMethod( RootOfDynamicTree,
    "for dynamic trees in default representation",
    true,
    [IsDynamicTree and IsDefaultDynamicTreeRep], 0,
    function(v)
        local w;

        SplayDynamicTree(v);
        w := v;
        while w!.right <> 0 do
            w := w!.right;
        od;
        SplayDynamicTree(w);
        return w;
    end );
InstallMethod( LinkDynamicTrees,
    "for dynamic trees in default representation",
    true,
    [IsDynamicTree and IsDefaultDynamicTreeRep,
     IsDynamicTree and IsDefaultDynamicTreeRep], 0,
    function(v, w)
        SplayDynamicTree(v);
        SplayDynamicTree(w);
        v!.parent := w;
    end );
InstallMethod( CutDynamicTree,
    "for dynamic trees in default representation",
    true,
    [IsDynamicTree and IsDefaultDynamicTreeRep], 0,
    function(v)
        SplayDynamicTree(v);
        if v!.right <> 0 then
            v!.right!.parent := 0;
            v!.right := 0;
        fi;
    end );
