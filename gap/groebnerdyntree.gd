# GAP Declarations
# This file was generated from
# $Id: dyntree.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
DeclareCategory( "IsDynamicTree", IsObject);
DeclareRepresentation( "IsDefaultDynamicTreeRep",
    IsComponentObjectRep,
    [ "item", "parent", "left", "right" ] );
BindGlobal( "DynamicTreeFamily", 
    NewFamily("DynamicTreeFamily", IsDynamicTree) );
DeclareGlobalFunction( "CreateDynamicTree" );
DeclareOperation("RightRotateDynamicTree", [IsDynamicTree]);
DeclareOperation("LeftRotateDynamicTree", [IsDynamicTree]);
DeclareOperation("SpliceDynamicTree", [IsDynamicTree]);
DeclareOperation("SplayDynamicTree", [IsDynamicTree]);
DeclareOperation( "RootOfDynamicTree", [IsDynamicTree] );
DeclareOperation( "LinkDynamicTrees", [IsDynamicTree, IsDynamicTree] );
DeclareOperation( "CutDynamicTree", [IsDynamicTree] );
