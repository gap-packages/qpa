# GAP Declarations
# This file was generated from
# $Id: dyndict.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
DeclareInfoClass( "InfoDynamicDictionary" );
DeclareCategory( "IsSuffixTree", IsObject );
DeclareRepresentation( "IsDefaultSuffixTreeRep", 
    IsComponentObjectRep,
    ["patterns", "root", "deadLength", "family"] );
DeclareCategory( "IsSuffixTreeElement", IsObject );
DeclareCategory( "IsSuffixTreeNode", IsSuffixTreeElement );
DeclareRepresentation( "IsDefaultSuffixTreeNodeRep", 
    IsComponentObjectRep,
    ["parentEdge", "label", "dynTree"] );
DeclareCategory( "IsInternalSuffixTreeNode", IsSuffixTreeNode );
DeclareRepresentation( "IsDefaultInternalSuffixTreeNodeRep",
        IsDefaultSuffixTreeNodeRep and IsComponentObjectRep,
        ["edges", "suffixLink", "marked"] );
DeclareCategory( "IsExternalSuffixTreeNode", IsSuffixTreeNode );
DeclareCategory( "IsSuffixTreeEdge", IsSuffixTreeElement );
DeclareRepresentation( "IsDefaultSuffixTreeEdgeRep",
    IsComponentObjectRep,
    [ "label", "parent", "child" ] );
DeclareGlobalFunction( "CreateSuffixTree" );
DeclareGlobalFunction( "CreateInternalSuffixTreeNode" );
DeclareGlobalFunction( "CreateExternalSuffixTreeNode" );
DeclareGlobalFunction( "CreateSuffixTreeEdge" );
DeclareOperation( "ParentNode", [IsSuffixTreeNode] );
DeclareOperation( "AddInternalSuffixTreeNode", 
    [IsSuffixTreeEdge, IsPosInt, IsInternalSuffixTreeNode] );
DeclareOperation( "AddExternalSuffixTreeNode",
    [IsInternalSuffixTreeNode, IsExternalSuffixTreeNode] );
DeclareOperation( "LabelOfSuffixTreeNode", [IsSuffixTreeNode] );
DeclareOperation( "LabelOfSuffixTreeEdge", [IsSuffixTreeEdge] );
DeclareOperation( "SuffixTreeEdgeStartingAt", 
    [IsInternalSuffixTreeNode, IsObject] );
DeclareOperation( "InsertPatternIntoSuffixTree", [IsSuffixTree, IsList] );
DeclareOperation("DeletePatternFromSuffixTree",
    [IsSuffixTree, IsPosInt]);
DeclareOperation( "AllPatternsInSequence", [IsSuffixTree, IsList] );
DeclareOperation( "SequenceInPatterns", [IsSuffixTree, IsList] );
