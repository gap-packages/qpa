DeclareCategory( "IsPoset", IsObject and IsAttributeStoringRep );
DeclareRepresentation( "IsPosetRep", IsPoset, [ "set", "binary_relation", "minimal_relations", "lessthan" ] );
DeclareAttribute( "Size", IsPoset );
DeclareOperation( "Poset", [ IsList, IsList ] );
DeclareOperation( "UnderlyingSet", [ IsPoset ] );
DeclareOperation( "PartialOrderOfPoset", [ IsPoset ] );
BindGlobal( "TheFamilyOfPosets", NewFamily( "TheFamilyOfPosets" ) );
BindGlobal( "TheTypePoset", NewType( TheFamilyOfPosets, IsPosetRep ) );