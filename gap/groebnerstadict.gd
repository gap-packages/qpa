# GAP Declarations
# This file was generated from
# $Id: stadict.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
DeclareInfoClass( "InfoStaticDictionary" );
DeclareCategory("IsStaticDictionary", IsObject);
DeclareRepresentation("IsStaticDictionaryDefaultRep", 
                      IsAttributeStoringRep, ["tree"]);
DeclareAttribute( "Patterns", IsStaticDictionary );
DeclareGlobalFunction("StaticDictionary");
DeclareOperation( "Search", [IsStaticDictionary, IsDenseList] );
DeclareCategory("IsQuiverStaticDictionary", IsStaticDictionary);
DeclareRepresentation("IsQuiverStaticDictionaryDefaultRep", 
                      IsStaticDictionaryDefaultRep, ["tsorted"]);
DeclareProperty( "IsFiniteDifference", IsQuiverStaticDictionary );
DeclareAttribute( "DifferenceSize", IsQuiverStaticDictionary );
DeclareAttribute( "DifferenceWords", IsQuiverStaticDictionary );
DeclareGlobalFunction("QuiverStaticDictionary");
DeclareOperation( "PrefixSearch", [IsQuiverStaticDictionary, IsObject] );
