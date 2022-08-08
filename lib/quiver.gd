# GAP Declarations
# This file was generated from
# $Id: quiver.gd,v 1.4 2012/02/27 12:26:34 sunnyquiver Exp $
DeclareInfoClass( "InfoQuiver" );

DeclareCategory("IsPath", IsMultiplicativeElement);
DeclareCategory("IsQuiverVertex", IsPath);
DeclareRepresentation("IsQuiverVertexRep", IsComponentObjectRep, ["vertex_name","gen_pos"]);
DeclareCategory("IsArrow", IsPath);
DeclareRepresentation("IsArrowRep", IsComponentObjectRep, ["arrow_name","gen_pos"]);
DeclareCategoryFamily("IsPath");

DeclareAttribute("SourceOfPath", IsPath);
DeclareAttribute("TargetOfPath", IsPath);
DeclareAttribute("LengthOfPath", IsPath);
DeclareAttribute("WalkOfPath", IsPath);
DeclareProperty("IsZeroPath", IsPath);
InstallTrueMethod( IsPath, IsZeroPath );

DeclareAttribute("IncomingArrowsOfVertex", IsQuiverVertex, "mutable");
DeclareAttribute("OutgoingArrowsOfVertex", IsQuiverVertex, "mutable");
DeclareAttribute("InDegreeOfVertex", IsQuiverVertex, "mutable");
DeclareAttribute("OutDegreeOfVertex", IsQuiverVertex, "mutable");
DeclareAttribute("NeighborsOfVertex", IsQuiverVertex, "mutable" );

DeclareGlobalFunction( "Path" );
DeclareCategory("IsQuiver", IsSemigroup and IsRecord);
DeclareRepresentation("IsQuiverRep", IsQuiver and IsComponentObjectRep,["pieces"]);
DeclareCategory("IsQuiverSA", IsQuiver);
DeclareRepresentation("IsQuiverIteratorRep", IsComponentObjectRep, ["quiver","position"]);

# DeclareRepresentation( "IsQuiverEnumerator",
#    IsDomainEnumerator and IsComponentObjectRep,
#    [ "quiver" ] ) ;
#
DeclareRepresentation("IsQuiverEnumerator", IsComponentObjectRep, ["quiver"]);
DeclareAttribute("VerticesOfQuiver", IsQuiver);
DeclareAttribute("ArrowsOfQuiver", IsQuiver);
DeclareAttribute("AdjacencyMatrixOfQuiver", IsQuiver);
DeclareSynonymAttr( "GeneratorsOfQuiver", GeneratorsOfMagma);


DeclareAttribute("NumberOfVertices", IsQuiver);
DeclareAttribute("NumberOfArrows", IsQuiver);
DeclareAttribute( "OrderingOfQuiver", IsQuiver);

DeclareGlobalFunction("Quiver");
DeclareGlobalFunction("UniqueQuiverName");
DeclareOperation("OrderedBy",[IsQuiver, IsQuiverOrdering]);
DeclareOperation("NextPath", [IsQuiver, IsObject]); 
DeclareOperation( "\[\]", [ IsQuiverEnumerator, IsPosInt ]);

DeclareOperation( "QuiverContainingPath", [ IsPath ] );
DeclareOperation( "VertexIndex", [ IsQuiverVertex ] );
DeclareOperation( "ArrowIndex", [ IsArrow ] );
DeclareOperation( "GeneratorIndex", [ IsPath ] );

DeclareAttribute( "OppositeQuiver", IsQuiver );
DeclareOperation( "OppositePath", [ IsPath ] );
DeclareAttribute( "OppositeQuiverNameMap", IsQuiver );

DeclareOperation( "SeparatedQuiver", [ IsQuiver ] );
DeclareOperation( "DynkinQuiverAn", [ IS_INT, IsList ] );
DeclareOperation( "DynkinQuiverEn", [ IS_INT, IsList ] );
DeclareOperation( "DynkinQuiverDn", [ IS_INT, IsList ] );
DeclareOperation( "DynkinQuiver", [ IsString, IS_INT, IsList ] );
DeclareAttribute( "DoubleQuiver", IsQuiver );
