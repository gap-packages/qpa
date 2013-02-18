# GAP Declarations
# This file was generated from
# $Id: quiver.gd,v 1.4 2012/02/27 12:26:34 sunnyquiver Exp $
DeclareInfoClass( "InfoQuiver" );

DeclareCategory("IsPath", IsMultiplicativeElement);
DeclareCategory("IsVertex", IsPath);
DeclareRepresentation("IsVertexRep", IsComponentObjectRep, ["vertex_name","gen_pos"]);
DeclareCategory("IsArrow", IsPath);
DeclareRepresentation("IsArrowRep", IsComponentObjectRep, ["arrow_name","gen_pos"]);
DeclareCategoryFamily("IsPath");

DeclareAttribute("SourceOfPath", IsPath);
DeclareAttribute("TargetOfPath", IsPath);
DeclareAttribute("LengthOfPath", IsPath);
DeclareAttribute("WalkOfPath", IsPath);
DeclareProperty("IsZeroPath", IsPath);

DeclareAttribute("IncomingArrowsOfVertex", IsVertex, "mutable");
DeclareAttribute("OutgoingArrowsOfVertex", IsVertex, "mutable");
DeclareAttribute("InDegreeOfVertex", IsVertex, "mutable");
DeclareAttribute("OutDegreeOfVertex", IsVertex, "mutable");
DeclareAttribute("NeighborsOfVertex", IsVertex, "mutable" );

DeclareGlobalFunction( "Path" );
DeclareCategory("IsQuiver", IsSemigroup and IsRecord);
DeclareRepresentation("IsQuiverRep", IsQuiver and IsComponentObjectRep,["pieces"]);
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
DeclareProperty("IsAcyclicQuiver", IsQuiver);
DeclareProperty("IsConnectedQuiver", IsQuiver);
DeclareProperty("IsTreeQuiver", IsQuiver); # better name??
DeclareProperty("IsUAcyclicQuiver", IsQuiver); # better name??
InstallImmediateMethod(IsFinite, IsQuiver and HasIsAcyclicQuiver, 0, IsAcyclicQuiver);
DeclareAttribute("NumberOfVertices", IsQuiver);
DeclareAttribute("NumberOfArrows", IsQuiver);
DeclareAttribute( "OrderingOfQuiver", IsQuiver);

DeclareGlobalFunction("Quiver");
DeclareGlobalFunction("UniqueQuiverName");
DeclareOperation("OrderedBy",[IsQuiver, IsQuiverOrdering]);
DeclareOperation("NextPath", [IsQuiver, IsObject]); 
DeclareOperation( "\[\]", [ IsQuiverEnumerator, IsPosInt ]);

DeclareOperation( "QuiverContainingPath", [ IsPath ] );
DeclareOperation( "VertexIndex", [ IsVertex ] );
DeclareOperation( "ArrowIndex", [ IsArrow ] );
DeclareOperation( "GeneratorIndex", [ IsPath ] );

DeclareAttribute( "OppositeQuiver", IsQuiver );
DeclareOperation( "OppositePath", [ IsPath ] );
DeclareAttribute( "OppositeQuiverNameMap", IsQuiver );