# GAP Implementation
# This contains various graph algorithms for testing quiver's properties
# (created A. Mroz, 21.02.2013)

DeclareProperty("IsAcyclicQuiver", IsQuiver);
InstallImmediateMethod(IsFinite, IsQuiver and HasIsAcyclicQuiver, 0, IsAcyclicQuiver);
DeclareProperty("IsConnectedQuiver", IsQuiver);
DeclareProperty("IsTreeQuiver", IsQuiver); # better name??
DeclareProperty("IsUAcyclicQuiver", IsQuiver); # better name??
DeclareProperty("IsDynkinQuiver", IsQuiver);

DeclareOperation("FullSubquiver", [IsQuiver, IsList]);
DeclareOperation("ConnectedComponentsOfQuiver", [IsQuiver]);

