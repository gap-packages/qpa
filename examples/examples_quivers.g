q1 := Quiver(["u","v"],[["u","u","a"],["u","v","b"], 
              ["v","u","c"],["v","v","d"]]);
VerticesOfQuiver(q1);
ArrowsOfQuiver(q1);
q2 := Quiver(2,[[1,1],[2,1],[1,2]]);
ArrowsOfQuiver(q2);
VerticesOfQuiver(q2);
q3 := Quiver(2,[[1,1,"a"],[2,1,"b"],[1,2,"c"]]);
ArrowsOfQuiver(q3);
q4 := Quiver([[1,1],[2,1]]);
VerticesOfQuiver(q4);
ArrowsOfQuiver(q4);
SourceOfPath(q4.a2);
TargetOfPath(q4.a2);
#######################
q1 := Quiver(2,[[1,2]]);
IsQuiver("v1");
IsQuiver(q1);
IsAcyclicQuiver(q1); IsUAcyclicQuiver(q1); 
IsConnectedQuiver(q1); IsTreeQuiver(q1);
q2 := Quiver(["u","v"],[["u","v"],["v","u"]]);
IsAcyclicQuiver(q2); IsUAcyclicQuiver(q2); 
IsConnectedQuiver(q2); IsTreeQuiver(q2);
IsFinite(q1); IsFinite(q2);
q3 := Quiver(["u","v"],[["u","v"],["u","v"]]);
IsAcyclicQuiver(q3); IsUAcyclicQuiver(q3); 
IsConnectedQuiver(q3); IsTreeQuiver(q3); 
q4 := Quiver(2, []);
IsAcyclicQuiver(q4); IsUAcyclicQuiver(q4); 
IsConnectedQuiver(q4); IsTreeQuiver(q4);
######################
q1 := Quiver(4,[[1,4],[4,2],[3,4]]);
IsDynkinQuiver(q1);
q2 := Quiver(2,[[1,2],[1,2]]);
IsDynkinQuiver(q2);
q3 := Quiver(5,[[1,5],[2,5],[3,5],[4,5]]);
IsDynkinQuiver(q3);
#############
q1 := Quiver(["u","v"],[["u","u","a"],["u","v","b"],
              ["v","u","c"],["v","v","d"]]);
q1.a;
q1.v;
VerticesOfQuiver(q1);
ArrowsOfQuiver(q1);
AdjacencyMatrixOfQuiver(q1);
GeneratorsOfQuiver(q1);
NumberOfVertices(q1);
NumberOfArrows(q1);
OrderingOfQuiver(q1);
q1_op := OppositeQuiver(q1);
VerticesOfQuiver(q1_op);
ArrowsOfQuiver(q1_op);
##############
Q := Quiver(6, [ [1,2],[1,1],[3,2],[4,5],[4,5] ]);
VerticesOfQuiver(Q);
FullSubquiver(Q, [Q.v1, Q.v2]);
ConnectedComponentsOfQuiver(Q);
###############
q1 := Quiver(["u","v"],[["u","u","a"],["u","v","b"],
              ["v","u","c"],["v","v","d"]]);
IsPath(q1.b);
IsPath(q1.u);
IsQuiverVertex(q1.c);
IsZeroPath(q1.d);
###############
q1 := Quiver(["u","v"],[["u","u","a"],["u","v","b"],
              ["v","u","c"],["v","v","d"]]);
SourceOfPath(q1.v);
p1:=q1.a*q1.b*q1.d*q1.d;
TargetOfPath(p1);
p2:=q1.b*q1.b;
WalkOfPath(p1);
WalkOfPath(q1.a);
LengthOfPath(p1);
LengthOfPath(q1.v);
###############
q1 := Quiver(["u","v"],[["u","u","a"],["u","v","b"],
              ["v","u","c"],["v","v","d"]]);
OutgoingArrowsOfVertex(q1.u);
InDegreeOfVertex(q1.u);
NeighborsOfVertex(q1.v);