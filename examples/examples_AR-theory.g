A := KroneckerAlgebra(GF(4),2);       
S := SimpleModules(A)[1];             
ass := AlmostSplitSequence(S);   
DecomposeModule(Range(ass[1]));
PredecessorsOfModule(S,5);   
A:=NakayamaAlgebra([5,4,3,2,1],GF(4));
S := SimpleModules(A)[1];             
PredecessorsOfModule(S,5);    
