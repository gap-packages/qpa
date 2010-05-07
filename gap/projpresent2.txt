# Projective presentation of RightModuleOverPathAlgebra

# Example for testing
q := Quiver(3,[[1,2,"a"],[1,2,"b"],[2,2,"c"],[2,3,"d"],[3,1,"e"]]);
K := Rationals;
pa := PathAlgebra(K,q);
M1:=RightModuleOverPathAlgebra(pa,[["a",[[1,0],[2,4],[3,5]]],["b",[[0,1],[0,1],[4,3]]],["c",[[0,1],[0,0]]],["d",[[0,0],[0,0]]],["e",[[2,3,0],[7,5,6]]]]);

DeclareOperation("ProjectivePresentationOfRightModuleOverPathAlgebra",[IsPathAlgebraMatModule, IsQuiver]);

InstallMethod(ProjectivePresentationOfRightModuleOverPathAlgebra, "for projective presentation of a right module over a path algebra", [IsPathAlgebraMatModule, IsQuiver], 0,

function(M,q)

	local B, num_v, d, i, v, num_basis, cnt, j, gens, pa, P0_v, P0, a, num_a, e, l, sum_d, sum_ae, m_ij_a, alpha, w, w_ij, zerolist, h, W, P1;

	B:=Basis(M);
	pa:=RightActingAlgebra(M);
	gens:=GeneratorsOfAlgebra(pa);
	num_v:=Length(VerticesOfQuiver(q));
	num_a:=Length(ArrowsOfQuiver(q));

	# list of vertices v1,...,vn of the quiver
	v:=[];
	for i in [1..num_v] do
		Add(v,gens[i]);
	od;

	# list of arrows a1,...,am of the quiver
	a:=[];
	for i in [1..num_a] do
		Add(a,gens[i+num_v]);
	od;

	# list of dimensions of M^v1,...,M^vn
	num_basis:=Length(B);
	d:=[];
	for i in [1..num_v] do
		cnt:=0;
		for j in [1..num_basis] do
			if B[j]^v[i] <> Zero(M) then
				cnt:=cnt+1;
			fi;
		od;
		Add(d,cnt);
	od;

	# e_ij
	e:=[];
	for i in [1..num_v] do
		e[i]:=[];
		for j in [1..d[i]] do
			e[i][j]:=[];
			sum_d:=0;
			for l in [1..i-1] do
				sum_d:=sum_d+d[l];
			od;
			for l in [1..sum_d+j-1] do
				e[i][j][l]:=0;
			od;
			e[i][j][sum_d+j]:=v[i];
			for l in [sum_d+j+1..num_basis] do
				e[i][j][l]:=0;
			od;
		od;
	od;
#	Print(e);
	
	# w_ij?
	zerolist:=[];
	for i in [1..num_basis] do
		Add(zerolist,Zero(pa));
	od;
	w:=[];	
	for h in [1..num_a] do
		for i in [1..num_v] do
			for j in [1..d[i]] do
				sum_d:=0;
				for l in [1..i-1] do
					sum_d:=sum_d+d[l];
				od;
				m_ij_a:= B[sum_d+j]^a[h];
#				Print(m_ij_a);
				alpha:= Coefficients(B,m_ij_a);
#				Print(alpha);
				l:=1;
				while (l <= num_v) and (m_ij_a^v[l] = Zero(M)) do
					l:=l+1;
				od;
				if l>3 then l:=3; fi;
#				Print(alpha*v[l]);
				w_ij:= e[i][j]*a[h]-alpha*v[l];
				if w_ij <> zerolist then
					Add(w,w_ij);
				fi;
			od;
		od;
	od;
	
	# create P0
	P0_v:=[];
	for i in [1..num_v] do
		for j in [1..d[i]] do
			Add(P0_v,v[i]);
		od;
	od;
	P0:=RightProjectiveModule(pa,P0_v);
	
	# create P1
	W:=[];
	for i in [1..Length(w)] do
		Add(W,Vectorize(P0,w[i]));
	od;
	
	return [W,P0];
	

end);
