# BG 14, 9
# BG 13 dokl. przetestowany
# BG 12, potest, DOP2 usuniety! testy, ok, usuniete Has_Ad_zero_col
# DOPWD, DOSWD DOP < DOP2 i wyszuk direct pred trwa prawie nic (ulamki, tylko 2 razy sekunda z kawalkiem, ale w por. do kilkadziesiat sek. szukanie pred)  na BG 7
# BG 7 wklepane, sprawdzone dimy kilka
# BG 5 rzetestowane, tez DOP < DOP2 i wyszuk direct pred trwa prawie nic
# DegOrderDirectPredecessors, DegOrderDirectSuccessors, zrobione, potestowane!
# testcodimconj zrobiona, 
# DegOrderPredecessors2, potestowane na wszystkich, DegOrderPredecessors lepszy!!!
# DegOrderPredecessors2, potestowane na wszystkich, DegOrderPredecessors lepszy!!!
# DegOrderPredecessors, DegOrderSuccessors,  reszta algebr, przetestowane!!!

# GAP Implementation
# This contains various tools for computing degenerations of modules in finite type
# (created A. Mroz, 01.07.2013)

 

######################################################
##
#F ARQuiverNumerical( <ind> , <proj>, <list> )
#F ARQuiverNumerical( <name> )
#F ARQuiverNumerical( <name>, <param1> )
#F ARQuiverNumerical( <name>, <param1>, <param2> )
##
## This function "initializes" Auslander-Reiten quiver and 
## performs all necessary preliminary computations concerning
## mainly determining the matrix of dimensions of all Hom-spaces between
## indecomposables.
## It returns <AR>, an object from the category IsARQuiverNumerical
## which represents Auslander-Reiten quiver and additionally contains all data 
## necessary for further computations in this subpackage:
## - attributes: NumberOfIndecomposables - number of indecomposable modules in our category
##               NumberOfProjectives - number of indecomposable projective modules in our category
## - components (access by <AR>!.<component>):
##   ARdesc = numerical description of AR quiver (= <list>)
##   DimHomMat = matrix [dim Hom (i,j)] (=> rows 1..p contain dim. vectors of all indecomposables)
##   Simples = list of numbers of simple modules
##
## Arguments: <ind> - number of indecomposable modules in our category
##            <proj> - number of indecomposable projective modules in our category
##            <list> - list of lists containing description of meshes in AR quiver
##            defined as follows: <list>[i] = description of mesh ending in vertex
##            (indec. mod.) no. i having the shape
##            [a1,...,an,t]
##             where a1,...,an = numbers of direct predecessors of i in AR quiver
##             t = number of tau i, or 0 if tau i does not exist (<=> i is projective).
##             In particular if i is projective <list>[i]=[a1,...,an,0]
##             where a1,...,an are indec. summands of rad i. 
##   REMARK: we ALWAYS assume that indecomposables with numbers
##           1..<proj> are projectives and the only projectives
##           (further dimension vectors are interpreted according
##            to this order of projectives!) 
## Alternative arguments:
##            <name> = string with the name of predefined AR quiver - see description of 
##                     PredefARQuivers (put "what" as <name> to obtain the description of 
##                     currently available quivers).
##            <param1> = (optional) parameter for <name>
##            <param2> = (optional) second parameter for <name>
##        

InstallGlobalFunction( ARQuiverNumerical,
  function( arg )
    local AR, fam, T, Ct, i, j, k, found, LL, n, p, t, data, simples,
          col_sum, are_proj;
    
    if (Length(arg) = 3) and (IsPosInt(arg[1])) and (IsPosInt(arg[2]))
       and (IsList(arg[3]))  then
      data := arg;      
    elif (Length(arg) = 3) and IsString(arg[1]) then
      data := PredefARQuivers(arg[1], arg[2], arg[3]);
    elif (Length(arg) = 2) and IsString(arg[1]) then
      data := PredefARQuivers(arg[1], arg[2]);
    elif (Length(arg) = 1) and IsString(arg[1]) then
      data := PredefARQuivers(arg[1]); 
      if data = false then # case when the user only asks about what predefined
                           # algebras are available (parameter "what")
        return false;
      fi;        
    else      
      Error("wrong number or type of arguments!");
    fi;
    
    n := data[1]; # number of indecomposables
    p := data[2]; # number of projectives
    LL := data[3]; # description of AR quiver (see a note before function)
      
  	# Check if LL[1],...,LL[p] define projectives  
  	for i in [1..p] do
	    if LL[i][Length(LL[i])] <> 0 then
	      Error("wrong AR description, indecomposables with numbers 1,...,",p," should be projectives!");
      fi;
	  od;
	
    # Computing Auslander-Reiten matrix T
    T := NullMat(n,n);
    for i in [1..n] do 
      T[i][i] := 1;
    od;
    for i in [1..Length(LL)] do
      for j in [1..Length(LL[i])-1] do
	      T[i][LL[i][j]] := -1;
	    od;
	    t := LL[i][Length(LL[i])];
	    if t <> 0 then
	      T[i][t]:=T[i][t]+1;
	    fi;  
    od; 
	
	# for i in [1..Length(T)] do
	  # for j in [1..Length(T[1])] do
	    # if T[i][j]=-1 then Print("x"); elif T[i][j]=0 then Print(" "); else Print(T[i][j]); fi;
		
	  # od;Print("|",i,"\n");
	# od;
      
    # Computing matrix Ct = [dim Hom (i,j)] = T^{-1}^tr
    Ct := TransposedMat(Inverse(T));
    if not ForAll(Ct, s -> ForAll(s, l -> (IsPosInt(l) or (l=0)))) then
      Error("wrong AR description, Ct is not a noneg. int. matrix!");
    fi;
      
      
    # Determining positions of simples
    simples := [];
    for i in [1..p] do
      for j in [1..n] do
        if Ct[i][j] = 1 then 
          found := true;
          for k in [1..p] do
            if (k <> i) and (Ct[k][j] <> 0) then
              found := false;
            fi;
          od;          
          if found then
            Append(simples, [j]);
          fi;
        fi;
      od;
    od;
    if Length(simples) <> p then
      Error("wrong AR description - cannot find simples!");
    fi;      
      
	# #### no longer current ###
    # Checking if Ad := [p+1..n,1..n]-submatrix of Ct
    # has zero column (if yes, we change further algorithms a little!)
    # has_Ad_zero_col := false;
    # i := 1;
    # while (i <= n) and (not has_Ad_zero_col) do
      # col_sum := 0;
      # for j in [p+1..n] do
        # col_sum := col_sum + Ct[j][i];
      # od;
      # if col_sum = 0 then
        # has_Ad_zero_col := true;
      # fi;
      # i := i + 1;      
    # od;
    
    
    
    fam := NewFamily( "ARQuiverNumericalFamily", IsARQuiverNumerical );
    AR := Objectify( NewType( fam, IsARQuiverNumerical and IsAttributeStoringRep),
                     rec( ARdesc := data[3],
                          DimHomMat := Ct,
                          Simples := simples ) );
    SetNumberOfIndecomposables(AR, data[1]);
    SetNumberOfProjectives(AR, data[2]);                       
    
    return AR;

  end
); # ARQuiverNumerical



######################################################
##
#F HowManyTimes( <aa>, <bb> )
##
## Auxiliary local function. Determines how many times vector aa
## fits in vector bb, i.e. it returns max l s.t. l*aa<=bb

HowManyTimes := function(aa, bb)
  local i, ile, q, l;
   l := 1;
   while (l <= Length(aa)) and (aa[l] = 0) do
     l := l + 1;
   od;
   if l > Length(aa) then # for correct data it is impossible!
     Error("HowManyTimes: vector aa is zero?!");
   fi;
   ile := QuoInt(bb[l],aa[l]);
   for i in [l+1..Length(bb)] do
     if aa[i] <> 0 then
       q := QuoInt(bb[i], aa[i]);
       if q < ile then
         ile := q;
       fi;
     fi;
   od;
   return ile;
end; # HowManyTimes






######################################################
##
#F NextNonsimple( <simples>,  <k> )
##
## Auxiliary local function. Determines the least
## number l > k s.t. l is a number of nonsimple
## (i.e. not in the list <simples>)

NextNonsimple := function(ss, k)
  local l;
  l := k + 1;
  while l in ss do
    l := l + 1;
  od;
  return l;  
end; # NextNonsimple





######################################################
##
#O ModulesOfDimVect( <AR> , <which>)
##
## This function generates all modules (= multiplicity vectors)
## with dimension vector equal to vv, where
## vv = dimension vector of module <which>.
## <which> can be number of indec. module or dimension vector.
## (i.e. nonneg. integer solutions of Ax=vv, where A is
##  [1..p,1..n]-submatrix of <AR>!.DimHomMat; we use here the fact
##  that simple modules form an "identity submatrix" in A -
##  this allows to use efficient algorithm, much better than "naive").
##

InstallMethod( ModulesOfDimVect,
               
  [ IsARQuiverNumerical, IsObject ],
  function( AR, which )
    local solutions, vv, y, collected, Generate,
          Ct, simples, p, n, i, j, ttime, colsofA, vec;
    
    # Auxiliary recursive function generating all solutions
    Generate := function(yy, current, lcollected)
      local ii, howmany, veccurrent, next, yyy, jj, kk, res;
      
       veccurrent := colsofA[current];       
       howmany := HowManyTimes(veccurrent, vv - lcollected);
 
       next := NextNonsimple(simples, current);
       
       for ii in [0..howmany] do
         
         yy[current] := ii;  # generating all possible values at position "current"
 
         if next > n then # now we have whole vector generated
            yyy := ShallowCopy(yy);
            for jj in [1..p] do
              res := 0;
              for kk in [1..n] do
                res := res + (Ct[jj][kk] * yyy[kk]);
              od;
              yyy[simples[jj]] := vv[jj] - res;
            od;
            Append(solutions, [ShallowCopy(yyy)]);
         fi;  
 
         if next <= n then
           Generate(yy, next, lcollected + (ii * veccurrent));
         fi;
         
       od;
      
    end; # Generate local 
    
    # Read data from AR
    n := NumberOfIndecomposables(AR);
    p := NumberOfProjectives(AR);
    Ct := AR!.DimHomMat;
    simples := AR!.Simples;
    
    ttime := Runtime();
    
    # Prepare list of columns of "matrix Ad"
    # to speed up further computations
    colsofA := [];
    for i in [1..n] do
      vec := [];
      for j in [1..p] do
        Append(vec, [Ct[j][i]]);
      od;
      Append(colsofA, [ShallowCopy(vec)]);
    od;
    
    # Prepare vector vv
    if IsPosInt(which) and (which in [1..n]) then # number of indecomposable
      vv := colsofA[which];
    elif IsList(which) and (Length(which) = p) then # dimension vector
      vv := which;
    else
      Error("second argument should be a number of indecomposable or multiplicity vector!");
    fi;
    
    solutions := [];
    y := NullMat(1, n)[1];
    collected := NullMat(1, p)[1];
    Generate(y, NextNonsimple(simples, 0), collected);
   
    #Print("Time: ", Float((Runtime()-ttime)/1000), "\n");
    
    return solutions;
  
 end
); # ModulesOfDimVect



######################################################
##
#O DegOrderPredecessors( <AR> , <which>)
##
## This function generates all modules (= multiplicity vectors)
## which are the predecessors of module <which> in deg order,
## that is all M's such that M degenerates to <which> (M < <which>).
## <which> can be a number of indecomposable or the multiplicity vector.
## 

InstallMethod( DegOrderPredecessors,
               
  [ IsARQuiverNumerical, IsObject ],
  function( AR, which )
    local solutions, vv, vd, y, collected, Generate,
          Ct, simples, p, n, i, j, ttime, colsofA, vec, 
		      pos_of_which, found, vecwhich;
    
    # Auxiliary recursive function generating all solutions
    Generate := function(yy, current, lcollected)
      local ii, howmany, veccurrent, next, yyy, jj, kk, res, isleq;
      
       veccurrent := colsofA[current];       
       howmany := HowManyTimes(veccurrent, vv - lcollected);
 
       next := NextNonsimple(simples, current);
       
       for ii in [0..howmany] do
         
         yy[current] := ii;  # generating all possible values at position "current"
 
         if next > n then # now we have whole vector generated
            yyy := ShallowCopy(yy);
            for jj in [1..p] do
              res := 0;
              for kk in [1..n] do
                res := res + (Ct[jj][kk] * yyy[kk]);
              od;
              yyy[simples[jj]] := vv[jj] - res;
            od;
			
		   	    kk := p + 1;
			      isleq := true;
			      while isleq and (kk <= n) do # here we check if  yyy DegOrderLEQ which
              if Ct[kk] * yyy > vd[kk-p] then
                isleq := false;
              fi;
              kk := kk + 1;
            od;
		
		        if isleq then 
                Append(solutions, [ShallowCopy(yyy)]);
			      fi;
         fi;  
 
         if next <= n then
           Generate(yy, next, lcollected + (ii * veccurrent));
         fi;
         
       od;
      
    end; # Generate local 
    
    # Read data from AR
    n := NumberOfIndecomposables(AR);
    p := NumberOfProjectives(AR);
    Ct := AR!.DimHomMat;
    simples := AR!.Simples;
    
    ttime := Runtime();
    
    # Prepare list of columns of "matrix Ad"
    # to speed up further computations
    colsofA := [];
    for i in [1..n] do
      vec := [];
      for j in [1..p] do
        Append(vec, [Ct[j][i]]);
      od;
      Append(colsofA, [ShallowCopy(vec)]);
    od;
    
    # Prepare vector vv = dim vect (which)
    vv := DimensionVector(AR, which); 
	
	  # Prepare the vector vd = [ [X,which] ] for all indec. nonproj. X
	  if IsPosInt(which) then
	    vd := [];
      for j in [p+1..n] do
        Append(vd, [Ct[j][which]]);
      od;
	  else # which is a multiplicity vector
	    vd := Ct{[p+1..n]} * which;
	  fi;
    
    solutions := [];
    y := NullMat(1, n)[1];
    collected := NullMat(1, p)[1];
    Generate(y, NextNonsimple(simples, 0), collected);
	
	  # throw out <which> from solutions
	  if IsPosInt(which) then
	    vecwhich := NullMat(1, n)[1];
	    vecwhich[which] := 1;
	  else
	    vecwhich := which;
	  fi;
    pos_of_which := Position(solutions, vecwhich);
	  solutions := Concatenation(solutions{[1..pos_of_which-1]}, 
		  					  solutions{[pos_of_which+1..Length(solutions)]});
   
    #Print("Time (=, <=): ", Float((Runtime()-ttime)/1000), "\n");
    
    return solutions;
  
 end
); # DegOrderPredecessors




######################################################
##
#O DegOrderDirectPredecessors( <AR> , <which>)
##
## This function generates all modules (= multiplicity vectors)
## which are the DIRECT predecessors of module <which> in deg order,
## that is all M's such that M degenerates to <which> (M < <which>)
## and there is no M' s.t. M < M' < <which>
##
## <which> can be a number of indecomposable or the multiplicity vector.
## (Naive method of choosing direct ones, but enough fast.)

InstallMethod( DegOrderDirectPredecessors,
               
  [ IsARQuiverNumerical, IsObject ],
  function( AR, which )
    local preds, factorizes, i, j, r, ttime, results;
    
    preds := DegOrderPredecessors(AR, which);
    
    ttime := Runtime();
    r := Length(preds);
    results := [];
    for i in [1..r] do # check if preds[i] factorizes by sth
      factorizes := false;
      j := 1;
      while (not factorizes) and (j <= r) do
        if (i <> j) and (DegOrderLEQNC(AR, preds[i], preds[j])) then
          factorizes := true;
        fi;
        j := j + 1;
      od;
      if not factorizes then
        Append(results, [preds[i]]);
      fi;
    od;
    
    #Print("Time (Direct): ", Float((Runtime()-ttime)/1000), "\n");
    
    return results;
  
 end
); # DegOrderDirectPredecessors



######################################################
##
#O DegOrderPredecessorsWithDirect( <AR> , <which>)
##
## This function returns a pair (2-element list) [p, dp] where
## p = result of function DegOrderPredecessors;
## dp = result of function DegOrderDirectPredecessors;
##
## The function generates predecessors only once,
## so the runtime is exactly the same as  DegOrderDirectPredecessors
##
## <which> can be a number of indecomposable or the multiplicity vector.
## Naive method!!

InstallMethod( DegOrderPredecessorsWithDirect,
               
  [ IsARQuiverNumerical, IsObject ],
  function( AR, which )
    local preds, factorizes, i, j, r, ttime, results;
    
    preds := DegOrderPredecessors(AR, which);
    
    ttime := Runtime();
    r := Length(preds);
    results := [];
    for i in [1..r] do # check if preds[i] factorizes by sth
      factorizes := false;
      j := 1;
      while (not factorizes) and (j <= r) do
        if (i <> j) and (DegOrderLEQNC(AR, preds[i], preds[j])) then
          factorizes := true;
        fi;
        j := j + 1;
      od;
      if not factorizes then
        Append(results, [preds[i]]);
      fi;
    od;
    
    #Print("Time (Direct): ", Float((Runtime()-ttime)/1000), "\n");
    
    return [preds, results];
  
 end
); # DegOrderPredecessorsWithDirect










######################################################
##
#O DegOrderSuccessors( <AR> , <which>)
##
## This function generates all modules (= multiplicity vectors)
## which are the successors of module <which> in deg order,
## that is all M's such that <which> degenerates to M (<which> < M).
## <which> can be a number of indecomposable or the multiplicity vector.
## 

InstallMethod( DegOrderSuccessors,
               
  [ IsARQuiverNumerical, IsObject ],
  function( AR, which )
    local solutions, vv, vd, y, collected, Generate,
          Ct, simples, p, n, i, j, ttime, colsofA, vec, 
		  pos_of_which, found, vecwhich;
    
    # Auxiliary recursive function generating all solutions
    Generate := function(yy, current, lcollected)
      local ii, howmany, veccurrent, next, yyy, jj, kk, res, isgeq;
      
       veccurrent := colsofA[current];       
       howmany := HowManyTimes(veccurrent, vv - lcollected);
 
       next := NextNonsimple(simples, current);
       
       for ii in [0..howmany] do
         
         yy[current] := ii;  # generating all possible values at position "current"
 
         if next > n then # now we have whole vector generated
            yyy := ShallowCopy(yy);
            for jj in [1..p] do
              res := 0;
              for kk in [1..n] do
                res := res + (Ct[jj][kk] * yyy[kk]);
              od;
              yyy[simples[jj]] := vv[jj] - res;
            od;
			
		  	    kk := p + 1;
			      isgeq := true;
			      while isgeq and (kk <= n) do # here we check if  which DegOrderLEQ yyy
              if Ct[kk] * yyy < vd[kk-p] then
                isgeq := false;
              fi;
              kk := kk + 1;
            od;
		
		        if isgeq then 
              Append(solutions, [ShallowCopy(yyy)]);
			      fi;
         fi;  
 
         if next <= n then
           Generate(yy, next, lcollected + (ii * veccurrent));
         fi;
         
       od;
      
    end; # Generate local 
    
    # Read data from AR
    n := NumberOfIndecomposables(AR);
    p := NumberOfProjectives(AR);
    Ct := AR!.DimHomMat;
    simples := AR!.Simples;
    
    ttime := Runtime();
    
    # Prepare list of columns of "matrix Ad"
    # to speed up further computations
    colsofA := [];
    for i in [1..n] do
      vec := [];
      for j in [1..p] do
        Append(vec, [Ct[j][i]]);
      od;
      Append(colsofA, [ShallowCopy(vec)]);
    od;
    
    # Prepare vector vv = dim vect (which)
    vv := DimensionVector(AR, which); 
	
	  # Prepare the vector vd = [ [X,which] ] for all indec. nonproj. X
	  if IsPosInt(which) then
	    vd := [];
      for j in [p+1..n] do
        Append(vd, [Ct[j][which]]);
      od;
	  else # which is a multiplicity vector
	    vd := Ct{[p+1..n]} * which;
	  fi;
    
    solutions := [];
    y := NullMat(1, n)[1];
    collected := NullMat(1, p)[1];
    Generate(y, NextNonsimple(simples, 0), collected);
	
	  # throw out <which> from solutions
	  if IsPosInt(which) then
	    vecwhich := NullMat(1, n)[1];
	    vecwhich[which] := 1;
	  else
	    vecwhich := which;
	  fi;
    pos_of_which := Position(solutions, vecwhich);
	  solutions := Concatenation(solutions{[1..pos_of_which-1]}, 
							  solutions{[pos_of_which+1..Length(solutions)]});
   
    #Print("Time (DegOrderSuccessors =, >=): ", Float((Runtime()-ttime)/1000), "\n");
    
    return solutions;
  
 end
); # DegOrderSuccessors




######################################################
##
#O DegOrderDirectSuccessors( <AR> , <which>)
##
## This function generates all modules (= multiplicity vectors)
## which are the DIRECT successors of module <which> in deg order,
## that is all M's such that <which> degenerates to M (<which> < M)
## and there is no M' s.t. <which> < M' < M
##
## <which> can be a number of indecomposable or the multiplicity vector.
## Naive method!!

InstallMethod( DegOrderDirectSuccessors,
               
  [ IsARQuiverNumerical, IsObject ],
  function( AR, which )
    local succs, factorizes, i, j, r, ttime, results;
    
    succs := DegOrderSuccessors(AR, which);
    
    ttime := Runtime();
    r := Length(succs);
    results := [];
    for i in [1..r] do # check if succs[i] factorizes by sth
      factorizes := false;
      j := 1;
      while (not factorizes) and (j <= r) do
        if (i <> j) and (DegOrderLEQNC(AR, succs[j], succs[i])) then
          factorizes := true;
        fi;
        j := j + 1;
      od;
      if not factorizes then
        Append(results, [succs[i]]);
      fi;
    od;
    
    #Print("Time (Direct): ", Float((Runtime()-ttime)/1000), "\n");
    
    return results;
  
 end
); # DegOrderDirectSuccessors


######################################################
##
#O DegOrderSuccessorsWithDirect( <AR> , <which>)
##
## This function returns a pair (2-element list) [s, ds] where
## s = result of function DegOrderSuccessors;
## ds = result of function DegOrderDirectSuccessors;
##
## The function generates successors only once,
## so the runtime is exactly the same as  DegOrderDirectSuccessors
##
## <which> can be a number of indecomposable or the multiplicity vector.
## Naive method!!

InstallMethod( DegOrderSuccessorsWithDirect,
               
  [ IsARQuiverNumerical, IsObject ],
  function( AR, which )
    local succs, factorizes, i, j, r, ttime, results;
    
    succs := DegOrderSuccessors(AR, which);
    
    ttime := Runtime();
    r := Length(succs);
    results := [];
    for i in [1..r] do # check if succs[i] factorizes by sth
      factorizes := false;
      j := 1;
      while (not factorizes) and (j <= r) do
        if (i <> j) and (DegOrderLEQNC(AR, succs[j], succs[i])) then
          factorizes := true;
        fi;
        j := j + 1;
      od;
      if not factorizes then
        Append(results, [succs[i]]);
      fi;
    od;
    
    #Print("Time (Direct): ", Float((Runtime()-ttime)/1000), "\n");
    
    return [succs, results];
  
 end
); # DegOrderSuccessorsWithDirect





######################################################
##
#O DimensionVector( <AR> , <which>)
##
## This function returns the dimension vector of module
## indicated by <which>. The order of dimension is exactly
## the same as the order of corresponding projective indecomposables
## in ARQuiverNumerical <AR>.
## <which> can be a number of indecomposable or the multiplicity vector.
##

InstallMethod( DimensionVector,
               
  [ IsARQuiverNumerical, IsObject ],
  function( AR, which )
    local j, p, n, vec;
    
    n := NumberOfIndecomposables(AR);
    p := NumberOfProjectives(AR);
    
    if IsList(which) and (Length(which) = n) then
      return AR!.DimHomMat{[1..p]} * which;
    elif IsPosInt(which) and (which in [1..n]) then
      vec := [];
      for j in [1..p] do
        Append(vec, [AR!.DimHomMat[j][which]]);
      od;
      return vec;
    else
      Error("second parameter should be a number if indecomposable or the multiplicity vector!");    
    fi;
    
 end
); # DimensionVector



######################################################
##
#O DimHom( <AR> , <M>, <N>)
##
## This function returns the dimension of Homspace
## between modules <M> and <N>.
## <M> and <N> can be either numbers of indecomposables or the multiplicity vectors.
##

InstallMethod( DimHom,
               
  [ IsARQuiverNumerical, IsObject, IsObject ],
  function( AR, M, N )
    local i, j, n, vec, sum;
    
    n := NumberOfIndecomposables(AR);
    
    if IsPosInt(M) and IsPosInt(N) and (M in [1..n]) and (N in [1..n]) then
      # number, number
      return AR!.DimHomMat[M][N];
    elif IsPosInt(M) and (M in [1..n]) and IsList(N) and (Length(N) = n) then
      # number, vector
      sum := 0;
      for j in [1..n] do
        sum := sum + (AR!.DimHomMat[M][j] * N[j]);
      od;
      return sum;
    elif IsPosInt(N) and (N in [1..n]) and IsList(M) and (Length(M) = n) then
      # vector, number
      sum := 0;
      for j in [1..n] do
        sum := sum + (M[j] * AR!.DimHomMat[j][N]);
      od;
      return sum;
    elif IsList(M) and IsList(N) and (Length(M) = n) and (Length(N) = n) then
      # vector, vector
      sum := 0;
      for i in [1..n] do
        for j in [1..n] do
          sum := sum + (M[i] * AR!.DimHomMat[i][j] * N[j]);
        od;
      od;
      return sum;
    else    
      Error("second and third parameter should be a number if indecomposable or the multiplicity vector!");    
    fi;
    
 end
); # DimHom


######################################################
##
#O DimEnd( <AR> , <M> )
##
## This function returns the dimension of endomorphism algebra
## of module <M>.
## <M> can be either number of indecomposables or the multiplicity vector.
##

InstallMethod( DimEnd,
               
  [ IsARQuiverNumerical, IsObject ],
  function( AR, M )
    return DimHom(AR, M, M);  
 end
); # DimEnd


######################################################
##
#O OrbitCodim( <AR> , <M> , <N>)
##
## This function returns the codimension of orbits of modules <M> and <N>
## (= dim End(<N>) - dim End(<M>)).
## <M> and <N> can be either numbers of indecomposables or the multiplicity vectors.
##
## NOTE: Function does not check if it makes sense, i.e. if <M> and <N> are in the same 
## variety ( = dimension vectors coincide)!
##

InstallMethod( OrbitCodim,
               
  [ IsARQuiverNumerical, IsObject, IsObject ],
  function( AR, M, N )
    return DimEnd(AR, N) - DimEnd(AR, M);  
 end
); # OrbitCodim


######################################################
##
#O DegOrderLEQ( <AR> , <M> , <N>)
##
## This function returns true if <M> <= <N> in degeneration order
## i.e. <N> is a degeneration of <M>
## (<=> <M> <= <N> in Hom-order; we are in finite type, cf. Zwara's thm).
##
## NOTE: Function checks if it makes sense, i.e. if <M> and <N> are in the same 
## variety ( = dimension vectors coincide). If not, it returns false
## and prints warning. 
##

InstallMethod( DegOrderLEQ,
               
  [ IsARQuiverNumerical, IsObject, IsObject ],
  function( AR, M, N )
    local i, n, isleq;
    
    if DimensionVector(AR, M) = DimensionVector(AR, N) then
      isleq := true;
      n := NumberOfIndecomposables(AR);
	    i := NumberOfProjectives(AR) + 1;
	  
      if IsList(M) and IsPosInt(N) then
        while isleq and (i <= n) do
          if AR!.DimHomMat[i] * M > AR!.DimHomMat[i][N] then
            isleq := false;
          fi;
          i := i + 1;
        od;
	    elif IsList(M) and IsList(N) then
        while isleq and (i <= n) do
          if AR!.DimHomMat[i] * M > AR!.DimHomMat[i] * N then
            isleq := false;
          fi;
          i := i + 1;
        od;
      elif IsPosInt(M) and IsList(N) then
        while isleq and (i <= n) do
          if AR!.DimHomMat[i][M] > AR!.DimHomMat[i] * N then
            isleq := false;
          fi;
          i := i + 1;
        od;	  
	    elif IsPosInt(M) and IsPosInt(N) then
        while isleq and (i <= n) do
          if AR!.DimHomMat[i][M] > AR!.DimHomMat[i][N] then
            isleq := false;
          fi;
          i := i + 1;
        od;
	    fi;
	  
      return isleq;
    else
      Print("Modules are not in the same variety (their dimension vectors are different)!\n");
      return false;
    fi;      
 end
); # DegOrderLEQ

######################################################
##
#O DegOrderLEQNC( <AR> , <M> , <N>)
##
## This function returns true if <M> <= <N> in degeneration order
## i.e. <N> is a degeneration of <M>
## (<=> <M> <= <N> in Hom-order; we are in finite type, cf. Zwara's thm).
##
## NOTE: Function does not check if it makes sense, i.e. if <M> and <N> are in the same 
## variety ( = dimension vectors coincide). If not, the result doesn't make sense!
##

InstallMethod( DegOrderLEQNC,
               
  [ IsARQuiverNumerical, IsObject, IsObject ],
  function( AR, M, N )
    local i, n, isleq;
    
      isleq := true;
      n := NumberOfIndecomposables(AR);
	    i := NumberOfProjectives(AR) + 1;
	  
      if IsList(M) and IsList(N) then
        while isleq and (i <= n) do
          if AR!.DimHomMat[i] * M > AR!.DimHomMat[i] * N then
            isleq := false;
          fi;
          i := i + 1;
        od;
      elif IsList(M) and IsPosInt(N) then
        while isleq and (i <= n) do
          if AR!.DimHomMat[i] * M > AR!.DimHomMat[i][N] then
            isleq := false;
          fi;
          i := i + 1;
        od;
	    
      elif IsPosInt(M) and IsList(N) then
        while isleq and (i <= n) do
          if AR!.DimHomMat[i][M] > AR!.DimHomMat[i] * N then
            isleq := false;
          fi;
          i := i + 1;
        od;	  
	    elif IsPosInt(M) and IsPosInt(N) then
        while isleq and (i <= n) do
          if AR!.DimHomMat[i][M] > AR!.DimHomMat[i][N] then
            isleq := false;
          fi;
          i := i + 1;
        od;
	    fi;
	  
      return isleq;
    
 end
); # DegOrderLEQNC




######################################################
##
#O PrintMultiplicityVector( <mv>)
##
## This function prints the multiplicity vector <mv>
## in a more "readible" way (especially useful if
## <mv> is long and sparse). It prints a "sum" of nonzero multiplicities
## in the form "multiplicity * (no.-of-indecomposable)".
##

InstallMethod( PrintMultiplicityVector,
               
  [ IsList ],
  function( mv )
    local i, firstprinted;
    
	firstprinted := false;
	
	for i in [1..Length(mv)] do
	  if mv[i] <> 0 then
	    if not firstprinted then
		  Print(mv[i],"*(",i,")");
		  firstprinted := true;
		else
		  Print(" + ",mv[i],"*(",i,")");
		fi;
	  fi;
	od;
	if not firstprinted then 
	  Print("0");
	fi;
	Print("\n");
          
 end
); # PrintMultiplicityVector


######################################################
##
#O PrintMultiplicityVectors( <list>)
##
## This function prints all the multiplicity vectors from the <list>
## in a more "readible" way, as PrintMultiplicityVector. 
##

InstallMethod( PrintMultiplicityVectors,
               
  [ IsList ],
  function( l )
    local mv;
    
	  for mv in l do
      PrintMultiplicityVector(mv);
    od;    
          
 end
); # PrintMultiplicityVectors




InstallMethod( PrintObj,
  "for ARQuiverNumerical",
  true,
  [ IsARQuiverNumerical ], 0,
  function( G )
    
    Print( "<ARQuiverNumerical with " );
    Print( NumberOfIndecomposables(G) );
    Print( " indecomposables and " );
    Print( NumberOfProjectives(G) );
    Print( " projectives>" );
  
  end
);

InstallMethod( ViewObj,
  "for ARQuiverNumerical",
  true,
  [ IsARQuiverNumerical ], NICE_FLAGS + 1,
  function( G )

    Print( "<ARQuiverNumerical with " );
    Print( NumberOfIndecomposables(G) );
    Print( " indecomposables and " );
    Print( NumberOfProjectives(G) );
    Print( " projectives>" );
end
);  


######################################################
##
#O TestCodimConjecture( <AR>, <printbool>, <breakbool> )
##
## This function tests the A-R quiver <AR> for the following
## hypothesis:
## Given M,N such that M < N and N is indecomposable.
## M is a direct deg-predecessor of N (i.e. there is no M' s.t. M < M' < N)
## <=> Codim(M, N) = 1.
## The function tests case by case all indecomposables.
## Additional arguments:
## <printbool> = true - function prints more info (consecutive tested indesomposables).
## <printbool> = false - function prints less info (only a summary at the end).
## <breakbool> = true - function breaks when a counterexample is found.
## <breakbool> = false - function doesn't break when a counterexample is found
##               (it collects all possible counterexamples).
## Function prints a summary.

InstallMethod( TestCodimConjecture,
               
  [ IsARQuiverNumerical, IsBool, IsBool ],
  function( AR, printbool, breakbool )
    local ind, n, ttime, cexs, preddpred, vec, 
          no_of_dpred, no_of_pred, no_of_cexs, no_of_cexs_loc;
    
    ttime := Runtime();
    
    n := NumberOfIndecomposables(AR);
    no_of_dpred := 0;
	  no_of_pred := 0;
    no_of_cexs := 0;
    cexs := [];
    
    for ind in [1..n] do
      preddpred := DegOrderPredecessorsWithDirect(AR, ind); 
      no_of_pred := no_of_pred + Length(preddpred[1]);
	    no_of_dpred := no_of_dpred + Length(preddpred[2]);
      no_of_cexs_loc := 0;
      for vec in preddpred[2] do 
         if OrbitCodim(AR, vec, ind) <> 1 then
           no_of_cexs_loc := no_of_cexs_loc + 1;
           if breakbool then
             Print("Counterexample founded!\n");
             Print("Total runtime: ", Float((Runtime()-ttime)/1000), " sec.\n");
             return [vec, ind];
           fi;
           Append(cexs, [[vec, ind]]);           
         fi;         
      od;
      no_of_cexs := no_of_cexs + no_of_cexs_loc;
      if printbool then
        Print("Ind.: ", ind, "/",n,"; deg-predecessors: ", Length(preddpred[1]),
		      "; direct: ", Length(preddpred[2]),"\n");
        if no_of_cexs_loc > 0 then
          Print("\nCounterexamples: ",no_of_cexs_loc,"!\n");
        fi;
      fi;
    od;
    
    if printbool then
      Print("\n");
    fi;
    Print("Summary:\n", n," indecomposables,\n", no_of_pred, " deg-predecessors of indecomposables,\n");
    Print(no_of_dpred, " direct deg-predecessors of indecomposables.\n");
  	Print(no_of_cexs, " counterexample(s) founded.\n");
    if no_of_cexs > 0 then
      Print("Conjecture disproved!");
    else
      Print("Conjecture confirmed.");
    fi;
    Print("\nTotal runtime: ", Float((Runtime()-ttime)/1000), " sec.\n");
    
    return cexs;
  
 end
); # TestCodimConjecture






######################################################
##
#F PredefARQuivers(<name>, [<param1> [,<param2>]])
##
## An auxiliary function containing predefined AR quivers.
## it returns a list of three elements <ind>, <proj>, <list>
## as in the description of ARQuiverNumerical.
## So far we have the following AR quivers:
## <name> = 
##
## <name> = "BG" = 
##        an algebra from Bongartz-Gabriel list of maximal 
##        algebras of finite type with 2 simples. Then
##        <param1> = number of an algebra from B-G list 
##                  (so far <param1> = 1,2)
##        <param2> = additional parameter defining algebra no. <param1>
##                   <param2> is necessary for B-G algebras with numbers:
##                   1
##        
## 

InstallGlobalFunction( PredefARQuivers,
  function( arg )
    local data, i;
    
    data := [];
    
    if (Length(arg) = 1) and (arg[1] = "what") then
      Print("1) (\"BG\",i), for i=2,5,7,9,12,13,14 algebra no. i from Bongartz-Gabriel list ");
      Print("of maximal finite type with 2 simples;\n");
      Print("2) (\"BG\",i,j), for i=1, j=1,2,3, algebra no. i with parameter j ");
      Print("from\n Bongartz-Gabriel list of maximal finite type with 2 simples;\n");
      
      Print("3) (\"D4 subspace\"), path algebra of Dynkin quiver D4 with ");
      Print("subspace\n orientation of arrows;\n");
      Print("4) (\"E6 subspace\"), path algebra of Dynkin quiver E6 with ");
      Print("subspace\n orientation of arrows;\n");
  	  Print("5) (\"A3 zero\"), path algebra of equioriented Dynkin quiver A3  ");
      Print("modulo unique zero relation of length 2;\n");
	    Print("6) (\"A1,2 zero\"), path algebra of Euclidean quiver A1,2 modulo zero relation  ");
      Print("(i.e. Quiver(3,[[3,1,\"a\"],[3,2,\"b\"],[2,1,\"c\"]])/(bc) );\n");
	    Print("7) (\"A2,2 comm\"), path algebra of Euclidean quiver A2,2 modulo unique commutativity relation;\n");
	    Print("8) (\"R nilp\"), one arrow + one loop modulo nilpotency deg. 2  ");
      Print("(i.e. Quiver(2,[[2,1,\"a\"],[1,1,\"b\"]])/(b^2) );\n");
      return false;
    fi;
    
    if Length(arg) = 3 then
    
      if (arg[1] = "BG") and (arg[2] = 1) and (arg[3] = 1) then
        data[1] := 8;
        data[2] := 2;
        data[3] := [[4,0],[3,0],[1,6,4],[2,5,3],[3,8,6],[4,7,5],[5,8],[6,7]];
      fi;
    
      if (arg[1] = "BG") and (arg[2] = 1) and (arg[3] = 2) then
        data[1] := 16;
        data[2] := 2;
        data[3] := [ [4,0],[3,0],[1,6,4],[2,5,3],[3,8,6],[4,7,5],[5,10,8],[6,9,7],
				  	[7,12,10],[8,11,9],[9,14,12],[10,13,11],[11,16,14],[12,15,13],[13,16],[14,15]
           	];
      fi;
     
      if (arg[1] = "BG") and (arg[2] = 1) and (arg[3] = 3) then
        data[1] := 24;
        data[2] := 2;
        data[3] := [ [4,0],[3,0],[1,6,4],[2,5,3],[3,8,6],[4,7,5],[5,10,8],[6,9,7],
					[7,12,10],[8,11,9],[9,14,12],[10,13,11],[11,16,14],[12,15,13],[13,18,16],[14,17,15],
					[15,20,18],[16,19,17],[17,22,20],[18,21,19],[19,24,22],[20,23,21],[21,24],[22,23]
         	];
      fi;
    
    fi;
    
    if Length(arg) = 2 then
    
      if (arg[1] = "BG") and (arg[2] = 2)  then # Bongartz-Gabriel list no. 2
        data[1] := 66;
        data[2] := 2;
        data[3] := [ 
          [11,0],[19,0],[12,1],[13,3],#4
          [14,4],[15,5],[16,6],[17,7],#8
          [18,8],[20,2],[21,10],[22,1,11],#12
          [23,3,12],[24,4,13],[25,5,14],[26,6,15],#16
          [27,7,16],[28,8,17],[61,66],[30,2,19],#20
          [31,10,20],[32,11,21],[33,12,22],[34,13,23],#24
          [35,14,24],[36,15,25],[37,16,26],[38,17,27],#28
          [39,18,28],[19,54,61],[40,20,30],[42,21,31],#32
          [44,22,32],[46,23,33],[48,24,34],[50,25,35],#36
          [52,26,36],[54,27,37],[40,28,38],[30,55,38,54],#40
          [40,55],[39,31,41,40],[42,41],[56,32,43,42],#44
          [44,43],[45,57,33,44],[46,45],[47,58,34,46],#48
          [48,47],[49,59,35,48],[50,49],[51,60,36,50],#52
          [52,51],[53,61,37,52],[54,53],[29,42,39],#56
          [62,44,56],[63,46,57],[64,48,58],[65,50,59],#60
          [66,52,60],[56,29],[57,62],[58,63],#64
          [59,64],[60,65]
         	];
      elif (arg[1] = "BG") and (arg[2] = 5)  then # Bongartz-Gabriel list no. 5
        data[1] := 72;
        data[2] := 2;
        data[3] := [ 
          [11,0],[17,0],[12,1],[13,3],#4
		      [14,4],[18,2],[19,6],[20,7],#8
		      [21,8],[22,9],[23,10],[24,1,11],#12
		      [25,3,12],[26,4,13],[27,5,14],[63,69],#16
	    	  [28,16],[29,2,17],[6,30,18],[7,31,19],#20
	    	  [8,32,20],[9,33,21],[10,34,22],[11,35,23],#24
	     	  [12,36,24],[13,37,25],[14,38,26],[47,16,63],#28
	    	  [17,39,28],[18,40,29],[19,41,30],[20,42,31],#32
	    	  [21,43,32],[22,44,33],[23,45,34],[24,46,35],#36
	    	  [25,47,36],[26,39,37],[37,48,28,47],[29,49,38,39],#40
	    	  [30,50,57,40],[31,51,58,41],[32,52,59,42],[33,53,60,43],#44
	    	  [34,54,61,44],[35,55,62,45],[36,56,63,46],[47,56],#48
	    	  [39,48],[40,49],[41,50],[42,51],[43,52],[44,53],#54
	    	  [45,54],[46,55],[40,27,38],[41,64,57],#58
	    	  [42,65,58],[43,66,59],[44,67,60],[45,68,61],#62
	    	  [46,69,62],[57,15,27],[58,70,64],[59,71,65],#66
    		  [60,72,66],[61,67],[62,68],[64,15],#70
	    	  [65,70],[66,71] #72
		 
         	];
      elif (arg[1] = "BG") and (arg[2] = 7)  then # Bongartz-Gabriel list no. 7
        data[1] := 143;
        data[2] := 2;
        data[3] := [ 
          [16,0],[10,0],[11,2],[17,1],[18,4],#5
          [19,5],[20,6],[21,7],[22,8],[23,9],#10
          [2,24,10],[3,25,11],[26,12],[27,13],[28,14],#15
          [32,31],[1,33,16],[4,34,17],[5,35,18],[6,36,19],#20
          [7,37,20],[8,38,21],[9,39,22],[10,40,23],[11,41,24],#25
          [12,42,25],[13,43,26],[14,44,27],[15,45,28],[46,29],#30
          [47,30],[31,48,47],[16,49,32],[17,50,33],[18,51,34],#35
          [19,52,35],[20,53,36],[21,54,37],[22,55,38],[23,56,39],#40
          [24,57,40],[25,58,41],[26,59,42],[27,60,43],[28,61,44],#45
		      [29,62,45],[30,63,46],[47,64,63],[32,65,48],[33,66,49],#50
          [34,67,50],[35,68,51],[36,69,52],[37,70,53],[38,71,54],#55
          [39,72,55],[40,73,56],[41,74,57],[42,75,58],[43,76,59],#60
          [44,77,60],[45,78,61],[46,79,62],[63,96,79],[48,97,64],#65
          [49,98,65],[50,99,66],[51,100,67],[52,101,68],[53,102,69],#70
          [54,103,70],[55,104,71],[56,105,72],[57,106,73],[58,107,74],#75
          [59,108,75],[60,109,76],[61,110,77],[62,111,78],[111,95],#80
          [96,80],[97,81],[98,82],[99,83],[100,84],#85
          [101,85],[102,86],[103,87],[104,88],[105,89],#90
          [106,90],[107,91],[108,92],[109,93],[110,94],#95
          [79,80,112,111],[64,81,113,96],[65,82,114,97],[66,83,115,98],[67,84,116,99],#100
          [68,85,117,100],[69,86,118,101],[70,87,119,102],[71,88,120,103],[72,89,121,104],#105
          [73,90,122,105],[74,91,123,106],[75,92,124,107],[76,93,125,108],[77,94,126,109],#110
          [78,95,127,110],[111,128,127],[96,129,112],[97,130,113],[98,131,114],#115
          [99,132,115],[100,133,116],[101,134,117],[102,135,118],[103,136,119],#120
          [104,137,120],[105,138,121],[106,139,122],[107,140,123],[108,141,124],#125
          [109,142,125],[110,143,126],[127,143],[112,128],[113,129],#130
          [114,130],[115,131],[116,132],[117,133],[118,134],#135
          [119,135],[120,136],[121,137],[122,138],[123,139],#140
          [124,140],[125,141],[126,142] #143
         	];
	  elif (arg[1] = "BG") and (arg[2] = 9)  then # Bongartz-Gabriel list no. 9
        data[1] := 28;
        data[2] := 2;
        data[3] := [ 
          [5,0],[21,26,0],[6,1],[11,7],#4
	    	  [8,4],[1,9,5],[3,10,6],[4,12,19,11],#8
	    	  [5,13,16,8],[6,14,17,9],[7,15,18,10],[11,15],#12
	    	  [8,12],[9,13],[10,14],[8,23,19],#16
	    	  [9,20,16],[10,21,17],[11,22,18],[16,25,23],#20
	    	  [17,20],[2,18,21],[19,24,22],[22,27,2],#24
	    	  [23,28,24],[25,28],[2,26],[24,27]#28
         ];
      elif (arg[1] = "BG") and (arg[2] = 12)  then # Bongartz-Gabriel list no. 12
        data[1] := 153;
        data[2] := 2;
        data[3] := [ 
          [21,0],[6,0],[7,2],[8,3],[22,1],#5
	    	  [23,5],[2,24,6],[3,25,7],[4,26,8],[27,9],#10
         ];
	      for i in [1..7] do Add(data[3],[27+i,9+i]); od; #17
		    Append(data[3],[ 
		     [17,35,34],[36,18],[37,19],#20
		     [38,20],[1,39,21],[5,40,22] #23
		     ]);
		    for i in [1..11] do Add(data[3],[5+i,40+i,22+i]); od; #34
		    Append(data[3], [[34,52,51],[18,53,35]]); #36
		    for i in [1..15] do Add(data[3],[18+i,53+i,35+i]); od; #51
		    Append(data[3], [[51,69,68],[35,70,52]]); #53
		    for i in [1..15] do Add(data[3],[35+i,70+i,52+i]); od; #68
		    Append(data[3], [[68,103,85],[52,104,69]]); #70
		    for i in [1..15] do Add(data[3],[52+i,104+i,69+i]); od; #85
		    Append(data[3], [[119,102],[103,86]]); #87
		    for i in [1..15] do Add(data[3],[103+i,86+i]); od; #102
		    Append(data[3], [[85,86,120,119],[69,87,121,103]]); #104
		    for i in [1..15] do Add(data[3],[69+i,87+i,121+i,103+i]); od; #119
		    Append(data[3], [[119,137,136],[103,138,120]]); #121
		    for i in [1..15] do Add(data[3],[103+i,138+i,120+i]); od; #136
		    Append(data[3], [[136,153],[120,137]]); #138
		    for i in [1..15] do Add(data[3],[120+i,137+i]); od; #153
	    elif (arg[1] = "BG") and (arg[2] = 13)  then # Bongartz-Gabriel list no. 13
        data[1] := 26;
        data[2] := 2;
        data[3] := [ 
          [2,0],[8,0],[4,1],[1,9,2],#4
		      [3,10,4],[11,5],[12,6],[13,20],#8
		      [2,15,8],[4,17,9],[5,19,10],[6,13,11],#12
		      [11,20,23,19],[13,26,23],[8,12,14,13],[15,14],#16
		      [9,16,21,15],[17,16],[10,18,22,17],[19,18],#20
		      [7,15,12],[17,24,21],[19,25,22],[21,7],#24
		      [22,24],[23,25] #26
         ];
	    elif (arg[1] = "BG") and (arg[2] = 14)  then # Bongartz-Gabriel list no. 14
        data[1] := 20;
        data[2] := 2;
        data[3] := [ 
          [5,0],[20,0],[1,6,5],[7,3],#4
		      [8,4],[5,9,8],[3,10,6],[4,11,7],#8
		      [8,14,11],[6,12,9],[7,13,10],[9,15,20,14],#12
		      [10,16,18,12],[11,17,19,13],[14,17],[12,15],#16
		      [13,16],[2,12,20],[13,18],[14,19] # 20
         ];
      fi;
      
      
    fi;
    
    
    if Length(arg) = 1 then
      if (arg[1] = "D4 subspace")  then
        data[1] := 12;
        data[2] := 4;
        data[3] := [
           [0], [1,0], [1,0], [1,0],
           [2,3,4,1],[5,2],[5,3],[5,4],
           [6,7,8,5],[9,6],[9,7],[9,8]
         	];
      fi;
      if (arg[1] = "E6 subspace")  then
        data[1] := 36;
        data[2] := 6;
        data[3] := [ [0], [1,0], [2,0], [1,0], [1,0], [5,0],
           [2,4,5,1], [3,7,2], [8,3], [7,4],
           [7,6,5], [11,6], [8,10,11,7], [9,13,8], 
           [14,9], [13,10], [13,12,11], [17,12], 
           [14,16,17,13], [15,19,14], [20,15], [19,16],
           [19,18,17], [23,18],[20,22,23,19],[21,25,20],
           [26,21],[25,22],[25,24,23],[29,24],
           [26,28,29,25],[27,31,26],[32,27],[31,28],
           [31,30,29],[35,30]
         	];
      fi;
	  if (arg[1] = "A3 zero")  then
        data[1] := 5;
        data[2] := 3;
        data[3] := [
           [0], [1,0], [4,0], [2,1], [3,4]
         	];
      fi;
	  if (arg[1] = "A1,2 zero")  then
        data[1] := 9;
        data[2] := 3;
        data[3] := [
           [0], [1,0], [4,1,0], [7,5], 
		   [3,4],[3,2,1],[5,6,3],[6,2],[7,8,6]
         	];
      fi;
	  if (arg[1] = "A2,2 comm")  then
        data[1] := 11;
        data[2] := 4;
        data[3] := [
           [0], [1,0], [1,0], [5,0],
		   [2,3,1], [5,2], [5,3], [6,7,4,5],
		   [8,6], [8,7], [9,10,8]
         	];
      fi;
	  if (arg[1] = "R nilp")  then
        data[1] := 7;
        data[2] := 2;
        data[3] := [
           [4,0], [1,0], [1,6,4], [3,6], 
		   [3,2,1], [4,5,3], [5,2]
         	];
      fi;
    fi;
    
    
    if Length(data) <> 3 then
      Error("wrong number or type of arguments!");
    fi;
    
    return data;
  end
 ); # PredefARQuivers
 
 
 
 