InstallMethod( GBNPGroebnerBasis,
  "compute a Groebner Basis (no longer using GBNP)",
  [ IsList, IsPathAlgebra ],
  function(els, A)
    local gb, el, el_tip,
          n, i, j, x, y, k, l, r, b, c,
          overlap, remainder;

    if not QPA_InArrowIdeal(els, A) then
      Error("elements do not belong to the arrow ideal of the path algebra");
    fi;

    els := ReducedList(MakeUniform(els), A);

    gb := [];

    while Length(els) > 0 do
      for el in els do
        el_tip := Tip(el);
        Add(gb, el/TipCoefficient(el_tip));
      od;

      n := Length(gb);
      els := [];

      for i in [1..n] do
        x := TipWalk(gb[i]);
        k := Length(x);

        for j in [1..n] do
          y := TipWalk(gb[j]);
          l := Length(y);

          for r in [Maximum(0, k-l)..k-1] do
            if x{[r+1..k]} = y{[1..k-r]} then
              b := x{[1..r]};
              c := y{[k-r+1..l]};

              overlap := gb[i]*Product(c, One(A)) - Product(b, One(A))*gb[j];
              remainder := RemainderOfDivision(overlap, gb, A);

              if not IsZero(remainder) then
                AddSet(els, remainder);
              fi;
            fi;
          od;
        od;
      od;
    od;

    gb := TipReducedList(gb, A);
    gb := ReducedList(gb, A);

    return gb;
  end
);


InstallMethod( ReducedList,
  "for a list of path-algebra elements",
  [ IsList, IsPathAlgebra ],
  function(els, A)
    local res, i, r;

    res := Filtered(els, el -> not IsZero(el));

    i := Length(res);
    while i > 0 do
      r := RemainderOfDivision(res[i], res{Concatenation([1..i-1], [i+1..Length(res)])}, A);

      if IsZero(r) then
        Remove(res, i);
      else
        res[i] := r;
      fi;

      i := i-1;
    od;

    return res;
  end
);


InstallMethod( TipReducedList,
  "for a list of path-algebra elements",
  [ IsList, IsPathAlgebra ],
  function(els, A)
    local res, el, i;

    res := [];

    for el in els do
      if not IsZero(el) then
        AddSet(res, el);
      fi;
    od;

    i := Length(res);
    while i > 0 do
      if ForAny([1..i-1], j -> LeftmostOccurrence(TipWalk(res[i]), TipWalk(res[j])) <> fail) then
        Remove(res, i);
      fi;
      i := i-1;
    od;

    return res;
  end
);


InstallMethod( RemainderOfDivision,
  "for a path-algebra element and a list of path-algebra elements",
  [ IsElementOfMagmaRingModuloRelations, IsList, IsPathAlgebra ],
  function(y, X, A)
    local r, n, y_tip, y_wtip, divided, i, p, u, v;

    r := Zero(A);
    n := Length(X);

    while not IsZero(y) do
      y_tip := Tip(y);
      y_wtip := TipWalk(y_tip);

      divided := false;

      for i in [1..n] do
        p := LeftmostOccurrence(y_wtip, TipWalk(X[i]));

        if p <> fail then
          u := Product(y_wtip{[1..p[1]-1]}, One(A));
          v := Product(y_wtip{[p[2]+1..Length(y_wtip)]}, One(A));

          y := y - TipCoefficient(y_tip)/TipCoefficient(X[i]) * u*X[i]*v;

          divided := true;
          break;
        fi;
      od;

      if not divided then
        r := r + y_tip;
        y := y - y_tip;
      fi;
    od;

    return r;
  end
);


InstallMethod( LeftmostOccurrence,
  "find second list as sublist of first list",
  [ IsList, IsList ],
  function(b, c)
    local lb, lc, i;

    lb := Length(b);
    lc := Length(c);

    for i in [1..lb-lc+1] do
      if b{[i..i+lc-1]} = c then
        return [i, i+lc-1];
      fi;
    od;

    return fail;
  end
);


InstallMethod( MakeUniform,
  "return array of uniform elements",
  true,
  [ IsElementOfMagmaRingModuloRelations ],
  0,
  function( el )
    local l,s,t,v1,v2,terms,uterms,newel;

    uterms := [];

    # get coefficients and monomials for all
    #  terms for el:
    terms := CoefficientsAndMagmaElements(el);

    # Get source vertices for all terms:
    l := List(terms{[1,3..Length(terms)-1]}, x -> SourceOfPath(x));

    # Make the list unique:
    s := Unique(l);

    # Get terminus vertices for all terms:
    l := List(terms{[1,3..Length(terms)-1]}, x -> TargetOfPath(x));

    # Make the list unique:
    t := Unique(l);

    # Create uniformized array based on el:
    for v1 in s do
      for v2 in t do
        newel := v1*el*v2;
	if not IsZero(newel) then
	  Add(uterms,newel);
	fi;
      od;
    od;

    return uterms;

  end
);


InstallMethod( MakeUniform,
  "return array of uniform elements",
  true,
  [ IsList ],
  0,
  function( l )
    return Unique(Flat(List(l,x -> MakeUniform(x))));
  end
);


InstallMethod( QPA_InArrowIdeal,
  "return array of uniform elements",
  true,
  [ IsList, IsPathAlgebra ],
  0,
  function( els, pa )
    local q,numv,retval,e,extrep,term;

    retval := false;

    # first check all elements are in given path algebra:
    if ForAll(els, x -> \in(x,pa) ) then

      # Get quiver:
      q := QuiverOfPathAlgebra(pa);

      # Get number vertices:
      numv := NumberOfVertices(q);

      # The following checks that each element in
      #  given list is an element of the given
      #  path algebra and that each element 
      #  consists of only non-vertex-only paths:

      retval := true;

      for e in els do

        # Get exterior (list) representation of element:
        extrep := ExtRepOfObj(e)[2];

	# Now iterate over only the monomials of this representation
	#  (the odd indices):
        for term in extrep{[1,3..Length(extrep)-1]} do

	  # Check to see if single generator word, and if that
	  # generator is a vertex:
	  if Length(term) = 1 then
	    if term[1] <= numv then
	      # bail out this is a vertex:
	      return false;
	    fi;
          fi;

	od;
      od;

    fi; 

    return retval;

  end
);


# The function returns a list of relations necessary to
# enforce structure from a given path algebra in the
# ideal for a quotient of a free algebra.
# This returned list has format that GBNP expects.
InstallMethod( QPA_RelationsForPathAlgebra,
  "create relations of path algebra for free algebra",
  true,
  [ IsPathAlgebra ],
  0,
  function( pa )
    local i,j,quiv,newrel,rels,gens,numa,numg,numv,terms,coefs,
          ar,arep,s,t,zero,algarrows;

    rels := [];
    terms := [];
    coefs := [];

    # Get quiver for given path algbra:
    quiv := QuiverOfPathAlgebra(pa);   

    # Get number of total generators for algebra:
    gens := GeneratorsOfAlgebra(pa);
    numg := Length(gens);

    # Get number of vertices:
    numv := NumberOfVertices(quiv); 

    # Calculate number of arrows:
    numa := numg - numv;

    ######
    #
    # Add relation that sum of primitive orthogonal
    #  idempotents (vertices) equals one:
    #  Example relation: u+v+w-1  
    #
    for i in [1..numv] do
      Add(terms,[i]);
      Add(coefs,1);
    od;

    # Add (empty) term for identity element:
    Add(terms,[]);

    # Add -1 as coefficient of identity element:
    Add(coefs,-1);

    # Construct this relation in GBNP format:
    newrel := [terms,coefs];

    # Add this new relation to our set of all relations:
    Add(rels,newrel);
    #
    ######


    ######
    #
    # Add relation for each primitive idempotent:
    #  Example relations: u*u-u  
    #
    for i in [1..numv] do
       newrel := [[[i,i],[i]],[1,-1]];
       Add(rels,newrel);
    od;
    #
    ######


    ######
    #
    # Add orthogonality relations for each pair of vertices:
    #  Example relations: u*v,v*u  (where v<>u) 
    #
    for i in [1..numv] do
      for j in [i+1..numv] do
        Add(rels,[[[i,j]],[1]]);
        Add(rels,[[[j,i]],[1]]);
      od;
    od;
    #
    ######


    ######
    #
    # For each arrow "a" in the quiver, add a relation
    #  that collapses u*a to a and a*v to a
    #  where u is the starting vertex for a, and
    #  v is the ending vertex for a.
    #  Example relations:  u*a-a, a*v-a
    #

    # Iterate over the set of all arrows:
    for ar in ArrowsOfQuiver(quiv) do

      # Get representation for arrow:
      arep := ExtRepOfObj(ar)[1];

      # Get starting and ending vertex for arrow:
      s := ExtRepOfObj(SourceOfPath(ar))[1];
      t := ExtRepOfObj(TargetOfPath(ar))[1];

      # Create and add relation for arrow:
      Add(rels,[[[s,arep],[arep]],[1,-1]]);
      Add(rels,[[[arep,t],[arep]],[1,-1]]);

    od;
    #
    ######


    ######
    #
    # Create arrow orthogonality relations, i.e. if
    # a*b = 0 then we need to add the representation
    # for the relation a*b to our set 'rels'.

    # Get zero for algebra:
    zero := Zero(pa);

    # Get arrows within the path algebra so that
    #  we can multiply them:
    algarrows := gens{[numv+1..numg]};

    for i in [numv+1..numg] do
      for j in [numv+1..numg] do
        # Check if relation is zero, if so, add
	#  appropriate relation:
	if (gens[i]*gens[j] = zero) then
          Add(rels,[[[i,j]],[1]]);
	fi;
      od;
    od;
    #
    ######

    return rels;

  end
);
