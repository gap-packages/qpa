InstallMethod( GBNPGroebnerBasis,
  "call GBNP for Groebner Basis",
  true,
  [ IsFLMLOR ],
  0,
  function( I )
    local A, GB, gens;

#    if not IsFamilyElementOfPathRing(ElementsFamily(FamilyObj(I))) then
#      TryNextMethod();
#    fi;

    # Get parent algebra:
    A := LeftActingRingOfIdeal(I);

    # Get ideal generators:
    gens := GeneratorsOfIdeal(I);

    return GBNPGroebnerBasis( gens, A );

  end
);


InstallMethod( GBNPGroebnerBasis,
  "call GBNP for Groebner Basis",
  true,
  [ IsList,
    IsPathAlgebra ],
  0,
  function( els, pa )
    local q,ord,creps,grob,pgrob;

    creps := [];

    q := QuiverOfPathAlgebra(pa);   
    ord := OrderingOfAlgebra(pa);

    #  Check that all elements are in 
    #  the given path algebra 'pa', and that they
    #  are in the Arrow Ideal J,
    #  (i.e. are not vertices):
    if (QPA_InArrowIdeal(els,pa)) then

      # Should convert all given elements 
      #  to their uniform components:
      els := MakeUniform(els);

      # Convert list of path algebra elements
      #  to Cohen's format:
      creps := QPA_Path2Cohen(els);

      # Call Cohen to get Groebner basis:
      grob := SGrobner(creps);

      # Convert results back to path algebra
      #  elements:
      pgrob :=  QPA_Cohen2Path(grob,pa);

    else
      Print("Please make sure all elements are in the given path algebra,",
            "\nand each summand of each element is not (only) a constant",
	    "\ntimes a vertex.\n");
      pgrob := false;
    fi;

    return pgrob;

  end
);


InstallMethod( GBNPGroebnerBasisNC,
  "call GBNP for Groebner Basis",
  true,
  [ IsList,
    IsPathAlgebra ],
  0,
  function( els, pa )
    local q,ord,creps,grob,pgrob,numv,parels;

    creps := [];

    q := QuiverOfPathAlgebra(pa);   
    ord := OrderingOfAlgebra(pa);

    # Check that all elements are elements in
    #  the given path algebra 'pa':
    if ForAll( els, x -> \in(x,pa) ) then

      # Should convert all given elements 
      #  to their uniform components:
      els := MakeUniform(els);

      # Add relations to preserve structure of path
      #  algebra (first we get number of vertices,
      #  if there's only one, we have a free algebra):     
      numv := NumberOfVertices(q); 

      if numv > 1 then

        Print("The given path algebra is not a free algebra,\n",
	      " adding relations to preserve path algebra structure.\n");

        # Convert list of path algebra elements
        #  to Cohen's format:
        creps := QPA_Path2Cohen(els);

        parels := QPA_RelationsForPathAlgebra(pa);
#        Print(parels,"\n");        
#        Print(QPA_Cohen2Path(parels,pa),"\n");        

        # Combine the lists:
	Append(creps,parels);

      else

        Print("The given path algebra is isomorphic to a free algebra.\n");

        # Convert list of path algebra elements
        #  to Cohen's format for this special
	#  case where our path algebra is isomorphic
	#  to a free algebra, i.e. pa has only one
	#  vertex:
        creps := QPA_Path2CohenFree(els);

      fi;

      # Call Cohen to get Groebner basis:
      grob := SGrobner(creps);

#      Print(grob,"\n");        

      # Convert results back to path algebra
      #  elements:
      pgrob :=  QPA_Cohen2Path(grob,pa);

      # Remove any zeroes we may've gathered from adding relations:
      pgrob := Filtered(pgrob, x -> x <> Zero(pa));

    else
      Print("\nPlease make sure all elements are in the path algebra ",pa,
            ".\n");
      pgrob := false;
    fi;

    return pgrob;

  end
);


# Convert to format Cohen expects.
InstallMethod( QPA_Path2Cohen,
  "convert path algebra elements to Cohen reps",
  true,
  [ IsList ],
  0,
  function( els )
    local e,i,oldrep,coefs,mons,creps;

    creps := [];

    for e in els do
      coefs := [];
      mons := [];

      oldrep := ExtRepOfObj(e)[2];

      for i in [ 2,4 .. Length(oldrep) ] do
        Add(coefs,oldrep[i]);
        Add(mons,oldrep[i-1]);
      od;

      Add(creps, [mons,coefs]);
    od;

    return creps;
  end
);


# Convert to format Cohen expects.
#
# NOTE: we are assuming that the elements in given
#  list are from a path algebra isomorphic to a free
#  algebra, i.e. there is only one vertex.
#
InstallMethod( QPA_Path2CohenFree,
  "convert path (special case: free) algebra elements to Cohen reps",
  true,
  [ IsList ],
  0,
  function( els )
    local e,i,oldrep,coefs,mons,creps;

    creps := [];

    for e in els do
      coefs := [];
      mons := [];

      oldrep := ExtRepOfObj(e)[2];

      for i in [ 2,4 .. Length(oldrep) ] do

        Add(coefs,oldrep[i]);

	# Check to see if this is our identity
	# (same as our only vertex, our first
	#  generator), if so correctly convert:
	if ( oldrep[i-1] = [1] ) then
          Add(mons,[]);
	else
          Add(mons,oldrep[i-1]);
	fi;

      od;

      Add(creps, [mons,coefs]);
    od;

    return creps;
  end
);


# Convert to GAP path algebra element from a Cohen representation.
InstallMethod( QPA_Cohen2Path,
  "convert Cohen reps to path algebra elements",
  true,
  [ IsList, IsPathAlgebra ],
  0,
  function( reps, pa )
    local e,i,j,mons,coefs,gens,els,word,poly,zero,one;

    gens := GeneratorsOfAlgebra(pa);
    zero := Zero(pa);
    one := One(pa);
    els := [];

    for e in reps do

#      Print("Rep: ",e,"\n");
      poly := zero;

      mons := e[1]; 
      coefs := e[2];

      # For each term in rep we build a word:
      for i in [1 .. Length(mons)] do
#        Print("\tMonomial: ",mons[i],"\n");

        # Build word:
	# NOTE:
	#   this might be dangerous in the general case:
        word := coefs[i]*one;
        for j in [ 1 .. Length(mons[i]) ] do
          word := word*gens[mons[i][j]];
        od;

#        Print("\t\tWord is: ",word,"\n");

        poly := poly + word;

      od; 

      # Add new element to list for return:
      Add(els,poly);

    od;

    return els;

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
