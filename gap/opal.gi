# GAP Implementation
# This file was generated from 
# $Id: opal.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $

Opal := rec();

InstallMethod( AlgebraToOpal,
    "for quotients of path algebras and output streams",
    true,
    [IsFLMLOR, IsOutputStream, IsList], 0,
    function(A, opalFile, verts)
        local Q, ordering, lexTable, vertices, arrows, v, a, first,
              elementFam, zero, one, nextOrder, weights;
        
        Q := QuiverOfPathAlgebra(A);
        ordering := OrderingOfQuiver(Q);
        if not IsAdmissibleOrdering(ordering) then
            Error("The ordering of elements must be admissible.");
        fi;
        elementFam := ElementsFamily(FamilyObj(A));
        zero := Zero(LeftActingDomain(A));
        one := One(LeftActingDomain(A));

        PrintTo(opalFile, "let g = Graph( [ ");
        if HasLexicographicTable(ordering) then
            lexTable := LexicographicTable(ordering);
        else
            lexTable := GeneratorsOfQuiver(Q);
        fi;
        # First just the vertices    
        vertices := Reversed(Filtered(lexTable, IsVertex));

        first := true;
        for v in vertices do
            if not first then
                AppendTo(opalFile, ", ");
            else
                first := false;
            fi;
            
            Opal.(String(v!.gen_pos)) := ObjByExtRep(elementFam, 
                                                     [zero, [ExtRepOfObj(v), one]]);
            AppendTo(opalFile, "(Opal'", v!.gen_pos, ")");
            Add(verts, v!.gen_pos);
        od;
        AppendTo(opalFile, "],\n");
        # Now the arrows
        AppendTo(opalFile, "[ ");
        arrows := Reversed(Filtered(lexTable, IsArrow));

        first := true;
        for a in arrows do
            if not first then
                AppendTo(opalFile, ",\n");
            else
                first := false;
            fi;
            Opal.(String(a!.gen_pos)) := ObjByExtRep(elementFam, 
                                                     [zero, [ExtRepOfObj(a), one]]);
            AppendTo(opalFile, "Opal'", a!.gen_pos, " : ");
            AppendTo(opalFile, "(Opal'", SourceOfPath(a)!.gen_pos, ") -> ");
            AppendTo(opalFile, "(Opal'", TargetOfPath(a)!.gen_pos, ")");
        od;
        AppendTo(opalFile, " ] );\n");
        nextOrder := ordering;
        while not IsWeightOrdering(nextOrder) and HasNextOrdering(nextOrder) do
            nextOrder := NextOrdering(nextOrder);
        od;
        if IsWeightOrdering(nextOrder) then
            weights := nextOrder!.weight;
            AppendTo(opalFile, "let g = ArcWeight(g, [\n");

            first := true;
            for a in arrows do
                if not first then
                    AppendTo(opalFile, ",\n");
                else
                    first := false;
                fi;
                AppendTo(opalFile, "Opal'", a!.gen_pos, ":", 
                         weights[a!.gen_pos - Length(vertices)]);
            od;
            AppendTo(opalFile, "]);\n");

            while HasNextOrdering(nextOrder) do
                nextOrder := NextOrdering(nextOrder);
                if IsWeightOrdering(nextOrder) then
                    Error("Two weight orderings encountered. Opal only supports one.");
                fi;
            od;
        fi;
end );


InstallMethod( AssumePathAlgebra,
    "for arbitrary fields",
    true,
    [IsField, IsOutputStream], 0,
    function( F, outputStream )
        Error("AssumePathAlgebra is not supported for the field: ", F);
end );


InstallMethod( AssumePathAlgebra,
    "for finite prime fields",
    true,
    [IsPrimeField and IsFinite, IsOutputStream], 0,
    function( F, outputStream)
        # Make sure this is a prime field for now
        AppendTo(outputStream,
                 "assume pathalgebra(integer(", Characteristic(F),"),g);\n");
end );


InstallMethod( AssumePathAlgebra,
    "for the rationals",
    true,
    [IsRationals, IsOutputStream], 0,
    function( F, outputStream )
        AppendTo(outputStream,
                 "assume pathalgebra( Rationals, g );\n");
end );


InstallMethod( WriteCoefficientToOpal,
    "for finite prime fields",
    true,
    [IsFFE, IsOutputStream], 0,
    function( coeff, opalFile )
        if DegreeFFE(coeff) <> 1 then
            TryNextMethod();
        fi;
        AppendTo(opalFile, IntFFE(coeff));
end);


InstallMethod( WriteCoefficientToOpal,
    "for the rationals",
    true,
    [IsRat, IsOutputStream], 0,
    function( coeff, opalFile )
        AppendTo(opalFile, coeff);
end);


InstallMethod(WriteRelatorsToOpal, 
    "for quotients of path algebras and output streams",
    true,
    [IsFLMLOR, IsOutputStream, IsList], 0,
    function(I, opalFile, verts)
        local relators, rel, extRel, coeff, term, generator, 
              i, firstTerm, firstGen, firstRel, F;
        AppendTo(opalFile, "let F = {\n");

        relators := GeneratorsOfIdeal(I);
        firstRel := true;
        F := LeftActingDomain(I);

        for rel in relators do
            if not firstRel then
                AppendTo(opalFile, ",\n");
            else
                firstRel := false;
            fi;
            extRel := ExtRepOfObj(rel)[2];

            # Make sure the first term is positive
            # because opal doesn't handle unary -
            if extRel[2] < Zero(F) then
                rel := -One(F)*rel;
                extRel := ExtRepOfObj(rel)[2];
            fi;
            
            firstTerm := true;
            for i in [1,3..Length(extRel) - 1] do
                coeff := extRel[i+1];
                term := extRel[i];
                if not firstTerm then
                    if coeff < Zero(F) then
                        coeff := -coeff;
                        AppendTo(opalFile, " - ");
                    else
                        AppendTo(opalFile, " + ");
                    fi;
                else
                    firstTerm := false;
                fi;
                firstGen := true;
                WriteCoefficientToOpal(coeff, opalFile);
                AppendTo(opalFile,"*");
                for generator in term do
                    if not firstGen then
                        AppendTo(opalFile, "*");
                    else
                        firstGen := false;
                    fi;
                    if generator in verts then
                        AppendTo(opalFile, "(");
                    fi; 
                    AppendTo(opalFile, "Opal'", generator);
                    if generator in verts then
                        AppendTo(opalFile, ")");
                    fi; 
                od;
            od;
        od;
        AppendTo(opalFile, "\n};\n");
end);


InstallGlobalFunction( WriteFinalOpalCommands, function( I, opalFile )
    AppendTo(opalFile, "let B = Basis(F);\n");
    AppendTo(opalFile, "B;\n");
    AppendTo(opalFile, "IsGrobner(B);\n");
end );


InstallGlobalFunction( WriteFinalOpalCommandsWithBound, function( I, n, opalFile )
    AppendTo(opalFile, "let B = Basis(F,",n,");\n");
    AppendTo(opalFile, "B;\n");
    AppendTo(opalFile, "IsGrobner(B);\n");
end );


InstallGlobalFunction( WriteFinalOpalCommandsFinite, function( I, opalFile )
    AppendTo(opalFile, "let B = FinBasis(F);\n");
    AppendTo(opalFile, "B;\n");
    AppendTo(opalFile, "IsGrobner(B);\n");
end );


InstallMethod( OrderingToOpal,
    "for left lexicographic orderings",
    true,
    [IsLeftLexicographicOrdering], 0,
    function(ordering)
        local retval;

        retval := "lex";
        if HasNextOrdering(ordering) then
            retval := Concatenation(retval, " ", 
                      OrderingToOpal(NextOrdering(ordering)));
        fi;
        return retval;     
end );


InstallMethod( OrderingToOpal,
    "for right lexicographic orderings",
    true,
    [IsRightLexicographicOrdering], 0,
    function(ordering)
        local retval;

        retval := "rlex";
        if HasNextOrdering(ordering) then
            retval := Concatenation(retval, " ", 
                      OrderingToOpal(NextOrdering(ordering)));
        fi;
        return retval;     
end );


InstallMethod( OrderingToOpal,
    "for length orderings",
    true,
    [IsLengthOrdering], 0,
    function(ordering)
        local retval;

        retval := "length";
        if HasNextOrdering(ordering) then
            retval := Concatenation(retval, " ", 
                      OrderingToOpal(NextOrdering(ordering)));
        fi;
        return retval;     
end );


InstallMethod( OrderingToOpal,
    "for reverse orderings",
    true,
    [IsReverseOrdering], 0,
    function(ordering)
        local retval;

        retval := "rev";
        if HasNextOrdering(ordering) then
            retval := Concatenation(retval, " ", 
                      OrderingToOpal(NextOrdering(ordering)));
        fi;
        return retval;     
end );


InstallMethod( OrderingToOpal,
    "for left vector ordering",
    true,
    [IsLeftVectorOrdering], 0,
    function(ordering)
        local retval;

        retval := "vec";
        if HasNextOrdering(ordering) then
            retval := Concatenation(retval, " ", 
                      OrderingToOpal(NextOrdering(ordering)));
        fi;
        return retval;     
end );


InstallMethod( OrderingToOpal,
    "for right vector ordering",
    true,
    [IsRightVectorOrdering], 0,
    function(ordering)
        local retval;

        retval := "rvec";
        if HasNextOrdering(ordering) then
            retval := Concatenation(retval, " ", 
                      OrderingToOpal(NextOrdering(ordering)));
        fi;
        return retval;     
end );


InstallMethod( OrderingToOpal,
    "for weight ordering",
    true,
    [IsWeightOrdering], 0,
    function(ordering)
        local retval;

        retval := "weight";
        if HasNextOrdering(ordering) then
            retval := Concatenation(retval, " ", 
                      OrderingToOpal(NextOrdering(ordering)));
        fi;
        return retval;     
end );


InstallGlobalFunction( ExecuteOpal,
    function(ordering, opalFilename, gapFilename)
    local packageDir, gapopal, opal, orderingName;

    packageDir := DirectoriesPackagePrograms("hopf");
    gapopal := Filename(packageDir, "gapopal");
    opal := Filename(packageDir, "opal");
    orderingName := Concatenation("\"", OrderingToOpal(ordering), "\"");
    if gapopal = fail then
        Error("The gapopal script does not exist. Cannot execute opal.");
    fi;
    Info(InfoOpal, 1, "The gapopal script path is: ", gapopal);
    Exec(gapopal, opal, orderingName, opalFilename, gapFilename);
end );


InstallMethod( OpalGroebnerBasis,
    "for a path algebra",
    true,
    [ IsFLMLOR ], 0,
    function( I )
        local tempdir, gapFilename, opalFilename, opalFile, GB, A, verts;

        if not IsFamilyElementOfPathRing(ElementsFamily(FamilyObj(I))) then
            TryNextMethod();
        fi;

        # Setup the temp directory and input and output files.
        tempdir := DirectoryTemporary();
        gapFilename := Filename(tempdir, "opalOutput");
        opalFilename := Filename(tempdir, "opalInput");
        Info(InfoOpal, 1, "Opal input file name is: ", opalFilename);
        Info(InfoOpal, 1, "GAP input file name is: ", gapFilename);
        opalFile := OutputTextFile(opalFilename, false);

        # Print out the information to an output file
        A := LeftActingRingOfIdeal(I);
        verts := [];

        AlgebraToOpal(A, opalFile, verts);
        AssumePathAlgebra(LeftActingDomain(I), opalFile);
        WriteRelatorsToOpal(I, opalFile, verts);
        WriteFinalOpalCommands(I, opalFile);

        # Close the output stream and execute Opal      
        CloseStream(opalFile);
        ExecuteOpal(OrderingOfAlgebra(A), opalFilename, gapFilename);
        Read(gapFilename); # Sets Opal.groebnerBasis

        GB := GroebnerBasis(I, Opal.groebnerBasis);
        SetIsCompleteGroebnerBasis(GB, Opal.isCompleteGroebnerBasis);

        return GB;
end );        


InstallMethod( OpalBoundedGroebnerBasis,
    "for a path algebra and a bound",
    true,
    [ IsFLMLOR , IsPosInt], 0,
    function( I , n )
        local tempdir, gapFilename, opalFilename, opalFile, GB, A, verts;

        if not IsFamilyElementOfPathRing(ElementsFamily(FamilyObj(I))) then
            TryNextMethod();
        fi;

        # Setup the temp directory and input and output files.
        tempdir := DirectoryTemporary();
        gapFilename := Filename(tempdir, "opalOutput");
        opalFilename := Filename(tempdir, "opalInput");
        Info(InfoOpal, 1, "Opal input file name is: ", opalFilename);
        Info(InfoOpal, 1, "GAP input file name is: ", gapFilename);
        opalFile := OutputTextFile(opalFilename, false);

        # Print out the information to an output file
        A := LeftActingRingOfIdeal(I);
        verts := [];
        AlgebraToOpal(A, opalFile, verts);
        AssumePathAlgebra(LeftActingDomain(I), opalFile);
        WriteRelatorsToOpal(I, opalFile, verts);
        WriteFinalOpalCommandsWithBound(I, n, opalFile);

        # Close the output stream and execute Opal      
        CloseStream(opalFile);
        ExecuteOpal(OrderingOfAlgebra(A), opalFilename, gapFilename);
        Read(gapFilename); # Sets Opal.groebnerBasis

        GB := GroebnerBasis(I, Opal.groebnerBasis);
        SetIsCompleteGroebnerBasis(GB, Opal.isCompleteGroebnerBasis);

        return GB;
end );


InstallMethod( OpalFiniteGroebnerBasis,
    "for a path algebra and a finite basis",
    true,
    [ IsFLMLOR ], 0,
    function( I )
        local tempdir, gapFilename, opalFilename, opalFile, GB, A, verts;

        if not IsFamilyElementOfPathRing(ElementsFamily(FamilyObj(I))) then
            TryNextMethod();
        fi;

        # Setup the temp directory and input and output files.
        tempdir := DirectoryTemporary();
        gapFilename := Filename(tempdir, "opalOutput");
        opalFilename := Filename(tempdir, "opalInput");
        Info(InfoOpal, 1, "Opal input file name is: ", opalFilename);
        Info(InfoOpal, 1, "GAP input file name is: ", gapFilename);
        opalFile := OutputTextFile(opalFilename, false);

        # Print out the information to an output file
        A := LeftActingRingOfIdeal(I);
        AlgebraToOpal(A, opalFile, verts);
        AssumePathAlgebra(LeftActingDomain(I), opalFile);
        WriteRelatorsToOpal(I, opalFile, verts);
        WriteFinalOpalCommandsFinite(I, opalFile);

        # Close the output stream and execute Opal      
        CloseStream(opalFile);
        ExecuteOpal(OrderingOfAlgebra(A), opalFilename, gapFilename);
        Read(gapFilename); # Sets Opal.groebnerBasis

        GB := GroebnerBasis(I, Opal.groebnerBasis);
        SetIsCompleteGroebnerBasis(GB, Opal.isCompleteGroebnerBasis);

        return GB;
end );  
