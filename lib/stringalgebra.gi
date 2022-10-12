InstallGlobalFunction (computeLPSarray,
  function(pat)
    local lps, i, len, m;
    m := Length(pat);
    lps := [];
    for i in [1..m] do
        Append(lps, [1]);
    od;
    len := 1;
    i := 2;
    while i <= m do
        if pat[i] = pat[len] then
            len := len + 1;
            lps[i] := len;
            i := i+1;
        else
            if len <> 1 then
                len := lps[len - 1];
            else
                lps[i] := 1;
                i := i + 1;
            fi;
        fi;
    od;
    for i in [1..Length(lps)] do
        lps[i] := lps[i] - 1;
    od;
    return lps;
end
);

InstallGlobalFunction (SourceOfArrow,
  function(Q, ch)
    local i, ch1;
    ch1 := ch;
    if SIntChar(ch) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
        ch := CharInt(SIntChar(ch) + 32);
    fi;
    if SIntChar(ch) - SIntChar('a') + 1 > NumberOfArrows(Q) then
        Error("The input string is invalid");
        return 0;
    fi;
    for i in [1..NumberOfArrows(Q)] do
        if String(ArrowsOfQuiver(Q)[i])[1] = ch then
            break;
        fi;
    od;
    if SIntChar(ch1) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
        return TargetOfPath(ArrowsOfQuiver(Q)[i]);
    else
        return SourceOfPath(ArrowsOfQuiver(Q)[i]);
    fi;
end
);

InstallGlobalFunction (TargetOfArrow,
  function(Q, ch)
    local i, ch1;
    ch1 := ch;
    if SIntChar(ch) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
        ch := CharInt(SIntChar(ch) + 32);
    fi;
    if SIntChar(ch) - SIntChar('a') + 1 > NumberOfArrows(Q) then
        Error("The input string is invalid");
        return 0;
    fi;
    for i in [1..NumberOfArrows(Q)] do
        if String(ArrowsOfQuiver(Q)[i])[1] = ch then
            break;
        fi;
    od;
    if SIntChar(ch1) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
        return SourceOfPath(ArrowsOfQuiver(Q)[i]);
    else
        return TargetOfPath(ArrowsOfQuiver(Q)[i]);
    fi;
end
);

InstallGlobalFunction (CyclicPermutationOfABand,
  function(band)
    local array, temp, i, j;
    array := [];
    for i in [1..Length(band)] do
        Append(array, [0]);
    od;
    for i in [1..Length(band)] do
        temp := band[1];
        for j in [1..Length(band) - 1] do
            band[j] := band[j + 1];
        od;
        band[j + 1] := temp;
        array[i] := Concatenation(band,"");
    od;
    return array;
end
);

InstallGlobalFunction (MaximumLengthOfDirectString,
  function(A)
    local Q, arr, i, max, x, y, array, sigma, eps;
    Q := QuiverOfPathAlgebra(A);
    array := QPAStringSigmaEps(A);
    sigma := array[1];
    eps := array[2];
    arr := [];
    for i in [1..NumberOfArrows(Q)] do
        x := [CharInt(i + 96)];
        while true do
            y := QPAStringDirectLeft(A, x, sigma, eps);
            if y = "Cannot Perform The Operation" then
                Append(arr, [Length(x)]);
                break;
            else x := y;
            fi;
        od;
    od;
    max := arr[1];
    for i in [2..Length(arr)] do
        if arr[i] > max then max := arr[i];
        fi;
    od;
    return max;
end
);

InstallGlobalFunction (NumberOfJoints,
  function(A)
    local Q, arr, i, temp, array, sigma, eps;
    Q := QuiverOfPathAlgebra(A);
    arr := [];
    array := QPAStringSigmaEps(A);
    sigma := array[1];
    eps := array[2];
    for i in [1..NumberOfArrows(Q)] do
      temp := QPAStringInverseRight(A, [CharInt(i+96)], sigma, eps);
      if temp <> "Cannot Perform The Operation" then
        Append(arr, [temp]);
      fi;
    od;
    return Length(arr);
end
);

InstallMethod( IsABand,
"for stringalgebras and strings",
true,
[ IsAlgebra, IsString ], 0,
function ( A, input_str )
    local Q, rho1, rho, vertices, arrows, onerel, tworel, str, output, w1, w2,
    n, i, j, k, arr, m_arr, small, capital, count, temp;

    if IsStringAlgebra(A) = false then
      Error("The input algebra is not a string algebra");
      return 0;
    fi;

    Q := QuiverOfPathAlgebra(A);
    rho1 := RelationsOfAlgebra(A);

    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);

    for i in [1..Length(vertices)] do
      if Length(String(vertices[i])) <> 2 or String(vertices[i])[1] <> 'v' or
        String(vertices[i])[2] <> String(i)[1] then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    for i in [1..Length(arrows)] do
      if Length(String(arrows[i])) <> 1 or String(arrows[i])[1] <> CharInt(i+96) then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    rho := [];
    for i in [1..Length(rho1)] do
        onerel := CoefficientsAndMagmaElements(rho1[i]);
        tworel := WalkOfPath(onerel[1]);
        str := Concatenation(List(tworel, a -> String(a)));
        Append(rho, [str]);
    od;

    output := IsValidString(A, input_str);
    if output = false then
        Error("The input string is invalid");
        return 0;
    fi;

    if output = true and input_str[1] = '(' then return false;
    fi;
    if Length(input_str) = 1 then return false;
    fi;

    n := Length(input_str);
    w1 := SourceOfArrow(Q, input_str[n]);
    w2 := TargetOfArrow(Q, input_str[1]);

    if w1 <> w2 then return false;
    fi;

    capital := 0;
    small := 0;

    for i in [1..Length(input_str)] do
        if SIntChar(input_str[i]) >= SIntChar('A') and
          SIntChar(input_str[i]) <= SIntChar('Z') then
            capital := 1;
        elif SIntChar(input_str[i]) >= SIntChar('a') and
          SIntChar(input_str[i]) <= SIntChar('z') then
            small := 1;
        fi;
    od;
    if capital = 0 or small = 0 then return false;
    fi;
    temp := Concatenation(input_str, input_str);
    if IsValidString(A, temp) = false then
        return false;
    fi;
    arr := computeLPSarray(input_str);
    if arr[Length(arr)] <> 0 and Length(input_str) mod (Length(input_str) - arr[Length(arr)]) = 0 then
        return false;
    else
        return true;
    fi;
end
);

InstallMethod( StringsLessThan,
"for stringalgebras and levels",
true,
[ IsAlgebra, IsPosInt ], 0,
function ( A, level )
    local Q, rho1, rho, vertices, arrows, onerel, tworel, str,making_tree, k,
    tree, array, i, vertex, kQ, gens, rel, j, temp, arr1, arr;

    if IsStringAlgebra(A) = false then
      Error("The input algebra is not a string algebra");
      return 0;
    fi;

    Q := QuiverOfPathAlgebra(A);
    rho1 := RelationsOfAlgebra(A);

    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);

    for i in [1..Length(vertices)] do
      if Length(String(vertices[i])) <> 2 or String(vertices[i])[1] <> 'v' or
        String(vertices[i])[2] <> String(i)[1] then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    for i in [1..Length(arrows)] do
      if Length(String(arrows[i])) <> 1 or String(arrows[i])[1] <> CharInt(i+96) then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    rho := [];
    for i in [1..Length(rho1)] do
        onerel := CoefficientsAndMagmaElements(rho1[i]);
        tworel := WalkOfPath(onerel[1]);
        str := Concatenation(List(tworel, a -> String(a)));
        Append(rho, [str]);
    od;

    making_tree := function(A, vertex, level)
        local tree, x, i, count, j, array, sigma,eps;
        array := QPAStringSigmaEps(A);
        sigma := array[1];
        eps := array[2];
        tree := [];
        tree[1] := [1, vertex];
        i := 1;
        count := 0;
        x := QPAStringDirectLeft(A, vertex, sigma, eps);
        if x <> "Cannot Perform The Operation" then
            Add(tree, [2 * 1, x]);
            j := 2;
        else
            count := count + 1;
        fi;

        x := QPAStringInverseLeft(A, vertex, sigma, eps);
        if x <> "Cannot Perform The Operation" then
            Add(tree, [2 * 1 + 1, x]);
            j := 3;
        else
            count := count + 1;
        fi;

        if count = 2 then return [];
        fi;
        i := 2;
        while true do
            if 2*tree[i][1] > 2^(level+1) - 1 then break;
            fi;
            x := QPAStringDirectLeft(A, tree[i][2], sigma, eps);
            if x <> "Cannot Perform The Operation" then
                Add(tree, [2 * tree[i][1], x]);
                j := 2*tree[i][1];
            fi;

            x := QPAStringInverseLeft(A, tree[i][2], sigma, eps);
            if x <> "Cannot Perform The Operation" then
      			Add(tree, [2 * tree[i][1] + 1, x]);
      			j := 2 * tree[i][1] + 1;
            fi;
            i := i + 1;
        od;
        return tree;
    end;

    array := [];
    if level < 0 then return array;
    elif level = 0 then
        for k in [1..NumberOfVertices(Q)] do
            vertex := Concatenation("(", String(k), ",", String(1), ")");
            Add(array, vertex);
            vertex := Concatenation("(", String(k), ",", String(-1), ")");
            Add(array, vertex);
        od;
        return array;
    fi;
    for k in [1..NumberOfVertices(Q)] do
        vertex := Concatenation("(", String(k), ",", String(1), ")");
        tree := making_tree(A, vertex, level);
        for i in [1..Length(tree)] do
            Add(array, tree[i][2]);
        od;
        vertex := Concatenation("(", String(k), ",", String(-1), ")");
        tree := making_tree(A, vertex, level);
        for i in [1..Length(tree)] do
            Add(array, tree[i][2]);
        od;
    od;
    return array;
end
);

InstallMethod( BandsLessThan,
"for stringalgebras and levels",
true,
[ IsAlgebra, IsPosInt ], 0,
function ( A, level )
    local Q, rho1, rho, vertices, arrows, i, onerel, tworel, str, q, b, j, k,
    kQ, gens, rel, arr, arr1, temp;

    if IsStringAlgebra(A) = false then
      Error("The input algebra is not a string algebra");
      return 0;
    fi;

    Q := QuiverOfPathAlgebra(A);
    rho1 := RelationsOfAlgebra(A);

    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);

    for i in [1..Length(vertices)] do
      if Length(String(vertices[i])) <> 2 or String(vertices[i])[1] <> 'v' or
        String(vertices[i])[2] <> String(i)[1] then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    for i in [1..Length(arrows)] do
      if Length(String(arrows[i])) <> 1 or String(arrows[i])[1] <> CharInt(i+96) then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    #rho := [];
    #for i in [1..Length(rho1)] do
    #    onerel := CoefficientsAndMagmaElements(rho1[i]);
    #    tworel := WalkOfPath(onerel[1]);
    #    str := Concatenation(List(tworel, a -> String(a)));
    #    Append(rho, [str]);
    #od;

    q := StringsLessThan(A, level);
    b := [];
    for i in [1..Length(q)] do
        if IsABand(A, q[i]) = true then
            Append(b, [q[i]]);
        fi;
    od;
    return b;
end
);

InstallMethod( BandRepresentativesLessThan,
"for stringalgebras and levels",
true,
[ IsAlgebra, IsPosInt ], 0,
function ( A, level )
    local Q, rho1, rho, vertices, arrows, onerel, tworel, str, b, tempo, br, i,
    j, k, arr, arr1, kQ, temp, gens, rel;

    if IsStringAlgebra(A) = false then
      Error("The input algebra is not a string algebra");
      return 0;
    fi;

    Q := QuiverOfPathAlgebra(A);
    rho1 := RelationsOfAlgebra(A);

    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);

    for i in [1..Length(vertices)] do
      if Length(String(vertices[i])) <> 2 or String(vertices[i])[1] <> 'v' or
        String(vertices[i])[2] <> String(i)[1] then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    for i in [1..Length(arrows)] do
      if Length(String(arrows[i])) <> 1 or String(arrows[i])[1] <> CharInt(i+96) then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    #rho := [];
    #for i in [1..Length(rho1)] do
    #    onerel := CoefficientsAndMagmaElements(rho1[i]);
    #    tworel := WalkOfPath(onerel[1]);
    #    str := Concatenation(List(tworel, a -> String(a)));
    #    Append(rho, [str]);
    #od;

    b := BandsLessThan(A,level);
    tempo := [];
    br := [];
    for i in [1..Length(b)] do
        if \in(b[i], tempo) = false and SIntChar(b[i][1]) >= 97 and
          SIntChar(b[i][Length(b[i])]) < 97 then
            Append(br, [b[i]]);
            Append(tempo, CyclicPermutationOfABand(b[i]));
        fi;
    od;
    return br;
end
);

InstallMethod( IsDomesticStringAlgebra,
"for stringalgebras",
true,
[ IsAlgebra], 0,
function( A )
    local Q, rho1, rho, vertices, arrows, onerel, tworel, str, k, kQ, gens, rel,
    arr1, arr, i, j, N, L, level, b, count, arr_quiver, len, temp;

    if IsStringAlgebra(A) = false then
      Error("The input algebra is not a string algebra");
      return 0;
    fi;

    Q := QuiverOfPathAlgebra(A);
    rho1 := RelationsOfAlgebra(A);

    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);

    for i in [1..Length(vertices)] do
      if Length(String(vertices[i])) <> 2 or String(vertices[i])[1] <> 'v' or
        String(vertices[i])[2] <> String(i)[1] then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    for i in [1..Length(arrows)] do
      if Length(String(arrows[i])) <> 1 or String(arrows[i])[1] <> CharInt(i+96) then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    rho := [];
    for i in [1..Length(rho1)] do
        onerel := CoefficientsAndMagmaElements(rho1[i]);
        tworel := WalkOfPath(onerel[1]);
        str := Concatenation(List(tworel, a -> String(a)));
        Append(rho, [str]);
    od;

    N := NumberOfJoints(A);
    L := MaximumLengthOfDirectString(A);
    level := 2 * L * N;
    b := BandsLessThan(A, level);

    count := 0;
    arr_quiver := [];
    len := Length(b);
    for i in [1..len] do
        for j in [1..len] do
            if i <> j then
                temp := Concatenation(b[i], b[j]);
                if IsValidString(A, temp) = true and temp[1] <> '(' then
                    count := count + 1;
                    Append(arr_quiver, [[j, i, Concatenation("a", String(count))]]);
                fi;
            fi;
        od;
    od;
    return IsAcyclicQuiver(Quiver(len, arr_quiver));
end
);

InstallMethod( BridgeQuiver,
"for stringalgebras",
true,
[ IsAlgebra], 0,
function( A )
    local vertices, arrows, Q, rho_temp, onerel, tworel, rho, arr, arr1, kQ,
    gens, rel, KMPSearch, Q2, WeakBridgeQuiver, BridgeQuiver1,
    BandFreeStrings, IsBandFreeString, wb, wb1, lambda, i, j, k,
    str, temp, array, array1, perm, Q1, num, l, N, L, level, q, b, br, tempo,
    count, arr_quiver, len;

    if IsStringAlgebra(A) = false then
      Error("The input algebra is not a string algebra");
      return 0;
    fi;

    Q := QuiverOfPathAlgebra(A);
    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);

    for i in [1..Length(vertices)] do
      if Length(String(vertices[i])) <> 2 or String(vertices[i])[1] <> 'v' or
        String(vertices[i])[2] <> String(i)[1] then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    for i in [1..Length(arrows)] do
      if Length(String(arrows[i])) <> 1 or String(arrows[i])[1] <> CharInt(i+96) then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    rho_temp := RelationsOfAlgebra(A);
    rho := [];
    for i in [1..Length(rho_temp)] do
        onerel := CoefficientsAndMagmaElements(rho_temp[i]);
        tworel := WalkOfPath(onerel[1]);
        str := Concatenation(List(tworel, a -> String(a)));
        Append(rho, [str]);
    od;

    IsBandFreeString := function(Q, rho, input_str)
        local i, j, k, temp;
        #var := IsValidString(Q,rho,input_str);
        #if var <> "Valid Positive Length String" then
        #  return var;
        #fi;
        if input_str[1] = '(' then return 1;
        fi;

        for i in [1..Length(input_str)] do
            for j in [i + 1..Length(input_str)] do
                if TargetOfArrow(Q, input_str[i]) = SourceOfArrow(Q, input_str[j]) then
                    temp := [];
                    for k in [i..j] do Add(temp,input_str[k]);
                    od;
                    if IsABand(A, temp) = true then return 0;
                    fi;
                fi;
            od;
        od;

        return 1;
    end;

    BandFreeStrings := function(Q, rho, v1, v2, q)
        local tree, tree1, i, j, k, array, count, vertex, x;
        array := [];
        if v1 > NumberOfVertices(Q) or v2 > NumberOfVertices(Q) then
            return array;
        fi;

        for i in [1..Length(q)] do
            if q[i][1] = '(' then
                if v1 = v2 and q[i][2] = String(v1)[1] then
                    Add(array, q[i]);
                fi;
            else
                if String(TargetOfArrow(Q, q[i][1])) = Concatenation("v", String(v2))
                  and String(SourceOfArrow(Q, q[i][Length(q[i])])) =
                  Concatenation("v", String(v1)) then
                    if IsBandFreeString(Q, rho, q[i]) = 1 then
                        Add(array, q[i]);
                    fi;
                fi;
            fi;
        od;
        return array;
    end;

    WeakBridgeQuiver := function(Q, rho, br, q)
        local i, j, k, arr_quiver, arr_quiver1, len, j1, j2, temp, Q1, wb, wb1;
        arr_quiver := [];
        arr_quiver1 := [];
        len := Length(br);
        for i in [1..len] do
            for j in [1..len] do
                #if i <> j then
                j1 := SIntChar(String(TargetOfArrow(Q,br[j][1]))[2]) - SIntChar('0');
                j2 := SIntChar(String(SourceOfArrow(Q,br[i][Length(br[i])]))[2]) - SIntChar('0');
                temp := BandFreeStrings(Q, rho, j1, j2, q);
                for k in [1..Length(temp)] do
                    if IsValidString(A, Concatenation(br[i], temp[k], br[j]))
                      = true and Concatenation(br[i], temp[k], br[j])[1] <> '(' then
                        Add(arr_quiver1, [j, i, temp[k]]);
                        Add(arr_quiver, [br[j], br[i], temp[k]]);
                    fi;
                od;
                #fi;
            od;
        od;
        wb := arr_quiver;
        wb1 := arr_quiver1;
        Q1 := Quiver(len, arr_quiver1);
        return [Q1, wb, wb1];
    end;

    KMPSearch := function(pat, txt)
        local M, N, lps, i, j;
        M := Length(pat);
        N := Length(txt);
        lps := computeLPSarray(pat);
        for i in [1..Length(lps)] do
            lps[i] := lps[i] + 1;
        od;
        i := 1;
        j := 1;
        while i <= N do
            if pat[j] = txt[i] then
                j := j + 1;
                i := i + 1;
            fi;
            if j=M + 1 then
                return i - j + 1;
            elif i <= N and pat[j] <> txt[i] then
                if j <> 1 then
                    j := lps[j - 1];
                else
                    i := i + 1;
                fi;
            fi;
        od;
        return 0;
    end;

    BridgeQuiver1 := function(Q,rho,br,wb,wb1)
        local lambda, i, j, k, str, temp, array, array1, perm, Q1, num, l;
        temp := [];
        lambda := [];
        for k in [1..Length(wb)] do
            lambda[k] := 1;
        od;
        for i in [1..Length(wb)] do
            for j in [1..Length(wb)] do
                if wb[i][1] = wb[j][2] then
                    str := Concatenation(wb[i][3], wb[j][3]);
                    if IsValidString(A, str) = true and str[1] = '(' then
                        if IsBandFreeString(Q, rho, str) = 1 then
                            for k in [1..Length(wb)] do
                                if str = wb[k][3] then
                                    lambda[k] := 0;
                                    break;
                                fi;
                            od;
                        else
                            perm := CyclicPermutationOfABand(wb[i][1]);
                            for k in [1..Length(perm)] do
                                temp := [];
                                num := KMPSearch(perm[k], str);
                                if num <> 0 then
                                    for l in [1..num - 1] do
                                        Add(temp, str[l]);
                                    od;
                                    for l in [num + Length(perm[k])..Length(str)] do
                                        Add(temp, str[l]);
                                    od;
                                    break;
                                fi;
                            od;

                            for k in [1..Length(wb)] do
                                if temp = wb[k][3] then
                                    lambda[k] := 0;
                                    break;
                                fi;
                            od;
                        fi;
                    fi;
                fi;
            od;
        od;

        array := [];
        array1 := [];
        for k in [1..Length(lambda)] do
            if lambda[k] = 1 then
                Add(array, wb[k]);
                Add(array1, wb1[k]);
            fi;
        od;
        Q1 := Quiver(Length(br), array1);
        return Q1;
    end;

    N := NumberOfJoints(A);
    L := MaximumLengthOfDirectString(A);
    level := 2 * L * N;
    q := StringsLessThan(A, level);

    b := [];
    for i in [1..Length(q)] do
        if IsABand(A, q[i]) = true then
            Append(b, [q[i]]);
        fi;
    od;

    tempo := [];
    br := [];
    for i in [1..Length(b)] do
        if \in(b[i], tempo) = false and SIntChar(b[i][1]) >= 97 and
          SIntChar(b[i][Length(b[i])]) < 97 then
            Append(br, [b[i]]);
            Append(tempo, CyclicPermutationOfABand(b[i]));
        fi;
    od;

    count := 0;
    arr_quiver := [];
    len := Length(b);
    for i in [1..len] do
        for j in [1..len] do
            if i <> j then
                temp := Concatenation(b[i], b[j]);
                if IsValidString(A, temp) = true and temp[1] <> '(' then
                    count := count + 1;
                    Append(arr_quiver, [[j, i, Concatenation("a", String(count))]]);
                fi;
            fi;
        od;
    od;

    if IsAcyclicQuiver(Quiver(len, arr_quiver)) = false then
        Error("The String Algebra represented by the input is not a Domestic String Algebra.");
        return 0;
    fi;

    Q1 := WeakBridgeQuiver(Q, rho, br, q);
    Q2 := BridgeQuiver1(Q, rho, br, Q1[2], Q1[3]);
    return Q2;
end
);
