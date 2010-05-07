# GAP Implementation
# This file was generated from 
# $Id: algtens.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallMethod(TensorProductOfAlgebras,
   "for a list of algebras",
   true,
   [IsDenseList],
   0,
   function(list)

   local basis,domain,fam,fin,g,genlist,gens,i,j,newgens,newlist,nextgens,O,one,tone,type;

   domain:=LeftActingDomain(list[1]);
   if not ForAll(list,x->LeftActingDomain(x)=domain) then
      Error("All algebras must be defined over the same domain.");
   fi;
   if not ForAll(list,x->IsAlgebra(x)) then
      Error("All parameters must be algebras.");
   fi;
   # flatten algebra tensors to component algebras

   newlist:=[];
   for i in list do
      if IsAlgebraTensor(i) then
         Append(newlist,ConstituentAlgebras(i));
      else
         Add(newlist,i);
      fi;
   od;
   list:=newlist;
   fin:=ForAll(list,x->IsFiniteDimensional(x));
   fam:=NewFamily("TensorElementsFam",IsAlgebraTensorElement);
   type:=NewType(fam,IsMonomialElementRep);
   tone:=Objectify(type,[[],not fin]);
   fam!.monomialElementDefaultType:=type;
   fam!.zeroCoefficient:=Zero(domain);
   # compute generating sets of each algebra or a basis if available

   if fin then
      fam!.constituentBases:=List(list,Basis);
      tone![1]:=[List(list,x->One(x)),One(domain)];

      # construct Cartesian product of generating sets

      gens:=Cartesian(fam!.constituentBases);
      Sort(gens);
      gens:=List(gens,x->Objectify(type,[[x,One(domain)],true]));
   else
      fam!.constituentBases:=[];
      tone![1]:=[[],One(domain)];
      for i in [1..Length(list)] do
         # compute one of space
         one:=One(list[i]);
         if not HasOne(list[i]) or one=fail then
            Error("Unable to compute one of algebra.");
         fi;
         Add(tone![1]![1],one);
         # build generating set
         if IsFiniteDimensional(list[i]) then
            fam!.constituentBases[i]:=BasisVectors(Basis(list[i]));
         else
            fam!.constituentBases[i]:=GeneratorsOfAlgebra(list[i]);
         fi;
         # ensure one in generating set
         if not one in fam!.constituentBases[i] then
            fam!.constituentBases[i]:=ShallowCopy(fam!.constituentBases[i]);
            Add(fam!.constituentBases[i],one);
         fi;
      od;

      # construct generating set using one generator and 1's

      gens:=[tone];
      for i in [1..Length(list)] do
         for g in fam!.constituentBases[i] do
            if not g=One(list[i]) then
               newgens:=ShallowCopy(tone![1]![1]);
               newgens![i]:=g;
               Add(gens,Objectify(type,[[newgens,One(domain)],true]));
            fi;
         od;
      od;
   fi;
   # construct algebra tensor

   O:=Algebra(domain,gens);
   SetFilterObj(O,IsAlgebraTensor);
   SetOne(O,tone);
   SetIsFiniteDimensional(O,fin);
   O!.constituentAlgebras:=list;
   SetElementsFamily(O,fam);
   fam!.algebraTensor:=O;

   # compute and store basis for finite algebra tensors

   if fin then
      basis:=BasisOfMonomialSpace(O,gens);
      basis!.echelonBasis:=gens;
      basis!.baseChange:=List([1..Length(gens)],x->[[x,One(domain)]]);
      basis!.heads:=List(gens,x->ExtRepOfObj(x)[1]);
      basis!.zeroCoefficient:=Zero(domain);
      SetBasis(O,basis);
   fi;
   return O;
end);
InstallOtherMethod(TensorProductOfAlgebras,
   "for two algebras",
   true,
   [IsAlgebra,IsAlgebra],
   0,
   function(X,Y)
 
   return TensorProductOfAlgebras([X,Y]);
end);
InstallMethod(\in,
   "for tensor products of algebras",
   IsElmsColls,
   [ IsAlgebraTensorElement, IsAlgebraTensor], 0,
   function( t, T )
       local algebras, rep, tElms, aElms, i;

       algebras := ConstituentAlgebras(T);
       rep := t![1];
       # match against 0 element
       if Length(rep) = 0 then
           return true;
       fi;

       # get the constituent elements
       tElms := rep{[1,3..Length(rep)-1]};
       for i in [1..Length(algebras)] do
           # select elements from ith algebra
           aElms := List(tElms, x -> x[i]);
           if not IsSubset(algebras[i], aElms) then
               return false;
           fi;
       od;
       return true;
   end );
InstallMethod(TensorElement,
   "for an algebra tensor and list of algebra elements or lists of algebra elements",
   true,
   [IsAlgebraTensor,IsDenseList],
   0,
   function(T,list)

   local algebras,i,elem,j,k,len,newlist,nextlist,O,templist;

   algebras:=ConstituentAlgebras(T);
   if IsList(list[1]) then
      elem:=Zero(T);
      # apply to each summand
      for i in [1..Length(list)] do
         elem:=elem+TensorElement(T,list[i]);
      od;
      return elem;
   fi;
   # flatten algebra tensor elements to component algebra elements

   newlist:=[[]];
   for i in [1..Length(list)] do
      len:=Length(newlist);
      if IsAlgebraTensorElement(list[i]) then
         # take apart existing tensor element
         elem:=ConstituentElements(list[i]);
         if elem=fail then
            Error("Unable to decompose existing tensor element.");
         fi;
         nextlist:=[];
         # add new element to each part of existing tensor element
         for j in [1..Length(newlist)] do
            for k in [1..Length(elem)] do 
               templist:=ShallowCopy(newlist[j]);
               Append(templist,elem[k]);
               Add(nextlist,templist);
            od;             
         od;
         newlist:=nextlist;
      else
         # add new element to each part
         for j in [1..Length(newlist)] do
            Add(newlist[j],list[i]);
         od;
      fi;
   od;
   # construct tensor element as a sum of tensor products of algebra elements

   elem:=Zero(T);
   for i in [1..Length(newlist)] do
      list:=newlist[i];
      len:=Length(list);
      if not len=Length(algebras) then
         Error("Must have the same number of components as algebra tensor.");
      fi;
      if not ForAll([1..len],i->list[i] in algebras [i]) then
         Error("Algebra elements must be contained in tensor domain.");
      fi;
      O:=Objectify(ElementsFamily(T)!.monomialElementDefaultType,[[list,One(LeftActingDomain(T))],false]);
      # attach summand
      elem:=elem+O;
   od;
   return elem;
end);
InstallOtherMethod(TensorElement,
   "for an algebra tensor and two algebra elements",
   true,
   [IsAlgebraTensor,IsRingElement,IsRingElement],
   0,
   function(T,e1,e2)

   return TensorElement(T,[e1,e2]);
end);
InstallMethod(ConstituentElements,
   "for an algebra tensor element",
   true,
   [IsAlgebraTensorElement],
   0,
   function(elem)

   local coeff,decomp,i,part;

   decomp:=[];
   # walk along non-coefficients
   for i in [1,3..Length(elem![1])-1] do
      part:=ShallowCopy(elem![1]![i]);
      # replace missing terms with 0
      if Length(part)=0 then
         part:=List(ConstituentAlgebras(FamilyObj(elem)!.algebraTensor),x->Zero(x));
      fi;
      part![1]:=part![1]*elem![1]![i+1];
      Add(decomp,part);
   od;
   return decomp;
end);
InstallMethod(ConstituentAlgebras,
   "for an algebra tensor",
   true,
   [IsAlgebraTensor],
   0,
   function(T)

   return ShallowCopy(T!.constituentAlgebras);
end);
InstallMethod(\*,
   "for algebra tensor elements",
   IsIdenticalObj,
   [IsAlgebraTensorElement,IsAlgebraTensorElement],
   0,
   function(X,Y)

   local f,fam,prod,zero;

   # function for computing inner product of two vectors

   f:=function(X,Y)
      local i,Z;

      Z:=[];
      for i in [1..Length(X)] do
         Z[i]:=X[i]*Y[i];
      od;
      return Z;
   end;
   fam:=FamilyObj(X);
   zero:=fam!.zeroCoefficient;
   # use internal GAP magic
   prod:=ZippedProduct(X![1],Y![1],zero,[f,\<,\+,\*]);
   if prod=[] then
      # replace empty list with 0
      return Objectify(fam!.monomialElementDefaultType,[[[],zero],true]);
   fi;
   return SimplifyTensorElement(Objectify(fam!.monomialElementDefaultType,[prod,false]));
end);
InstallMethod(TensorProductOfAlgebraMaps,
   "for a list of algebra maps",
   true,
   [IsDenseList],
   0,
   function(list)

   local algebra,elem,elems,gens,i,j,k,len,range,ranget,source,sourcet;

   if not ForAll(list,x->IsAlgebraGeneralMapping(x)) then
      Error("List of algebra maps required.");
   fi;
   source:=List([1..Length(list)],x->Source(list[x]));
   range:=List([1..Length(list)],x->Range(list[x]));
   # source is tensor product of sources
   sourcet:=TensorProductOfAlgebras(source);
   # range is tensor product of ranges
   ranget:=TensorProductOfAlgebras(range);
   # follow maps componentwise by generators
   gens:=GeneratorsOfAlgebra(sourcet);
   elems:=[];
   for i in [1..Length(gens)] do
      Add(elems,[]);
      k:=1;
      # apply map to each component
      for j in [1..Length(list)] do
         if IsAlgebraTensor(source[j]) then
            # follow down every map in tensor product
            len:=Length(ConstituentAlgebras(source[i]));
            elem:=TensorElement(source[j],gens[i]![1]![1]{[k..k+len]});
         else
            # follow down map
            len:=1;
            elem:=gens[i]![1]![1]![k];
         fi;
         # add new image
         Add(elems[i],Image(list[j],elem));
         k:=k+len;
      od;
      # put together images
      elems[i]:=TensorElement(ranget,elems[i]);
   od;
   if ForAll(list,x->IsAlgebraHomomorphism(x)) then
      return AlgebraHomomorphismByImagesNC(sourcet,ranget,gens,elems);
   else
      return AlgebraGeneralMappingByImages(sourcet,ranget,gens,elems);
   fi;
end);
InstallOtherMethod(TensorProductOfAlgebraMaps,
   "for two algebra maps",
   true,
   [IsAlgebraGeneralMapping,IsAlgebraGeneralMapping],
   0,
   function(X,Y)

   return TensorProductOfAlgebraMaps([X,Y]);
end);
InstallMethod(SimplifyTensorElement,
   "for an algebra tensor element",
   true,
   [IsAlgebraTensorElement],
   0,
   function(elem)

   local azeros,fam,i,simple;

   fam:=FamilyObj(elem);
   # construct zero
   azeros:=List(ConstituentAlgebras(fam!.algebraTensor),x->Zero(x));
   simple:=[];
   # walk along non-coefficients
   for i in [1,3..Length(elem![1])-1] do
      # toss zero summands
      if ForAll([1..Length(elem![1]![i])],x->not elem![1]![i]![x]=azeros[x]) then
         Append(simple,[ShallowCopy(elem![1]![i]),elem![1]![i+1]]);
      fi;      
   od;
   if simple=[] then
      # replace empty list with 0
      return Objectify(fam!.monomialElementDefaultType,[[[],fam!.zeroCoefficient],true]);
   fi;
   return Objectify(fam!.monomialElementDefaultType,[simple,false]);
end);
InstallMethod(ViewObj,
   "for tensor products of algebras",
   true,
   [IsAlgebraTensor],
   1000,
   function(a)

   Print("<tensor product of ",Length(ConstituentAlgebras(a))," algebras, with ",Length(GeneratorsOfAlgebra(a))," generators>");   
end);
