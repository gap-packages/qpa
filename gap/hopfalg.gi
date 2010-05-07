# GAP Implementation
# This file was generated from 
# $Id: hopfalg.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
InstallMethod(Comultiply,
      "for a Hopf algebra element",
      true,
      [IsRingElement],
      0,
      function(e)
   return e^ComultiplicationMap(FamilyObj(e)!.hopfalgebra);
end);
InstallMethod(Counit,
      "for a Hopf algebra element",
      true,
      [IsRingElement],
      0,
      function(e)
   return e^CounitMap(FamilyObj(e)!.hopfalgebra);
end);
InstallMethod(Antipode,
      "for a Hopf algebra element",
      true,
      [IsRingElement],
      0,
      function(e)
   return e^AntipodeMap(FamilyObj(e)!.hopfalgebra);
end);
InstallMethod(HopfAlgebra,
   "for a group ring",
   true,
   [IsGroupRing],
   0,
   function(gr)

   if IsField(LeftActingDomain(gr)) then
      return HopfAlgebraNC(gr);
   else
      Error("Left acting domain of group ring is not a field.");
   fi;
end);

InstallMethod(HopfAlgebraNC,
   "for a group ring",
   true,
   [IsGroupRing],
   0,
   function(gr)

   local b,d,o,t;

   t:=TensorProductOfAlgebras(gr,gr);
   b:=BasisVectors(Basis(gr));
   d:=LeftActingDomain(gr);
   o:=One(d);
   # comultiplication: elem goes to elem tensor elem
   # counit: elem goes to 1
   # antipode: elem goes to group inverse
   return HopfAlgebraNC(gr,AlgebraHomomorphismByImagesNC(gr,t,b,List(b,x->TensorElement(t,x,x))),AlgebraHomomorphismByImagesNC(gr,d,b,List(b,x->o)),AlgebraHomomorphismByImagesNC(gr,gr,b,List(b,x->x^-1)));
end);


InstallMethod(HopfAlgebraNC,
   "for an algebra and algebra maps",
   true,
   [IsAlgebra,IsAlgebraGeneralMapping,IsAlgebraGeneralMapping,IsAlgebraGeneralMapping],
   function(a,comult,counit,antipode)

   local fam,hopf,type;

   type:=NewType(NewFamily("HopfAlgebraElementsFamily",IsHopfAlgebraElement),IsPackedElementDefaultRep);
   if IsAlgebraWithOne(a) then
      hopf:=AlgebraWithOneByGenerators(LeftActingDomain(a),List(GeneratorsOfAlgebraWithOne(a),x->Objectify(type,[x])));
      SetOne(hopf,Objectify(type,[One(a)]));
   else
      hopf:=AlgebraByGenerators(LeftActingDomain(a),List(GeneratorsOfAlgebra(a),x->Objectify(type,[x])));
   fi;
   SetComultiplicationMap(hopf,comult);
   SetCounitMap(hopf,counit);
   SetAntipodeMap(hopf,antipode);
   SetUnderlyingAlgebra(hopf,a);
   SetFilterObj(hopf,IsHopfAlgebra);
   fam:=ElementsFamily(FamilyObj(hopf));
   fam!.packedType:=type;
   fam!.underlyingAlgebraEltsFam:=ElementsFamily(FamilyObj(a));
   fam!.hopfalgebra:=hopf;
   return hopf;
end);


InstallMethod(ViewObj,
   "for a Hopf algebra",
   true,
   [IsHopfAlgebra],
   1000,
   function(a)

   Print("<Hopf algebra from ");
   ViewObj(UnderlyingAlgebra(a));
   Print(">");
end);

InstallMethod(\in,
   "for a Hopf algebra element and Hopf algebra",
   true,
   [IsHopfAlgebraElement,IsHopfAlgebra],
   function(e,a)

   return e![1] in UnderlyingAlgebra(a);
end);

InstallMethod(Dimension,
      "for a Hopf algebra",
      true,
      [IsHopfAlgebra],
      function(a)
   return Dimension(UnderlyingAlgebra(a));
end);

InstallMethod(IsFiniteDimensional,
      "for a Hopf algebra",
      true,
      [IsHopfAlgebra],
      function(a)
   return IsFiniteDimensional(UnderlyingAlgebra(a));
end);

InstallMethod(Basis,
         "for a Hopf algebra",
         true,
         [IsHopfAlgebra],
         function(a)
      local b;

      b:=Objectify(NewType(FamilyObj(a),IsHopfAlgebraBasis and IsAttributeStoringRep),rec());;
      SetUnderlyingLeftModule(b,a);
   b!.underlyingBasis:=Basis(UnderlyingAlgebra(a));
   return b;
end);

InstallMethod(CanonicalBasis,
         "for a Hopf algebra",
         true,
         [IsHopfAlgebra],
         function(a)
      local b;

      b:=Objectify(NewType(FamilyObj(a),IsHopfAlgebraBasis and IsAttributeStoringRep),rec());;
      SetUnderlyingLeftModule(b,a);
   b!.underlyingBasis:=CanonicalBasis(UnderlyingAlgebra(a));
   if b!.underlyingBasis=fail then
      return fail;
   fi;
   SetIsCanonicalBasis(b,true);
   return b;
end);

InstallMethod( Coefficients,
   "for a Hopf algebra basis and a Hopf algebra element",
   true,
   [IsHopfAlgebraBasis,IsHopfAlgebraElement and IsPackedElementDefaultRep],
   function(b,e)

   return Coefficients(b!.underlyingBasis,e![1]);
end);

InstallMethod(BasisVectors,
   "for a Hopf algebra basis",
   true,
   [IsHopfAlgebraBasis],
   function(b)

   local type,vectors;

   vectors:=BasisVectors(b!.underlyingBasis);
   type:=ElementsFamily(FamilyObj(b))!.packedType;
   return List(vectors,x->Objectify(type,[x]));
end);
InstallMethod(\=,
      "for two Hopf algebra elements",
      IsIdenticalObj,
      [IsHopfAlgebraElement and IsPackedElementDefaultRep,IsHopfAlgebraElement and IsPackedElementDefaultRep],
      function(e1,e2)
   return e1![1]=e2![1];
end);

InstallMethod(\<,
      "for two Hopf algebra elements",
      IsIdenticalObj,
      [IsHopfAlgebraElement and IsPackedElementDefaultRep,IsHopfAlgebraElement and IsPackedElementDefaultRep],
      function(e1,e2)
   return e1![1]<e2![1];
end);

InstallMethod(\+,
      "for two Hopf algebra elements",
      IsIdenticalObj,
      [IsHopfAlgebraElement and IsPackedElementDefaultRep,IsHopfAlgebraElement and IsPackedElementDefaultRep],
      function(e1,e2)
   return Objectify(TypeObj(e1),[e1![1]+e2![1]]);
end);

InstallMethod(\*,
      "for two Hopf algebra elements",
      IsIdenticalObj,
      [IsHopfAlgebraElement and IsPackedElementDefaultRep,IsHopfAlgebraElement and IsPackedElementDefaultRep],
      function(e1,e2)
   return Objectify(TypeObj(e1),[e1![1]*e2![1]]);
end);

InstallMethod(\*,
   "for a Hopf algebra element and scalar",
   true,
   [IsHopfAlgebraElement and IsPackedElementDefaultRep,IsScalar],
   function(e,s)
   
   return Objectify(TypeObj(e),[e![1]*s]);
end);

InstallMethod(\*,
   "for a scalar and Hopf algebra element",
   true,
   [IsScalar,IsHopfAlgebraElement and IsPackedElementDefaultRep],
   function(s,e)

   return Objectify(TypeObj(e),[s*e![1]]);
end);

InstallMethod(AINV,
      "for a Hopf algebra element",
      true,
      [IsHopfAlgebraElement],
      function(e)
   return Objectify(TypeObj(e),[-e![1]]);
end);

InstallMethod(\^,
   "for a Hopf algebra element and general algebra mapping",
   true,
   [IsHopfAlgebraElement,IsAlgebraGeneralMapping],
   function(e,map)

   return e![1]^map;
end);

InstallMethod(ZeroOp,
      "for a Hopf algebra element",
      true,
      [IsHopfAlgebraElement],
      function(e)
   return Objectify(TypeObj(e),[0*e![1]]);
end);

InstallMethod(PrintObj,
   "for a Hopf algebra element",
   true,
   [IsHopfAlgebraElement],
   1000,
   function(e)

   Print("Hopf(",e![1],")");
end);
   
InstallMethod(ViewObj,
   "for a Hopf algebra element",
   true,
   [IsHopfAlgebraElement],
   1000,
   function(e)

   ViewObj(e![1]);
end);
InstallMethod(AlgebraAntihomomorphismByImages,
   "for an algebra, an algebra, a list of generators, and list of images",
   [IsAlgebra,IsAlgebra,IsDenseList,IsDenseList],
   0,
   function(A,B,gens,imgs)

   local hom,opB,opElemFamily,opImgs,map;

   # treat as homomorphism with opposite algebra
   opB:=OppositeAlgebra(B);
   opElemFamily:=ElementsFamily(FamilyObj(opB));
   opImgs:=List(imgs,x->Objectify(opElemFamily!.packedType,[x]));
   hom:=AlgebraHomomorphismByImages(A,opB,gens,opImgs);
   if hom=fail then
      Error("Unable to construct underlying homomorphism");
   fi;
   # take map to be unwrapped version of opposite element
   map:=function(e)
      local opImg;

      opImg:=Image(hom,e);
      return opImg![1];
   end;
   map:=MappingByFunction(A,B,map);
   SetFilterObj(map,IsAlgebraAntihomomorphism);
   return map;
end);
InstallMethod(AlgebraAntihomomorphismByImagesNC,
   "for an algebra, an algebra, a list of generators, and list of images",
   [IsAlgebra,IsAlgebra,IsDenseList,IsDenseList],
   0,
   function(A,B,gens,imgs)

   local hom,opB,opElemFamily,opImgs,map;

   # treat as homomorphism with opposite algebra
   opB:=OppositeAlgebra(B);
   opElemFamily:=ElementsFamily(FamilyObj(opB));
   opImgs:=List(imgs,x->Objectify(opElemFamily!.packedType,[x]));
   hom:=AlgebraHomomorphismByImagesNC(A,opB,gens,opImgs);
   if hom=fail then
      Error("Unable to construct underlying homomorphism");
   fi;
   # take map to be unwrapped version of opposite element
   map:=function(e)
      local opImg;

      opImg:=Image(hom,e);
      return opImg![1];
   end;
   map:=MappingByFunction(A,B,map);
   SetFilterObj(map,IsAlgebraAntihomomorphism);
   return map;
end);
InstallMethod(AlgebraWithOneAntihomomorphismByImages,
   "for an algebra-with-one, an algebra-with-one, a list of generators, and list of images",
   [IsAlgebraWithOne,IsAlgebraWithOne,IsDenseList,IsDenseList],
   0,
   function(A,B,gens,imgs)

   local hom,opB,opElemFamily,opImgs,map;

   # treat as homomorphism with opposite algebra
   opB:=OppositeAlgebra(B);
   opElemFamily:=ElementsFamily(FamilyObj(opB));
   opImgs:=List(imgs,x->Objectify(opElemFamily!.packedType,[x]));
   hom:=AlgebraWithOneHomomorphismByImages(A,opB,gens,opImgs);
   if hom=fail then
      Error("Unable to construct underlying homomorphism");
   fi;
   # take map to be unwrapped version of opposite element
   map:=function(e)
      local opImg;

      opImg:=Image(hom,e);
      return opImg![1];
   end;
   map:=MappingByFunction(A,B,map);
   SetFilterObj(map,IsAlgebraWithOneAntihomomorphism);
   return map;
end);
InstallMethod(AlgebraWithOneAntihomomorphismByImagesNC,
   "for an algebra-with-one, an algebra-with-one, a list of generators, and list of images",
   [IsAlgebraWithOne,IsAlgebraWithOne,IsDenseList,IsDenseList],
   0,
   function(A,B,gens,imgs)

   local hom,opB,opElemFamily,opImgs,map;

   # treat as homomorphism with opposite algebra
   opB:=OppositeAlgebra(B);
   opElemFamily:=ElementsFamily(FamilyObj(opB));
   opImgs:=List(imgs,x->Objectify(opElemFamily!.packedType,[x]));
   hom:=AlgebraWithOneHomomorphismByImages(A,opB,gens,opImgs);
   if hom=fail then
      Error("Unable to construct underlying homomorphism");
   fi;
   # take map to be unwrapped version of opposite element
   map:=function(e)
      local opImg;

      opImg:=Image(hom,e);
      return opImg![1];
   end;
   map:=MappingByFunction(A,B,map);
   SetFilterObj(map,IsAlgebraWithOneAntihomomorphism);
   return map;
end);
