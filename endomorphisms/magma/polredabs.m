intrinsic Polred(f::RngUPolElt : Best:="dynamic") -> RngUPolElt, SeqEnum, BoolElt
  { A smallest generating polynomial of the number field, using polredabs/polredbest in pari. }

  require &or[ Best cmpeq elt : elt in [* "dynamic", false, true *]]: "Best must be in [\"dynamic\", false, true]";
  if Best cmpeq "dynamic" then
    Best := Degree(f) gt 8;
  end if;

  vprint EndoFind, 3 : "Starting Polred...";
  if Best then
    cmdp := "polredbest";
  else
    cmdp := "polredabs";
  end if;

  try
    cmd := Sprintf(
       "{u = %o(Pol(Vecrev(%o)),1); print(Vecrev(Vec(u[1])),Vecrev(Vec(lift(u[2]))))}",
       cmdp, Coefficients(f));
    s := Pipe("gp -q -D timer=0", cmd);
    c := Index(s,"][");
    spol := s[1..c];
    sroot := s[c+1..#s-1];
    sspol := [ StringToInteger(x) : x in Split(spol, ", []\n") | x ne "" ];
    ssroot := [ StringToRational(x) : x in Split(sroot, ", []\n") | x ne "" ];
    assert #ssroot le Degree(f);
    ssroot := ssroot cat [0 : i in [1..Degree(f)-#ssroot]];
  catch e
    print("WARNING: need gp at command-line for polredabs, without this many examples become intractable\n");
    vprint EndoFind, 3 : "done.";
    assert Type(BaseRing(Parent(f))) eq FldRat;
    lcm := LCM([ Denominator(c) : c in Coefficients(f) ]);
    f0 := Evaluate(f, (Parent(f).1)/lcm);
    f0 /:= LeadingCoefficient(f0);
    return f0, [0,1] cat [0: i in [1..Degree(f)-2]], false;
  end try;
  vprint EndoFind, 3 : "Polred done.";
  return Parent(f) ! sspol, ssroot, true;
end intrinsic

intrinsic Polredabs(f::RngUPolElt) -> RngUPolElt, SeqEnum, BoolElt
  { A smallest generating polynomial of the number field, using polredabs in pari. }
  return Polred(f : Best:=false);
end intrinsic;


intrinsic Polredbest(f::RngUPolElt) -> RngUPolElt, SeqEnum, BoolElt
  { A smallest generating polynomial of the number field, using polredbest in pari. }
  return Polred(f : Best:=true);
end intrinsic;


intrinsic Polredbestabs(f::RngUPolElt) -> RngUPolElt, SeqEnum, BoolElt
  {A smallest generating polynomial of the number field, using pari.  First polredbest, then polredabs.}

  K := (Degree(f) ne 1) select NumberField(f) else RationalsAsNumberField();
  fbest, fbest_root := Polredbest(f);
  fredabs, fredabs_root, bl := Polredabs(fbest);
  Kbest := (Degree(fbest) ne 1) select NumberField(fbest) else RationalsAsNumberField();
  iotabest := hom<K -> Kbest | fbest_root>;
  Kredabs := (Degree(fredabs) ne 1) select NumberField(fredabs) else RationalsAsNumberField();
  iotaredabs := hom<Kbest -> Kredabs | fredabs_root>;
  iota := iotabest*iotaredabs;  // functional composition is backwards in Magma, for some reason
  return fredabs, Eltseq(iota(K.1)), bl;
end intrinsic;

intrinsic Polred(K::Fld : Best:="dynamic") -> FldNum, Map, BoolElt
  { A smallest generating polynomial of the number field, using polredabs/polredbest pari. }

  if Type(K) eq FldRat then
    return K, hom< K -> K | >;
  else
    fredabs, fredabs_root, bl := Polredabs(DefiningPolynomial(K) : Best:=Best);
    if IsZero(fredabs) then
        return K, hom<K -> K | K.1>, false;
    end if;
    Kredabs := NumberField(fredabs);
    iotaredabs := hom<K -> Kredabs | fredabs_root>;
    return Kredabs, iotaredabs, bl;
  end if;
end intrinsic;

intrinsic Polredbest(K::Fld) -> FldNum, Map, BoolElt
  { A smallest generating polynomial of the number field, using polredabs in pari. }
  return Polred(K : Best:=false);
end intrinsic;


intrinsic Polredbest(K::Fld) -> FldNum, Map, BoolElt
  { A smallest generating polynomial of the number field, using polredbest in pari. }
  return Polred(K : Best:=true);
end intrinsic;

intrinsic Polredbestabs(K::Fld) -> RngUPolElt, Map, BoolElt
  {A smallest generating polynomial of the number field, using pari.  First polredbest, then polredabs.}

  if Type(K) eq FldRat then
    return K, hom< K -> K | >;
  else
    f := DefiningPolynomial(K);
    fbest, fbest_root, bl0 := Polredbest(f);
    if IsZero(fbest) then
        return K, hom<K -> K | K.1>, false;
    end if;
    fredabs, fredabs_root, bl1 := Polredabs(fbest);
    if IsZero(fredabs) then
        return K, hom<K -> K | K.1>, false;
    end if;
    assert bl0 eq bl1;

    Kbest := NumberField(fbest);
    iotabest := hom<K -> Kbest | fbest_root>;
    Kredabs := NumberField(fredabs);
    iotaredabs := hom<Kbest -> Kredabs | fredabs_root>;
    iota := iotabest*iotaredabs;  // functional composition is backwards in Magma, for some reason
    return Kredabs, iota, bl0;
  end if;
end intrinsic;

intrinsic Eval(alpha::FldNumElt, v::PlcNumElt, is_conjugated::BoolElt) -> FldComElt
  {}
  z := Evaluate(alpha, v);
  if is_conjugated then
    return ComplexConjugate(z);
  else
    return z;
  end if;
end intrinsic;

intrinsic Eval(l::List) -> FldComElt
  {}
  return Eval(l[1].1, l[2], l[3]);
end intrinsic;

intrinsic MachineZero(z::FldReElt) -> BoolElt
  {}
  prec := Precision(Parent(z));
  if AbsoluteValue(z) le 10^(-prec/2) then
    return true;
  else
    return false;
  end if;
end intrinsic;

intrinsic MachineZero(z::FldComElt) -> BoolElt
  {}
  prec := Precision(Parent(z));
  if AbsoluteValue(z) le 10^(-prec/2) then
    return true;
  else
    return false;
  end if;
end intrinsic;
