function theta_int(z,M,n,tau)
e1:=M[1];
f1:=M[2];
C:=BaseRing(tau);
I:=C.1;
pi:=Pi(C);
g:=#z;
S:=0;
e:=ChangeRing(e1,C);
f:=ChangeRing(f1,C);
for i:=1 to g do
for j:=1 to g do
S:=S+ (n[i]+e[i]/2)*tau[i,j]*(n[j]+e[j]/2);
end for;
end for;
S:=1/2*S;
for i:=1 to g do
S:=S+(n[i]+e[i]/2)*(z[i]+f[i]/2);
end for;
return Exp(2*I*pi*S);
end function;

function theta(z,M,tau,prec)
g:=#z;
L:=CartesianPower([-prec..prec],g);
S:=0;
for n in L do
S:=S+theta_int(z,M,n,tau);
end for;
return S;
end function;

function InitializeOddChar()
	OddChar := [];
	for m11 in [0..1] do
	for m12 in [0..1] do
	for m21 in [0..1] do
	for m22 in [0..1] do
		parity := m11*m21+m12*m22;
		if parity mod 2 eq 1 then
			Append(~OddChar, Matrix(Integers(),2,2,[[m11, m12], [m21, m22]]) );
		end if;
	end for;
	end for;
	end for;
	end for;
	return OddChar;
end function;

function BracketProduct(tau, m1, m2)
	OddChar := InitializeOddChar();
	print OddChar;
	prec := 10;
	result := 1;
	for m in OddChar do
		if m ne m1 and m ne m2 then
			print m1+m2-m;
			result := result * theta([0,0],m1+m2-m,tau,prec);
		end if;
	end for;
	return result;
end function;

function SymmetricModel(tau)
	OddChar := InitializeOddChar();
	l123 := BracketProduct(tau, OddChar[1], OddChar[3]) / BracketProduct(tau, OddChar[2], OddChar[3]);
	l124 := BracketProduct(tau, OddChar[1], OddChar[4]) / BracketProduct(tau, OddChar[2], OddChar[4]);
	l125 := BracketProduct(tau, OddChar[1], OddChar[5]) / BracketProduct(tau, OddChar[2], OddChar[5]);
	l126 := BracketProduct(tau, OddChar[1], OddChar[6]) / BracketProduct(tau, OddChar[2], OddChar[6]);

	MinimalPolynomial(l123,12);

	R<x> := PolynomialRing(Parent(l123));
	f := x*(x-l123)*(x-l124)*(x-l125)*(x-l126);
	C := HyperellipticCurve(f);
	return C;
end function;

SymmetricModel(tau);




> CReconstructed;
Hyperelliptic Curve defined by y^2 = x^5 +
    (4.051840519856159834215128450945191339334516733980945729971762139882586061\
    081228080566873366982299857 + 0.0948029265796288672597342519306711595299127\
    3389642605668199782354732783655411642418587038781265321353*I)*x^4 +
    (6.098268488322249148437103256271907084928773091944057430685561297229004204\
    277907510412786234022325713 + 0.3666227623388291106884159617041902017518818\
    441910635919498633630002339810936747191596854391576472919*Cpx.1)*x^3 +
    (4.072078333452894616947527776379381897801661118022906776343356732143685106\
    391708053637954621532826944 + 0.2588690433225026031026683800495259605243749\
    633892199668053707404288728095688717380207969392653544270*Cpx.1)*x^2 +
    (0.999999999999999999999999999999999999999999999999999999999999999999999999\
    9999999999999999999999999999 + 2.567304091762875015929360357865547987598465\
    956139521966525183645088307529898631117365753121153544740E-100*Cpx.1)*x over
Cpx
>



function test();
Cpx := ComplexField(200);
R<x> := PolynomialRing(Cpx);
f := x^5 + 3*x^2 + 4*x +1;
C := HyperellipticCurve(f);
J := AnalyticJacobian(C);
tau := SmallPeriodMatrix(J);
roots := RosenhainInvariants(tau);
fRec := &* [(x-r) : r in roots];
CRec := HyperellipticCurve(x*(x-1)*fRec);
CRecSym := SymmetricModel(tau);
g2 := IgusaInvariants(C);
g2Rec := IgusaInvariants(CRec);
g2Sym := IgusaInvariants(CRecSym);
return g2, g2Rec, g2Sym;
end function;





CRec := SymmetricModel(tau);
CRec;





for d in [1..3] do
	MinimalPolynomial(Coefficients(fRec)[d], 20);
end for;




5.9512249446856273722296865167965987197060027075881794366728021624556881775440967098598091890523874380862132869609180773086541874833454498809800403530221973665733718856982380112523884110165782491077135 -1.4626333142152073125803412408881664408669752264548521323325303506154163555058763453688998497771110862332325882955482618686778859390689382594650638668436675044078849621248278257937313085696395513018241*I