#  Copyright (C) 2016-2017 Nicolas Mascot, Edgar Costa
#  See LICENSE file for license details.
#
#
#To use it, first load the code, then do
#
#SetPariPrec(500) // Set the decimal precision used in Pari (replace 500 with what you like). /!\ NOTE : 500 is the number of significant DECIMAL digits, not bits.
#
#PariIntNumInit(3); // Initialise an integration mesh in Pari. This takes a bit of time but only needs to be done once (after the decimal precision has been set of course). Don't forget the semicolon in the end, or Sage's interpreter will spit out the whole mess (sorry I haven't found how to avoid that). (Seems to be OK now; JRS)
#The argument (3 here) is optional. Giving n means multiply the number of integration points by 2^n (and thus make the integrations 2^n times slower).  Giving nothing means n=0. In practise, n=2 is numerically stable up to thousands of decimals I think.
#
#None of these commands produce any return, they just initialise Pari. Then to evaluate AJ, do
#
#AJ2(f,(0,1),(1,sqrt(2)))
#
#that is to say, the polynomial f(x) such that y^2=f(x) (must have degree 6 for now), and a couple of points. The output is a sage couple of Pari complexes, representing the integral of dx/y and x dx/y from infty_+ + infty_- to P+Q.
#

from sage.all import ComplexField, RR
from sage.all import ceil, gp, log, vector

#def load_gp():
#    for i in range(10):
#        try:
#            assert 2 == gp("1+1");
#            load_gp_raw();
#            assert 2 == gp("1+1");
#            print "gp code loaded";
#            break;
#        except TypeError:
#            print "failed, retrying"
#            pass;



def load_gp():
    if  not all([gp("RNEG") == 101, gp("RPOS") == 102, gp("INEG") == 103, gp("IPOS") == 104]):
        #print "loading gp code"
        gp._reset_expect();
        # Branch values:
        gp("RNEG=101; RPOS=102; INEG=103; IPOS=104;");

        # Print branch label given a branch value.
        # The label corresponds to the branch cut.
        gp(
            """
            printbranch(branch)=
            {
                if( branch == RNEG, print("RNEG"));
                if( branch == RPOS, print("RPOS"));
                if( branch == INEG, print("INEG"));
                if( branch == IPOS, print("IPOS"));
                return();
                error("Unknown sqrt branch");
            }
            """.replace('\n','')
        );

        # Square root determined by branch value, which should be thought of as a
        # continuation to a half-plane.
        gp(
            """
            rac2(z, branch)=
            {
                if( branch == RNEG, return( sqrt(z) ) );
                if( branch == RPOS, return( I * sqrt(-z) ) );
                if( branch == INEG, return( sqrt(I) * sqrt(-I*z) ) );
                if( branch == IPOS, return( sqrt(-I) * sqrt(I * z) ) );
                error("Unknown sqrt branch");
            }
            """.replace('\n','')
        );

        # Best branch is the one that stays farthest away from the cut corresponding to
        # it.
        gp(
            """
            best_branch(z)=
            {
                my( x = real(z), y = imag(z) );
                if( abs(y) > abs(x),
                    if( y > 0, return(INEG), return(IPOS))
                    ,
                    if( x > 0, return(RNEG), return(RPOS))
                );
            }
            """.replace('\n','')
        );

        # Evaluation:
        gp("EvPol(f, x) = subst(f, variable(f), x);");

        # Heuristic zero test:
        gp("IsApprox0(x) = abs(x) < 10^(-default(realprecision) / 2);");

        # Return real positive roots of a complex polynomial
        gp(
        """
        {
        polposroots( f ) = 
            if( poldegree( f ) <= 0, return( [] ) );

            my( Z = polroots( f ) );
            RZ = List();
            for(i = 1,#Z,
                if( IsApprox0( imag( Z[i] ) ) && real( Z[i] ) > 10^(-default(realprecision) / 2), listput( RZ, real(Z[i]) ))
            );
            Vec(RZ);
        }
        """.replace('\n',''));

        # Return distinct roots in [0, 1] of a complex polynomial
        gp(
            """
            {
            pol01roots( f ) =
                if( poldegree( f ) <= 0, return([]));
                my( Z = polroots(f), i, z, eps = 10^(-default(realprecision) / 2));
                RZ = List();
                for(i = 1,#Z,
                    if( IsApprox0( imag(Z[i]) ),
                        z = real( Z[i] );
                        if( z > eps && z < 1-eps && (#RZ==0 || abs( z - RZ[#RZ]) > eps),
                            listput(RZ, z)
                        )
                    )
                );
                Vec(RZ);
            }
            """.replace('\n','')
        );

        # Normalises a polynomial by removing the null leading terms (should fix the bug reported by Jeroen)
        gp(
            """
            {
            normalise_lc_pol(f) =
                my( d = poldegree(f), d2 = -1, var = variable(f), eps=10^(-default(realprecision)/2) );
                for(i=0,d,
                    if(abs(polcoeff(f,i))>eps,d2=i)
                );
                sum(i=0,d2,polcoeff(f,i)*var^i);
            }
            """.replace('\n','')
        );

        # Returns values in [0, 1] where f crosses one of the branch cuts
        gp(
            """
            {
            polcrossaxes01(f) =
                my( u = real(f), v = imag(f), Z, Zf, eps = 10^(-default(realprecision)/2) );
                Z = pol01roots(f);
                u = normalise_lc_pol(u);
                v = normalise_lc_pol(v);
                if( u == 0 || v == 0, return(Z));
                Zf=prod(i=1,#Z,'t-Z[i]);
                Z=concat(Z,pol01roots(u\Zf)); Z=concat(Z,pol01roots(v\Zf));
                Z=vecsort(Z); /* TODO fix can be multiple roots there */
            }
            """.replace('\n','')
        );

        # Integrates division by sqrt(f) from a to b (complexes), assuming no branch
        # cut and f square free. Detects zeros at ends and simplifies them out by changing
        # variable.
        gp("""{
        AJ_fin_nocut(f, a, b, branch, mesh)= 

         /*print("AJfin");*/
         my(t, g, h, S, z0, z1);

         g = EvPol(f, a + (b - a) * 't);

         z0 = IsApprox0( polcoeff(g, 0) );
         z1 = IsApprox0( EvPol(g, 1) );

         if(z0,

          if(z1,
           /*print("we have zeros at both ends");*/
           h = g\('t-'t^2);
           S = 6 * intnum( t = 0, 1, [1, a + (b - a) * t^2 * (3 - 2*t)] / rac2((3 - 2*t)*(2*t + 1) * EvPol(h, t^2 * (3 - 2*t)) ,branch),mesh)
          ,
           /*print("a is a zero, and b is not");*/
           h = g\ 't;
           S = 2 * intnum( t = 0, 1, [1 , a + (b - a) * t^2 ] / rac2( EvPol( h, t^2 ), branch), mesh)
          )
         ,
          if(z1,
           /*print("b is zero, and a is not");*/
           h = g\(1-'t);
           S = 2 * intnum(t = 0, 1, [1, b + (a - b) * t^2 ] / rac2( EvPol( h, 1 - t^2 ), branch), mesh)
          ,
           /*print("neither a or b are zeros");*/
           S = intnum(t = 0, 1, [1, a + (b - a) * t ] / rac2(EvPol(g, t), branch), mesh)
          )
         );
         (b - a) * S;
        }""".replace('\n',''));





        # Abelian integral from x=a straight to x=b, starting at y=ya. Also returns yb.
        # Technique: Keep track of a branch switch throughout.
        # (The correct signed root is chosen by explicit continuation.)
        gp(
        """
        AJ_fin(f, a, b, ya, mesh) = 
        {
         my(g, stops, S = [0, 0], s, yprev, z, i, w, branch, ynew);

         /*print(["Int",a,b]);*/
         g = EvPol(f, a + (b - a) * 't);

         /* Checks for which 0<t<1 the image by f of the integration path crosses the axes */
         stops = polcrossaxes01( g ); 

         /*print(Str(#stops, "obstacles"));*/
         /* Turn these t into z */
         stops = apply( t-> a + t * (b-a), stops);
         /* Add endpoints */
         stops = concat([a], concat(stops,[b])); 
         /* print(stops); */
         
         /* Keep track of the sheet */
         s = 1;
         yprev = ya;
         z = a;

         for(i = 1,#stops - 1,

          w = EvPol(f,z);
          /* w is in the area we're starting from now. */
          if( IsApprox0(w), w = EvPol(deriv(f),z));
          branch = best_branch(w);

          /* Choose a sqrt accordingly */
          /* print(Str("New value ",w)); */
          /* printbranch(branch); */
          if( !IsApprox0( EvPol(f, z) ),
          /* adjust sign so that we stay in the same sheet */
           ynew = rac2( EvPol(f, z), branch);
           if( s * real(ynew/yprev) < 0, s = -s );
          );
          S += s * AJ_fin_nocut(f, stops[i], stops[i+1], branch, mesh);

          /* integrate 1 step */

          /* compute new start and new y */
          z = stops[i + 1];
          yprev = s * rac2(EvPol(f, z) , branch);
         );
         [S, yprev];
        }
        """.replace('\n',''));

        # Image of P+Q-OO by Abel-Jacobi
        gp(
        """
        {
        AJ2(f, P, Q, mesh) =
         my();
         Z = polroots(f);
         /* Integrate P+Q-2(z,0) instead, it is the same up to a period */
         z = Z[1];
         AJ_fin(f, P[1], z, P[2], mesh)[1] + AJ_fin(f, Q[1], z, Q[2], mesh)[1];
        }
        """.replace('\n',''));

        # Image of P-OO by Abel-Jacobi
        gp(
        """
        {
        AJ1(f,P,mesh)= 
            AJ1(f, P, mesh, 1);
        }
        """.replace('\n',''));
        gp(
        """
        {
        AJ1(f,P,mesh,i)= 
         my();
         Z = polroots(f);
         z = Z[i]; /* Integrate P-(z,0) instead, it is the same up to a period */
         AJ_fin(f,P[1],z,P[2],mesh)[1];
        }
        """.replace('\n',''));
        gp(
        """
        {
        whatZ(f)=
         Z = polroots(f);
         z = Z[1];
         z;
        }
        """.replace('\n',''));
load_gp();


# Initializers:
def SetPariPrec(prec):
    gp.set_default('realprecision',prec)

def PariIntNumInit(n=1):
    gp("my_mesh = intnuminit(0,1,{})".format(str(n)));

# Wrappers of functions above:
def AJ2(f,P,Q):
    J = gp("AJ2({},[{},{}],[{},{}],my_mesh)".format(str(f),str(P[0]),str(P[1]),str(Q[0]),str(Q[1])))
    return (J[1],J[2])


def AJ1(R, f,P, i = 1):
    #if f.degree() == 5:
    #    return AJ1odd(R, f, P, i)
    assert f.degree() == 6
    J = gp("AJ1({},[{},{}],my_mesh,{})".format( str(f),str(P[0]),str(P[1]), str(i) ))
    return vector(R, (J[1],J[2]))


def Z(f):
    print gp("whatZ({})".format( str(f) ));


# Dealing with hyperelliptic polynomials of odd degree by inversion:
def AJ1odd(R, f, P, i = 1):
    x0, y0 = P
    if x0 == 0:
        x = (f.parent()).gen()
        f = f(x - 1);
        x0 += 1;
    g = f.parent()(list(reversed(f.padded_list(7))))
    print gp.polroots(g)[1]
    v = R(1/x0)
    w = R(y0/(x0**3))
    total = ( vector(R, AJ1(R, g,(0,0), i)) - vector(R, AJ1(R, g,(v,w), i)) ).list();
    total.reverse();
    return vector(R, total);


def AJ1_digits(f, P, desired_digits, mesh = 4):
    for k in range(1, f.degree()+1):
        assert mesh >= 2;
        # the guess should take less than 2s
        working_digits = 300;
        SetPariPrec(working_digits);
        CF = ComplexField(ceil(RR(log(10)/log(2)) * working_digits));
        aj = {};
        correct_digits = {};
        sum = 0
        total = 0;
        guess = 10000;
        CF = ComplexField(log(10)/log(2)*working_digits)
        for i in range(mesh, max(0, mesh -3) ,-1):
            PariIntNumInit(i);
            aj[i] = AJ1(CF, f, P, k);
            if i < mesh:
                correct_digits[i] = RR((-log(max([ abs(CF((aj[mesh][j] - x)/aj[mesh][j])) for j,x in enumerate(aj[i]) ]))/ log(10.))/working_digits);
            if i + 1 < mesh:
                sum +=  correct_digits[i+1]/correct_digits[i]
                total += 1;
        for i in sorted(correct_digits.keys()):
            if correct_digits[i] > 1:
                guess = i;
                break;
        avg = sum/total
        if guess == 0 and avg > 1.1:
            guess = mesh  - 1 + ceil( log( (1.1/correct_digits[mesh - 1]) )/log( avg) );

        if guess != 0 and 2**guess * desired_digits**2 * log(RR(desired_digits))**3 < ( 300**2 * 2**11 * log(RR(300))**3 ):
            SetPariPrec( ceil(desired_digits) + 10*guess );
            PariIntNumInit(guess);
            CF = ComplexField(log(10)/log(2)*desired_digits);
            result = vector(CF, AJ1(CF, f, P, k));
            # avoiding memory leaks
            gp._reset_expect();
            load_gp();
            return result; 
    else:
        gp._reset_expect();
        load_gp();
        raise OverflowError;


