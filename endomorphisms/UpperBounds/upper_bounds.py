#  Copyright (C) 2018, Edgar Costa
#  See LICENSE file for license details.

# exposes some of the functionality mentioned in Section 7.3 and Section 7.4

from endomorphisms.OverFiniteField import endomorphism_frob
from .utils import field_intersection_matrix, RR_representation, polredabs
from sage.all import Set, NumberField


def endomorphisms_upper_bound(frob_list, eta_char0 = None):
    r"""

    INPUT:

    - ``frob_list`` -- a list of Frobenius polynomials of A
    - ``eta_char0`` -- a putative value for eta(A)

    OUTPUT:

    - a boolean = did we manage to provide a putative upper bound on the centers?

    - a message, explaining what we have achieved to do.

    - None or a putative value for eta(A) (if this matches the input value, then we proved that this value is correct)

    - None or the number of factors assuming that eta(A) computed above is correct

    - None or a list of tuples
            [ (e_j n_j, n_j * dim A_j, L_j, RR_j ) for j in range(1, t + 1)],
      where we are assuming that the eta(A) and t above are correct and therefore we have
            A^{al} =  A^{al} = (A_1)^n_1 x ... x (A_k)^n_t
      where e_j ^2 is the dimension of End(A_j ^{al}) over its center and
      its center is a subfield of L_j.
      Assuming that L_j is indeed the center, then RR_j = End( A_j ^{n_j} ) \otimes RR.

    - Assuming eta(A) and t are correct, an upper bound for dim_Q End(A)

    EXAMPLE:

        sage: from endomorphisms import endomorphisms_upper_bound
        sage: ZZT.<T> = ZZ[]
        sage: F3 = 1 - T^2 + 9*T^4;
        sage: F7 = 1 + 4*T^2 + 49*T^4;
        sage: F13 = 1 - 8*T^2 + 169*T^4;
        sage: endomorphisms_upper_bound([[3, F3], [7, F7], [13, F13]])
        (True,
         'We have putatively computed eta and t. Under this assumption, we bounded the corresponding centers.',
         4,
         1,
         [(2, 2, [T, [T]], ['M_2(RR)'])],
         4)

    """
    g = int(frob_list[0][1].degree()/2);
    if eta_char0 is None:
        eta = 4*g*g; # max value found for eta_p
    else:
        eta = 2*eta_char0;

    t = g;
    eta_lower = [];
    for p, f in frob_list:
        dimtotal, fieldext, endo = endomorphism_frob(f);
        if dimtotal < eta:
            eta = dimtotal
            t = len(endo);
            eta_lower = [];

        if dimtotal == eta:
            if len(endo) < t:
                t = len(endo)
                eta_lower = []
            if len(endo) == t:
                eta_lower.append([p, f, len(endo), dimtotal, fieldext, endo])

    if len(eta_lower) == 0:
        return False, "We  did not manage to find any prime where eta(A_p) = 2 * eta(A)", None, None, None, None

    eta_char0 = eta/2;

    multiset_char0 = None
    frob_factors = [[None]*len(eta_lower) for _ in range(t)];
    for i, (_, _, tp, etap, _, endo) in enumerate(eta_lower):
        assert t == tp
        assert eta == etap
        # endo[j] = mpj, mpj*deg(hpj), hpj
        # the multiset_char0 in the paper has y divided by 2
        multiset = sorted([(x, y) for x, y, _ in endo])
        if multiset_char0 is None:
            multiset_char0 = multiset
            frob_factors = {}
            for pair in Set(multiset):
                frob_factors[pair] = [ []  for _ in eta_lower ];
        if multiset_char0 != multiset:
            # we only managed to bound eta
            message = "We only managed to find an upper bound for eta.";
            message += " If the upper bound for eta indeed is eta, then the number of factors is a strict upper bound";
            return False, message, eta_char0, t, None, None
        for x, y, hpj in endo:
            # endo[j] = mpj, mpj*deg(hpj), hpj
            frob_factors[(x,y)][i].append( polredabs(hpj) );


    # it looks like we have a consistent upper bound for eta and t
    message = "We have putatively computed eta and t.";
    message += " Under this assumption, we bounded the corresponding centers."

    #We can try to bound the center of each factor
    output = [];
    total_dim = 0;
    for pair, frob_matrix in frob_factors.iteritems():
        L = field_intersection_matrix( frob_matrix );
        ejnj, njdimAj = pair
        njdimAj = njdimAj//2;
        for Lj in L:
            # the only real functionality of the NumberField that we use is
            # Ljmax.is_CM()
            # so it doesn't matter which field of maximal degree we take
            Ljmax = NumberField(Lj[-1][-1], 'a');
            RRj = RR_representation(njdimAj, Ljmax, ejnj)
            output.append( (ejnj, njdimAj, Lj, RRj) )
            total_dim += ejnj**2 * Ljmax.degree()

    return True, message,  eta_char0, t, output, total_dim





def RR_upper_bound(frob_list):
    endo_upper_bound = endomorphisms_upper_bound(frob_list);

    if endo_upper_bound[0] == False:
        return None
    else:
        # endo_upper_bound[4][j] = (ejnj, njdimAj, Lj, RRj)
        return [ elt[3] for elt in endo_upper_bound[4] ]




