"""
 *  Some examples of endomorphism rings of genus 2 hyperelliptic curves
 *
 *  Copyright (C) 2016  J.R. Sijsling (sijsling@gmail.com)
 *
 *  Distributed under the terms of the GNU General License (GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc., 51
 *  Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

# Add if no initilization script set:
load('Initialize.sage')

# Ambient ring:
R.<x> = PolynomialRing(QQ)

# Curve input: specify g and h in its equation y^2 + h y = g.
# Hyperelliptic:
f = -15*x^8 + 420*x^6 - 3780*x^4 + 8400*x^2 - 23580*x + 1680
h = 0
f = x^7 - 14*x^6 + 210*x^5 - 658*x^4 + 245*x^3 + 588*x^2 + 637*x - 686
h = 0

End = EndomorphismData(f, h, prec = 200)
AsAlg, As, Rs = End.geometric_representations()
print AsAlg
print [ magma.Discriminant(magma.MinimalPolynomial(A)) for A in AsAlg ]

dim = len(AsAlg)
B = 1
D = [-B..B]
CP = cartesian_product([ D for i in range(dim) ])
deg_min = 10^6
for tup in CP:
    print tup
    R = sum([ tup[i]*Rs[i+1] for i in range(dim) ])
    A = sum([ tup[i]*AsAlg[i+1] for i in range(dim) ])
    p = magma.MinimalPolynomial(R)
    deg = End.degree_estimate(A)
    if magma.Degree(p) == dim:
        if deg < deg_min:
            deg_min = deg
            R_min = R
            A_min = A
print deg_min
print R_min
print A_min
print magma.Parent(A_min[1,1])

