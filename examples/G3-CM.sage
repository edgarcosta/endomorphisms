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
load('../Initialize.sage')

# Ambient ring:
F = QQ
R.<x> = PolynomialRing(F)

# Curve input: specify g and h in its equation y^2 + h y = g.
# Hyperelliptic:
f = -15*x^8 + 420*x^6 - 3780*x^4 + 8400*x^2 - 23580*x + 1680
h = 0
f = x^7 - 14*x^6 + 210*x^5 - 658*x^4 + 245*x^3 + 588*x^2 + 637*x - 686
h = 0
X = mHyperellipticCurve(f, h)

print X
# The main functionality
Endo = EndomorphismData(X, prec = 300, have_oldenburg = True)

print "Field of definition:"
print Endo.endomorphism_field()

print "Geometric representation:"
print Endo.geometric().representation()

print "Lattice:"
print Endo.lattice().pretty_print()
