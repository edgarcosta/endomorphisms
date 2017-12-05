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
load('../../Initialize.sage')

# Ambient ring:
F = QQ
P2.<x,y,z> = ProjectiveSpace(F, 2)

# Curve:
f = x^4 + 8*x^3*z + 2*x^2*y*z + 25*x^2*z^2 - x*y^3 + 2*x*y^2*z + 8*x*y*z^2 \
    + 36*x*z^3 + y^4 - 2*y^3*z + 5*y^2*z^2 + 9*y*z^3 + 20*z^4
X = mPlaneCurve(f)

print X
Endo = EndomorphismData(X, prec = 100, have_oldenburg = True)

print "Geometric representation:"
overK = Endo.over_field("geometric")
endodict = overK.full()
rep = overK.representation()
print rep

#print "Verifying endomorphisms:"
#test = Endo.verify()
#print test
