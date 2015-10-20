/*
  Implementation of Lambert W function

  Copyright (C) 2009 Darko Veberic, darko.veberic@ung.si

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _utl_LambertW_h_
#define _utl_LambertW_h_

/**
  \author Darko Veberic
  \version $Id: LambertW.h 14958 2009-10-19 11:11:19Z darko $
  \date 25 Jun 2009
*/


  /** Approximate Lambert W function
    Accuracy at least 5 decimal places in all definition range.
    See LambertW() for details.

    \param branch: valid values are 0 and -1
    \param x: real-valued argument \f$\geq-1/e\f$
    \ingroup math
  */

  template<int branch>
  double LambertWApproximation(const double x);


  /** Lambert W function
    \image html utl_LambertW.png

    Lambert function \f$y={\rm W}(x)\f$ is defined as a solution
    to the \f$x=ye^y\f$ expression and is also known as
    "product logarithm". Since the inverse of \f$ye^y\f$ is not
    single-valued, the Lambert function has two real branches
    \f${\rm W}_0\f$ and \f${\rm W}_{-1}\f$.

    \f${\rm W}_0(x)\f$ has real values in the interval
    \f$[-1/e,\infty]\f$ and \f${\rm W}_{-1}(x)\f$ has real values
    in the interval \f$[-1/e,0]\f$.
    Accuracy is the nominal double type resolution
    (16 decimal places).

    \param branch: valid values are 0 and -1
    \param x: real-valued argument \f$\geq-1/e\f$ (range depends on branch)
    \ingroup math
  */

  template<int branch>
  double LambertW(const double x);

#endif
