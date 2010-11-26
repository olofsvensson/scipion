/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#ifndef _PROG_SYMMETRIZE_HH
#  define _PROG_SYMMETRIZE_HH

#include <data/funcs.h>
#include <data/image.h>
#include <data/mask.h>
//#include <data/matrix3d.h >

#include <data/symmetries.h>

/**@defgroup SymmetrizeProgram symmetrize (Symmetrize a volume)
   @ingroup ReconsLibrary */
//@{
/* Test parameters --------------------------------------------------------- */
/// Symmetrize Parameters
class ProgSymmetrize : public XmippProgram
{
public:
    /// input file
    FileName        fn_in;
    /// output file
    FileName        fn_out;
    /// symmetry file
    FileName        fn_sym;
    /// Do not generate subgroup
    bool            do_not_generate_subgroup;
    /// wrap or don't wrap input file during symmetrisation
    bool            wrap;
    /// use Bsplines for interpolating
    bool            useBsplines;
public:
    /** Read parameters from command line. */
    void readParams();

    /** Define Parameters */
    void defineParams();

    /** Run */
    void run();

    /** Show parameters */
    //friend std::ostream & operator << (std::ostream &out, const Symmetrize_Parameters &prm);
};

/** Really symmetrize.*/
void symmetrize(const SymList &SL, MultidimArray<double> &V_in, MultidimArray<double> &V_out,
                int Splinedegree=LINEAR, bool wrap=true, bool show_progress=false,
                bool do_outside_avg=false);

//@}

#endif
