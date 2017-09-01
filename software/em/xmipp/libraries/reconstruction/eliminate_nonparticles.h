/***************************************************************************
 * Authors:     Tomas Majtner (tmajtner@cnb.csic.es)
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

#ifndef _PROG_ELIMINATE_NONPARTICLES
#define _PROG_ELIMINATE_NONPARTICLES

#include <data/metadata.h>
#include <data/xmipp_program.h>


class ProgEliminateNonParticles: public XmippProgram
{
public:
	/// Name of the input metadata
    FileName fnIn;

    /// Name of the output metadata
    FileName fnOut;

    /// Threshold for variance of variance
    float threshold;

    /// SelFile containing input data
    MetaData SF;

    /// Input image data
    Image<double> Iref;

public:
    /// Read input parameters
    void readParams();

    /// Define input parameters
    void defineParams();

    /// Show
    void show();

    /// Execute
    void run();

    /// Function for recognizing non-noise particles
    bool isParticle();
};

#endif