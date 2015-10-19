/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
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

#ifndef _PROG_ANGULAR_DISCRIMINATE
#define _PROG_ANGULAR_DISCRIMINATE

#include <data/xmipp_program.h>
#include <data/ctf.h>
#include "fourier_projection.h"
#include "fourier_filter.h"

/**@defgroup AngularDiscriminate Discriminate between several models
   @ingroup ReconsLibrary */
//@{

/** Predict Continuous Parameters. */
class ProgAngularDiscriminate: public XmippMetadataProgram
{
public:
    /** Filename of the reference volumes */
    FileName fnVols;
    /** Maximum gray scale change */
    double maxA;
    /** Maximum gray shift change */
    double maxB;
    // Phase Flipped
    bool phaseFlipped;
    // Sampling rate
    double Ts;
    /** Maximum frequency (A) */
    double maxResol;
    /** Padding factor */
    int pad;
    /** Maximum radius */
    int Rmax;
public:
    // 2D mask in real space
    MultidimArray<int> mask2D;
    // Inverse of the sum of Mask2D
    double iMask2Dsum;
    // Fourier projector
    std::vector<FourierProjector *> projector;
    // Reference volume
    std::vector<Image<double> *> V;
    // Volume size
    size_t Xdim;
    // Input image
	Image<double> I;
	// Theoretical projection
	std::vector<Projection *> P;
	// Filter
    FourierFilter filter;
	// Has CTF
	bool hasCTF;
	// CTF
	CTFDescription ctf;
	// CTF image
	MultidimArray<double> *ctfImage;
public:
    /// Empty constructor
    ProgAngularDiscriminate();

    /// Destructor
    ~ProgAngularDiscriminate();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /** Produce side info.
        An exception is thrown if any of the files is not found*/
    void preProcess();

    /** Predict angles and shift.
        At the input the pose parameters must have an initial guess of the
        parameters. At the output they have the estimated pose.*/
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    /** Update CTF image */
    void updateCTFImage(double defocusU, double defocusV, double angle);

    /** Post process */
    void postProcess();
};
//@}
#endif
