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

#include "angular_discriminate.h"
#include <data/mask.h>
#include <data/numerical_tools.h>

// Empty constructor =======================================================
ProgAngularDiscriminate::ProgAngularDiscriminate()
{
    produces_a_metadata = true;
}

ProgAngularDiscriminate::~ProgAngularDiscriminate()
{
	for (size_t i=0; i<V.size(); ++i)
	{
		delete projector[i];
		delete V[i];
		delete P[i];
	}
}

// Read arguments ==========================================================
void ProgAngularDiscriminate::readParams()
{
	XmippMetadataProgram::readParams();
    fnVols = getParam("--ref");
    maxA = getDoubleParam("--max_gray_scale");
    maxB = getDoubleParam("--max_gray_shift");
    Ts = getDoubleParam("--sampling");
    phaseFlipped = checkParam("--phaseFlipped");
    maxResol = getDoubleParam("--max_resolution");
    pad = getIntParam("--padding");
    Rmax = getIntParam("--Rmax");
}

// Show ====================================================================
void ProgAngularDiscriminate::show()
{
    if (!verbose)
        return;
	XmippMetadataProgram::show();
    std::cout
    << "Reference volumes:   " << fnVols             << std::endl
	<< "Max. Gray Scale:     " << maxA               << std::endl
	<< "Max. Gray Shift:     " << maxB               << std::endl
    << "Sampling:            " << Ts                 << std::endl
    << "Phase flipped:       " << phaseFlipped       << std::endl
    << "Max. Resolution:     " << maxResol           << std::endl
    << "Padding factor:      " << pad                << std::endl
    << "Max. Radius:         " << Rmax               << std::endl
    ;
}

// usage ===================================================================
void ProgAngularDiscriminate::defineParams()
{
    addUsageLine("Discriminate an image amongst several volumes");
	defaultComments["-i"].clear();
	defaultComments["-i"].addComment("Metadata with initial alignment");
	defaultComments["-o"].clear();
	defaultComments["-o"].addComment("Metadata with output alignment");
    XmippMetadataProgram::defineParams();
    addParamsLine("   --ref <volumes>             : Metadata with several reference volumes");
    addParamsLine("  [--max_gray_scale <a=0.05>]  : Maximum gray scale change");
    addParamsLine("  [--max_gray_shift <b=0.05>]  : Maximum gray shift change as a factor of the image standard deviation");
    addParamsLine("  [--sampling <Ts=1>]          : Sampling rate (A/pixel)");
    addParamsLine("  [--phaseFlipped]             : Input images have been phase flipped");
    addParamsLine("  [--max_resolution <f=8>]     : Maximum resolution (A)");
    addParamsLine("  [--Rmax <R=-1>]              : Maximum radius (px). -1=Half of volume size");
    addParamsLine("  [--padding <p=2>]            : Padding factor");
}

// Produce side information ================================================
void ProgAngularDiscriminate::preProcess()
{
    // Read the reference volume
	MetaData mdRefs;
	mdRefs.read(fnVols);
	if (mdRefs.size()==0)
		REPORT_ERROR(ERR_MD_NOOBJ,"The input metadata of reference volumes is empty");
	FileName fnVol;
	FOR_ALL_OBJECTS_IN_METADATA(mdRefs)
	{
		mdRefs.getValue(MDL_IMAGE,fnVol,__iter.objId);
		Image<double> *ptrV=new Image<double>();
		ptrV->read(fnVol);
		(*ptrV)().setXmippOrigin();
		V.push_back(ptrV);

	    projector.push_back(new FourierProjector((*ptrV)(),pad,Ts/maxResol,BSPLINE3));
	}
    Xdim=XSIZE((*V[0])());

    // Construct mask
    if (Rmax<0)
    	Rmax=Xdim/2;
    Mask mask;
    mask.type = BINARY_CIRCULAR_MASK;
    mask.mode = INNER_MASK;
    mask.R1 = Rmax;
    mask.generate_mask(Xdim,Xdim);
    mask2D=mask.get_binary_mask();
    iMask2Dsum=1.0/mask2D.sum();

    // Low pass filter
    filter.FilterBand=LOWPASS;
    filter.w1=Ts/maxResol;
    filter.raised_w=0.02;
}

void ProgAngularDiscriminate::updateCTFImage(double defocusU, double defocusV, double angle)
{
	ctf.K=1; // get pure CTF with no envelope
	ctf.DeltafU=defocusU;
	ctf.DeltafV=defocusV;
	ctf.azimuthal_angle=angle;
	ctf.produceSideInfo();
	if (ctfImage==NULL)
	{
		ctfImage = new MultidimArray<double>();
		ctfImage->resizeNoCopy(projector[0]->projection());
		STARTINGY(*ctfImage)=STARTINGX(*ctfImage)=0;
	}
	ctf.generateCTF(YSIZE(projector[0]->projection()),XSIZE(projector[0]->projection()),*ctfImage,Ts);
	if (phaseFlipped)
		FOR_ALL_ELEMENTS_IN_ARRAY2D(*ctfImage)
			A2D_ELEM(*ctfImage,i,j)=fabs(A2D_ELEM(*ctfImage,i,j));
}

// Predict =================================================================
//#define DEBUG
void ProgAngularDiscriminate::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
#ifdef NEVERDEFINED
    rowOut=rowIn;

    // Read input image and initial parameters
//  ApplyGeoParams geoParams;
//	geoParams.only_apply_shifts=false;
//	geoParams.wrap=DONT_WRAP;

	rowIn.getValue(MDL_ANGLE_ROT,old_rot);
	rowIn.getValue(MDL_ANGLE_TILT,old_tilt);
	rowIn.getValue(MDL_ANGLE_PSI,old_psi);
	rowIn.getValue(MDL_SHIFT_X,old_shiftX);
	rowIn.getValue(MDL_SHIFT_Y,old_shiftY);
	rowIn.getValue(MDL_FLIP,old_flip);
	double old_scaleX=0, old_scaleY=0;
	old_grayA=1;
	old_grayB=0;
	if (rowIn.containsLabel(MDL_CONTINUOUS_GRAY_A))
	{
		rowIn.getValue(MDL_CONTINUOUS_GRAY_A,old_grayA);
		rowIn.getValue(MDL_CONTINUOUS_GRAY_B,old_grayB);
		rowIn.getValue(MDL_CONTINUOUS_SCALE_X,old_scaleX);
		rowIn.getValue(MDL_CONTINUOUS_SCALE_Y,old_scaleY);
		rowIn.getValue(MDL_CONTINUOUS_X,old_shiftX);
		rowIn.getValue(MDL_CONTINUOUS_Y,old_shiftY);
		rowIn.getValue(MDL_CONTINUOUS_FLIP,old_flip);
	}

	if (rowIn.containsLabel(MDL_CTF_DEFOCUSU) || rowIn.containsLabel(MDL_CTF_MODEL))
	{
		hasCTF=true;
		ctf.readFromMdRow(rowIn);
		ctf.produceSideInfo();
		old_defocusU=ctf.DeltafU;
		old_defocusV=ctf.DeltafV;
		old_defocusAngle=ctf.azimuthal_angle;
	}
	else
		hasCTF=false;

	if (verbose>=2)
		std::cout << "Processing " << fnImg << std::endl;
	I.read(fnImg);
	I().setXmippOrigin();
	Istddev=I().computeStddev();

    Ifiltered()=I();
    filter.applyMaskSpace(Ifiltered());

    Matrix1D<double> p(12), steps(12);
    p(0)=old_grayA; // a in I'=a*I+b
    p(1)=old_grayB; // b in I'=a*I+b
    p(4)=old_scaleX;
    p(5)=old_scaleY;

    // Optimize
	double cost=-1;
	if (fabs(old_scaleX)>maxScale || fabs(old_scaleY)>maxScale)
    	rowOut.setValue(MDL_ENABLED,-1);
	else
	{
		try
		{
			cost=1e38;
			int iter;
			steps.initZeros();
			if (optimizeGrayValues)
				steps(0)=steps(1)=1.;
			if (optimizeShift)
				steps(2)=steps(3)=1.;
			if (optimizeScale)
				steps(4)=steps(5)=1.;
			if (optimizeAngles)
				steps(6)=steps(7)=steps(8)=1.;
			if (optimizeDefocus)
				steps(9)=steps(10)=steps(11)=1.;
			powellOptimizer(p, 1, 12, &continuous2cost, this, 0.01, cost, iter, steps, verbose>=2);
			if (cost>1e30 || (cost>0 && contCost==CONTCOST_CORR))
			{
				rowOut.setValue(MDL_ENABLED,-1);
				p.initZeros();
			    p(0)=old_grayA; // a in I'=a*I+b
			    p(1)=old_grayB; // b in I'=a*I+b
			    p(4)=old_scaleX;
			    p(5)=old_scaleY;
			}
			else
			{
				if (fnResiduals!="")
				{
					FileName fnResidual;
					fnResidual.compose(fnImgOut.getPrefixNumber(),fnResiduals);
					E.write(fnResidual);
					rowOut.setValue(MDL_IMAGE_RESIDUAL,fnResidual);
				}
			}
			if (contCost==CONTCOST_CORR)
				cost=-cost;
			if (verbose>=2)
				std::cout << "I'=" << p(0) << "*I" << "+" << p(1) << " Dshift=(" << p(2) << "," << p(3) << ") "
				          << "scale=(" << 1+p(4) << "," << 1+p(5) << ") Drot=" << p(6) << " Dtilt=" << p(7)
				          << " Dpsi=" << p(8) << " DU=" << p(9) << " DV=" << p(10) << " Dalpha=" << p(11) << std::endl;
			// Apply
			FileName fnOrig;
			rowIn.getValue(MDL::str2Label(originalImageLabel),fnOrig);
			I.read(fnImg);
			if (XSIZE(Ip())!=XSIZE(I()))
			{
				scaleToSize(BSPLINE3,Ip(),I(),XSIZE(Ip()),YSIZE(Ip()));
				I()=Ip();
			}
			A(0,2)=p(2)+old_shiftX;
			A(1,2)=p(3)+old_shiftY;
			A(0,0)=1+p(4);
			A(1,1)=1+p(5);

			if (old_flip)
			{
				MAT_ELEM(A,0,0)*=-1;
				MAT_ELEM(A,0,1)*=-1;
				MAT_ELEM(A,0,2)*=-1;
			}
			applyGeometry(BSPLINE3,Ip(),I(),A,IS_NOT_INV,DONT_WRAP);
			MultidimArray<double> &mIp=Ip();
			double ia=1.0/p(0);
			double b=p(1);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mIp)
			{
				if (DIRECT_MULTIDIM_ELEM(mask2D,n))
					DIRECT_MULTIDIM_ELEM(mIp,n)=ia*(DIRECT_MULTIDIM_ELEM(mIp,n)-b);
				else
					DIRECT_MULTIDIM_ELEM(mIp,n)=0.0;
			}
			Ip.write(fnImgOut);
		}
		catch (XmippError XE)
		{
			std::cerr << XE << std::endl;
			std::cerr << "Warning: Cannot refine " << fnImg << std::endl;
			rowOut.setValue(MDL_ENABLED,-1);
		}
	}
    rowOut.setValue(MDL_IMAGE_ORIGINAL, fnImg);
    rowOut.setValue(MDL_IMAGE, fnImgOut);
    rowOut.setValue(MDL_ANGLE_ROT,  old_rot+p(6));
    rowOut.setValue(MDL_ANGLE_TILT, old_tilt+p(7));
    rowOut.setValue(MDL_ANGLE_PSI,  old_psi+p(8));
    rowOut.setValue(MDL_SHIFT_X,    0.);
    rowOut.setValue(MDL_SHIFT_Y,    0.);
    rowOut.setValue(MDL_FLIP,       false);
    rowOut.setValue(MDL_COST,       cost);
    rowOut.setValue(MDL_CONTINUOUS_GRAY_A,p(0));
    rowOut.setValue(MDL_CONTINUOUS_GRAY_B,p(1));
    rowOut.setValue(MDL_CONTINUOUS_SCALE_X,p(4));
    rowOut.setValue(MDL_CONTINUOUS_SCALE_Y,p(5));
    rowOut.setValue(MDL_CONTINUOUS_X,p(2)+old_shiftX);
    rowOut.setValue(MDL_CONTINUOUS_Y,p(3)+old_shiftY);
    rowOut.setValue(MDL_CONTINUOUS_FLIP,old_flip);
    if (hasCTF)
    {
    	rowOut.setValue(MDL_CTF_DEFOCUSU,old_defocusU+p(9));
    	rowOut.setValue(MDL_CTF_DEFOCUSV,old_defocusV+p(10));
    	rowOut.setValue(MDL_CTF_DEFOCUS_ANGLE,old_defocusAngle+p(11));
    	if (old_defocusU+p(9)<0 || old_defocusU+p(10)<0)
    		rowOut.setValue(MDL_ENABLED,-1);
    }
#endif

#ifdef DEBUG
    std::cout << "p=" << p << std::endl;
    MetaData MDaux;
    MDaux.addRow(rowOut);
    MDaux.write("PPPmd.xmd");
    Image<double> save;
    save()=P();
    save.write("PPPprojection.xmp");
    save()=I();
    save.write("PPPexperimental.xmp");
    //save()=C;
    //save.write("PPPC.xmp");
    Ip.write("PPPexperimentalp.xmp");
    Ifiltered.write("PPPexperimentalFiltered.xmp");
    Ifilteredp.write("PPPexperimentalFilteredp.xmp");
    E.write("PPPresidual.xmp");
    std::cout << A << std::endl;
    std::cout << fnImgOut << " rewritten\n";
    std::cout << "Press any key" << std::endl;
    char c; std::cin >> c;
#endif
}
#undef DEBUG

void ProgAngularDiscriminate::postProcess()
{
#ifdef NEVERDEFINED
	MetaData &ptrMdOut=*getOutputMd();
	ptrMdOut.removeDisabled();
	if (contCost==CONTCOST_L1)
	{
		double minCost=1e38;
		FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
		{
			double cost;
			ptrMdOut.getValue(MDL_COST,cost,__iter.objId);
			if (cost<minCost)
				minCost=cost;
		}
		FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
		{
			double cost;
			ptrMdOut.getValue(MDL_COST,cost,__iter.objId);
			ptrMdOut.setValue(MDL_WEIGHT_CONTINUOUS2,minCost/cost,__iter.objId);
		}
	}
	else
	{
		double maxCost=-1e38;
		FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
		{
			double cost;
			ptrMdOut.getValue(MDL_COST,cost,__iter.objId);
			if (cost>maxCost)
				maxCost=cost;
		}
		FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
		{
			double cost;
			ptrMdOut.getValue(MDL_COST,cost,__iter.objId);
			ptrMdOut.setValue(MDL_WEIGHT_CONTINUOUS2,cost/maxCost,__iter.objId);
		}
	}

	ptrMdOut.write(fn_out.replaceExtension("xmd"));
#endif
}
