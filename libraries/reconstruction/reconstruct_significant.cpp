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

#include "reconstruct_significant.h"
#include <algorithm>

// Define params
void ProgReconstructSignificant::defineParams()
{
    //usage
    addUsageLine("Generate 3D reconstructions from projections using random orientations");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    addParamsLine("  --initvolumes <md_file>      : Set of initial volumes");
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--iter <N=10>]              : Number of iterations");
    addParamsLine("  [--alpha0 <N=0.05>]          : Initial significance");
    addParamsLine("  [--alphaF <N=0.005>]         : Final significance");
    addParamsLine("  [--keepIntermediateVolumes]  : Keep the volume of each iteration");
    addParamsLine("  [--thr <n=1>]                : Number of threads");
    addParamsLine("  [--angularSampling <a=5>]    : Angular sampling in degrees for generating the projection gallery");
}

// Read arguments ==========================================================
void ProgReconstructSignificant::readParams()
{
    fnIn = getParam("-i");
    fnDir = getParam("--odir");
    fnSym = getParam("--sym");
    fnInit = getParam("--initvolumes");
    alpha0 = getDoubleParam("--alpha0");
    alphaF = getDoubleParam("--alphaF");
    Niter = getIntParam("--iter");
    Nthr = getIntParam("--thr");
    keepIntermediateVolumes = checkParam("--keepIntermediateVolumes");
    angularSampling=getDoubleParam("--angularSampling");
}

// Show ====================================================================
void ProgReconstructSignificant::show()
{
    if (verbose > 0)
    {
        std::cout << "Input metadata              : "  << fnIn        << std::endl;
        std::cout << "Output directory            : "  << fnDir       << std::endl;
        std::cout << "Initial significance        : "  << alpha0      << std::endl;
        std::cout << "Final significance          : "  << alphaF      << std::endl;
        std::cout << "Number of iterations        : "  << Niter       << std::endl;
        std::cout << "Number of threads           : "  << Nthr        << std::endl;
        std::cout << "Keep intermediate volumes   : "  << keepIntermediateVolumes << std::endl;
        std::cout << "Angular sampling            : "  << angularSampling << std::endl;
        if (fnSym != "")
            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
        if (fnInit !="")
            std::cout << "Initial volume              : "  << fnInit << std::endl;
    }
}

// Image alignment ========================================================
//#define DEBUG
void progReconstructSignificantThreadAlign(ThreadArgument &thArg)
{
	ProgReconstructSignificant &prm=*((ProgReconstructSignificant *)thArg.workClass);
	MultidimArray<double> mGalleryProjection, mCurrentImage, mCurrentImageAligned;
	Matrix2D<double> M;
	std::vector< Matrix2D<double> > allM;

	int Nimgs=ZSIZE(prm.cc);
	int Nvols=YSIZE(prm.cc);
	int Ndirs=XSIZE(prm.cc);
	int Nbest=std::ceil(prm.currentAlpha*Ndirs*Nvols);
	MultidimArray<double> imgcc(Nvols*Ndirs);
	MultidimArray<double> cdf;
	double one_alpha=1-prm.currentAlpha;

	for (size_t nVolume=0; nVolume<Nvols; ++nVolume)
	{
		prm.mdReconstructionPartial[nVolume*prm.Nthr+thArg.thread_id].clear();
		prm.mdProjectionMatching[nVolume*prm.Nthr+thArg.thread_id].clear();
	}

	FileName fnImg;
	for (int nImg=0; nImg<Nimgs; ++nImg)
	{
		if ((nImg+1)%prm.Nthr==thArg.thread_id)
		{
#ifdef DEBUG
			std::cout << "Processing: " << prm.mdInp[nImg] << std::endl;
#endif
			mCurrentImage.aliasImageInStack(prm.inputImages(),nImg);
			mCurrentImage.setXmippOrigin();
			fnImg=prm.mdInp[nImg];
			allM.clear();

			double bestCorr=-2, bestRot, bestTilt;
			Matrix2D<double> bestM;
			int bestVolume=-1;

			// Compute all correlations
	    	for (size_t nVolume=0; nVolume<Nvols; ++nVolume)
		    	for (size_t nDir=0; nDir<Ndirs; ++nDir)
				{
					mCurrentImageAligned=mCurrentImage;
					mGalleryProjection.aliasImageInStack(prm.gallery[nVolume](),nDir);
					mGalleryProjection.setXmippOrigin();
					double corr=alignImagesConsideringMirrors(mGalleryProjection,mCurrentImageAligned,M,DONT_WRAP);

//					std::cout << prm.mdGallery[nVolume][nDir].fnImg << " corr=" << corr << std::endl;
//					Image<double> save;
//					save()=mGalleryProjection;
//					save.write("PPPgallery.xmp");
//					save()=mCurrentImage;
//					save.write("PPPcurrentImage.xmp");
//					save()=mCurrentImageAligned;
//					save.write("PPPcurrentImageAligned.xmp");
//					char c; std::cin >> c;

					DIRECT_A3D_ELEM(prm.cc,nImg,nVolume,nDir)=corr;
					DIRECT_A1D_ELEM(imgcc,nVolume*Ndirs+nDir)=corr;
					allM.push_back(M);

					if (corr>bestCorr)
					{
						bestM=M;
						bestCorr=corr;
						bestVolume=(int)nVolume;
						bestRot=prm.mdGallery[nVolume][nDir].rot;
						bestTilt=prm.mdGallery[nVolume][nDir].tilt;
					}
				}

	    	// Keep the best assignment for the projection matching
			MetaData &mdProjectionMatching=prm.mdProjectionMatching[bestVolume*prm.Nthr+thArg.thread_id];
			double scale, shiftX, shiftY, anglePsi;
			bool flip;
			transformationMatrix2Parameters2D(bestM,flip,scale,shiftX,shiftY,anglePsi);
			anglePsi*=-1;

			size_t recId=mdProjectionMatching.addObject();
			mdProjectionMatching.setValue(MDL_IMAGE,fnImg,recId);
			mdProjectionMatching.setValue(MDL_ENABLED,1,recId);
			mdProjectionMatching.setValue(MDL_MAXCC,bestCorr,recId);
			mdProjectionMatching.setValue(MDL_ANGLE_ROT,bestRot,recId);
			mdProjectionMatching.setValue(MDL_ANGLE_TILT,bestTilt,recId);
			mdProjectionMatching.setValue(MDL_ANGLE_PSI,anglePsi,recId);
			mdProjectionMatching.setValue(MDL_SHIFT_X,shiftX,recId);
			mdProjectionMatching.setValue(MDL_SHIFT_Y,shiftY,recId);
			mdProjectionMatching.setValue(MDL_FLIP,flip,recId);

	    	// Get the best images
			imgcc.cumlativeDensityFunction(cdf);
			for (size_t nVolume=0; nVolume<Nvols; ++nVolume)
			{
				MetaData &mdPartial=prm.mdReconstructionPartial[nVolume*prm.Nthr+thArg.thread_id];
				for (size_t nDir=0; nDir<Ndirs; ++nDir)
				{
					double cdfthis=DIRECT_A1D_ELEM(cdf,nVolume*Ndirs+nDir);
					if (cdfthis>=one_alpha)
					{
						double cc=DIRECT_A1D_ELEM(imgcc,nVolume*Ndirs+nDir);
						transformationMatrix2Parameters2D(allM[nVolume*Ndirs+nDir],flip,scale,shiftX,shiftY,anglePsi);

						anglePsi*=-1;
						double weight=DIRECT_A1D_ELEM(cdf,nVolume*Ndirs+nDir);
						double angleRot=prm.mdGallery[nVolume][nDir].rot;
						double angleTilt=prm.mdGallery[nVolume][nDir].tilt;
			#ifdef DEBUG
						std::cout << "   Getting Gallery: " << prm.mdGallery[nVolume][nDir].fnImg
								  << " corr=" << cc << " weight=" << weight << " rot=" << angleRot
								  << " tilt=" << angleTilt << std::endl;
			#endif

						recId=mdPartial.addObject();
						mdPartial.setValue(MDL_IMAGE,fnImg,recId);
						mdPartial.setValue(MDL_ENABLED,1,recId);
						mdPartial.setValue(MDL_MAXCC,cc,recId);
						mdPartial.setValue(MDL_ANGLE_ROT,angleRot,recId);
						mdPartial.setValue(MDL_ANGLE_TILT,angleTilt,recId);
						mdPartial.setValue(MDL_ANGLE_PSI,anglePsi,recId);
						mdPartial.setValue(MDL_SHIFT_X,shiftX,recId);
						mdPartial.setValue(MDL_SHIFT_Y,shiftY,recId);
						mdPartial.setValue(MDL_FLIP,flip,recId);
						mdPartial.setValue(MDL_WEIGHT,weight,recId);
					}
				}
			}
#ifdef DEBUG
			std::cout << "Press any key" << std::endl;
			char c; std::cin >> c;
#endif
		}

		if (thArg.thread_id==0)
			progress_bar(nImg+1);
	}
}
#undef DEBUG

// Main routine ------------------------------------------------------------
void ProgReconstructSignificant::run()
{
    show();
    produceSideinfo();

    ThreadManager thMgr(Nthr,this);
    currentAlpha=alpha0;
    double deltaAlpha;
    if (Niter>1)
    	deltaAlpha=(alphaF-alpha0)/(Niter-1);
    else
    	deltaAlpha=0;

    MetaData mdReconstruction, mdPM;
    for (iter=1; iter<=Niter; iter++)
    {
    	// Generate projections from the different volumes
    	generateProjections();

    	cc.initZeros(mdInp.size(),mdGallery.size(),mdGallery[0].size()); // Nimgs, Nvols, Ndirections

    	// Align the input images to the projections
    	std::cout << "Aligning images ...\n";
    	init_progress_bar(mdIn.size());
    	thMgr.run(progReconstructSignificantThreadAlign);
    	progress_bar(mdIn.size());

    	// Write the corresponding angular metadata
		for (size_t nVolume=0; nVolume<mdInit.size(); ++nVolume)
		{
			mdReconstruction.clear();
			mdPM.clear();
			for (int thr=0; thr<Nthr; thr++)
			{
				mdReconstruction.unionAll(mdReconstructionPartial[nVolume*Nthr+thr]);
				mdPM.unionAll(mdProjectionMatching[nVolume*Nthr+thr]);
			}
			if (mdReconstruction.size()>0)
				mdReconstruction.write(formatString("%s/angles_iter%02d_%02d.xmd",fnDir.c_str(),iter,nVolume));
			else
				REPORT_ERROR(ERR_UNCLASSIFIED,formatString("%s/angles_iter%02d_%02d.xmd is empty. Not written.",fnDir.c_str(),iter,nVolume));
			mdPM.write(formatString("%s/images_iter%02d_%02d.xmd",fnDir.c_str(),iter,nVolume));
	    	deleteFile(formatString("%s/gallery_iter%02d_%02d_sampling.xmd",fnDir.c_str(),iter,nVolume));
	    	deleteFile(formatString("%s/gallery_iter%02d_%02d.doc",fnDir.c_str(),iter,nVolume));
	    	deleteFile(formatString("%s/gallery_iter%02d_%02d.stk",fnDir.c_str(),iter,nVolume));
		}

    	// Reconstruct
    	reconstructCurrent();

    	currentAlpha+=deltaAlpha;
    }
}

void ProgReconstructSignificant::reconstructCurrent()
{
	std::cout << "Reconstructing volumes ..." << std::endl;
	for (size_t nVolume=0; nVolume<mdInit.size(); ++nVolume)
	{
		FileName fnAngles=formatString("%s/angles_iter%02d_%02d.xmd",fnDir.c_str(),iter,nVolume);
		FileName fnVolume=formatString("%s/volume_iter%02d_%02d.vol",fnDir.c_str(),iter,nVolume);
		String args=formatString("-i %s -o %s --sym %s --weight --thr %d -v 0",fnAngles.c_str(),fnVolume.c_str(),fnSym.c_str(),Nthr);
		String cmd=(String)"xmipp_reconstruct_fourier "+args;
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

		args=formatString("-i %s --mask circular %d -v 0",fnVolume.c_str(),-XSIZE(inputImages())/2);
		cmd=(String)"xmipp_transform_mask "+args;
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
	}
}

void ProgReconstructSignificant::generateProjections()
{
	int Nvol=(int)mdInit.size();
	FileName fnVol, fnGallery, fnAngles, fnGalleryMetaData;
	std::vector<GalleryImage> galleryNames;
	mdGallery.clear();
	for (int n=0; n<Nvol; n++)
	{
		fnVol=formatString("%s/volume_iter%02d_%02d.vol",fnDir.c_str(),iter-1,n);
		fnGallery=formatString("%s/gallery_iter%02d_%02d.stk",fnDir.c_str(),iter,n);
		fnAngles=formatString("%s/angles_iter%02d_%02d.xmd",fnDir.c_str(),iter-1,n);
		fnGalleryMetaData=formatString("%s/gallery_iter%02d_%02d.doc",fnDir.c_str(),iter,n);
		String args=formatString("-i %s -o %s --sampling_rate %f --sym %s --compute_neighbors --angular_distance -1 --experimental_images %s -v 0",
				fnVol.c_str(),fnGallery.c_str(),angularSampling,fnSym.c_str(),fnAngles.c_str());
		String cmd=(String)"xmipp_angular_project_library "+args;
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

		mdGallery.push_back(galleryNames);
		MetaData mdAux(fnGalleryMetaData);
		galleryNames.clear();
		FOR_ALL_OBJECTS_IN_METADATA(mdAux)
		{
			GalleryImage I;
			mdAux.getValue(MDL_IMAGE,I.fnImg,__iter.objId);
			mdAux.getValue(MDL_ANGLE_ROT,I.rot,__iter.objId);
			mdAux.getValue(MDL_ANGLE_TILT,I.tilt,__iter.objId);
			mdGallery[n].push_back(I);
		}
		gallery[n].read(fnGallery);
	}
}

void ProgReconstructSignificant::produceSideinfo()
{
	mdIn.read(fnIn);
	mdIn.removeDisabled();

	mdIn.fillConstant(MDL_MAXCC,"0.0");
	mdIn.fillConstant(MDL_ANGLE_ROT,"0.0");
	mdIn.fillConstant(MDL_ANGLE_TILT,"0.0");
	mdIn.fillConstant(MDL_ANGLE_PSI,"0.0");
	mdIn.fillConstant(MDL_SHIFT_X,"0.0");
	mdIn.fillConstant(MDL_SHIFT_Y,"0.0");
	mdIn.fillConstant(MDL_WEIGHT,"1.0");

	// Read all input images in memory
	size_t xdim, ydim, zdim, ndim;
	getImageSize(mdIn,xdim, ydim, zdim, ndim);
	inputImages().resizeNoCopy(mdIn.size(), 1, ydim, xdim);

	Image<double> I;
	size_t n=0;
	MultidimArray<double> mCurrentImage;
	FileName fnImg;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE,fnImg,__iter.objId);
		mdInp.push_back(fnImg);
		I.read(fnImg);
		mCurrentImage.aliasImageInStack(inputImages(),n++);
		memcpy(MULTIDIM_ARRAY(mCurrentImage),MULTIDIM_ARRAY(I()),MULTIDIM_SIZE(mCurrentImage)*sizeof(double));
	}

	// Copy all input values as iteration 0 volumes
	mdInit.read(fnInit);
	FileName fnVol, fnAngles;
	Image<double> V, galleryDummy;
	int idx=0;
	FOR_ALL_OBJECTS_IN_METADATA(mdInit)
	{
		mdInit.getValue(MDL_IMAGE,fnVol,__iter.objId);
		V.read(fnVol);
		fnVol=formatString("%s/volume_iter00_%02d.vol",fnDir.c_str(),idx);
		V.write(fnVol);
		mdInit.setValue(MDL_IMAGE,fnVol,__iter.objId);
		fnAngles=formatString("%s/angles_iter00_%02d.xmd",fnDir.c_str(),idx);
		mdIn.write(fnAngles);
		gallery.push_back(galleryDummy);

		// Create a partial metadata for each thread and each volume
		for (int thr=0; thr<Nthr; ++thr)
		{
			mdReconstructionPartial.push_back(MetaData());
			mdProjectionMatching.push_back(MetaData());
		}
		idx++;
	}

	iter=0;
}
