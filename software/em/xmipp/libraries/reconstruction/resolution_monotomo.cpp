/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2016)
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

#include "resolution_monotomo.h"

//#define DEBUG
//#define DEBUG_MASK

void ProgMonoTomoRes::readParams()
{
	fnOdd = getParam("--odd_volume");
	fnEven = getParam("--even_volume");
	fnTomo = getParam("--volume");
	fnOut = getParam("-o");
	fnchim = getParam("--chimera_volume");
	sampling = getDoubleParam("--sampling_rate");
	minRes = getDoubleParam("--minRes");
	maxRes = getDoubleParam("--maxRes");
	fnSym = getParam("--sym");
	N_freq = getDoubleParam("--number_frequencies");
	trimBound = getDoubleParam("--trimmed");
	significance = getDoubleParam("--significance");
	fnMd = getParam("--md_outputdata");
}


void ProgMonoTomoRes::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --odd_volume <vol_file=\"\">   : Input volume");
	addParamsLine("  [--even_volume <vol_file=\"\">]: Half volume 2");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--volume <vol_file=\"\">]: Mean volume of half1 and half2 (only it is neccesary the two haves are used)");
	addParamsLine("  --sym <symmetry>: Symmetry (c1, c2, c3,..d1, d2, d3,...)");
	addParamsLine("  [--chimera_volume <output=\"Chimera_resolution_volume.vol\">]: Local resolution volume for chimera viewer (in Angstroms)");
	addParamsLine("  [--sampling_rate <s=1>]   : Sampling rate (A/px)");
	addParamsLine("                            : Use -1 to disable this option");
	addParamsLine("  [--number_frequencies <w=50>]       : The resolution is computed at a number of frequencies between mininum and");
	addParamsLine("                            : maximum resolution px/A. This parameter determines that number");
	addParamsLine("  [--minRes <s=30>]         : Minimum resolution (A)");
	addParamsLine("  [--maxRes <s=1>]          : Maximum resolution (A)");
	addParamsLine("  [--trimmed <s=0.5>]       : Trimming percentile");
	addParamsLine("  [--significance <s=0.95>]    : The level of confidence for the hypothesis test.");
	addParamsLine("  [--md_outputdata <file=\".\">]  : It is a control file. The provided mask can contain voxels of noise.");
	addParamsLine("                                  : Moreover, voxels inside the mask cannot be measured due to an unsignificant");
	addParamsLine("                                  : SNR. Thus, a new mask is created. This metadata file, shows, the number of");
	addParamsLine("                                  : voxels of the original mask, and the created mask");
}


void ProgMonoTomoRes::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;
	Image<double> V;
	Image<double> V1, V2;
	V1.read(fnOdd);
	V2.read(fnEven);

	if (fnTomo == ""){
		V()=0.5*(V1()+V2());

		V.write("meanVolume.vol");
	}
	else{
	    V.read(fnTomo);
	}
	//V().setXmippOrigin();



	//Defining the volume
	MultidimArray<double> &inputVol = V();
	MultidimArray<double> noiseVolume;



	//Defining the noise
	V1()-=V2();
	noiseVolume = V1()/2;
	Image<double> savenoise = noiseVolume;
	savenoise.write("noise_volume.vol");

	V1.clear();
	V2.clear();

	//Fourier transform slice by slice
	size_t Ydim, Xdim, Zdim, Ndim, Ydim_aux, Xdim_aux, Zdim_aux, Ndim_aux;
	inputVol.getDimensions(Xdim, Ydim, Zdim, Ndim);
	inputVol.printShape();

	MultidimArray<double> slice, sliceNoise, aux_slice, aux_noise;
	MultidimArray< std::complex<double> > fftSlice, fftNoise, fftSlice_aux, fftNoise_aux;

	const MultidimArray<double> &myvolume = inputVol;
	const MultidimArray<double> &mynoise = noiseVolume;
	myvolume.getSlice(0, slice);
	mynoise.getSlice(0, sliceNoise);

	VRiesz.resizeNoCopy(slice);

	FourierTransformer transformer, transformer2;
	transformer.FourierTransform(slice, fftSlice_aux);
	transformer2.FourierTransform(sliceNoise, fftNoise_aux);
	fftSlice_aux.getDimensions(Xdim_aux, Ydim_aux, Zdim_aux, Ndim_aux);
	fftVol.initZeros(1, Zdim, Ydim_aux, Xdim_aux);
	fftNoiseVol = fftVol;
	fftVol.setSlice(0, fftSlice_aux);
	fftNoiseVol.setSlice(0, fftNoise_aux);
	std::cout << "Xdim = " << Xdim << "Ydim = " << Ydim << "Zdim = " << Zdim << std::endl;
	std::cout << "Zdim_aux = " << Zdim_aux << std::endl;
	fftVol.printShape();
	fftNoiseVol.printShape();
	for (size_t j = 1; j< Zdim; j++)
	{
		myvolume.getSlice(j, slice);
		mynoise.getSlice(j, sliceNoise);
		transformer.FourierTransform(slice, fftSlice_aux);
		transformer2.FourierTransform(sliceNoise, fftNoise_aux);
		fftVol.setSlice(j, fftSlice_aux);
		fftNoiseVol.setSlice(j, fftNoise_aux);
	}


	// Calculate u and first component of Riesz vector
	double uy, uy2, ux, uz2, u2;
	long n=0;
	iu.initZeros(Ydim, Xdim);
	for(size_t i=0; i<YSIZE(fftSlice); ++i)
	{
		FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
		uy2=uy*uy;
		for(size_t j=0; j<XSIZE(fftSlice); ++j)
		{
			FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
			u2=uy2+ux*ux;
			if ((i != 0) || (j != 0))
				DIRECT_A2D_ELEM(iu,i,j) = 1.0/sqrt(u2);
			else
				DIRECT_A2D_ELEM(iu,i,j) = 1e38;
			++n;
		}
	}
	#ifdef DEBUG
	Image<double> savedebug;
	savedebug = iu;
	savedebug.write("frequencies.xmp");
	#endif


	// Prepare low pass filter
	lowPassFilter.FilterShape = RAISED_COSINE;
	lowPassFilter.raised_w = 0.01;
	lowPassFilter.do_generate_3dmask = false;
	lowPassFilter.FilterBand = LOWPASS;

	// Prepare mask
	MultidimArray<int> &pMask=mask();

	mask().initZeros(inputVol);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pMask)
	{
		NVoxelsOriginalMask++;
		DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
	}

	V.clear();

}


void ProgMonoTomoRes::amplitudeMonogenicSignal3D(MultidimArray< std::complex<double> > &myfftV,
		double w1, double w1h, double w1l, MultidimArray<double> &amplitudeVol, int count,
		FileName fnDebug)
{

	//fftVRiesz.initZeros(myfftV);
	MultidimArray<double>  amplitude;

	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude

	double iw=1.0/w1;
	double iwl=1.0/w1l;
	double ideltal=PI/(w1-w1l);

	MultidimArray< std::complex<double> > fftSlice;
	size_t Xdim, Ydim, Zdim, Ndim, Zdim_aux;
	myfftV.getDimensions(Xdim, Ydim, Zdim_aux, Ndim);
	VRiesz.getDimensions(Xdim, Ydim, Zdim, Ndim);
	amplitudeVol.resizeNoCopy(Ndim, Zdim_aux, Ydim, Xdim);
	amplitude.initZeros(VRiesz);
	//////////////////////////////
	std::cout << Zdim_aux << std::endl;
	for (size_t ss = 0; ss < Zdim_aux; ss++)
	{
		//std::cout << ss << std::endl;
		myfftV.getSlice(ss, fftSlice);
		long n=0;
		fftVRiesz.initZeros(fftSlice);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftSlice)
		{
			double iun=DIRECT_MULTIDIM_ELEM(iu,n);
			double un=1.0/iun;
			if (w1l<=un && un<=w1)
			{
				//double H=0.5*(1+cos((un-w1)*ideltal));
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(fftSlice, n);
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-w1)*ideltal));//H;
			} else if (un>w1)
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(fftSlice, n);
		}

		transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

		#ifdef DEBUG
		Image<double> filteredvolume;
		filteredvolume = VRiesz;
		filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
		#endif

		#ifdef DEBUG
		FileName iternumber;
		iternumber = formatString("_Volume_%i.vol", count);
		Image<double> saveImg2;
		saveImg2() = VRiesz;
		  if (fnDebug.c_str() != "")
		  {
			saveImg2.write(fnDebug+iternumber);
		  }
		saveImg2.clear();
		#endif


		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

		// Calculate first component of Riesz vector
		fftVRiesz.initZeros(fftSlice);
		double uz, uy, ux;
		n=0;

		for(size_t j=0; j<XSIZE(fftSlice); ++j)
		{
			FFT_IDX2DIGFREQ(j,YSIZE(amplitude),ux);
				for(size_t i=0; i<YSIZE(fftSlice); ++i)
				{
					double iun=DIRECT_MULTIDIM_ELEM(iu,n);
					double un=1.0/iun;
					if (w1l<=un && un<=w1)
					{
						//double H=0.5*(1+cos((un-w1)*ideltal));
						//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-J*ux*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
						//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *=0.5*(1+cos((un-w1)*ideltal));//H;
						//Next lines are an optimization of the commented ones
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(fftSlice, n);
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= -ux*iun*0.5*(1+cos((un-w1)*ideltal));//H;
					} else if (un>w1)
					{
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-ux*iun)*DIRECT_MULTIDIM_ELEM(fftSlice, n);
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
					}
					++n;
				}
		}

		transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

		// Calculate second component of Riesz vector
		fftVRiesz.initZeros(fftSlice);
		n=0;

		for(size_t i=0; i<YSIZE(fftSlice); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(amplitude),uy);
				for(size_t j=0; j<XSIZE(fftSlice); ++j)
				{
					double iun=DIRECT_MULTIDIM_ELEM(iu,n);
					double un=1.0/iun;
					if (w1l<=un && un<=w1)
					{
						//double H=0.5*(1+cos((un-w1)*ideltal));
						//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-J*uy*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
						//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *=0.5*(1+cos((un-w1)*ideltal));//H;
						//Next lines are an optimization of the commented ones
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= -uy*iun*0.5*(1+cos((un-w1)*ideltal));//H;
					} else if (un>w1)
					{
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-uy*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
						DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
					}
					++n;
				}
		}

		transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		{
			DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
			DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
		}

		#ifdef DEBUG
		if (fnDebug.c_str() != "")
		{
		Image<double> saveImg;
		saveImg = amplitude;
		iternumber = formatString("_Amplitude_%i.vol", count);
		saveImg.write(fnDebug+iternumber);
		saveImg.clear();
		}
		#endif // DEBUG
	//
		// Low pass filter the monogenic amplitude
		lowPassFilter.w1 = w1;
		amplitude.setXmippOrigin();
		lowPassFilter.applyMaskSpace(amplitude);

		#ifdef DEBUG
		saveImg2 = amplitude;
		FileName fnSaveImg2;
		if (fnDebug.c_str() != "")
		{
			iternumber = formatString("_Filtered_Amplitude_%i.vol", count);
			saveImg2.write(fnDebug+iternumber);
		}
		saveImg2.clear();
		#endif // DEBUG

		amplitudeVol.setSlice(ss, amplitude);

	}
}

/*
void ProgMonoTomoRes::postProcessingLocalResolutions(MultidimArray<double> &resolutionVol,
		std::vector<double> &list, MultidimArray<double> &resolutionChimera,
		double &cut_value, MultidimArray<int> &pMask)
{
	std::cout << "post processing...." << std::endl;
	MultidimArray<double> resolutionVol_aux = resolutionVol;
	double last_resolution_2 = sampling/list[(list.size()-1)];
	std::cout << "last_res = " << last_resolution_2 << std::endl;

	// Count number of voxels with resolution
	size_t N=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n)>(last_resolution_2-0.001)) //the value 0.001 is a tolerance
			++N;

	std::cout << "post processing...." << std::endl;
	// Get all resolution values
	MultidimArray<double> resolutions(N);
	std::cout << "N = " << N << std::endl;
	size_t N_iter=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n)>(last_resolution_2-0.001))
			DIRECT_MULTIDIM_ELEM(resolutions,N_iter++)=DIRECT_MULTIDIM_ELEM(resolutionVol, n);
	// Sort value and get threshold
	std::cout << "post processing...." << std::endl;
	std::sort(&A1D_ELEM(resolutions,0),&A1D_ELEM(resolutions,N));
	double filling_value = A1D_ELEM(resolutions, (int)(0.5*N)); //median value
	double trimming_value = A1D_ELEM(resolutions, (int)((1-cut_value)*N));

	double freq, res, init_res, last_res;
	std::cout << "post processing...." << std::endl;
	init_res = sampling/list[0];
	last_res = sampling/list[(list.size()-1)];
	
	std::cout << "--------------------------" << std::endl;
	std::cout << "last_res = " << last_res << std::endl;

	resolutionChimera = resolutionVol;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
	{
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n) < last_res)
		{
			DIRECT_MULTIDIM_ELEM(resolutionChimera, n) = filling_value;
			DIRECT_MULTIDIM_ELEM(resolutionVol, n) = 0;
			DIRECT_MULTIDIM_ELEM(pMask,n) = 0;
		}
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n) >= trimming_value)
		{
		  DIRECT_MULTIDIM_ELEM(pMask,n) = 0;
		  DIRECT_MULTIDIM_ELEM(resolutionVol, n) = 0;
		}
	}
	//#ifdef DEBUG_MASK
	Image<int> imgMask;
	imgMask = pMask;
	imgMask.write(fnMaskOut);
	//#endif
}

*/
void ProgMonoTomoRes::run()
{
	produceSideInfo();

	Image<double> outputResolution;
	outputResolution().initZeros(mask());

	MultidimArray<int> &pMask = mask();
	MultidimArray<double> &pOutputResolution = outputResolution();
	MultidimArray<double> &pVfiltered = Vfiltered();
	MultidimArray<double> &pVresolutionFiltered = VresolutionFiltered();
	MultidimArray<double> amplitudeMS, amplitudeMN;

	std::cout << "Looking for maximum frequency ..." << std::endl;
	double criticalZ=icdf_gauss(significance);
	double criticalW=-1;
	double resolution, last_resolution = 10000;  //A huge value for achieving last_resolution < resolution
	double freq, freqH, freqL, resVal, counter;
	double max_meanS = -1e38;
	double cut_value = 0.025;

	double range = maxRes-minRes;

	double R_ = range/N_freq;

	if (R_<0.1)
		R_=0.1;

	double w0 = sampling/maxRes;
	double wF = sampling/minRes;
	double w=w0;
	bool doNextIteration=true;
	bool lefttrimming = false;
	int iter=0;
	int count_res = 0;
	std::vector<double> list;

	std::cout << "Analyzing frequencies" << std::endl;
	std::vector<double> noiseValues;
	FileName fnDebug;

	do
	{
		resolution = maxRes - count_res*R_;
		freqL = sampling/(resolution+R_);
		freq = sampling/resolution;
		freqH = sampling/(resolution-R_);
		if (freq > 0.5)
		{
		  std::cout << "search stopped due to Nyquist limit has been reached" << std::endl;
		  break;
		}
		++count_res;
		if (count_res<=2)
			counter = 0; //maxRes/R_;
		else
			counter = 2;//count_res-2;

		std::cout << "Iteration " << iter << " Freq = " << freq << " Resolution = " << resolution << " (A)" << std::endl;
		//std::cout << "             " << " FreqLOW = " << freqL << " FreqHIGH = " << freqH << std::endl;

		fnDebug = "Signal";

		std::cout << "amplitude signal " << std::endl;
		amplitudeMonogenicSignal3D(fftVol, freq, freqH, freqL, amplitudeMS, iter, fnDebug);
		amplitudeMonogenicSignal3D(fftNoiseVol, freq, freqH, freqL, amplitudeMN, iter, fnDebug);
		std::cout << "amplitude noise finished " << std::endl;


		list.push_back(freq);

		double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
		noiseValues.clear();

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
			double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
			if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
			{
				sumS  += amplitudeValue;
				sumS2 += amplitudeValue*amplitudeValue;
				++NS;
			}
			if (DIRECT_MULTIDIM_ELEM(pMask, n)>=0)
			{
				sumN  += amplitudeValueN;
				sumN2 += amplitudeValueN*amplitudeValueN;
				++NN;
			}
		}
	
		#ifdef DEBUG
		std::cout << "NS" << NS << std::endl;
		std::cout << "NVoxelsOriginalMask" << NVoxelsOriginalMask << std::endl;
		std::cout << "NS/NVoxelsOriginalMask = " << NS/NVoxelsOriginalMask << std::endl;
		#endif
		
		if (NS == 0)
		{
			std::cout << "There are no points to compute inside the mask" << std::endl;
			std::cout << "If the number of computed frequencies is low, perhaps the provided"
					"mask is not enough tight to the volume, in that case please try another mask" << std::endl;
			break;
		}

		double meanS=sumS/NS;
		double sigma2S=sumS2/NS-meanS*meanS;
		double meanN=sumN/NN;
		double sigma2N=sumN2/NN-meanN*meanN;

		if (meanS>max_meanS)
			max_meanS = meanS;

		if (meanS<0.001*max_meanS)
		{
			std::cout << "Search of resolutions stopped due to too low signal" << std::endl;
			break;
		}

		// Check local resolution
		double thresholdNoise;
		thresholdNoise = meanN+criticalZ*sqrt(sigma2N);

		#ifdef DEBUG
		  std::cout << "Iteration = " << iter << ",   Resolution= " << resolution << ",   Signal = " << meanS << ",   Noise = " << meanN << ",  Threshold = " << thresholdNoise <<std::endl;
		#endif

		std::cout << "calculating is is higher than threshold" << std::endl;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
				if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
				{
					DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
					DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution;//sampling/freq;
				}
				else{

					DIRECT_MULTIDIM_ELEM(pMask, n) = DIRECT_MULTIDIM_ELEM(pMask, n) + 1;
					if (DIRECT_MULTIDIM_ELEM(pMask, n) >2)
					{
						DIRECT_MULTIDIM_ELEM(pMask, n) = -1;
						DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution + counter*R_;//maxRes - counter*R_;
					}
				}
		}
		std::cout << "after if it is higher than threshold" << std::endl;
		#ifdef DEBUG_MASK
		FileName fnmask_debug;
		fnmask_debug = formatString("maske_%i.vol", iter);
		mask.write(fnmask_debug);
		#endif

		// Is the mean inside the signal significantly different from the noise?
		double z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);
		#ifdef DEBUG
			std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
			std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
			std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
			std::cout << "  z=" << z << " (" << criticalZ << ")" << std::endl;
		#endif
		if (z<criticalZ)
		{
			criticalW = freq;
			doNextIteration=false;
		}
		if (doNextIteration)
		{
			if (resolution <= (minRes-0.001))
				doNextIteration = false;
		}

		iter++;
	} while (doNextIteration);


	amplitudeMN.clear();
	amplitudeMS.clear();

	double last_resolution_2 = resolution;
	if (fnSym!="c1")
	{
		SymList SL;
		SL.readSymmetryFile(fnSym);
		MultidimArray<double> VSimetrized;
		symmetrizeVolume(SL, pOutputResolution, VSimetrized, LINEAR, DONT_WRAP);
		outputResolution() = VSimetrized;
		VSimetrized.clear();
	}

	#ifdef DEBUG
		outputResolution.write("resolution_simple_simmetrized.vol");
	#endif

	//MultidimArray<double> resolutionFiltered, resolutionChimera;
	//postProcessingLocalResolutions(pOutputResolution, list, resolutionChimera, cut_value, pMask);


	Image<double> outputResolutionImage;
	outputResolutionImage() = pOutputResolution;//resolutionFiltered;
	outputResolutionImage.write(fnOut);
	outputResolutionImage() = pOutputResolution;
	outputResolutionImage.write(fnchim);


	#ifdef DEBUG
		outputResolution.write("resolution_simple.vol");
	#endif


	Nvoxels = NVoxelsOriginalMask;
	
	
	MetaData md;
	size_t objId;
	objId = md.addObject();
	md.setValue(MDL_IMAGE, fnOut, objId);
	md.setValue(MDL_COUNT, (size_t) NVoxelsOriginalMask, objId);
	md.setValue(MDL_COUNT2, (size_t) Nvoxels, objId);

	md.write(fnMd);

}
