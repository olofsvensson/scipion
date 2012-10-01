/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#include "xmippmodule.h"

/***************************************************************/
/*                    MDLabels constants                       */
/**************************************************************/

void addIntConstant(PyObject * dict, const char * name, const long &value)
{
    PyObject * pyValue = PyInt_FromLong(value);
    PyDict_SetItemString(dict, name, pyValue);
    Py_DECREF(pyValue);
}

//Macro to add constants to xmipp module with the same name
#define ADD_CONST(label) addIntConstant(dict, #label, (long)label)
// and with different name
#define ADD_CONST2(name, label) addIntConstant(dict, name, (long)label)

void addLabels(PyObject * dict)
{
  //Add constants
  ADD_CONST(AGGR_COUNT);
  ADD_CONST(AGGR_MAX);
  ADD_CONST(AGGR_SUM);
  ADD_CONST(AGGR_AVG);

  ADD_CONST(UNION);
  ADD_CONST(UNION_DISTINCT);
  ADD_CONST(INTERSECTION);
  ADD_CONST(SUBSTRACTION);
  ADD_CONST(INNER_JOIN);
  ADD_CONST(LEFT_JOIN);
  ADD_CONST(NATURAL_JOIN);
  ADD_CONST(OUTER_JOIN);
  ADD_CONST(INNER);
  ADD_CONST(LEFT);
  ADD_CONST(OUTER);
  ADD_CONST(NATURAL);
  ADD_CONST(EQ);
  ADD_CONST(NE);
  ADD_CONST(GT);
  ADD_CONST(LT);
  ADD_CONST(GE);
  ADD_CONST(LE);
  ADD_CONST(MD_OVERWRITE);
  ADD_CONST(MD_APPEND);
  ADD_CONST(MDL_UNDEFINED);
  ADD_CONST(MDL_FIRST_LABEL);
  //Metadata Labels
  ADD_CONST(MDL_OBJID);
  ADD_CONST(MDL_ANGLE_PSI2);
  ADD_CONST(MDL_ANGLE_PSI);
  ADD_CONST(MDL_ANGLE_PSI_DIFF);
  ADD_CONST(MDL_ANGLE_ROT2);
  ADD_CONST(MDL_ANGLE_ROT);
  ADD_CONST(MDL_ANGLE_ROT_DIFF);
  ADD_CONST(MDL_ANGLE_TILT2);
  ADD_CONST(MDL_ANGLE_TILT);
  ADD_CONST(MDL_ANGLE_TILT_DIFF);
  ADD_CONST(MDL_ANGLE_DIFF);
  ADD_CONST(MDL_ANGLE_Y);
  ADD_CONST(MDL_ANGLE_Y2);
  ADD_CONST(MDL_IMAGE1);
  ADD_CONST(MDL_IMAGE2);
  ADD_CONST(MDL_IMAGE3);
  ADD_CONST(MDL_IMAGE4);
  ADD_CONST(MDL_IMAGE5);
  ADD_CONST(MDL_AVG);
  ADD_CONST(MDL_BGMEAN);
  ADD_CONST(MDL_BLOCK_NUMBER);
  ADD_CONST(MDL_CLASSIFICATION_DATA);
  ADD_CONST(MDL_CRYSTAL_CELLX);
  ADD_CONST(MDL_CRYSTAL_CELLY);
  ADD_CONST(MDL_COMMENT);
  ADD_CONST(MDL_COST);
  ADD_CONST(MDL_COUNT);
  ADD_CONST(MDL_CTF_INPUTPARAMS);
  ADD_CONST(MDL_CTF_MODEL);
  ADD_CONST(MDL_CTF_MODEL2);
  ADD_CONST(MDL_CTF_SAMPLING_RATE);
  ADD_CONST(MDL_CTF_VOLTAGE);
  ADD_CONST(MDL_CTF_DEFOCUSU);
  ADD_CONST(MDL_CTF_DEFOCUSV);
  ADD_CONST(MDL_CTF_DEFOCUS_ANGLE);
  ADD_CONST(MDL_CTF_CS);
  ADD_CONST(MDL_CTF_CA);
  ADD_CONST(MDL_CTF_GROUP);
  ADD_CONST(MDL_CTF_ENERGY_LOSS);
  ADD_CONST(MDL_CTF_LENS_STABILITY);
  ADD_CONST(MDL_CTF_CONVERGENCE_CONE);
  ADD_CONST(MDL_CTF_LONGITUDINAL_DISPLACEMENT);
  ADD_CONST(MDL_CTF_TRANSVERSAL_DISPLACEMENT);
  ADD_CONST(MDL_CTF_Q0);
  ADD_CONST(MDL_CTF_K);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN_K);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN_SIGMAU);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN_SIGMAV);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN_CU);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN_CV);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN_ANGLE);
  ADD_CONST(MDL_CTF_BG_SQRT_K);
  ADD_CONST(MDL_CTF_BG_SQRT_U);
  ADD_CONST(MDL_CTF_BG_SQRT_V);
  ADD_CONST(MDL_CTF_BG_SQRT_ANGLE);
  ADD_CONST(MDL_CTF_BG_BASELINE);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN2_K);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN2_SIGMAU);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN2_SIGMAV);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN2_CU);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN2_CV);
  ADD_CONST(MDL_CTF_BG_GAUSSIAN2_ANGLE);
  ADD_CONST(MDL_CTF_CRIT_PSDCORRELATION90);
  ADD_CONST(MDL_CTF_CRIT_FIRSTZERORATIO);
  ADD_CONST(MDL_CTF_CRIT_FIRSTZEROAVG);
  ADD_CONST(MDL_CTF_CRIT_FIRSTZERODISAGREEMENT);
  ADD_CONST(MDL_CTF_CRIT_DAMPING);
  ADD_CONST(MDL_CTF_CRIT_PSDRADIALINTEGRAL);
  ADD_CONST(MDL_CTF_CRIT_FITTINGSCORE);
  ADD_CONST(MDL_CTF_CRIT_FITTINGCORR13);
  ADD_CONST(MDL_CTF_CRIT_PSDVARIANCE);
  ADD_CONST(MDL_CTF_CRIT_PSDPCA1VARIANCE);
  ADD_CONST(MDL_CTF_CRIT_PSDPCARUNSTEST);
  ADD_CONST(MDL_DATATYPE);
  ADD_CONST(MDL_DEFGROUP);
  ADD_CONST(MDL_DM3_IDTAG);
  ADD_CONST(MDL_DM3_NODEID);
  ADD_CONST(MDL_DM3_NUMBER_TYPE);
  ADD_CONST(MDL_DM3_PARENTID);
  ADD_CONST(MDL_DM3_TAGCLASS);
  ADD_CONST(MDL_DM3_TAGNAME);
  ADD_CONST(MDL_DM3_SIZE);
  ADD_CONST(MDL_DM3_VALUE);
  ADD_CONST(MDL_ENABLED);
  ADD_CONST(MDL_EMX_MICROGRAPH_URL);
  ADD_CONST(MDL_EMX_MICROGRAPH_SAMPLING);
  ADD_CONST(MDL_EMX_MICROGRAPH_DEFOCUSU);
  ADD_CONST(MDL_EMX_MICROGRAPH_DEFOCUSV);
  ADD_CONST(MDL_EMX_MICROGRAPH_ASTIGMATISM_ANGLE);
  ADD_CONST(MDL_EMX_MICROGRAPH_VOLTAGE);
  ADD_CONST(MDL_EMX_MICROGRAPH_CS);
  ADD_CONST(MDL_EMX_MICROGRAPH_AMPLITUDE_CONTRAST);
  ADD_CONST(MDL_EMX_MICROGRAPH_FOM);
  ADD_CONST(MDL_EMX_PARTICLE_COORDINATE_X);
  ADD_CONST(MDL_EMX_PARTICLE_COORDINATE_Y);
  ADD_CONST(MDL_EMX_PARTICLE_URL);
  ADD_CONST(MDL_EMX_PARTICLE_MICROGRAPH_URL);
  ADD_CONST(MDL_EMX_PARTICLE_DEFOCUSU);
  ADD_CONST(MDL_EMX_PARTICLE_DEFOCUSV);
  ADD_CONST(MDL_EMX_PARTICLE_ASTIGMATISM_ANGLE);
  ADD_CONST(MDL_EMX_P_PARTICLE_CLASS_ID);
  ADD_CONST(MDL_EMX_P_PARTICLE_ACTIVE_FLAG);
  ADD_CONST(MDL_EMX_P_PARTICLE_FOM);
  ADD_CONST(MDL_EMX_P_PARTICLE_PARTICLE_URL);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_1_1);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_2_1);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_3_1);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_4_1);

  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_1_2);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_2_2);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_3_2);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_4_2);

  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_1_3);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_2_3);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_3_3);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_4_3);

  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_1_4);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_2_4);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_3_4);
  ADD_CONST(MDL_EMX_P_PARTICLE_TRANSFORMATION_MATRIX_4_4);

  ADD_CONST(MDL_FLIP);
  ADD_CONST(MDL_CLASS_COUNT);
  ADD_CONST(MDL_IMAGE);
  ADD_CONST(MDL_IMAGE_ORIGINAL);
  ADD_CONST(MDL_IMAGE_TILTED);
  ADD_CONST(MDL_IMGMD);
  ADD_CONST(MDL_INTSCALE);
  ADD_CONST(MDL_ITER);
  ADD_CONST(MDL_KSTEST);
  ADD_CONST(MDL_LL);
  ADD_CONST(MDL_MAGNIFICATION);
  ADD_CONST(MDL_MASK);
  ADD_CONST(MDL_MAXCC);
  ADD_CONST(MDL_MAX);
  ADD_CONST(MDL_MICROGRAPH);
  ADD_CONST(MDL_MICROGRAPH_ORIGINAL);
  ADD_CONST(MDL_MICROGRAPH_TILTED);
  ADD_CONST(MDL_MICROGRAPH_TILTED_ORIGINAL);
  ADD_CONST(MDL_MIN);
  ADD_CONST(MDL_MIRRORFRAC);
  ADD_CONST(MDL_MISSINGREGION_NR);
  ADD_CONST(MDL_MISSINGREGION_TYPE);
  ADD_CONST(MDL_MISSINGREGION_THY0);
  ADD_CONST(MDL_MISSINGREGION_THYF);
  ADD_CONST(MDL_MISSINGREGION_THX0);
  ADD_CONST(MDL_MISSINGREGION_THXF);
  ADD_CONST(MDL_MODELFRAC);
  ADD_CONST(MDL_NEIGHBORS);
  ADD_CONST(MDL_NEIGHBOR);
  ADD_CONST(MDL_NEIGHBORHOOD_RADIUS);
  ADD_CONST(MDL_NMA);

  ADD_CONST(MDL_NMA_MODEFILE);
  ADD_CONST(MDL_NOISE_ANGLES);
  ADD_CONST(MDL_NOISE_PARTICLE_COORD);
  ADD_CONST(MDL_NOISE_COORD);
  ADD_CONST(MDL_NOISE_PIXEL_LEVEL);

  ADD_CONST(MDL_CRYSTAL_LATTICE_A);
  ADD_CONST(MDL_CRYSTAL_LATTICE_B);
  ADD_CONST(MDL_CRYSTAL_DISAPPEAR_THRE);
  ADD_CONST(MDL_CRYSTAL_SHFILE);
  ADD_CONST(MDL_CRYSTAL_ORTHO_PRJ);
  ADD_CONST(MDL_CRYSTAL_PROJ);

  ADD_CONST(MDL_ORDER);
  ADD_CONST(MDL_ORIGIN_X);
  ADD_CONST(MDL_ORIGIN_Y);
  ADD_CONST(MDL_ORIGIN_Z);
  ADD_CONST(MDL_PICKING_FAMILY);
  ADD_CONST(MDL_COLOR);
  ADD_CONST(MDL_PICKING_PARTICLE_SIZE);
  ADD_CONST(MDL_PICKING_FAMILY_STATE);
  ADD_CONST(MDL_PICKING_MICROGRAPH_FAMILY_STATE);
  ADD_CONST(MDL_PMAX);

  ADD_CONST(MDL_PRJ_DIMENSIONS);
  ADD_CONST(MDL_PRJ_ANGFILE);
  ADD_CONST(MDL_PRJ_PSI_NOISE);
  ADD_CONST(MDL_PRJ_PSI_RANDSTR);
  ADD_CONST(MDL_PRJ_PSI_RANGE);
  ADD_CONST(MDL_PRJ_ROT_NOISE);
  ADD_CONST(MDL_PRJ_ROT_RANDSTR);
  ADD_CONST(MDL_PRJ_ROT_RANGE);
  ADD_CONST(MDL_PRJ_TILT_NOISE);
  ADD_CONST(MDL_PRJ_TILT_RANDSTR);
  ADD_CONST(MDL_PRJ_TILT_RANGE);
  ADD_CONST(MDL_PRJ_VOL);

  ADD_CONST(MDL_DIMENSIONS_3D);
  ADD_CONST(MDL_DIMENSIONS_2D);

  ADD_CONST(MDL_PSD);
  ADD_CONST(MDL_PSD_ENHANCED);
  ADD_CONST(MDL_RANDOMSEED);
  ADD_CONST(MDL_REF3D);
  ADD_CONST(MDL_REF);
  ADD_CONST(MDL_REFMD);
  ADD_CONST(MDL_RESOLUTION_DPR);
  ADD_CONST(MDL_RESOLUTION_ERRORL2);
  ADD_CONST(MDL_RESOLUTION_FRC);
  ADD_CONST(MDL_RESOLUTION_FRCRANDOMNOISE);
  ADD_CONST(MDL_RESOLUTION_FREQ);
  ADD_CONST(MDL_RESOLUTION_FREQREAL);
  ADD_CONST(MDL_SAMPLINGRATE);
  ADD_CONST(MDL_SAMPLINGRATE_ORIGINAL);
  ADD_CONST(MDL_SAMPLINGRATE_X);
  ADD_CONST(MDL_SAMPLINGRATE_Y);
  ADD_CONST(MDL_SAMPLINGRATE_Z);
  ADD_CONST(MDL_SCALE);
  ADD_CONST(MDL_SELFILE);
  ADD_CONST(MDL_SERIE);
  ADD_CONST(MDL_SHIFT_X);
  ADD_CONST(MDL_SHIFT_Y);
  ADD_CONST(MDL_SHIFT_Z);
  ADD_CONST(MDL_SHIFT_X2);
  ADD_CONST(MDL_SHIFT_Y2);
  ADD_CONST(MDL_SHIFT_X_DIFF);
  ADD_CONST(MDL_SHIFT_Y_DIFF);
  ADD_CONST(MDL_SHIFT_DIFF);
  ADD_CONST(MDL_CRYSTAL_SHIFTX);
  ADD_CONST(MDL_CRYSTAL_SHIFTY);
  ADD_CONST(MDL_CRYSTAL_SHIFTZ);
  ADD_CONST(MDL_SIGMANOISE);
  ADD_CONST(MDL_SIGMAOFFSET);
  ADD_CONST(MDL_SIGNALCHANGE);
  ADD_CONST(MDL_STDDEV);
  ADD_CONST(MDL_SUM);
  ADD_CONST(MDL_SUMWEIGHT);
  ADD_CONST(MDL_SYMNO);
  ADD_CONST(MDL_TRANSFORMATIONMTRIX);
  ADD_CONST(MDL_VOLTAGE);
  ADD_CONST(MDL_WEIGHT);
  ADD_CONST(MDL_WROBUST);
  ADD_CONST(MDL_XCOOR);
  ADD_CONST(MDL_XCOOR_TILT);
  ADD_CONST(MDL_XSIZE);
  ADD_CONST(MDL_X);
  ADD_CONST(MDL_YCOOR);
  ADD_CONST(MDL_YCOOR_TILT);
  ADD_CONST(MDL_Y);
  ADD_CONST(MDL_YSIZE);
  ADD_CONST(MDL_ZCOOR);
  ADD_CONST(MDL_Z);
  ADD_CONST(MDL_ZSCORE);
  ADD_CONST(MDL_LAST_LABEL);
  ADD_CONST(LABEL_NOTYPE);
  ADD_CONST(LABEL_INT);
  ADD_CONST(LABEL_BOOL);
  ADD_CONST(LABEL_DOUBLE);
  ADD_CONST(LABEL_VECTOR_DOUBLE);
  ADD_CONST(LABEL_STRING);
  ADD_CONST(LABEL_SIZET);
  ADD_CONST(LABEL_VECTOR_SIZET);
  ADD_CONST(_NONE);
  ADD_CONST(HEADER);
  ADD_CONST2("HEADER_ALL", _HEADER_ALL);
  ADD_CONST(DATA);
  ADD_CONST2("DATA_ALL", _DATA_ALL);
  ADD_CONST(WRAP);
  ADD_CONST(ALL_IMAGES);
  ADD_CONST(FILENAMENUMBERLENGTH);
  ADD_CONST2("XMIPP_BLACK", BLACK);
  ADD_CONST2("XMIPP_RED", RED);
  ADD_CONST2("XMIPP_GREEN", GREEN);
  ADD_CONST2("XMIPP_YELLOW", YELLOW);
  ADD_CONST2("XMIPP_BLUE", BLUE);
  ADD_CONST2("XMIPP_MAGENTA", MAGENTA);
  ADD_CONST2("XMIPP_CYAN", CYAN);
  ADD_CONST2("XMIPP_WHITE", WHITE);
  ADD_CONST2("XMIPP_RND_UNIFORM", RND_UNIFORM);
  ADD_CONST2("XMIPP_RND_GAUSSIAN", RND_GAUSSIAN);

  ADD_CONST2("DT_DEFAULT", DT_Default);
  ADD_CONST2("DT_UNKNOWN", DT_Unknown);
  ADD_CONST2("DT_UCHAR", DT_UChar);
  ADD_CONST2("DT_SCHAR", DT_SChar);
  ADD_CONST2("DT_USHORT", DT_UShort);
  ADD_CONST2("DT_SHORT", DT_Short);
  ADD_CONST2("DT_UINT", DT_UInt);
  ADD_CONST2("DT_INT", DT_Int);
  ADD_CONST2("DT_LONG", DT_Long);
  ADD_CONST2("DT_FLOAT", DT_Float);
  ADD_CONST2("DT_DOUBLE", DT_Double);
  ADD_CONST2("DT_COMPLEXSHORT", DT_CShort);
  ADD_CONST2("DT_COMPLEXINT", DT_CInt);
  ADD_CONST2("DT_COMPLEXFLOAT", DT_CFloat);
  ADD_CONST2("DT_COMPLEXDOUBLE", DT_CDouble);
  ADD_CONST2("DT_BOOL", DT_Bool);
  ADD_CONST2("DT_LASTENTRY", DT_LastEntry);

}
