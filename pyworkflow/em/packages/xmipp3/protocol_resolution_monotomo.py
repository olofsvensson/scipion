# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import (PointerParam, StringParam, BooleanParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from convert import readSetOfVolumes
from pyworkflow.object import Float
from pyworkflow.em import ImageHandler
from pyworkflow.utils import getExt
import numpy as np
import pyworkflow.em.metadata as md


MONORES_METHOD_URL = 'http://github.com/I2PC/scipion/wiki/XmippProtMonoRes'

OUTPUT_RESOLUTION_FILE = 'mgresolution.vol'
FN_FILTERED_MAP = 'filteredMap.vol'
OUTPUT_RESOLUTION_FILE_CHIMERA = 'MG_Chimera_resolution.vol'
OUTPUT_MASK_FILE = 'output_Mask.vol'
FN_MEAN_VOL = 'mean_volume.vol'
METADATA_MASK_FILE = 'mask_data.xmd'


class XmippProtMonoTomo(ProtAnalysis3D):
    """    
    Given a map the protocol assigns local resolutions to each voxel of the map.
    """
    _label = 'local MonoTomo'
    _version = VERSION_1_1
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.min_res_init = Float() 
        self.max_res_init = Float() 
        
    
    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('tomogram', PointerParam, pointerClass='Volume',
                      label="Tomogram", important=True, allowsNull = True,
                      help='Select a tomogram for determining its local resolution.')

        form.addParam('oddVolume', PointerParam, pointerClass='Volume',
                      label="Odd Tomogram", important=True,
                      help='Select a volume for determining its local resolution.')

        form.addParam('evenVolume', PointerParam, pointerClass='Volume',
                      label="Even Tomogram", important=True,
                      help='Select a second volume for determining a local resolution.')

        group = form.addGroup('Extra parameters')
        group.addParam('symmetry', StringParam, default='c1',
                      label="Symmetry",
                      help='Symmetry group. By default = c1.'
                      'See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]'
                      'for a description of the symmetry groups format, If no symmetry is present, give c1.')

        line = group.addLine('Resolution Range (A)',
                            help="If the user knows the range of resolutions or only a"
                                 " range of frequency needs to be analysed")
        
        group.addParam('significance', FloatParam, default=0.95, expertLevel=LEVEL_ADVANCED,
                      label="Significance",
                      help='Relution is computed using hipothesis tests, this value determines'
                      'the significance of that test')
        
        line.addParam('minRes', FloatParam, default=1, label='High')
        line.addParam('maxRes', FloatParam, default=30, label='Low')
        line.addParam('stepSize', FloatParam, allowsNull=True,
                      expertLevel=LEVEL_ADVANCED, label='Step')


    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        
        self.micsFn = self._getPath()
        self.vol1Fn = self.oddVolume.get().getFileName()
        self.vol2Fn = self.evenVolume.get().getFileName()
        
        if self.tomogram.hasValue():    
            self.vol0Fn = self.tomogram.get().getFileName()
        else:
            self.tomogram.set(None)


            # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep', )

        MS = self._insertFunctionStep('resolutionMonoTomoStep',
                                      prerequisites=[convertId])

        self._insertFunctionStep('createOutputStep', prerequisites=[MS])

        self._insertFunctionStep("createHistrogram")

    
    def convertInputStep(self):
        """ Read the input volume.
        """
        if self.tomogram.hasValue():    
            extVol0 = getExt(self.vol0Fn)
            if (extVol0 == '.mrc') or (extVol0 == '.map'):
                self.vol0Fn = self.vol0Fn + ':mrc'

        extVol1 = getExt(self.vol1Fn)
        extVol2 = getExt(self.vol2Fn)
        if (extVol1 == '.mrc') or (extVol1 == '.map'):
            self.vol1Fn = self.vol1Fn + ':mrc'
        if (extVol2 == '.mrc') or (extVol2 == '.map'):
            self.vol2Fn = self.vol2Fn + ':mrc'


    def resolutionMonoTomoStep(self):

#         #Number of frequencies
        if (self.stepSize.hasValue()):
            Nfreqs = round((self.maxRes.get() - self.minRes.get())/self.stepSize.get())
        else:
            Nfreqs = 50
  
              
        params =  ' --odd_volume %s' % self.vol1Fn
        params += ' --even_volume %s' % self.vol2Fn
        if self.tomogram.hasValue():
            params += ' --volume %s' % self.tomogram.get().getFileName()
        else:
            params += ' --volume '

        params += ' -o %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        params += ' --sampling_rate %f' % self.oddVolume.get().getSamplingRate()
        params += ' --number_frequencies %f' % Nfreqs
        params += ' --minRes %f' % self.minRes.get()
        params += ' --maxRes %f' % self.maxRes.get()
        params += ' --chimera_volume %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE_CHIMERA)
        params += ' --sym %s' % self.symmetry.get()
        params += ' --significance %f' % self.significance.get()
        params += ' --md_outputdata %s' % self._getExtraPath('mask_data.xmd')  

        self.runJob('xmipp_resolution_monotomo', params)


    def createHistrogram(self):

        params = ' -i %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        params += ' --steps %f' % 30
        params += ' --range %f %f' % (self.min_res_init, self.max_res_init)#(self.minRes.get(), self.maxRes.get())
        params += ' -o %s' % self._getExtraPath('hist.xmd')

        self.runJob('xmipp_image_histogram', params)
        
        
    def readMetaDataOutput(self):
        mData = md.MetaData(self._getExtraPath(METADATA_MASK_FILE))
        NvoxelsOriginalMask = float(mData.getValue(md.MDL_COUNT, mData.firstObject()))
        NvoxelsOutputMask = float(mData.getValue(md.MDL_COUNT2, mData.firstObject()))
        nvox = int(round(((NvoxelsOriginalMask-NvoxelsOutputMask)/NvoxelsOriginalMask)*100))
        return nvox

    def getMinMax(self, imageFile):
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        min_res = round(np.amin(imgData) * 100) / 100
        max_res = round(np.amax(imgData) * 100) / 100
        return min_res, max_res

    def createOutputStep(self):
        volume_path = self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        self.volumesSet = self._createSetOfVolumes('resolutionVol')
        if self.tomogram.hasValue():
            self.volumesSet.setSamplingRate(self.tomogram.get().getSamplingRate())
        else:
            self.volumesSet.setSamplingRate(self.oddVolume.get().getSamplingRate())
        readSetOfVolumes(volume_path, self.volumesSet)
        self._defineOutputs(outputVolume=self.volumesSet)
        if (self.tomogram.hasValue()):
            self._defineSourceRelation(self.tomogram, self.volumesSet)
        else:
            self._defineSourceRelation(self.oddVolume, self.volumesSet)
            
        #Setting the min max for the summary
        imageFile = self._getExtraPath(OUTPUT_RESOLUTION_FILE_CHIMERA)
        min_, max_ = self.getMinMax(imageFile)
        self.min_res_init.set(round(min_*100)/100)
        self.max_res_init.set(round(max_*100)/100)
        self._store(self.min_res_init)
        self._store(self.max_res_init)


    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'outputVolume'):
            messages.append(
                'Information about the method/article in ' + MONORES_METHOD_URL)
        return messages
    
    def _summary(self):
        summary = []
        summary.append("Highest resolution %.2f Å,   "
                       "Lowest resolution %.2f Å. \n" % (self.min_res_init,
                                                         self.max_res_init))
        Nvox = self.readMetaDataOutput()

        if (Nvox>10):
            summary.append("The resolution of %i %% of the mask voxels could not be computed. Maybe the mask was"
            "not correctly created, it is too wide or the resolution range does not cover the resolution at those voxels. "
            "If it is not the problem, decrease the significance in the advaced parameters can be an alternative" % Nvox)

        return summary

    def _citations(self):
        return ['Vilas2017']

