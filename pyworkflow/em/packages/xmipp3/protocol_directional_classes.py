# **************************************************************************
# *
# * Authors:     Javier Vargas (jvargas@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from os.path import join
from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        StringParam, EnumParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.utils.path import moveFile, makePath
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
import pyworkflow.em.metadata as md


PROJECTION_MATCHING = 0
SIGNIFICANT = 1


class XmippProtDirectionalClasses(ProtAnalysis3D):
    """    
    Analyze 2D classes as assigned to the different directions
    """

    _label = 'directional_classes'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",  
                      help='Select the input volume.')     
        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles', 
                      label="Input particles",  
                      help='Select the input projection images with an angular assignment.') 
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('targetResolution', FloatParam, default=8,
                      label='Target resolution (A)')
        
        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):        
        self._insertFunctionStep('convertInputStep', self.inputParticles.get().getObjId())
#         self._insertFunctionStep('constructGroups', self.inputParticles.get().getObjId())
#         self._insertFunctionStep('classifyGroups', self.inputParticles.get().getObjId(), self.inputVolume.get().getObjId())
#         self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        writeSetOfParticles(self.inputParticles.get(), self._getPath('input_particles.xmd'))
        Xdim = self.inputParticles.get().getDimensions()[0]
        Ts = self.inputParticles.get().getSamplingRate()
        newTs = self.targetResolution.get()*0.4
        newXdim = Xdim*Ts/newTs
        self.runJob("xmipp_image_resize","-i %s -o %s --save_metadata_stack %s --dim %d"%
                    self._getPath('input_particles.xmd'),
                    self._getExtraPath('scaled_particles.stk'),
                    self._getExtraPath('scaled_particles.xmd'),
                    newXdim)
    
   
    def createOutputStep(self):
        outputVols = self._createSetOfVolumes()
        
        for vol in self._iterInputVols():
            volume = vol.clone()
            volDir = self._getVolDir(vol.getObjId())
            volPrefix = 'vol%03d_' % (vol.getObjId())
            validationMd = self._getExtraPath(volPrefix + 'validation.xmd')
            moveFile(join(volDir, 'validation.xmd'), 
                     validationMd)
            clusterMd = self._getExtraPath(volPrefix + 'clusteringTendency.xmd')
            moveFile(join(volDir, 'clusteringTendency.xmd'), clusterMd)
            
            mData = md.MetaData(validationMd)
            weight = mData.getValue(md.MDL_WEIGHT, mData.firstObject())
            volume._xmipp_weight = Float(weight)
            volume.clusterMd = String(clusterMd)
            volume.cleanObjId() # clean objects id to assign new ones inside the set
            outputVols.append(volume)
        
        outputVols.setSamplingRate(self.partSet.getSamplingRate())
        self._defineOutputs(outputVolumes=outputVols)
        self._defineTransformRelation(self.inputVolumes, outputVols)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if self.inputVolumes.get() and not self.inputVolumes.hasValue():
            validateMsgs.append('Please provide an input reference volume.')
        if self.inputParticles.get() and not self.inputParticles.hasValue():
            validateMsgs.append('Please provide input particles.')            
        return validateMsgs
    
    def _summary(self):
        summary = []
        return summary
    
    #--------------------------- UTILS functions --------------------------------------------

    
    