# coding: utf-8
# **************************************************************************
# *
# * Author:     Olof Svensson (svensson@esrf.fr) [1]
# *
# * [1] European Synchrotron Radiation Facility
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

# This code is based on the "protocol_monitor_ispyb.py" written by
# J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es) [1] and
# Kevin Savage (kevin.savage@diamond.ac.uk) [2]
# [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# [2] Diamond Light Source, Ltd

import os
import collections
import ConfigParser

from suds.client import Client
from suds.transport.http import HttpAuthenticated

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_1
from pyworkflow.em.protocol import ProtMonitor, Monitor, PrintNotifier
from pyworkflow.em.protocol import ProtImportMovies, ProtAlignMovies, ProtCTFMicrographs

from ispyb_esrf_utils import ISPyB_ESRF_Utils

class ProtMonitorISPyB_ESRF(ProtMonitor):
    """ 
    Monitor to communicated with ISPyB system at ESRF.
    """
    _label = 'monitor to ISPyB at the ESRF'
    _lastUpdateVersion = VERSION_1_1

    def _defineParams(self, form):
        ProtMonitor._defineParams(self, form)

        group = form.addGroup('Experiment')
        group.addParam('proposal', params.StringParam,
                      label="Proposal",
                      help="Proposal")

        form.addParam('db', params.EnumParam,
                      choices=["production", "valid"],
                      label="Database",
                      help="Select which ISPyB database you want to use.")

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    #--------------------------- STEPS functions -------------------------------
    def monitorStep(self):
        config = ConfigParser.ConfigParser()
        credentialsConfig = ConfigParser.ConfigParser()
    
        # Configuration files
        config.read(os.path.join(os.path.dirname(__file__), 'ispyb.properties'))    
        credentialsConfig.read(os.path.join(os.path.dirname(__file__), 'credentials.properties'))
    
    
        username = str(credentialsConfig.get('Credential', 'user'))
        password = str(credentialsConfig.get('Credential', 'password'))
        url = str(config.get('Connection', 'url'))
    
        # Authentication
        httpAuthenticatedToolsForAutoprocessingWebService = HttpAuthenticated(username = username, password = password ) 
        client = Client( url, transport = httpAuthenticatedToolsForAutoprocessingWebService, cache = None, timeout = 15 )  
              
        proposalCode = config.get('Proposal', 'type')
        proposalNumber = config.get('Proposal', 'number')

        sampleAcronym = "ACRONYM"

        monitor = MonitorISPyB_ESRF(self, workingDir=self._getPath(),
                                        samplingInterval=self.samplingInterval.get(),
                                        monitorTime=100, client=client, 
                                        proposalCode=proposalCode, proposalNumber=proposalNumber,
                                        sampleAcronym=sampleAcronym)
    
        monitor.addNotifier(PrintNotifier())
        monitor.loop()


class MonitorISPyB_ESRF(Monitor):
    """ This will will be monitoring a CTF estimation protocol.
    It will internally handle a database to store produced
    CTF values.
    """
    def __init__(self, protocol, client=None, 
                 proposalCode=None, proposalNumber=None,
                 sampleAcronym=None, **kwargs):
        Monitor.__init__(self, **kwargs)
        self.protocol = protocol
        self.allIds = collections.OrderedDict()
        self.numberOfFrames = None
        self.imageGenerator = None
        self.proposal = self.protocol.proposal.get()
        self.project = self.protocol.getProject()
        self.client = client
        self.proposalCode = proposalCode
        self.proposalNumber = proposalNumber
        self.sampleAcronym = sampleAcronym
        self.movieDirectory = None
        self.currentDir = os.getcwd()

    def step(self):
        self.info("MonitorISPyB: only one step")

        prot = self.protocol

#         proxy = ISPyBProxy(["prod", "dev", "test"][prot.db.get()],
#                            experimentParams={'proposal': prot.proposal.get()})


        runs = [p.get() for p in self.protocol.inputProtocols]
        g = self.project.getGraphFromRuns(runs)

        nodes = g.getRoot().iterChildsBreadth()

        allParams = collections.OrderedDict()

        for n in nodes:
            prot = n.run
            self.info("Protocol name: {0}".format(prot.getRunName()))
            self.info("Protocol: {0}".format(type(prot)))

            if isinstance(prot, ProtImportMovies):
                self.uploadImportMovies(prot, allParams)
            elif isinstance(prot, ProtAlignMovies) and hasattr(prot, 'outputMicrographs'):
                self.uploadAlignMovies(prot, allParams)
            elif isinstance(prot, ProtCTFMicrographs):
                self.uploadCTFMicrographs(prot, allParams)

        for itemId, params in allParams.iteritems():
            self.info("Params: {0}".format(params))
#             ispybId = proxy.sendMovieParams(params)
            # Use -1 as a trick when ispyb is not really used and id is None
#             self.allIds[itemId] = ispybId or -1

        self.info("Closing proxy")
#         proxy.close()

        return False

    def iter_updated_set(self, objSet):
        objSet.load()
        objSet.loadAllProperties()
        for obj in objSet:
            yield obj
        objSet.close()


    def uploadImportMovies(self, prot, allParams):
        for movie in self.iter_updated_set(prot.outputMovies):
            movieFilePath = movie.getFileName()
            movieId = movie.getObjId()
            
            filesPath = prot.filesPath.get('').strip()
            self.movieDirectory = os.path.dirname(filesPath)
            
            jpeg, mrc, xml = ISPyB_ESRF_Utils.getMovieJpegMrcXml(filesPath)
            
            self.info("ESRF ISPyB upload import movies:")
            self.info("proposal: {0}".format(self.proposalCode+self.proposalNumber))
            self.info("sampleAcronym: {0}".format(self.sampleAcronym))
            self.info("imageDirectory: {0}".format(self.movieDirectory))
            self.info("jpeg: {0}".format(jpeg))
            self.info("mrc: {0}".format(mrc))
            self.info("xml: {0}".format(xml))
            
            self.client.service.addMovie(proposal=self.proposalCode+self.proposalNumber, 
                                        sampleAcronym=self.sampleAcronym, 
                                        imageDirectory=self.movieDirectory,
                                        jpeg=jpeg,
                                        mrc=mrc,
                                        xml=xml)

            allParams[movieId] = {
                'id': self.allIds.get(movieId, None),
                'imgdir': os.path.dirname(movieFilePath),
                'imgprefix': pwutils.removeBaseExt(movieFilePath),
                'imgsuffix': pwutils.getExt(movieFilePath),
                'file_template': movieFilePath,
                'n_images': self.numberOfFrames
             }
            self.info("allParams: {0}".format(allParams))

    def uploadAlignMovies(self, prot, allParams):
        self.info("allParams: {0}".format(dict(allParams)))
        self.info("prot.outputMicrographs: {0}".format(prot.outputMicrographs))
        for micrograph in self.iter_updated_set(prot.outputMicrographs):
            mrcFilePath = os.path.join(self.currentDir, micrograph.getFileName())
            png, logFilePath = ISPyB_ESRF_Utils.getAlignMoviesPngLogFilePath(mrcFilePath)
            jpeg = "/data/motionCorrJpegePath"
            self.info("ESRF ISPyB upload align movies:")
            self.info("proposal: {0}".format(self.proposalCode+self.proposalNumber))
            self.info("imageDirectory: {0}".format(self.movieDirectory))
            self.info("jpeg: {0}".format(jpeg))
            self.info("mrc: {0}".format(mrcFilePath))
            self.info("logFilePath: {0}".format(logFilePath))
            self.client.service.addMotionCorrection(proposal=self.proposalCode+self.proposalNumber, 
                                                    imageDirectory=self.movieDirectory,
                                                    jpeg=jpeg,
                                                    png=png,
                                                    mrc=mrcFilePath,
                                                    logFilePath=logFilePath)

    def uploadCTFMicrographs(self, prot, allParams):
        self.info("ESRF ISPyB upload ctf micrographs:")
        self.info("allParams: {0}".format(dict(allParams)))
        workingDir = os.path.join(self.currentDir, str(prot.workingDir))
        self.info("workingDir: {0}".format(workingDir))
        for ctf in self.iter_updated_set(prot.outputCTF):
            jpeg = "/data/ctfJpeg"
            mrcFilePath = ctf.getMicrograph().getFileName()
            self.info("mrcFilePath: {0}".format(mrcFilePath))
            ctfMrcFilePath, outputOne, outputTwo, logFilePath = ISPyB_ESRF_Utils.getCtfMetaData(workingDir, mrcFilePath)
            self.info("ESRF ISPyB upload align movies:")
            self.info("proposal: {0}".format(self.proposalCode+self.proposalNumber))
            self.info("imageDirectory: {0}".format(self.movieDirectory))
            self.info("jpeg: {0}".format(jpeg))
            self.info("mrc: {0}".format(ctfMrcFilePath))
            self.info("outputOne: {0}".format(outputOne))
            self.info("outputTwo: {0}".format(outputTwo))
            self.info("logFilePath: {0}".format(logFilePath))
            self.client.service.addCTF(proposal=self.proposalCode+self.proposalNumber, 
                                       imageDirectory=self.movieDirectory,
                                       jpeg=jpeg,
                                       mrc=ctfMrcFilePath,
                                       outputOne=outputOne,
                                       outputTwo=outputTwo,
                                       logFilePath=logFilePath)
