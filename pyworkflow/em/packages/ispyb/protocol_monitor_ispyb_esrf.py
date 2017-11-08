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
import sys
import pprint
import collections
import ConfigParser

sys.path.insert(0, "/opt/pxsoft/EDNA/vMX/edna/libraries/suds-0.4")

from suds.client import Client
from suds.transport.http import HttpAuthenticated

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_1
from pyworkflow.em.protocol import ProtMonitor, Monitor, PrintNotifier
from pyworkflow.em.protocol import ProtImportMovies, ProtAlignMovies, ProtCTFMicrographs
from pyworkflow.protocol import getProtocolFromDb

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
        self.beamlineName = "cm01"
        self.allParams = collections.OrderedDict()

    def getUpdatedProtocol(self, protocol):
        """ Retrieve the updated protocol and close db connections
            """
        prot2 = getProtocolFromDb(self.currentDir,
                                  protocol.getDbPath(),
                                  protocol.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def step(self):
        self.info("MonitorISPyB: start step")
                
        runs = [self.getUpdatedProtocol(p.get()) for p in self.protocol.inputProtocols] 
        
        g = self.project.getGraphFromRuns(runs)

        nodes = g.getRoot().iterChildsBreadth()

        for n in nodes:
            prot = n.run
            self.info("Protocol name: {0}".format(prot.getRunName()))
            self.info("Protocol: {0}".format(type(prot)))

            if isinstance(prot, ProtImportMovies):
                self.uploadImportMovies(prot)
            elif isinstance(prot, ProtAlignMovies) and hasattr(prot, 'outputMicrographs'):
                self.uploadAlignMovies(prot)
            elif isinstance(prot, ProtCTFMicrographs) and hasattr(prot, 'outputCTF'):
                self.uploadCTFMicrographs(prot)

        self.info("MonitorISPyB: end step")

        return False

    def iter_updated_set(self, objSet):
        objSet.load()
        objSet.loadAllProperties()
        for obj in objSet:
            yield obj
        objSet.close()


    def uploadImportMovies(self, prot):
        self.info("ESRF ISPyB upload import movies:")
        self.info("prot.getMatchFiles(): {0}".format(prot.getMatchFiles()))
        for movieFullPath in prot.getMatchFiles():
            listMovieFullPath = [ self.allParams[movieNumber]["movieFullPath"] for movieNumber in self.allParams if "movieFullPath" in self.allParams[movieNumber]]
            self.info("listMovieFullPath: {0}".format(listMovieFullPath))
            if not movieFullPath in listMovieFullPath:                
                self.info("movieFullPath: {0}".format(movieFullPath))
                dictFileNameParameters = ISPyB_ESRF_Utils.getMovieFileNameParameters(movieFullPath)
                self.movieDirectory = dictFileNameParameters["directory"]
                prefix = dictFileNameParameters["prefix"]   
                id1 = dictFileNameParameters["id1"]   
                id2 = dictFileNameParameters["id2"]   
                id3 = dictFileNameParameters["id3"]   
                date = dictFileNameParameters["date"]   
                hour = dictFileNameParameters["hour"]   
                movieNumber = dictFileNameParameters["movieNumber"]   
                self.info("movieNumber: {0}".format(movieNumber))
                suffix = dictFileNameParameters["suffix"]               
                
                self.movieDirectory = os.path.dirname(movieFullPath)
                
                micrographSnapshotFullPath, micrographFullPath, xmlMetaDataFullPath, gridSquareSnapshotFullPath = \
                   ISPyB_ESRF_Utils.getMovieJpegMrcXml(movieFullPath)
    
                self.info("proposal: {0}".format(self.proposalCode+self.proposalNumber))
                self.info("sampleAcronym: {0}".format(self.sampleAcronym))
                self.info("imageDirectory: {0}".format(self.movieDirectory))
                self.info("micrographSnapshotFullPath: {0}".format(micrographSnapshotFullPath))
    
                micrographSnapshotPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(micrographSnapshotFullPath)
                self.info("micrographSnapshotPyarchPath: {0}".format(micrographSnapshotPyarchPath))
                micrographPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(micrographFullPath)
                self.info("micrographPyarchPath: {0}".format(micrographPyarchPath))
                xmlMetaDataPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(xmlMetaDataFullPath)
                self.info("xmlMetaDataPyarchPath: {0}".format(xmlMetaDataPyarchPath))
                gridSquareSnapshotPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(gridSquareSnapshotFullPath)
                self.info("gridSquareSnapshotPyarchPath: {0}".format(gridSquareSnapshotPyarchPath))
                
                dictMetaData = ISPyB_ESRF_Utils.getXmlMetaData(xmlMetaDataFullPath)
                voltage = dictMetaData["accelerationVoltage"]
                magnification = dictMetaData["nominalMagnification"]
                imagesCount = dictMetaData["numberOffractions"]
                positionX = dictMetaData["positionX"]
                positionY = dictMetaData["positionY"]
                dosePerImage = dictMetaData["dose"]
                sphericalAberration = None
                amplitudeContrast = None
                scannedPixelSize = None
    
                self.info("voltage: {0}".format(voltage))
                self.info("sphericalAberration: {0}".format(sphericalAberration))
                self.info("amplitudeContrast: {0}".format(amplitudeContrast))
                self.info("magnification: {0}".format(magnification))
                self.info("scannedPixelSize: {0}".format(scannedPixelSize))
                self.info("imagesCount: {0}".format(imagesCount))
                self.info("dosePerImage: {0}".format(dosePerImage))
                self.info("positionX: {0}".format(positionX))
                self.info("positionY: {0}".format(positionY))
                self.info("beamlineName: {0}".format(self.beamlineName))
                
                
                movieObject = self.client.service.addMovie(proposal=self.proposalCode+self.proposalNumber, 
                                sampleAcronym=self.sampleAcronym, 
                                movieDirectory=self.movieDirectory,
                                movieFullPath=movieFullPath,
                                movieNumber=movieNumber,
                                micrographFullPath=micrographPyarchPath,
                                micrographSnapshotFullPath=micrographSnapshotPyarchPath,
                                xmlMetaDataFullPath=xmlMetaDataPyarchPath,
                                voltage=voltage,
                                sphericalAberration=sphericalAberration,
                                amplitudeContrast=amplitudeContrast,
                                magnification=magnification,
                                scannedPixelSize=scannedPixelSize,
                                imagesCount=imagesCount,
                                dosePerImage=dosePerImage,
                                positionX=positionX,
                                positionY=positionY,
                                beamlineName=self.beamlineName,
                                gridSquareSnapshotFullPath=gridSquareSnapshotPyarchPath,
                                )
                if movieObject is not None:
                    self.info("movieObject: {0}".format(pprint.pformat(dict(movieObject))))
                    movieId = movieObject.movieId
                else:
                    self.info("ERROR: movieObject is None!")
                    movieId = None
    
                self.allParams[movieNumber] = {
                    "movieFullPath": movieFullPath,
                    "prefix": prefix,   
                    "id1": id1,   
                    "id2": id2,   
                    "id3": id3,   
                    "date": date,   
                    "hour": hour,   
                    "movieId": movieId,   
                    "suffix": suffix,               
                 }
            self.info("self.allParams: {0}".format(self.allParams))

    def uploadAlignMovies(self, prot):
        self.info("ESRF ISPyB upload align movies:")
        self.info("prot.outputMicrographs: {0}".format(prot.outputMicrographs))
        for micrograph in self.iter_updated_set(prot.outputMicrographs):
            micrographFullPath = os.path.join(self.currentDir, micrograph.getFileName())
            self.info("micrographFullPath: {0}".format(micrographFullPath))
            dictFileNameParameters = ISPyB_ESRF_Utils.getMovieFileNameParameters(micrographFullPath)
            movieNumber = dictFileNameParameters["movieNumber"]
            if movieNumber in self.allParams and not "motionCorrectionId" in self.allParams[movieNumber]:
                movieFullPath = self.allParams[movieNumber]["movieFullPath"]
                dictResult = ISPyB_ESRF_Utils.getAlignMoviesPngLogFilePath(micrographFullPath)
                driftPlotFullPath = dictResult["globalShiftPng"]
                if "doseWeightMrc" in dictResult:
                    correctedDoseMicrographFullPath = dictResult["doseWeightMrc"]
                else:
                    correctedDoseMicrographFullPath = None
                if "thumbnailPng" in dictResult:
                    micrographSnapshotFullPath = dictResult["thumbnailPng"]
                else:
                    micrographSnapshotFullPath = None
                logFileFullPath = dictResult["logFileFullPath"]
                firstFrame = 1
                lastFrame = 28
                dosePerFrame = 1.0 
                doseWeight = 2.0
                totalMotion = 3.0
                averageMotionPerFrame = 4.0 
                driftPlotPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(driftPlotFullPath)
                micrographPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(micrographFullPath)
                correctedDoseMicrographPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(correctedDoseMicrographFullPath)
                micrographSnapshotPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(micrographSnapshotFullPath)
                logFilePyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(logFileFullPath)
                self.info("ESRF ISPyB upload align movies:")
                self.info("proposal: {0}".format(self.proposalCode+self.proposalNumber))
                self.info("movieFullPath: {0}".format(movieFullPath))
                self.info("firstFrame: {0}".format(firstFrame))
                self.info("lastFrame: {0}".format(lastFrame))
                self.info("dosePerFrame: {0}".format(dosePerFrame))
                self.info("doseWeight: {0}".format(doseWeight))
                self.info("totalMotion: {0}".format(totalMotion))
                self.info("averageMotionPerFrame: {0}".format(averageMotionPerFrame))
                self.info("driftPlotFullPath: {0}".format(driftPlotPyarchPath))
                self.info("micrographFullPath: {0}".format(micrographPyarchPath))
                self.info("correctedDoseMicrographFullPath: {0}".format(correctedDoseMicrographPyarchPath))
                self.info("micrographSnapshotFullPath: {0}".format(micrographSnapshotPyarchPath))
                self.info("logFilePath: {0}".format(logFilePyarchPath))
                motionCorrectionObject = self.client.service.addMotionCorrection(proposal=self.proposalCode+self.proposalNumber, 
                                                movieFullPath=movieFullPath,
                                                firstFrame=firstFrame,
                                                lastFrame=lastFrame,
                                                dosePerFrame=dosePerFrame,
                                                doseWeight=doseWeight,
                                                totalMotion=totalMotion,
                                                averageMotionPerFrame=averageMotionPerFrame,
                                                driftPlotFullPath=driftPlotPyarchPath,
                                                micrographFullPath=micrographPyarchPath,
                                                correctedDoseMicrographFullPath=correctedDoseMicrographPyarchPath,
                                                micrographSnapshotFullPath=micrographSnapshotPyarchPath,
                                                logFileFullPath=logFilePyarchPath)
                if motionCorrectionObject is not None:
                    self.info("motionCorrectionObject: {0}".format(pprint.pformat(dict(motionCorrectionObject))))
                    motionCorrectionId = motionCorrectionObject.motionCorrectionId
                else:
                    self.info("ERROR: motionCorrectionObject is None!")
                    motionCorrectionId = None
                self.allParams[movieNumber]["motionCorrectionId"] = motionCorrectionId

    def uploadCTFMicrographs(self, prot):
        self.info("ESRF ISPyB upload ctf micrographs:")
        self.info("prot.outputCTF: {0}".format(prot.outputCTF))
        workingDir = os.path.join(self.currentDir, str(prot.workingDir))
        for ctf in self.iter_updated_set(prot.outputCTF):
            micrographFullPath = ctf.getMicrograph().getFileName()
            dictFileNameParameters = ISPyB_ESRF_Utils.getMovieFileNameParameters(micrographFullPath)
            movieNumber = dictFileNameParameters["movieNumber"]
            if movieNumber in self.allParams and not "CTFid" in self.allParams[movieNumber]:
                movieFullPath = self.allParams[movieNumber]["movieFullPath"]
                dictResults = ISPyB_ESRF_Utils.getCtfMetaData(workingDir, micrographFullPath)
                spectraImageSnapshotFullPath = dictResults["spectraImageSnapshotFullPath"]
                spectraImageSnapshotPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(spectraImageSnapshotFullPath)
                spectraImageFullPath = dictResults["spectraImageFullPath"]
                spectraImagePyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(spectraImageFullPath)
                defocusU = dictResults["defocusU"]
                defocusV = dictResults["defocusV"]
                angle = dictResults["angle"]
                crossCorrelationCoefficient = dictResults["crossCorrelationCoefficient"]
                resolutionLimit = dictResults["resolutionLimit"]
                estimatedBfactor = dictResults["estimatedBfactor"]
                logFilePath = ISPyB_ESRF_Utils.copyToPyarchPath(dictResults["logFilePath"])
                self.info("micrographFullPath: {0}".format(micrographFullPath))
                self.info("ESRF ISPyB upload align movies:")
                self.info("proposal: {0}".format(self.proposalCode+self.proposalNumber))
                self.info("imageDirectory: {0}".format(self.movieDirectory))
                self.info("defocusU: {0}".format(dictResults["defocusU"]))
                self.info("defocusV: {0}".format(dictResults["defocusV"]))
                self.info("crossCorrelationCoefficient: {0}".format(dictResults["crossCorrelationCoefficient"]))
                self.info("resolutionLimit: {0}".format(dictResults["resolutionLimit"]))
                self.info("estimatedBfactor: {0}".format(dictResults["estimatedBfactor"]))
                self.info("spectraImageFullPath: {0}".format(spectraImagePyarchPath))
                self.info("spectraImageSnapshotFullPath: {0}".format(spectraImageSnapshotPyarchPath))
                self.info("logFilePath: {0}".format(dictResults["logFilePath"]))
                ctfObject = self.client.service.addCTF(proposal=self.proposalCode+self.proposalNumber, 
                                        movieFullPath=movieFullPath,
                                        spectraImageSnapshotFullPath=spectraImageSnapshotPyarchPath,
                                        spectraImageFullPath=spectraImagePyarchPath,
                                        defocusU=defocusU,
                                        defocusV=defocusV,
                                        angle=angle,
                                        crossCorrelationCoefficient=crossCorrelationCoefficient,
                                        resolutionLimit=resolutionLimit,
                                        estimatedBfactor=estimatedBfactor,
                                        logFilePath=logFilePath)
                if ctfObject is not None:
                    self.info("ctfObject: {0}".format(pprint.pformat(dict(ctfObject))))
                    CTFid = ctfObject.CTFid
                else:
                    self.info("ERROR: ctfObject is None!")
                    CTFid = None
                self.allParams[movieNumber]["CTFid"] = CTFid

