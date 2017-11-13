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
import time
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

        group.addParam('sampleAcronym', params.StringParam,
                      label="Sample acronym",
                      help="Name of the sample acronym")

        group.addParam('proteinAcronym', params.StringParam,
                      label="Protein acronym",
                      help="Name of the protein acronym")

        form.addParam('db', params.EnumParam,
                      choices=["production", "valid"],
                      label="Database",
                      help="Select which ISPyB database you want to use.")

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')
        self._params = {}
        self._params['db'] = "valid"

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
        self.client = Client( url, transport = httpAuthenticatedToolsForAutoprocessingWebService, cache = None, timeout = 15 )  
              
        monitor = MonitorISPyB_ESRF(self, workingDir=self._getPath(),
                                        samplingInterval=self.samplingInterval.get(),
                                        monitorTime=100)
    
        monitor.addNotifier(PrintNotifier())
        monitor.loop()


class MonitorISPyB_ESRF(Monitor):
    """ This will will be monitoring a CTF estimation protocol.
    It will internally handle a database to store produced
    CTF values.
    """
    def __init__(self, protocol, **kwargs):
        Monitor.__init__(self, **kwargs)
        self.protocol = protocol
        self.allIds = collections.OrderedDict()
        self.numberOfFrames = None
        self.imageGenerator = None
        self.project = self.protocol.getProject()
        self.client = protocol.client
        self.proposal = protocol.proposal.get()
        self.proteinAcronym = protocol.proteinAcronym.get()
        self.sampleAcronym = protocol.sampleAcronym.get()
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
        self.info("MonitorISPyB: start step ------------------------")
                
        runs = [self.getUpdatedProtocol(p.get()) for p in self.protocol.inputProtocols] 
        
        g = self.project.getGraphFromRuns(runs)

        nodes = g.getRoot().iterChildsBreadth()

        for n in nodes:
            prot = n.run
            self.info("Protocol name: {0}".format(prot.getRunName()))

            if isinstance(prot, ProtImportMovies):
                self.uploadImportMovies(prot)
            elif isinstance(prot, ProtAlignMovies) and hasattr(prot, 'outputMicrographs'):
                self.uploadAlignMovies(prot)
            elif isinstance(prot, ProtCTFMicrographs) and hasattr(prot, 'outputCTF'):
                self.uploadCTFMicrographs(prot)

        self.info("MonitorISPyB: end step --------------------------")

        return False

    def iter_updated_set(self, objSet):
        objSet.load()
        objSet.loadAllProperties()
        for obj in objSet:
            yield obj
        objSet.close()


    def uploadImportMovies(self, prot):
        for movieFullPath in prot.getMatchFiles():
            listMovieFullPath = [ self.allParams[movieNumber]["movieFullPath"] for movieNumber in self.allParams if "movieFullPath" in self.allParams[movieNumber]]
            # self.info("listMovieFullPath: {0}".format(listMovieFullPath))
            if not movieFullPath in listMovieFullPath:                
                self.info("Import movies: movieFullPath: {0}".format(movieFullPath))
                dictFileNameParameters = ISPyB_ESRF_Utils.getMovieFileNameParameters(movieFullPath)
                self.movieDirectory = dictFileNameParameters["directory"]
                prefix = dictFileNameParameters["prefix"]   
                date = dictFileNameParameters["date"]   
                hour = dictFileNameParameters["hour"]   
                movieNumber = dictFileNameParameters["movieNumber"]   
                
                self.movieDirectory = os.path.dirname(movieFullPath)
                
                micrographSnapshotFullPath, micrographFullPath, xmlMetaDataFullPath, gridSquareSnapshotFullPath = \
                   ISPyB_ESRF_Utils.getMovieJpegMrcXml(movieFullPath)
                
                time.sleep(1)
                while micrographSnapshotFullPath is None or micrographFullPath is None or xmlMetaDataFullPath is None or gridSquareSnapshotFullPath is None:
                    self.info("Import movies: waiting for meta-data files to appear on disk...")
                    time.sleep(5)
                    micrographSnapshotFullPath, micrographFullPath, xmlMetaDataFullPath, gridSquareSnapshotFullPath = \
                       ISPyB_ESRF_Utils.getMovieJpegMrcXml(movieFullPath)
                time.sleep(1)
    
                self.info("Import movies: micrographSnapshotFullPath: {0}".format(micrographSnapshotFullPath))
    
                micrographSnapshotPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(micrographSnapshotFullPath)
                micrographPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(micrographFullPath)
                xmlMetaDataPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(xmlMetaDataFullPath)
                gridSquareSnapshotPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(gridSquareSnapshotFullPath)
                
                dictMetaData = ISPyB_ESRF_Utils.getXmlMetaData(xmlMetaDataFullPath)
                voltage = dictMetaData["accelerationVoltage"]
                magnification = dictMetaData["nominalMagnification"]
                imagesCount = dictMetaData["numberOffractions"]
                positionX = dictMetaData["positionX"]
                positionY = dictMetaData["positionY"]
                dosePerImage = round( float(dictMetaData["dose"]) / 10.0**20 / float(imagesCount), 2)
                sphericalAberration = None
                amplitudeContrast = None
                scannedPixelSize = None
                    
                movieObject = self.client.service.addMovie(proposal=self.proposal,
                                proteinAcronym=self.proteinAcronym, 
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
                    movieId = movieObject.movieId
                else:
                    self.info("ERROR: movieObject is None!")
                    movieId = None
    
                self.allParams[movieNumber] = {
                    "movieFullPath": movieFullPath,
                    "prefix": prefix,   
                    "date": date,   
                    "hour": hour,   
                    "movieId": movieId,   
                    "imagesCount": imagesCount,     
                    "dosePerFrame": prot.dosePerFrame.get(),
                 }
                self.info("Import movies done, movieId = {0}".format(self.allParams[movieNumber]["movieId"]))

    def uploadAlignMovies(self, prot):
        self.info("ESRF ISPyB upload ")
        for micrograph in self.iter_updated_set(prot.outputMicrographs):
            micrographFullPath = os.path.join(self.currentDir, micrograph.getFileName())
            dictFileNameParameters = ISPyB_ESRF_Utils.getMovieFileNameParameters(micrographFullPath)
            movieNumber = dictFileNameParameters["movieNumber"]
            if movieNumber in self.allParams and not "motionCorrectionId" in self.allParams[movieNumber]:
                self.info("Align movies: movie {0}".format(os.path.basename(self.allParams[movieNumber]["movieFullPath"])))
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
                lastFrame = self.allParams[movieNumber]["imagesCount"]
                dosePerFrame = self.allParams[movieNumber]["dosePerFrame"]
                doseWeight = "Dummy value: 2.0"
                totalMotion = "Dummy value: 3.0"
                averageMotionPerFrame = "Dummy value: 4.0" 
                driftPlotPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(driftPlotFullPath)
                micrographPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(micrographFullPath)
                correctedDoseMicrographPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(correctedDoseMicrographFullPath)
                micrographSnapshotPyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(micrographSnapshotFullPath)
                logFilePyarchPath = ISPyB_ESRF_Utils.copyToPyarchPath(logFileFullPath)
                motionCorrectionObject = self.client.service.addMotionCorrection(proposal=self.proposal, 
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
                    motionCorrectionId = motionCorrectionObject.motionCorrectionId
                else:
                    self.info("ERROR: motionCorrectionObject is None!")
                    motionCorrectionId = None
                self.allParams[movieNumber]["motionCorrectionId"] = motionCorrectionId
                self.info("Align movies done, motionCorrectionId = {0}".format(motionCorrectionId))

    def uploadCTFMicrographs(self, prot):
        workingDir = os.path.join(self.currentDir, str(prot.workingDir))
        for ctf in self.iter_updated_set(prot.outputCTF):
            micrographFullPath = ctf.getMicrograph().getFileName()
            dictFileNameParameters = ISPyB_ESRF_Utils.getMovieFileNameParameters(micrographFullPath)
            movieNumber = dictFileNameParameters["movieNumber"]
            if movieNumber in self.allParams and not "CTFid" in self.allParams[movieNumber]:
                self.info("CTF: movie {0}".format(os.path.basename(self.allParams[movieNumber]["movieFullPath"])))
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
                ctfObject = self.client.service.addCTF(proposal=self.proposal, 
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
                    CTFid = ctfObject.CTFid
                else:
                    self.info("ERROR: ctfObject is None!")
                    CTFid = None
                self.allParams[movieNumber]["CTFid"] = CTFid
                self.info("CTF done, CTFid = {0}".format(CTFid))

