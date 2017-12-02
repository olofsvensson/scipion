#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Olof Svensson (svensson@esrf.fr)
# *
# * European Synchrotron Radiation Facility
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

import os
import sys
import glob
import time
import pprint
import datetime
import tempfile
from pyworkflow.manager import Manager
import pyworkflow.utils as pwutils
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.packages.ispyb.ispyb_esrf_utils import ISPyB_ESRF_Utils

def usage(error):
    print """
    ERROR: %s

    Usage: scipion python scripts/esrf_launch_workflow.py
        dataDirectory="folder" location of raw data (movie files)
        filesPattern="pattern" template for movie files
        name="Scipion project name (must be unique)"
        proteinAcronym="Protein acronym name"
        sampleAcronym="Sample acronym name"
        doseInitial="Total dose"
        dosePerFrame="Dose per frame"
        This script will create and run a project 
    """ % error
    sys.exit(1)

def getUpdatedProtocol(protocol):
    """ Retrieve the updated protocol and close db connections
        """
    prot2 = getProtocolFromDb(os.getcwd(),
                              protocol.getDbPath(),
                              protocol.getObjId())
    # Close DB connections
    prot2.getProject().closeMapper()
    prot2.closeMappers()
    return prot2

n = len(sys.argv)

if n != 8:
    usage("Incorrect number of input parameters")

dataDirectory = sys.argv[1]
filesPattern = sys.argv[2]
projName = sys.argv[3]
proteinAcronym = sys.argv[4]
sampleAcronym = sys.argv[5]
doseInitial = float(sys.argv[6])
dosePerFrame = float(sys.argv[7])

path = os.path.join(os.environ['SCIPION_HOME'], 'pyworkflow', 'gui', 'no-tkinter')
sys.path.insert(1, path)

# Set up location
if "RAW_DATA" in dataDirectory:
    location = dataDirectory.replace("RAW_DATA", "PROCESSED_DATA")
    if not os.path.exists(location):
        os.makedirs(location, 0755)
else:
    location = tempfile.mkdtemp(prefix="ScipionUserData_")

# All param json file
allParamsJsonFile = os.path.join(location, "{0}_{1}.json".format(proteinAcronym,sampleAcronym))
#index = 1
#while os.path.exists(allParamsJsonFile):
#    allParamsJsonFile = os.path.join(location, "{0}_{1}_{2}.json".format(proteinAcronym,sampleAcronym, index))
#    index += 1

print("All param json file: {0}".format(allParamsJsonFile))

# Get meta data like phasePlateUsed

doPhaseShiftEstimation = "false"
listMovies = glob.glob(os.path.join(dataDirectory, filesPattern))

if len(listMovies) == 0:
    print("ERROR! No movies acqured yet.")
    sys.exit(1)

firstMovieFullPath = listMovies[0]

# Check proposal
# db=0: production
# db=1: valid
# db=2: lindemaria
# db=3: localhost
proposal = ISPyB_ESRF_Utils.getProposal(firstMovieFullPath)
if proposal is None:
    print("WARNING! No valid proposal could be found for movie {0}.".format(firstMovieFullPath))
    print("No data will be uploaded to ISPyB.")
    db = 3
else:
    print("Proposal: {0}".format(proposal))
    if proposal == "mx415":
        # Use valid data base
        db = 1
    else:
        # Use productiond data base
        db = 0

jpeg, mrc, xml, gridSquareThumbNail =  ISPyB_ESRF_Utils.getMovieJpegMrcXml(firstMovieFullPath)

dictResults = ISPyB_ESRF_Utils.getXmlMetaData(xml)
doPhaseShiftEstimation = dictResults["phasePlateUsed"]
nominalMagnification = int(dictResults["nominalMagnification"])
superResolutionFactor = int(dictResults["superResolutionFactor"])

samplingRate = 1.1 / float(superResolutionFactor) 
# Create json file

jsonString = """[
    {
        "object.className": "ProtImportMovies",
        "object.id": "2",
        "object.label": "scipion - import movies",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "importFrom": 0,
        "filesPath": "%s",
        "filesPattern": "%s",
        "copyFiles": false,
        "haveDataBeenPhaseFlipped": false,
        "acquisitionWizard": null,
        "voltage": 300.0,
        "sphericalAberration": 2.7,
        "amplitudeContrast": 0.1,
        "magnification": %d,
        "samplingRateMode": 0,
        "samplingRate": %f,
        "scannedPixelSize": 5.0,
        "doseInitial": %f,
        "dosePerFrame": %f,
        "gainFile": null,
        "darkFile": null,
        "dataStreaming": true,
        "timeout": 7200,
        "fileTimeout": 30,
        "inputIndividualFrames": false,
        "numberOfIndividualFrames": null,
        "stackFrames": false,
        "writeMoviesInProject": false,
        "movieSuffix": "_frames.mrcs",
        "deleteFrames": false,
        "streamingSocket": false,
        "socketPort": 5000
    },
    {
        "object.className": "ProtMotionCorr",
        "object.id": "77",
        "object.label": "motioncorr - motioncorr alignment",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "gpuMsg": "True",
        "GPUIDs": "0",
        "alignFrame0": 1,
        "alignFrameN": 0,
        "useAlignToSum": true,
        "sumFrame0": 1,
        "sumFrameN": 0,
        "binFactor": 1.0,
        "cropOffsetX": 0,
        "cropOffsetY": 0,
        "cropDimX": 0,
        "cropDimY": 0,
        "doSaveAveMic": true,
        "doSaveMovie": false,
        "doComputePSD": false,
        "doComputeMicThumbnail": true,
        "computeAllFramesAvg": false,
        "extraParams": "",
        "useMotioncor2": true,
        "doApplyDoseFilter": true,
        "patchX": 5,
        "patchY": 5,
        "group": 1,
        "tol": 0.5,
        "doMagCor": false,
        "useEst": true,
        "scaleMaj": 1.0,
        "scaleMin": 1.0,
        "angDist": 0.0,
        "extraParams2": "",
        "doSaveUnweightedMic": true,
        "hostName": "localhost",
        "numberOfThreads": 1,
        "numberOfMpi": 1,
        "inputMovies": "2.outputMovies"
    },
    {
        "object.className": "ProtGctf",
        "object.id": "195",
        "object.label": "gctf - CTF estimation on GPU",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "recalculate": false,
        "sqliteFile": null,
        "ctfDownFactor": 1.0,
        "lowRes": 0.05,
        "highRes": 0.35,
        "minDefocus": 0.25,
        "maxDefocus": 4.0,
        "astigmatism": 100.0,
        "windowSize": 512,
        "plotResRing": true,
        "GPUCore": 0,
        "doEPA": true,
        "EPAsmp": 4,
        "doBasicRotave": false,
        "bfactor": 150,
        "overlap": 0.5,
        "convsize": 85,
        "doHighRes": true,
        "HighResL": 30.0,
        "HighResH": 5.0,
        "HighResBf": 50,
        "doValidate": false,
        "doPhShEst": %s,
        "phaseShiftL": 0.0,
        "phaseShiftH": 180.0,
        "phaseShiftS": 10.0,
        "phaseShiftT": 0,
        "inputMicrographs": "77.outputMicrographs"
    },
    {
        "object.className": "ProtMonitorISPyB_ESRF",
        "object.id": "259",
        "object.label": "ispyb - monitor to ISPyB at the ESRF",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "inputProtocols": ["2", "77", "195"],
        "samplingInterval": 30,
        "proposal": "%s",
        "proteinAcronym": "%s",
        "sampleAcronym": "%s",
        "db": %d,
        "allParamsJsonFile": "%s"
    }
]""" % (dataDirectory, filesPattern, nominalMagnification, samplingRate, \
        doseInitial, dosePerFrame, doPhaseShiftEstimation, proposal, \
        proteinAcronym, sampleAcronym, db, allParamsJsonFile)

# Write json file
fd, jsonFile = tempfile.mkstemp(suffix=".json", prefix="scipion_workflow_")
os.write(fd, jsonString)
os.close(fd)
os.chmod(jsonFile, 0644)
print("jsonFile: {0}".format(jsonFile))


# Create a new project
manager = Manager()

if manager.hasProject(projName):
    usage("There is already a project with this name: %s"
          % pwutils.red(projName))

if jsonFile is not None and not os.path.exists(jsonFile):
    usage("Inexistent json file: %s" % pwutils.red(jsonFile))

project = manager.createProject(projName, location=location)

if jsonFile is not None:
    protDict = project.loadProtocols(jsonFile)
    
# Start the project
runs = project.getRuns()

# Now assuming that there is no dependencies between runs
# and the graph is lineal
#for prot in runs:
#    project.scheduleProtocol(prot)


# Monitor the execution:
doContinue = False
while doContinue:
    doContinue = False
    updatedRuns = [getUpdatedProtocol(p) for p in runs]
    print("") 
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')) 
    for prot in updatedRuns:
        print("{0} status: {1}".format(prot.getRunName(), prot.getStatusMessage()))
        if prot.isActive():
            doContinue = True
    time.sleep(15)
        
    