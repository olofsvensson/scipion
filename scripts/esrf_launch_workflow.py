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
import tempfile
from pyworkflow.manager import Manager
import pyworkflow.utils as pwutils


def usage(error):
    print """
    ERROR: %s

    Usage: scipion python scripts/esrf_launch_workflow.py
        dataDirectory="folder" location of raw data (movie files)
        filesPattern="pattern" template for movie files
        name="project name"
        sampleAcronym="Sample acronym name"
        This script will create and run a project 
    """ % error
    sys.exit(1)


n = len(sys.argv)

if n != 6:
    usage("Incorrect number of input parameters")

dataDirectory = sys.argv[1]
filesPattern = sys.argv[2]
projName = sys.argv[3]
proteinAcronym = sys.argv[4]
sampleAcronym = sys.argv[5]

path = os.path.join(os.environ['SCIPION_HOME'], 'pyworkflow', 'gui', 'no-tkinter')
sys.path.insert(1, path)

# Set up location
if "RAW_DATA" in dataDirectory:
    location = dataDirectory.replace("RAW_DATA", "PROCESSED_DATA")
    if not os.path.exists(location):
        os.makedirs(location, 0755)
else:
    location = tempfile.mkdtemp(prefix="ScipionUserData_")


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
        "voltage": 200.0,
        "sphericalAberration": 2.0,
        "amplitudeContrast": 0.1,
        "magnification": 50000,
        "samplingRateMode": 0,
        "samplingRate": 1.0,
        "scannedPixelSize": 7.0,
        "doseInitial": 0.0,
        "dosePerFrame": 1.0,
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
        "HighResL": 15.0,
        "HighResH": 4.0,
        "HighResBf": 50,
        "doValidate": false,
        "doPhShEst": false,
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
        "proposal": "mx415",
        "proteinAcronym": "%s",
        "sampleAcronym": "%s",
        "db": 1
    }
]""" % (dataDirectory, filesPattern, proteinAcronym, sampleAcronym)

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