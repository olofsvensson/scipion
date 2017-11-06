# coding: utf-8
# **************************************************************************
# *
# * Authors:     Alejandro de Maria (demariaa@esrf.fr) [1]
# *              Olof Svensson (svensson@esrf.fr) [1]
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

import datetime

from suds.client import Client
from suds.transport.http import HttpAuthenticated
from suds.sax.date import DateTime


import ConfigParser


if __name__ == "__main__":

	config = ConfigParser.ConfigParser()
	credentialsConfig = ConfigParser.ConfigParser()

	# Configuration files
	config.read('ispyb.properties')	
	credentialsConfig.read('credentials.properties')

	username = str(credentialsConfig.get('Credential', 'user'))
	password = str(credentialsConfig.get('Credential', 'password'))
	url = str(config.get('Connection', 'url'))

	# Authentication
	httpAuthenticatedToolsForAutoprocessingWebService = HttpAuthenticated(username = username, password = password ) 
	client = Client( url, transport = httpAuthenticatedToolsForAutoprocessingWebService, cache = None, timeout = 15 )
	
	# Proposal parameters
	proposalCode = config.get('Proposal', 'type')
	proposalNumber = config.get('Proposal', 'number')


	sampleAcronym = "ACRONYM"
	movieDirectory = "/data/imageDirectory2"
	movieFullPath = "/users/svensson/cryoem/CWAT_ESRF_RawData_K2/170619_bGal1/Images-Disc1/GridSquare_19141127/Data/FoilHole_19150795_Data_19148847_19148848_20170619_2101-0344.mrc"
	movieNumber = "301"
	micrographFullPath = "/data/pyarch/fileName.mrc"
	micrographSnapshotFullPath = "/data/pyarch/fileName.png"
	xmlMetaDataFullPath = "/data/pyarch/fileName.png"
	voltage = "200.0"
	sphericalAberration = "2.0"
	amplitudeContrast = "0.1"
	magnification = "5000"
	scannedPixelSize = "7.0"
	imagesCount = "28"
	dosePerImage = "1.0"
	positionX = "101"
	positionY = "102"
	beamlineName = "cm01"			
	if True:								
		movieObject = client.service.addMovie(proposal=proposalCode+proposalNumber, 
							sampleAcronym=sampleAcronym, 
							movieDirectory=movieDirectory,
							movieFullPath=movieFullPath,
							movieNumber=movieNumber,
							micrographFullPath=micrographFullPath,
							micrographSnapshotFullPath=micrographSnapshotFullPath,
							xmlMetaDataFullPath=xmlMetaDataFullPath,
							voltage=voltage,
							sphericalAberration=sphericalAberration,
							amplitudeContrast=amplitudeContrast,
							magnification=magnification,
							scannedPixelSize=scannedPixelSize,
							imagesCount=imagesCount,
							dosePerImage=dosePerImage,
							positionX=positionX,
							positionY=positionY,
							beamlineName=beamlineName,
							)
	
		print(movieObject)
	firstFrame = 1
	lastFrame = 28
	dosePerFrame = 1.0 
	doseWeight = 2.0
	totalMotion = 3.0
	averageMotionPerFrame = 4.0 
	driftPlotFullPath = "/data/pyarch/motionCorr.jpeg"
	micrographFullPath = "/data/pyarch/motionCorr.mrc"
	correctedDoseMicrographFullPath = "/data/pyarch/motionCorr_DW.mrc"
	micrographSnapshotFullPath = "/data/pyarch/motionCorr.png"
	logFileFullPath = "/data/pyarch/motionCorrLog.txt"

	
	if True:								
		motionCorrectionObject = client.service.addMotionCorrection(proposal=proposalCode+proposalNumber, 
											movieFullPath=movieFullPath,
											firstFrame=firstFrame,
											lastFrame=lastFrame,
											dosePerFrame=dosePerFrame,
											doseWeight=doseWeight,
											totalMotion=totalMotion,
											averageMotionPerFrame=averageMotionPerFrame,
											driftPlotFullPath=driftPlotFullPath,
											micrographFullPath=micrographFullPath,
											correctedDoseMicrographFullPath=correctedDoseMicrographFullPath,
											micrographSnapshotFullPath=micrographSnapshotFullPath,
											logFileFullPath=logFileFullPath)
		print(motionCorrectionObject)
		
	
	spectraImageSnapshotFullPath = "/data/pyarch/ctf.jpeg"
	spectraImageFullPath = "/data/pyarch/ctf.mrc"
	defocusU = "1"
	defocusV = "2"
	angle = "22"
	crossCorrelationCoefficient = "3"
	resolutionLimit = "4"
	estimatedBfactor = "5"
	logFilePath = "/data/pyarch/ctfLogFile.txt"
	if True:
		ctfObject = client.service.addCTF(proposal=proposalCode+proposalNumber, 
							movieFullPath=movieFullPath,
							spectraImageSnapshotFullPath=spectraImageSnapshotFullPath,
							spectraImageFullPath=spectraImageFullPath,
							defocusU=defocusU,
							defocusV=defocusV,
							angle=angle,
							crossCorrelationCoefficient=crossCorrelationCoefficient,
							resolutionLimit=resolutionLimit,
							estimatedBfactor=estimatedBfactor,
							logFilePath=logFilePath)
		print(ctfObject)
