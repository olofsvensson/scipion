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

from suds.client import Client
from suds.transport.http import HttpAuthenticated
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
	imageDirectory = "/data/imageDirectory"
	jpeg = "/data/jpeg"
	mrc = "/data/mrc"
	xml = "/data/xml"

	client.service.addMovie(proposal=proposalCode+proposalNumber, 
							sampleAcronym=sampleAcronym, 
							imageDirectory=imageDirectory,
							jpeg=jpeg,
							mrc=mrc,
							xml=xml)
	
	jpeg = "/data/motionCorrJpeg"
	mrc = "/data/motionCorrMrc"
	png = "/data/motionCorrPng"
	logFilePath = "/data/motionCorrLogFilePath"
	client.service.addMotionCorrection(proposal=proposalCode+proposalNumber, 
									imageDirectory=imageDirectory,
									jpeg=jpeg,
									png=png,
									mrc=mrc,
									logFilePath=logFilePath)
	
	jpeg = "/data/ctfJpeg"
	mrc = "/data/ctfMrc"
	outputOne = "/data/ctfOutputOne"
	outputTwo = "/data/ctfOutputTwo"
	logFilePath = "/data/ctfLogFilePath"
	client.service.addCTF(proposal=proposalCode+proposalNumber, 
							imageDirectory=imageDirectory,
							jpeg=jpeg,
							mrc=mrc,
							outputOne=outputOne,
							outputTwo=outputTwo,
							logFilePath=logFilePath)
	