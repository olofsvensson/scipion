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


import os
import glob

class ISPyB_ESRF_Utils(object):

    @staticmethod
    def getMovieJpegMrcXml(movieFilePath):
        jpeg = None
        mrc = None
        xml = None
        doContinue = True
        currentDirectory = os.path.dirname(movieFilePath)
        movieFileName = os.path.basename(movieFilePath)
        listDirectory = []
        index = 0
        while doContinue:
            # Go up one directory level
            listDirectory.append(currentDirectory)
            #print(listDirectory)
            currentDirectory = os.path.dirname(currentDirectory)
            # Go through all sub directories
            for root, dirs, files in os.walk(currentDirectory):
                print("root: ", root)
                if not root in listDirectory:
                    for fileName in files:
                        fileNameWithoutSuffix = os.path.splitext(fileName)
                        if movieFileName.startswith(fileNameWithoutSuffix):
                            if fileName.endswith("jpg"):
                                jpeg = os.path.join(root, fileName)
                            elif fileName.endswith("mrc"):
                                mrc = os.path.join(root, fileName)
                            elif fileName.endswith("xml"):
                                xml = os.path.join(root, fileName)
                        if jpeg is not None and mrc is not None and xml is not None:
                            break
                listDirectory.append(root)
            index += 1
            if index > 4 or (jpeg is not None and mrc is not None and xml is not None):
                doContinue = False
            
        return jpeg, mrc, xml


    @staticmethod
    def getAlignMoviesPngLogFilePath(mrcFilePath):
        png = None
        logFilePath = None
        # Locate png file in same directory
        mrcDirectory = os.path.dirname(mrcFilePath)
        listPng = glob.glob(os.path.join(mrcDirectory, "*.png"))
        if len(listPng) > 1:
            print("WARNING! More than one PNG file found in ISPyB_ESRF_Utils.getAlignMoviesJpegMrcXml")
        elif len(listPng) > 0:
            png = listPng[0]
        # Find log file
        logFilePath = os.path.join(os.path.dirname(mrcDirectory), "logs", "run.log")
        return png, logFilePath
    
    @ staticmethod
    def getCtfMetaData(workingDir, mrcFilePath):
        ctfMrcFilePath = None
        outputOne = None
        outputTwo = None
        logFilePath = None
        # Find MRC directory
        mrcFileName = os.path.splitext(os.path.basename(mrcFilePath))[0]
        mrcDirectory = os.path.join(workingDir, "extra", mrcFileName)
        if os.path.exists(mrcDirectory):
            ctfMrcFilePath = os.path.join(mrcDirectory, "ctfEstimation.mrc")
            outputOne = os.path.join(mrcDirectory, "ctfEstimation.txt")
            outputTwo = os.path.join(mrcDirectory, "ctfEstimation_EPA.txt")
        # Find log file
        logFilePath = os.path.join(workingDir, "logs", "run.log")
        return ctfMrcFilePath, outputOne, outputTwo, logFilePath