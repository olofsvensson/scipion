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
import re
import glob
import time
import socket
import shutil
import pprint
import xml.etree.ElementTree


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
        dictResult = {}
        # Locate png file in same directory
        mrcDirectory = os.path.dirname(mrcFilePath)
        listPng = glob.glob(os.path.join(mrcDirectory, "*.png"))
        for pngFile in listPng:
            dictFileNameParameters = ISPyB_ESRF_Utils.getMovieFileNameParameters(pngFile)
            if dictFileNameParameters["extra"] == "_global_shifts":
                dictResult["globalShiftPng"] = pngFile
            elif dictFileNameParameters["extra"] == "_thumbnail":
                dictResult["thumbnailPng"] = pngFile
        listMrc = glob.glob(os.path.join(mrcDirectory, "*.mrc"))
        for mrcFile in listMrc:
            dictFileNameParameters = ISPyB_ESRF_Utils.getMovieFileNameParameters(mrcFile)
            if "DW" in dictFileNameParameters["extra"]:
                dictResult["doseWeightMrc"] = mrcFile            
        # Find log file
        dictResult["logFileFullPath"] = os.path.join(os.path.dirname(mrcDirectory), "logs", "run.log")
        return dictResult
    
    @staticmethod
    def etree_to_dict(t):
        p = re.compile("^\{(.*)\}")
        m = p.match(t.tag)
        if m is not None:
            t.tag = t.tag[m.span()[1]:]
        listTmp = map(ISPyB_ESRF_Utils.etree_to_dict, t.getchildren())
        if len(listTmp) > 0:
            d = {t.tag : listTmp}
        else:
            d = {t.tag : t.text}
        return d
 
    @staticmethod
    def get_recursively(search_dict, field):
        """
        Takes a dict with nested lists and dicts,
        and searches all dicts for a key of the field
        provided.
        See: https://stackoverflow.com/questions/14962485/finding-a-key-recursively-in-a-dictionary
        """
        fields_found = []
    
        for key, value in search_dict.iteritems():
    
            if key == field:
                fields_found.append(value)
    
            elif isinstance(value, dict):
                results = ISPyB_ESRF_Utils.get_recursively(value, field)
                for result in results:
                    fields_found.append(result)
    
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        more_results = ISPyB_ESRF_Utils.get_recursively(item, field)
                        for another_result in more_results:
                            fields_found.append(another_result)
    
        return fields_found
    
    @staticmethod
    def getXmlMetaData(xmlMetaDataFullPath):
        dictResults = {}
        root = xml.etree.ElementTree.parse(xmlMetaDataFullPath).getroot()
        dictXML = ISPyB_ESRF_Utils.etree_to_dict(root)
        listKeyValue = ISPyB_ESRF_Utils.get_recursively(dictXML, "KeyValueOfstringanyType")
        for dictKey, dictValue in listKeyValue:
            if dictKey["Key"] == "Dose":
                dictResults["dose"] = dictValue["Value"]
        dictResults["numberOffractions"] = ISPyB_ESRF_Utils.get_recursively(dictXML, "NumberOffractions")[0]
        dictResults["nominalMagnification"] = ISPyB_ESRF_Utils.get_recursively(dictXML, "NominalMagnification")[0]
        dictResults["positionX"] = ISPyB_ESRF_Utils.get_recursively(dictXML, "X")[0]
        dictResults["positionY"] = ISPyB_ESRF_Utils.get_recursively(dictXML, "Y")[0]
        dictResults["accelerationVoltage"] = ISPyB_ESRF_Utils.get_recursively(dictXML, "AccelerationVoltage")[0]
        dictResults["acquisitionDateTime"] = ISPyB_ESRF_Utils.get_recursively(dictXML, "acquisitionDateTime")[0]
        return dictResults

    @staticmethod
    def _getKeyValue(root, key):
        value = None
        keyFound = False
        for child in root:
            for subChild in child:
                if subChild.tag.endswith("Key") and subChild.text == key:
                    keyFound = True
                if keyFound and subChild.tag.endswith("Value"):
                    value = subChild.text
                    keyFound = False
                    break
        return value
        
    
    
    @ staticmethod
    def getCtfMetaData(workingDir, mrcFilePath):
        dictResults = {
            "spectraImageSnapshotFullPath": None,
            "spectraImageFullPath": None,
            "defocusU": None,
            "defocusV": None,
            "angle": None,
            "crossCorrelationCoefficient": None,
            "resolutionLimit": None,
            "estimatedBfactor": None,
            "logFilePath": None,
        }
        # Find MRC directory
        mrcFileName = os.path.splitext(os.path.basename(mrcFilePath))[0]
        mrcDirectory = os.path.join(workingDir, "extra", mrcFileName)
        if os.path.exists(mrcDirectory):
            spectraImageFullPath = os.path.join(mrcDirectory, "ctfEstimation.mrc")
            if os.path.exists(spectraImageFullPath):
                dictResults["spectraImageFullPath"] = spectraImageFullPath
                spectraImageSnapshotFullPath = os.path.join(mrcDirectory, "ctfEstimation.jpeg")
                os.system("bimg {0} {1}".format(spectraImageFullPath, spectraImageSnapshotFullPath))
                if os.path.exists(spectraImageSnapshotFullPath):
                    dictResults["spectraImageSnapshotFullPath"] = spectraImageSnapshotFullPath
            ctfEstimationPath = os.path.join(mrcDirectory, "ctfEstimation.txt")
            if os.path.exists(ctfEstimationPath):
                f = open(ctfEstimationPath)
                lines = f.readlines()
                f.close()
                index = 0
                for index in range(len(lines)):
                    if "Defocus_U" in lines[index]:
                        index += 1
                        listValues = lines[index].split()
                        dictResults["defocusU"] = listValues[0]
                        dictResults["defocusV"] = listValues[1]
                        dictResults["angle"] = listValues[2]
                        dictResults["crossCorrelationCoefficient"] = listValues[3]
                    elif "Resolution limit" in lines[index]:
                        listValues = lines[index].split()
                        print(listValues)
                        dictResults["resolutionLimit"] = listValues[-2]
                    elif "Estimated Bfactor" in lines[index]:
                        listValues = lines[index].split()
                        dictResults["estimatedBfactor"] = listValues[-2]
        # Find log file
        logFilePath = os.path.join(workingDir, "logs", "run.log")
        if os.path.exists(logFilePath):
            dictResults["logFilePath"] = logFilePath
        return dictResults

    @ staticmethod
    def getMovieFileNameParameters(mrcFilePath):
        """
        FoilHole_19150795_Data_19148847_19148848_20170619_2101-0344.mrc
        """
        dictResult = {}
        p = re.compile("^(.*)/(.*)_([0-9]*)_Data_([0-9]*)_([0-9]*)_([0-9]*)_([0-9]*)-([0-9]*)(_?.*)\.(.*)")
        m = p.match(mrcFilePath)
        dictResult["directory"] = m.group(1)   
        dictResult["prefix"] = m.group(2)   
        dictResult["id1"] = m.group(3)   
        dictResult["id2"] = m.group(4)   
        dictResult["id3"] = m.group(5)   
        dictResult["date"] = m.group(6)   
        dictResult["hour"] = m.group(7)   
        dictResult["movieNumber"] = m.group(8)   
        dictResult["extra"] = m.group(9)   
        dictResult["suffix"] = m.group(10)   
        return dictResult

    @ staticmethod
    def copyToPyarchPath(filePath):
        # Test path:
        testPath = "/data/pyarch/2017/cm01/test"
        # Add date
        datePath = os.path.join(testPath, time.strftime("%Y%m%d", time.localtime(time.time())))
        # Loop until done
        isDone = False
        fileName = os.path.basename(filePath)
        pyarchFilePath = None
        while not isDone:
            timePath = os.path.join(datePath, time.strftime("%H%M%S", time.localtime(time.time())))
            print(timePath)
            if os.path.exists(timePath):
                time.sleep(1)
            else:
                pyarchFilePath = os.path.join(timePath, fileName)
                if "linsvensson" in socket.gethostname() and os.path.getsize(filePath) < 1e6:
                    # For the moment, only copy file if smaller than 1 MB
                    os.system("ssh mxhpc2-1705 'mkdir -p {0}'".format(timePath))
                    os.system("scp {0} mxhpc2-1705:{1}".format(filePath, pyarchFilePath))
                else:
                    os.makedirs(timePath, 0755)
                    shutil.copy(filePath, pyarchFilePath)
                isDone = True
        else:
            pyarchFilePath = filePath
        return pyarchFilePath
        
                