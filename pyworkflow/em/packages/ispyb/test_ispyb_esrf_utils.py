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

import pprint
import unittest
from ispyb_esrf_utils import ISPyB_ESRF_Utils

class Test(unittest.TestCase):


    def tes_getMovieJpegMrcXml(self):
        movieFilename = "/users/svensson/cryoem/CWAT_ESRF_RawData_K2/170619_bGal1/Images-Disc1/GridSquare_19141127/Data/FoilHole_19150795_Data_19148847_19148848_20170619_2101-0344.mrc"
        jpeg, mrc, xml =  ISPyB_ESRF_Utils.getMovieJpegMrcXml(movieFilename)
        print(jpeg, mrc, xml)

    def tes_getAlignMoviesPngLogFilePath(self):
        mrcFilePath = "/users/svensson/cryoem/testRunData/000580_ProtMotionCorr/extra/FoilHole_19150795_Data_19148847_19148848_20170619_2101-0344_aligned_mic.mrc"
        png, logFilePath = ISPyB_ESRF_Utils.getAlignMoviesPngLogFilePath(mrcFilePath)
        print(png, logFilePath)
        
    def test_getXmlMetaData(self):
        xmlMetaDataFullPath = "/scisoft/data/cryoem/CWAT_ESRF_MicroscopePC/170619_bGal1/Images-Disc1/GridSquare_19141127/Data/FoilHole_19150716_Data_19148705_19148706_20170619_1744.xml"
        dictResults = ISPyB_ESRF_Utils.getXmlMetaData(xmlMetaDataFullPath)
        pprint.pprint(dictResults)
    

    def test_getCtfMetaData(self):
        workingDir = "/scisoft/pxsoft/data/cryoem/testRunData/20171017/000977_ProtGctf"
        mrcFilePath = "/scisoft/pxsoft/data/cryoem/testRunData/20171017/000859_ProtMotionCorr/extra/FoilHole_19150795_Data_19148847_19148848_20170619_2101-0344_aligned_mic.mrc"
        dictResults = ISPyB_ESRF_Utils.getCtfMetaData(workingDir, mrcFilePath)
        pprint.pprint(dictResults)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()