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


    def tes_getMovieFileNameParameters(self):
        movieFullPath = "/users/svensson/cryoem/CWAT_ESRF_RawData_K2/170619_bGal1/Images-Disc1/GridSquare_19141127/Data/FoilHole_19150795_Data_19148847_19148848_20170619_2101-0344.mrc"
        dictResult = ISPyB_ESRF_Utils.getMovieFileNameParameters(movieFullPath)
        pprint.pprint(dictResult)
        self.assertEqual(dictResult["movieNumber"], "0344")
        mrcFullPath = "/users/svensson/ScipionUserData/projects/TestPipeLine/Runs/000859_ProtMotionCorr/extra/FoilHole_19150795_Data_19148847_19148848_20170619_2101-0344_aligned_mic.mrc"
        dictResult = ISPyB_ESRF_Utils.getMovieFileNameParameters(mrcFullPath)
        pprint.pprint(dictResult)
        
        

    def tes_getMovieJpegMrcXml(self):
#        movieFullPath = "/data/visitor/mx415/cm01/20171103/RAW_DATA/test4/CWAT_ESRF_RawData_K2/170619_bGal1/Images-Disc1/FoilHole_19150716_Data_19148705_19148706_20170619_1744-0055.mrc"
#        jpeg, mrc, xml, gridSquareThumbNail =  ISPyB_ESRF_Utils.getMovieJpegMrcXml(movieFullPath)
#        print(jpeg, mrc, xml, gridSquareThumbNail)
#        movieFullPath = "/data/visitor/mx415/cm01/20171108/RAW_DATA/test1/CWAT_ESRF_RawData_K2/170620_TMV_1/Images-Disc1/GridSquare_20174003/Data/FoilHole_20182354_Data_20179605_20179606_20170620_1523-1198.mrc"
#        jpeg, mrc, xml, gridSquareThumbNail =  ISPyB_ESRF_Utils.getMovieJpegMrcXml(movieFullPath)
#        print(jpeg, mrc, xml, gridSquareThumbNail)
        movieFullPath = "/data/cm01/inhouse/nicetest/supervisor_20171109_143730/Images-Disc1/GridSquare_9191979/Data/FoilHole_9208897_Data_9209280_9209281_20171109_1611-0556.mrc"
        jpeg, mrc, xml, gridSquareThumbNail =  ISPyB_ESRF_Utils.getMovieJpegMrcXml(movieFullPath)
        print(jpeg, mrc, xml, gridSquareThumbNail)
 
    def tes_getAlignMoviesPngLogFilePath(self):
        mrcFilePath = "/scisoft/pxsoft/data/cryoem/testRunData/20171113/Runs/000056_ProtMotionCorr/extra/FoilHole_9208892_Data_9209286_9209287_20171109_1540-0539_aligned_mic.mrc"
        dictResult = ISPyB_ESRF_Utils.getAlignMoviesPngLogFilePath(mrcFilePath)
        pprint.pprint(dictResult)
        
    def test_getShiftData(self):
        mrcFilePath = "/scisoft/pxsoft/data/cryoem/testRunData/20171113/Runs/000056_ProtMotionCorr/extra/FoilHole_9208892_Data_9209286_9209287_20171109_1540-0539_aligned_mic.mrc"
        dictResult = ISPyB_ESRF_Utils.getShiftData(mrcFilePath)
        pprint.pprint(dictResult)

    def tes_getXmlMetaData(self):
        xmlMetaDataFullPath = "/scisoft/data/cryoem/CWAT_ESRF_MicroscopePC/170619_bGal1/Images-Disc1/GridSquare_19141127/Data/FoilHole_19150716_Data_19148705_19148706_20170619_1744.xml"
        dictResults = ISPyB_ESRF_Utils.getXmlMetaData(xmlMetaDataFullPath)
        pprint.pprint(dictResults)
    

    def tes_getCtfMetaData(self):
        workingDir = "/scisoft/pxsoft/data/cryoem/testRunData/20171017/000977_ProtGctf"
        mrcFilePath = "/scisoft/pxsoft/data/cryoem/testRunData/20171017/000859_ProtMotionCorr/extra/FoilHole_19150795_Data_19148847_19148848_20170619_2101-0344_aligned_mic.mrc"
        dictResults = ISPyB_ESRF_Utils.getCtfMetaData(workingDir, mrcFilePath)
        pprint.pprint(dictResults)
        
    def tes_getPyarchFilePath(self):
        mrcFilePath = "/data/visitor/mx415/cm01/20171108/RAW_DATA/test2/CWAT_ESRF_RawData_K2/170620_TMV_1/Images-Disc1/GridSquare_20174003/Data/FoilHole_20182354_Data_20179605_20179606_20170620_1523-1198.mrc"
        pyarchFilePath = ISPyB_ESRF_Utils.getPyarchFilePath(mrcFilePath)
        self.assertEqual("/data/pyarch/2017/cm01/mx415/20171108/RAW_DATA/test2/CWAT_ESRF_RawData_K2/170620_TMV_1/Images-Disc1/GridSquare_20174003/Data/FoilHole_20182354_Data_20179605_20179606_20170620_1523-1198.mrc", pyarchFilePath)
        mrcFilePath = "/mntdirect/_data_visitor/mx415/cm01/20171108/PROCESSED_DATA/test2/CWAT_ESRF_RawData_K2/170620_TMV_1/Images-Disc1/GridSquare_20174003/Data/GridSquare_20174003_test2/Runs/000056_ProtMotionCorr/extra/FoilHole_20182354_Data_20179605_20179606_20170620_1523-1198_aligned_mic.mrc"
        pyarchFilePath = ISPyB_ESRF_Utils.getPyarchFilePath(mrcFilePath)
        self.assertEqual("/data/pyarch/2017/cm01/mx415/20171108/PROCESSED_DATA/test2/CWAT_ESRF_RawData_K2/170620_TMV_1/Images-Disc1/GridSquare_20174003/Data/GridSquare_20174003_test2/Runs/000056_ProtMotionCorr/extra/FoilHole_20182354_Data_20179605_20179606_20170620_1523-1198_aligned_mic.mrc", pyarchFilePath)
#        mrcFilePath = "/mntdirect/_data_cm01_inhouse/opcm01/20171108/RAW_DATA/nicetest/Frame.mrc"
#        pyarchFilePath = ISPyB_ESRF_Utils.getPyarchFilePath(mrcFilePath)
#        print(pyarchFilePath)
#        self.assertEqual("/data/pyarch/2017/cm01/mx415/20171108/PROCESSED_DATA/test2/CWAT_ESRF_RawData_K2/170620_TMV_1/Images-Disc1/GridSquare_20174003/Data/GridSquare_20174003_test2/Runs/000056_ProtMotionCorr/extra/FoilHole_20182354_Data_20179605_20179606_20170620_1523-1198_aligned_mic.mrc", pyarchFilePath)

    def tes_copyToPyarchPath(self):
        mrcFilePath = "/data/visitor/mx415/cm01/20171108/RAW_DATA/test2/CWAT_ESRF_RawData_K2/170620_TMV_1/Images-Disc1/GridSquare_20174003/Data/FoilHole_20182354_Data_20179605_20179606_20170620_1523-1198.mrc"
        pyarchFilePath = ISPyB_ESRF_Utils.copyToPyarchPath(mrcFilePath)
        print(pyarchFilePath)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()