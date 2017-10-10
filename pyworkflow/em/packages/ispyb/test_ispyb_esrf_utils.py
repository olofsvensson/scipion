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

    def test_getCtfMetaData(self):
        workingDir = "/users/svensson/cryoem/testRunData/000680_ProtGctf"
        mrcFilePath = "/users/svensson/cryoem/testRunData/000580_ProtMotionCorr/extra/FoilHole_19150795_Data_19148847_19148848_20170619_2101-0344_aligned_mic.mrc"
        ctfMrcFilePath, outputOne, outputTwo, logFilePath = ISPyB_ESRF_Utils.getCtfMetaData(workingDir, mrcFilePath)
        print(ctfMrcFilePath, outputOne, outputTwo, logFilePath)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()