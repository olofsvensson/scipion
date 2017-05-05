# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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


from protocol_refmac import CCP4ProtRunRefmac
from pyworkflow.protocol.params import LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer, Viewer
from pyworkflow.gui.text import _open_cmd
from pyworkflow.em.data import EMSet, EMObject
from pyworkflow.object import Float, String
from pyworkflow.em.viewer import ObjectView, TableView
from tkMessageBox import showerror

from pyworkflow.em.viewer import ImageView, ChimeraView
import os


def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
    except:
        print("Error:", msg)


class ParseFile():
    """class that parse a refmac log file"""
    LASTITERATIONRESULTS="lastIteration"
    FINALRESULTS="finalResults"
    FOMPLOT="fomplot"
    def __init__(self, fileName, tkParent=None, lastIteration=0):
        self.headerDict = {}#parsed headers goes here
        self.dataDict = {}#parsed data goes here
        self.msgDict = {}#error messages goes here
        self.tkParent = tkParent
        self.fileName = fileName
        self._parsefile(lastIteration)

    def _parseFinalResults(self, filePointer):
        headerList = []
        dataList = []
        stop = False
        msg=""
        while 1:
            line = filePointer.readline()
            if line.find('$TEXT:Result: $$ Final results $$') != -1:  # detect final results
                break
            if not line:
                stop = True
                break
        if stop:
            msg = 'Can not find "Final result" information in log file: %s'% self.fileName
        else:
            # finalResultsDict={'header':}
            #parse header
            headerList.append(" ")
            line = filePointer.readline()
            words = line.strip().split()
            headerList.extend([words[0],words[1]])
            #parse data: another 4 lines
            for i in range(4):
                row = []
                line = filePointer.readline()
                words = line.strip().split()
                #the first column has 2 words
                row.extend([words[0]+" "+ words[1], words[2], words[3]])
                dataList.append(tuple(row))
            #TODO: remove debug lines
        return headerList, dataList, msg

    def retrievefinalResults(self):
        return self.headerDict[self.FINALRESULTS], \
               self.dataDict[self.FINALRESULTS], \
               self.msgDict[self.FINALRESULTS]

    def _parseLastIteration(self, filePointer, iteration):
        headerList = ["variable","value"]
        dataList = []
        stop = False
        msg=""
        while 1:
            line = filePointer.readline()
            #TODO: double check space after Cycle
            if line.find('$GRAPHS:Cycle    %d. M(Fom) v. resln :N:1,3,5,7,8,9,10:'%iteration)!=-1:  # detect final results
                break
            if not line:
                stop = True
                break
        if stop:
            msg = 'Can not find "Last Iteration" information in log file: %s' % self.fileName
        else:
            # find three lines with $$
            counter=4
            while counter!=0:
                line = filePointer.readline()
                if line.find("$$")!=-1:
                    counter -= 1
                if not line:
                    stop = True
                    break

            for i in range(14):
                row = []
                line = filePointer.readline()
                words = line.strip().split("=")
                # the first column has 2 words
                row.extend([words[0].strip(), words[1].strip()])
                dataList.append(tuple(row))
        return headerList, dataList, msg

    def retrievelastIteration(self):
        return self.headerDict[self.LASTITERATIONRESULTS],\
               self.dataDict[self.LASTITERATIONRESULTS],\
               self.msgDict[self.LASTITERATIONRESULTS]

    def _parsefile(self, lastIteration=0):
        """ call the different functions that parse the data in the right order"""
        with open(self.fileName,"r") as filePointer:
            headerList,dataList, msg = self._parseLastIteration(filePointer, lastIteration)
            self.headerDict[self.LASTITERATIONRESULTS] = headerList
            self.dataDict[self.LASTITERATIONRESULTS]   = dataList
            self.msgDict[self.LASTITERATIONRESULTS]    =  msg
            headerList,dataList, msg = self._parseLastIteration(filePointer, lastIteration)
            self.headerDict[self.LASTITERATIONRESULTS] = headerList
            self.dataDict[self.LASTITERATIONRESULTS]   = dataList
            self.msgDict[self.LASTITERATIONRESULTS]    =  msg

class CCP4ProtRunRefmacViewer(ProtocolViewer):
    """ Viewer for CCP4 program refmac
    """
    _label = 'Refmac Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [CCP4ProtRunRefmac]

    # ROB: do we need this memory for something?
    # _memory = False
    # temporary metadata file with ctf that has some resolution greathan than X
    # tmpMetadataFile = 'viewersTmp.sqlite'
    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        self.parseFile = ParseFile(self.protocol._getExtraPath(self.protocol.refineLogFileName),
                                   self.getTkRoot(),
                                   self.protocol.nRefCycle.get() + 1)
    def _defineParams(self, form):
        form.addSection(label='Visualization of Refmac results')
        # group = form.addGroup('Overall results')
        form.addParam('displayMask', LabelParam,
                      label="PDB based Mask",
                      help="Display Masked map")
        form.addParam('showFinalResults', LabelParam,
                      label="Final Results Table",
                      help="Table of Final Results from refine.log file")
        form.addParam('showLogFile', LabelParam,
                      label="Show log file",
                      help="open refmac log file in a text editor")
        form.addParam('showLastIteration', LabelParam,
                      label="Results Table (last iteration)",
                      help="Table stored in log file summarizing the last iteration")
        form.addParam('displayRFactorPlot', LabelParam,
                      label="R-factor vs. iteration",
                      help="Plot R-factor as a function of the iteration")
        form.addParam('displayFOMPlot', LabelParam,
                      label="FOM vs. iteration",
                      help="Plot Figure Of Merit as a function of the iteration")
        form.addParam('displayLLPlot', LabelParam,
                      label="-LL vs. iteration",
                      help="Plot Log likelihood as a function of the iteration")
        form.addParam('displayLLfreePlot', LabelParam,
                      label="-LLfree vs. iteration",
                      help="Plot Log likelihood as a function of the iteration")
        form.addParam('displayGeometryPlot', LabelParam,
                      label="Geometry vs. iteration",
                      help="""Plot Geometry as a function of the iteration:
Geometry includes rmsBOND (root mean square bond lengths)
zBOND (zscore of the deviation of bond lengths)
rmsANGL (root mean square bond angles)
zANGL (zscore of the deviation of bond angles)
and rmsCHIRAL (root mean square of chiral index""")

    def _getVisualizeDict(self):
        return {
            'showFinalResults': self._visualizeFinalResults,
            'showLastIteration': self._visualizeLastIteration,
            'displayMask': self._visualizeMask,
            'displayRFactorPlot': self._visualizeRFactorPlot,
            'displayFOMPlot': self._visualizeFOMPlot,
            'displayLLPlot': self._visualizeLLPlot,
            'displayLLfreePlot': self._visualizeLLfreePlot,
            'displayGeometryPlot': self._visualizeGeometryPlot,
            'showLogFile': self._visualizeLogFile
        }

    def _visualizeMask(self):
        pass

    def _visualizeFinalResults(self, e=None):


        """
        views = []
        labels = '_1 _2'
        emSet = EMSet(filename="/tmp/kk.sqlite")
        emObject = EMObject()
        emObject._1 = String('first parameter')
        emObject._2 = Float(12.)
        emSet.append(emObject)
        emObject = EMObject()
        emObject._1 = String('second parameter')
        emObject._2 = Float(22.)
        emSet.append(emObject)
        emSet.write()
        views.append(ObjectView(self._project,
                                self.protocol.strId(),
                                "/tmp/kk.sqlite",
                                viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels}))
        return views
"""
        #Selection of lines from 'refine.log' file that include Refmac final results.
        #These lines will be saved in outputLines list

        headerList, dataList, msg = self.parseFile.retrievefinalResults()
        if not dataList:
            errorWindow(self.getTkRoot(), msg)
            return

        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg="Values for a good fitted 3D map. R factor ~ 0.3, Rms BondLength ~ 0.02.",
                  title= "Refmac: Final Results Summary",
                  height=len(dataList), width=200,padding=40)

    def _visualizeLogFile(self, e=None):
        """Show refmac log file."""
        refineLogFileName = self.protocol._getExtraPath(self.protocol.refineLogFileName)
        _open_cmd(refineLogFileName, self.getTkRoot())

    def _visualizeLastIteration(self, e=None):
        # Selection of lines from 'refine.log' file that include last iteration .
        # characteristics
        headerList, dataList, msg = self.parseFile.retrievelastIteration()
        if not dataList:
            errorWindow(self.getTkRoot(), msg)
            return

        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg=" ",
                  title= "Refmac: Last Iteration summary",
                  height=len(dataList), width=200,padding=40)

    def _visualizeRFactorPlot(self, e=None):
        """
        xplotter = Plotter(windowTitle="SAXS Curves")
        a = xplotter.createSubPlot('SAXS curves', 'Angstrongs^-1', 'log(SAXS)', yformat=False)
        a.plot(x[:, 0], numpy.log(x[:, 1]))
        a.plot(x[:, 0], numpy.log(x[:, 2]))
        if obj.experimentalSAXS.empty():
            xplotter.showLegend(['SAXS in solution', 'SAXS in vacuo'])
        else:
            xplotter.showLegend(['Experimental SAXS', 'SAXS from volume'])
        xplotter.show()
        """
        pass

    def _visualizeFOMPlot(self, e=None):
        """ Plot FOM
        """
        headerList, dataList, msg = self.parseFile.retrieveFomPlot()
        if not dataList:
            errorWindow(self.getTkRoot(), msg)
            return
        #make plotter here
        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg=" ",
                  title= "Refmac: Last Iteration summary",
                  height=len(dataList), width=200,padding=40)

    def _visualizeLLPlot(self, e=None):
        pass

    def _visualizeLLfreePlot(self, e=None):
        pass

    def _visualizeGeometryPlot(self, e=None):
        pass