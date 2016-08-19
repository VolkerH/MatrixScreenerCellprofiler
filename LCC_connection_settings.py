####################################################################
#  LCConnect
#  Leica CAM Connect Module
#  
#  handles/stores communication settings between Cellprofiler and 
#  Leica CAM server for use in the LCCWaitForImage and 
#  LCCImageObjectWithMicroscope modules
#
######################################################################  
#  requires cam_communicator_class.py
#
######################################################################  
#
#  AUTHOR: Volker Hilsenstein, EMBL, 
#          volker.hilsenstein at embl.de
#  Latest changes: 1. September 2014
#
#
#  LICENSE:
#  This software is distributed under an EMBLEM academic use software 
#  license, see the file LCC_license.txt for details.   
#
#  Note that the "LCC Module"-suite is developed at EMBL and all 
#  support requests relating to these modules should be related to the
#  author, not to the CellProfiler team.
#
######################################################################  
#
# If you want to acknowledge the use of these modules in your research 
# publications please cite:
#
# C.Tischer, V.Hilsenstein, K.Hanson, R.Pepperkok
# "Adaptive Fluorescence Microscopy by Online Feedback Image Analysis", 
# In "Quantitative Imaging in Cell Biology"
# Methods Cell Biol. 2014;123:489-503. 
# doi: 10.1016/B978-0-12-420138-5.00026-4.
# 
######################################################################  


import cellprofiler.cpmodule as cpm
import cellprofiler.pipeline as cpp
import cellprofiler.settings as cps
import cellprofiler.cpimage as cpi

import cam_communicator_class as cc

global CAMC 
CAMC = cc.CAMcommunicator()

IP_ADDRESS_TEXT = "IP Address of CAM server"
BASEPATH_TEXT = """Path to images (on machine running CellProfiler, excluding "Subfolder")"""

class LCConnect(cpm.CPModule):
    
    ################### Name ##########################
    variable_revision_number = 1
    module_name = "LCConnect"
    category = "MicroscopeAutomation"


    CAMC = cc.CAMcommunicator()

    ################# GUI Settings ########################
    def create_settings(self):

        self.IP_address = cps.Text(IP_ADDRESS_TEXT, "127.0.0.1",
                                   metadata = False,
                                   doc="""
                                   IP address of the leica cam server""")

        self.basepath = cps.Text(BASEPATH_TEXT, "Z:\\MatrixExportFolder",
                                 metadata = False,
                                 doc="""
                                 Basepath to images on the computer running cellprofiler. Use the path seperators for the operating system that Cellprofiler is running on.""")

        self.sysID = cps.Integer("Leica /sys value (typically 0)", value = 0, minval = 0, doc = """some Matrix screener CAM commands require passing in a system identifier. On most microscopes I have seen this ID is zero, but in some rare cases you may have to use a value of 1 (or something else)""")
        
        
    def do_connect(self):
        print "Connecting"
        CAMC.setIP(self.IP_address.value)
        CAMC.setSysID(self.sysID.value)
        CAMC.open()

    def do_disconnect(self):
        print "Disconnecting"
        CAMC.close()

    def do_check_status(self):
        # TODO: this connection status doesn't actually check whether the connection is still alive.
        # maybe remove as it isn't terribly useful in its current state
        print "Socket is", ("disconnected","connected")[CAMC.isConnected()]

    def settings(self):
        return [ self.IP_address,  self.basepath, self.sysID] 
    
    def visible_settings(self):
        return self.settings()
    
    def getCAMCommunicator(self):
        return CAMC

        
    def run(self, workspace):
        print "setting basepath to ", self.basepath.value
        CAMC.basepath=self.basepath.value
        print "setting IP address ", self.IP_address.value
        CAMC.setIP(self.IP_address.value)
        CAMC.setSysID(self.sysID.value)
