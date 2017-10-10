####################################################################
#  CAMcommunicator 
#  Python class for talking to a Leica CAM server 
#  
#  A simple python library for talking to a Leica CAM server, can
#  be used interactively in ipython.
#
#
######################################################################  
#
#  NOTE: this class is required by the "LCC Module" suite for inter- 
#  facing CellProfiler with a Leica CAM server.
#
######################################################################  
#
#  AUTHOR: Volker Hilsenstein, EMBL, 
#          volker.hilsenstein at embl.de
#  Latest changes: 13. November 2014
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


"""
CAM Communicator Class.

Volker.Hilsenstein@embl.de

Talks to a Leica CAM (computer aided microscopy) over TCP/IP link.
Provides a number of convenience functions for functions that are used often.

In terms of software architecture there are a few issues, the major one is that this software is not multi-threaded.
Ideally, there would be a worker thread that constantly talks to the microscopes and places received communication on a queue
for parsing. Certain events (e.g. new image received) should then trigger callbacks. Currently, we poll the microscope whenever we require some information. Communication that doesn't relate to the event we're look
ing for gets discarded or ignored. So we rely on information from the microscope arriving in the order expected by us (but this is not guaranteed).


TODO:

* Use python logging module for diagnostic output

* Make sure each exception block displays the actual error that occured. At the moment, the error is often too unspecific as one of several commands in the block may have caused the error.

* added flag connected and corresponding isConnected() method. This doesn't check whether the connection still exists, it just relies on the open and close functions to set the flag

* make sure the methods affecting/querying stage position can work with the StagePosition class in layout.py

""" 

import socket
import time
import sys
import os
import string
import pdb
import numpy as np
try:
    import layout
    nolayoutmodule=False
except:
    nolayoutmodule=True
import re

iplocal = "127.0.0.1"
ipSP5A = "10.11.112.16" # these are convenient shorthands for internal use
ipSP5B = "10.11.112.18" 


class CAMcommunicator:
    def __init__(self):
        # Settings for TCP/IP communication
        # changed all of these from
        # class to object variables
        # which seems to make more sense
        # in case one has several instances
        # (which is unlikely in practice, though)
        self.IP_address = "127.0.0.1"
        self.sysID='1'
        self.port = 8895
        self.timeout = 120
        self.delay=0.8 # changed from default
        self.buffersize = 4096
        self.leicasocket = None
        self.verbose=True
        self.basepath = "" # Set this to the data exporter base path
        self.connected=False
        self.cmdlist = [] # CAM command list
        self.previous_file = ""
        self.previous_cmd = ""
        self.stage_settle_time = 6.5
        sequence_counter=0

    def setSysID(self, newsysID):
        self.sysID = str(newsysID)

    def getSysID(self):
        return self.sysID

    
    ###############################################
    #  network connection settings and functions
    ###############################################

    def setIP(self,ip):
        self.IP_address= string.strip(ip)

    def getIP(self):
        return self.IP_address
        



    def printSettings(self):
        print "CAMCommunicator settings"
        print "not yet implemented"

    def open(self):
        """Open connection to CAM server. Returns True if successful, False otherwise"""
        try:
            if self.verbose:
                print "Trying to open connection to Leica at ",self.IP_address,":",str(self.port)
            self.leicasocket = socket.socket()
            self.leicasocket.connect((self.IP_address,self.port))
            if self.verbose:
                print("Connected.")
            self.connected=True
            return True
        except:
            if self.verbose:
                print "Error opening connection to ", self.IP_address
            self.connected=False
            return False

    def close(self):
        """ Close connection to CAM server. Returns True if successful, False otherwise"""
        if self.verbose:
            print "Disconnecting from ", self.IP_address
        if self.leicasocket is not None:
            try:
                self.leicasocket.close()
                self.connected=False
                return True
            except:
                print "Error closing connection. Did you pass in the correct socket object ?"
                return False

    def isConnected(self):
        """ returns true if connected to a CAM server, false otherwise. This does not actually test the connectivity, it simply
        returns a flag that is modified by the open() and close() methods. If the connection has been dropped by the server or
        lost due to a network error this function will not give the proper result"""
        return self.connected
    
    def flushCAMreceivebuffer(self):
        """ reads and discards all data waiting at socket """
        self.leicasocket.setblocking(False)
        try:
            while(True):
                self.leicasocket.recv(self.buffersize)
        except:
            pass

    def FixLineEndingsForWindows(self,str):
               """Helper function to make the line ending of a string windows-compatible.
               Depending on the input string it leaves the string as is, replaces a newline with a CR+newline, or adds a CR+newline. """
               # TODO: this should not really be part of this class
               if str[-2:]=="\r\n":
                   return str
               if str[-1:]=="\n":
                   return str[:-1]+"\r\n"
               else:
                   return str + "\r\n"

    ###############################################
    #  receiving and parsing CAM notifications
    ###############################################

    def parseCAMcmd(self,cmdstr):
        '''Parses a single line received from the Leica cam server and returns a dictionary where the keys
        correspond to the part behind the slash eg: app, cli, relpath etc.'''
        tmp = cmdstr.split('/')
        cmds = [c.strip() for c in tmp if c!=''] # remove empty results and strip trailing and leading whitespaces
        result_dict = {}
        for c in cmds:
            cmdname, cmdvalue = c.split(':')
            result_dict[cmdname] = cmdvalue
        return result_dict

    def receiveCAMNotification(self):
        """ waits for a notification from the Leica Cam Server telling us that image acquisition has finished"""
        self.leicasocket.setblocking(True)
        self.leicasocket.settimeout(self.timeout)
        try:
            while True:
                fromCAMServer=self.leicasocket.recv(self.buffersize)
                # TODO check for duplicates
                lines=fromCAMServer.splitlines()
                files = []
                for line in lines:
                    parsed = self.parseCAMcmd(line)
                    if 'relpath' in parsed.keys():
                        print("new image received")
                        filename = parsed['relpath']
                        files.append(filename)
                    else:
                        print("command incomplete or doesn't contain new image")
                        return(files)
        except:
            exctype, value = sys.exc_info()[:2]
            print "exception type: ", exctype, "value: ", value
            if self.verbose:
                print("Didn't receive anything from CAM Server for " + str(self.timeout) + " seconds. Timed out.") 
            pass


    def waitforimage(self,jobnr=None, jobname=None, ignoreduplicates = True, timeout=None, stopcallback=None, processGUIEvents=None): # timeeout option ?
        """Waits for an image from the CAM server.
        Sometimes matrix screener will notify us about an image twice. If ignoreduplicates is True, such duplicates will be ignored on the second notification.
        If jobnr is not None, images which do not match the job number are ignored.
        If jobname is not None, images which do not match the job name are ignored.
        If both jobnr and jobname are not None both have to match.

        the function returns a tuple (fullfilename, metadata)

        fullfilename is the full path assembled from self.basepath and the relpath given by the CAM server. Path separators are adjusted to match the operating system,
        backslashes on Windows, forward slashes on Mac/Unix.
        metadata is a dict of metadata fields extracted from the filename (TODO- Test with CAM images)
        """

        re_pattern = "(?P<Prefix>.*)(?P<Loop>--[Ll][0-9]*)(?P<Slide>--S[0-9]*)(?P<U>--[Uu][0-9]*)(?P<V>--[Vv][0-9]*)(?P<Job>--J[0-9]*)(?P<E>--[Ee].*)(?P<O>--O.*)(?P<X>--[Xx][0-9]*)(?P<Y>--[Yy][0-9]*)(?P<T>--[Tt][0-9]*)(?P<Zpos>--[Zz][0-9]*)(?P<Channel>--[Cc][0-9]*)(?P<Suffix>.*)(\.ome.tif$)"
        # on some systems, the file name includes an additional --M field, in this case we need to use a different regular
        # expressen
        re_pattern_m = """(?P<Prefix>.*)(?P<Loop>--[Ll][0-9]*)(?P<Slide>--S[0-9]*)(?P<M>--[Mm][0-9]*)(?P<U>--[Uu][0-9]*)(?P<V>--[Vv][0-9]*)(?P<Job>--J[0-9]*)(?P<E>--[Ee].*)(?P<O>--O.*)(?P<X>--[Xx][0-9]*)(?P<Y>--[Yy][0-9]*)(?P<T>--[Tt][0-9]*)(?P<Zpos>--[Zz][0-9]*)(?P<Channel>--[Cc][0-9]*)(?P<Suffix>.*)\.ome.tif$"""

        try:
            starttime = time.time()
            while True:
                if timeout is not None:
                    if time.time()-starttime > timeout:
                        return None
                if  processGUIEvents is not None:
                    processGUIEvents()
                if stopcallback is not None:
                    if not stopcallback():
                        return None
                    

                msgs=self.readandparseCAM(stopcallback=stopcallback, processGUIEvents=processGUIEvents)
                # when timed out msgs will be None. Need to make sure this is not the case
                if msgs is not None: 
                    metadata = {}
                    if jobnr is not None:
                        jobstr = "J" + str(jobnr).zfill(2)
                    else:
                        jobstr = ""

                    # go through messages in reverse order
                    #
                    # WARNING: if several images were received this function will only return the last one as we analyze the different
                    # messages from the CAM server in reverse order.
                    #
                    # For our purposes this is good enough but be aware of the limitation

                    # TODO sometimes msgs appear to be None

                    msgs.reverse()
                    for m in msgs:
                        if 'relpath' in m.keys():
                            #print m['relpath']
                            #fname = self.basepath + os.sep + m['relpath'].replace("\\",os.sep)
                            fname = self.basepath + m['relpath'].replace("\\",os.sep)
                            if ignoreduplicates is False or fname != self.previous_file:
                                self.previous_file = fname # workaround as some files are reported twice
                                print "New file " + fname
                                # TODO: check whether filename contains the substring CAM and use different re pattern if necessary
                                if ("--m" in fname) or ("--M" in fname):
                                    pattern = re_pattern_m
                                    withM = True
                                    print("filename contains --M")
                                else:
                                    pattern = re_pattern
                                    withM = False
                                    print("regular filename")
                                re_m = re.match(pattern, fname)
                                if re_m is not None:
                                    metadata['prefix'] = (re_m.group('Prefix'))
                                    metadata['loop'] = (re_m.group('Loop'))
                                    metadata['slide'] = (re_m.group('Slide'))
                                    if withM:
                                        metadata['M'] = (re_m.group('M'))
                                    else:
                                        metadata['M'] = ''
                                    metadata['U'] = (re_m.group('U'))
                                    metadata['V'] = (re_m.group('V'))
                                    metadata['job'] = (re_m.group('Job'))
                                    metadata['E'] = (re_m.group('E'))
                                    metadata['other'] = (re_m.group('O'))
                                    metadata['X'] = (re_m.group('X'))
                                    metadata['Y'] = (re_m.group('Y'))
                                    metadata['tpoint'] = (re_m.group('T'))
                                    metadata['zpos'] = (re_m.group('Zpos'))
                                    metadata['channel'] = (re_m.group('Channel'))
                                    metadata['suffix'] = (re_m.group('Suffix'))
                                else:
                                    print "Error  extracting metadata from filename ", fname
                                    raise

                                if jobname is None or ('jobname' in m.keys() and m['jobname'].lower()==jobname.lower()):
                                    if jobstr in metadata['job']:
                                        # our job matches all criteria
                                        print "file matches selection criteria. breaking out of loop"
                                        return (fname, metadata)
                                    else:
                                        print "job number is different from requested"
                                else:
                                    print "job name does not match required name"
                            else:
                                print "Ignoring Duplicate! Cam server reported file twice."
        except:
            print "Unexpected error:", sys.exc_info()
            # TODO ... make sure we actually catch the correct exception
            return None
            

                                                
    def readandparseCAM(self, stopcallback=None, processGUIEvents=None):
        """ reads pending (or waits for incoming until timeout) CAM responses from the server.
        The responses are parsed into python dictionaries and a list with the parsed responses (each list entry is a dictionary) is returned. """
        assert self.leicasocket is not None
        self.leicasocket.setblocking(True)
        self.leicasocket.settimeout(self.timeout)
        try:
            while True:
                if stopcallback is not None:
                    if not stopcallback():
                        return None
                if  processGUIEvents is not None:
                    processGUIEvents()
                fromCAMServer=self.leicasocket.recv(self.buffersize)
                lines=fromCAMServer.splitlines()
                responses = []
                for line in lines:
                    parsed = self.parseCAMcmd(line)
                    responses.append(parsed)
                return(responses)
        except:
            print("Didn't receive anything from CAM Server for " + str(self.timeout) + " seconds. Timed out.") 
            return None


    #############################################################################
    #  Scanfield enable/disable/assign job, 
    #############################################################################

    def enableScanField(self, allfields=False, wellx=1, welly=1, fieldx=1, fieldy=1, value=True, slide=1):
        '''Enables either the scanfield at the specified position or all scanfields if allfields is True'''
        c = "/cli:python /app:matrix "
        print ["disabling","enabling"][value], " fields"
        if allfields is False:
            c += "/cmd:enable /slide:" + str(slide) + " /wellx:"+str(wellx) + " /welly:" + str(welly)
            c += " /fieldx:"+str(fieldx) + " /fieldy:" + str(fieldy)
        else:
            print "allfields"
            c +=  "/cmd:enableall"

        c += " /value:" + ("false","true")[value]
        self.sendCMDstring(c)
        time.sleep(0.2)

    def enableScanFields(self, fields, wellx=1, welly=1, value=True, slide=1):
        '''Enables each scanfield provided in the list of tuples "fields" in well wellx, welly'''
        print ["disabling","enabling"][value], " fields", fields
        basec = "/cli:python /app:matrix "
        basec += "/cmd:enable /slide:" + str(slide) + " /wellx:"+str(wellx) + " /welly:" + str(welly)
        for fieldx, fieldy in fields:
            c = basec + " /fieldx:"+str(fieldx) + " /fieldy:" + str(fieldy)
            c += " /value:" + ("false","true")[value]
            self.sendCMDstring(c)
            time.sleep(0.2)

    def disableScanField(self, allfields=False, wellx=1, welly=1, fieldx=1, fieldy=1, slide=1):
        '''disableScanField is just a convenience wrapper that calls enableScanField mit value=False'''
        self.enableScanField(allfields, wellx, welly, fieldx, fieldy, value=False, slide=1)

    def disableScanFields(self, fields, wellx=1, welly=1, slide=1):
        '''disableScanFields is just a convenience wrapper that calls enableScanField mit'''
        self.enableScanField( fields, wellx, welly, value=False, slide=1)

    def selectScanField(self, allfields=False, wellx=1, welly=1, fieldx=1, fieldy=1):
        """Select scan field in wellx, welly, fieldx, fieldy. If allfields==True, all scan fields are selected"""
        # TODO ... allow to pass lists with fields
      
        if allfields is False:
            c = "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:selectfield /wellx:"+str(wellx) + "/welly:"+str(welly) + " /fieldx:"+str(fieldx) + " /fieldy:" + str(fieldy)
            print "selecting scanfield ", wellx, " ", welly, " ", fieldx, " ", fieldy
        else:
            print "selecting all scanfields"
            c =  "/cli:python /app:matrix /cmd:selectallfields"
        self.sendCMDstring(c)
        time.sleep(0.5)

    def assignJob(self, jobname):
        """assign a Job to the the currently selected positions"""
        c = "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:assignjob /job:"+jobname.lower()  # convert jobname to lowercase as workaround
        print "Assigning ", jobname
        self.sendCMDstring(c)
        time.sleep(0.5)

    def assignJobToScanFieldDisableAllOthers(self, jobname, wellx=1, welly=1, fieldx=1, fieldy=1):
        """Convenience function for assigning a job to a scanfield, and disabling all other scanfields"""
        self.selectScanField(False, wellx, welly, fieldx, fieldy)
        self.assignJob(jobname)
        self.disableScanField(allfields=True)
        self.enableScanField(False, wellx, welly, fieldx, fieldy)

    def assignJobToScanFields(self, jobname, fields,  wellx=1,  welly=1):
        """Convenience function for assigning a job to a scanfield, and disabling all other scanfields"""
        time.sleep(0.2)
        first=True
        for fieldx, fieldy in fields:
            self.selectScanField(False, wellx, welly, fieldx, fieldy)
            time.sleep(0.2)
            self.assignJob(jobname)
            if first:
                # need an extra long wait after assigning the first job.
                # this is because Matrix Screener somehow takes that long to switch
                # between jobs and as we don't know which job was previously selected we want
                # to be on the save side
                time.sleep(5)
                first=False
            else:
                time.sleep(2)
        time.sleep(0.2)
        #self.disableScanField(allfields=True)
        #for fieldx, fieldy in tmpfields:
         #   self.enableScanField(False, wellx, welly, fieldx, fieldy)


    #############################################################################
    #  Drift compensation
    #############################################################################

    def enableDriftCompensation(self, allfields=False, wellx=1, welly=1, fieldx=1, fieldy=1, value=True, slide=1):
        '''Enables the scanfield at the specified position or all scanfields if allfields is True'''
        c = "/cli:python /app:matrix /cmd:enableattribute "

        if allfields is False:
            c += " /slide:" + str(slide) + " /wellx:"+str(wellx) + " /welly:" + str(welly)
            c += " /fieldx:"+str(fieldx) + " /fieldy:" + str(fieldy)

        c += " /drift:" + ("false","true")[value]
        c += " /track:false /pump:false" 
        self.sendCMDstring(c)

    def disableDriftCompensation(self, allfields=False, wellx=1, welly=1, fieldx=1, fieldy=1, slide=1):
        '''Diables the scanfield at the specified position or all scanfields if allfields is True'''
        self.enableDriftCompensation(allfields, wellx, welly, fieldx, fieldy, value=False, slide=slide)
    ####################################################
    #  Starting,stopping, dealing with CAMlist
    #############################################################################

    # normal experiment
        
    def startScan(self):
        self.sendCMDstring("/cli:python /app:matrix /cmd:startscan")

    def pauseScan(self):
        self.sendCMDstring("/cli:python /app:matrix /cmd:pausescan")

    def stopScan(self):
        self.sendCMDstring("/cli:python /app:matrix /cmd:stopscan")

    # AF Scan
    def startAFScan(self):
        self.sendCMDstring("/cli:python /app:matrix /cmd:autofocusscan")

    # CAM list 
    def addJobToCAMlist(self, jobname, dxpos, dypos, slide=1, wellx=1, welly=1, fieldx=1, fieldy=1, ext=None):
        #TODO
        prefix = "/cli:python /app:matrix /cmd:add"
        tar = "/tar:camlist"
        exp = "/exp:" + jobname
        s = "/slide:" + `int(slide)`
        wx = "/wellx:" + `int(wellx)`
        wy = "/welly:" + `int(welly)`
        fx = "/fieldx:" + `int(fieldx)`
        fy = "/fieldy:" + `int(fieldy)`
        dx = "/dxpos:" + `int(dxpos)` 
        dy = "/dypos:" + `int(dypos)`
        if ext is None:
            e= ""
        else:
            e="/ext:"+ext
        cmd = " ".join((prefix,tar,exp,s,wx,wy,fx,fy,dx,dy,e))
        self.sendCMDstring(cmd)
    
    def deleteCAMList(self):
        self.sendCMDstring("/cli:python /app:matrix /cmd:deletelist")

    def startCAMScan(self, runtime=None, repeattime=None, afinterval=None, trackinterval=None, pumpinterval=None):
        # TODO add /runtime /repeattime /afinterval /trackinterval /pumpinterval options
        self.sendCMDstring("/cli:python /app:matrix /cmd:startcamscan")

    def stopCAMScan(self):
        self.sendCMDstring("/cli:python /app:matrix /cmd:stopcamscan")

    def stopWaitingForCAM(self):
        self.sendCMDstring("/cli:python /app:matrix /cmd:stopwaitingforcam")

    def waitForScanToFinish(self):
        """"loop indefinitely until we receive scanfinished"""
        while True:
            answers = self.readandparseCAM()
            if answers is not None:
                for a in answers:
                    if 'inf' in a.keys():
                        if a['inf'] == "scanfinished":
                            return
    
    ###############################################
    # CMD list - handling
    # TODO: explain concept
    ###############################################

    def emptyCMDlist(self):
        self.cmdlist=[]

    def getCMDlist(self):
        return self.cmdlist

    def addtoCMDlist(self,cmd):
        self.cmdlist.append(cmd)

    def augmentCMDlist(self):
        # these should be deprecated
        self.cmdlist.insert(0,"/cli:python /app:matrix /cmd:deletelist")
        self.addtoCMDlist("/cli:python /app:matrix /cmd:startcamscan")

    def saveCMDlist(self, filename):
        """Given a list of CAM commands in cmdlist, save them into file with filename. """
        f = open(filename,'w')
        for cmd in self.cmdlist:
            try:
                tmp = self.FixLineEndingsForWindows(cmd)
                f.write(tmp)
            except:
                return False
        f.close()
        return True

    def sendCMDlist(self):
        """ This function sends each string in cmdlist to the CAMserver.
        A delay between successive commands can be specified  in self.delay (default is 0.2s).
        Line endings are fixed to be Windows-compatible, i.e. CR+LF.
        After successful completion the list is emptied""" 

        if self.cmdlist:
            for cmd in self.cmdlist:
                try:
                    tmp = self.FixLineEndingsForWindows(cmd)
                    charssent= self.leicasocket.send(tmp)
                    # we actually need to make sure
                    # we sent the whole string by comparing charssent.
                    if charssent != len(tmp):
                        print "Error sending commands"
                        raise CAMSendCharsError
                except:
                    print "error sending command", cmd
                    return False
                time.sleep(self.delay) # wait some time between sending each line
            self.emptyCMDlist()
            time.sleep(self.delay)
            
    def sendCMDstring(self, cmdstr, seq_counter=False):
        """Sends cmdstr to the leica. Internally the CMDlist is emptied, the string is added and the list is cleared"""
        
        self.flushCAMreceivebuffer() # TODO ... can we really flush here ?
        self.emptyCMDlist()
        self.addtoCMDlist(cmdstr)

        # The following lines are a workaround around a bug in CAM server.
        # Sometimes (not always) the CAM server echos received commands twice.
        # In some of the client applications we filter out duplicate commands, as that often makes sense,
        # e.g. we can assume that the same image isn't saved twice. However, there are situations, such
        # as when using CAM server as a communication interface with external hardware, where a command
        # is genuinely sent twice in immediate succession. We can't distinguish between erroneous duplicate 
        # echos and genuine commands being sent twice. Therefore the option seq_counter was added. If set to True,
        # we send a second command with a unique sequence number to disambiguate between these two cases.

        # TODO test this
        if seq_counter:
            self.addtoCMDlist("/cli:python /app:external /name:disambiguation /cmd:seq /par1:"+str(self.sequence_counter))
            self.sequence_counter +=1 
        
        self.sendCMDlist()


    ###############################################
    #  load / save scanning templates
    ###############################################

    #load

    #save

    ###############################################
    #  stage positioning
    ###############################################

    def setStagePositionSafely(self, pos, high_stage_position=-0.00075):
        """ This function will first raise the stage, to make more space to clear the objective, then
        travel to commanded X,Y Position and then move the stage to the commanded Z position. """
        # TODO: disabled moving stage up and down. Test and revisit.
        #

        if len(pos) == 2: # if no Z-position is given, we remember the current Z position
            finalZ= self.getCurrentStagePosition()[2]
        else:
            finalZ= pos[2]

        # First, raise the stage to the maximum setting (NEED TO ACTIVATE Hi-Z OPTION IN LASAF CONFIGURATION !!!)
        # This will give us maximum clearance above the objective
        #self.setStageZPosition(high_stage_position)
        # Now we move to the target X,Y position ...
        self.setStageXYPosition(pos[0:2])
        # ... and wait a while so we can be sure the stage has arrived
        time.sleep(self.stage_settle_time)
        # finally we lower the stage to the commanded Z position (or the Z position before the move, if no
        # Z position was specified)
        #if not np.isnan(finalZ) and finalZ is not None:
        #    self.setStageZPosition(finalZ)
        #time.sleep(self.stage_settle_time)
        

    def setStageXYPosition(self,pos, relative=False):
        c= "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:setposition /typ:" + ("absolute","relative")[relative] +" /dev:stage /unit:meter /xpos:" + "{:.12f}".format(pos[0]) + " /ypos:" + "{:.12f}".format(pos[1])
        self.sendCMDstring(c)

    def setStagePosition(self,pos, relative=False):
        """Sets the stage position. Sets either X,Y only if input has less than 3 coordinates or third coordinate is nan,
        or X,Y,Z otherwise"""
        if len(pos) == 3:
            if not np.isnan(pos[2]):
                self.setStageZPosition(pos[2], relative)
        self.setStageXYPosition(pos[0:2],relative)

    def setStageZPosition(self,z, relative=False):
        c= "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:setposition /typ:" + ("absolute","relative")[relative] + " /dev:zdrive /unit:meter /zpos:" + "{:.12f}".format(z)
        self.sendCMDstring(c)
        
    def getCurrentStagePosition(self):
        """Queries the CAM server for the current stage position and returns that position as a 3-tuple of float values"""
        
        # send query
        self.sendCMDstring("/cli:python /app:matrix /sys:"+self.sysID+" /cmd:getinfo /dev:stage")
        # wait for and parse response
        resp=self.readandparseCAM()[0]
        if resp['dev']=='stage':
            if nolayoutmodule:
                sp = (float(resp['xpos']),float(resp['ypos']),float(resp['zpos']))
            else:
                sp = layout.StagePosition()
                sp[:] = (float(resp['xpos']),float(resp['ypos']),float(resp['zpos']))
            return sp
        else:
            return None

    def saveCurrentStagePosition(self):
        """Saves the current stage position (X,Y) - use returnToStagePosition to revisit"""
        c= "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:savecurrentposition /dev:stage"
        self.sendCMDstring(c)

    def returnToStagePosition(self):
        """returns to saved stage position (X,Y) - use saveCurrentStagePosition to set """ 
        c= "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:returntosavedposition /dev:stage"
        self.sendCMDstring(c)

    def saveCurrentZPosition(self):
        """Saves the current zdrive position (Z) - use returnToZPosition to revisit"""
        c= "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:savecurrentposition /dev:zdrive"
        self.sendCMDstring(c)

    def returnToZPosition(self):
        """returns to saved Z position (Z) - use saveCurrentZPosition to set """ 
        c= "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:returntosavedposition /dev:zdrive"
        self.sendCMDstring(c)
        
    def moveStageToWell(self, u, v):
        """moves stage to well u,v """ 
        c= "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:moveteowell "
        c += " /upos:"+str(u)
        c += " /vpos:"+str(v)
        self.sendCMDstring(c)
        
    ##################################################
    #  Scan job settings (laser, PMT, scan pattern (not implemented))
    ###############################################

    def setLaserPower(self,wavelength, power, job, seq=1):
        """Set the laser power for a job. Laser is identified by wavelength (int or str)"""
        self.emptyCMDlist()

        c = "/cli:python /app:matrix /sys:"+self.sysID+" /cmd:adjustls /seq:" + str(seq) + " /exp:" + job
        if int(wavelength) == 405:
            c += " /lsq:uv /lid:405"
        else:
            c += " /lsq:vis /lid:" + str(wavelength)
        c+= " /tar:laser /value:" + str(float(power))
        self.addtoCMDlist(c)
        self.sendCMDlist()

    def setPMT(self,job, n, gain=None, offset=None):
        """Set the photomultiplier parameters for PMT number n and the specified job. Both gain and offset can be specified individually or together"""
        if gain==None and offset==None:
            return
        self.emptyCMDlist()
        c = "/cli:python /app:matrix /cmd:adjust /tar:pmt /exp:" + job +" /num:"+ str(int(n))
        # Set gain
        if gain is not None:
            tmp= c + " /prop:gain /value:" + str(float(gain))
            self.addtoCMDlist(tmp)

        # Set offset
        if offset is not None:
            tmp= c + " /prop:offset /value:" + str(float(offset))
            self.addtoCMDlist(tmp)

        # send to CAM server
        self.sendCMDlist()

    def setPinhole(self, job, value):
        """Set the pinhole in job <job> to <value> ... Not sure about units"""
        c ="/cli:python /app:matrix /cmd:adjust /tar:pinhole /exp:"+ job + " /value:" +str(value)
        self.addtoCMDlist(c)
        self.sendCMDlist()

    def setPump(self, time, wait):
        """Not tested, presumably sets the time that the water immersion pump is active <time> and the time to wait after pumping."""
        c ="/cli:python /app:matrix /cmd:pump /time:"+str(time)+ " /value:"+ str(value)
        self.addtoCMDlist(c)
        self.sendCMDlist()

    ###############################################
    #  Mosaic
    ###############################################

    # there are a whole lot of mosaic options that we have never tested

    ###############################################
    #  FCS functions
    #  legacy FCS functions.
    #  These commands are probably only supported
    #  on Malte's machine and not on any other SP5.
    #  This special matrix
    #  screener version pauses after every prescan
    #  image until it receives a FCS startscan command.
    #
    #  The whole procedure is very similar to CAM list,
    #  however, the job is only specified in pixel coordinates
    #  for the previous prescan image. There is no
    #  way to send well or field coordinates.
    #  currently zpos is always zero
    ###############################################

    def addFCSPoint(self, dxpos, dypos, dzpos=None):
        """adds an FCS measurement point at coordinates dxpos, dypos, dzpos to the list"""
        prefix = "/cli:python /app:fcs /cmd:add"
        dx = "/xpos:" + `int(dxpos)` 
        dy = "/ypos:" + `int(dypos)`
        if dzpos is None:
            dz = "/zpos:0"
        else:
            dz = "/zpos:"+dzpos
    
        
        cmd = " ".join((prefix,dx,dy,dz))
        #print "FCS command string"
        #print cmd
        self.sendCMDstring(cmd)

    def removeFCSPoints(self):
        """clears the list of FCS measurement point"""
        c = "/cli:python /app:fcs /cmd:removeall"
        self.sendCMDstring(c)

    def startFCSscan(self):
        """starts FCS measurements for the points that have been added to the list"""
        c = "/cli:python /app:fcs /cmd:startscan"
        self.sendCMDstring(c)

    def stopFCSscan(self):
        """no documentation. presumably stops an ongoing FCS scan"""
        c = "/cli:python /app:fcs /cmd:stopscan"
        self.sendCMDstring(c)
    
    ###############################################
    #  external devices (talk to immersion oil
    #  dispenser, drug dispenser, laser etc.)
    #  via the /app:external command
    #  need one program running on the Leica side
    #  to listen to the communication
    ###############################################


    ### liquid dispenser
    ld_base_cmd = "/cli:python /app:external /name:ld "  
    # client side commands
    def ld_home(self):
        c = self.ld_base_cmd + " /cmd:home"
        self.sendCMDstring(c, True)
        
    def ld_selectSolvent(self, reservoir_nr):
        c = self.ld_base_cmd + " /cmd:selectSolvent /par1:" + str(reservoir_nr)
        # TODO check whether number is in range
        self.sendCMDstring(c, True)

    def ld_dispense(self, seconds):
        c = self.ld_base_cmd + " /cmd:dispense " + "/par1:" + str(seconds)
        self.sendCMDstring(c, True)

    def ld_purge(self, seconds):
        c = self.ld_base_cmd + " /cmd:purge " + "/par1:" + str(seconds)
        self.sendCMDstring(c, True)

    def ld_wash(self, seconds):
        c = self.ld_base_cmd + " /cmd:wash " + "/par1:" + str(seconds)
        self.sendCMDstring(c, True)

    # server side-commands
    def ld_response(self, success, text):
        c = self.ld_base_cmd + " /cmd:response " + "/par1:" + ("false", "true")[success] + " /par2:" + text 

    ### immersion oil dispenser
    od_base_cmd = "/cli:python /app:external /name:od "  
    # client side commands
    def od_pump_on(self):
        c = self.od_base_cmd + " /cmd:pumpon"
        self.sendCMDstring(c, True)
        
    def od_pump_off(self): 
        c = self.od_base_cmd + " /cmd:pumpoff"
        self.sendCMDstring(c, True)
        
    def od_pump_speed(self, value):
        c = self.od_base_cmd + " /cmd:pumpspeed " + "/par1:" + str(value)
        self.sendCMDstring(c, True)
        
    def od_pump_for_x_seconds(self, seconds):
        c = self.od_base_cmd + " /cmd:pumpxseconds " + "/par1:" + str(seconds)
        self.sendCMDstring(c, True)

    # server side commands        
    def od_response(self, success, text):
        c = self.od_base_cmd + " /cmd:response " + "/par1:" + ("false", "true")[success] + " /par2:" + text 

    ### bleaching laser
    bl_base_cmd = "/cli:python /app:external /name:bl "  
    # client side commands
    def bl_on(self):
        c = self.bl_base_cmd + " /cmd:on"
        self.sendCMDstring(c, True)

    def bl_off(self):
        c = self.bl_base_cmd + " /cmd:off"
        self.sendCMDstring(c, True)

    def bl_power_level(self, value):
        c = self.bl_base_cmd + " /cmd:power /par1:" + str(value)
        self.sendCMDstring(c, True)

    # server side commands
    def bl_response(self, success, text):
        c = self.bl_base_cmd + " /cmd:response " + "/par1:" + ("false", "true")[success] + " /par2:" + text 

    ###############################################
    #  Utility/various
    ###############################################

    def loop(self,runhours, runmin, runsec, repmin, repsec):
        """not sure what this one does, presumably adjusts time lapse parameters ... untested """
        c ="/cli:python /app:matrix /cmd:loop /runh:"+str(runhours)+ " /runm:"+str(runmin)
        c += " /runs:"+str(runsec) + " /repm:"+ str(repmin) + " /reps:"+str(reps)
        self.addtoCMDlist(c)
        self.sendCMDlist()
        
    def createAbsPath(self,filename):
        """ given a base path and a filename received from Leica CAM server, create the absolute filepath on the local machine. Changes path delimeters from \\ to / depending on operating system"""
        if "Subfolder" in self.basepath:
            print "Warning !!!\nYou provided baspath: "+self.basepath +"\nThis includes /Subfolder/. You probably need to specify the path without Subfolder." 
        return self.basepath + os.sep + filename.replace("\\", os.sep)

        # TODO: currently this returns only files (relpath), everything else is ignored although messages are fully parsed
        # change such that all parsed commands are returned
        # filtering for "relpaths" or other commands should happen outside

    def corrRing(self, angle):
        """this should set the correction ring"""
        c ="/cli:python /app:matrix /cmd:setcorring" + str(float(angle))
        self.sendCMDstring(c)

    def getJobDict(self):
        """gets the list of defined jobs from matrix screener"""
        c = "/cli:python /app:matrix /cmd:getinfo /dev:joblist"
        self.sendCMDstring(c)
        time.sleep(self.delay)
        answers = self.readandparseCAM()
        joblist = {}
        for a in answers:
            if a['dev']=='joblist':
                for i in range(int(a['count'])):
                    nr = a['jobid' +str(i+1)]
                    name = a['jobname' +str(i+1)].lower()
                    joblist[name]=nr
            else:
                print "no joblist in answers"
        return joblist

    def getPatternDict(self):
        """gets the list of defined patterns from matrix screener"""
        c = "/cli:python /app:matrix /cmd:getinfo /dev:patternlist"
        self.sendCMDstring(c)
        time.sleep(self.delay)
        answers = self.readandparseCAM()
        patternlist = {}
        for a in answers:
            if a['dev']=='patternlist':
                for i in range(int(a['count'])):
                    nr = a['patternid' +str(i+1)]
                    name = a['patternname' +str(i+1)].lower()
                    patternlist[name]=nr
            else:
                print "no patternlist in answers"
        return patternlist
