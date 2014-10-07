####################################################################
#  LCCimageObject
#  Leica CAM Connect image object module
#  
#  adds object coordinates of CellProfiler objects to the Leica 
#  Matrix Screener CAM list and invokes a CAM scan.
#  Several options provide the flexibility to correct for 
#  flipped stage coordinate systems, add fixed offsets for compensating
#  misalignments between different objectives. Multiple LCCimage 
#  modules can be used in the same pipeline, e.g. to invoke different
#  imaging jobs for different object classes. 
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

'''<b>...LCCimageObject</b> - start a high-resolution imaging job for the given object. 
<hr>
This module is  used to invoke a scanning job or pattern on the Leica SP5/SP8 for the given object set by sending the object coordinates to the Leica CAM server.

'''

# TODO
#
# reset object counter in prepare_run !

#################################
#
# Imports from CellProfiler
#
##################################
#import pdb
import cellprofiler.cpimage as cpi
import cellprofiler.cpmodule as cpm
import cellprofiler.measurements as cpmeas
import cellprofiler.objects as cpo
import cellprofiler.settings as cps

from cellprofiler.modules.identify import M_LOCATION_CENTER_X, M_LOCATION_CENTER_Y
import cam_communicator_class as cc


from LCC_connection_settings import CAMC


##################################
# Constants
###################################

# TODO - move into text file and import to avoid boilerplate code
MIC_WAITFORIMAGE = "Wait for image"
MIC_DONOTHING = "Do nothing"

CHOICE_STARTCAMJOB_DEFAULT = "Send start CAM job"
CHOICE_STARTCAMJOB_NONE = "Do net send start CAM job"
CHOICE_STARTCAMJOB_POSTRUN = "Send start CAM job - postrun"

CHOICE_STOPWAITING_DEFAULT = "Send 'stopwait' command"
CHOICE_STOPWAITING_NONE = "Don't send 'stopwait' command"
CHOICE_STOPWAITING_POSTRUN = "Send 'stopwait' command - postrun"

CENTER_LOCATION = "Location_X/Y (default)"
CENTER_AREASHAPECENTER = "AreaShape_Center_X/Y (if available)"
CENTER_BB = "Bounding Box Center (if available)"

'''This is the measurement template category'''
C_MEASUREMENT_TEMPLATE = "MT"

#########
# Globals
#########

# this is a dictionary that will keep track how many objects have been imaged in a particular well.
# indexing is via a string of time points, well coordinates and object name, e.g. "T0000,U04,V06,nuclei"
#nr_objs_in_well={}


###################################
## The actual ImageObjectWithMicroscope module
####################################

class LCCimageObject(cpm.CPModule):

    ################### Name ##########################
    module_name = "LCCimageObject"

    category = "MicroscopeAutomation"
    variable_revision_number = 1
    
    ################# GUI Settings ########################
    def create_settings(self):

        self.input_object_name = cps.ObjectNameSubscriber(
            "Input object name",
            doc = """Pick the objects you want to image.""")


        self.CAMJob = cps.Text("Name of CAM job", "Hiresjob", metadata = False,
                               doc="""Name of the CAM job to be called for each object""")

        self.advancedOptions = cps.Binary("Advanced options ", False, doc = """set advanced options""")

        self.maxNrObjsPerWell = cps.Integer(
            "Maximum nr of objects to image per well (-1 is no limit)",
            -1, doc="""This will allow you to limit the number of hires scans for this object for each well. If the value is -1 there is no limit""")


        self.stopWaitingForCAM = cps.Choice(
            "send StopWaitingForCAM command ?",
            [CHOICE_STOPWAITING_DEFAULT,CHOICE_STOPWAITING_NONE, CHOICE_STOPWAITING_POSTRUN], doc="Choose if and when to send the stop waiting for CAM job command to the CAM server.")

        self.startCAMJob = cps.Choice(
            "Send startcamjob command ?",
            [CHOICE_STARTCAMJOB_DEFAULT,
             CHOICE_STARTCAMJOB_NONE,
             CHOICE_STARTCAMJOB_POSTRUN], doc="""The microscope won't start working through the list of CAM jobs (e.g. hires scans) that this module is sending
            until it receives the startcamjob command. You can choose whether you want to <ul>
            <li>send this command (default),
            <li>not send it (because you may have additional instances of this module, in which case you probably only want to tick this
            option in the last instance), or
            <li>send it after all images have been analysed (post-run). The post-run option is useful if you are processing previously
            captured images offline, i.e. without a LeicaWaitForImage module.
            </ul>""")
        
        
        self.centerMeasurement = cps.Choice(
            "Which measurement to use for object center ?",
            [CENTER_LOCATION, CENTER_AREASHAPECENTER, CENTER_BB], doc="Choose which measurement to use for the object location.")


        doc_offset = """Sometimes it is required to specify a pixel offset, e.g. when the objective lens is changed for the highres job and there is some lateral displacement. Normally you want to leave this at zero."""


        doc_flip='''You can flip the direction (sign) of the coordinate axes here, this might be necessary depending on the orientation of your microscope stage'''
        
        self.deleteCAMList = cps.Binary(
            "delete existing CAM list",
            True, doc=doc_flip)

        #####
        #  The following options could potentially be moved to the CAM connect module.
        #  Flipx,flipy, swapxy shouldn't be different for different objects
        #  offsetx,offsety may be different though.
        #

        
        self.flipx = cps.Binary(
            "Invert direction of x-coordinate axis",
            False, doc=doc_flip)

        self.flipy = cps.Binary(
            "Invert direction of y-coordinate axis",
            False, doc=doc_flip)
        
        self.swapxy = cps.Binary(
            "Swap the x and y coordinates",
            False, doc="Swap the x and y coordinates.")

        self.offsetX = cps.Integer(
            "Pixel offset along X axis",
            0, doc=doc_offset)

        self.offsetY = cps.Integer(
            "Pixel offset along Y axis",
            0, doc=doc_offset)
        
    #
    def settings(self):
        self.base_settings = [self.input_object_name, self.CAMJob,self.advancedOptions]
        self.advanced_settings = [ self.stopWaitingForCAM,   self.startCAMJob, self.deleteCAMList, self.maxNrObjsPerWell, self.offsetX, self.offsetY,self.flipx, self.flipy, self.swapxy, self.centerMeasurement ]
        return  self.base_settings+self.advanced_settings

    def visible_settings(self):
        self.settings() # need to call, otherwise self.base_settings can be None
        if self.advancedOptions.value is False:
            return self.base_settings
        else:
            return self.base_settings+self.advanced_settings
    
    def prepare_run(self, pipeline, image_set_list, frame):
        self.nr_objs_in_well = {}
        return True
    # Main
    def run(self, workspace):
        global CAMC

        
        if CAMC is None:
            print "No CAMC Object. Exiting"
            return 
    	if not CAMC.isConnected():
    		print "CAM server not connected. Trying to connect."
    		CAMC.open()	
        ###
        # Find centre coordinates of object to image at high_res
        ## 
        # Get the measurements object 

        measurements = workspace.measurements
        assert isinstance(measurements, cpmeas.Measurements)

        print "current image set number ", measurements.image_set_number
        print "is first image ", measurements.is_first_image
    
    	features =  measurements.get_feature_names(self.input_object_name.value)
    	
    	if self.centerMeasurement.value==CENTER_AREASHAPECENTER and 'AreaShape_Center_X' in features:
    	    print 'Using AreaShape_Center_X/Y'
    	    xcentres = measurements.get_current_measurement(self.input_object_name.value, 'AreaShape_Center_X')
    	    ycentres = measurements.get_current_measurement(self.input_object_name.value, 'AreaShape_Center_Y')
        elif self.centerMeasurement.value==CENTER_BB and 'AreaShape_Boundingbox_X_Centre' in features:
    	    print 'Using AreaShape_Boundingbox_X/Y_Centre'
            xcentres = measurements.get_current_measurement(self.input_object_name.value, 'AreaShape_Boundingbox_X_Centre')
            ycentres = measurements.get_current_measurement(self.input_object_name.value, 'AreaShape_Boundingbox_Y_Centre')
        else:
    	    print 'Using Location_X/Y'        
            xcentres = measurements.get_current_measurement(self.input_object_name.value, M_LOCATION_CENTER_X)
            ycentres = measurements.get_current_measurement(self.input_object_name.value, M_LOCATION_CENTER_Y)
        
        
        # x and y will be arrays if there are multiple objects.
        # if no object is present they will be empty lists.

        # delete existing CAM list even when there are no objects found, otherwise
        # we end up imaging the "old" cam list again.
        if self.deleteCAMList.value:
            CAMC.deleteCAMList()
        
        if len(xcentres)==0:
            print "No object found in current image."
            # TODO how to handle this
        else:
            try: # this will only work with LeicaWaitForImage module which sets the appropriate measurements
                nx = workspace.measurements.get_current_image_measurement("Metadata_image_width")
                ny = workspace.measurements.get_current_image_measurement("Metadata_image_height")
            except: # not sure whether the following will always work
                print "No Image measurements, trying to  figure out image dimensions"
                obj = workspace.get_objects(self.input_object_name.value)
                if obj.has_parent_image:
                	parent_image = obj.get_parent_image()
                	nx, ny = parent_image.pixel_data.shape
                else:
                	nx, ny = obj.segmented.shape

            # Debug output
            #print "Image dimensions width,height: " , nx, ny
            #print "Flipping x coordinates: ", ("no","yes")[self.flipx.value]
            #print "Flipping y coordinates: ", ("no","yes")[self.flipy.value]
            #print "Swap x and y coordinates: ", ("no","yes")[self.swapxy.value]
            #print "Offset x: ", self.offsetX.value
            #print "Offset y: ", self.offsetY.value
                
            for x,y in zip(xcentres, ycentres): 
                # print some debugging info to the console:
               
                print "Object centre X coordinate :",  x
                print "Object centre Y coordinate :",  y

                if self.flipx.value:
                    x=(nx-1)-x
                if self.flipy.value:
                    y=(ny-1)-y

                if self.swapxy.value:
                    x,y = y,x # isn't python great ! we can swap variables without a temporary variable ! 

                U = str(measurements.get_current_image_measurement("Metadata_ChamberU"))
                V = str(measurements.get_current_image_measurement("Metadata_ChamberV"))
                PosX = str(measurements.get_current_image_measurement("Metadata_PosX"))
                PosY = str(measurements.get_current_image_measurement("Metadata_PosY"))
                Slide = str(measurements.get_current_image_measurement("Metadata_Slide"))
                TimePoint = str(measurements.get_current_image_measurement("Metadata_T"))

                # TODO make sure the coordinates are still within range after adding offset, either by clipping
                # the values or sending an error message

                dxpos =  `int(round(float(x)-(float(nx))/2.0)+self.offsetX.value)` # might have to use nx-1, ny-1
                dypos = `int(round(float(y)-(float(ny))/2.0)+self.offsetY.value)` # depending on whether Leica starts counting at 1 or 0 

                # Increase counter for this type of object in this well.
                # This feature can be used to limit the number of CAM jobs 
                # to be called for a particular type of object in each well.
                # To keep track of the count we create a unique string index from 
                # object name and well coordinates. 
                time_well_index =   ",".join((TimePoint, U, V, self.input_object_name.value))
                if time_well_index in self.nr_objs_in_well.keys():
                    self.nr_objs_in_well[time_well_index] +=1
                else: 
                    self.nr_objs_in_well[time_well_index] = 1
                print time_well_index, " count is ", self.nr_objs_in_well[time_well_index]

                if self.maxNrObjsPerWell.value==-1 or self.nr_objs_in_well[time_well_index] <= self.maxNrObjsPerWell.value:
                    CAMC.addJobToCAMlist(jobname=self.CAMJob.value, dxpos=dxpos, dypos=dypos, slide=`int(Slide)+1`, wellx= `int(U)+1`, welly= `int(V)+1`, fieldx=`int(PosX)+1`, fieldy=`int(PosY)+1` )
                else:
                    print "Max nr of objects to image for this well reached. Not adding to CAM list"
                    
        if self.startCAMJob.value in (CHOICE_STARTCAMJOB_DEFAULT):
            CAMC.startCAMScan()

        if self.stopWaitingForCAM.value in (CHOICE_STOPWAITING_DEFAULT):
            CAMC.stopWaitingForCAM()

    def is_interactive(self):
        return False
    
    def post_run(self, workspace):
        '''Send the startcamjob command after all images have been analysed'''
        #
        # Do nothing in test mode
        #
        if workspace.pipeline.test_mode:
            return

        global CAMC
        if CAMC is None or not CAMC.isConnected():
            print "No Connection to CAM Server. Connecting"    
        else:
             CAMC.close()

        # send Post-Run commands if selected:
        if self.startCAMJob.value in (CHOICE_STARTCAMJOB_POSTRUN):
            CAMC.startCAMScan()
        if self.stopWaitingForCAM.value in (CHOICE_STOPWAITING_POSTRUN):
            CANC.stopWaitingForCAM()
