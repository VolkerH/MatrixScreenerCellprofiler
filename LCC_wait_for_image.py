####################################################################
#  LCCwaitForImage
#  Leica CAM Connect wait for image module 
#  
#  This module continuously listens to messages from the CAM server,
#  
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


'''<b>LCCwaitForImage</b> - Leica CAM Connect WaitForImage 
<hr>
This module receives notifications about new images from a Leica CAM server (running on e.g. a Leica SP5,SP8 microscope).
The images are provided to the CP pipeline. You can read in several input images that you can assign to different imaging channels.
'''

#TODO:
# * the current loop is structure is not very nice. Unless the number of images is known exactly in advance, the   post_run() functions won't be called for any module. This means that data can't be exported to CSV.




#### General imports #############################
import numpy as np
import os.path
import os.path
import sys
import pdb  

import cellprofiler.measurements as cpmeas
from cellprofiler.modules.loadimages import load_using_bioformats 
import cam_communicator_class as cc

from LCC_connection_settings import CAMC

import re
try:
    import bioformats.formatreader as formatreader
    import bioformats.metadatatools as metadatatools
    formatreader.jutil.attach()
    try:
        FormatTools = formatreader.make_format_tools_class()
        ImageReader = formatreader.make_image_reader_class()
        ChannelSeparator = formatreader.make_reader_wrapper_class(
            "loci/formats/ChannelSeparator")
        has_bioformats = True
    finally:
        formatreader.jutil.detach()
except:
    traceback.print_exc()
    has_bioformats = False


#### CP imports  #############################


import cellprofiler.cpimage as cpi
import cellprofiler.cpmodule as cpm
import cellprofiler.measurements as cpmeas
import cellprofiler.objects as cpo
import cellprofiler.settings as cps


#### Constants ##############################


# define a few global variables for this module 
leica_scanjob_nr = -1    # this will hold the number of the scanjob we're interested in - determined by first received image
leica_previous_filename = "" # keeps track of the previous file - required as CAM module sometimes sends the filenames twice 

# TODO - move into text file and import to avoid boilerplate code
MIC_WAITFORIMAGE = "Wait for image"
MIC_DONOTHING = "Do nothing"

STACK_NONE = "None"
STACK_MEAN = "Mean projection"
STACK_MAX  = "Max projection"
STACK_BEST_FOCUS = "Best focused slice"

# these are copied from loadimages
'''The FileName measurement category'''
C_FILE_NAME = "FileName"
'''The PathName measurement category'''
C_PATH_NAME = "PathName"

IP_ADDRESS_TEXT = "IP Address of CAM server"
BASEPATH_TEXT = """Path to images (on local machine, excluding "Subfolder")"""

##### Documentation #########################


doc_channel = """which channel - warning - there's no check on the actual channel number"""

doc_nrimages =  """Because Cellprofiler is originally designed as a batch processing application that works on image files, it first creates a python list structure with as many elements as the image sets to be batch processed. When images come from the microscope, Cellprofiler cannot know how many images to expect and when to terminate. With this parameter you can tell CellProfiler after how many image sets to stop processing and call post_run() events such as ExportToSpreadsheet. If you want to run indefintely, put in a very large number here but be aware that you will have to terminate manually. """

doc_stack =  """Choose a method for handling Z-Stacks. None - read single image reported by CAM. Mean/Max Projection - check for filename variations with different Z - then read all slices while performing a mean or maximum projection."""

doc_flush = """this will flush any outstanding messages from the CAM server before starting. In most cases you want to enable this option to discard any outstanding messages from previous experiments"""

doc_jobofinterest = """The Leica CAM server notifies us about newly captured images from the microscope. This incluses images from all jobs, e.g. autofocus job, hires job, lowres job, etc. Typically we are only interested in images from a particular job and want to ignore all other images. This parameter specifies the job number for the job of interest (just the number without 'J'). This job number will change with your experiment settings and you can find it in the image filename after the 'J'. If you leave this value at -1, the module will assume that the first image received is from the lowres job and take that number."""

doc_outputimage = """Specify a name with which you want to refer to this image in subsequent modules."""

doc_addchannel = """"tick if you want to read in additional image channels"""

class LCCwaitForImage(cpm.CPModule):

    module_name = "LCCwaitForImage"
    category = "MicroscopeAutomation"
    variable_revision_number = 1
    
    filename = ""

    
    ###############################################
    # create_settings 
    ################################################
    def create_settings(self):

        channels = [ str(x) for x in range(1, 10) ]

        self.channel = cps.Choice("Channel number", channels, doc = doc_channel)        

        self.nr_of_images = cps.Integer("Number of image sets to process", value = 100000, minval = -1, doc = doc_nrimages)

        self.stackOption = cps.Choice("Z Stack handling:", [STACK_NONE, STACK_MEAN, STACK_MAX], STACK_NONE,  doc = doc_stack)

        self.flush_input = cps.Binary("Flush pending messages from CAM server:", False, doc=doc_flush)

        self.job_of_interest = cps.Integer("Job number of lowres job", value = -1, minval = -1, doc = doc_jobofinterest)

        self.output_image_name = cps.ImageNameProvider("Output image name:","OutputImage",doc = doc_outputimage)

        self.ch2_active = cps.Binary("Read additional channel:", False, doc=doc_addchannel)
        self.channel2 = cps.Choice("Channel number", channels, doc = doc_channel)
        self.output_image_name_ch2 = cps.ImageNameProvider("Output image name (additional channel)", "OutputImageCh2", doc = doc_outputimage)

        self.ch3_active = cps.Binary("Read additional channel:", False, doc=doc_addchannel)
        self.channel3 = cps.Choice("Channel number", channels, doc = doc_channel)
        self.output_image_name_ch3 = cps.ImageNameProvider("Output image name (additional channel)","OutputImageCh3", doc = doc_outputimage)

        self.ch4_active = cps.Binary("Read additional channel:", False, doc=doc_addchannel)
        self.channel4 = cps.Choice("Channel number", channels, doc = doc_channel)
        self.output_image_name_ch4 = cps.ImageNameProvider("Output image name (additional channel)","OutputImageCh4", doc = doc_outputimage)

        self.ch5_active = cps.Binary("Read additional channel:", False, doc=doc_addchannel)
        self.channel5 = cps.Choice("Channel number", channels, doc = doc_channel)
        self.output_image_name_ch5 = cps.ImageNameProvider("Output image name (additional channel)","OutputImageCh5", doc = doc_outputimage)



    def settings(self):
        self.base_settings = [self.job_of_interest, self.flush_input,self.channel,  self.output_image_name, self.stackOption, self.ch2_active]
        self.ch2_settings = [self.channel2, self.output_image_name_ch2, self.ch3_active]
        self.ch3_settings = [self.channel3, self.output_image_name_ch3, self.ch4_active]
        self.ch4_settings = [self.channel4, self.output_image_name_ch4, self.ch5_active]
        self.ch5_settings = [self.channel5, self.output_image_name_ch5]
        
        #return [self.zeissaction_choice, self.microscope_choice,self.channel, self.output_image_name]
        return self.base_settings + self.ch2_settings + self.ch3_settings + self.ch4_settings + self.ch5_settings

    def visible_settings(self):
        # TODO boilerplate code very similar to settings(). How to avoid ?
        self.base_settings = [self.job_of_interest, self.flush_input, self.nr_of_images, self.channel, self.output_image_name, self.stackOption, self.ch2_active]
        self.ch2_settings = [self.channel2, self.output_image_name_ch2, self.ch3_active]
        self.ch3_settings = [self.channel3, self.output_image_name_ch3, self.ch4_active]
        self.ch4_settings = [self.channel4, self.output_image_name_ch4, self.ch5_active]
        self.ch5_settings = [self.channel5, self.output_image_name_ch5]

        if self.ch5_active.value:
            return self.base_settings + self.ch2_settings + self.ch3_settings + self.ch4_settings + self.ch5_settings
        if self.ch4_active.value:
            return self.base_settings + self.ch2_settings + self.ch3_settings + self.ch4_settings
        if self.ch3_active.value:
            return self.base_settings + self.ch2_settings + self.ch3_settings
        if self.ch2_active.value:
            return self.base_settings + self.ch2_settings
        return self.base_settings

    def prepare_run(self, pipeline, image_set_list, frame):
        """ prepare_run gets called by the cellprofiler framework. This is where you populate
        the list of images that Analyze Images will batch process. We initialize a very large image list as
        a dirty hack so that the pipeline doesn't stop before the microscope is finished with the low-resolution
        scan. You may have to increase this number for very large experiments."""

        for i in range(self.nr_of_images.value):
            image_set = image_set_list.get_image_set(i)
        return True

    def get_measurement_columns(self, pipeline):
        # Publish the metadata fields here, so we can use them in export and
        # 
        return [(cpmeas.IMAGE, "Metadata_image_width", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_image_height", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_PosX", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_PosY", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_Slide", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_ChamberU", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_ChamberV", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_Loop", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_Zpos", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_PosX", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_PosY", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_Other", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_Job", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_Channel", cpmeas.COLTYPE_INTEGER),
                (cpmeas.IMAGE, "Metadata_T", cpmeas.COLTYPE_INTEGER)
                #,
                #(cpmeas.IMAGE, "_".join((C_FILE_NAME,self.output_image_name.value), cpmeas.COLTYPE_VARCHAR_FILE_NAME)),
                #(cpmeas.IMAGE, "_".join((C_PATH_NAME,self.output_image_name.value), cpmeas.COLTYPE_VARCHAR_FILE_NAME))
                # publishing these like this appears to break everything to do with measurements in the pipeline - unclear why
                # maybe try changing COLTYPE_VARCHAR_FILE_NAME to COLTYPE_VARCHAR
                ]

    # Main 
    def run(self, workspace):
      
        output_image_name = self.output_image_name.value
        image_set = workspace.image_set
        assert isinstance(image_set, cpi.ImageSet)


        global CAMC

        # some debugging output on the console
        print "entering run in leica_interface.py"
        print "Ip address of CAM server ", CAMC.getIP()
        print "Baspath for files", CAMC.basepath
        #print self.flush_input.value

        if workspace.measurements.is_first_image: 
            print "First image of analysis run, closeing/reopening connection resetting previous_filename"
            print "closing"
            CAMC.close()
            
            print "reopening"
            CAMC.open()
            CAMC.previous_file = ""  
            CAMC.timeout=400

        # make sure we are connected
        if CAMC.isConnected():
          if self.flush_input.value:
            print("flushing input buffer on CAM server connection")
            CAMC.flushCAMreceivebuffer()
        else:
          print "No connection to cam server. Establish connection using LeicaCAMConnect module first"
          raise Exception("No Connection")
          


        # wait for an image from the microscope 
        imageresponse = CAMC.waitforimage(jobnr=self.job_of_interest.value)
        if imageresponse is not None:
            self.fullfilename, md = imageresponse # get file name and metadata
        else:
            print "No image received ... Timeout ?"
            raise Exception("Timeout")
        print "New input file ", self.fullfilename
        print "Metadata "
        print md

        #self.path, self.filename = os.path.split(self.fullfilename)

        # sample pattern we want to match: image--L0003--S00--U00--V00--J07--E00--O00--X00--Y00--T0003--Z00--C01.ome.tif
        # This is the default Leica CAM module file name pattern
        # m = re.match("(?P<Prefix>.*)(?P<Loop>--L[0-9]*)(?P<Slide>--S[0-9]*)(?P<U>--U[0-9]*)(?P<V>--V[0-9]*)(?P<Job>--J[0-9]*)(?P<E>--E.*)(?P<O>--O.*)(?P<X>--X[0-9]*)(?P<Y>--Y[0-9]*)(?P<T>--T[0-9]*)(?P<Zpos>--Z[0-9]*)(?P<Channel>--C[0-9]*)(?P<Suffix>.*)(\.ome.tif$)", self.fullfilename)

        if md is None: 
            print "Error - no metadata for file"
            raise Exception("no metadata")

        # reassemble base name
        base= md['prefix']+md['loop']+md['slide']+md['M']+md['U']+md['V']+md['job']+md['E']+md['other']+md['X']+md['Y']+md['tpoint']+md['zpos']
        

        
        def read_stack(filename):
            slicedigits = len(md['zpos'][3:])
            lastslice=int(md['zpos'][3:])
            #pdb.set_trace()

            if  self.stackOption.value == STACK_NONE:
                # Read single image
                tmpimg, scale = load_using_bioformats(
                    filename,
                    rescale = False,
                    wants_max_intensity = True)
                img = tmpimg.astype(np.float64)
            else:
                # Read stack and calculate projection on the fly
                slices = range(0,lastslice+1) 
                        
                for z in slices:
                    # cobble together filename for the current slice
                    currentslice_string = '--Z' + str(z).zfill(slicedigits)
                    slicefile = filename.replace(md['zpos'], currentslice_string)
                    #print "slicefile ", slicefile
                    # now read as usual 
                    tmpimg, tmpscale = load_using_bioformats(
                        slicefile, 
                        rescale = False,
                        wants_max_intensity = True)
                    if z == 0:
                        # copy first image and change type
                        img = tmpimg.astype(np.float64)
                        # store scale in this module
                        scale = tmpscale
                        print "data type:", tmpimg.dtype, "slice:", z, "scale:", scale
                    else:
                        if  self.stackOption.value == STACK_MEAN:
                            img += tmpimg
                            scale += tmpscale # increase scale with each slice
                        elif  self.stackOption.value == STACK_MAX:
                            img = np.maximum(img, tmpimg)
                        else:
                            print "stack option not implemented"
                            pass

            print "maxpix ", img.max()
            img /= scale
            print "maxpix after rescaling", img.max()
            return img

        # TODO: turn this into a loop or similar to avoid boilerplate code
        self.filech1 = base + "--C0"+str(int(self.channel.value)-1)+md['suffix']+".ome.tif"
        print "Reading " + self.filech1

        pixel_data=read_stack(self.filech1)
       
        if self.ch2_active.value:
            self.filech2 = base+"--C0"+str(int(self.channel2.value)-1)+md['suffix']+".ome.tif"
            print "Reading " + self.filech2
            pixel_data_ch2=read_stack(self.filech2)

        if self.ch3_active.value:
            self.filech3 = base+"--C0"+str(int(self.channel3.value)-1)+md['suffix']+".ome.tif"
            print "Reading " + self.filech3
            pixel_data_ch3=read_stack(self.filech3)

        if self.ch4_active.value:
            self.filech4 = base+"--C0"+str(int(self.channel4.value)-1)+md['suffix']+".ome.tif"
            print "Reading " + self.filech4
            pixel_data_ch4=read_stack(self.filech4)

        if self.ch5_active.value:
            self.filech5 = base+"--C0"+str(int(self.channel5.value)-1)+md['suffix']+".ome.tif"
            print "Reading " + self.filech5
            pixel_data_ch5=read_stack(self.filech5)



        # create the cellprofiler image objects from the pixel data and add the image objects to the set of images
        # Lots of boilerplate code for the differnt channels - TODO make this neater with a loop

        tmppath, tmpfile= os.path.split(self.filech1)
        output_image = cpi.Image(pixel_data, path_name=tmppath,file_name = tmpfile, scale=255)
        image_set.add(self.output_image_name.value, output_image)
        workspace.measurements.add_measurement("Image","_".join((C_FILE_NAME,self.output_image_name.value)), tmpfile, can_overwrite=True)
        workspace.measurements.add_measurement("Image","_".join((C_PATH_NAME,self.output_image_name.value)), tmppath, can_overwrite=True)
        print "added ", self.output_image_name.value, " to image_set."
        
        if self.ch2_active.value:
            tmppath, tmpfile= os.path.split(self.filech2)
            output_image_ch2 = cpi.Image(pixel_data_ch2, path_name=tmppath,file_name = tmpfile)
            image_set.add(self.output_image_name_ch2.value, output_image_ch2)
            workspace.measurements.add_measurement("Image","_".join((C_FILE_NAME,self.output_image_name_ch2.value)), tmpfile, can_overwrite=True)
            workspace.measurements.add_measurement("Image","_".join((C_PATH_NAME,self.output_image_name_ch2.value)), tmppath, can_overwrite=True)
            print "added ", self.output_image_name_ch2.value, " to image_set."

        if self.ch3_active.value:
            tmppath, tmpfile= os.path.split(self.filech2)
            output_image_ch3 = cpi.Image(pixel_data_ch3, path_name=tmppath,file_name = tmpfile)
            image_set.add(self.output_image_name_ch3.value, output_image_ch3)
            workspace.measurements.add_measurement("Image","_".join((C_FILE_NAME,self.output_image_name_ch3.value)), tmpfile, can_overwrite=True)
            workspace.measurements.add_measurement("Image","_".join((C_PATH_NAME,self.output_image_name_ch3.value)), tmppath, can_overwrite=True)
            print "added ", self.output_image_name_ch3.value, " to image_set."

        if self.ch4_active.value:
            tmppath, tmpfile= os.path.split(self.filech4)
            output_image_ch4 = cpi.Image(pixel_data_ch4, path_name=tmppath,file_name = tmpfile)
            image_set.add(self.output_image_name_ch4.value, output_image_ch4)
            workspace.measurements.add_measurement("Image","_".join((C_FILE_NAME,self.output_image_name_ch4.value)), tmpfile, can_overwrite=True)
            workspace.measurements.add_measurement("Image","_".join((C_PATH_NAME,self.output_image_name_ch4.value)), tmppath, can_overwrite=True)
            print "added ", self.output_image_name_ch4.value, " to image_set."

        if self.ch5_active.value:
            tmppath, tmpfile= os.path.split(self.filech5)
            output_image_ch5 = cpi.Image(pixel_data_ch5, path_name=tmppath,file_name = tmpfile)
            image_set.add(self.output_image_name_ch5.value, output_image_ch5)
            workspace.measurements.add_measurement("Image","_".join((C_FILE_NAME,self.output_image_name_ch5.value)), tmpfile, can_overwrite=True)
            workspace.measurements.add_measurement("Image","_".join((C_PATH_NAME,self.output_image_name_ch5.value)), tmppath, can_overwrite=True)
            print "added ", self.output_image_name_ch5.value, " to image_set."

        width = pixel_data.shape[0]
        height = pixel_data.shape[1]

        #
        #  Add Image Measurements
        #
        # Add filename and pathname to measurements, so that we can use them when saving the images later
        
        # Metadata measurements
        workspace.measurements.add_measurement("Image","Metadata_image_width", np.array(width), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_image_height", np.array(height), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_PosX", np.array(int(md['X'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_PosY", np.array(int(md['Y'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_Slide", np.array(int(md['slide'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_Other", np.array(int(md['other'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_Job", np.array(int(md['job'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_Channel", np.array(int(md['channel'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_T", np.array(int(md['tpoint'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_M", np.array(int(md['M'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_ChamberU", np.array(int(md['U'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_ChamberV", np.array(int(md['V'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_Loop", np.array(int(md['loop'][3:])), can_overwrite=True)
        workspace.measurements.add_measurement("Image","Metadata_Zpos", np.array(int(md['zpos'][3:])), can_overwrite=True)
    
    def is_interactive(self):
        return False

    def other_providers(self, group):
        if group == 'imagegroup':
            return [self.output_image_name.value]
        return []

    def is_load_module(self):
        # this is necessary for checking the pipeline, 
        # otherwise LoadSingleImage won't work
        # in the same pipeline
        return True

    def post_run(self, workspace):
        #
        # No post-run needed
        #
        return
