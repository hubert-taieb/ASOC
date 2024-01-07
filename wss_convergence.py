# Import modules

# Import files that are part of the package
import parameters
import create_files
import stl_manipulation
import methods
#Import other useful modules
import time #Modules to know the current time and format it
from datetime import datetime
import os #Module to send commands to the shell
from os import listdir #List directories
from os.path import isfile,isdir, join #Checks if path is a file or a directory
import shutil #Used to copy/paste files, trees
import os
import statistics #module to do statistics on lists

# Use the trivalues from create_files to store the wss data
# The idea is to use this program during the runtime of OpenFOAM
# For each case, wss avg will be calculated for each cell and stored
# Each iteration is compared to the last one and if the relative difference between them is lower than 0.4% for all of them it's ok

# Reads the wss file passed as input and returns the number of faces associated with the target cell, the mean wss, the sdev on the value and the complete list of vectors (trivalues)
def read_wss(input_file="",cell_name=""):
    # List of trivalues (defined in create_files, basically vectors)
    wss_list=[]
    # Replica containing only the magnitudes of the wss vectors so it can be handled by the stats module
    wss_mag_list=[]

    cell_exists=False

	#Check that the file exists and is not being modified anymore (closed at least 1 second ago)
    while not os.path.exists(input_file) and (float(datetime.now().strftime('%s'))-os.path.getmtime(input_file))>1.0:
        time.sleep(0.2)

    f_handle=open(input_file)
    f_content=f_handle.read() # Read the content so it's easier to go through the lines by index
    f_lines=f_content.split("\n")
    f_handle.close()

    # Identify the block that stores the values linked to the cell of interest
    for line_index in range(len(f_lines)):
        if f_lines[line_index]=="    "+cell_name:
            N_faces=int(f_lines[line_index+4])
            start_index=line_index+6
            cell_exists=True

    if cell_exists:

        #Add each value to the list
        for line_index in range(start_index,start_index+N_faces):
            xyz=f_lines[line_index].split()
            #Retrieve the values of the kinematic shear stress (m²/s²) and convert it directly into Pascal using the density stored in parameters
            x,y,z=float(xyz[0][1:])*parameters.rho,float(xyz[1])*parameters.rho,float(xyz[2][:-1])*parameters.rho
            wss_list.append(create_files.trivalue(x,y,z))
            wss_mag_list.append(wss_list[-1].calc_mag())

        # Number of faces on the boundary, avg wss, stdev wss, list of wss values on the different faces
        return N_faces,statistics.mean(wss_mag_list),statistics.stdev(wss_mag_list),wss_list

    else:
        return None,None,None,None

#Identifies which iteration was the last one to be saved (with the WSS)
def latest_saved_it(wd=""):
    it_list=["100","200","300","400","500","600","700","800","900","1000","1100","1200","1300","1400","1500","1600","1700","1800","1900","2000"]
    latest_it,previous_it="0","0" #boolean tells if the wss has been computed at least twice (so the convergence can be evaluated)
    for ind in range(len(it_list)-1,-1,-1):
        previous_file = wd+it_list[ind-1]+"/wallShearStress"
        current_file = wd+it_list[ind]+"/wallShearStress"
        #Check the relative difference between the 2 last file to estimate if the convergence is alright
        if isfile(current_file):
            latest_it=it_list[ind]
            if isfile(previous_file): previous_it=it_list[ind-1]
            return latest_it,previous_it
    return latest_it,previous_it

#Returns a boolean that tells if the sim has reached an acceptable level of precision for the cell of interest
def wss_convergence_bool(wd="",cell_name="",reldiff_target=parameters.reldiff_target):
    latest_it,previous_it=latest_saved_it(wd=wd)
    previous_file = wd+previous_it+"/wallShearStress"
    current_file = wd+latest_it+"/wallShearStress"
    if previous_it!="0":
        wss_previous_mean=read_wss(input_file=previous_file,cell_name=cell_name)[1]
        wss_current_mean=read_wss(input_file=current_file,cell_name=cell_name)[1]
        wss_reldiff=abs((wss_current_mean-wss_previous_mean)/wss_current_mean)
        #If the criteria is respected, the sim is considered to have converged (wss wise)
        if wss_reldiff<reldiff_target:
            return True
    return False

#For a given cell in a given simulation, retrieve the number of faces of the cell, the associated WSS avg and stdev, and the complete wss list
def transfer_latest_wss(sim,Cell):
    latest_it,previous_it=latest_saved_it(wd=sim.buffer_folder)
    sim.latest_saved_it=latest_it #store it in the object
    input_file=sim.buffer_folder+str(latest_it)+"/wallShearStress"
    cell_name=sim.FOV+"_"+Cell.name+"_"+parameters.cell_channel+"_"+sim.Time
    #return read_wss(input_file=input_file,cell_name=cell_name)
    Cell.N_faces,Cell.wss_mean,Cell.wss_stdev,Cell.wss_list=read_wss(input_file=input_file,cell_name=cell_name)
