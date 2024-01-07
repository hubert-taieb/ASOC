#Import general libraries

import os
import shutil
import time
from datetime import datetime
import wss_convergence
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from os.path import isfile,isdir, join
import numpy as np
import statistics

#Import local libraries

import parameters

### Class that can represent any 3D-vector, usually either a position or a velocity vector
class cell_class():
    #Name="Cell_NNN", x, y, x in um
    def __init__(self,name,x,y,z):
        self.name=name
        self.x=float(x) #position of the center of the bounding box
        self.y=float(y)
        self.z=float(z)
        self.size_x=None #size of the bounding box
        self.size_y=None
        self.size_z=None
        self.input_file=None
        self.N_faces=None
        self.wss_mean=None
        self.wss_stdev=None
        self.wss_list=None
    def store_wss(self):
        return

### Class that can represent any 3D-vector, usually either a position or a velocity vector
class sim_class():
    def __init__(self,FOV,Time):
        self.FOV=FOV
        self.Time=Time
        self.name=FOV+"_"+Time
        self.Cell_list=[]
        self.buffer_folder=None
        self.output_folder=None
        self.mesh_start_time=""
        self.sim_start_time=""
        self.sim_end_time=""
        self.transfer_time=""
        self.job_id=None
        self.time_it="Unknown"
        self.latest_saved_it="Unknown"
        self.continuity_error="Unknown"
        self.wss_convergence=False
        self.continuity_convergence=True
        self.transferred=False
        self.transfer_paths=[]
        print("Simulation "+FOV+"_"+Time+" has been created.")

    # Adds a cell with its coordinate to the Cell_list
    def add_Cell(self,Cell,x,y,z):
        self.Cell_list.append(cell_class(Cell,x,y,z))

    # Writes the recap file for each cell
    def transfer(self,output_folder):
        log_string=double_print("\nStarting to transfer "+self.name+"\n")
        self.output_folder=output_folder+self.FOV+"/openfoam_files/"+self.Time+"/"
        sim_folder=self.buffer_folder

        #Go through all the cells present in the sim (if it converged)
        if self.wss_convergence and self.continuity_convergence:
            for Cell_ind in range(len(self.Cell_list)):
                Cell=self.Cell_list[Cell_ind]
                wss_convergence.transfer_latest_wss(self,Cell)
                #Write the file with the wss values
                cell_summary_file=self.FOV+"_"+Cell.name+"_"+parameters.cell_channel+"_"+self.Time+".txt"
                f_handle=open(sim_folder+cell_summary_file,"w+")
                f_handle.write("Faces = "+str(Cell.N_faces)+"\n")
                f_handle.write("WSS_mean = "+str(Cell.wss_mean)+" Pa\n")
                f_handle.write("WSS_stdev = "+str(Cell.wss_stdev)+" Pa\n")
                f_handle.write("\n")
                f_handle.write("WSS_x(Pa) WSS_y(Pa) WSS_z(Pa) WSS_mag(Pa)\n")
                for wss in Cell.wss_list:
                    f_handle.write(str(wss.x)+" "+str(wss.y)+" "+str(wss.z)+" "+str(wss.mag)+"\n")
                f_handle.close()
                #Add the summary to the list of files to be transferred
                self.transfer_paths.append(cell_summary_file)

        #List all the files to be transferred (not yet added to the list)
        self.transfer_paths.append("0") #initial conditions
        self.transfer_paths.append(str(self.latest_saved_it)) #latest values of the sim
        self.transfer_paths.append("constant")
        self.transfer_paths.append("postProcessing")
        self.transfer_paths.append("system")
        self.transfer_paths.append(self.name+"_mesh")
        self.transfer_paths.append(self.name+"_mesh.err")
        self.transfer_paths.append(self.name+"_mesh.out")
        self.transfer_paths.append(self.name+"_sim")
        self.transfer_paths.append(self.name+"_sim.out")
        self.transfer_paths.append(self.name+"_sim.err")
        self.transfer_paths.append("model.fms")
        self.transfer_paths.append("model.stl")
        self.transfer_paths.append("open_"+self.name+".foam")
        self.transfer_paths.append("trap.stl")

        #Actually transfer everything
        #Create the output tree of folders
        if not os.path.exists(self.output_folder): os.makedirs(self.output_folder)
        for path in self.transfer_paths:
            input_path,output_path=sim_folder+path,self.output_folder+path
            #should add a condition on the size of the file
            if not os.path.exists(output_path):
                if isdir(input_path): shutil.copytree(input_path,output_path)
                elif isfile(input_path): shutil.copy(input_path,output_path)
        self.transferred,self.transfer_time=True,datetime.now()
        log_string+=double_print("Transfer finished\n")

        return log_string

# Method that checks that the FOV subfolder has the right structure
def check_subfolder_FOV(folder,subfolder_name):
    #Check that the item is a folder, that it has
    if isdir(join(folder,subfolder_name)):
        subfolder_split=subfolder_name.split("_")
        #Checks that there 3 parts separated by "_", that the second part has 3 elements and the third one has 4
        if len(subfolder_split)==3 and len(subfolder_split[1])==3 and len(subfolder_split[2])==4:
            return True
    return False

# Retrieve the index of a sim within a list of simulations based on FOV, time
def sim_index(sim_list,name):
    for sim_ind in range(len(sim_list)):
        sim=sim_list[sim_ind]
        if sim.name==name:
            return sim_ind

# Sort the sim_list (list of objects of type sim_class) based on their name, alphabetically
def sort_sim_list(sim_list):
    sim_name_list=[sim.name for sim in sim_list]
    sim_name_list.sort()
    sim_list_size=len(sim_list)
    sorted_sim_list=[None]*sim_list_size
    for new_sim_ind in range(sim_list_size):
        old_sim_ind=sim_index(sim_list,sim_name_list[new_sim_ind])
        sorted_sim_list[new_sim_ind]=sim_list[old_sim_ind]
    return sorted_sim_list

# Reads the ..._sim.out file in the buffer folder and retrieves the job id and time step continuity error
def read_sim_out(sim):
    sim_out_file=sim.buffer_folder+sim.name+"_sim.out"
    #File might not exist the first time this method is called
    if os.path.exists(sim_out_file):
        sim_out_handle=open(sim_out_file)
        iteration_bool=False

        for line in sim_out_handle:
            #words=line.split()
            if line[:12]=="SLURM_JOBID=":
                sim.job_id=line.split("=")[-1].split("\n")[0]
                end=time.time()
            elif line[:7]=="Time = ":
                sim.time_it=line.split(" = ")[-1].split("\n")[0]
                iteration_bool=True
            elif iteration_bool and line[0]=="t":
                sim.continuity_error=line.split()[-1].split("\n")[0]
                if abs(float(sim.continuity_error))>parameters.continuity_limit:
                    sim.continuity_convergence=False

        sim_out_handle.close()

#Returns a string from 2 datetime.datetime() objects
def time_diff(start,end):
    diff=(end-start).total_seconds() #converts the objects in seconds
    diff_gmt=time.gmtime(diff) #converts it to the proper format
    diff_str=time.strftime("%H:%M:%S",diff_gmt)
    return diff_str

#Double print used to both print the output of the console and store it in a string
def double_print(string):
    print(string)
    return string+"\n"

#Write the definitive log file using log_string and recapping how all the sims ended (success or convergence fail)
def write_log(start_time,ASOC_name,output_folder,sim_list,log_string):
    log_string+=double_print("\n#####          ASOC global recap            #####\n")
    log_string+=double_print("ASOC ran for "+time_diff(start_time,datetime.now())+"\n")
    success_count=0
    for sim_ind in range(len(sim_list)):
        sim=sim_list[sim_ind]

        if not sim.continuity_convergence:
            log_string+=double_print("\n"+sim.name+" failed due to continuity error")
        elif not sim.wss_convergence:
            log_string+=double_print("\n"+sim.name+" did not converge to a satisfying value of WSS on one of the cell boundaries")
        else:
            success_count+=1
            if sim.transferred: log_string+=double_print(sim.name+" ran and was transferred successfully")
            else:  log_string+=double_print(sim.name+" ran successfully but was not transferred, check buffer folder for investigation")
            log_string+=double_print("Meshing ran for "+time_diff(sim.mesh_start_time,sim.sim_start_time))
            log_string+=double_print("Simulation ran for "+time_diff(sim.sim_start_time,sim.sim_end_time))

            for Cell in sim.Cell_list:
                log_string+=double_print(Cell.name+" WSS="+str(Cell.wss_mean)+"Pa, stdev="+str(Cell.wss_stdev)+"Pa")

    log_string+=double_print("Overall "+str(success_count)+"/"+str(len(sim_list))+" simulations were successful")
    output_file=output_folder+ASOC_name+"_log.txt"
    output_file_handle=open(output_file,"w+")
    output_file_handle.write(log_string)
    output_file_handle.close()
