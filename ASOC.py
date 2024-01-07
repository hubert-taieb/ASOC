################################################################################################
############################### Automated Simulations On Cluster ###############################
################################################################################################
### Author: Guillaume Herment
### Python version: 3.8.2
### Script version: 1.0 - 2020-12-07
### Runs automatically the simulations based on the variables set in parameters.py
### This script sets up the simulations in a buffer folder (V not being accessible from hot)
### It fetches the stl files obtained from the confocal, copies them in the buffer
### Based on a ..._positions.txt file it shifts them to the right position
### The mesh is created after merging them with the existing trap.stl file
### The script launches the simulations right after the meshing is over
### Once they're all launched it checks whether there is a continuity divergence
### And also whether the WSS on the cell boundaries in converging
### If everything is alright, it transfers the files of interest back to V
### Also creates a log file per FOV with relevant information
### And a log file recapping the whole process
################################################################################################

# Import files that are part of the package
import parameters
import create_files
import stl_manipulation
import methods
from methods import double_print #Import this method directly since it's used extensively
import wss_convergence
#Import other useful modules
import time #Modules to know the current time and format it
from datetime import datetime
import pickle #Module to save python objects
from pexpect import pxssh #Module to connect to the cluster in ssh
import os #Module to send commands to the shell
from os import listdir #List directories
from os.path import isfile,isdir, join #Checks if path is a file or a directory
import shutil #Used to copy/paste files, trees


#Saves all global variables in a file
def save_backup():
    with open(buffer_sd+"backup.dat","wb") as backup_file:
        pickle.dump([start_time,log_string,exp_date,exp_type,exp_year,exp_month,exp_day,input_folder,output_folder,FOV_subs,sim_list,ASOC_name,buffer_sd,sim_list],backup_file)

def load_backup(wd):
    with open(wd+"backup.dat","rb") as backup_file:
        objects=pickle.load(backup_file)
        start_time,log_string,exp_date,exp_type,exp_year,exp_month,exp_day,input_folder,output_folder,FOV_subs,sim_list,ASOC_name,buffer_sd,sim_list=objects[0],objects[1],objects[2],objects[3],objects[4],objects[5],objects[6],objects[7],objects[8],objects[9],objects[10],objects[11],objects[12],objects[13]
        return start_time,log_string,exp_date,exp_type,exp_year,exp_month,exp_day,input_folder,output_folder,FOV_subs,sim_list,ASOC_name,buffer_sd,sim_list

#If the script crashes, you can load the backup by commenting everything up until the last time the backup was saved
#Then write start_time,log_string,exp_date,exp_type,exp_year,exp_month,exp_day,input_folder,output_folder,FOV_subs,sim_list,ASOC_name,buffer_sd,sim_list=load_backup(buffer_sd) by actually replacing buffer_sd
#The variables are recharged and the script will be executed identically from then on (with the error that made it crash if you didn't solve it!)


########## Initialize global variables
start_time=datetime.now()
log_string=""

########## Create useful variables
exp_date,exp_type=parameters.exp_date,parameters.exp_type
exp_year,exp_month,exp_day=exp_date[0:4],exp_date[5:7],exp_date[8:10] #just splits the exp_date in years, month, day
input_folder=parameters.input_folder_prefix+exp_year+"/"+exp_type+"/"+exp_date+"/"
output_folder=input_folder.replace("02 Experiments","03 Simulations")

########## Setup before launching sims on the cluster
FOV_subs=[f for f in listdir(input_folder) if methods.check_subfolder_FOV(input_folder,f)] # list all the subfolders with the right structure

# Create a list of cell (class in method)
sim_list=[]
# Go through each FOV
for FOV_ind in range(len(FOV_subs)):
    FOV=FOV_subs[FOV_ind]
    FOV_folder=input_folder+FOV+"/"
    pos_file=FOV_folder+FOV+"_positions.txt"
    pos_handle=open(pos_file)

    # Go through the positions file
    for line in pos_handle:
        Cell,Time,x,y,z=line.split()
        sim_ind=methods.sim_index(sim_list,FOV+"_"+Time)
        # Create a new sim object if necessary
        if sim_ind==None:
            sim_list.append(methods.sim_class(FOV,Time))
            sim_ind=len(sim_list)-1
        # Add the cell to the sim object
        sim_list[sim_ind].add_Cell(Cell,x,y,z)
    pos_handle.close()

# Create buffer subfolder and copy the basecase
ASOC_name="Exp_"+exp_date+"_Sim_"+start_time.strftime("%Y-%m-%d_%H%M")
buffer_sd=parameters.buffer_folder+ASOC_name+"/"
if not os.path.exists(buffer_sd):
    os.mkdir(buffer_sd)
    shutil.copytree(parameters.buffer_folder+"basecase",buffer_sd+"basecase")

#Sort the sim_list based on the name of the sim
sim_list=methods.sort_sim_list(sim_list)

save_backup()

#start_time,log_string,exp_date,exp_type,exp_year,exp_month,exp_day,input_folder,output_folder,FOV_subs,sim_list,ASOC_name,buffer_sd,sim_list=load_backup("/usr/data/bgfs1/herment/01_Simulations/02_Chip_simulations/02_Simulations/Exp_2020-11-31_Sim_2020-12-10_1840/")

########## Go through all the sims that have to be launched
log_string+=double_print("\n##### Setting up the simulations parameters #####\n")
for sim_ind in range(len(sim_list)):
    sim=sim_list[sim_ind]
    sim.buffer_folder=create_files.setup_case(buffer_sd,sim.name)
    sim_folder=sim.buffer_folder
    #Go through all the cells present in the sim
    for Cell_ind in range(len(sim.Cell_list)):
        Cell=sim.Cell_list[Cell_ind]
        input_folder_cell=input_folder+sim.FOV+"/"+sim.FOV+"_amira_files/"
        stl_name=sim.FOV+"_"+Cell.name+"_"+parameters.cell_channel+"_"+sim.Time+".stl"
        #Copy the stl file of the cell in the sim folder
        shutil.copy(input_folder_cell+stl_name,sim_folder+stl_name)
        #Add the stl file to the list of files to be transferred at the end
        sim.transfer_paths.append(stl_name)
        #Move, rename the solid to the position given by ...positions.txt, retrieve the size of the solid
        sim.Cell_list[Cell_ind].size_x,sim.Cell_list[Cell_ind].size_y,sim.Cell_list[Cell_ind].size_z=stl_manipulation.change_stl(sim_folder,stl_name,Cell.x,Cell.y,Cell.z)[3:]
        log_string+=double_print(stl_name+" copied to the sim (buffer) folder, solid renamed and moved to the good position")
    #Create the meshDict based on the position, size of the cells
    create_files.setup_mesh(sim_folder,sim.Cell_list)
    #Set the propre value of flow
    create_files.input_flow(sim_folder,parameters.exp_flow)
    log_string+=double_print("Mesh dictionary created for "+sim.name+", flow will be set to "+str(parameters.exp_flow)+" uL/min\n")

save_backup()

########## Connect to the cluster, launch meshing + sim
ssh_hot = pxssh.pxssh()
if not ssh_hot.login (parameters.localhost_hot,parameters.username_hot,parameters.password_hot):
    log_string+=double_print("SSH session failed on login.")
    log_string+=double_print(str(ssh_hot))
else:
    log_string+=double_print("\n#####      SSH session login successful     #####\n")

    ########## Go through all the sims that have to be launched to start the meshing
    log_string+=double_print("\n#####     Starting to create the meshes     #####\n")
    for sim_ind in range(len(sim_list)):
        sim=sim_list[sim_ind]
        sim_folder=sim.buffer_folder
        # Pass the names of the stl file to the create_mesh function
        stl_name_list=[sim.FOV+"_"+Cell.name+"_"+parameters.cell_channel+"_"+sim.Time+".stl" for Cell in sim.Cell_list]
        create_files.create_mesh(sim_folder,sim.name,stl_name_list,ssh_hot)
        sim.mesh_start_time=datetime.now()
        log_string+=double_print("Meshing of "+sim.name+" has been launched on the cluster\n")

    save_backup()

    log_string+=double_print("\n#####       Meshing jobs all running        #####\n")

    ########## Launch the simulation only once the meshing has been done
    log_string+=double_print("\n#####  Starting to launch the simulations   #####\n")
    mesh_end_list=[]
    #The loop stops when all the simulations have finished to meshes
    #What if it fails ? It's not checked (it doesn't happen usually)
    #The sim will fail but this is handled so the script shouldn't crash anyway
    while len(mesh_end_list)!=len(sim_list):
        #A dummy file ..._mesh_end should be created once the meshing job has ended
        #The creation of this file is detected, it is deleted after the sim has been added to end_list
        end_list=create_files.retrieve_end_list(buffer_sd,sim_list)
        for sim_ind in range(len(sim_list)):
            sim=sim_list[sim_ind]
            if sim.name in end_list:
                sim_folder=sim.buffer_folder
                casename=sim.name
                create_files.launch_sim(sim_folder,casename,ssh_hot)
                mesh_end_list.append(casename)
                sim.sim_start_time=datetime.now()
                log_string+=double_print("Sim of "+casename+" has been launched on the cluster, "+str(len(sim_list)-len(mesh_end_list))+" to go")
            if len(end_list)==0: log_string+=double_print("Waiting for meshing jobs to complete, "+str(len(sim_list)-len(mesh_end_list))+" to go")
            log_string+=double_print("ASOC launched "+methods.time_diff(start_time,datetime.now())+" ago")
            time.sleep(10)

    save_backup()

    log_string+=double_print("\n#####      Simulation jobs all running      #####\n")

    ########## Go through all the sims that have been launched and check their convergence
    sim_end_list=[]
    log_string+=double_print("\n#####         Checking convergence          #####\n")
    #The loop stops once all simulations have converged or crashed
    while len(sim_end_list)!=len(sim_list):
        time.sleep(30) #better to put it at the beginning since most likely none has finished by the first iteration
        for sim_ind in range(len(sim_list)):
            sim=sim_list[sim_ind]
            if sim.name not in sim_end_list:
                sim_folder=sim.buffer_folder
                convergence_bool=True #Initialize the bool as True
                #Go through all the cells present in the sim
                for Cell_ind in range(len(sim.Cell_list)):
                    Cell=sim.Cell_list[Cell_ind]
                    #All the cell wss must converge
                    convergence_bool=convergence_bool and wss_convergence.wss_convergence_bool(wd=sim_folder,cell_name=sim.FOV+"_"+Cell.name+"_"+parameters.cell_channel+"_"+sim.Time)
                #Read the ...sim_out file to retrieve job_id, current iteration, continuity
                methods.read_sim_out(sim)
                #Add the sim to the list if it has either converged or ended (without converging)
                if convergence_bool:
                    sim.wss_convergence=convergence_bool
                    ssh_hot.sendline("scancel "+sim.job_id)
                    sim_end_list.append(sim.name)
                    log_string+=double_print("\nSimulation "+sim.name+" has converged by iteration "+sim.time_it+", Continuity error "+sim.continuity_error)
                elif os.path.exists(sim_folder+sim.name+"_sim_end"):
                    sim_end_list.append(sim.name)
                    log_string+=double_print("\nSimulation "+sim.name+" did not converge by iteration "+sim.time_it+", but ended")
                elif not sim.continuity_convergence:
                    ssh_hot.sendline("scancel "+sim.job_id)
                    sim_end_list.append(sim.name)
                    log_string+=double_print("\nSimulation "+sim.name+" was terminated because the continuity error was too large, continuity error: "+sim.continuity_error)
                else: log_string+=double_print("Simulation "+sim.name+": Current iteration "+sim.time_it+", Continuity error "+sim.continuity_error+", has not converged yet")
            #Transfer the finished sim if not done yet
            if sim.name in sim_end_list and not sim.transferred:
                sim.sim_end_time=datetime.now()
                log_string+=sim.transfer(output_folder)
        log_string+=double_print(str(len(sim_list)-len(sim_end_list))+" simulations still running, ASOC launched "+methods.time_diff(start_time,datetime.now())+" ago\n")

    save_backup()

    ssh_hot.logout() #End the session

save_backup()

#Write the log file recapping everything
methods.write_log(start_time,ASOC_name,output_folder,sim_list,log_string)
