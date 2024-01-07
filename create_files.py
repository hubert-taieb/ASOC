#!/usr/bin/python -tt

### Import libraries

#Import general libraries

import os
import shutil
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import statistics

#Import local libraries

import parameters

### Class that can represent any 3D-vector, usually either a position or a velocity vector
class trivalue():
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
        self.found=False
        self.mag=None
    #In the case of a probe, tells if the probe has been found
    def found_status(self,bool):
        self.found=bool
    #Calculates and returns the magnitude of a vector
    def calc_mag(self):
        if self.x!="NaN":
            self.mag=np.sqrt(self.x**2+self.y**2+self.z**2)
        else: self.mag="NaN"
        return self.mag
    #Returns each component in the right order to create a matrix with the right dimensions
    def current_axis(self,x_fixed,y_fixed,z_fixed):
        if x_fixed: return self.z,self.y,self.x
        elif y_fixed: return self.z,self.x,self.y
        else: return self.y,self.x,self.z

### Creates a bash script that can be launched from the terminal or via os.system
### when using python
class cluster_job():
    #initialize a new cluster job by giving it a name etc.
    def __init__(self,name,mailtype,mailuser,wd):
        self.name = name #give a name to the job
        self.mailtype = mailtype #specify if a mail should be sent
        self.mailuser = mailuser #to which address
        self.wd = wd
        self.tasks = []

    #method to add the shell command to be executed
    def add_task(self,task):
        self.tasks.append(task)

    #method to actually writer the file
    def write_file(self):
        fhandle = open(self.wd+self.name,"w+")
        #some parts of the file are stored in the parameters file (cjpX)
        fhandle.write(parameters.cjp1+self.name+"\n\
\n")
        #write the mailtype and the mail user in case there is one
        if self.mailtype != "none":
            fhandle.write("# User notify on events (begin,end,fail,requeue,all)\n\
#SBATCH --mail-type="+self.mailtype+"\n\
#SBATCH --mail-user="+self.mailuser+"\n\
\n")
        fhandle.write("# redirect output to file.out instead to default slurm-%j.out.\n\
#SBATCH --output "+self.wd+self.name+".out\n\
#SBATCH --error "+self.wd+self.name+".err\n"+parameters.cjp2)
        #write the tasks
        for task in self.tasks:
            fhandle.write(task+"\n")
        #Create a dummy file to know when the last command has been executed
        fhandle.write("touch "+self.wd+self.name+"_end\n\n"+parameters.cjp3)
        fhandle.close()
    """
    #Launch the job on the cluster
    def launch_job(self,ssh_obj=""):
		#The first option makes sense only if the code is launched from the terminal «hile being connected to hot
		#commented 2020-12-14 error compiling
		'''
		if ssh_obj=="":
			os.system("sbatch "+self.wd+self.name)
		else:
			print("sbatch "+self.wd+self.name)
			ssh_obj.sendline("sbatch "+self.wd+self.name)
		'''
		ssh_obj.sendline("sbatch "+self.wd+self.name)
	"""
    def launch_job(self,ssh_obj=""):
        print("sbatch "+self.wd+self.name)
        while not os.path.exists(self.wd+self.name):
            time.sleep(0.2)
        ssh_obj.sendline("sbatch "+self.wd+self.name)

    def wait_end(self,killtime):
        wait_begintime = time.time()
        wait_endtime = wait_begintime + killtime
        #Wait until the job on the cluster ends, or the runtime is exceeded
        #before continuing
        while time.time() < wait_endtime:
            if os.path.exists(self.wd+self.name+"_end"):
                os.remove(self.wd+self.name+"_end")
                return True
            else: time.sleep(10)
        return False


#Wait until all tasks are done before continuing to run the python script
def wait_end(wd,casename,killtime):
    wait_begintime = time.time()
    wait_endtime = wait_begintime + killtime
    #Wait until the job on the cluster ends, or the runtime is exceeded
    #before continuing
    while time.time() < wait_endtime:
        if os.path.exists(wd+casename+"_end"):
            os.remove(wd+casename+"_end")
            return True
        else: time.sleep(10)
    return False

#Returns a list of all cases for which the mesh has been createdi
#Allows to launch simulations in parallel right after the mesh is good
def retrieve_end_list(sd,sim_list):
    end_list=[]
    for sim in sim_list:
        casename=sim.name
        #Identify when the dummy end file has been created
        if os.path.exists(sd+casename+"/"+casename+"_mesh_end"):

            os.remove(sd+casename+"/"+casename+"_mesh_end")
            end_list.append(casename)
    return end_list

def retrieve_end_list_2(wd,casename_list):
    end_list=[]
    for casename in casename_list:
        #Identify when the dummy end file has been created
        if os.path.exists(wd+casename+"_end"):
            os.remove(wd+casename+"_end")
            end_list.append(casename)
    return end_list

###Various methods

# Change the value of the fluid flow at the input
def input_flow(wd,Q):
    #Convert to uL/min to m^3/s, and divide by the number of channels
    Q_met=Q*1e-9/(8*60)
    # Define the current working directory
    cwd = wd+"0/"
    # Create the 0/U file
    fhandle = open(cwd+"U","w+")
    #Write the file from parameters and input
    fhandle.write(parameters.U_init1+str(Q_met)+parameters.U_init2)
    #Close the files
    fhandle.close()

# Change the mesh setups
def setup_mesh(wd,Cell_list):
    # Define the current working directory
    cwd = wd+"system/"
    # Create meshDict
    fhandle = open(cwd+"meshDict","w+")
    # Write the first part of the file
    fhandle.write(parameters.mDf1)
    # Write the boxcell for each cell
    Cell_ind=0
    for Cell in Cell_list:
        fhandle.write("    box_"+Cell.name+"\n")
        fhandle.write("    {\n")
        fhandle.write("        type box;\n")
        fhandle.write("        cellSize "+str(parameters.cell_CS)+";\n")
        fhandle.write("        centre ("+str(Cell.x)+" "+str(Cell.y)+" "+str(parameters.T_channel/2)+");\n")
        fhandle.write("        lengthX "+str(Cell.size_x+parameters.mesh_margin)+";\n")
        fhandle.write("        lengthY "+str(Cell.size_y+parameters.mesh_margin)+";\n")
        fhandle.write("        lengthZ "+str(parameters.T_channel)+";\n")
        fhandle.write("    }\n")
        Cell_ind+=1
    # Write the second part of the file
    fhandle.write(parameters.mDf2)
    #Close the file
    fhandle.close()

# Setup the case using the basecase as a reference
def setup_case(sd,casename):
    #os.system("cp -r "+sd+"basecase/",sd+casename+"/")
    shutil.copytree(sd+"basecase/",sd+casename+"/")
    os.system("touch "+sd+casename+"/open_"+casename+".foam")
    return sd+casename+"/"

# Change the type of boundary for the input/output (wall -> patch)
# 2020-12-01 function rendered useless by the renaming done automatically in meshDict
def patch_IO(wd):
    # Define the current working directory
    cwd = wd+"constant/polyMesh/"
    # Open the existing boundary file
    fhandle = open(cwd+"boundary")
    # Create the new file
    fhandle_temp = open(cwd+"boundary_temp","w+")
    # Create a boolean to know if the current entry is within the inlet/outlet block
    IObool = False
    # Go through the entir file line by line
    for line in fhandle:
        # Identify the inlet/outlet blocks
        if line == "inlet\n" or line == "outlet\n":
            fhandle_temp.write(line)#write the line normally
            IObool = True #change the value of the boolean
        # Identify the line defining the type of the boundary when in the block
        elif IObool and line == "    type wall;\n":
            fhandle_temp.write("    type patch;\n") #replace wall by patch
            IObool = False #turn off the identifier
        # Write the line normally in all the other cases
        else:
            fhandle_temp.write(line)

    #Close the files
    fhandle_temp.close()
    fhandle.close()

    os.remove(cwd+"boundary") #delete old file
    os.rename(cwd+"boundary_temp",cwd+"boundary")

# Create a job to merge the 3D models of the chip + cell, and create the mesh

def create_mesh(wd,casename,stl_list,ssh_obj):
    #Create a new cluster job object
    meshing = cluster_job(casename+"_mesh",parameters.mailtype,parameters.mailuser,wd)
    #Concatenate trap stl file with cell stl file
    #If there is no stl file for the cells, the model is the empty trap alone
    if len(stl_list)==0:
        meshing.add_task("cat "+wd+"trap.stl >> "+wd+"model.stl")
    else:
        cat_str="cat "+wd+"trap.stl "
        for stl_name in stl_list:
            cat_str+=wd+stl_name+" "
        cat_str+=">> "+wd+"model.stl"
        meshing.add_task(cat_str)
    #Convert to FMS (openfoam file specific file format)
    meshing.add_task("surfaceFeatureEdges "+wd+"model.stl "+wd+"model.fms")
    #Create a mesh from the result
    meshing.add_task("cartesianMesh -case "+wd)
    #Scale it down as all the original files are in meters
    meshing.add_task("transformPoints -scale 1e-6 -case "+wd)
    meshing.write_file()
    #Launch the job
    meshing.launch_job(ssh_obj=ssh_obj)

### Create a job for the sim itself

def launch_sim(wd,casename,ssh_obj):
    #Launch the sim
    sim = cluster_job(casename+"_sim",parameters.mailtype,parameters.mailuser,wd)
    sim.add_task("simpleFoam -case "+wd)
    #sim.add_task("simpleFoam -postProcess -func wallShearStress -case "+wd)
    sim.write_file()
    sim.launch_job(ssh_obj)
    #sim.wait_end(parameters.sim_killtime)

### Post-process the data to retrieve the velocity fields

### Give a list of trivalue for the positions to be probed
def create_probes(wd,pos_list):
    fhandle = open(wd+"system/probes","w+")
    #some parts of the file are stored in the parameters file (pfX)
    fhandle.write(parameters.pf1)

    #Go through the list of positions
    for pos in pos_list:
        fhandle.write("  ("+pos[0]+" "+pos[1]+" "+pos[2]+")\n")

    #Write the last part of the file
    fhandle.write(parameters.pf2)

    fhandle.close()

#Creates a list of position from a file (obtained following fluo beads)
def create_pos_list_from_file(file_path,z_list):

    f_handle_trajectories=open(file_path)
    line_count=0
    pos_list=[]

    #Reads the file and create a list of probes with the defined z-coordinate
    for line in f_handle_trajectories:
        if line_count!=0:
            words=line.split()
            for z in z_list:
                pos_list.append([str(round(float(words[2])*1e-6,7)),str(round(float(words[3])*1e-6,7)),str(z)])
        line_count+=1
    return pos_list

#Creates the same file as before but adds the sim output
def create_file_output(file_path,probe_list,probe_content_list,z_list):
    f_handle_trajectories=open(file_path)
    f_handle_output=open(file_path[:-4]+"_sim.txt","w+")
    line_count=0
    pos_list=[]
    for line in f_handle_trajectories:
        if line_count!=0:
            probe_ind=(line_count-1)*len(z_list)
            words=line.split()
            loc_x,loc_y=round(float(words[2])*1e-6,7),round(float(words[3])*1e-6,7)
            ux_list,uy_list,umag_list=[],[],[]
            if probe_list[probe_ind].x==loc_x and probe_list[probe_ind].y==loc_y:
                z_ind=0
                for z in z_list:
                    if probe_list[probe_ind+z_ind].z==z:
                        ux_list.append(probe_content_list[0][probe_ind+z_ind].x*1e6)
                        uy_list.append(probe_content_list[0][probe_ind+z_ind].y*1e6)
                        umag_list.append(probe_content_list[0][probe_ind+z_ind].calc_mag()*1e6)
                    z_ind+=1
                ux_sim,sdev_ux_sim,uy_sim,sdev_uy_sim,u_sim,sdev_u_sim=statistics.mean(ux_list),statistics.stdev(ux_list),statistics.mean(uy_list),statistics.stdev(uy_list),statistics.mean(umag_list),statistics.stdev(umag_list)
                line_str=line[:-1]+" "+str(ux_sim)+" "+str(sdev_ux_sim)+" "+str(uy_sim)+" "+str(sdev_uy_sim)+" "+str(u_sim)+" "+str(sdev_u_sim)+"\n"
                print(words[0],words[4],str(ux_sim),words[5],str(uy_sim),float(words[4])/ux_sim,float(words[5])/uy_sim)

            else:
                print("No correspondance between input file and probed locations")
        else:
            line_str="idx_id times(s) x(µm) y(µm) ux(µm/s) uy(µm/s) u(µm/s) ux_sim(µm/s) sdev_ux_sim(µm/s) uy_sim(µm/s) sdev_uy_sim(µm/s) u_sim(µm/s) sdev_u_sim(µm/s)\n"
        #print(line_str)
        line_count+=1
        f_handle_output.write(line_str)
    f_handle_output.close()
    #return pos_list

#Creates a list of coordinates from a probe_list
def create_probelist_file(wd,probe_list,probe_content_list):
    f_handle_output=open(wd+"/positions.txt","w+")
    line_str="x(µm) y(µm) z(µm)\n"
    for probe_ind in range(len(probe_list)):
        line_str+=str(round(probe_list[probe_ind].x*1e6,1))+" "+str(round(probe_list[probe_ind].y*1e6,1))+" "+str(round(probe_list[probe_ind].z*1e6,1))+"\n"
    f_handle_output.write(line_str)
    f_handle_output.close()

#Creates a list of only velocities
def create_velocitylist_file(wd,flow,probe_list,probe_content_list):
    f_handle_output=open(wd+"/"+flow+".txt","w+")
    line_str="ux(µm/s) uy(µm/s) uz(µm/s) umag(µm/s)\n"
    for probe_ind in range(len(probe_list)):
        #line_str+=str(round(probe_list[probe_ind].x*1e6,1))+" "+str(round(probe_list[probe_ind].y*1e6,1))+"\n"
        ux,uy,uz,umag=probe_content_list[0][probe_ind].x,probe_content_list[0][probe_ind].y,probe_content_list[0][probe_ind].z,probe_content_list[0][probe_ind].calc_mag()
        if ux!="NaN": ux,uy,uz,umag=str(round(ux*1e6,3)),str(round(uy*1e6,3)),str(round(uz*1e6,3)),str(round(umag*1e6,3))
        else: ux,uy,uz,umag="NaN","NaN","NaN","NaN"
        line_str+=ux+" "+uy+" "+uz+" "+umag+"\n"
    f_handle_output.write(line_str)
    f_handle_output.close()

#Create an uniform mesh of probes within a box, with a given resolution
def create_probes_box(wd,x_min,x_max,y_min,y_max,z_min,z_max,res):
    pos_list=[]
    x_it=int((x_max-x_min)/res)+1
    y_it=int((y_max-y_min)/res)+1
    z_it=int((z_max-z_min)/res)+1

    for z_ind in range(z_it):
        for y_ind in range(y_it):
            for x_ind in range(x_it):
                x,y,z=str(round(x_ind*res+x_min,7)),str(round(y_ind*res+y_min,7)),str(round(z_ind*res+z_min,7))
                pos_list.append([x,y,z])
    """
    #Comment previous block and uncomment this one to fix z at 15um (middle of the channel)
    for y_ind in range(y_it):
        for x_ind in range(x_it):
            x,y,z=str(round(x_ind*res+x_min,7)),str(round(y_ind*res+y_min,7)),str(1.5e-5)
            pos_list.append([x,y,z])
    """
    create_probes(wd,pos_list)

#Create an uniform mesh of probes within a box, with a given resolution in each direction
def create_probes_box_2(wd,x_min,x_max,y_min,y_max,z_min,z_max,res_x,res_y,res_z):
    pos_list=[]
    x_it=int((x_max-x_min)/res_x)+1
    y_it=int((y_max-y_min)/res_y)+1
    z_it=int((z_max-z_min)/res_z)+1

    for z_ind in range(z_it):
        for y_ind in range(y_it):
            for x_ind in range(x_it):
                x,y,z=str(round(x_ind*res+x_min,7)),str(round(y_ind*res+y_min,7)),str(round(z_ind*res+z_min,7))
                pos_list.append([x,y,z])
    """
    #Comment previous block and uncomment this one to fix z at 15um (middle of the channel)
    for y_ind in range(y_it):
        for x_ind in range(x_it):
            x,y,z=str(round(x_ind*res+x_min,7)),str(round(y_ind*res+y_min,7)),str(1.5e-5)
            pos_list.append([x,y,z])
    """
    create_probes(wd,pos_list)

#Post-process the probes on the cluster - cause why not ?
def post_process_probes(wd,casename="post_process_probes"):
    ppp = cluster_job(casename,parameters.mailtype,parameters.mailuser,wd)
    ppp.add_task("postProcess -func probes -case "+wd)
    ppp.write_file()
    ppp.launch_job()
    #ppp.wait_end(parameters.sim_killtime)

#Read the output file obtained after post-processing using probes
#Still some things do deal with regarding the "time"
def read_probes(wd,field,time):

    #Create a list for the positions of each probe
    probe_list=[]
    #Create a counter to identify the current probe
    #probe_index=0
    #Create a list of lists containing the field values
    probe_content_list=[]
    #Create a list to store the timestamps
    time_list=[]
    #Create a boolean to identify when the actual data starts being red
    time_bool=False
    fhandle = open(wd+"postProcessing/probes/"+time+"/"+field)
    #Go through the lines of the file

    for line in fhandle:
        #Split them into words
        words=line.split()
        #Use the keyword "Time" to identify when the list of probes stops and the data starts

        if not time_bool:

            if words[1]=="Probe" and words[3][0]=="(":
                #List the probes only if they were found within the mesh
                probe_list.append(trivalue(float(words[3][1:]),float(words[4]),float(words[5][:-1])))
                #Change the "found" id of the current probe if it's been found

                if len(words)==6:
                    probe_list[-1].found_status(True)
                #probe_index+=1
            #Identify when the probed data starts

            elif words[1]=="Time":
                time_bool=True

        elif time_bool:
            probe_content_list.append([])

            for probe_index in range(len(probe_list)):

                #will need to be field specific (scalar vs vector)
                if probe_list[probe_index].found:
                    probe_content_list[len(time_list)].append(trivalue(float(words[3*probe_index+1][1:]),float(words[3*probe_index+2]),float(words[3*probe_index+3][:-1])))

                #add a null value in case the probe hasn't been found
                else:
                    probe_content_list[len(time_list)].append(trivalue("NaN","NaN","NaN"))

            time_list.append(int(words[0]))
            time_bool=False

    fhandle.close()

    return probe_list,probe_content_list


#Re-organizes the data into a matrix (list of list) that can be used as an entry for a colormap
def create_mag_list(probe_list,probe_content_list,x_min,x_max,y_min,y_max,z_min,z_max):

    #Identify the iso-plane and use booleans and a common variable to store the info
    if x_min==x_max:
        x_fixed=True
        y_fixed=False
        z_fixed=False
        min_1,max_1=z_min,z_max
        min_2,max_2=y_min,y_max
        fixed_3=x_min
    elif y_min==y_max:
        x_fixed=False
        y_fixed=True
        z_fixed=False
        min_1,max_1=z_min,z_max
        min_2,max_2=x_min,x_max
        fixed_3=y_min
    elif z_min==z_max:
        x_fixed=False
        y_fixed=False
        z_fixed=True
        min_1,max_1=y_min,y_max
        min_2,max_2=x_min,x_max
        fixed_3=z_min
    else:
        print("At least one of the coordinates must be fixed to generate a 2D-matrix from 3D-data")

    #Define an absurd initial value to avoid conflicts
    last_1=1e300
    mag_list=[]

    # Generate some data to plot
    for probe_ind in range(len(probe_list)):
        #Retrieve the values for the 2 non-fixed axis 1 and 2
        #Axis (1, 2) are either (y,x),(z,x) or (z,y)
        #Axis 3 is the fixed one
        current_1,current_2,current_3=probe_list[probe_ind].current_axis(x_fixed,y_fixed,z_fixed)
        #Check that the current values are within the desired bounds
        if current_3==fixed_3:
            if min_1 <= current_1 <= max_1:
                if last_1==current_1:
                    if min_2 <= current_2 <= max_2:
                        mag_list[-1].append(probe_content_list[-1][probe_ind].calc_mag())
                #add a new vector to the matrix when the coordinate of axis  has changed
                else:
                    last_1=current_1
                    mag_list.append([])

    return mag_list

#Does the same thing as the first but in a more versatile way
def create_mag_list_2(probe_list,probe_content_list,x_min,x_max,y_min,y_max,z_min,z_max,res):

    #Identify the iso-plane and use booleans and a common variable to store the info
    if x_min==x_max:
        x_fixed=True
        y_fixed=False
        z_fixed=False
        min_1,max_1=z_min,z_max
        min_2,max_2=y_min,y_max
        fixed_3=x_min
        amp_1=3e-5
        amp_2=6e-4
    elif y_min==y_max:
        x_fixed=False
        y_fixed=True
        z_fixed=False
        min_1,max_1=z_min,z_max
        min_2,max_2=x_min,x_max
        fixed_3=y_min
        amp_1=3e-5
        amp_2=9e-4
    elif z_min==z_max:
        x_fixed=False
        y_fixed=False
        z_fixed=True
        min_1,max_1=y_min,y_max
        min_2,max_2=x_min,x_max
        fixed_3=z_min
        amp_1=6e-4
        amp_2=9e-4
    else:
        print("At least one of the coordinates must be fixed to generate a 2D-matrix from 3D-data")

    N_1=int((max_1-min_1)/res)
    N_2=int((max_2-min_2)/res)

    mag_list=[]

    #Create the matrix beforhand
    for coord_1 in range(N_1):
        mag_list.append([])
        for coord_2 in range(N_2):
            mag_list[-1].append(0.0)

    # Generate some data to plot
    for probe_ind in range(len(probe_list)):
        #Retrieve the values for the 2 non-fixed axis 1 and 2
        #Axis (1, 2) are either (y,x),(z,x) or (z,y)
        #Axis 3 is the fixed one
        current_1,current_2,current_3=probe_list[probe_ind].current_axis(x_fixed,y_fixed,z_fixed)

        #Retrieve the values for the 2 non-fixed axis 1 and 2
        #Axis (1, 2) are either (y,x),(z,x) or (z,y)
        #Axis 3 is the fixed one
        current_1,current_2,current_3=probe_list[probe_ind].current_axis(x_fixed,y_fixed,z_fixed)
        ind_1=int((current_1+amp_1/2)/res)
        ind_2=int((current_2+amp_2/2)/res)
        #print(ind_1,ind_2)
        mag_list[ind_1][ind_2]=probe_content_list[-1][probe_ind].calc_mag()

    #print(mag_list)

    return mag_list



#Using a list of probes and their associated content, plus a plane, this function creates a colormap of the field
#This function rests on create_mag_list, so in one direction the min and the max must be equal
def create_colormap(wd,probe_list,probe_content_list,x_min,x_max,y_min,y_max,z_min,z_max,res,save_fig=True):

    mag_list=create_mag_list_2(probe_list,probe_content_list,x_min,x_max,y_min,y_max,z_min,z_max,res)

    # Create figure object
    fig = plt.figure()

    # make a color map of fixed colors
    cmap = mpl.cm.jet


    # tell imshow about color map so that only set colors are used
    img = plt.imshow(mag_list, origin='lower',cmap=cmap)

    # make a color bar
    cbar=plt.colorbar(img, cmap=cmap)
    cbar.set_label("Velocity (m/s)")

    if save_fig:
        plt.savefig(wd+"X_"+str(x_min)+"_"+str(x_max)+"Y_"+str(y_min)+"_"+str(y_max)+"Z_"+str(z_min)+"_"+str(z_max),format='png',dpi=300,quality=95)
    else:
        plt.show()


### Some pieces of code that relate to PIV / probing the velocity field

######## Create a csv with x,y,z,ux,uy,uz,umag,p

"""
wd="/usr/data/bgfs1/herment/01_Simulations/02_Chip_simulations/01_Testbed/automatization_testbed/20201104-flowtest/2020-11-04_1925_Q1.00/"
file="C3T1_1ul.min_10x_data_speed_to_fetch.txt"
z_list=[1.5e-6,3e-6,4.5e-6,6e-6,7.5e-6,9e-6,10.5e-6,12e-6,13.5e-6,15e-6]
x_min,x_max,y_min,y_max,z_min,z_max,res=-4.5e-4,4.5e-4,-3e-4,3e-4,1.5e-5,1.5e-5,1e-6

pos_list=create_pos_list_file(wd+file,z_list)

#create_probes(wd,pos_list)
#post_process_probes(wd)

probe_list,probe_content_list=read_probes(wd,"U","2000")
#create_file_output(wd+file,probe_list,probe_content_list,z_list)

#probe_list,probe_content_list=read_probes(wd,"U","2000")
create_colormap(wd,probe_list,probe_content_list,x_min,x_max,y_min,y_max,z_min,z_max,res)
"""
#create_colormap(probe_list,probe_content_list,3e-6,-3e-4,-4.5e-4)

#should manage the pressure as well as the velocity
#also automatize the post-processing (job etc) - maybe create some probes in advance in the 3 directions, crossing the center of the cell
'''
wd='/usr/data/bgfs1/herment/01_Simulations/02_Chip_simulations/01_Testbed/automatization_testbed/20201104-flowtest/2020-11-04_1925_Q7.00/'
for z_ind in range(10):
    for y_ind in range(200):
        for x_ind in range(300):
            x,y,z=str(round(x_ind*3e-6-4.5e-4,6)),str(round(y_ind*3e-6-3e-4,6)),str(round(z_ind*3e-6,6))
            #x,y,z="{:.2e}".format(round(x_ind*1e-6-4.5e-4,6)),"{:.2e}".format(round(y_ind*1e-6-3e-4,6)),"{:.2e}".format(round(z_ind*1e-6,6))
            #x,y,z=remove_useless_zeros(x),remove_useless_zeros(y),remove_useless_zeros(z)
            pos_list.append([x,y,z])
create_probes(wd,pos_list)
'''
