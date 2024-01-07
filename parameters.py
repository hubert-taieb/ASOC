### Variables

# Sim specific variables (should be modified each time)

exp_date="2020-11-33" #date of the experiment, formatted yyyy-mm-dd
exp_flow=0.1 #syringe flow (uL/min)
exp_type="12 3D_Chips_SiR_Flow_TimeLapse_MDA" #experiment type, given by the folder number within ".../02 Experiments/yyyy/"
input_folder_prefix="/run/user/10327/manuallymounted/V/Group_Cipitria/Hubert Taieb (V)/02 PhD/02 Experiments/" #folder that contains the subfolders with the stl (date, exp type dependent)

buffer_folder = "/usr/data/bgfs1/herment/01_Simulations/02_Chip_simulations/02_Simulations/" #buffer folder to perform the simulations
cell_channel="SiR" #channel used for the cells, in the name of the stl
mailtype = "none" #can be "begin","end","fail","requeue","all" or "none"
mailuser = "hubert.taieb@mpikg.mpg.de"

### SSH connection to the hot cluster
localhost_hot="hot"
username_hot="herment"
password_hot=""

####### Sim parameters

### Geometric parameters

W_channel = 40 #width of the entry channel is 40um (not used at the moment)
T_channel = 30 #thickness of the channel is 30um

### Fluid properties
rho=994.8 #kg/m3 density of water @37C --> MEASURE IT
mu=8.775e-4 #Pa.s dynamic viscosity of the media used for the MDA cells
nu=mu/rho #m2/s kinematic viscosity of the media (approx)

#Mesh parameters in um - if minCS=maxCS the mesh is uniform (except for the other refinements)
#For many different sims in parallel, minCS = 4, maxCS = 4, post_CS = 2, cell_CS = 1.5 works well
#If you have an infinite amount of time, run them at minCS = maxCS = post_CS = cell_CS = 1.5
#post_CS is the loc refinement around the posts, cell_CS is the refinement around the cells
minCS,maxCS,post_CS,cell_CS = 4,4,2,1.5 #local refinement around the cell (um)
mesh_margin = 2 #defines how much larger than the cell the refined mesh around the cell should me (in um)
mesh_killtime = 300 #300 seconds before killing the process
sim_killtime = 86400 #a whole day before killing the process
reldiff_target = 0.004 # Relative wss difference between 2 consecutive 100 iterations (convergence target)
continuity_limit = 1 # Limit value on the continuity error

################################################################################
################################ Constants #####################################
################################################################################

############################## Cluster job #####################################

#Part 1 of the cluster job file content

cjp1="#!/bin/bash -l\n\
\n\
#########################################################################################\n\
# SLURM job script\n\
##########################################################################################\n\
\n\
#SBATCH --partition=cpu,cpu2\n\
\n\
#SBATCH --job-name="

#Part 2 of the cluster job file content

cjp2='\n\
# Specify number of parallel tasks\n\
#SBATCH --ntasks=1\n\
#SBATCH --nodes=1-1\n\
\n\
# clean environment loaded at submit time\n\
module purge\n\
\n\
# load environment\n\
module add python/3.8.2\n\
module add openfoam/v2006.first\n\
module add openfoam/v2006.second\n\
\n\
### application config\n\
PROGRAM="foam"\n\
PROGRAM_OPTIONS=""\n\
\n\
### application call\n\
set\n'

#Part 3 of the cluster job file content

cjp3="# submit: sbatch <this_file_name>\n\
# job state: squeue -u <username>\n\
# queue state: sinfo\n\
# view fair-share-priority: sprio -ll\n\
# delete a job: scancel <jobnumber>\n"


############################## Probes file #####################################

pf1='/*--------------------------------*- C++ -*----------------------------------*\\\n\
| =========                 |                                                 |\n\
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
|  \\\\    /   O peration     | Version:  v2006                                 |\n\
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n\
|    \\\\/     M anipulation  |                                                 |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
\n\
// Where to load it from\n\
libs ( "libsampling.so" );\n\
\n\
type        probes;\n\
timeStart 2000;\n\
timeEnd 2000;\n\
\n\
// Name of the directory for probe data\n\
name        probes;\n\
\n\
// Fields to be probed\n\
fields\n\
(\n\
  U\n\
);\n\
\n\
probeLocations\n\
(\n'

pf2=");\n\
\n\
\n\
\n\
// ************************************************************************* //\n"

################################ meshDict file ##################################

# First part of meshDict
# Defines the mesh mostly everywhere with the parameteres specified at the top of this file
# Creates the local refinement around the posts
mDf1='/*--------------------------------*- C++ -*----------------------------------*\\\n\
| =========                 |                                                 |\n\
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
|  \\\\    /   O peration     | Version:  v2006                                 |\n\
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n\
|    \\\\/     M anipulation  |                                                 |\n\
\\*---------------------------------------------------------------------------*/\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    location  "system";\n\
    object    meshDict;\n\
}\n\
\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\
\n\
surfaceFile "model.fms";\n\
\n\
minCellSize '+str(minCS)+';\n\
\n\
maxCellSize '+str(maxCS)+';\n\
\n\
objectRefinements\n\
{\n\
    boxcenter\n\
    {\n\
        type box;\n\
        cellSize '+str(post_CS)+';\n\
        centre (1 0 15);\n\
        lengthX 30;\n\
        lengthY 58.1;\n\
        lengthZ 30;\n\
    }\n' #lengthX,lengthY,lengthZ are the size of the post

# End of the file
# Renames automatically the inlet, outlet and changes their type to patch
mDf2='}\n\
\n\
renameBoundary\n\
{\n\
    defaultName     fixedWalls;\n\
    defaultType     wall;\n\
\n\
    newPatchNames\n\
    {\n\
        "inlet.*"\n\
        {\n\
            type    patch;\n\
            newName inlet;\n\
        }\n\
\n\
        "outlet.*"\n\
        {\n\
            type    patch;\n\
            newName outlet;\n\
        }\n\
    }\n\
}\n\
\n\
// ************************************************************************* //'


################################## 0/U file ####################################

#Beginning of the files, defines the flow at the inlet (added in input_flow)
U_init1='/*--------------------------------*- C++ -*----------------------------------*\\\n\
| =========                 |                                                 |\n\
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
|  \\\\    /   O peration     | Version:  v2006                                 |\n\
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n\
|    \\\\/     M anipulation  |                                                 |\n\
\\*---------------------------------------------------------------------------*/\n\
FoamFile\n\
{\n\
    version     2.0;\n\
    format      ascii;\n\
    class       volVectorField;\n\
    object      U;\n\
}\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\
\n\
dimensions      [0 1 -1 0 0 0 0];\n\
\n\
internalField   uniform (0 0 0);\n\
\n\
boundaryField\n\
{\n\
    inlet\n\
    {\n\
        type flowRateInletVelocity;\n\
        //float must be replaced by a float value in m3/s\n\
        volumetricFlowRate '

#End of the file, boundary coundition at the outlet and no slip condition on the walls
U_init2=';\n\
        extrapolateProfile yes;\n\
    }\n\
\n\
    outlet\n\
    {\n\
        type            zeroGradient;\n\
    }\n\
\n\
    wall\n\
    {\n\
        type            noSlip;\n\
    }\n\
\n\
}\n\
\n\
// ************************************************************************* //'
