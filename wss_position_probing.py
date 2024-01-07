#As of 20201112 - script unfinished
#Should be able to link a value of wss on a boundary to its position in the 3D model
#superdirectory
sd="/usr/data/bgfs1/herment/01_Simulations/02_Chip_simulations/01_Testbed/automatization_testbed/20201111-positiontest/2020-11-11_1322_x035_y042/"
#boundary of interest
int_boundary="walls"

class face_wss():
    def __init__(self,face):
        self.face=face
        self.points=[]
        self.points_coord=[]
        self.wss=[]
    def add_points(self,line):
        xyz=line.split()
        x,y,z=float(xyz[0][1:]),float(xyz[1]),float(xyz[2][:-1])
        mag=np.sqrt(x**2+y**2+z**2)
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)
        self.mag.append(mag)
    def mean(self):
        return statistics.mean(self.x),statistics.mean(self.y),statistics.mean(self.z),statistics.mean(self.mag)
    def stdev(self):
        return statistics.stdev(self.x),statistics.stdev(self.y),statistics.stdev(self.z),statistics.stdev(self.mag)

boundary_start_face=0 #start face of the boundary of interest (as in the boundar file)
boundary_N_faces=0 #number of faces in the boudary of interest

polyMesh_path=sd+"/constant/polyMesh/"
boundary_path=polyMesh_path+"boundary"

#Open the boundary file
f_handle_boundary=open(boundary_path)
#Identifier of the int bound block start
int_boundary_bool=False

#Go through the boundary file to find the one of interest
for line in f_handle_boundary:
    words=line.split()
    if len(words)!=0 and words[0]==int_boundary: int_boundary_bool=True
    elif int_boundary_bool:
        if words[0]=='nFaces': boundary_N_faces=int(words[1][:-1])
        elif words[0]=='startFace':
            boundary_start_face=int(words[1][:-1])
            int_boundary_bool=False

#Close the boundary file
f_handle_boundary.close()



### Now let's go check the faces file
faces_path=polyMesh_path+"faces"

#The first face is always written at the line 21
line_face_0=21
line_face_counter=0

#Create the list that stores the faces and the associates points, wss
face_wss_list=[]

#Open the faces file
f_handle_faces=open(faces_path)

for line in f_handle_faces:
    line_face_counter+=1
    loc_face_index=line_face_counter-line_face_0-boundary_start_face
    #print(loc_face_index)
    #Check that the indexes fit
    if 0 <= loc_face_index < boundary_N_faces:

        words=line.split()

        face_wss_list.append(face_wss(loc_face_index))

        for word in words:
            #print(word)
            if word[1]=="(": face_wss_list[-1].points.append(int(word[2:]))
            elif word[-1]==")": face_wss_list[-1].points.append(int(word[:-1]))
            else: face_wss_list[-1].points.append(int(word))

#Close the faces file
f_handle_faces.close()




### Now let's go check the points file
points_path=polyMesh_path+"points"

#Open the points file
f_handle_points=open(boundary_path)

#The first point is always written at the line 21
line_N_points=19
line_point_0=21
line_point_counter=0
boundary_N_points=0

#Create the list that stores the faces and the associates points, wss
points_list=[]

#Open the faces file
f_handle_points=open(points_path)


for line in f_handle_points:
    line_point_counter+=1

    #Extract the number of points
    if line_point_counter == line_N_points: boundary_N_points=int(line)
    #Go the points
    elif line_point_0 <= line_point_counter < line_point_0 + boundary_N_points:

        points_list.append([])

        words=line.split()
        #Create the list containing the coordinates of the points
        points_list[-1].append(float(words[0][1:]))
        points_list[-1].append(float(words[1]))
        points_list[-1].append(float(words[2][:-1]))


#Close the points file
f_handle_points.close()


#Add the coordinates of the points of the face
for face_wss in face_wss_list:
    for point in face_wss.points:
        face_wss.points_coord.append(points_list[point])




###All that's left to do is to retrieve the associated wss value, run a double/triple for loop on coordinates while checking all the probes to identify those we want, list the wss, average the shit out of it and boom create the colormap
