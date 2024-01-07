#Script that allows to scale and translate stl files - doesn't allow for rotation

import parameters

#Given a scaling factor sf, scale a stl file
def scale_stl(vertices_x,vertices_y,vertices_z,sf):
    vertices_x_new=[v_x*sf for v_x in vertices_x]
    vertices_y_new=[v_y*sf for v_y in vertices_y]
    vertices_z_new=[v_z*sf for v_z in vertices_z]
    return vertices_x_new,vertices_y_new,vertices_z_new

#Given a value in each direction, translate a stl file
def translate_stl(vertices_x,vertices_y,vertices_z,tx,ty,tz):
    vertices_x_new=[v_x+tx for v_x in vertices_x]
    vertices_y_new=[v_y+ty for v_y in vertices_y]
    vertices_z_new=[v_z+tz for v_z in vertices_z]
    return vertices_x_new,vertices_y_new,vertices_z_new

#Read all the vertices values present in a stl file
def read_stl(stl_file,vertices_list):
    f_in=open(stl_file)

    vertices_x,vertices_y,vertices_z=[],[],[]

    for line in f_in:
        words = line.split()
        if words[0] == 'vertex':
            vertices_x.append(float(words[1]))
            vertices_y.append(float(words[2]))
            vertices_z.append(float(words[3]))

    f_in.close()

    min_x,min_y,min_z=min(vertices_x),min(vertices_y),min(vertices_z)
    max_x,max_y,max_z=max(vertices_x),max(vertices_y),max(vertices_z)
    center_x,center_y,center_z=(max_x+min_x)/2,(max_y+min_y)/2,(max_z+min_z)/2
    size_x,size_y,size_z=max_x-min_x,max_y-min_y,max_z-min_z
    if vertices_list: return center_x,center_y,center_z,size_x,size_y,size_z,vertices_x,vertices_y,vertices_z
    else: return center_x,center_y,center_z,size_x,size_y,size_z


#Given a stl file and the new coordinates, will translate a stl file to the appropriate positions
#Also modifies the name of the solid
#Returns the new coordinates and size of the object so it can be used to create a local refinement during the meshing
def change_stl(stl_folder,stl_name,center_x_new,center_y_new,center_z_new):

    stl_file=stl_folder+stl_name
    solid_name=stl_name[:-4]
    vertices_list=True #tell the read_stl function to return the vertices list
    center_x,center_y,center_z,size_x,size_y,size_z,vertices_x,vertices_y,vertices_z=read_stl(stl_file,vertices_list)
    T_channel=parameters.T_channel #thickness of the channel in um

    #move the cell to the floor or the ceiling depending on floor_bool
    if center_z_new<T_channel/2: center_z_new=size_z/2
    else: center_z_new=T_channel-size_z/2

    tx,ty,tz=center_x_new-center_x,center_y_new-center_y,center_z_new-center_z

    vertices_x,vertices_y,vertices_z=translate_stl(vertices_x,vertices_y,vertices_z,tx,ty,tz)

    string='' #string variable that will store the content of the output file
    vertex_count=0
    f_in=open(stl_file)
    for line in f_in:
        words = line.split()
        if words[0] == 'vertex':
            string+="      vertex "+str(vertices_x[vertex_count])+" "+str(vertices_y[vertex_count])+" "+str(vertices_z[vertex_count])+"\n"
            vertex_count+=1
        #rename the solid within the file
        elif words[0] == 'solid':
            string+="solid "+solid_name+"\n"
        elif words[0] == 'endsolid':
            string+="endsolid "+solid_name+"\n"
        else:
            string+=line
    f_in.close()

    f_out=open(stl_folder+solid_name+".stl", 'w+')
    f_out.write(string)
    f_out.close()

    return center_x_new,center_y_new,center_z_new,size_x,size_y,size_z
