import sympy as sym
import numpy as np
import math
No_of_element = int(input('Give number of element = '))
No_of_nodes = No_of_element+2
List_of_matrices = []
List_of_force_atrices = []
support_list = ['H','NS','NS','NS','R','NS']
# a = input("\nWrite 'PS' for plane stress conditions\nAnd for plane strain condition write 'PST'\n")
# E = float(input("\nGive value of elasticity = "))
# µ = float(input("\nGive value of poisson's ratio = "))
t=float(input('\ngive the value of thickness of element='))
for index in range (0,No_of_element):
    def coordinates(n=None):
        global cordx,cordy
        cordx=[]
        cordy=[]
        if n:
            print("Please give the value of x & y as shown below, at particular node in sequence wise from 1\n i.e. 10,20")
            for index in range(0,n):
                x, y = map(float,input('\n(X{0},Y{1})='.format(index+1,index+1)).split(",")) 
                #print("\n(X{0},Y{1})=({2},{3})".format(index+1,index+1,x,y))
                cordx.append(x)
                cordy.append(y)
        else:
            A=int(input('Please give the number of total nodes\n'))
            B=print("Please give the value of x & y as shown below, at particular node in sequence wise from 1\n i.e. 10,20")
            for index in range(0,A):
                x, y = map(float,input('\n(X{0},Y{1})='.format(index+1,index+1)).split(","))
                #print("\n(X{0},Y{1})=({2},{3})".format(index+1,index+1,x,y))
                cordx.append(x)
                cordy.append(y)
        return
    coordinates(3)
    X,Y=sym.symbols("X,Y")
    U1=sym.Matrix([[1,X,Y]])
    # print("U1={0}".format(U1))
    U=np.array([[1,cordx[0],cordy[0]], [1,cordx[1],cordy[1]],[1,cordx[2],cordy[2]]])
    print("\nU={0}".format(U))
    C=np.linalg.inv(U)
    det=np.around(np.linalg.det(U),decimals=3)
    print("det={0}".format(det))
    A=np.around((det*C),decimals=2)
    print("\na1={0}\na2={1},\na3={2},\nb1={3},\nb2={4},\nb3={5},\nc1={6},\nc2={7},\nc3={8}".format(A[0,0],A[0,1],A[0,2],A[1,0],A[1,1],A[1,2],A[2,0],A[2,1],A[2,2]))
    N1=U1*A
    N=sym.matrix2numpy(N1)
    shape_Function=np.array([[N[0,0],0,N[0,1],0,N[0,2],0],[0,N[0,0],0,N[0,1],0,N[0,2]]])
    # print("Shape function\n{0}".format(shape_Function))
    print("\nN1={0}\nN2={1}\nN3={2}".format(N[0,0],N[0,1],N[0,2]))

    def strain_displacement_matrix(a=None,E=None,µ=None):
        # a=input("\nWrite 'PS' for plane stress conditions\nAnd for plane strain condition write 'PST'\n")
        # E=float(input("\nGive value of elasticity = "))
        # µ=float(input("\nGive value of poisson's ratio = "))
        global B,D 
        B = np.array([[A[1,0],0,A[1,1],0,A[1,2],0],[0,A[2,0],0,A[2,1],0,A[2,2]],[A[2,0],A[1,0],A[2,1],A[1,1],A[2,2],A[1,2]]])
        print("B matrix =\n{0}".format(B))
        if a == 'PST':
            G=1-(µ**2)
            D=np.around((E/G)*np.array([[1,µ,0], [µ,1,0],[0,0,(1-µ)/2]]),decimals=3)
            print("D matrix =\n{0}".format(D))
        else:
            G=1-(2*µ)
            G1=1+µ
            D=np.around((E/G/G1)*np.array([[1-µ,µ,0], [µ,1-µ,0],[0,0,(1-2*µ)/2]]),decimals=3)
            print("D matrix =\n{0}".format(D))
        return
    strain_displacement_matrix(a='PST',E=100000,µ=25)
    BtD=np.transpose(B).dot(D)
    BtDB=np.transpose(B).dot(D).dot(B)
    K=np.around((2*t/(4*det))*np.transpose(B).dot(D).dot(B),decimals=3)
    List_of_matrices.append(K)
    def loadvector(c=None):
        if c:
            for index in range (0,c):
                print("\nGive type of loading as shown below\nWrite 'P' for point load\nWrite 'E' for edge loading\nWrite 'S' for surface loading" )
                f=input('\nGive type of loading as shown above\ncase {0}\n'.format(index+1))
                if f=='P':
                    print("\nGive value of point load and θ as asked below\n Important note that θ value must be w.r.t the +ve x-axis ")
                    P,θ = map(float,input('\ni.e. Point load , θ\n(P,0)=').split(","))
                    x,y = map(float,input('\ni.e. give value of x and y co-ordinates, where load is applied, θ\n(x,y)=').split(","))
                    p=np.around(np.array([[P*math.cos(math.radians(θ))],[P*math.sin(math.radians(θ))]]),decimals=3)
                    print("\nP matrix =\n{0}".format(p))
                    U3=sym.Matrix([[1,x,y]])
                    U2=sym.matrix2numpy(U3)
                    sp=np.array(U2).dot(A)/det
                    print("\nN1={0}\nN2={1}\nN3={2}".format(sp[0,0],sp[0,1],sp[0,2]))
                    spf=np.array([[sp[0,0],0,sp[0,1],0,sp[0,2],0],[0,sp[0,0],0,sp[0,1],0,sp[0,2]]])
                    F=np.transpose(spf).dot(p)
                    # print("load vector = \n{0}".format(F))
                    List_of_force_atrices.append(F)
                elif f=='E':
                    print("\nGive following details for edge loading as asked below\n Write 'H' for only horizontal "
                          "loading\nWrite 'V' for only verticle loading\nWrite 'I' for inclined loading\n" )
                    H=input('\nGive details as asked above\n')
                    if H=='H':
                        P1,P2,alpha,member_number = map(float,input('\ni.e. small value , large value,alpha,member '
                                                                    'number\n(P1,P2,alpha,member number)=').split(","))
                        if member_number==1:
                            l=math.sqrt((cordx[1]-cordx[0])**2 + (cordy[1]-cordy[0])**2)
                        elif member_number==2:
                            l=math.sqrt((cordx[2]-cordx[1])**2 + (cordy[2]-cordy[1])**2)
                        elif member_number==3:
                            l=math.sqrt((cordx[2]-cordx[0])**2 + (cordy[2]-cordy[0])**2)
                        Bar_shape=sym.Matrix([[1-(X/l),0,X/l,0],[0,1-(X/l),0,X/l]])
                        Px=P1+((P2-P1)*((X/l)**alpha))
                        p=sym.Matrix([[Px],[0]])
                        F2=sym.Matrix(Bar_shape.T)
                        F1=F2*p
                        F=(sym.integrate(F1,(X,0,l)))
                        # print("load vector = \n{0}".format(sym.matrix2numpy(F)))
                        List_of_force_atrices.append(sym.matrix2numpy(F))
                    if H=='V':
                        P1,P2,alpha,member_number = map(float,input('\ni.e. small value , large value,alpha,member '
                                                                    'number\n(P1,P2,alpha,member number)=').split(","))
                        if member_number==1:
                            l=math.sqrt((cordx[1]-cordx[0])**2 + (cordy[1]-cordy[0])**2)
                            print("l={0}".format(l))
                        elif member_number==2:
                            l=math.sqrt((cordx[2]-cordx[1])**2 + (cordy[2]-cordy[1])**2)
                            print("l={0}".format(l))
                        elif member_number==3:
                            l=math.sqrt((cordx[2]-cordx[0])**2 + (cordy[2]-cordy[0])**2)
                            print("l={0}".format(l))
                        Bar_shape=sym.Matrix([[1-(X/l),0,X/l,0],[0,1-(X/l),0,X/l]])
                        Py=P1+((P2-P1)*((X/l)**alpha))
                        p=sym.Matrix([[0],[Py]])
                        F2=sym.Matrix(Bar_shape.T)
                        F1=F2*p
                        F=sym.integrate(F1,(X,0,l))
                        print("load vector = \n{0}".format(sym.matrix2numpy(F)))
                    if H=='I':
                        P1,P2,θ,alpha,member_number = map(float,input('\ni.e. small value , large value,θ w.r.t +ve '
                                                                      'x-axis,alpha,member number\n(P1,P2,θ,alpha,member'
                                                                      ' number)=').split(","))
                        if member_number==1:
                            l=math.sqrt((cordx[1]-cordx[0])**2 + (cordy[1]-cordy[0])**2)
                        elif member_number==2:
                            l=math.sqrt((cordx[2]-cordx[1])**2 + (cordy[2]-cordy[1])**2)
                        elif member_number==3:
                            l=math.sqrt((cordx[2]-cordx[0])**2 + (cordy[2]-cordy[0])**2)
                        Bar_shape=sym.Matrix([[1-(X/l),0,X/l,0],[0,1-(X/l),0,X/l]])
                        P1x=P1*math.cos(math.radians(θ))
                        P2x=P2*math.cos(math.radians(θ))
                        P1y=P1*math.sin(math.radians(θ))
                        P2y=P2*math.sin(math.radians(θ))
                        Px=P1x+((P2x-P1x)*((X/l)**alpha))
                        Py=P1y+((P2y-P1y)*((X/l)**alpha))
                        p=sym.Matrix([[Px],[Py]])
                        F2=sym.Matrix(Bar_shape.T)
                        F1=F2*p
                        F=sym.integrate(F1,(X,0,l))
                        print("load vector = \n{0}".format(sym.matrix2numpy(F)))
                elif f=='S':
                    print("\nGive following details for Surface loading as asked below\n Write 'H' for only "
                          "horizontal loading\nWrite 'V' for only verticle loading\nWrite 'I' for inclined loading\n" )
                    H=input('\nGive details as asked above\n')
                    if H=='H':
                        Sx = float(input("Give value of Sx ="))
                        p=sym.Matrix([[Sx],[0],[Sx],[0],[Sx],[0]])
                        F=p*det*t/6
                        print("load vector = \n{0}".format(sym.matrix2numpy(F)))
                    if H=='V':
                        Sy = float(input("Give value of Sy ="))
                        p=sym.Matrix([[0],[Sy],[0],[Sy],[0],[Sy]])
                        F=p*det*t/6
                        print("load vector = \n{0}".format(sym.matrix2numpy(F)))
                    if H=='I':
                        Sx1,Sy1,θ = map(float,input('\ni.e. Surface load in X-direction,Surface load '
                                                    'in Y-direction,θ\n(Sx,Sy,θ)=').split(","))
                        Sx=Sx1*math.cos(math.radians(θ))
                        Sy=Sy1*math.sin(math.radians(θ))
                        p=sym.Matrix([[Sx],[Sy],[Sx],[Sy],[Sx],[Sy]])
                        F=p*det*t/6
                        print("load vector = \n{0}".format(sym.matrix2numpy(F)))
        return
    def stresses():
        a=input("\nWrite 'PS' for plane stress conditions\nAnd for plane strain condition write 'PST'\n")
        E=float(input("\nGive value of elasticity = "))
        µ=float(input("\nGive value of poisson's ratio µ = "))
        u1,v1,u2,v2,u3,v3 = map(float,input('\nGive value of nodal displacements \n(u1,v1,u2,v2,u3,v3)=').split(","))
        global B,D,u 
        B = np.array([[A[1,0],0,A[1,1],0,A[1,2],0],[0,A[2,0],0,A[2,1],0,A[2,2]],[A[2,0],A[1,0],A[2,1],A[1,1],A[2,2],A[1,2]]])
        d = np.array([[u1],[v1],[u2],[v2],[u3],[v3]])
        DB=np.around(np.array(D).dot(B),decimals=3)
        print("D*B={0}".format(DB))
        if a == 'PST':
            G=1-(µ**2)
            D=(E/G)*np.array([[1,µ,0], [µ,1,0],[0,0,(1-µ)/2]])
            σ=np.around(np.array(D).dot(B).dot(d),decimals=3)
            print("\nσx={0}\nσy={1}\nτxy={2}".format(σ[0,0],σ[1,0],σ[2,0]))
            #principle stresses
            σmax=((σ[0,0]+σ[1,0])/2)+math.sqrt(((σ[0,0]-σ[1,0])/2)**2+(σ[2,0])**2)
            σmin=((σ[0,0]+σ[1,0])/2)-math.sqrt(((σ[0,0]-σ[1,0])/2)**2+(σ[2,0])**2)
            θp=0.5*math.degrees(math.atan(((2*σ[2,0])/(σ[0,0]-σ[1,0]))))
            print("\nσmax={0}\nσmin={1}\nθp={2}".format(σmax,σmin,θp))
        else:
            G=1-(2*µ)
            G1=1+µ
            D=(E/G/G1)*np.array([[1-µ,µ,0], [µ,1-µ,0],[0,0,(1-2*µ)/2]])
            σ=np.around(np.array(D).dot(B).dot(d),decimals=3)
            print("\nσx={0}\nσy={1}\nτxy={2}".format(σ[0,0],σ[1,0],σ[2,0]))
            #principle stresses
            σmax=((σ[0,0]+σ[1,0])/2)+math.sqrt(((σ[0,0]-σ[1,0])/2)**2+(σ[2,0])**2)
            σmin=((σ[0,0]+σ[1,0])/2)-math.sqrt(((σ[0,0]-σ[1,0])/2)**2+(σ[2,0])**2)
            θp=0.5*math.degrees(math.atan(((2*σ[2,0])/(σ[0,0]-σ[1,0]))))
            print("\nσmax={0}\nσmin={1}\nθp={2}".format(σmax,σmin,θp))
        return
    loadvector(1)
    # stresses()
Joint_stiffnessmatrix = np.around(np.zeros((2 * No_of_nodes, 2 * No_of_nodes)), decimals=3)
initial_index = 0
size_1 = 2
size_2 = 4
index_i = 2
index1_step = 4
index_j = 6
for index in range(0, len(List_of_matrices)):
    matrix = List_of_matrices[index]
    print("\n K{0} = {1}".format(index + 1, matrix))
    index_next_matrices = index + 1
    if index != len(List_of_matrices) - 1:
        next_matrices = List_of_matrices[index_next_matrices]
    if index == 0:
        Joint_stiffnessmatrix[initial_index:size_1, initial_index:size_2+2] = matrix[initial_index:size_1,
                                                                            initial_index:size_2+2]
        Joint_stiffnessmatrix[size_1:size_2+2, initial_index:size_1] = matrix[size_1:size_2+2, initial_index:size_1]
        if index != len(List_of_matrices) - 1:
            Joint_stiffnessmatrix[size_1:size_2+2, size_1:size_2+2] = matrix[size_1:size_2+2,
                                                                  size_1:size_2+2] + next_matrices[
                                                                                   initial_index:size_1+2,
                                                                                   initial_index:size_1+2]
        else:
            Joint_stiffnessmatrix[size_1:size_2+2, size_1:size_2+2] = matrix[size_1:size_2+2, size_1:size_2+2]
    else:
        Joint_stiffnessmatrix[index_i:index1_step, index1_step+2:index_j+2] = matrix[initial_index:size_1,
                                                                          size_1+2:size_2+2]
        Joint_stiffnessmatrix[index1_step+2:index_j+2, index_i:index1_step] = matrix[size_1+2:size_2+2,
                                                                          initial_index:size_1]
        if index != len(List_of_matrices) - 1:
            Joint_stiffnessmatrix[index1_step:index_j+2, index1_step:index_j+2] = matrix[size_1:size_2+2,
                                                                              size_1:size_2+2] + next_matrices[
                                                                                               initial_index:size_1+2,
                                                                                               initial_index:size_1+2]
        else:
            Joint_stiffnessmatrix[index1_step:index_j+2, index1_step+2:index_j+2] = matrix[size_1:size_2+2,
                                                                                    size_1+2:size_2+2]
            Joint_stiffnessmatrix[index1_step+2:index_j + 2, index1_step:index_j] = matrix[size_1+2:size_2+2,
                                                                                          size_1:size_2]
        index_i += 2
        index1_step += 2
        index_j += 2

print("\nSj=\n{0}".format(Joint_stiffnessmatrix))
index = 0
List_of_type_of_support = []
counter_r = 0
counter_h = 0
counter_ns = 0
for index in range(0, No_of_nodes):
    # B = input('write support type as describe above\n')
    B = support_list[index]
    intput_user = str(B).strip()
    if intput_user == 'R' or intput_user == 'r':
        counter_r += 1
    elif intput_user == 'H' or intput_user == 'h':
        counter_h += 1
    elif intput_user == 'ns' or intput_user == 'NS':
        counter_ns += 1
    List_of_type_of_support.append(B)
SFF = np.zeros(((counter_r + counter_h + (2 * counter_ns)), (counter_r + counter_h + (2 * counter_ns))))
key_pair_columns = {}
for index in range(0, len(List_of_type_of_support)):
    user_input = str(List_of_type_of_support[index]).strip().lower()
    if user_input == 'r' or user_input == 'h':
        key_pair_columns[index + index + 1] = Joint_stiffnessmatrix[:, index + index + 1]
    elif user_input == 'NS' or user_input == 'ns':
        key_pair_columns[index + index] = Joint_stiffnessmatrix[:, index + index]
        key_pair_columns[index + index + 1] = Joint_stiffnessmatrix[:, index + index + 1]
# print(">>>>>>>>>>>>>>>>3{0}".format(key_pair_columns))
index_counter = 0
for key, column in key_pair_columns.items():
    # print("{0}".format(key))
    # print("{0}".format(column))
    temp_list = []
    for key_index in key_pair_columns.keys():
        temp_list.append(column[key_index])
    SFF[:, index_counter] = temp_list
    index_counter += 1
print("\n\n SFF\n:{0}".format(SFF))