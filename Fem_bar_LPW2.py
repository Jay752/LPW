##########################################################################################################################################################
## developed by: Jaykumar Viradiya(JV)
##########################################################################################################################################################
import sympy
from sympy import *
import sympy as sym
from sympy import pprint
import numpy as np
from numpy import *
import math
from scipy.misc import derivative
from scipy.integrate import quad
from Fem_bar_LPW2_function import *
init_printing(num_columns=100)
print("write '2' for two node bar element\nwrite '3' for three node bar element")
D=int(input('write number of nodes for bar element\n'))
NO_OF_ELEMENT=int(input('write number of members\n'))
length_of_int_element = []
length_of_element = []
length_of_element_1 = []
X,L = sym.symbols('X,L')
List_of_Stiffness_matrices = []
List_of_Force_vector = []
if D==2:
    for index in range (0,NO_OF_ELEMENT):
        print("****************for member = {0}*************\n\n ".format(index+1))
        J=float(input('Give length of member {0} = '.format(index+1)))
        length_of_int_element.append(J)
        length_of_element.append(index*L)
        length_of_element_1.append((index+1)*L)
        U_x=sym.Matrix([[1,X]])
        U_xx=sym.Matrix([[1,0],[1,L]])
        In = U_xx.inv()
        U=U_x*In
        print("Shape Function for member {0}= \n".format(index+1))
        pprint(U, use_unicode=false, wrap_line=false)
        B=derive_by_array(U,X)
        print("B matrix for member {0}= \n".format(index+1))
        pprint(B, use_unicode=false, wrap_line=false)
        print("\nWrite 'E' if elasticity is not same throughout the bar\nWrite 'A' if area is not same throughout the bar\nIf all three properties are same throughout the bar then write 'SAP'\nIf there are temprature variation write 'T'")
        F=input('Write name of property which has not same value along the member {0} as given in above list\n'.format(index+1))
        if F in ('E','e'):
            α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
            EE1=float(input('give Value of Elasticity at left end (E1) in N/mm^2\n'))
            EE2=float(input('give Value of Elasticity at right end (E2) in N/mm^2\n')) 
            AA =float(input('give Value of area in mm^2\n'))
            E1,E2,A = sym.symbols('E1,E2,A')
            Ex = E1+((E2-E1)*((X/L)))
            j=sym.Matrix([B])
            N=j.T*j*Ex*A
            K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
            K1=sym.lambdify([L,E1,E2,A],K)
            K2=K1(J,EE1,EE2,AA)
            List_of_Stiffness_matrices.append(K2)
            #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
            print("Stiffness matrix K{0} = \n".format(index+1))
            pprint(K, use_unicode=false, wrap_line=false)
        if F in ('A','a'):
            α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
            AA1=float(input('give Value of area at left end (A1) = \n'))
            AA2=float(input('give Value of area at right end (A2) =\n')) 
            EE = float(input('give Value of Elasticity E =\n'))
            A1,A2,E = sym.symbols('A1,A2,E')
            Ax = A1+((A2-A1)*((X/L)**α))
            j=sym.Matrix([B])
            N=j.T*j*E*Ax
            K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
            K1=sym.lambdify([L,A1,A2,E],K)
            K2=K1(J,AA1,AA2,EE)
            List_of_Stiffness_matrices.append(K2)
            #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
            print("Stiffness matrix K{0} = \n".format(index+1))
            pprint(K, use_unicode=false, wrap_line=false)
        if F in ('SAP','sap'):
            EEE = float(input('give Value of Elasticity E =\n'))
            AAA = float(input('give Value of area in mm^2\n'))
            A,E = sym.symbols('A,E')
            j=sym.Matrix([B])
            N=j.T*j*E*A
            K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
            K1=sym.lambdify([L,A,E],K)
            K2=K1(J,AAA,EEE)
            List_of_Stiffness_matrices.append(K2)
            #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
            print("Stiffness matrix K{0} = \n".format(index+1))
            pprint(K, use_unicode=false, wrap_line=false)
        
        if F in ('T','t'):
            print("\nWrite 'E' if elasticity is not same throughout the bar\nWrite 'A' if area is not same throughout the bar\nIf all three properties are same throughout the bar then write 'SAP'\n")
            Z=input('Write name of property which has not same value along the member {0} as given in above list\n'.format(index+1))
            ΔΔT=float(input('Temprature variation '))
            α_00 = float(input('Temprature Co-efficient α_0 =  '))
            if Z in ('E','e'):
                α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
                EE1=float(input('give Value of Elasticity at left end (E1) in N/mm^2\n'))
                EE2=float(input('give Value of Elasticity at right end (E2) in N/mm^2\n')) 
                AA =float(input('give Value of area in mm^2\n'))
                E1,E2,A = sym.symbols('E1,E2,A')
                Ex = E1+((E2-E1)*((X/L)))
                j=sym.Matrix([B])
                N=j.T*j*Ex*A
                K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
                K1=sym.lambdify([L,E1,E2,A],K)
                K2=K1(J,EE1,EE2,AA)
                List_of_Stiffness_matrices.append(K2)
                #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
                print("Stiffness matrix K{0} = \n".format(index+1))
                pprint(K, use_unicode=false, wrap_line=false)
                ΔT,α_0 = sym.symbols('ΔT,α_0')
                F = integrate(j.T*A*Ex*ΔT*α_0,(X,(length_of_element[index]),length_of_element_1[index]))
                F1 = sym.lambdify([L,E1,E2,A,ΔT,α_0],F)
                F2 = F1(J,EE1,EE2,AA,ΔΔT,α_00)
                print("Force matrix F{0}= \n".format(index+1))
                pprint(F,use_unicode=false,wrap_line=false)
                List_of_Force_vector.append(F2)
            if Z in ('A','a'):
                α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
                AA1=float(input('give Value of area at left end (A1) = \n'))
                AA2=float(input('give Value of area at right end (A2) =\n')) 
                EE = float(input('give Value of Elasticity E =\n'))
                A1,A2,E = sym.symbols('A1,A2,E')
                Ax = A1+((A2-A1)*((X/L)**α))
                j=sym.Matrix([B])
                N=j.T*j*E*Ax
                K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
                K1=sym.lambdify([L,A1,A2,E],K)
                K2=K1(J,AA1,AA2,EE)
                List_of_Stiffness_matrices.append(K2)
                #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
                print("Stiffness matrix K{0} = \n".format(index+1))
                pprint(K, use_unicode=false, wrap_line=false)
                ΔT,α_0 = sym.symbols('ΔT,α_0')
                F = integrate(j.T*Ax*E*ΔT*α_0,(X,(length_of_element[index]),length_of_element_1[index]))
                F1 = sym.lambdify([L,A1,A2,E,ΔT,α_0],F)
                F2 = F1(J,AA1,AA2,EE,ΔΔT,α_00)
                print("Force matrix F{0}= \n".format(index+1))
                pprint(F,use_unicode=false,wrap_line=false)
                List_of_Force_vector.append(F2)
            if Z in ('SAP','sap'):
                α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
                EEE = float(input('give Value of Elasticity E =\n'))
                AAA = float(input('give Value of area in mm^2\n'))
                A,E = sym.symbols('A,E')
                j=sym.Matrix([B])
                N=j.T*j*E*A
                K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
                K1=sym.lambdify([L,A,E],K)
                K2=K1(J,AAA,EEE)
                List_of_Stiffness_matrices.append(K2)
                #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
                print("Stiffness matrix K{0} = \n".format(index+1))
                pprint(K, use_unicode=false, wrap_line=false)
                ΔT,α_0 = sym.symbols('ΔT,α_0')
                F = integrate(j.T*A*E*ΔT*α_0,(X,(length_of_element[index]),length_of_element_1[index]))
                F1 = sym.lambdify([L,A,E,ΔT,α_0],F)
                F2 = F1(J,AAA,EEE,ΔΔT,α_00)
                print("Force matrix F{0}= \n".format(index+1))
                pprint(F,use_unicode=false,wrap_line=false)
                List_of_Force_vector.append(F2)  
        Force = input('Give value of loading as asked below\nWrite "P" for point load\nWrite "S" for same load throughout the member or for UDL on member throughout length\nWrite "D" for unsymetrical loading\n')
        if Force in ('s,S'):
            W1 = float(input('Give value of load = '))
            W = sym.symbols('W')
            F = integrate(U.T*W,(X,(length_of_element[index]),length_of_element_1[index]))
            F1 = sym.lambdify([L,W],F)
            F2=F1(J,W1)
            print("Force matrix F{0}= \n".format(index+1))
            pprint(F,use_unicode=false,wrap_line=false)
            List_of_Force_vector.append(F2)
        elif Force in ('D,d'):
            WW1 = float(input('Give value of load at left end W1= '))
            WW2 = float(input('Give value of load at right end W2= '))
            W1,W2 = sym.symbols('W1,W2')
            Wx = W1+((W2-W1)*((X/L)**α))
            F = integrate(U.T*Wx,(X,(length_of_element[index]),length_of_element_1[index]))
            F1 = sym.lambdify([L,W1,W2],F)
            F2=F1(J,WW1,WW2)
            print("Force matrix F{0}= \n".format(index+1))
            pprint(F,use_unicode=false,wrap_line=false)
            List_of_Force_vector.append(F2)
        elif Force in ('P,p'):
            W1 = float(input('\nGive value of point load = '))
            dist=float(input('\nGive value of distance from left end within member = '))
            W = sym.symbols('W')
            F = integrate(U.T*W,(X,(length_of_element[index]),length_of_element_1[index]))
            F1 = sym.lambdify([L,W],F)
            F2=F1(dist,W1)
            print("Force matrix F{0}= \n".format(index+1))
            pprint(F,use_unicode=false,wrap_line=false)
            List_of_Force_vector.append(F2)

    def twonodebar_jointstiffnessmatrix():
        global Joint_stiffnessmatrix
        Joint_stiffnessmatrix=np.zeros((NO_OF_ELEMENT+1,NO_OF_ELEMENT+1))
        initial_index= 0
        size_1 = 1
        size_2 = 2
        index_i = 1
        index1_step = 2
        index_j = 3
        for index  in range(0, len(List_of_Stiffness_matrices)):
            matrix = List_of_Stiffness_matrices[index]
            print("\n K{0}{1} = {2}".format(index+1,index+1,matrix))
            index_next_matrices = index+1
            if index != len(List_of_Stiffness_matrices)-1:
                next_matrices = List_of_Stiffness_matrices[index_next_matrices]
            if index == 0 :
                Joint_stiffnessmatrix[initial_index:size_1,initial_index:size_2] = matrix[initial_index:size_1,initial_index:size_2]
                Joint_stiffnessmatrix[size_1:size_2,initial_index:size_1] = matrix[size_1:size_2,initial_index:size_1]
                if index != len(List_of_Stiffness_matrices)-1:
                    Joint_stiffnessmatrix[size_1:size_2,size_1:size_2] = matrix[size_1:size_2,size_1:size_2]+next_matrices[initial_index:size_1,initial_index:size_1]
                else:
                    Joint_stiffnessmatrix[size_1:size_2,size_1:size_2] = matrix[size_1:size_2,size_1:size_2]
            else:
                Joint_stiffnessmatrix[index_i:index1_step,index1_step:index_j] = matrix[initial_index:size_1,size_1:size_2]
                Joint_stiffnessmatrix[index1_step:index_j,index_i:index1_step] = matrix[size_1:size_2,initial_index:size_1]
                if index != len(List_of_Stiffness_matrices)-1:
                    Joint_stiffnessmatrix[index1_step:index_j,index1_step:index_j] = matrix[size_1:size_2,size_1:size_2]+next_matrices[initial_index:size_1,initial_index:size_1]
                else:
                    Joint_stiffnessmatrix[index1_step:index_j,index1_step:index_j] = matrix[size_1:size_2,size_1:size_2]
                index_i +=1
                index1_step += 1
                index_j += 1
        print("\nSj=\n{0}".format(np.around(Joint_stiffnessmatrix,decimals=3)))
    def joint_force_vector():
        global Joint_Force_vector
        Joint_Force_vector=np.zeros((NO_OF_ELEMENT+1,1))
        initial_index= 0
        size_1 = 1
        size_2 = 2
        index1_step = 2
        index_j = 3
        for index  in range(0, len(List_of_Force_vector)):
            matrix = List_of_Force_vector[index]
            print("\n F{0}{1} = {2}".format(index+1,index+1,matrix))
            index_next_matrices = index+1
            if index != len(List_of_Force_vector)-1:
                next_matrices = List_of_Force_vector[index_next_matrices]
            if index == 0 :
                Joint_Force_vector[initial_index:size_1,initial_index:size_1] = matrix[initial_index:size_1,initial_index:size_1]
                if index != len(List_of_Force_vector)-1:
                    Joint_Force_vector[size_1:size_2,initial_index:size_1] = matrix[size_1:size_2,initial_index:size_1]+next_matrices[initial_index:size_1,initial_index:size_1]
                else:
                    Joint_Force_vector[size_1:size_2,initial_index:size_1] = matrix[size_1:size_2,initial_index:size_1]
            else:
                if index != len(List_of_Force_vector)-1:
                    Joint_Force_vector[index1_step:index_j,initial_index:size_1] = matrix[size_1:size_2,initial_index:size_1]+next_matrices[initial_index:size_1,initial_index:size_1]
                else:
                    Joint_Force_vector[index1_step:index_j,initial_index:size_1] = matrix[size_1:size_2,initial_index:size_1]
                index1_step += 1
                index_j += 1
        print("\nF=\n{0}".format(np.around(Joint_Force_vector,decimals=3)))
    def displacement_vector_for_2node():
        global Final_displacement_matrix,Final_displacement_matrix1
        print("\nPlease give following details to calculate stresses at nodes\n")
        r = input('\nWrite "BOTH" for both end fixed\n Write "SINGLE" fore any one end fixed\n')
        if r in ('BOTH,both'):
            displacement_matrix_row=np.delete(Joint_stiffnessmatrix,0,0)
            displacement_matrix_column=np.delete(displacement_matrix_row,0,1)
            displacement_matrix_row1=np.delete(displacement_matrix_column,-1,0)
            displacement_matrix_column1=np.delete(displacement_matrix_row1,-1,1)
            force_row_first=np.delete(Joint_Force_vector,0,0)
            force_row_last=np.delete(force_row_first,-1,0)
            m=np.linalg.inv(displacement_matrix_column1)
            Final_displacement_matrix= np.array(m).dot(force_row_last)
            Final_displacement_matrix2=np.insert(Final_displacement_matrix,0,0,axis=0) 
            Final_displacement_matrix1=np.insert(Final_displacement_matrix2,len(Final_displacement_matrix2),0,axis=0)
            # print("K = {0}".format(displacement_matrix_column1))
            # print("f = {0}".format(force_row_last))
            print("\nDisplacement vector=\n{0}".format(np.around(Final_displacement_matrix1,decimals=3)))
        elif r in ('SINGLE,single'):
            s = input('\nGive at which end support is fixed as asked below\nWrite "LEFT" for left end fixed support\nWrite "RIGHT" for right end fixed support\n')
            if s in ('LEFT,left'):
                displacement_matrix_row=np.delete(Joint_stiffnessmatrix,0,0)
                displacement_matrix_column=np.delete(displacement_matrix_row,0,1)
                force_row_first=np.delete(Joint_Force_vector,0,0)
                m=np.linalg.inv(displacement_matrix_column)
                Final_displacement_matrix= np.array(m).dot(force_row_first)
                Final_displacement_matrix1=np.insert(Final_displacement_matrix,0,0,axis=0)
                # print("K = {0}".format(displacement_matrix_column))
                # print("f = {0}".format(force_row_first))
                print("Displacement vector=\n{0}".format(np.around(Final_displacement_matrix1,decimals=3)))
            elif s in ('RIGHT,right'):
                displacement_matrix_row1=np.delete(Joint_stiffnessmatrix,-1,0)
                displacement_matrix_column1=np.delete(displacement_matrix_row1,-1,1)
                force_row_last=np.delete(Joint_Force_vector,-1,0)
                m=np.linalg.inv(displacement_matrix_column1)
                Final_displacement_matrix= np.array(m).dot(force_row_last)
                Final_displacement_matrix1=np.insert(Final_displacement_matrix,len(Final_displacement_matrix),0,axis=0)
                # print("K = {0}".format(displacement_matrix_column1))
                # print("f = {0}".format(force_row_last))
                print("\nDisplacement Vector=\n{0}".format(np.around(Final_displacement_matrix1,decimals=3)))
        print("\nPlease give following details if there is lack of fit at any nodes")
        h = input('\nWrite "BOTH" for both end having lack of fit with fixed ends\n Write "SINGLE" fore any one end having lack of fit\nWrite "NONE" for no ends having lack of fit')
        if h in ('BOTH,both'):
            d = float(input('\nGive lack of fit Value at left end = '))
            q = float(input('\nGive lack of fit Value at right end = '))
            Final_displacement_matrix_lack_of_fit=np.insert(Final_displacement_matrix,0,d,axis=0)
            Final_displacement_matrix_lack_of_fit_final=np.insert(Final_displacement_matrix_lack_of_fit,len(Final_displacement_matrix_lack_of_fit),q,axis=0)
            print("\nDisplacement Vector in case of lack of fit=\n{0}".format(np.around(Final_displacement_matrix_lack_of_fit_final,decimals=3)))
        elif h in ('SINGLE,single'):
            hh = input('\nGive at which end support is fixed as asked below\nWrite "LEFT" for left end fixed support\nWrite "RIGHT" for right end fixed support')
            if hh in ('LEFT,left'):
                d = float(input('\nGive lack of fit Value at left end = '))
                Final_displacement_matrix_lack_of_fit=np.insert(Final_displacement_matrix,0,d,axis=0)
                print("\nDisplacement Vector in case of lack of fit=\n{0}".format(np.around(Final_displacement_matrix_lack_of_fit,decimals=3)))
            elif hh in ('RIGHT,right'):
                q = float(input('\nGive lack of fit Value at right end = '))
                Final_displacement_matrix_lack_of_fit_final=np.insert(Final_displacement_matrix,len(Final_displacement_matrix),q,axis=0)
                print("\nDisplacement Vector in case of lack of fit=\n{0}".format(np.around(Final_displacement_matrix_lack_of_fit_final,decimals=3)))
        elif h in ('NONE,none'):
            print("\n")
    def Stresses_for_2node():
        give_member_number=int(input('\nPlease Give member number, where stress value is required\n'))
        Elasticity= float(input('\nGive value of elasticity = '))
        t = float(input('\nGive total value from left node\nGive value at which stress is required within member = '))
        if give_member_number>NO_OF_ELEMENT:
            print('\nGiven member number is greater than the given number of element')
        else:
            B2 = sym.lambdify([X,L],B)
            length = length_of_int_element[0:give_member_number]
            suma=sum(length)
            BB=B2(t,suma)
            # print("F2={0}".format(BB))
            c = Final_displacement_matrix1[give_member_number-1:give_member_number+1]
            # print("c={0}".format(c))
            stress = np.array(BB).dot(c).dot(Elasticity)
            print("Stress=\n{0}".format(np.around(stress,decimals=3)))
    twonodebar_jointstiffnessmatrix()
    joint_force_vector()
    displacement_vector_for_2node()
    Stresses_for_2node()
#########################################################################################################################################################

    #####################################################################################################################################################               
            
    
##########################################################################################################################################################
if D==3:
    for index in range (0,NO_OF_ELEMENT):
        print("****************for member = {0}*************\n\n ".format(index+1))
        J=float(input('Give length of member {0} = '.format(index+1)))
        length_of_int_element.append(J)
        length_of_element.append(index*L)
        length_of_element_1.append((index+1)*L)
        U_x=sym.Matrix([[1,X,X**2]])
        U_xx=sym.Matrix([[1,0,0],[1,L/2,L**2/4],[1,L,L*L]])
        In = U_xx.inv()
        U=U_x*In
        print("Shape Function for member {0}= \n".format(index+1))
        pprint(U, use_unicode=false, wrap_line=false)
        B=derive_by_array(U,X)
        print("B matrix for member {0}= \n".format(index+1))
        pprint(B, use_unicode=false, wrap_line=false)
        print("\nWrite 'E' if elasticity is not same throughout the bar\nWrite 'A' if area is not same throughout the bar\nIf all three properties are same throughout the bar then write 'SAP'")
        D=input('Write name of property which has not same value along the member {0} as given in above list\n'.format(index+1))
        if D in ('E','e'):
            α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
            EE1=float(input('give Value of Elasticity at left end (E1) in N/mm^2\n'))
            EE2=float(input('give Value of Elasticity at right end (E2) in N/mm^2\n')) 
            AA =float(input('give Value of area in mm^2\n'))
            E1,E2,A = sym.symbols('E1,E2,A')
            Ex = E1+((E2-E1)*((X/L)))
            j=sym.Matrix([B])
            N=j.T*j*Ex*A
            K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
            K1=sym.lambdify([L,E1,E2,A],K)
            K2=K1(J,EE1,EE2,AA)
            List_of_Stiffness_matrices.append(K2)
            #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
            print("Stiffness matrix K{0} = \n".format(index+1))
            pprint(K, use_unicode=false, wrap_line=false)
        if D in ('A','a'):
            α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
            AA1=float(input('give Value of area at left end (A1) = \n'))
            AA2=float(input('give Value of area at right end (A2) =\n')) 
            EE = float(input('give Value of Elasticity E =\n'))
            A1,A2,E = sym.symbols('A1,A2,E')
            Ax = A1+((A2-A1)*((X/L)**α))
            j=sym.Matrix([B])
            # print("BtB = ***********\n")
            # BtB=j.T*j
            # print("N = ***********\n")
            N=j.T*j*E*Ax
            # pprint(BtB, use_unicode=false, wrap_line=false)
            # pprint(N, use_unicode=false, wrap_line=false)
            K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
            K1=sym.lambdify([L,A1,A2,E],K)
            K2=K1(J,AA1,AA2,EE)
            List_of_Stiffness_matrices.append(K2)
            #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
            print("Stiffness matrix K{0} = \n".format(index+1))
            pprint(K, use_unicode=false, wrap_line=false)
            print("Stiffness matrix K{0} = \n{1}".format(index+1,K))
        if D in ('SAP','sap'):
            EEE = float(input('give Value of Elasticity E =\n'))
            AAA = float(input('give Value of area in mm^2\n'))
            A,E = sym.symbols('A,E')
            j=sym.Matrix([B])
            N=j.T*j*E*A
            K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
            K1=sym.lambdify([L,A,E],K)
            K2=K1(J,AAA,EEE)
            List_of_Stiffness_matrices.append(K2)
            #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
            print("Stiffness matrix K{0} = \n".format(index+1))
            pprint(K, use_unicode=false, wrap_line=false)
        if D in ('T','t'):
            print("\nWrite 'E' if elasticity is not same throughout the bar\nWrite 'A' if area is not same throughout the bar\nIf all three properties are same throughout the bar then write 'SAP'\n")
            Z=input('Write name of property which has not same value along the member {0} as given in above list\n'.format(index+1))
            ΔΔT=float(input('Temprature variation '))
            α_00 = float(input('Temprature Co-efficient α_0 =  '))
            if Z in ('E','e'):
                α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
                EE1=float(input('give Value of Elasticity at left end (E1) in N/mm^2\n'))
                EE2=float(input('give Value of Elasticity at right end (E2) in N/mm^2\n')) 
                AA =float(input('give Value of area in mm^2\n'))
                E1,E2,A = sym.symbols('E1,E2,A')
                Ex = E1+((E2-E1)*((X/L)))
                j=sym.Matrix([B])
                N=j.T*j*Ex*A
                K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
                K1=sym.lambdify([L,E1,E2,A],K)
                K2=K1(J,EE1,EE2,AA)
                List_of_Stiffness_matrices.append(K2)
                #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
                print("Stiffness matrix K{0} = \n".format(index+1))
                pprint(K, use_unicode=false, wrap_line=false)
                ΔT,α_0 = sym.symbols('ΔT,α_0')
                F = integrate(j.T*A*Ex*ΔT*α_0,(X,(length_of_element[index]),length_of_element_1[index]))
                F1 = sym.lambdify([L,E1,E2,A,ΔT,α_0],F)
                F2 = F1(J,EE1,EE2,AA,ΔΔT,α_00)
                print("Force matrix F{0}= \n".format(index+1))
                pprint(F,use_unicode=false,wrap_line=false)
                List_of_Force_vector.append(F2)
            if Z in ('A','a'):
                α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
                AA1=float(input('give Value of area at left end (A1) = \n'))
                AA2=float(input('give Value of area at right end (A2) =\n')) 
                EE = float(input('give Value of Elasticity E =\n'))
                A1,A2,E = sym.symbols('A1,A2,E')
                Ax = A1+((A2-A1)*((X/L)**α))
                j=sym.Matrix([B])
                N=j.T*j*E*Ax
                K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
                K1=sym.lambdify([L,A1,A2,E],K)
                K2=K1(J,AA1,AA2,EE)
                List_of_Stiffness_matrices.append(K2)
                #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
                print("Stiffness matrix K{0} = \n".format(index+1))
                pprint(K, use_unicode=false, wrap_line=false)
                ΔT,α_0 = sym.symbols('ΔT,α_0')
                F = integrate(j.T*Ax*E*ΔT*α_0,(X,(length_of_element[index]),length_of_element_1[index]))
                F1 = sym.lambdify([L,A1,A2,E,ΔT,α_0],F)
                F2 = F1(J,AA1,AA2,EE,ΔΔT,α_00)
                print("Force matrix F{0}= \n".format(index+1))
                pprint(F,use_unicode=false,wrap_line=false)
                List_of_Force_vector.append(F2)
            if Z in ('SAP','sap'):
                α = int(input('Give value of α, if linearly varying then write "1"\nif Parabolically varying then write "2"\n'))
                EEE = float(input('give Value of Elasticity E =\n'))
                AAA = float(input('give Value of area in mm^2\n'))
                A,E = sym.symbols('A,E')
                j=sym.Matrix([B])
                N=j.T*j*E*A
                K=integrate(N,(X,(length_of_element[index]),length_of_element_1[index]))
                K1=sym.lambdify([L,A,E],K)
                K2=K1(J,AAA,EEE)
                List_of_Stiffness_matrices.append(K2)
                #print("Stiffness matrix K{0} = {1}\n".format(index+1,K2))
                print("Stiffness matrix K{0} = \n".format(index+1))
                pprint(K, use_unicode=false, wrap_line=false)
                ΔT,α_0 = sym.symbols('ΔT,α_0')
                F = integrate(j.T*A*E*ΔT*α_0,(X,(length_of_element[index]),length_of_element_1[index]))
                F1 = sym.lambdify([L,A,E,ΔT,α_0],F)
                F2 = F1(J,AAA,EEE,ΔΔT,α_00)
                print("Force matrix F{0}= \n".format(index+1))
                pprint(F,use_unicode=false,wrap_line=false)
                List_of_Force_vector.append(F2)  
        Force = input('Give value of loading as asked below\nWrite "P" for point load\nWrite "S" for same load throughout the member or for UDL on member throughout length\nWrite "D" for unsymetrical loading\n')
        if Force in ('s,S'):
            W1 = float(input('Give value of load = '))
            W = sym.symbols('W')
            F = integrate(U.T*W,(X,(length_of_element[index]),length_of_element_1[index]))
            F1 = sym.lambdify([L,W],F)
            F2=F1(J,W1)
            print("Force matrix F{0}= \n".format(index+1))
            pprint(F,use_unicode=false,wrap_line=false)
            List_of_Force_vector.append(F2)
        elif Force in ('D,d'):
            WW1 = float(input('Give value of load at left end W1= '))
            WW2 = float(input('Give value of load at right end W2= '))
            W1,W2 = sym.symbols('W1,W2')
            Wx = W1+((W2-W1)*((X/L)**α))
            F = integrate(U.T*Wx,(X,(length_of_element[index]),length_of_element_1[index]))
            F1 = sym.lambdify([L,W1,W2],F)
            F2=F1(J,WW1,WW2)
            print("Force matrix F{0}= \n".format(index+1))
            pprint(F,use_unicode=false,wrap_line=false)
            List_of_Force_vector.append(F2)
        elif Force in ('P,p'):
            W1 = float(input('\nGive value of point load = '))
            dist=float(input('\nGive value of distance from left end within member = '))
            W = sym.symbols('W')
            F = integrate(U.T*W,(X,(length_of_element[index]),length_of_element_1[index]))
            F1 = sym.lambdify([L,W],F)
            F2=F1(dist,W1)
            print("Force matrix F{0}= \n".format(index+1))
            pprint(F,use_unicode=false,wrap_line=false)
            List_of_Force_vector.append(F2)        
    def threenodebar_jointstiffnessmatrix():
        global Joint_stiffnessmatrix
        Joint_stiffnessmatrix=np.zeros((2*NO_OF_ELEMENT+1,2*NO_OF_ELEMENT+1))
        initial_index= 0
        size_1 = 2
        size1=1
        #size_2 = 4
        size_2 = 3
        index_i = 2
        index1_step = 4
        index_j = 6
        for index  in range(0, len(List_of_Stiffness_matrices)):
            matrix = List_of_Stiffness_matrices[index]
            print("\n K{0}{1} = {2}".format(index+1,index+1,matrix))
            index_next_matrices = index+1
            if index != len(List_of_Stiffness_matrices)-1:
                next_matrices = List_of_Stiffness_matrices[index_next_matrices]
            if index == 0 :
                Joint_stiffnessmatrix[initial_index:size_1,initial_index:size_2] = matrix[initial_index:size_1,initial_index:size_2]
                Joint_stiffnessmatrix[size_1:size_2,initial_index:size_1] = matrix[size_1:size_2,initial_index:size_1]
                if index != len(List_of_Stiffness_matrices)-1:
                    Joint_stiffnessmatrix[size_1:size_2,size_1:size_2] = matrix[size_1:size_2,size_1:size_2]+next_matrices[initial_index:size1,initial_index:size1]
                else:
                    Joint_stiffnessmatrix[size_1:size_2,size_1:size_2] = matrix[size_1:size_2,size_1:size_2]
            else:
                Joint_stiffnessmatrix[index_i:index1_step,index1_step-1:index_j-1] = matrix[initial_index:size_1,size1:size_2]
                Joint_stiffnessmatrix[index1_step-1:index_j-1,index_i:index1_step-1] = matrix[size1:size_2,initial_index:size1]
                Joint_stiffnessmatrix[index1_step:index_j-1,index_i+1:index1_step] = matrix[size_1:size_2,size1:size_1]
                if index != len(List_of_Stiffness_matrices)-1:
                    Joint_stiffnessmatrix[index1_step:index_j-1,index1_step:index_j-1] = matrix[size_1:size_2,size_1:size_2]+next_matrices[initial_index:size1,initial_index:size1]
                else:
                    Joint_stiffnessmatrix[index1_step:index_j-1,index1_step:index_j-1] = matrix[size_1:size_2,size_1:size_2]
                index_i +=2
                index1_step += 2
                index_j += 2
        print("\nSj=\n{0}".format(Joint_stiffnessmatrix))
    def joint_force_vector():
        global Joint_Force_vector
        Joint_Force_vector=np.zeros((2*NO_OF_ELEMENT+1,1))
        initial_index= 0
        size_1 = 1
        size_2 = 2
        index_i =3
        index1_step = 4
        index_j = 5
        for index  in range(0, len(List_of_Force_vector)):
            matrix = List_of_Force_vector[index]
            print("\n F{0}{1} = {2}".format(index+1,index+1,matrix))
            index_next_matrices = index+1
            if index != len(List_of_Force_vector)-1:
                next_matrices = List_of_Force_vector[index_next_matrices]
            if index == 0 :
                Joint_Force_vector[initial_index:size_2,initial_index:size_1] = matrix[initial_index:size_2,initial_index:size_1]
                if index != len(List_of_Force_vector)-1:
                    Joint_Force_vector[size_2:size_2+1,initial_index:size_1] = matrix[size_2:size_2+1,initial_index:size_1]+next_matrices[initial_index:size_1,initial_index:size_1]
                else:
                    Joint_Force_vector[size_2:size_2+1,initial_index:size_1] = matrix[size_2:size_2+1,initial_index:size_1]
            else:
                Joint_Force_vector[index_i:index1_step,initial_index:size_1] = matrix[size_1:size_2,initial_index:size_1]
                if index != len(List_of_Force_vector)-1:
                    Joint_Force_vector[index1_step:index_j,initial_index:size_1] = matrix[size_2:size_2+1,initial_index:size_1]+next_matrices[initial_index:size_1,initial_index:size_1]
                else:
                    Joint_Force_vector[index1_step:index_j,initial_index:size_1] = matrix[size_2:size_2+1,initial_index:size_1]
                index1_step += 2
                index_j += 2
                index_i +=2
        print("\nF=\n{0}".format(Joint_Force_vector))
    def displacement_vector_for_3node():
        global Final_displacement_matrix,Final_displacement_matrix1
        print("\nPlease give following details to calculate stresses at nodes\n")
        r = input('\nWrite "BOTH" for both end fixed\n Write "SINGLE" fore any one end fixed\n')
        if r in ('BOTH,both'):
            displacement_matrix_row=np.delete(Joint_stiffnessmatrix,0,0)
            displacement_matrix_column=np.delete(displacement_matrix_row,0,1)
            displacement_matrix_row1=np.delete(displacement_matrix_column,-1,0)
            displacement_matrix_column1=np.delete(displacement_matrix_row1,-1,1)
            force_row_first=np.delete(Joint_Force_vector,0,0)
            force_row_last=np.delete(force_row_first,-1,0)
            m=np.linalg.inv(displacement_matrix_column1)
            Final_displacement_matrix= np.array(m).dot(force_row_last)
            Final_displacement_matrix2=np.insert(Final_displacement_matrix,0,0,axis=0) 
            Final_displacement_matrix1=np.insert(Final_displacement_matrix2,len(Final_displacement_matrix2),0,axis=0) 
            # print("K = {0}".format(displacement_matrix_column1))
            # print("f = {0}".format(force_row_last))
            print("\nDisplacement vector=\n{0}".format(Final_displacement_matrix1))
        elif r in ('SINGLE,single'):
            s = input('\nGive at which end support is fixed as asked below\nWrite "LEFT" for left end fixed support\nWrite "RIGHT" for right end fixed support\n')
            if s in ('LEFT,left'):
                displacement_matrix_row=np.delete(Joint_stiffnessmatrix,0,0)
                displacement_matrix_column=np.delete(displacement_matrix_row,0,1)
                force_row_first=np.delete(Joint_Force_vector,0,0)
                m=np.linalg.inv(displacement_matrix_column)
                Final_displacement_matrix= np.array(m).dot(force_row_first)
                Final_displacement_matrix1=np.insert(Final_displacement_matrix,0,0,axis=0) 
                # print("K = {0}".format(displacement_matrix_column))
                # print("f = {0}".format(force_row_first))
                print("Displacement vector=\n{0}".format(Final_displacement_matrix1))
            elif s in ('RIGHT,right'):
                displacement_matrix_row1=np.delete(Joint_stiffnessmatrix,-1,0)
                displacement_matrix_column1=np.delete(displacement_matrix_row1,-1,1)
                force_row_last=np.delete(Joint_Force_vector,-1,0)
                m=np.linalg.inv(displacement_matrix_column1)
                Final_displacement_matrix= np.array(m).dot(force_row_last)
                Final_displacement_matrix1=np.insert(Final_displacement_matrix,len(Final_displacement_matrix),0,axis=0) 
                # print("K = {0}".format(displacement_matrix_column1))
                # print("f = {0}".format(force_row_last))
                print("\nDisplacement Vector=\n{0}".format(Final_displacement_matrix1))
        print("\nPlease give following details if there is lack of fit at any nodes")
        h = input('\nWrite "BOTH" for both end having lack of fit with fixed ends\n Write "SINGLE" fore any one end having lack of fit\nWrite "NONE" for no ends having lack of fit')
        if h in ('BOTH,both'):
            d = float(input('\nGive lack of fit Value at left end = '))
            q = float(input('\nGive lack of fit Value at right end = '))
            Final_displacement_matrix_lack_of_fit=np.insert(Final_displacement_matrix,0,d,axis=0)
            Final_displacement_matrix_lack_of_fit_final=np.insert(Final_displacement_matrix_lack_of_fit,len(Final_displacement_matrix_lack_of_fit),q,axis=0)
            print("\nDisplacement Vector in case of lack of fit=\n{0}".format(Final_displacement_matrix_lack_of_fit_final))
        elif h in ('SINGLE,single'):
            hh = input('\nGive at which end support is fixed as asked below\nWrite "LEFT" for left end fixed support\nWrite "RIGHT" for right end fixed support')
            if hh in ('LEFT,left'):
                d = float(input('\nGive lack of fit Value at left end = '))
                Final_displacement_matrix_lack_of_fit=np.insert(Final_displacement_matrix,0,d,axis=0)
                print("\nDisplacement Vector in case of lack of fit=\n{0}".format(Final_displacement_matrix_lack_of_fit))
            elif hh in ('RIGHT,right'):
                q = float(input('\nGive lack of fit Value at right end = '))
                Final_displacement_matrix_lack_of_fit_final=np.insert(Final_displacement_matrix,len(Final_displacement_matrix),q,axis=0)
                print("\nDisplacement Vector in case of lack of fit=\n{0}".format(Final_displacement_matrix_lack_of_fit_final))
        elif h in ('NONE,none'):
            print("\n")
    def Stresses_for_3node():
        give_member_number=int(input('\nPlease Give member number, where stress value is required\n'))
        Elasticity= float(input('\nGive value of elasticity = '))
        t = float(input('\nGive total value from left node\nGive value at which stress is required within member = '))
        if give_member_number>NO_OF_ELEMENT:
            print('\nGiven member number is greater than the given number of element')
        else:
            B2 = sym.lambdify([X,L],B)
            length = length_of_int_element[0:give_member_number]
            suma=sum(length)
            BB=B2(t,suma)
            #print("F2={0}".format(BB))
            if give_member_number == 1:
                c = Final_displacement_matrix1[0:3]
                stress = np.array(BB).dot(c).dot(Elasticity)
                print("Stress=\n{0}".format(stress))
            else:
                c = Final_displacement_matrix1[2*give_member_number-2:2*give_member_number+1]
                #print("c={0}".format(c))
                stress = np.array(BB).dot(c).dot(Elasticity)
                print("Stress=\n{0}".format(stress))
    threenodebar_jointstiffnessmatrix()
    joint_force_vector()
    displacement_vector_for_3node()
    Stresses_for_3node()
    
##########################################################################################################################################################
 ## developed by: JV
##########################################################################################################################################################
