import warnings

warnings.filterwarnings("ignore")
import numpy as np
from numpy import *
import math


def listsublistsum(input_list1, size_of_innerlist):
    global final_list
    list1 = input_list1
    a = size_of_innerlist
    final_list = []
    for _ in range(0, a):
        final_list.append(0)
    for ele in range(0, len(list1)):
        sum_list = []
        list2 = list1[ele]
        for (item1, item2) in zip(list2, final_list):
            # if ele == len(list1)-1:
            #     fin_list = []
            # fin_list.append(item2+item1)
            # print(">>>>>>>>>>>>>>>>>>>{0}".format(sum_list))
            sum_list.append(item1 + item2)
        final_list = sum_list


print(
    "Write 'E' if all the members having same elasticity\nWrite 'I' if all the members have same C/S\nWrite 'L' if all the mrmbers have same length\nIf all three properties are same for all the members then write 'SAP'\nIf not a single property is same then write 'D'")
No_of_element = int(input("Please enter the number of elements for the beam: "))
A = input('Write name of property which has same value for all elemnt as given in above list\n')
length_list = []
if A in ("E", "e"):
    E = float(input("Please give the elasticity of elements for the beam: "))
elif A in ("I", "i"):
    I = float(input("Please give the inertia of elements for the beam: "))
elif A in ("L", "l"):
    L = float(input("Please give the Length of elements for the beam: "))
    length_list.append(L)
elif A in ("SAP", "sap"):
    E = float(input("Please give the elasticity of elements for the beam: "))
    I = float(input("Please give the inertia of elements for the beam: "))
    L = float(input("Please give the Length of elements for the beam: "))
    length_list.append(L)
No_of_nodes = No_of_element + 1
List_of_matrices = []

for index in range(0, No_of_element):
    if A in ("E", "e"):
        I = float(input("Please give the inertia of elements for the beam: "))
        L = float(input("Please give the Length of elements for the beam: "))
        k = np.around(np.array(
            [[(12 * E * I) / (L * L * L), (6 * E * I) / (L * L), (-12 * E * I) / (L * L * L), (6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (4 * E * I) / (L), (-6 * E * I) / (L * L), (2 * E * I) / (L)],
             [(-12 * E * I) / (L * L * L), (-6 * E * I) / (L * L), (12 * E * I) / (L * L * L), (-6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (2 * E * I) / (L), (-6 * E * I) / (L * L), (4 * E * I) / (L)]]), decimals=2)
        length_list.append(L)
    elif A in ("I", "i"):
        E = float(input("Please give the elasticity of elements for the beam: "))
        L = float(input("Please give the Length of elements for the beam: "))
        k = np.around(np.array(
            [[(12 * E * I) / (L * L * L), (6 * E * I) / (L * L), (-12 * E * I) / (L * L * L), (6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (4 * E * I) / (L), (-6 * E * I) / (L * L), (2 * E * I) / (L)],
             [(-12 * E * I) / (L * L * L), (-6 * E * I) / (L * L), (12 * E * I) / (L * L * L), (-6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (2 * E * I) / (L), (-6 * E * I) / (L * L), (4 * E * I) / (L)]]), decimals=2)
        length_list.append(L)
    elif A in ("L", "l"):
        E = float(input("Please give the elasticity of elements for the beam: "))
        I = float(input("Please give the inertia of elements for the beam: "))
        k = np.around(np.array(
            [[(12 * E * I) / (L * L * L), (6 * E * I) / (L * L), (-12 * E * I) / (L * L * L), (6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (4 * E * I) / (L), (-6 * E * I) / (L * L), (2 * E * I) / (L)],
             [(-12 * E * I) / (L * L * L), (-6 * E * I) / (L * L), (12 * E * I) / (L * L * L), (-6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (2 * E * I) / (L), (-6 * E * I) / (L * L), (4 * E * I) / (L)]]), decimals=2)
    elif A in ("SAP", "sap"):
        k = np.around(np.array(
            [[(12 * E * I) / (L * L * L), (6 * E * I) / (L * L), (-12 * E * I) / (L * L * L), (6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (4 * E * I) / (L), (-6 * E * I) / (L * L), (2 * E * I) / (L)],
             [(-12 * E * I) / (L * L * L), (-6 * E * I) / (L * L), (12 * E * I) / (L * L * L), (-6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (2 * E * I) / (L), (-6 * E * I) / (L * L), (4 * E * I) / (L)]]), decimals=2)
        length_list.append(L)
    elif A in ("d", "D"):
        E = float(input("Please give the elasticity of elements for the beam: "))
        I = float(input("Please give the inertia of elements for the beam: "))
        L = float(input("Please give the Length of elements for the beam: "))
        k = np.around(np.array(
            [[(12 * E * I) / (L * L * L), (6 * E * I) / (L * L), (-12 * E * I) / (L * L * L), (6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (4 * E * I) / (L), (-6 * E * I) / (L * L), (2 * E * I) / (L)],
             [(-12 * E * I) / (L * L * L), (-6 * E * I) / (L * L), (12 * E * I) / (L * L * L), (-6 * E * I) / (L * L)],
             [(6 * E * I) / (L * L), (2 * E * I) / (L), (-6 * E * I) / (L * L), (4 * E * I) / (L)]]), decimals=2)
        length_list.append(L)
    List_of_matrices.append(k)
Joint_stiffnessmatrix = np.around(np.zeros((2 * No_of_nodes, 2 * No_of_nodes)), decimals=3)
initial_index = 0
size_1 = 2
size_2 = 4
index_i = 2
index1_step = 4
index_j = 6
for index in range(0, len(List_of_matrices)):
    matrix = List_of_matrices[index]
    print("\n K{0}{1} = {2}".format(index + 1, index + 1, matrix))
    index_next_matrices = index + 1
    if index != len(List_of_matrices) - 1:
        next_matrices = List_of_matrices[index_next_matrices]
    if index == 0:
        Joint_stiffnessmatrix[initial_index:size_1, initial_index:size_2] = matrix[initial_index:size_1,
                                                                            initial_index:size_2]
        Joint_stiffnessmatrix[size_1:size_2, initial_index:size_1] = matrix[size_1:size_2, initial_index:size_1]
        if index != len(List_of_matrices) - 1:
            Joint_stiffnessmatrix[size_1:size_2, size_1:size_2] = matrix[size_1:size_2, size_1:size_2] + next_matrices[
                                                                                                         initial_index:size_1,
                                                                                                         initial_index:size_1]
        else:
            Joint_stiffnessmatrix[size_1:size_2, size_1:size_2] = matrix[size_1:size_2, size_1:size_2]
    else:
        Joint_stiffnessmatrix[index_i:index1_step, index1_step:index_j] = matrix[initial_index:size_1, size_1:size_2]
        Joint_stiffnessmatrix[index1_step:index_j, index_i:index1_step] = matrix[size_1:size_2, initial_index:size_1]
        if index != len(List_of_matrices) - 1:
            Joint_stiffnessmatrix[index1_step:index_j, index1_step:index_j] = matrix[size_1:size_2,
                                                                              size_1:size_2] + next_matrices[
                                                                                               initial_index:size_1,
                                                                                               initial_index:size_1]
        else:
            Joint_stiffnessmatrix[index1_step:index_j, index1_step:index_j] = matrix[size_1:size_2, size_1:size_2]
        index_i += 2
        index1_step += 2
        index_j += 2
print("\nSj=\n{0}".format(Joint_stiffnessmatrix))
print(
    "Write type of support and location of same as given in below detail\n give type of support as per sequence from "
    "left to right like \nFor Fixed Support type-'F'\nFor Roller Support type-'R'\nFor Hinge or Pin Support type-'H'\n"
    "The program will ask support types for No.of Elements+1time")
index = 0
List_of_type_of_support = []
counter_r = 0
counter_h = 0
counter_ns = 0
for index in range(0, No_of_nodes):
    B = input('write support type as describe above\n')
    # B = support_list[index]
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
#########################################################################################
print("\n\n SFF\n:{0}".format(SFF))
G=int(input('Give total number of loads, Which should not be greater than number of member\n'))

Member_end_actions = []
Member_end_actions_inter = []
# #print("*************{0}".format(length_list))
for idx in range(0, No_of_element):
    Member_NO, NO_Of_Lcases = map(int, input('\n(member number,Number of load case)=').split(","))
    for index in range(0, NO_Of_Lcases):
        print(" \nWrite P for point load \nWrite 'UDL' for uniformly distributed load\nWrite 'M' for moment")
        C = input('\n Give types of loading = ')
        if NO_Of_Lcases > 1:
            if C in ("P", "p"):
                W, a = map(float, input('\n(Load value,Distance from left node)=').split(","))
                AML = np.array([W * (length_list[idx] - a) / length_list[idx],
                                W * a * (length_list[idx] - a) * (length_list[idx] - a) / length_list[idx] /
                                length_list[idx], W * a / length_list[idx],
                                -W * (length_list[idx] - a) * a * a / length_list[idx] / length_list[idx]])
                Member_end_actions_inter.append(AML)
                # Member_end_actions.append(Member_end_actions_inter)
            elif C in ("UDL", "udl"):
                W = float(input('Give value of loading sequence wise from left to right in KN/m \n'))
                AML = np.array([W * length_list[idx] / 2, W * length_list[idx] * W * length_list[idx] / 12,
                                W * length_list[idx] / 2, -W * length_list[idx] * length_list[idx] / 12])
                Member_end_actions_inter.append(AML)
                # Member_end_actions.append(Member_end_actions_inter)
            elif C in ("M", "m"):
                W, a = map(float, input('\n(Moment value,Distance from left node)=').split(","))
                AML = np.array([W * a * (length_list[idx] - a) / length_list[idx] / length_list[idx],
                                -W * a * a * (length_list[idx] - a) * (length_list[idx] - a) / length_list[idx] /
                                length_list[idx] / length_list[idx],
                                W * a * (length_list[idx] - a) / length_list[idx] / length_list[idx],
                                -W * a * (length_list[idx] - a) * (length_list[idx] - a) / length_list[idx] /
                                length_list[idx] / length_list[idx]])
                Member_end_actions_inter.append(AML)
        listsublistsum(Member_end_actions_inter, 4)
        print("\n mEMBER_INTER:  {0}".format(Member_end_actions_inter))
        if NO_Of_Lcases == 1:
            if C in ("P", "p"):
                W, a = map(float, input('\n(Load value,Distance from left node)=').split(","))
                AML = np.array([W * (length_list[idx] - a) / length_list[idx],
                                W * a * (length_list[idx] - a) * (length_list[idx] - a) / length_list[idx] /
                                length_list[idx], W * a / length_list[idx],
                                -W * (length_list[idx] - a) * a * a / length_list[idx] / length_list[idx]])
                Member_end_actions.append(AML)
            elif C in ("UDL", "udl"):
                W = float(input('Give value of loading sequence wise from left to right in KN/m \n'))
                AML = np.array([W * length_list[idx] / 2, W * length_list[idx] * W * length_list[idx] / 12,
                                W * length_list[idx] / 2, -W * length_list[idx] * length_list[idx] / 12])
                Member_end_actions.append(AML)
            elif C in ("M", "m"):
                W, a = map(float, input('\n(Moment value,Distance from left node)=').split(","))
                AML = np.array([W * a * (length_list[idx] - a) / length_list[idx] / length_list[idx],
                                -W * a * a * (length_list[idx] - a) * (length_list[idx] - a) / length_list[idx] /
                                length_list[idx] / length_list[idx],
                                W * a * (length_list[idx] - a) / length_list[idx] / length_list[idx],
                                -W * a * (length_list[idx] - a) * (length_list[idx] - a) / length_list[idx] /
                                length_list[idx] / length_list[idx]])
                Member_end_actions.append(AML)

    Member_end_actions.append(final_list)

print("\n member end matrices for element:  {0}".format(Member_end_actions))
#     Final_AML = np.zeros((4,1))
#     Final_AML2 = np.zeros((4,1))
#     Final_AML[0:4,0:1] = AML_list[index6:index7,0:1]
#     Member_end_actions.append(AML)
# Final_AML = np.zeros((4,1))
# initial_index= 0
# size_1 = 1
# size_2 = 4
# indexl = 0
# indexk = 4
# for index  in range(0, len(Member_end_actions)):
#     matrix = Member_end_actions[index]
#     print("\n AML{0}{1} = {2}".format(index+1,index+1,memberendactions))
#     index_next_matrices = index+1
#     if index != len(Member_end_actions)-1:
#         next_matrices = Member_end_actions[index_next_matrices]
#     if index == 0 :
#         Final_AML[initial_index:size_2,initial_index:size_1] = Member_end_actions[initial_index:size_2,initial_index:size_1]
#     else:
#         Final_AML[initial_index:size_2,initial_index:size_1] = Member_end_actions[initial_index:size_2,initial_index:size_1]+next_matrices[indexl:indexk,initial_index:size_1]
#     indexl +=4
#     indexk += 4
# print("\nFinal Aml=\n{0}".format(Final_AML))
