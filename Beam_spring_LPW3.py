import sympy as sym
from sympy import *

sym.init_printing()
import math
import numpy as np

member_num = 2
l = 2500
temp_list = 0
print("\n**********************************************************\n\nTaken data for Group 2\nlength = 5000mm\n"
      "Elasticity = 205000N/mm^2\nI=300000000 mm^4\nr=2\nalpha=1\nW=100N/mm"
      "\nspring constant Ks = {0}\n************\n".format(temp_list))
No_of_nodes = member_num + 1
support_list = ['F', 'NS', 'F']
# support_list = ['F', 'F', 'NS', 'F', 'F']
# support_list = ['F', 'F','F','NS','F', 'F', 'F']
# support_list = ['F', 'F', 'F','F','NS','F','F', 'F', 'F']
List_of_matrices = []
List_of_force_atrices = []
for index in range(0, member_num):
    length_list = (index + 1) * l
    x = sym.symbols('x')
    U1 = sym.Matrix([[1, x, x * x, x * x * x]])
    inv = sym.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [1, l, l ** 2, l ** 3], [0, 1, 2 * l, 3 * l ** 2]])
    inversw = inv.inv()
    U = U1 * inversw
    print("shape function N{0}".format(index + 1))
    pprint(simplify(U.T), use_unicode=false, wrap_line=false)
    B1 = derive_by_array(U, x)
    B2 = derive_by_array(B1, x)
    j = sym.Matrix([B2])
    print("Bmatrix B{0}".format(index + 1))
    pprint(simplify(j.T), use_unicode=false, wrap_line=false)
    E = 205000
    I = 300000000
    alpha = 1
    Ix = I * (1 + (2 * ((x / l) ** alpha)))
    # btb=j.T*j
    # print("\n B^T*B = \n")
    # pprint(btb, use_unicode=false, wrap_line=false)
    N1 = j.T * j * E * Ix
    Kb = integrate(N1, (x, 0, length_list))
    spring_matrix = []
    Final_matrix = []
    inverse_matrix = []
    # temp_list = [0, 10, 50, 100, 150, 200]
    N2 = U.T * temp_list * U
    Ks = integrate(N2, (x, 0, length_list))
    K = np.add(Kb, Ks)
    List_of_matrices.append(K)
    # for index in range (0,len(temp_list)):
    #     print("Stiffness matrix Ks for spring value {0} = ".format(index+1))
    #     pprint(spring_matrix[index], use_unicode=false, wrap_line=false,mat_symbol_style="plain")
    #     print("Final Stiffness matrix for spring value {0}".format(index + 1))
    #     pprint(Final_matrix[index], use_unicode=false, wrap_line=false, mat_symbol_style="plain")
    W = 100
    inte = U.T * W
    F = integrate(U.T * W, (x, 0, length_list))
    List_of_force_atrices.append(F)
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
        Joint_stiffnessmatrix[initial_index:size_1, initial_index:size_2] = matrix[initial_index:size_1,
                                                                            initial_index:size_2]
        Joint_stiffnessmatrix[size_1:size_2, initial_index:size_1] = matrix[size_1:size_2, initial_index:size_1]
        if index != len(List_of_matrices) - 1:
            Joint_stiffnessmatrix[size_1:size_2, size_1:size_2] = matrix[size_1:size_2,
                                                                  size_1:size_2] + next_matrices[
                                                                                   initial_index:size_1,
                                                                                   initial_index:size_1]
        else:
            Joint_stiffnessmatrix[size_1:size_2, size_1:size_2] = matrix[size_1:size_2, size_1:size_2]
    else:
        Joint_stiffnessmatrix[index_i:index1_step, index1_step:index_j] = matrix[initial_index:size_1,
                                                                          size_1:size_2]
        Joint_stiffnessmatrix[index1_step:index_j, index_i:index1_step] = matrix[size_1:size_2,
                                                                          initial_index:size_1]
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
Joint_forcematrix = np.around(np.zeros((2 * No_of_nodes, 1)), decimals=3)
initial_index = 0
size_1 = 2
size_2 = 1
size_22 = 4
index_i = 2
index1_step = 4
index_j = 6
for index in range(0, len(List_of_force_atrices)):
    matrix = List_of_force_atrices[index]
    print("\n F{0} = {1}".format(index + 1, matrix))
    index_next_matrices = index + 1
    if index != len(List_of_force_atrices) - 1:
        next_matrices = List_of_force_atrices[index_next_matrices]
    if index == 0:
        Joint_forcematrix[initial_index:size_1, initial_index:size_2] = matrix[initial_index:size_1,
                                                                        initial_index:size_2]
        if index != len(List_of_force_atrices) - 1:
            Joint_forcematrix[size_1:size_22, initial_index:size_2] = matrix[size_1:size_22,
                                                                      initial_index:size_2] + next_matrices[
                                                                                              initial_index:size_1,
                                                                                              initial_index:size_2]
        else:
            Joint_forcematrix[size_1:size_22, initial_index:size_2] = matrix[size_1:size_22, initial_index:size_2]
    else:
        if index != len(List_of_force_atrices) - 1:
            Joint_forcematrix[index1_step:index_j, initial_index:size_2] = matrix[size_1:size_22,
                                                                           initial_index:size_2] + next_matrices[
                                                                                                   initial_index:size_1,
                                                                                                   initial_index:size_2]
        else:
            Joint_forcematrix[index1_step:index_j, initial_index:size_2] = matrix[size_1:size_22, initial_index:size_2]
        index_i += 2
        index1_step += 2
        index_j += 2
print("\njoint force matrix =\n{0}".format(Joint_forcematrix))
Ff = np.zeros(((counter_r + counter_h + (2 * counter_ns)), (1)))
key_pair_columns_force = {}
for index in range(0, len(List_of_type_of_support)):
    user_input = str(List_of_type_of_support[index]).strip().lower()
    if user_input == 'r' or user_input == 'h':
        key_pair_columns_force[index + index + 1] = Joint_forcematrix[index + index + 1, :]
    elif user_input == 'NS' or user_input == 'ns':
        key_pair_columns_force[index + index] = Joint_forcematrix[index + index,:]
        key_pair_columns_force[index + index + 1] = Joint_forcematrix[index + index + 1,:]
temp_list_force = []
for key, column in key_pair_columns_force.items():
    temp_list_force.append(column)
Ff[:, :] = temp_list_force
print("\n\n Ff\n:{0}".format(Ff))
K_inv = np.linalg.inv(SFF)
delta = np.dot(K_inv, Ff)
print("\nCentral slope & deflection for spring value {0} = \n{1}\n".format(temp_list, delta))
