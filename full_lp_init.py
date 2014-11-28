#!/usr/bin/env python
#complete pivoting for a standard Linear Program(given initial dictionary is feasible) +
# initializing using extra positive variable x0,( initialization using dual method to be added later)
import sys
from threading import Thread
import os
import math
import csv
from time import gmtime, strftime

def read_values():
	fp = open('/Users/pgaura1/Desktop/LP/part3TestCases/assignmentParts/part10.dict','r')
	# n is the number of equations
	# m is the number of constraints
	line = fp.readline()
	arr = line.split(' ')
	m = int(arr[0])
	n = int(arr[1])
	# now the basic and non basic integers
	# basic are B1...Bm which are not zero(mostly)
	# non basic are N1..Nn which are zero and reside in the objective function
	line = fp.readline()
	basic_ind = []
	non_basic_ind = []
	arr = line.split(' ')
	for elem in arr:
		if(elem=='' or elem == ' ' or elem == '\n'):
			continue;
		basic_ind.append(int(elem))
	line = fp.readline()
	arr = line.split(' ')
	for elem in arr:
		if(elem=='' or elem == ' ' or elem == '\n'):
			continue;
		non_basic_ind.append(int(elem))
	# now b1, b2..bm which are values for the constraints
	constraint_val  = []
	line = fp.readline()
	arr = line.split(' ')
	for elem in arr:
		if(elem=='' or elem == ' ' or elem == '\n'):
			continue;
		constraint_val.append(float(elem))
	# now get tge coefficient of the constraints
	# m x n matrix
	mat = [[]]
	i = 0
	while i < m:
		#each of the m rows contain n coefficients against the n non basic integers
		line = fp.readline()
		arr = line.split(' ')
		mat.append([])
		for elem in arr:
			if(elem=='' or elem == ' ' or elem == '\n'):
				continue
			mat[i].append(float(elem))
		i+=1;
	# now mat is an m x n matrix where each row consists of coefficients of the constraints for each of the N variables
	# now get the objective row coefficients
	# print the matrix

	line = fp.readline()
	arr = line.split(' ')
	objective_coff = []
	for elem in arr:
		if(elem=='' or elem == ' ' or elem == '\n'):
			continue
		objective_coff.append(float(elem))
	#m = number of constraints or the basic variables
	#n = number of non-basic variables or the variables which appear 
	#basic_ind (m variables)
	#non_basic_ind (n variables)
	#mat = m * n  matrix containing the constraints coefficients (will change with pivoting)
	#constraint_val = m values which are the constants(reals for the m constraints)
	#objective_coff = n coefficients for the n variables in the objective function 
	
	# make  a dictionary of dictionaries where primary key denotes the index of the basic variables in the current dictionary
	# the secondary keys would be n non basic variables along with their coefficients as values
	# an extra key indexed by -1 would denote the current constant/real value associated with that index
	pivot_dict = {}
	i = 0
	for var in basic_ind:
		pivot_dict[var] = {}
		#now add the constants and the variables coefficients associated with it
		pivot_dict[var][-1] = constraint_val[i]
		j = 0
		for nb_var in non_basic_ind:
			pivot_dict[var][nb_var] = mat[i][j]
			j += 1
		i += 1;

	objective_dict = {}
	objective_dict[-1] = objective_coff[0]
	ob_len = len(objective_coff)
	j=0;
	for i in range(1,ob_len):
		objective_dict[non_basic_ind[j]] = objective_coff[i]
		j += 1

	return(pivot_dict,objective_dict)

def perform_pivot(pivot_dict,objective_dict):
	#m = number of contraints
	# n = number of decision variables involved in the objective function

	m = len(pivot_dict)
	n = len(objective_dict) - 1
	# entering variables are in the objective
	# leaving variables are in the constraint
	# now find the candidates for entering variable
	leaving_var = -1;
	tup = [99999999999,-1]
	#tup[0] is the entering variable
	#tup[1] is the leaving variable
	#f = open('part4.ouput','w')

	#check if all the non_basic var in objective dict
	flag=1
	for key in objective_dict:
		if(key==-1):
			continue;
		if(objective_dict[key] > 0):
			flag = 0
			break

	if(flag==1):
		return(2,[objective_dict[-1],0])

	for key in objective_dict:
		if(key!=-1):
			if(objective_dict[key]>0):
				#possible candidate f
				#find the lowest possible indexed leaving variable for this entering variable, -1 for unbounded,
				leaving_var = find_leaving_variable(key,pivot_dict,m);
				if(leaving_var == -1):
					#f.write("UNBOUNDED")
					break;
				if(key < tup[0]):
					tup = [key,leaving_var]
	if(leaving_var == -1):
		return(0,tup)
	#now since we know the entering /leaving var ,find the objective function
	val = pivot_dict[tup[1]][-1]/pivot_dict[tup[1]][tup[0]]*-1.0
	ob_val = objective_dict[-1] + val*objective_dict[tup[0]]
	#print('entering variable is %s\n')%tup[0]
	#print('leaving variable is %s\n')%tup[1]
	#print('current objective value is %s\n')%ob_val
	#f.write('%s\n'%tup[0])
	#f.write('%s\n'%tup[1])
	#f.write('%.4f\n'%ob_val)
	return(1,tup)

def get_new_dict_for_init(pivot_dict,objective_dict):
	#add an extra variable x0 and try to minimize that
	#find the row with the maximum neg value
	basic_var = -1
	max_val = 0
	for key in pivot_dict:
		if(pivot_dict[key][-1] < 0):
			if(abs(pivot_dict[key][-1]) > max_val):
				max_val = abs(pivot_dict[key][-1])
				basic_var = key
	# if all the bi's are positive
	if(basic_var == -1):
		return(0,pivot_dict,objective_dict)
	#get a new dict for x0
	pivot_dict[0] = {}
	pivot_dict[0][basic_var] = 1
	for key in pivot_dict[basic_var]:
		pivot_dict[0][key] = -1*pivot_dict[basic_var][key]
	# now use this to replace X0 in the rest of the dictionary
	# remove basic var

	pivot_dict.pop(basic_var,None)
	# now replace the occurence of X(basic_var) in rest of the program
	for key in pivot_dict:
		if(key==0):
			continue
		for sec_key in pivot_dict[key]:
			pivot_dict[key][sec_key] = pivot_dict[0][sec_key] + pivot_dict[key][sec_key]
		pivot_dict[key][basic_var] = pivot_dict[0][basic_var]
	# update the objective dictionary
	new_objective_dict = {}
	for key in pivot_dict[0]:
		new_objective_dict[key] = -1*pivot_dict[0][key]
		#print("objective_dict[%d]  %d")%(key,new_objective_dict[key])

	#print the pivot_dict
	'''
	for key in pivot_dict:
		print("the current baisc var is %d <-->")%(key)
		for non_basic_var in pivot_dict[key]:
			print("pivot[%d] is %d")%(non_basic_var,pivot_dict[key][non_basic_var])
	'''
	# now return the dict
	# 1 means some extra work (initialization) must be done to retrive the original pivot and objective dict
	return(1,pivot_dict,new_objective_dict)


#tup[0] is the entering variable
#tup[1] is the leaving variable
def change_dict(pivot_dict,objective_dict,tup):
	#pivot_dict is primary indexed by leaving variable(one of the m basic variables)
	#pivot_dict[lv][-1] = constant
	#pivot_dict[lv][nb_var] = some shitty real constant
	#objective_dict is primary indexed by entering variables(one of the n non-basic variables)
	#entering variables(non-basic variables) constitute the objective function and rows of the constraints
	#leaving variables(basic variables) constitute the indices of the constraints
	common_div = abs(pivot_dict[tup[1]][tup[0]])
	#make a new key index by tup[0] which is the entering variable
	pivot_dict[tup[0]] = {}
	for key in pivot_dict[tup[1]]:
		if(int(key)==int(tup[0])):
			continue
		pivot_dict[tup[0]][key] = pivot_dict[tup[1]][key]/common_div;

	pivot_dict[tup[0]][tup[1]] = - 1.0/common_div;
	# remove the tup[1] key from the pivot dict
	pivot_dict.pop(tup[1],None)
	# now change the rest of the keys and remove tup[0] as the sceondary index and add tup[1] and secondary index in rest of the constraints
	for basic_var in pivot_dict:
		if(int(basic_var) == int(tup[0])):
			continue
		for non_basic_var in pivot_dict[basic_var]:
			if(int(non_basic_var) == int(tup[0])):
				continue
			pivot_dict[basic_var][non_basic_var] += pivot_dict[basic_var][tup[0]]*pivot_dict[tup[0]][non_basic_var]
		pivot_dict[basic_var][tup[1]] = pivot_dict[basic_var][tup[0]]*pivot_dict[tup[0]][tup[1]]
		pivot_dict[basic_var].pop(tup[0],None)
	#change objective dict
	#remove tup[0] from objective dict
	#add tup[1] to the objective dict
	for non_basic_var in objective_dict:
		if(non_basic_var == tup[0]):
			continue
		objective_dict[non_basic_var] += objective_dict[tup[0]]*pivot_dict[tup[0]][non_basic_var]
	objective_dict[tup[1]] = objective_dict[tup[0]]*pivot_dict[tup[0]][tup[1]]
	objective_dict.pop(tup[0],None)

	return(pivot_dict,objective_dict)

def find_leaving_variable(var,pivot_dict,m):
	#check if unbounded
	count = 0
	#find the lowest valued ,lowest indiced variable
	min_tup = [9999999999.0,-1];
	for key in pivot_dict:
		if(pivot_dict[key][var] >= 0):
			count += 1
		else:
			if(pivot_dict[key][-1]/pivot_dict[key][var]*-1.0 < min_tup[0]):
				min_tup = [pivot_dict[key][-1]/pivot_dict[key][var]*-1.0,int(key)]
	#unbounded
	if(count == m):
		return(-1)
	#find the lowest indexed key with val = min_tup
	for key in pivot_dict:
		if(pivot_dict[key][var] < 0):
			if(pivot_dict[key][-1]/pivot_dict[key][var]*-1.0 == min_tup[0]):
				if(key < min_tup[1]):
					min_tup[1]= key;
	return(min_tup[1])

def perform_actual_simplex(f,pivot_dict,objective_dict):
	cur_iter = 0
	#f = open('/Users/pgaura1/Desktop/LP/part3TestCases/assignmentParts/ans.dict','w')
	while(1):	
		val,tup = perform_pivot(pivot_dict,objective_dict)
		cur_iter += 1
		if(val == 0):
			print("UNBOUNDED\n");
			f.write('UNBOUNDED\n')
			break;
		if(val == 2):
			print("%.3f\n")%(tup[0])
			f.write('%.3f\n'%tup[0])
			#print("%d\n")%(cur_iter-1)
			#f.write('%d\n'%(cur_iter-1))
			break
		pivot_dict,objective_dict = change_dict(pivot_dict,objective_dict,tup)

def get_org_objective_dict(objective_dict,org_objective_dict,pivot_dict):
	for key in objective_dict:
		objective_dict[key] = 0

	#write the original objective dictionary in terms of the current non-basic variables which are there in the current dictionary
	for key in org_objective_dict:
		if(key in pivot_dict):
			#basic var
			for non_basic_var in pivot_dict[key]:
				objective_dict[non_basic_var] += org_objective_dict[key]*pivot_dict[key][non_basic_var]

		else:
			#non basic var,
			objective_dict[key] += org_objective_dict[key]
	return(objective_dict)


pivot_dict,org_objective_dict = read_values()
objective_dict = {}
for key in org_objective_dict:
	objective_dict[key] = org_objective_dict[key]

#new objective dict would max -x0
val,pivot_dict,objective_dict = get_new_dict_for_init(pivot_dict,objective_dict)
flag = 0
f = open('/Users/pgaura1/Desktop/LP/part3TestCases/ans.dict','w')
if(val==1):
	# solve this linear program
	#print("The problem needs initialization")
	#print("start object value is %f")%(objective_dict[-1])
	cur_iter =0

	while(1):
		val,tup =  perform_pivot(pivot_dict,objective_dict)
		cur_iter += 1
		#val  = 0 means unbounded
		#val == 2 means iteration is done
		if(val == 2 or abs(objective_dict[-1]-0) <  0.0000000000000001):
			if(abs(objective_dict[-1]-0)>0.0000000000000001):
				#print("The problem is infeasible\n")
				flag = 1
			break
		pivot_dict,objective_dict = change_dict(pivot_dict,objective_dict,tup)
	if(flag==0):
		#retrive the original pivot and objective dict
		#find the n non-basic vars(variables which are zero,in case of multiple choose anyone)
		#atleast n+1 variables will be zero , if we remove x0 , then atleast n variables , take any N
		
		#two cases might happen, either x0 is the non-basic variable ,which case we are at a non-degenerate point ,otherwise we are at a degenerate point
		if(0 in pivot_dict):
			#here x0 is the basic variable, at a degenerate feasible point
			#get any random key from the constraint since all will be zero
			for key in pivot_dict[0]:
				if(key == -1):
					continue
				else:
					rand_key = key
					break
			#rand_key is one of the basic variables here
			pivot_dict[rand_key] = {}
			for key in pivot_dict[0]:
				if(key==rand_key):
					continue
				pivot_dict[rand_key][key] = (-1.0*pivot_dict[0][key])/pivot_dict[0][rand_key]

			pivot_dict.pop(0,None)
			# now replace rand_key variable everywhere
			for key in pivot_dict:
				if(key==rand_key):
					continue
				for non_basic_var in pivot_dict[key]:
					if(non_basic_var == rand_key):
						continue
					pivot_dict[key][non_basic_var] += pivot_dict[key][rand_key]*pivot_dict[rand_key][non_basic_var]
				pivot_dict[key].pop(rand_key,None)

			#remove rand_key from the objective dictionnary
			objective_dict.pop(rand_key,None)
			# now get the original objective dictionary
			objective_dict = get_org_objective_dict(objective_dict,org_objective_dict,pivot_dict)

		else:
			# remove the x0 from the contrainst
			for key in pivot_dict:
				pivot_dict[key].pop(0,None)
			# remove the x0 form the constraint
			objective_dict.pop(0,None)
			#simply change the objective dict
			objective_dict = get_org_objective_dict(objective_dict,org_objective_dict,pivot_dict)


		#now our objective and pivot dictionaries are ready
		#now we have the desired objective and pivot dictionaries, so start the main simplex from here
		#make it a function
		perform_actual_simplex(f,pivot_dict,objective_dict)
	else:
		# do nothing as the problem is infeasible
		print("INFEASIBLE")
		f.write("INFEASIBLE\n")
else:
	#simply perform the normal the pivoting scheme
	perform_actual_simplex(f,pivot_dict,objective_dict)

