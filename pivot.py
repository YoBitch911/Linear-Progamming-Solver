#!/usr/bin/env python
#Single iteration og pivoting for a standard Linear Program
import sys
from threading import Thread
import os
import csv
import numpy
from time import gmtime, strftime

def read_values():
	fp = open('part1.dict','r')
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
	tup = (99999999999,-1)
	#tup[0] is the entering variable
	#tup[1] is the leaving variable
	f = open('part4.ouput','w')
	for key in objective_dict:
		if(key!=-1):
			if(objective_dict[key]>0):
				#possible candidate f
				# find the lowest possible indexed leaving variable for this entering variable, -1 for unbounded,
				leaving_var = find_leaving_variable(key,pivot_dict,m);
				if(leaving_var == -1):
					f.write("UNBOUNDED")
					break;
				if(key < tup[0]):
					tup = (key,leaving_var)
	if(leaving_var == -1):
		return
	#now since we know the entering /leaving var ,find the objective function
	val = pivot_dict[tup[1]][-1]/pivot_dict[tup[1]][tup[0]]*-1.0
	ob_val = objective_dict[-1] + val*objective_dict[tup[0]]
	f.write('%s\n'%tup[0])
	f.write('%s\n'%tup[1])
	f.write('%.4f\n'%ob_val)
	
def find_leaving_variable(var,pivot_dict,m):
	#check if unbounded
	count = 0
	#find the lowest valued ,lowest indiced variable
	min_tup = (9999999999.0,-1);
	for key in pivot_dict:
		if(pivot_dict[key][var] >= 0):
			count += 1
		else:
			if(pivot_dict[key][-1]/pivot_dict[key][var]*-1.0 < min_tup[0]):
				min_tup = (pivot_dict[key][-1]/pivot_dict[key][var]*-1.0,int(key))
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

pivot_dict,objective_dict = read_values()
perform_pivot(pivot_dict,objective_dict)
