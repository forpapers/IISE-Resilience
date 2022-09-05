'''
NB:
This code can help you get the objectives under disruptions
Also, you can use the values of variables to get information such as transmission, generation and cost
'''

import os
print(os.getcwd())
import pandas as pd
import numpy as np
from pulp import *
import pulp as pl
from gurobipy import *

First_city = 'Guangdong'
Last_city = 'Yunnan'
df_Param = pd.read_csv('param.csv', index_col=0)
df_Cost = pd.read_csv('cost.csv', index_col=0)
df_Demand = pd.read_csv('demand.csv', index_col=0)
df_OperationHour = pd.read_csv('operation hour.csv', index_col=0)
df_TechCapacity = pd.read_csv('tech capacity.csv', index_col=0)
# 逐行读取csv文件
K = int(df_Param.loc['K', 'A'])
tf = {0:'Guangdong', 1:'Guangxi', 2:'Hainan', 3:'Guizhou', 4:'Yunnan'}
tran_dict = {}
import csv
with open('transmission2-24.csv','r',encoding="utf-8") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row['original'] not in tran_dict.keys(): tran_dict[row['original']] = {}
        if row['destination'] not in tran_dict[row['original']].keys(): tran_dict[row['original']][row['destination']] = []
        for iter in range(K):
            if iter < int(row['Quantity']):
                tran_dict[row['original']][row['destination']].append(float(row['Capacity(GW)']))
for i in range(5):
    if tf[i] not in tran_dict.keys(): tran_dict[tf[i]] = {}
    for j in range(5):
        if tf[j] not in tran_dict[tf[i]].keys(): tran_dict[tf[i]][tf[j]] = []
        list_len = len(tran_dict[tf[i]][tf[j]])
        if list_len != 9:
            for q in range(list_len,9):
                tran_dict[tf[i]][tf[j]].append(0)
ber_h = tran_dict
Cost = df_Cost.loc[First_city:Last_city,'Wind':'Oil'].values
Demand = df_Demand.loc[First_city:Last_city,'GWh'].values
OperationHour = df_OperationHour.loc[First_city:Last_city,'Wind':'Oil'].values
TechCapacity = df_TechCapacity.loc[First_city:Last_city,'Wind':'Oil'].values

ber_x = TechCapacity
alpha = Cost
y = Demand
va = OperationHour
number_J = len(alpha[0,:])
number_R = int(df_Param.loc['number_R', 'A'])
# beta: transmission cost
beta = df_Param.loc['c1':'c5', 'A':'E'].values
v = df_Param.loc['v', 'A']
CO2 = df_Param.loc['CO2', 'A']
penalty_cost = df_Param.loc['penalty_cost', 'A':'E'].values
rate_penalty_cost = df_Param.loc['rate_penalty_cost', 'A']
penalty_cost = penalty_cost*rate_penalty_cost
gamma = penalty_cost
b = df_Param.loc['b', 'A']
fff = df_Param.loc['fff', 'A':'G'].values
rate_fff = df_Param.loc['rate_fff', 'A']
print(fff)

fff /= rate_fff
d = df_Param.loc['c6':'c10', 'A':'E'].values
sigma = df_Param.loc['sigma', 'A']


BIG_PI, ber_theta  = 0,0
for r in range(number_R):
    for r_prime in range(number_R):
        for k in range(K):
            BIG_PI += d[r][r_prime] * ber_h[tf[r]][tf[r_prime]][k]
            ber_theta += ber_h[tf[r]][tf[r_prime]][k]


number_L = 100
buget_relaxation = 0.3
print(buget_relaxation)
BIG_GAMMA = 379720 + buget_relaxation*(439540-379720)
print(BIG_GAMMA)


# Situation 1
import random
AttackerCost_epsilon = np.zeros((number_R,number_R,K,number_L))
num_epsilon = 0
flag= True
while flag:
    tmp_epsilon = np.zeros((number_R, number_R, K))
    for r in range(number_R):
        for r_prime in range(number_R):
            for k in range(K):
                if ber_h[tf[r]][tf[r_prime]][k]!=0:
                    tmp_epsilon[r][r_prime][k] = random.randint(0,1)  #（a <= N <= b)
                    # print(r,r_prime,k,tmp_epsilon[r][r_prime][k])
    for r in range(number_R):
        for r_prime in range(number_R):
            for k in range(K):
                AttackerCost_epsilon[r][r_prime][k][num_epsilon] = tmp_epsilon[r][r_prime][k]
    num_epsilon = num_epsilon+1
    if num_epsilon==number_L:
        flag = False



SAVE_OBJ = np.zeros(number_L)


epsilon = AttackerCost_epsilon
for l in range(number_L):
    prob = LpProblem("problem",LpMinimize)

    var_x = LpVariable.dicts("var_x", (range(number_R),range(number_J)), 0) # Generation output of tech j in region r
    var_s = LpVariable.dicts("var_s", (range(number_R)), 0) # Elec deficit in region r
    var_h = LpVariable.dicts("var_h", (range(number_R),range(number_R),range(K)), 0) # Elec transmission, r*r*k

    # formula(10)
    prob += lpSum(alpha[r,j] * var_x[r][j] for j in range(number_J) for r in range(number_R)) + lpSum(beta[r,r_prime]*var_h[r][r_prime][k] for r_prime in range(number_R) for r in range(number_R) for k in range(K)) + lpSum(gamma[r] * var_s[r] for r in range(number_R))

    # formula(11)
    for r in range(number_R):
        prob += y[r] - var_s[r] <= lpSum(var_x[r][j] for j in range(number_J)) + lpSum(lpSum(var_h[r_prime][r][k] for k in range(K) for r_prime in range(number_R))) - lpSum(lpSum(var_h[r][r_prime][k] for k in range(K) for r_prime in range(number_R)))

    # formula(12)
    prob += lpSum(lpSum(fff[j]*var_x[r][j] for r in range(number_R) for j in range(number_J))) <=CO2

    # formula(13)
    for r in range(number_R):
        for j in range(number_J):
            prob += var_x[r][j] <= va[r][j] * ber_x[r][j]

    # formula(14)
    for r in range(number_R):
        for r_prime in range(number_R):
            for k in range(K):
                prob += var_h[r][r_prime][k] + var_h[r_prime][r][k] <= v*b*(1-epsilon[r][r_prime][k][l]*(1-pi[r][r_prime][k]))*ber_h[tf[r]][tf[r_prime]][k]

    solver = pl.GUROBI()
    prob.solve(solver)
    print("Status:", LpStatus[prob.status])
    print("Objective = ", value(prob.objective))

    SAVE_OBJ[l] = value(prob.objective)


print('worse')
print(max(SAVE_OBJ))

print('mean')
print(np.mean(SAVE_OBJ))


