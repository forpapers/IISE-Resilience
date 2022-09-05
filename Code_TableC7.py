'''
NB:
This code can get the breakdown information under given disruption
Also you can check all variable information according to your need
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
df_Param = pd.read_csv('param.csv', index_col=0)  # import parameters
df_Cost = pd.read_csv('cost.csv', index_col=0)
df_Demand = pd.read_csv('demand.csv', index_col=0)
df_OperationHour = pd.read_csv('operation hour.csv', index_col=0)
df_TechCapacity = pd.read_csv('tech capacity.csv', index_col=0)
K = int(df_Param.loc['K', 'A'])
tf = {0:'Guangdong', 1:'Guangxi', 2:'Hainan', 3:'Guizhou', 4:'Yunnan'}
tran_dict = {}
import csv
with open('transmission2-24.csv','r',encoding="utf-8") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row['original'] not in tran_dict.keys(): tran_dict[row['original']] = {}   # check 'original' exit
        if row['destination'] not in tran_dict[row['original']].keys(): tran_dict[row['original']][row['destination']] = []  # check 'destination' exit
        for iter in range(K):
            if iter < int(row['Quantity']):
                tran_dict[row['original']][row['destination']].append(float(row['Capacity(GW)']))
for i in range(5):   # number_R
    if tf[i] not in tran_dict.keys(): tran_dict[tf[i]] = {}
    for j in range(5):   # number_R
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
number_J = len(alpha[0,:])
number_R = int(df_Param.loc['number_R', 'A'])
y = Demand
va = OperationHour
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
fff /= rate_fff


prob = LpProblem("problem",LpMinimize)

var_x = LpVariable.dicts("Generation output", (range(number_R),range(number_J)), 0) # Generation output of tech j in region r
var_s = LpVariable.dicts("Elec deficit", (range(number_R)), 0) # Elec deficit in region r
var_h = LpVariable.dicts("Elec transmission", (range(number_R),range(number_R),range(K)), 0) # Elec transmission, r*r*k

prob += lpSum([alpha[r,j] * var_x[r][j]] for j in range(number_J) for r in range(number_R)) \
        + lpSum([beta[r_prime,r]*var_h[r_prime][r][k]] for r_prime in range(number_R) for r in range(number_R) for k in range(K)) \
        + lpSum([gamma[r] * var_s[r]] for r in range(number_R))

################# Constraints ###################

# sampling approach disruption scenario, Utilization ranking approach, Loss ranking approach
# prob += var_h[2][0][0] <= 0
# prob += var_h[4][0][1] <= 0
# prob += var_h[4][0][3] <= 0
# prob += var_h[4][1][1] <= 0
# prob += var_h[3][1][0] <= 0

# #  proposed approach disruption scenario
# prob += var_h[4][0][1] <= 0
# prob += var_h[4][0][3] <= 0
# prob += var_h[3][1][0] <= 0

# formula(1.2)
for r in range(number_R):
    prob += y[r] - var_s[r] <= lpSum([var_x[r][j]] for j in range(number_J)) + lpSum([lpSum([var_h[r_prime][r][k]] for k in range(K))] for r_prime in range(number_R)) - lpSum([lpSum([var_h[r][r_prime][k]] for k in range(K))] for r_prime in range(number_R))

# formula(1.3)
prob += lpSum([lpSum([fff[j]*var_x[r][j]] for r in range(number_R))] for j in range(number_J)) <= CO2

# formula(1.4)
for r in range(number_R):
    for j in range(number_J):
        prob += var_x[r][j] <= va[r][j] * ber_x[r][j]

# formula(1.5)
for r in range(number_R):
    for r_prime in range(number_R):
        for k in range(K):
            prob += var_h[r][r_prime][k] + var_h[r_prime][r][k] <= v*b*ber_h[tf[r]][tf[r_prime]][k]  # 可以用常数测试一下，eg:8000


solver = pl.GUROBI()
prob.solve(solver)
print("Status:", LpStatus[prob.status])
print("Objective = ", value(prob.objective))


# for v in prob.variables():
#     print(v.name, "=", v.varValue)


# for r in range(number_R):
#     for r_prime in range(number_R):
#         for k in range(K):
#             if value(var_h[r][r_prime][k]) != 0:
#                 print('var_h[',r,'][',r_prime,'][',k,']:',value(var_h[r][r_prime][k]))
#                 # print('var_h[',r,'][',r_prime,'][',k,']:',ber_h[tf[r]][tf[r_prime]][k])
#                 # print('var_h[',r,'][',r_prime,'][',k,']:',d[r][r_prime])


# for r in range(number_R):
#     for r_prime in range(number_R):
#         for k in range(K):
#             if ber_h[tf[r]][tf[r_prime]][k] != 0:
#                 if value(var_h[r][r_prime][k]) == 0 and value(var_h[r_prime][r][k]) == 0:
#                     print('var_h[',r,'][',r_prime,'][',k,']:',ber_h[tf[r]][tf[r_prime]][k],':',d[r][r_prime])


# for r in range(number_R):
#     # for j in range(number_R):
#     print('var_x:[',r,']',sum(value(var_x[r][j]) for j in range(number_R)))
#

print()
print('generation')
for r in range(number_R):
    print(r,sum(var_x[r][j].value() for j in range(number_J)))


print()
print('trans')
for r in range(number_R):
    for r_prime in range(number_R):
        for k in range(K):
            if var_h[r][r_prime][k].value()!=0:
                print(r,r_prime,k,var_h[r][r_prime][k].value())



print()
print('deficit')
for r in range(number_R):
    print('var_s[',r,']',value(var_s[r]))



print('gen cost')
print(sum(alpha[r,j] * var_x[r][j].value() for j in range(number_J) for r in range(number_R)))
print('tran cost')
print(sum(beta[r_prime,r]*var_h[r_prime][r][k].value() for r_prime in range(number_R) for r in range(number_R) for k in range(K)))
print('demand deficit')
print(sum(gamma[r] * var_s[r].value() for r in range(number_R)))


