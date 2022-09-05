'''
NB:
This code can get the fortification strategies
Also, all values of variables can be gained which can help you check information such as generation, transmission, cost
'''

ITER_N = 3
LIST = []
for i in range(ITER_N):
    line1 = list(input().split())
    line = [int(i) for i in line1]
    LIST.append(line)

NEED = []


for iter_n in range(ITER_N):

    UB = 48
    LB = 0

    PP=0.05
    buget_relaxation = 0.3
    print('buget_relaxation = ',buget_relaxation)
    BIG_GAMMA = 367805 + buget_relaxation*(439540-367805)
    # BIG_GAMMA = 379720 + buget_relaxation*(439540-379720)
    print(BIG_GAMMA)




    import time
    time_start = time.time()

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
    K = int(df_Param.loc['K', 'A']) # max number of tran lines
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
    va = OperationHour  # v*a
    number_J = len(alpha[0,:])
    number_R = int(df_Param.loc['number_R', 'A'])
    # beta: transmission cost
    beta = df_Param.loc['c1':'c5', 'A':'E'].values
    v = df_Param.loc['v', 'A']  # operation hour
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

    pi_final = []
    hat_theta_list = []
    threshold = 0.05
    ITERATION = 500000000
    result = []

    U = [[[0.96 for k in range(9)] for r in range(5)] for r_prime in range(5)]
    # UB = df_Param.loc['UB', 'A']

    # UB = 13
    # LB = 8
    # # LB = df_Param.loc['LB', 'A']
    # # buget_relaxation = df_Param.loc['buget_relaxation', 'A']
    #
    # buget_relaxation = 0.3
    # print('buget_relaxation = ',buget_relaxation)
    # BIG_GAMMA = 367805 + buget_relaxation*(439540-367805)
    # # BIG_GAMMA = 379720 + buget_relaxation*(439540-379720)
    # print(BIG_GAMMA)
    # # PP = df_Param.loc['PP', 'A']
    # PP=0.15


    # Ui = gamma - beta
    # for r in range(number_R):
    #     for r_prime in range(number_R):
    #         for k in range(K):
    #             if r == r_prime or ber_h[tf[r]][tf[r_prime]][k] == 0:
    #                 continue
    #             U[r][r_prime][k] = Ui[r][r_prime]


    for kk in range(ITERATION):
        if UB-LB < threshold:
            break

        epsilon_list = [[[[] for k in range(K)] for r_prime in range(number_R)] for r in range(number_R)]
        pi = [[[0 for k in range(K)] for r_prime in range(number_R)] for r in range(number_R)]  # set
        print('LB,UB:',LB,UB)

        for ll in range(ITERATION):
            hat_theta = (UB+LB)/2
            # M6; input: updated pi, hat_theta

            # <editor-fold desc="### M6 modeling ">
            prob6 = LpProblem("problem6", LpMaximize)
            var_epsilon = LpVariable.dicts("var_epsilon", (range(number_R), range(number_R), range(K)), lowBound=0,upBound=1, cat=LpInteger)
            var_lambda = LpVariable.dicts("var_lambda", (range(number_R)), lowBound=0, upBound=None, cat=LpContinuous)
            var_mu = LpVariable.dicts("var_mu", (range(number_R), range(number_J)), lowBound=0, upBound=None,cat=LpContinuous)
            var_omega = LpVariable("var_omega", lowBound=0, upBound=None, cat=LpContinuous)
            var_delta = LpVariable.dicts("var_delta", (range(number_R), range(number_R), range(K)), lowBound=0,upBound=None, cat=LpContinuous)
            var_eta = LpVariable.dicts("var_eta", (range(number_R), range(number_R), range(K)), lowBound=0, upBound=None,cat=LpContinuous)

            Objective, Obj_FirstItem, Obj_SecondItem, Obj_ThirdItem, Obj_fourthItem, Obj_fithItem = 0, 0, 0, 0, 0, 0
            for r in range(number_R):
                Obj_FirstItem += y[r] * var_lambda[r]
            for r in range(number_R):
                for j in range(number_J):
                    Obj_SecondItem -= va[r, j] * ber_x[r, j] * var_mu[r][j]
            for r in range(number_R):
                for r_prime in range(number_R):
                    for k in range(K):
                        Obj_ThirdItem -= v * b * ber_h[tf[r]][tf[r_prime]][k] * var_delta[r][r_prime][k]
            for r in range(number_R):
                for r_prime in range(number_R):
                    for k in range(K):
                        Obj_fourthItem += v * b * ber_h[tf[r]][tf[r_prime]][k] * var_eta[r][r_prime][k]
            Obj_fithItem -= CO2 * var_omega
            Objective = Obj_FirstItem + Obj_SecondItem + Obj_ThirdItem + Obj_fourthItem + Obj_fithItem
            prob6 += Objective

            # formula(6.2)
            for r in range(number_R):
                for j in range(number_J):
                    prob6 += var_lambda[r] - fff[j] * var_omega - var_mu[r][j] <= alpha[r][j]
            # formula(6.3)
            for r in range(number_R):
                for r_prime in range(number_R):
                    for k in range(K):
                        prob6 += var_lambda[r_prime] - var_lambda[r] - var_delta[r][r_prime][k] <= beta[r][r_prime]
            # formula(6.4)
            for r in range(number_R):
                prob6 += var_lambda[r] <= gamma[r]

            # # # formula(6.6)
            prob6 += lpSum([(1 - pi[r][r_prime][k]) * var_epsilon[r][r_prime][k] * ber_h[tf[r]][tf[r_prime]][k] +pi[r][r_prime][k] * var_epsilon[r][r_prime][k] * (hat_theta + sigma)] for k in range(K) for r_prime in range(number_R) for r in range(number_R)) <= hat_theta

            #  formula(6.7) —— (6.11)
            for r in range(number_R):
                for r_prime in range(number_R):
                    for k in range(K):
                        # formula(6.7)
                        prob6 += var_eta[r][r_prime][k] >= var_delta[r][r_prime][k] + var_epsilon[r][r_prime][k] * U[r][r_prime][k] - U[r][r_prime][k]
                        # formula(6.8)
                        prob6 += var_eta[r][r_prime][k] <= var_epsilon[r][r_prime][k] * U[r][r_prime][k]
                        # formula(6.9)
                        prob6 += var_eta[r][r_prime][k] <= var_delta[r][r_prime][k]
                        # formula(6.10)
                        prob6 += var_delta[r][r_prime][k] <= U[r][r_prime][k]

            solver = pl.GUROBI()
            prob6.solve(solver)

            print('C:', prob6.objective.value())
            print('LB,UB:', LB, UB)

            if prob6.objective.value() <= BIG_GAMMA:
                LB = hat_theta
                hat_theta_list.append(hat_theta)
                print('LB,UB:', LB, UB)
                break

            for r in range(number_R):
                for r_prime in range(number_R):
                    for k in range(K):
                        epsilon_list[r][r_prime][k].append(value(var_epsilon[r][r_prime][k]))

            prob8 = LpProblem("problem8", LpMaximize)

            var_phi = LpVariable("var_phi", lowBound=0, upBound=None, cat=LpContinuous)  # 标量
            var_pi = LpVariable.dicts("var_pi", (range(number_R), range(number_R), range(K)), 0, 1, LpInteger)

            prob8 += var_phi

            # formula(8.2)
            for iter in range(len(epsilon_list[0][0][0])):
                prob8 += lpSum(
                    (1 - var_pi[i][j][k]) * epsilon_list[i][j][k][iter] * ber_h[tf[i]][tf[j]][k] + var_pi[i][j][k] *epsilon_list[i][j][k][iter] * (hat_theta + sigma) for i in range(number_R) for j in range(number_R) for k in range(K)) >= var_phi

            # formula(8.3)
            prob8 += lpSum(var_pi[r][r_prime][k] * d[r][r_prime] * ber_h[tf[r]][tf[r_prime]][k] for r in range(number_R) for r_prime in range(number_R) for k in range(K)) <= PP * lpSum(d[r][r_prime] * ber_h[tf[r]][tf[r_prime]][k] for r in range(number_R) for r_prime in range(number_R) for k in range(K))

            # formula(8.4)
            for r in range(number_R):
                for r_prime in range(number_R):
                    for k in range(K):
                        if ber_h[tf[r]][tf[r_prime]][k] == 0.0:
                            prob8 += var_pi[r][r_prime][k] <= 1


            # # insert optimal pi from sampling approach
            # for r in range(number_R):
            #     for r_prime in range(number_R):
            #         for k in range(K):
            #             if r==0 and r_prime==2 and k==0:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             elif r==0 and r_prime==3 and k==1:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             elif r == 1 and r_prime == 3 and k == 4:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             # elif r==1 and r_prime==4 and k==1:
            #             #     prob8 += var_pi[r][r_prime][k] == 1
            #             # elif r==1 and r_prime==4 and k==0:
            #             #     prob8 += var_pi[r][r_prime][k] == 1
            #             else:
            #                 prob8 += var_pi[r][r_prime][k] == 0


            # # insert optimal pi from loss ranking approach
            # for r in range(number_R):
            #     for r_prime in range(number_R):
            #         for k in range(K):
            #             if r==4 and r_prime==0 and k==1:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             elif r==4 and r_prime==0 and k==2:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             elif r==4 and r_prime==0 and k==3:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             elif r==1 and r_prime==0 and k==1:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             elif r==4 and r_prime==1 and k==0:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             # elif r==3 and r_prime==0 and k==0:
            #             #     prob8 += var_pi[r][r_prime][k] == 1
            #             # elif r==4 and r_prime==1 and k==1:
            #             #     prob8 += var_pi[r][r_prime][k] == 1
            #             else:
            #                 prob8 += var_pi[r][r_prime][k] == 0


            # # insert optimal pi from utilization ranking approach
            # for r in range(number_R):
            #     for r_prime in range(number_R):
            #         for k in range(K):
            #             if r==4 and r_prime==0 and k==1:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             elif r==4 and r_prime==0 and k==2:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             elif r==4 and r_prime==0 and k==3:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             # elif r==4 and r_prime==1 and k==0:
            #             #     prob8 += var_pi[r][r_prime][k] == 1
            #             elif r==1 and r_prime==0 and k==0:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             elif r==3 and r_prime==1 and k==0:
            #                 prob8 += var_pi[r][r_prime][k] == 1
            #             # elif r==3 and r_prime==0 and k==0:
            #             #     prob8 += var_pi[r][r_prime][k] == 1
            #             else:
            #                 prob8 += var_pi[r][r_prime][k] == 0



            # insert optimal pi from approaches
            line = LIST[iter_n]  # PROTECT_LIST
            num_couple = int(len(line) / 3)
            LIS = []
            for num in range(num_couple):
                prob8 += var_pi[line[0 + 3 * num]][line[1 + 3 * num]][line[2 + 3 * num]] == 1
                LIS.append([line[0 + 3 * num], line[1 + 3 * num], line[2 + 3 * num]])
            print(LIS)
            for r in range(number_R):
                for r_prime in range(number_R):
                    for k in range(K):
                        if [r, r_prime, k] in LIS:
                            continue
                        else:
                            prob8 += var_pi[r][r_prime][k] == 0


            solver = pl.GUROBI()
            prob8.solve(solver)

            print('PHI:',prob8.objective.value())
            print('LB,UB:', LB, UB)

            if prob8.objective.value() <= hat_theta:
                UB = hat_theta
                hat_theta_list.append(hat_theta)
                print('LB,UB:', LB, UB)
                break

            # else: update pi
            number_pi_only1 = 0
            number_pi_total = 0
            for r in range(number_R):
                for r_prime in range(number_R):
                    for k in range(K):
                        pi[r][r_prime][k] = value(var_pi[r][r_prime][k])
                        if pi[r][r_prime][k]==1:
                            number_pi_total += 1
                        if pi[r][r_prime][k]==1 and ber_h[tf[r]][tf[r_prime]][k]!=0:
                            number_pi_only1 += 1
                            print(r,r_prime,k)

            print(number_pi_total)
            print(number_pi_only1)



    print(hat_theta_list)

    time_elapsed = (time.time() - time_start)
    hour = time_elapsed // 3600
    min = (time_elapsed - 3600*hour) // 60
    sec = time_elapsed - 3600*hour - 60*min

    print('computation time:', time_elapsed,)
    print(int(hour),':',int(min),':',int(sec))


    # ????????  SAVE INTO TEXT !!!!!!!
    print()
    print('pi:')
    for r in range(number_R):
        for r_prime in range(number_R):
            for k in range(K):
                if ber_h[tf[r]][tf[r_prime]][k]>0:
                    if var_pi[r][r_prime][k].value()>0:
                        print(r,r_prime,k)

    print('resilience:')
    print(round(hat_theta_list[-1]/44,2))


    NEED.append(round(hat_theta_list[-1]/44,2))


print(NEED)

print()
for i in NEED:
    print(i)

