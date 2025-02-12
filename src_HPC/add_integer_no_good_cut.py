import numpy as np
import gurobipy as gp
import global_vars

def add_INGC(self):
    coeff_list = []
    lhs_constant = 0
    for i in range(len(self.s_names)):
        if self.s_sol[i] <= 1e-4:
            coeff_list.append(1)
        else:
            coeff_list.append(-1)
            lhs_constant += 1
    
    a = []        
    for i in range(len(self.s_sol)):
        a.append(self.s_sol[i])
    print(a)        
    #print(coeff_list)
    
    #print(self.s_names[0])
    
    self.m.addConstr(gp.quicksum([coeff_list[i] * self.s_names[i] for i in range(len(coeff_list))]) >= 1 - lhs_constant)
    self.m.update()
    
    global_vars.number_of_INGC += 1
    
            