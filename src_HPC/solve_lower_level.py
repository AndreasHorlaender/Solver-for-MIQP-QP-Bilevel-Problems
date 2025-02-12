import numpy as np
import gurobipy as gp
from gurobipy import GRB, LinExpr
from refinement import refinement_procedure
import time
import global_vars

def solve_ll(self):
    #print("mofdify ll")
    ll_modify_start_time = time.time()
    for i in range(self.b_len):
        rhs_i = self.b[i] - np.dot(self.C[i], self.x_sol)
        self.ll.getConstrByName("ll_cnstr["+str(i)+"]").rhs=rhs_i
    
    global_vars.ll_modify_time += time.time() - ll_modify_start_time
    #self.ll.write("/home/horlaender/Schreibtisch/INDEF_LL_Python/LL.lp")  
    
    ll_solve_start_time = time.time()
    #print("solve LL")
    self.ll.reset()
    self.ll.optimize()
    global_vars.ll_solve_time += time.time() - ll_solve_start_time
        
    self.nr_of_nodes_LL += self.ll.NodeCount
    
    #print("LL_Status = ", self.ll.status)
    #print("y_ll = ", np.round(self.y_ll.x, 2))
    
    
    
    #print("AAAAAAAAAAAAAAAAAAAAAAAAAA", self.c_x @ self.x_sol + self.c_y @ self.y_ll.x)
    if self.ll.status == 9:
        add_cuts = False
        global_vars.time_limit = True
        self.ll_is_feas = False
    elif self.ll.status == 2:
        self.best_ll_obj_val = self.ll.ObjVal
    else:
        self.ll_is_feas = False
    
    #print("best LL obj Val = ", self.best_ll_obj_val)
    #print("Current LL obj Val = ", self.current_ll_obj_val)
    #print("")
    
    #return self.y_ll.x, self.best_ll_obj_val, self.current_ll_obj_val
    
    
