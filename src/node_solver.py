import numpy as np
import gurobipy as gp
from gurobipy import GRB, LinExpr
import time

from Subproblem_modeler import build_ll
from Subproblem_modeler import build_refinement_model
from solve_lower_level import solve_ll
from add_integer_no_good_cut import add_INGC
from refinement import refinement_procedure
import global_vars

class Node_Solver:
    
    def __init__(self, node_model, x, packed_data):
        self.m = node_model
        self.x = x
        self.packed_data = packed_data
        
        self.x_len = self.packed_data[0]
        self.y_len = self.packed_data[1]
        self.y_bounds = np.array(self.packed_data[3])
        self.number_of_x_binaries = np.array(self.packed_data[4])
        self.a_len = self.packed_data[6]
        self.b_len = self.packed_data[7]
        
        self.c_x = np.array(self.packed_data[8])
        self.c_y = np.array(self.packed_data[9])
       
        self.d_y = np.array(self.packed_data[10])
        self.Q_obj = self.packed_data[11]
        
        self.A = self.packed_data[13]
        self.B = self.packed_data[14] 
        self.C = self.packed_data[15]
        self.D = self.packed_data[16] 
        self.a = np.array(self.packed_data[17])
        self.b = np.array(self.packed_data[18])
        
        
        
        self.nr_of_nodes_KKT = 0
        self.nr_of_nodes_LL = 0
        
        global_vars.time_limit = False
        self.KKT_infeasible = False
        self.optimal_UL_obj = None
        
        
        self.incumbent = 1e+18
        self.gap = 1e+18
       
    def solve(self):
        
        ############### BUILD SUBPROBLEMS ################
        build_ll(self)
        build_refinement_model(self)
        
        self.x_bf = np.zeros(self.x_len)
        self.y_bf = np.zeros(self.y_len)
        
        self.x_global_opt = None
        self.y_global_opt = None
        
        self.best_lower_bound = - 1e+18
        
        add_cuts = True
        
        while add_cuts:
            self.new_incumbent = False
            
            KKT_solve_start_time = time.time()
            self.m.Params.TimeLimit = max(0,global_vars.finish_time - time.time())
            self.m.optimize()
            global_vars.KKT_solve_time += time.time() - KKT_solve_start_time
            
            self.nr_of_nodes_KKT += self.m.NodeCount ### ADDITIV DARSTELLEN ?
    
            
            #print(self.m.status)
            
            if self.m.status == 9: ##### gurobi runs into time limit
                add_cuts = False    
                global_vars.time_limit = True
            
            elif self.m.status == 2:
                self.x_sol = np.array([var.x for var in self.m.getVars() if "x" in var.VarName])
                self.y_sol = np.array([var.x for var in self.m.getVars() if "y" in var.VarName])
                self.s_sol = np.array([var.x for var in self.m.getVars() if "s" in var.VarName])
                #self.w_sol = np.array([var.x for var in self.m.getVars() if "w" in var.VarName])
                #self.lamb_sol = np.array([var.x for var in self.m.getVars() if "lambda" in var.VarName])
                
                self.x_names = [var for var in self.m.getVars() if "x" in var.VarName]
                self.y_names = [var for var in self.m.getVars() if "y" in var.VarName]
                self.s_names = [var for var in self.m.getVars() if "s" in var.VarName]
                
                #print("x_KKT = ", self.x_sol)
                #print("y_KKT = ", np.round(self.y_sol, 2))
                
                
                
                self.best_lower_bound = self.c_x @ self.x_sol + self.c_y @ self.y_sol #max(self.best_lower_bound, self.c_x @ self.x_sol + self.c_y @ self.y_sol)
                
                ################################ CHECK LL
                self.ll.Params.TimeLimit = max(0,global_vars.finish_time - time.time())
                self.ll_is_feas = True
                solve_ll(self)
                if self.ll_is_feas:
                    self.current_ll_obj_val = 0.5*(self.y_sol @ self.Q_obj @ self.y_sol) + self.d_y @ self.y_sol
                    if self.current_ll_obj_val - self.best_ll_obj_val <= 1e-4: ### bilevel-feasible
                        if global_vars.use_opt_cuts: ### THIS IS BILEVEL OPTIMAL
                            #print("Solution is bilevel-feasible!")
                            self.x_global_opt = self.x_sol
                            self.y_global_opt = self.y_sol
                            self.incumbent = min(self.incumbent, self.c_x @ self.x_sol + self.c_y @ self.y_sol)
                            add_cuts = False
                        elif self.best_lower_bound <= self.incumbent: ### THIS IS BILEVEL OPTIMAL
                            #print("Solution is bilevel-feasible!")
                            self.x_global_opt = self.x_sol
                            self.y_global_opt = self.y_sol
                            self.incumbent = min(self.incumbent, self.c_x @ self.x_sol + self.c_y @ self.y_sol)
                            add_cuts = False    
                        else:
                            add_cuts = False ### terminate with last bilevel-feasible point
                    else:
                        #print("Solution is NOT bilevel-feasible!")
                        self.ref.Params.TimeLimit = max(0,global_vars.finish_time - time.time())
                        refinement_procedure(self)
                        
                        
                        if self.ref_is_feasible:
                            #print("REF IS FEAS")
                            if self.c_x @ self.x_bf + self.c_y @ self.y_bf < self.incumbent:
                                self.new_incumbent = True
                                self.incumbent = min(self.incumbent, self.c_x @ self.x_bf + self.c_y @ self.y_bf)
                                self.x_global_opt = self.x_bf
                                self.y_global_opt = self.y_bf
                                #print("Best UL ZF = ", self.incumbent)
                                
                            if global_vars.use_opt_cuts and self.new_incumbent:
                                self.m.getConstrByName("opt_cut").rhs=self.incumbent
                                self.m.update()
                          
                        #else:
                            #print("REF IS INFEAS")
                        add_INGC(self)
                        #print("INC = ", self.incumbent)
                        
                       
                else:
                    add_INGC(self)
                
            else:
                #print("STATUS = ", self.m.status)
                add_cuts = False
                #print("KKT relaxation is now infeasible.")
                if global_vars.number_of_INGC == 0: ### then this was the first KKT solve
                    self.KKT_infeasible = True
                    #print("Bilevel problem is infeasible!")
                    #self.m.computeIIS()
                    #self.m.write("/home/horlaender/Schreibtisch/INDEF_LL_Python/model.ilp")
                    
                    
                    self.KKT_without_UL_feas = False
                    self.m.remove(self.m.getConstrs()[0:self.a_len])
                    self.m.update()
                    self.m.optimize()
                    if self.m.status == 2:
                        self.KKT_without_UL_feas = True
                    #else:
                    #    print("aaaaaaaaaaaa", self.m.status)
                    
                    
            ################## TERMINATE DUE TO GAP        
            if self.gap <= 1e-4 or self.best_lower_bound >= self.incumbent - 1e-6:
                add_cuts = False
            #else:
            #    add_cuts = False
                
                
                
        ######## HANDLE SOLUTION
        if not abs(self.incumbent) <= 1e-4:
            self.gap = abs(self.incumbent - self.best_lower_bound)/abs(self.incumbent)
            
            
        
            
        
        if global_vars.time_limit:
            self.optimal_UL_obj = "time_limit"
        elif self.KKT_infeasible:
            if self.KKT_without_UL_feas:
                self.optimal_UL_obj = "KKT_infeasible_but_feas_without_UL"
            else:    
                self.optimal_UL_obj = "KKT_infeasible"
        else:
            self.optimal_UL_obj = self.incumbent
            
        #print("BEST UL_obj = ", self.optimal_UL_obj)    
        #print("best_upper_bound = ", self.incumbent)
        #print("best lower bound = ", self.best_lower_bound)
        
        return self.optimal_UL_obj, self.gap, self.best_lower_bound, self.incumbent, self.nr_of_nodes_KKT, self.nr_of_nodes_LL, global_vars.number_of_INGC
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            #if self.m.status == 2:
            #    print(solution)
            #    return 1, rel_mip_gap, best_upper_bound, number_of_nodes
            #else:
            #    return self.m.status, rel_mip_gap, best_upper_bound, number_of_nodes
