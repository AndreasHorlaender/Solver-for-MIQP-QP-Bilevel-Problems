import numpy as np
import time
import global_vars

def refinement_procedure(self):
    def modify_ref_constrs(self):  
        ref_modify_start_time = time.time()
        # modify HPR constraints
        for i in range(self.a_len):
            rhs_i = self.a[i] - self.A[i] @ self.x_sol
            self.ref.getConstrByName("ref_ul_cnstr["+str(i)+"]").rhs=rhs_i
                                  
        for i in range(self.b_len):
            rhs_i = self.b[i] - self.C[i] @ self.x_sol
            self.ref.getConstrByName("ref_ll_cnstr["+str(i)+"]").rhs=rhs_i
                
        # add refinement constraint
        self.ref_cnstr = self.ref.addConstr(self.d_y @ self.y_ref +  0.5 * (self.y_ref @ self.Q_obj @ self.y_ref) <= self.best_ll_obj_val)
        
        global_vars.ref_modify_time += time.time() - ref_modify_start_time
        #self.ref.write("/home/horlaender/Schreibtisch/INDEF_LL_Python/REF.lp")  
        
    def solve_ref_model(self):
        self.ref_is_feasible = True
        ref_solve_start_time = time.time()
        self.ref.reset()
        self.ref.optimize()
        global_vars.ref_solve_time += time.time() - ref_solve_start_time
        
            
        if self.ref.status == 9:
            add_cuts = False
            global_vars.time_limit = True
            self.ref_is_feasible = False
        elif self.ref.status == 2:
            self.x_bf = self.x_sol
            for i in range(self.y_len):
                self.y_bf[i] = self.y_ref[i].x
        else:
            #print("The restricted HPR is infeasible.")
            self.ref_is_feasible = False
            #self.y_bf = None
            
        
        
    #print("Apply refinement.")    
    modify_ref_constrs(self)
    solve_ref_model(self)
    
    #print(self.d_y @ self.y_ref.x +  0.5 * (self.y_ref.x @ self.Q_obj @ self.y_ref.x))
    #print(self.y_bf)
    #err
    
    self.ref.remove(self.ref_cnstr)
    
    self.ref.update()
    
    #return self.y_bf
    
    
