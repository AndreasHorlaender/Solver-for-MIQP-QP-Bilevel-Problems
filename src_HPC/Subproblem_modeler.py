import numpy as np
import gurobipy as gp
from gurobipy import GRB, LinExpr
import time
import global_vars

def build_ll(self):
    ll_model_start_time = time.time()
    self.ll = gp.Model()
    self.ll.Params.LogToConsole = 0
    self.ll.Params.Threads = 1
    self.ll.Params.NonConvex = 2
    
    self.y_ll = self.ll.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=np.array(self.y_bounds)[:,0], ub=np.array(self.y_bounds)[:,1], name = "y")
    
    self.ll.addConstrs((LinExpr([self.D[i][j] for j in range(self.y_len)], [self.y_ll.tolist()[j] for j in range(self.y_len)]) >= 0 for i in range(self.b_len)), name="ll_cnstr")
    
    self.ll.setObjective(0.5*(self.y_ll @ self.Q_obj @ self.y_ll) + self.d_y @ self.y_ll, gp.GRB.MINIMIZE)
    
    self.ll.update()
    global_vars.ll_model_time = time.time() - ll_model_start_time
   
def build_refinement_model(self):
    ref_model_start_time = time.time()
    self.ref = gp.Model()
    
    self.ref.Params.LogToConsole = 0
    self.ref.Params.Threads = 1
    self.ref.Params.NonConvex = 2
    
    self.y_ref = self.ref.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=np.array(self.y_bounds)[:,0], ub=np.array(self.y_bounds)[:,1], name = "y_ref")
    
    # lin UL cnstrs
    self.ref.addConstrs((self.B[i] @ self.y_ref >= 0 for i in range(self.a_len)), name="ref_ul_cnstr")
    # lin LL cnstrs
    self.ref.addConstrs((self.D[i] @ self.y_ref >= 0 for i in range(self.b_len)), name="ref_ll_cnstr")
   
    # we min c_y @ y because x is fixed and, hence, the term c_x @ x would be constant
    self.ref.setObjective(self.c_y @ self.y_ref, gp.GRB.MINIMIZE) 
    
    self.ref.update()
    global_vars.ref_model_time = time.time() - ref_model_start_time
