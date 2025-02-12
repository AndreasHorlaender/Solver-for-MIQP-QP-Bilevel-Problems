import numpy as np
import gurobipy as gp
from gurobipy import GRB, LinExpr
import time
import global_vars

class Node_Builder:
    
    def __init__(self, packed_data):
        self.m = gp.Model()
        
        self.m.Params.LogToConsole = 0
        self.m.Params.Threads = 1
        
        if not global_vars.use_alphaBB:
            self.m.Params.NonConvex = 2
    
        #self.m.parameters.timelimit = 7200
        
        
        self.packed_data = packed_data
        # lenghts
        self.x_len = self.packed_data[0]
        self.y_len = self.packed_data[1]
        self.a_len = self.packed_data[6]
        self.b_len = self.packed_data[7]
        
        # bounds
        self.x_bounds = np.array(self.packed_data[2])
        self.y_bounds = np.array(self.packed_data[3])
        
        self.number_of_x_binaries = self.packed_data[4]
        self.number_of_y_binaries = self.packed_data[5]
        
          
        #print(self.number_of_x_binaries)
        # objective info
        self.c_x = np.array(self.packed_data[8])
        self.c_y = np.array(self.packed_data[9])
        self.d_y = np.array(self.packed_data[10])
        self.Q_obj = self.packed_data[11]
        self.Q_cnstr = self.packed_data[12]
        
        # LHS constraints
        self.A = self.packed_data[13]
        self.B = self.packed_data[14]
        self.C = self.packed_data[15]
        self.D = self.packed_data[16]
        
        self.D_bounds = np.vstack((self.D, np.eye(self.y_len)))
        #print(self.D_bounds)
        
        # RHS constraints
        self.a = np.array(self.packed_data[17])
        self.b = np.array(self.packed_data[18])
        self.b_bounds = np.concatenate((self.b, -np.array(self.y_bounds)[:,1]))
        #print(self.b_bounds)
        self.a_sense = np.array(self.packed_data[19])
        self.b_sense = np.array(self.packed_data[20])
        
        self.min_eigvals = self.packed_data[21]
        
        ###################################### alpha-underestimator
        self.alpha = (- self.min_eigvals/2)
        

    def build(self):
    
        def add_node_vars(self):  
            self.x = self.m.addMVar(self.x_len, vtype=GRB.INTEGER, lb=np.array(self.x_bounds)[:,0], ub=np.array(self.x_bounds)[:,1], name = "x")
            self.y = self.m.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=np.array(self.y_bounds)[:,0], ub=np.array(self.y_bounds)[:,1], name = "y")
            self.s = []
            #self.w = []
            
            if not global_vars.use_HPR:
                self.lamb = self.m.addMVar(self.b_len, vtype=GRB.CONTINUOUS, lb=np.zeros(self.b_len), name = "lambda")
                self.mu_lb = self.m.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=np.zeros(self.y_len), name = "mu_lb")
                self.mu_ub = self.m.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=np.zeros(self.y_len), name = "mu_ub")
                #self.lamb_bounds = self.m.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=np.zeros(self.y_len), name = "lambda_bounds")
                self.v = self.m.addMVar(self.b_len, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY)
                self.v_lb = self.m.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY)
                self.v_ub = self.m.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY)
                
                
            
            for j in range(self.x_len):
                self.s.append(self.m.addMVar(self.number_of_x_binaries[j], vtype=GRB.BINARY, name = "s_" + str(j)))
                
            
           
        def add_node_constrs(self):
            # lin UL cnstrs
            self.m.addConstrs((LinExpr([self.A[i][j] for j in range(self.x_len)], [self.x.tolist()[j] for j in range(self.x_len)]) + 
                                    LinExpr([self.B[i][j] for j in range(self.y_len)], [self.y.tolist()[j] for j in range(self.y_len)]) - 
                                    self.a[i] >= 0 for i in range(self.a_len)))
            
            # lin LL cnstrs
            self.m.addConstrs((LinExpr([self.C[i][j] for j in range(self.x_len)], [self.x.tolist()[j] for j in range(self.x_len)]) + 
                               LinExpr([self.D[i][j] for j in range(self.y_len)], [self.y.tolist()[j] for j in range(self.y_len)]) - 
                               self.b[i] >= 0 for i in range(self.b_len)))
            
            # binary extension
            self.m.addConstrs((self.x[j] == - 2**(self.number_of_x_binaries[j] -1) * self.s[j][self.number_of_x_binaries[j]-1]
                                               + LinExpr([2**k for k in range(self.number_of_x_binaries[j] - 1)], 
                                                         [self.s[j].tolist()[r] for r in range(self.number_of_x_binaries[j] - 1)]) for j in range(self.x_len)))
            
            
            if not global_vars.use_HPR:
                # dual feasibility
                self.m.addConstrs((LinExpr([self.Q_obj[i][j] for j in range(self.y_len)], [self.y.tolist()[j] for j in range(self.y_len)]) 
                                        + self.d_y[i] - LinExpr([self.D[j][i] for j in range(self.b_len)], 
                                                                [self.lamb.tolist()[j] for j in range(self.b_len)])
                                        - self.mu_lb[i] + self.mu_ub[i] == 0 for i in range(self.y_len)))
                
               
                
                # KKT complementarity
                self.m.addConstrs((self.v[i] == LinExpr([self.C[i][j] for j in range(self.x_len)], [self.x.tolist()[j] for j in range(self.x_len)]) + 
                                        LinExpr([self.D[i][j] for j in range(self.y_len)], [self.y.tolist()[j] for j in range(self.y_len)]) - 
                                        self.b[i] for i in range(self.b_len)))
                
                for i in range(self.y_len):
                    if not np.array(self.y_bounds)[i,0] == 0:
                        self.m.addConstr(self.v_lb[i] == np.array(self.y_bounds)[i,0] - self.y[i])
                    
                    if np.array(self.y_bounds)[i,1] >= 1e10:
                        self.m.addConstr(self.v_ub[i] == 0)
                    else:
                        self.m.addConstr(self.v_ub[i] == self.y[i] - np.array(self.y_bounds)[i,1])
                        
                        
                        
                        
                for i in range(self.b_len):
                    self.m.addSOS(GRB.SOS_TYPE1, [self.lamb[i], self.v[i]])
                    
                for i in range(self.y_len):
                    if np.array(self.y_bounds)[i,0] == 0:
                        self.m.addSOS(GRB.SOS_TYPE1, [self.mu_lb[i], self.y[i]]) 
                    else:    
                        self.m.addSOS(GRB.SOS_TYPE1, [self.mu_lb[i], self.v_lb[i]])       
                    
                for i in range(self.y_len):
                    self.m.addSOS(GRB.SOS_TYPE1, [self.mu_ub[i], self.v_ub[i]])     
                
                
            ############### Optimality Cut
            self.m.addConstr(self.c_x @ self.x + self.c_y @ self.y <= np.inf, name="opt_cut")
                
                
                
            
            
            
        def add_node_obj(self):
            self.m.setObjective(self.c_x @ self.x + self.c_y @ self.y, gp.GRB.MINIMIZE)
            
           
        KKT_model_start_time = time.time()    
        add_node_vars(self)
        add_node_constrs(self)
        add_node_obj(self)
        global_vars.KKT_model_time = time.time() - KKT_model_start_time
        #print("")
        #self.m.write("/home/horlaender/Schreibtisch/INDEF_LL_Python/model.lp")  
        return(self.m, self.x)
    
   


