import numpy as np
import gurobipy as gp
from gurobipy import GRB, LinExpr

class Node_Builder:
    
    def __init__(self, packed_data):
        self.m = gp.Model()
        
        self.m.Params.LogToConsole = 0
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
        
        self.D_bounds = np.vstack((self.D, -np.eye(self.y_len)))
        #print(self.D_bounds)
        
        # RHS constraints
        self.a = np.array(self.packed_data[17])
        self.b = np.array(self.packed_data[18])
        self.b_bounds = np.concatenate((self.b, -np.array(self.y_bounds)[:,1]))
        #print(self.b_bounds)
        self.a_sense = np.array(self.packed_data[19])
        self.b_sense = np.array(self.packed_data[20])
        
        
        #self.A = self.A * (-1)
        #self.B = self.B * (-1)
        #self.a = self.a * (-1)
        #self.C = self.C * (-1)
        #self.D = self.D * (-1)
        #self.b = self.b * (-1)
        
        

    def build(self):
    
        def add_node_vars(self):  
            self.x = self.m.addMVar(self.x_len, vtype=GRB.INTEGER, lb=np.array(self.x_bounds)[:,0], ub=np.array(self.x_bounds)[:,1], name = "x")
            self.y = self.m.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=np.array(self.y_bounds)[:,0], ub=np.array(self.y_bounds)[:,1], name = "y")
            self.lamb = self.m.addMVar(self.b_len + self.y_len, vtype=GRB.CONTINUOUS, lb=np.zeros(self.b_len + self.y_len), name = "lambda")
            #self.lamb_bounds = self.m.addMVar(self.y_len, vtype=GRB.CONTINUOUS, lb=np.zeros(self.y_len), name = "lambda_bounds")
            
            self.s = []
            self.w = []
            
            for j in range(self.x_len):
                self.s.append(self.m.addMVar(self.number_of_x_binaries[j], vtype=GRB.BINARY, name = "s_" + str(j)))
                
            for j in range(self.x_len):
                self.w.append(self.m.addMVar(self.number_of_x_binaries[j], vtype=GRB.BINARY, name = "w_" + str(j)))    
            
            
            
            self.lamb_ub = 1e10
            self.lamb_lb = -1e10
            
           
        def add_node_constrs(self):
            
            # lin UL cnstrs
            self.m.addConstrs((LinExpr([self.A[i][j] for j in range(self.x_len)], [self.x[j] for j in range(self.x_len)]) + 
                                    LinExpr([self.B[i][j] for j in range(self.y_len)], [self.y[j] for j in range(self.y_len)]) - 
                                    self.a[i] <= 0 for i in range(self.a_len)))
            
            # lin LL cnstrs
            self.m.addConstrs((LinExpr([self.C[i][j] for j in range(self.x_len)], [self.x[j] for j in range(self.x_len)]) + 
                               LinExpr([self.D[i][j] for j in range(self.y_len)], [self.y[j] for j in range(self.y_len)]) - 
                               self.b[i] <= 0 for i in range(self.b_len)))
            
            # feasibility of lower-level dual
            self.m.addConstrs((LinExpr([self.Q_obj[i][j] for j in range(self.y_len)], [self.y[j] for j in range(self.y_len)]) 
                                    + self.d_y[i] - LinExpr([self.D_bounds[j][i] for j in range(self.b_len + self.y_len)], 
                                                            [self.lamb[j] for j in range(self.b_len + self.y_len)]) <= 0 for i in range(self.y_len)))
            
            # binary extension
            self.m.addConstrs((self.x[j] == - 2**(self.number_of_x_binaries[j] -1) * self.s[j][self.number_of_x_binaries[j]-1]
                                               + LinExpr([2**k for k in range(self.number_of_x_binaries[j] - 1)], 
                                                         [self.s[j][r] for r in range(self.number_of_x_binaries[j] - 1)]) for j in range(self.x_len)))
            
            # add linearization constraints
            self.m.addConstrs((self.w[j][r] <= self.lamb_ub*self.s[j][r] for j in range(self.x_len) for r in range(self.number_of_x_binaries[j])))
            
            self.m.addConstrs((self.w[j][r] <= LinExpr([self.C[i][j] for i in range(self.b_len)], [self.lamb[i] for i in range(self.b_len)]) + 
                                    self.lamb_lb*(self.s[j][r]-1) for j in range(self.x_len) for r in range(self.number_of_x_binaries[j])))
            
            self.m.addConstrs((self.w[j][r] >= self.lamb_lb*self.s[j][r] for j in range(self.x_len) for r in range(self.number_of_x_binaries[j])))
            
            self.m.addConstrs((self.w[j][r] >= LinExpr([self.C[i][j] for i in range(self.b_len)], [self.lamb[i] for i in range(self.b_len)]) + 
                                    self.lamb_ub*(self.s[j][r]-1) for j in range(self.x_len) for r in range(self.number_of_x_binaries[j])))
            
            
        def add_node_obj(self):
            self.pi = 1e5
            
            #self.m.setObjective(1, GRB.MINIMIZE)
            
            self.m.setObjective(self.c_x @ self.x + self.c_y @ self.y + self.pi*(0.5*(self.y @ self.Q_obj @ self.y) + self.d_y @ self.y - self.b_bounds @ self.lamb
                                        + sum(LinExpr([2**k for k in range(self.number_of_x_binaries[j])], [self.w[j][r] for r in range(self.number_of_x_binaries[j])]) 
                                        for j in range(self.x_len) )), gp.GRB.MINIMIZE)
                                 

            
        add_node_vars(self)
        add_node_constrs(self)
        add_node_obj(self)
        print("")
        #self.m.export_as_lp("/home/horlaender/Schreibtisch/CPLEX_model/src_HPC/model.lp")  
        return(self.m, self.x)
    
   


