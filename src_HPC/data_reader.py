import numpy as np
from scipy import linalg
import cplex
import pandas as pd
import math
#import matlab.engine
        
class Data_Reader:       
        
    def store_mps_info(self):
        # save UL-objective coeff
        self.UO = self.mps_info.objective.get_linear()
       
        # save coeff of constraint matrix ABCD
        self.number_of_vars = self.mps_info.variables.get_num()
        self.number_of_cnstrs = self.mps_info.linear_constraints.get_num()
        
        # split the SparsePairs into 2 lists
        self.ABCD_ind = []
        self.ABCD_values = []
        for i in range(self.number_of_cnstrs):
            #print(self.mps_info.linear_constraints.get_rows()[i])
            self.ABCD_ind.append(self.mps_info.linear_constraints.get_rows()[i].unpack()[0])
            self.ABCD_values.append(self.mps_info.linear_constraints.get_rows()[i].unpack()[1])
            
        # save RHS of the constraints formed by ABCD
        self.RHS = self.mps_info.linear_constraints.get_rhs()
        self.sense = self.mps_info.linear_constraints.get_senses()
        #print("sense = ", self.sense)
        
        # save variable bounds for x and y
        self.bounds = []
        for i in range(len(self.mps_info.variables.get_lower_bounds())):
            self.bounds.append([self.mps_info.variables.get_lower_bounds(i), self.mps_info.variables.get_upper_bounds(i)])
            
            
    def store_aux_info(self):
        ### Store LL info
        self.UC = [] # UL Vars
        self.LC = [] # LL Vars
        self.UR = [] # UL Constraints LHS (coeff)
        self.LR = [] # LL Constraints LHS (coeff)
        self.a = []  # UL Constraints RHS
        self.b = []  # LL Constraints RHS
        self.a_sense = [] # <,>,=
        self.b_sense = [] # <,>,=
        
        self.d_y = [] # part of LL objective
        self.OS = 1 # objective scaling
    
        for i in range(self.LL_info.shape[0]):
            if self.LL_info[0][i] == "LC":
                self.LC.append(self.LL_info[1][i])
            elif self.LL_info[0][i] == "LR":
                self.LR.append(self.LL_info[1][i])
            elif self.LL_info[0][i] == "LO":
                self.d_y.append(self.LL_info[1][i]) 
            elif self.LL_info[0][i] == "OS": 
                self.OS = self.LL_info[1][i]
                
        self.d_y = np.multiply(self.d_y, self.OS)
                
        ### Store UL info        
        self.c_x = [] # part of UL objective
        self.c_y = [] # part of UL objective
        
        self.x_bounds = []
        self.y_bounds = []
        
        for i in range(self.number_of_vars): # für alle Vars
            if i not in self.LC:
                self.UC.append(i)
                self.c_x.append(self.UO[i])
                self.x_bounds.append(self.bounds[i])
            else:
                self.c_y.append(self.UO[i])
                self.y_bounds.append(self.bounds[i])
               
        for i in range(self.number_of_cnstrs): # für alle NB
            if i not in self.LR:
                self.UR.append(i)     
                self.a.append(self.RHS[i])
                self.a_sense.append(self.sense[i])
            else:
                self.b.append(self.RHS[i])
                self.b_sense.append(self.sense[i])
        
            
    def init_rest(self): 
        self.x_len = len(self.c_x)
        #print("x_len =", self.x_len)
        self.y_len = len(self.c_y)
        #print("y_len =", self.y_len)
        self.a_len = len(self.a)
        #print("a_len =", self.a_len)
        self.b_len = len(self.b)
        
        
        self.number_of_x_binaries = []
        for i in range(self.x_len):
            assert np.array(self.x_bounds)[i,1] >= 1
            self.number_of_x_binaries.append(int(math.log(np.array(self.x_bounds)[i,1], 2) + 2))
            
        self.number_of_y_binaries = []
        for i in range(self.y_len):
            assert np.array(self.y_bounds)[i,1] >= 1
            self.number_of_y_binaries.append(int(math.log(np.array(self.y_bounds)[i,1], 2) + 2))   
        
        self.A_ind = []
        self.A_values = []
        self.B_ind = []
        self.B_values = []
        self.C_ind = []
        self.C_values = []
        self.D_ind = []
        self.D_values = []
        
        self.I_init = []       
      
    def split_matrix_ABCD(self):
        self.ABCD = np.zeros((self.a_len + self.b_len, self.x_len + self.y_len))
        
        for row_ind in range(self.a_len + self.b_len):
            self.ABCD[row_ind][self.ABCD_ind[row_ind]] = self.ABCD_values[row_ind]
        #print(self.ABCD)    
        
        self.A = np.zeros((self.a_len, self.x_len))
        self.B = np.zeros((self.a_len, self.y_len))
        self.C = np.zeros((self.b_len, self.x_len))
        self.D = np.zeros((self.b_len, self.y_len))
        
        self.A += self.ABCD[list(map(int, self.UR)), :][:, list(map(int, self.UC))]
        self.B += self.ABCD[list(map(int, self.UR)), :][:, list(map(int, self.LC))]
        self.C += self.ABCD[list(map(int, self.LR)), :][:, list(map(int, self.UC))]
        self.D += self.ABCD[list(map(int, self.LR)), :][:, list(map(int, self.LC))]
        #print(self.D)
        
    def transform_problem_in_standard_form(self):
        for i in range(self.a_len):
            if self.a_sense[i] == 'L':
                self.a_sense[i] = 'G'
                self.a[i] = - self.a[i]
                for j in range(self.x_len):
                    self.A[i][j] = - self.A[i][j]
                for k in range(self.y_len):
                    self.B[i][k] = - self.B[i][k]
            elif self.a_sense[i] == 'E':
                self.a_sense[i] = 'G'
                    
                self.a_sense.append('G')
                self.a = np.hstack([self.a, -self.a[i]])
                self.A = np.vstack([self.A, -self.A[i]])
                self.B = np.vstack([self.B, -self.B[i]])
                self.a_len += 1
            else:
                assert self.a_sense[i] == 'G'
                
        #print("A = ", self.A)
        #print("B = ", self.B)
        #print("C = ", self.C)
        #print("D = ", self.D)
        #print("a = ", self.a)
        #print("b = ", self.b)
        #eror
         
        for i in range(self.b_len):
            if self.b_sense[i] == 'L':
                self.b_sense[i] = 'G'
                self.b[i] = - self.b[i]
                for j in range(self.x_len):
                    self.C[i][j] = - self.C[i][j]
                for k in range(self.y_len):
                    self.D[i][k] = - self.D[i][k]
            elif self.b_sense[i] == 'E':
                self.b_sense[i] = 'G'
                    
                self.b_sense.append('G')
                self.b = np.hstack([self.b, -self.b[i]])
                self.C = np.vstack([self.C, -self.C[i]])
                self.D = np.vstack([self.D, -self.D[i]])
                self.b_len += 1
            else:
                assert self.b_sense[i] == 'G'    
                
                
        #print("b = ", self.b)
        #print("d_y = ", self.d_y)
        
    def read_other_matrices(self, file):
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
      
        self.Q_obj, self.Q_cnstr = [], []
        
        for i in range(1, 1+self.y_len):
            self.Q_obj.append(list(map(int, lines[i].split())))
        self.Q_obj = np.array(self.Q_obj)
        
        eigenvalues = np.linalg.eigvals(self.Q_obj)
        #print("eigenvalues = ", eigenvalues)
        #print("eigenvalues_2 = ", linalg.eigvals(self.Q_obj))
        #print("Eigenvalue Q_obj MIN = ", min(eigenvalues))
        #print("Eigenvalue Q_obj MAX = ", max(eigenvalues))
        
        self.min_eigvals = min(eigenvalues)
            
        
        
    def read_from_files(self, current_mps_file, current_aux_file, current_txt_file): 
        
        with open(f'../data/{current_mps_file}'):
            with open(f'../data/{current_aux_file}'):     
                self.mps_info = cplex.Cplex("../data/" + current_mps_file)
                #self.mps_info.write("mps_info.lp")
                #self.mps_info = gp.read("../data/" + current_mps_file)   
                self.LL_info = pd.read_csv("../data/" + current_aux_file, sep=" ", header=None) 
                self.store_mps_info()
                self.store_aux_info()
                self.init_rest()
                self.split_matrix_ABCD()  
                self.transform_problem_in_standard_form()
                with open(f'../data/{current_txt_file}') as file:
                        self.read_other_matrices(file)  
                
        self.packed_data = [self.x_len,
                            self.y_len,
                            self.x_bounds,
                            self.y_bounds,
                            self.number_of_x_binaries,
                            self.number_of_y_binaries,
                            self.a_len,
                            self.b_len,
                            self.c_x,
                            self.c_y,
                            self.d_y,
                            self.Q_obj,
                            self.Q_cnstr,
                            self.A,
                            self.B,
                            self.C,
                            self.D,
                            self.a,
                            self.b,
                            self.a_sense,
                            self.b_sense,
                            self.min_eigvals]
            
        return self.packed_data           
