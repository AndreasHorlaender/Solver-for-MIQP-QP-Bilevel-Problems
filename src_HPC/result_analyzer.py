import csv 
import os
import global_vars
import numpy as np
from scipy.linalg import eigh

class Result_Analyzer():
    
    def __init__(self, packed_data, result_file, current_file_name, solution, solving_time, best_lower_bound, best_upper_bound, number_of_nodes_KKT, number_of_nodes_LL, number_of_INGC, rel_mip_gap):
        self.packed_data = packed_data
        self.x_len = self.packed_data[0]
        self.y_len = self.packed_data[1]
        self.a_len = self.packed_data[6]
        self.b_len = self.packed_data[7]
        self.c_x = np.array(self.packed_data[8])
        self.c_y = np.array(self.packed_data[9])
        self.d_y = np.array(self.packed_data[10])
        
        # compute min/max objective coeffs
        self.min_ulo = min(min(self.c_x), min(self.c_y))
        self.max_ulo = max(max(self.c_x), max(self.c_y))
        self.min_llo = min(self.d_y)
        self.max_llo = max(self.d_y)
        
        # compute min/max eigenvalues of Q and P
        self.Q_obj = self.packed_data[11]
        self.Q_cnstr = self.packed_data[12]
        self.Q_eig_min = min(eigh(self.Q_obj, eigvals_only=True))
        self.Q_eig_max = max(eigh(self.Q_obj, eigvals_only=True))
        
        self.result_file = result_file
        self.current_file_name = current_file_name
        
        
        self.solution = solution
        self.solving_time = solving_time
        self.best_lower_bound = best_lower_bound
        self.best_upper_bound = best_upper_bound
        self.number_of_nodes_KKT = number_of_nodes_KKT
        self.number_of_nodes_LL = number_of_nodes_LL
        self.number_of_INGC = number_of_INGC
        self.rel_mip_gap = rel_mip_gap
        
        self.KKT_model_time = global_vars.KKT_model_time 
        self.KKT_solve_time = global_vars.KKT_solve_time 
        self.ll_model_time = global_vars.ll_model_time 
        self.ll_modify_time = global_vars.ll_modify_time 
        self.ll_solve_time = global_vars.ll_solve_time 
        self.ref_model_time = global_vars.ref_model_time 
        self.ref_modify_time = global_vars.ref_modify_time 
        self.ref_solve_time = global_vars.ref_solve_time 
        
        self.spectrum = global_vars.matrix_spectrum
        self.method = None
        
        if global_vars.use_HPR:
            if global_vars.use_opt_cuts:
                self.method = "HPR_w/_opt_cuts"
            else:
                self.method = "HPR_w/o_opt_cuts"
        else:
            if global_vars.use_opt_cuts:
                self.method = "KKT_w/_opt_cuts"
            else:
                self.method = "KKT_w/o_opt_cuts"
                
        
        
        if global_vars.time_limit:
            self.solution = 'time_limit'
            self.number_of_nodes_KKT = 9999999999999
            self.number_of_INGC = 9999999999999
            self.solving_time = 9999999999999
            if self.rel_mip_gap == 'nan':
                self.rel_mip_gap = 9999999999
        if self.rel_mip_gap == 'nan':
            self.rel_mip_gap = 0   
        
        
        self.most_time_consumed = None
        if max(self.KKT_solve_time, self.ll_solve_time, self.ref_solve_time) == self.KKT_solve_time:
            self.most_time_consumed = 'KKT_solve_time'
        elif max(self.KKT_solve_time, self.ll_solve_time, self.ref_solve_time) == self.ll_solve_time:
            self.most_time_consumed = 'll_solve_time'    
        else:
            self.most_time_consumed = 'ref_solve_time'   
        
    def save_results(self):
        header = ['name', 'objective value', 'method', 'LL_spectrum', 'use_opt_cuts', 'use_alphaBB',
                  'sum_nodes_KKT', 'sum_nodes_LL', 'time', 'most_time_consumed', 'INGCs', 'KKT_model_time', 'KKT_solve_time',
                  'll_model_time', 'll_modify_time', 'll_solve_time', 'ref_model_time', 'ref_modify_time', 'ref_solve_time',
                  'mip_gap', 'best_lower_bound', 'best_upper_bound', 'nr_UL_vars',
                  'nr_LL_vars', 'nr_UL_cnstrs', 'nr_LL_cnstrs', 'min_UL_obj_coeff', 'max_UL_obj_coeff', 
                  'min_LL_obj_coeff', 'max_LL_obj_coeff', 'Q_eig_min', 'Q_eig_max']
        
        data = [self.current_file_name, self.solution, self.method, self.spectrum, global_vars.use_opt_cuts, global_vars.use_alphaBB,
                self.number_of_nodes_KKT, self.number_of_nodes_LL, self.solving_time, self.most_time_consumed, self.number_of_INGC, 
                self.KKT_model_time, self.KKT_solve_time, self.ll_model_time, 
                self.ll_modify_time, self.ll_solve_time, self.ref_model_time, 
                self.ref_modify_time, self.ref_solve_time, self.rel_mip_gap, self.best_lower_bound,
                self.best_upper_bound, self.x_len, self.y_len, self.a_len, self.b_len, 
                self.min_ulo, self.max_ulo, self.min_llo, self.max_llo, self.Q_eig_min, self.Q_eig_max]
            
        with open('../results/' + self.result_file + '.csv','a', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=';')
            
            if os.stat('../results/' + self.result_file + '.csv').st_size == 0: #wenn die Datei leer ist
                writer.writerow(header)
            writer.writerow(data)
