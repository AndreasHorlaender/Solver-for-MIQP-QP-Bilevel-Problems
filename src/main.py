import os
import sys
import time
from data_reader import Data_Reader
from model_builder import Node_Builder
from node_solver import Node_Solver
from result_analyzer import Result_Analyzer
import global_vars
import parameter_settings



'''
instance_name = "bmilplib_110_4"

global_vars.matrix_spectrum = 1001

global_vars.use_opt_cuts = True
global_vars.use_HPR = False




global_vars.use_alphaBB = False

'''



instance_name = sys.argv[1] 

global_vars.matrix_spectrum = 1000


if len(sys.argv) >= 3: # else: take settings from parameter_settings.py
    if sys.argv[2] == "KKT":
        global_vars.use_HPR = False
    elif sys.argv[2] == "HPR":
        global_vars.use_HPR = True  
    else:
        raise Exception("Invalid second argument. Should be 'KKT' or 'HPR'.")
else:
    global_vars.use_HPR = parameter_settings.use_HPR

if len(sys.argv) >= 4:
    if sys.argv[3] == "use_opt_cut":
        global_vars.use_opt_cuts = True
    else:
        global_vars.use_opt_cuts = False
else:
    global_vars.use_opt_cuts = parameter_settings.use_opt_cuts      






result_file = "Results"

#if global_vars.use_opt_cuts == "True":
#    result_file += "_with_opt_cuts"
    
#if global_vars.use_alphaBB == "True":
#    result_file += "_with_alphaBB"



time_limit = 7200
global_vars.finish_time = time.time() + time_limit

# preparing data parsing
data_directory_name = "../data/"
data_directory = os.fsencode(data_directory_name)

# sort files by size (smallest first)
list_of_files = filter( lambda x: os.path.isfile(os.path.join(data_directory, x)), os.listdir(data_directory) )
files_sorted_by_size = sorted( list_of_files, key =  lambda x: os.stat(os.path.join(data_directory, x)).st_size)
for file in files_sorted_by_size:
    data_filename = os.fsdecode(file)
    if data_filename.endswith(".aux"):
        if os.path.splitext(data_filename)[0] == instance_name:
            current_file_name = os.path.splitext(data_filename)[0]
            current_mps_file = current_file_name + ".mps" 
            current_aux_file = current_file_name + ".aux"
            current_txt_file = current_file_name + "_" + str(global_vars.matrix_spectrum) + ".txt"
                
            # for sepcific aux-files the name of the other files has to be adjusted
            if data_filename.endswith("_0_100.aux") or data_filename.endswith("_50_50.aux") or data_filename.endswith("_100_0.aux"):
                current_file_name = os.path.splitext(data_filename)[0]
                current_mps_file = current_file_name[:len(current_file_name)-6] + ".mps" 
                current_aux_file = current_file_name + ".aux"
                current_txt_file = current_file_name[:len(current_file_name)-6] + "_" + str(global_vars.matrix_spectrum) + ".txt"
                
            #print(current_aux_file)
            # read csv file
            data_reader = Data_Reader()
            packed_data = data_reader.read_from_files(current_mps_file, current_aux_file, current_txt_file)
            
            
            
            global_vars.number_of_INGC = 0
            
            global_vars.KKT_model_time = 0
            global_vars.KKT_solve_time = 0
            
            global_vars.ll_model_time = 0
            global_vars.ll_modify_time = 0
            global_vars.ll_solve_time = 0
            
            global_vars.ref_model_time = 0
            global_vars.ref_modify_time = 0
            global_vars.ref_solve_time = 0
            
            
            
            
            
            # create model
            #print("Building the node ...")
            node_builder = Node_Builder(packed_data)
            node_model, x = node_builder.build()
            
            #print("Solving the node ...")
            start_time = time.time()
            node_solver = Node_Solver(node_model, x, packed_data)
            solution, rel_mip_gap, best_lower_bound, best_upper_bound, number_of_nodes_KKT, number_of_nodes_LL, number_of_INGC = node_solver.solve()
            solving_time = time.time() - start_time
            
            # save results
            result_analyzer = Result_Analyzer(packed_data, result_file, current_file_name, solution, solving_time, best_lower_bound, best_upper_bound, number_of_nodes_KKT, number_of_nodes_LL, number_of_INGC, rel_mip_gap)
            result_analyzer.save_results()
            
