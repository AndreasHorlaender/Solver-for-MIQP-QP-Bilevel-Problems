To run the code, go to /src and type "python main.py {instance_name} {Parameter_1} {Parameter_2}" in the terminal where
instance_name can be one of the names in /src/Filenames.txt.
Parameter_1 can be "HPR", "KKT", or "" and controls if the KKT relaxation or the HPR is taken.
Parameter_2 can be "use_opt_cut" or "" and controls if an optimality cut is used.

If Parameter_1 or Parameter_2 are not specified, i.e., the are "", then the settings in /src/parameter_settings.py are taken.
The default setting uses the KKT relaxation and no optimality cut.

For example, running "python main.py bmilplib_10_1 KKT use_opt_cut" solves the instance bmilplib_10_1 using
the KKT relaxation and an optimality cut. The results are written on /results/Results.csv (the file will be created after 
solving the first instance). The name of the .csv file is defined in /src/parameter_settings.py.
