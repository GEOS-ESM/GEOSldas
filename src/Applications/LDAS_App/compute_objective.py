def compute_objective(exp_dir,num_ens):
    for ens in range(num_ens):
        output_dir = exp_dir+'/output/SMAP_EASEv2_M36/cat/ens'+ens
        for (dirpath,dirnames,filenames) in os.walk(output_dir):
            for this_dir in dirnames:
                curr_dir = output_dir+this_dir
