import numpy as np
dir_path='.'

d=2
for N in [5,10,20,100,1000]:
    for rho in [.90,.95]:
        initial_number_of_jobs = 2.8 if rho==.90 else 3.6
        oldFileName = '{}/traj/averageTraj_N{}_r{}_init{}.npz'.format(
            dir_path,N,int(rho*100),int(initial_number_of_jobs*10))
        newFileName = '{}/traj/averageTraj_N{}_d{}_r{}_init{}.npz'.format(
            dir_path,N,d,int(rho*100),int(initial_number_of_jobs*10))
        a = np.load(newFileName)
        np.savez_compressed(newFileName,x=a['x'],nb_samples=a['nb_samples'])
