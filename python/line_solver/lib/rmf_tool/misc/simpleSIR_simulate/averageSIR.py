from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
import os
import time

path = os.path.abspath(__file__)
dir_path = os.path.dirname(path)
os.system('cd {} && make'.format(dir_path))

def averageTraj(N,a,nbSamples=1000):
    fileName = '{}/traj/trajN{}a{}.npz'.format(dir_path,N,a)
    if os.path.exists(fileName):
        results = np.load(fileName)
        nbSamplesComputed = results['nbSamples']
        T = results['T']
        S = results['S']
        I = results['I']
        if nbSamplesComputed > nbSamples: return(T,S/nbSamplesComputed,I/nbSamplesComputed)
    else:
        nbSamplesComputed = 0
        S = np.zeros(1000)
        I = np.zeros(1000)
    T = np.linspace(0,10,1000)
    ti = time.time()
    for i in range(nbSamplesComputed,nbSamples):
        os.system('{0}/sir_simu N{1} a{2} t > {0}/traj/tmp_N{1}_a{2}'.format(dir_path,N,a))
        result = np.array(pd.read_csv('{0}/traj/tmp_N{1}_a{2}'.format(dir_path,N,a),sep=' ',header=None))
        f_S = interp1d(result[:,0],result[:,1],kind='zero',assume_sorted=True)
        f_I = interp1d(result[:,0],result[:,2],kind='zero',assume_sorted=True)
        S += f_S(T)
        I += f_I(T)
    np.savez(fileName,nbSamples=nbSamples,T=T,S=S,I=I)
    print(N,a,'computed in ',time.time()-ti,'seconds')
    results = np.load(fileName)
    return(results['T'],results['S']/nbSamples,results['I']/nbSamples)

def steadyState(N,a):
    fileName = '{}/traj/steadyStateN{}a{}.txt'.format(dir_path,N,a)
    if not os.path.exists(fileName):
        os.system('{}/sir_simu N{1} a{2} > {3}'.format(dir_path,N,a,fileName))
    return np.mean(np.loadtxt(fileName),0)

