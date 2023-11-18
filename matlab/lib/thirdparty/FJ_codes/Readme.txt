This set of scripts implements the approximation method proposed in the paper
'Beyond the Mean in Fork-Join Queues: Efficient Approximation for
Response-Time Tails', by Z. Qiu, J.F. Pérez, and P. Harrison, accepted in IFIP Performance 2015. 

An example of how to use the scripts is provided in the file main_example.m, 
where the model parameters are set and the main method (mainFJ.m) is executed. 
This method returns the approximated values of the response- time percentiles 
for a K-node FJ queue according to Section 6 of the paper. 
The approximation is based on the computation of the response-time 
percentiles for a 1-node and a 2-node FJ queue. 
The response-time percentiles for the 1-node queue are exact, while the results for a
2-node FJ queue are based on the approximation proposed in Section 4 of the paper. 

If you use the scripts available here in your work, please cite our paper. 



Copyright 2015 Imperial College London