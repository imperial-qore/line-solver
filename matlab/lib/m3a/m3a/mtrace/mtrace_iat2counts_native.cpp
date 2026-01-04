#include "mex.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>

using namespace std;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {

    if (nrhs != 3) {
        mexErrMsgIdAndTxt("mmap_toolbox:mtrace_iat2counts:nrhs", "Wrong number of input arguments");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("mmap_toolbox:mtrace_iat2counts:nlhs", "Wrong number of output arguments.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetN(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("mmap_toolbox:mtrace_iat2coutns:prhs0", "Invalid first argument: expected column vector of real number");
    }
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1]) > 1) {
        mexErrMsgIdAndTxt("mmap_toolbox:mtrace_iat2coutns:prhs1", "Invalid second argument: expected column vector of real number");
    }
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("mmap_toolbox:mtrace_iat2coutns:prhs2", "Invalid third argument: expected scalar real number");
    }
    if (mxGetM(prhs[0]) != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("mmap_toolbox:mtrace_iat2counts:prhs0_prhs1", "Mismatched vector length");
    }

    mwSize n = mxGetM(prhs[0]);
    double * T = mxGetPr(prhs[0]);
    double * A = mxGetPr(prhs[1]);
    double t = mxGetScalar(prhs[2]);

    // compute cumsum
    vector<double> C(n);
    C[0] = T[0];
    for (mwSize i = 1; i < n; ++i) {
        C[i] = C[i-1]+T[i];
    }

    // find class labels
    vector<double> L(A, A+n);
    sort(L.begin(), L.end());
    vector<double>::iterator Lend = std::unique(L.begin(), L.end());
    L.resize(distance(L.begin(),Lend));

    // get number of classes
    vector<double>::size_type m = L.size();

    // create counting process array
    mxArray * count = mxCreateDoubleMatrix(n, m, mxREAL);
    double * N = mxGetPr(count);

    // compute counting process
    mwSize i = 0;
    mwSize cur = 1;
    mwSize reached = 0;
    for (i = 0; i < (n-1); ++i) {
        mwSize prev_cur = std::max(cur,i+1);
        while (cur < n && C[cur] - C[i] <= t) {
            ++cur;
        }
        reached = i;
        // copy from previous observation
        if (i > 0) {
            for (mwSize j = 0; j < m; ++j) {
                N[n*j+i] = N[n*j+(i-1)];
            }
        } else {
            for (mwSize j = 0; j < m; ++j) {
                N[n*j+i] = 0;
            }
        }
        // discard current observation
        if (i > 0 && C[i] - C[i-1] <= t) {
            for (mwSize j = 0; j < m; ++j) {
                if (A[i] == L[j]) {
                    N[n*j+i] -= 1;
                    break;
                }
            }
        }
        // increase counters for new observations
        for (mwSize h = prev_cur; h < cur; ++h) {
            // increase counter corresponding to the arrival class
            for (mwSize j = 0; j < m; ++j) {
                if (A[h] == L[j]) {
                    N[n*j+i] += 1;
                    break;
                }
            }
        }
        // if C[cur-1] considered, exit
        if (cur == n) {
            break;
        }
    }

    // copy first k rows of count into the output argument
    mwSize k = reached+1;
    plhs[0] = mxCreateDoubleMatrix(k, m, mxREAL);
    double * Nout = mxGetPr(plhs[0]);
    for (mwSize col = 0; col < m; ++col) {
        memcpy(Nout + k*col, N + n*col, k*sizeof(double));
    }
    mxDestroyArray(count);
}
