#include "mex.h"
#include <string.h> // for memset
/*
 * Computing connected components of an undirected graph - assuming sA is symmetric
 *
 * Usage:
 *  [l c] = graph_conn_comp_mex(sA);
 *
 * Inputs:
 *  sA - sparse adjacency matrix (for directed graph - does not have to be symmetric)
 *
 * Outputs:
 *  l - components labels
 *  c - number of connected components
 *
 *
 * Compile using:
 * >> mex -O -largeArrayDims graph_conn_comp_mex.cpp
 *
 */
#line   __LINE__  "graph_conn_comp_mex"
#define     STR(s)      #s  
#define     ERR_CODE(a,b)   a ":" "line_" STR(b)
// inputs
enum {
    AIN = 0,
    NIN };
// outputs
enum {
    LOUT = 0,
    COUT,
    NOUT };
void
ConnComp(const mxArray* sA, unsigned int* pl, const double* pnc, mwIndex start_node);
mwIndex
FindUnLabeled(unsigned int* pl, mwIndex n);

void mexFunction( 
    int nout,
    mxArray* pout[],
    int nin,
    const mxArray* pin[])  {
    if ( nin != NIN )
         mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"must have %d inputs", NIN);
    if (nout==0)
        return;
    if (nout != NOUT )
         mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"must have exactly %d output", NOUT);

    if ( mxIsComplex(pin[AIN]) || !mxIsSparse(pin[AIN]) )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"sA must be real sparse matrix");
    if ( mxGetM(pin[AIN]) != mxGetN(pin[AIN]) )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"sA must be square matrix");

    mwIndex n = mxGetM(pin[AIN]); // number of variables
    mwIndex ii(0);
    // allocate outputs
    pout[LOUT] = mxCreateNumericMatrix(1, n, mxUINT32_CLASS, mxREAL);
    unsigned int* pl = (unsigned int*)mxGetData(pout[LOUT]);
    memset(pl, 0, n*sizeof(unsigned int)); // set initial labels to zero
    unsigned int cl = 0;
    // number of components
    pout[COUT] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* pnc = mxGetPr(pout[COUT]); 
    pnc[0] = 0; // number of components
    ii = 0;
    do {
        ConnComp(pin[AIN], pl, pnc, ii); // find conn comp
        pnc[0]++;        
        ii = FindUnLabeled(pl, n);
    } while ( ii < n );
}

mwIndex
FindUnLabeled(unsigned int* pl, mwIndex n) {
    for ( mwIndex ii(0); ii<n ; ii++ ){
        if ( pl[ii]==0 ) {
            return ii;
        }
    }
    return n+1;
}

void
ConnComp(const mxArray* sA, unsigned int* pl, const double* pnc, mwIndex start_node) {
    mwIndex n = mxGetM(sA);
    unsigned int curr_label = (unsigned int)(pnc[0]+1);
    mwIndex *stack = (mwIndex*)mxMalloc( (n)*sizeof(mwIndex) ); 
    memset(stack, 0, (n)*sizeof(mwIndex));
    mwIndex sp(0); // stack pointer
    // put start_node on top of stack after labeling it
    pl[start_node]=curr_label;
    stack[sp] = start_node;
    sp++;    
    mwIndex* pir = mxGetIr(sA);
    mwIndex* pjc = mxGetJc(sA);
    mwIndex  ii(0), neighbor(0);
    while ( sp > 0 ) {       
        // pop start_label from stack
        sp--;
        start_node = stack[sp];
        for ( ii = pjc[start_node] ; ii < pjc[start_node+1] ; ii++ ) {
            neighbor = pir[ii];
            if (pl[neighbor]==0) { // unlabeled
                pl[neighbor] = curr_label; // label it
                // push it into stack
                stack[sp] = neighbor;
                sp++;
            } else {
                if (pl[neighbor]!=curr_label)
                    mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"Got mixed labeling %d <-> %d",pl[neighbor], curr_label);
            }
        }
    }
    mxFree(stack);
}