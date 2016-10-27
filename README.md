# hierarchical_clustering

This function optimizes the following problem:

    \min_{W,Qi,F\in Ind} ||X'W+1b'-F||_F^2 + gamma*\sum_i^g||WQ_i||_{2,1}^2 + lambda*Tr(F'11'F)

Format of input:
    
    n: number of samples
    d: number of SNPs
    c: number of clusters
    X: d*n SNP data
    numC: number of clusters
    numG: number of groups of the clusters

Simply run the code in matlab as below:

    [outF, outQ, outW, outObj, outNumIter] = hiCluster(X, numC, numG);
