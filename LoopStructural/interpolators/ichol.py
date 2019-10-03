# import numpy as np
# from scipy.sparse import tril
# def ichol(A):
#     LA_row = []
#     LA_col = []
#     LA_val = []
#     n = A.shape[0]
#     #set diagonal elements to be equal to the sqrt of the A
#     LA_val.extend(np.sqrt(A.diagonal()).tolist())
#     LA_col.extend(np.arange(0,A.shape[0]).tolist())
#     LA_row.extend(np.arange(0,A.shape[0]).tolist())
# 
#     #return only the lower triangle
#     lower = tril(A,-1)
#     #set diagonal to sqrt of A
#     lower.setdiag(np.sqrt(A.diagonal()))
#     for k in range(0,n):
#         #divide column elements by diagonal element
#         lower.getcol(k)
# 
#     for k in range(0,n):
#         A[k,k] = np.sqrt(A[k,k])
#         for i in range(k+1,n):
#             if
