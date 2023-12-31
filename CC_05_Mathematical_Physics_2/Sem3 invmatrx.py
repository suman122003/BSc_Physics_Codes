
# Matrix Inversion ******************


def InvGaussElim(A):
    n=len(A)
    a=[[A[i][j] for j in range(n)] for i in range(n)]
    b=[[1.0 if i==j else 0.0 for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if j!=i:
                r=a[j][i]/a[i][i]
                for k in range(n):
                    a[j][k]=a[j][k]-r*a[i][k]
                    b[j][k]=b[j][k]-r*b[i][k]
    for i in range(n):
        for j in range(n):
            b[i][j]=b[i][j]/a[i][i]
    return b



A=[[2,5,6,8],[7,10,5,4],[1,3,4,8],[4,8,9,12]]
invA=InvGaussElim(A)
print("         Inverse of the Matrix A     ")
for i in invA:
    for ii in i:
        print("%0.3f,"%ii,end=" ")   # all ii in same line
    print("")
