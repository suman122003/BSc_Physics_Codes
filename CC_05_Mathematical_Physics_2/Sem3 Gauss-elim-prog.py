def GaussElim(A):
	n=len(A)
	a=[[A[i][j] for j in range(len(A[0]))] for i in range(n)]
	for i in range(n-1):
		mxa=a[i][i]
		m=i
		for j in range(i+1,n):
			if abs(a[j][i])>mxa:
				mxa=abs(a[j][i])
				m=j
		ta=a[i]
		a[i]=a[m]
		a[m]=ta

		for j in range(i+1,n):
			cf=a[j][i]/a[i][i]
			for k in range(n+1):
				a[j][k]=a[j][k]-cf*a[i][k]
# Back substitution
	X=[0.0 for i in range(n)]
	X[n-1]=a[n-1][n]/a[n-1][n-1]
	for i in range(n-2,-1,-1):
		sm=0.0
		for j in range(i+1,n):
			sm+=a[i][j]*X[j]	
		X[i]=1.0/a[i][i]*(a[i][n]-sm)
	return X
xyz=[[-5,16,-4,0],[10,-5,0,12],[0,-4,11,0]]
X=GaussElim(xyz)
print('coeff=',X)
 
