x=[1.1,2.2,3.3,4.0,5.5,6.6,7.7]
y=[2.3,3.4,4.5,5.2,6.7,7.8,8.9]
yfit=[]
n=len(x)
xav=sum(x)/n
yav=sum(y)/n
sxy=sum([(i-xav)*(j-yav) for i,j in zip(x,y)])
sxx=sum([(i-xav)**2 for i in x])
a1=sxy/sxx
a0=yav-a1*xav
for i in range(n):
    yfit1=a0+a1*x[i]
    yfit.append(yfit1)

# print
import matplotlib.pyplot as plt
plt.scatter(x,y)
plt.plot(x,yfit)
plt.title('Least Square Fit',size=18)
plt.xlabel('x',size=18)
plt.ylabel('y',size=18)
plt.show()
