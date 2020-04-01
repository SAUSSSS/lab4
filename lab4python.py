"""
Author: Ivanov Aleksandr 426g

"""





from numpy import linalg as LA

n = 20

def norma_matrix(a):
    max,summa = 0,0
    for j in range(n):
        for i in range(n):
            summa += a[i][j]
        if(max < summa):
            max = summa
        summa = 0
    return max


def vector_discrepancy(a,x,f):
        temp = []
        R = []
        otv = 0
        for i in range(n):
            R.append(0)
            temp.append(0)
        temp[0] = 1 * x[0]
        for i in range(1,n-1,1):
            temp[i] += (1 * x[i-1]) + (-2 * x[i]) + (1 * x[i+1])
        otv += 1 * x[0] + 1 * x[n-1]
        for i in range(1,n-2,1):
            otv += (2 * x[i])
        temp[n-1] = otv
        for i in range(n-1):
            R[i] = f[i] - temp[i]
        new_otv = 0
        max = 0
        for i in range(n-1):
            new_otv += R[i]
            if(max < new_otv):
                max = new_otv
                new_otv = 0
        return max


def gauss(a,f):
    x = []
    for i in range(n):
        x.append(0)
    k,index,maxx = 0,0,0
    eps = 0.00001
    while(k < n):
        maxx = abs(a[k][k])
        index = k
        for i in range(k+1,n,1):
            if(abs(a[i][k]) > maxx):
                maxx = abs(a[i][k])
                index = i
        if(maxx < eps):
            print("No solution")
        for j in range(n):
            temp1 = a[k][j]
            a[k][j] = a[index][j]
            a[index][j] = temp1
        temp2 = f[k]
        f[k] = f[index]
        f[index] = temp2
        for i in range(k,n):
            temp3 = a[i][k]
            if(abs(temp3) < eps): 
                continue
            for j in range(n):
                a[i][j] /= temp3
            f[i] /= temp3
            if(i == k):
                continue
            for j in range(n):
                a[i][j] -= a[k][j]
            f[i] -= f[k]
        k += 1
    for t in range(n-1,-1,-1):
        x[t] = f[t]
        for i in range(t):
            f[i] -= a[i][t] * x[t]
    return x

w = 1.15
iterr = 0 
                  
def ros(a,f):
    global iterr
    eps = 0.00001
    xn = []
    xx = []
    for i in range(n):
        xn.append(0)
        xx.append(0)
    norma = 0
    for t in range(1000):
        norma = 0
        for i in range(n):
            xx[i] = f[i]
            for j in range(n):
                if(i != j):
                    xx[i] = xx[i] - a[i][j] * xx[j]       
            xx[i] /= a[i][i]
            xx[i] = w * xx[i] + (1 - w) * xn[i]
            if(abs(xx[i] - xn[i]) > norma):
                norma = abs(xx[i] - xn[i])
            xn[i] = xx[i]
        iterr += 1
        if(norma < eps):
            break
    return xn
            
        



a = []
ca = []
x = []
f = []
y = []
cf = []

for i in range(n):
    a.append([])
    ca.append([])
    f.append(0)
    for j in range(n):
        a[i].append(0)
        ca[i].append(0)
        f.append(0)
        cf.append(0)
f[0] = 1
v = 1
for i in range(1,n-1,1):
    a[i][v-1] = 1
    a[i][v] = -2
    a[i][v+1] = 1
    f[i] = 2 / pow(i,2)
    v += 1
a[0][0] = 1     
a[n-1][0] = 1  
for i in range(1,n-1,1):
    a[n-1][i] = 2
a[n-1][n-1] = 1
f[n-1] =  -n/3

for i in range(n):
    cf[i] = f[i]
    for j in range(n):
        ca[i][j] = a[i][j]


for i in range(n):
    for j in range(n):
        print(a[i][j], end = " ")
    print(" = ",f[i])
print()
      
a_inv = LA.inv(a)
lambda_min = norma_matrix(a_inv)
lambda_max= norma_matrix(a)
cond = lambda_max * lambda_min
print("Lambda min:",lambda_min)
print("Lambda max:",lambda_max)
print("Condition:", cond)





x = gauss(a,f)
y = ros(ca,cf)


print("=====GAUSS=======")
for i in range(n):
    print(x[i])
vector_gauss = 0
vector_gauss = vector_discrepancy(a,x,cf)
print()
print("=====ROS===","w=", w ,"Iter:", iterr )
for i in range(n):
    print(y[i])
vector_ros = 0
vector_ros = vector_discrepancy(a,y,cf)