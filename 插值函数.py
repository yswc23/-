import numpy as np
from numpy.polynomial import Polynomial
from sympy import *

#计算多项式的值
def poly(lst, x):
    sum = 0
    for i in range(len(lst)):
        sum += lst[i] * x ** i
    return sum

#Lagrange插值函数
def Lagrange(x, y, x0):
    c = np.zeros(len(x))
    for i in range(len(x)):
        a = x[0]
        b = 1
        x.remove(a)
        L = np.poly(x)
        for j in range(len(x)):
            b *= a - x[j]
        x.append(a)
        L = y[i] * L / b
        c += L
    d = np.zeros(len(x))
    j = 0
    for i in c:
        d[j] = c[len(x)-j-1]
        j += 1
    L_x = Polynomial(coef=d)
    value = poly(d, x0)
    return L_x, value

#Newton插值函数
def Newton(x, y,x0):
    a = []
    a.append([])
    if len(x) <= 2:
        for i in range(len(x)-1):
            cs_1 = (y[i+1] - y[i]) / (x[i+1] - x[i])       #计算一阶差商
            a[0].append(cs_1)
    else:
        for i in range(len(x)-1):
            cs_1 = (y[i+1] - y[i]) / (x[i+1] - x[i])       #计算一阶差商
            a[0].append(cs_1)
        for j in range(len(x)-2):
            a.append([])
            for i in range(len(x)-j-2):
                cs = (a[j][i+1] - a[j][i]) / (x[i+2+j] - x[i])        #计算高阶差商
                a[j+1].append(cs)
    b = []
    for i in range(len(x)):
        if i == 0:
            b.append(y[0])
        else:
            b.append(a[i-1][0])
    c = np.zeros(len(x))
    d = []
    for i in range(len(x)):
        if i == 0:
            c[i] = b[i]
        else:
            d.append(x[i-1])
            e = np.poly(d)
            e = e.tolist()
            e.reverse()
            while len(e) < len(c):
                e.append(0)
            e = np.array(e)
            c += b[i] * e
    N_x = Polynomial(coef=c)
    value = poly(c, x0)
    return N_x, value

if __name__=='__main__':
    #输入计算数据：
    x = [2.0, 2.1, 2.2]
    y = [1.414214, 1.449138, 1.483240]
    x0 = 2.15

    fun_L, value_L = Lagrange(x, y, x0)
    print('Lagrange插值计算结果：\nf(x) = {}\nf({}) = {}'.format(fun_L, x0, value_L))      #输出插值函数
    fun_N, value_N = Newton(x, y, x0)
    print('Newton插值计算结果：\nf(x) = {}\nf({}) = {}'.format(fun_N, x0, value_N))  # 输出插值函数
    print('函数在x = 2.15处的真实值：\nsqrt(2.15) = ', format(sqrt(2.15)))