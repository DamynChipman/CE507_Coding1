
# Coding Assignment #1
### Damyn Chipman  -  CE 507  -  10.17.18

We start with the strong form of our problem. We wish to find the solution to:

\begin{equation}
\text{Find } u \text{ s.t.:} \\
u_{,xx} + f(x) = 0,\ \ \ x \in 0,1 \\
u(1) = 0, \ \ \ \ \ -u_{,x}(0) = 0
\end{equation}

Upon introduction of a trial solution space and weight function space, we can transition from the strong form to the weak form of our problem. From there, we introduce an approximation of both spaces and allow a summation of each to form the solution space of u. This allows us to write the Matrix form of our problem:

\begin{equation}
\sum_{B \in \eta \setminus \eta_g} d_B a(N_A(x),N_B(x)) = (N_A(x),f) + (N_A(x),h) - a(N_A(x), \sum_{B \in \eta_g} d_B N_B(x)) \\
\longrightarrow\ \ \sum_{B \in \eta \setminus \eta_g} d_B a(N_A(x),N_B(x)) = (N_A(x),f) - a(N_A(x), \sum_{B \in \eta_g} d_B N_B(x)) \\
\longrightarrow\ \ [K_{AB}] \{d_B\} = \{F_A\} \\
\text{Given } K \text{ and } \textbf{F} \text{, find } \textbf{d} \text{ s.t.:} \\
K \textbf{d} = \textbf{F}
\end{equation}

which is our Matrix form of our problem, where:

\begin{equation}
K_{AB} = a(N_A,N_B),\\
F_A = (N_A,f) + (N_A,h) - a(N_A,\sum_{B \in \eta_g} d_B N_B)
\end{equation}

(The derivation of this problem has already been explicitly expressed in previous homeworks.)

The goal is to solve for $\textbf{d}$ given global matrices $\textbf{K}$ and $\textbf{F}$. Both the global stiffness matrix $\textbf{K}$ and the global forcing vector $\textbf{F}$ can be built through an element stiffness matrix $\textbf{K}^e$ and an element forcing vector $\textbf{F}^e$ with appropiate integer mapping arrays.

To add another layer of robustness, we will create a routine that builds the global stiffness matrix and global forcing vector from inputs of the mapping arrays and element matrix and vector. In other words, we build a routine called *FE1D* such as:

\begin{equation}
FE1D(ID,\ IEN,\ LM,\ \textbf{K}^e,\ \textbf{F}^e,\ N_e) \longrightarrow \textbf{d}
\end{equation}

where *ID*, *IEN*, and *LM* are the mapping arrays and $\textbf{K}^e$ and $\textbf{F}^e$ are the element stiffness matrix and element forcing vector. This allows for a layer of robustness, as this routine will work for any 1D Linear finite element problem, given the specific inputs for the given problem. In addition, *FE1D* is built that it is easily adapted for higher order finite elements.

An outline of *FE1D* is given below:

\begin{align}
&\text{Given ID, IEN, LM, }\textbf{K}^e,\ \textbf{F}^e: \\
&\text{Initialize }\textbf{K} = \textbf{0}, \textbf{F} = \textbf{0} \\
&for\ \ e = 1,...,N_e \\
&\ \ \ \ for\ a = 1,...,N_{en} \\
&\ \ \ \ \ \ \ \ for\ b = 1,...,N{eq} \\
&\ \ \ \ \ \ \ \ \ \ \ \ calc\ k_{a,b}^e \\
&\ \ \ \ \ \ \ \ endfor \\
&\ \ \ \ calc\ f_a^e \\
&\ \ \ \ endfor \\
&\ \ \ \ for\ a = 1,...,N_{en} \\
&\ \ \ \ \ \ \ \ P = LM(a,e) \\
&\ \ \ \ \ \ \ \ if\ P \ne 0 \\
&\ \ \ \ \ \ \ \ \ \ \ \ for\ b = 1,...,N_{eq} \\
&\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ Q = LM(b,e) \\
&\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ if\ Q \ne 0 \\
&\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ K_{P,Q} = K_{P,Q} + k_{a,b}^e \\
&\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ endif \\
&\ \ \ \ \ \ \ \ \ \ \ \ endfor \\
&\ \ \ \ \ \ \ \ F_P = F_P + f_a^e \\
&\ \ \ \ \ \ \ \ endif \\
&\ \ \ \ endfor \\
&endfor \\
&\textbf{Solve }\textbf{K} \textbf{d} = \textbf{F} \text{ for }\textbf{d}
\end{align}

*FE1D* is implemented in C++ using Eigen for building the sparse matrices and solving the resulting equation. All numerical results are obtained in C++, and all visualization is in Python (this notebook).

We now turn our attention to the specific problem. We wish to solve the problem above using $N = 10,100,1000,10000$ (where N corresponds to nunmber of elements) and for $f = c,\ f = x,\ f = x^2$. The first task is to build the integer mapping arrays. The essential boundary conition at $x(1) = g = 0$ helps us build the ID array. Wherever we know the solution (i.e. @ $x=1$) there is no mapping (a $0$) between the node and global equation number. IEN is the mapping between the element node ID's to the global node ID's. LM is a mapping between element node ID's to the global equation number, or LM = ID(IEN). This gives us:

\begin{equation}
ID:
\begin{bmatrix}
\text{Node}     & 1 & 2 & ... & n & n+1 \\
\text{Equation} & 1 & 2 & ... & n & 0   \\
\end{bmatrix}
\end{equation}

\begin{equation}
IEN:
\begin{bmatrix}
(a,e) & e = 1 & e = 2 & ... & e = n \\
a = 1     & 1 & 2 & ... & n   \\
a = 2     & 2 & 3 & ... & n+1 \\
\end{bmatrix}
\end{equation}

\begin{equation}
LM:
\begin{bmatrix}
(a,e) & e = 1 & e = 2 & ... & e = n \\
a = 1     & 1 & 2 & ... & n \\
a = 2     & 2 & 3 & ... & 0 \\
\end{bmatrix}
\end{equation}

These are provided to *FE1D*. The element stiffness matrix and element forcing vectors can be found for specifically for linear elements. These are supplied as (derived in class and notes beginning from equations in Matrix form above):

\begin{equation}
\textbf{k}^e = \dfrac{1}{h^e} 
\begin{bmatrix}
1  & -1 \\
-1 & 1  \\
\end{bmatrix}
\ \ \ \textbf{f}^e = \dfrac{h^e}{6}
\begin{bmatrix}
2f_1^e + f_2^e \\
f_1^e + 2f_2^e \\
\end{bmatrix} + \textbf{BC}
\end{equation}

Our linear elements are defined as 

\begin{equation}
N_A(x) = \left\{
\begin{array}{ll}
      \dfrac{x - x_{A-1}}{x_A - x_{A-1}} & x_{A-1} \le x \le x_A \\
      \dfrac{x_{A+1} - x}{x_{A+1} - x_{A}} & x_A \le x \le x_{A+1} \\
      0 & otherwise
\end{array} 
\right.
\end{equation}

Once found using *FE1D*, the coefficients are used to form the solution, which is given as:

\begin{equation}
u^h = v^h + g^h \\
v^h = \sum_{B \in \eta \setminus \eta_g} d_B N_B(x) \\
g^h = \sum_{B \in \eta_g} d_B N_B(x) = 0 \\
\longrightarrow u(x) \approx d_1 N_1(x) + d_2 N_2(x) + ... + d_n N_n(x)
\end{equation}

There are three cases we wish to explore, for $f = c,\ f = x,\ f = x^2$. The analytical solutions for each are given by:

\begin{align}
u_A(x) &= \dfrac{1}{2} c (1 - x^2) \\
u_B(x) &= \dfrac{1}{6}(1 - x^3) \\
u_C(x) &= \dfrac{1}{12}(1 - x^4) \\
\end{align}


```python
import numpy as np
import matplotlib.pyplot as plt
import csv

N = [10,100,1000,10000]
var = 1
NPlot = [var*N[0]+1,var*N[1]+1,var*N[2]+1,var*N[3]+1]
file_path = 'CE507_Coding1/CE507_Coding1/outputs/'
Xfile_names = ['xN1.csv','xN2.csv','xN3.csv','xN4.csv']
U1file_names = ['u1N1.csv','u1N2.csv','u1N3.csv','u1N4.csv']
U2file_names = ['u2N1.csv','u2N2.csv','u2N3.csv','u2N4.csv']
U3file_names = ['u3N1.csv','u3N2.csv','u3N3.csv','u3N4.csv']
U1ACTfile_names = ['uAct1N1.csv','uAct1N2.csv','uAct1N3.csv','uAct1N4.csv']
U2ACTfile_names = ['uAct2N1.csv','uAct2N2.csv','uAct2N3.csv','uAct2N4.csv']
U3ACTfile_names = ['uAct3N1.csv','uAct3N2.csv','uAct3N3.csv','uAct3N4.csv']
U_file_names = [U1file_names,U2file_names,U3file_names]
UACT_file_names = [U1ACTfile_names,U2ACTfile_names,U3ACTfile_names]
titles = [['N = 10, f = c',
           'N = 10, f = x',
           'N = 10, f = x^2'],
          ['N = 100, f = c',
           'N = 100, f = x',
           'N = 100, f = x^2'],
          ['N = 1000, f = c',
           'N = 1000, f = x',
           'N = 1000, f = x^2'],
          ['N = 10000, f = c',
           'N = 10000, f = x',
           'N = 10000, f = x^2']]

plot_flag = 1
for n in [0,1,2,3]:
    # Create X to plot
    x = np.zeros(NPlot[n])
    with open(file_path+Xfile_names[n]) as csv_file:
        csv_reader = csv.reader(csv_file,delimiter=',')
        for row in csv_reader:
            for j in np.arange(0,NPlot[n]):
                x[j] = float(row[j])
            
    # Create U and UAct to plot
    for k in [0,1,2]:
        u = np.zeros(NPlot[n])
        with open(file_path+U_file_names[k][n]) as csv_file:
            csv_reader = csv.reader(csv_file,delimiter=',')
            for row in csv_reader:
                for j in np.arange(0,NPlot[n]):
                    u[j] = float(row[j])

        uAct = np.zeros(NPlot[n])
        with open(file_path+UACT_file_names[k][n]) as csv_file:
            csv_reader = csv.reader(csv_file,delimiter=',')
            for row in csv_reader:
                for j in np.arange(0,NPlot[n]):
                    uAct[j] = float(row[j])
                        
        # Plot results
        fig = plt.figure(1,figsize=(24,12))
        plt.subplot(4,3,plot_flag)
        plt.plot(x,uAct,'b',x,u,'r--')
        plt.xlim(0,1)
        plt.ylim(0,0.8)
        plt.grid(True)
        if plot_flag > 9: plt.xlabel("X")
        plt.ylabel("U")
        plt.title(titles[n][k])
        plt.legend(['Analytical','FE'])
        fig.suptitle("Overview of Solutions",fontsize=22)
        plot_flag = plot_flag + 1
```


![png](output_1_0.png)


As seen above, the FE solution agrees very well with the analytical solution. However, this trend fails after N = 1000 elements for reasons I am unable to determine. We can calculate the error between the analytical and FE solution through the norm of the difference between the functions. The norm for this function space is defined through an integral of the difference in functions over the domain, normalized for each value. We evaluate the integral using Gaussian Quadrature as shown below:

\begin{equation}
e = ||u - u^h|| = \left(\int_0^1 | u - u^h | ^2 dx \right)^{1/2} \\
= \left( \sum_e \int_{-1}^{1} | u(\xi) - u^h(\xi) |^2 \dfrac{\partial x(\xi)}{\partial \xi} d\xi \right)^{1/2} \\
\approx \left( \sum_e \sum_{i=1}^3 \int_{-1}^{1} | u(\xi_i) - u^h(\xi_i) |^2 \dfrac{\partial x(\xi_i)}{\partial \xi} w_i \right)^{1/2} \\
\text{where: } \xi_1 = -\sqrt{3/5},\ \ \xi_2 = 0,\ \ \xi_3 = \sqrt{3/5},\ \ w_1 = w_3 = 5/9,\ \ w_2 = 8/9
\end{equation}

These errors are plotted below in a log-log plot of the error versus the element spacing h.


```python
error_file_names = ['errorN1.csv','errorN2.csv','errorN3.csv','errorN4.csv']

del_x = np.zeros(4)
with open(file_path+'delX.csv') as csv_file:
    csv_reader = csv.reader(csv_file,delimiter=',')
    for row in csv_reader:
        for j in np.arange(0,4):
            del_x[j] = float(row[j])

error = np.zeros([4,3])
for n in [0,1,2,3]:
    with open(file_path+error_file_names[n]) as csv_file:
        csv_reader = csv.reader(csv_file,delimiter=',')
        for row in csv_reader:
            for j in np.arange(0,3):
                error[n][j] = float(row[j])


slope1 = (np.log(del_x[0]) - np.log(del_x[3]))/(np.log(error[0,0]) - np.log(error[3,0]))
slope2 = (np.log(del_x[0]) - np.log(del_x[3]))/(np.log(error[0,1]) - np.log(error[3,1]))
slope3 = (np.log(del_x[0]) - np.log(del_x[3]))/(np.log(error[0,2]) - np.log(error[3,2]))

# Plot error vs h
plt.figure(2,figsize=(10,8))
plt.plot(np.log(error[:,0]),np.log(del_x))
plt.plot(np.log(error[:,1]),np.log(del_x))
plt.plot(np.log(error[:,2]),np.log(del_x))
plt.legend(["f = constant","f = x","f = x^2"])
plt.title("Gaussin Quadrature Error for FE1D")
plt.xlabel("log(error)")
plt.ylabel("log(h)")
plt.text(-3,-4.4,"Slope of f = x: "+str(np.round(slope2,4)),fontsize=20)
plt.text(-3,-4.8,"Slope of f = x^2: "+str(np.round(slope3,4)),fontsize=20)
plt.grid(True)
```


![png](output_3_0.png)


Ideally, the error would decrease with the increase of N. The slope of the line would give the rate of convergance for the FE method to the exact solution. However, either through my implementation of the functions shown above or a miscalcualtion would cause for the FE solution to decay the larger that N gets, as opposed to converging to the analytical solution. I believe this has to do with how I calculate $\textbf{f^e}$ for each case of $f$. This is a place to start for debugging in the future.

### Appendix

Listing of the code for analysis is provided at https://github.com/camperD/CE507_Coding1

The speed of the program was greatly enhanced through the use of Eigen/Sparse rather than Eigen/Dense. Below are snapshots of terminal output for the two implementations:

<img src="evalTime.PNG" style="transform:scale(0.8)" />
