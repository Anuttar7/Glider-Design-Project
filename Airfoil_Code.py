import numpy as np
import matplotlib.pyplot as plt

#Main function where all the used variables are stored, and functions are called to
#compute required computations
def Main():
#variable declaration
    
    N = 1001    #Number of points on camber line
    M = 200    #Number of fourier coefficients of gamma considered for computation
    
    theta = np.linspace(0, np.pi, N)    #stores the parametric coordinate theta values
    x = 0.5*(1-np.cos(theta))           #stores the x-coordinates of camber-line
    y = 0                               #stores the y-coordinates of camber-line
    dy = 0
    
    alpha = 0                           #angle of attack
    v0 = 0                              #Free stream velocity
    m = 0                               #maximum camber
    p = 0                               #position of maximum camber
    q = 0                               #flag variable to tell airfoil is NACA or user defined
    
    An = 0                              #Fourier coefficients of gamma function
    gamma = 0                           #gamma function
    cl = 0                              #coefficient of lift
    
    
#Input variable values from user
    
    print("Enter 1 if airfoil is NACA, 2 if user defined: ")
    q = int(input())
    print("Enter max camber: ")
    m = float(input())
    if (q == 1):
        print("Enter position of max camber: ")
        p = float(input())
    else:
        print("Enter 1: parabolic, 2: Elliptical, 3: Hyperbolic, 4: Sinusoidal, 5: Circular Arc")
        qq = int(input())
        p = 0.5
    print("Enter Angle of Attack in degrees: ")
    alpha = float(input())
    print("Enter freestream velocity: ")
    v0 = float(input())
    
    
#Computations

    alpha = alpha*(np.pi/180)               #converting from deg to radians
    
    #to compute and store the y-coordinates if camber is NACA or user-defined
    if q == 1:
        y = f_NACA(x, m, p)
    else:
        y = f_user(x, m, qq)
    
    dy = derivative(x, y)                   #to store the slope at all camber-line points    
    An = an(dy, theta, N, M, alpha)         #An are computed and stored from n=0 to n=M-1
    gamma = gamma_func(theta, An, N, M, v0) #gamma is stored for all theta
    cl = Cl(An)                             #Cl is computed
    
    
#Outputs to user
    
    #camberLinePlotter(x, y)
    #camberLineSlopePlotter(x, dy)
    #slopeFindingFunc(x, dy)
    #print("Cl of the airfoil is: ", cl)
    #Cl_vs_Alpha(dy, theta)
    vectorPlot(theta, x, y, gamma, v0, alpha)
    #print("The total Circulation is: ", circulation_gammaIntegral(gamma, theta))
    #circulation_VelocityLineIntegral(x, y, gamma, theta, v0, alpha)
#______________________________________________________________________________

#function to compute the y-coordinates of NACA airfoil with m, p parameters
def f_NACA(x, m, p):
    n = len(x)                      #n = N
    y = np.zeros(n)                 #zero initialised array of size n
    
    #Computes the y coordinates
    for i in range(n):              
        if (x[i] <= p):
            y[i] = (m/p**2)*(2*p*x[i] - x[i]*x[i])
        else:
            y[i] = (m/(1-p)**2)*(1-2*p + 2*p*x[i] - x[i]*x[i])
    
    return y                        
#______________________________________________________________________________    

#User-defined function to compute y-coordinates of camber line
def f_user (x, m, q):
    if (q == 1):
        y = 4*m*x*(1-x)                                     #Parabolic
    elif (q == 2):
        y = 2*m*((x - x*x)**0.5)                            #Elliptical
    elif (q == 3):
        y = (m+1) - ( 1 + (m*m + 2*m)*(1 - 2*x)**2 )**0.5   #Hyperbolic
    elif (q == 4):
        y = m*np.sin(np.pi*x)                               #Sinusoidal
    else:
        k = (1 - 4*m*m)/(8*m)
        y = -k + (k*k + x*(1-x))**0.5                       #Circular Arc
    
    return y
#______________________________________________________________________________

#function to compute the slope of all camber-line points
def derivative(x, y):
    n = len(x)                      #n = N
    dy = np.zeros(n)                #zero initialised array of size n
    
    dy[0] = (y[1]-y[0])/(x[1]-x[0])
    for i in range(1, n):
        dy[i] = (y[i-1] - y[i])/(x[i-1] - x[i])
    return dy
#______________________________________________________________________________

#function to compute the An values for all n=0 to n=M-1
def an(dy, theta, N, M, alpha):
    An = np.zeros(M)
    for i in range(M):
        r = 0
        
        if (i == 0):
            for j in range(N):
                r = r + dy[j]
            r = alpha - (1/np.pi)*r*(np.pi/N)
            An[i] = r
        else:
            for j in range(N):
                r = r + dy[j]*np.cos(i*theta[j])
            r = (2/np.pi)*r*(np.pi/N)
            An[i] = r
    return An
#______________________________________________________________________________

#function to compute the Cl value of the airfoil at given angle of attack
def Cl(An):
    cl = np.pi*(2*An[0] + An[1])
    return cl
#______________________________________________________________________________

#function to compute the gamma function values for all theta
def gamma_func(theta, An, N, M, v0):
    gamma = np.zeros(N)
    for i in range(1, N-1):
        r = An[0]*(1+np.cos(theta[i]))/np.sin(theta[i])
        for j in range(1, M):
            r = r + An[j]*np.sin(j*theta[i])
        r = 2*v0*r
        gamma[i] = r
    return gamma
#______________________________________________________________________________

#Function to plot the velocity vectors in the region around the airfoil
def vectorPlot(theta, x, y, gamma, v0, alpha):
    N = len(x)
    
    #a grid of size m x n is created containing (x, y) points
    n = 20
    m = 30
    X, Y = np.meshgrid(np.linspace(0.5-2, 0.5+2, n), np.linspace(-1.5, 1.5, m))
    
    #u stores the x-component of induced velocity
    #v stores the y-component of induced velocity
    #u and v are zero-initialised below
    u = np.zeros_like(X)
    v = np.zeros_like(Y)
    
    #i, j are used to select different points on the grid
    for i in range(m):
        for j in range(n):
            #k is used to select a point on the chord
            for k in range(N):
                #in this loop, the velocity induced at point (X[i][j], Y[i][j])
                #by the entire chord is calculated via numerical integration
                
                #r gives the distance between 2 points
                r =  ((Y[i][j] - y[k])**2 + (X[i][j] - x[k])**2)**0.5
                
                #handling cases giving 0 denominators individually
                #Case 1: r = 0
                if r == 0:
                    continue
                
                K = gamma[k]/(2*np.pi*r)
                dzeta = 0.5*np.sin(theta[k])*(np.pi/N)
                
                #Case 2: delta = pi/2 (i.e slope not defined)
                if X[i][j] - x[k] == 0:
                    if Y[i][j] > y[k]:
                        u[i][j] = u[i][j] + K*dzeta
                    else:
                        u[i][j] = u[i][j] - K*dzeta
                    continue
                
                cos_delta = (X[i][j] - x[k])/r
                sin_delta = (Y[i][j] - y[k])/r
                
                u[i][j] = u[i][j] + K*sin_delta*dzeta
                v[i][j] = v[i][j] - K*cos_delta*dzeta
            
            #total velocity is computed by adding the freestream velocity
            u[i][j] = u[i][j] + v0*np.cos(alpha)
            v[i][j] = v[i][j] + v0*np.sin(alpha)
    
    camberLinePlotter(x,y)          #the camber line is plotted again for reference
    plt.quiver(X,Y,u,v)             #Plots the vector plot
#______________________________________________________________________________
    
#Function to find the total ciculation using integral of gamma along airfoil
def circulation_gammaIntegral (gamma, theta):
    r = 0
    n = len(gamma)
    for i in range(n):
        r += gamma[i]*np.sin(theta[i])
    r = r*(np.pi/n)/2
    return r
#______________________________________________________________________________

#Function to find the total circulation using velocity line integral around a
#circle containing the airfoil
def circulation_VelocityLineIntegral (x, y, gamma, theta, v0, alpha):
    n = len(gamma)                      #Number of points on camber line
    N = 10000                            #Number of points on curve
    beta = np.linspace(0,2*np.pi, N)    #parametric coordinate beta
    x1 = 0.5 + np.cos(beta)             #x-component of point on curve
    y1 = np.sin(beta)                   #y-component of point on curve
    u = np.zeros(N)                     #x-component of velocity on curve
    v = np.zeros(N)                     #y-component of velocity on curve
    circ = 0                            #the total circulation
    
    for i in range(N):
        #each i value represents a point on curve
        for j in range(n):
            #each j value represents a point on camber-line
            #In this loop, the induced velocity at point beta[i] on curve is calculated
            
            #distance of point on curve and point on camber line
            r =  ((y1[i] - y[j])**2 + (x1[i] - x[j])**2)**0.5
            
            #Cases to handle zero-division error
            #Case 1: r = 0
            if r == 0:
                continue
            
            dzeta = 0.5*np.sin(theta[j])*(np.pi/N)
            K = gamma[j]/(2*np.pi*r)
            
            #Case 2: delta = pi/2 (i.e slope is undefined)
            if x1[i] - x[j] == 0:
                if y1[i] > y[j]:
                    u[i] = u[i] + K*dzeta
                else:
                    u[i] = u[i] - K*dzeta
                continue
            
            cos_delta = (x1[i] - x[j])/r
            sin_delta = (y1[i] - y[j])/r
            
            u[i] = u[i] + K*sin_delta*dzeta
            v[i] = v[i] - K*cos_delta*dzeta
        
        #total velocity is calculated by summing the free stream velocity
        u[i] = u[i] + v0*np.cos(alpha)
        v[i] = v[i] + v0*np.sin(alpha)
            
        #an infinitely small length along the curve is calculated
        dsx = -1*(2*np.pi/N)*np.sin(beta[i])
        dsy = (2*np.pi/N)*np.cos(beta[i])
        
        #dot product is taken between total velocity and the infinitely small
        #curve length and the result gives a contribution to total circulation
        circ = circ + u[i]*dsx + v[i]*dsy
    
    # the integral is calculated opposite to the velocity direction and hence, 
    # will always return a negative result. Hence, -circ is returned
    circ = -1*circ
    print("Circulation around circle is: ", circ)
    return
#______________________________________________________________________________

#Function to plot the camber line of airfoil
def camberLinePlotter(x, y):
    plt.xlabel('x/c')
    plt.ylabel("Camber Line")
    plt.xlim(0.5-2, 0.5+2)
    plt.ylim(-1.5, 1.5)
    plt.title("Camber Line of Airfoil")
    plt.plot(x, y)
#______________________________________________________________________________

#Function to find the slope of a given input point on camber line
def slopeFindingFunc(x, dy):
    print("Enter position at which you want to know the slope: ")
    px = float(input())
    n = len(x)
    for i in range(n-1):
        if (x[i] <= px and x[i+1] > px):
            print ("Slope at the point is: ", dy[i])
            return
    if x[n-1]==px:
        print ("Slope at the point is: ", dy[i])
    else:
        print("Given point is not on airfoil")
    return
#______________________________________________________________________________

#Function to plot the slope vs x curve of camber line
def camberLineSlopePlotter(x, dy):
    plt.xlabel('x/c')
    plt.ylabel("Camber Line Slope")
    plt.title("Slope of the Camber Line of Airfoil")
    plt.plot(x, dy)
    plt.axis('equal')
#______________________________________________________________________________

#Function to plot the Cl vs alpha curve of the airfoil
def Cl_vs_Alpha(dy, theta):
    n = len(theta)
    alpha = np.linspace(-5, 20, 26)
    alpha = alpha*np.pi/180
    cl = 0
    a0 = 0
    a1 = 0
    
    r = 0
    for j in range(n):
        r += dy[j]
    r = r*(np.pi/n)
    a0 = np.sin(alpha) - (np.cos(alpha)/np.pi)*r
        
    r = 0
    for j in range(n):
        r += dy[j]*np.cos(theta[j])
    r = r*(np.pi/n)
    a1 = (np.cos(alpha)/np.pi)*r
    
    cl = np.pi*(2*a0 + a1)
    
    plt.title("Cl vs AoA")
    plt.xlabel("Angle of Attack (in deg)")
    plt.ylabel("Cl")
    plt.plot(alpha*180/np.pi, cl)
#______________________________________________________________________________

#Function to find the Cm0 value about the leading edge of the airfoil
def Cm (An):
    cm = (-1*np.pi/2)*(An[0] + An[1] - An[2]/2)
    return cm
#______________________________________________________________________________

#calling the main function
Main()