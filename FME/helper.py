import numpy as np

def rotation(axis,angle):
    c = np.cos(np.deg2rad(angle))
    s = np.sin((np.deg2rad(angle)))
    C = 1.0 - c
    x = axis[0]
    y = axis[1]
    z = axis[2]
    xs = x*s
    ys = y*s
    zs = z*s
    xC = x*C
    yC = y*C
    zC = z*C
    xyC = x*yC
    yzC = y*zC
    zxC = z*xC
    rotation_mat = np.zeros((3,3))
    rotation_mat[0][0] = x*xC+c
    rotation_mat[0][1] = xyC-zs
    rotation_mat[0][2] = zxC+ys
    
    rotation_mat[1][0] = xyC+zs
    rotation_mat[1][1] = y*yC+c
    rotation_mat[1][2] = yzC-xs

    rotation_mat[2][0] = zxC -ys
    rotation_mat[2][1] = yzC+xs
    rotation_mat[2][2] = z*zC+c
    return rotation_mat
    for i in range(3):
        for j in range(3):
            rv[i]+=rotation_mat[i][j]*vector[i]
        n+=rv[i]*rv[i]
    n = sqrt(n)
    for i in range(3):
        rv[i]/=n
    return rv
def mg2coords(X, Y, Z):
        return np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
def section2d(azi,dip,distance,origin=np.ones(3),scalev=np.ones(3),res=100):
    #create a unit square section
    xvalues = np.linspace(-1,1,res)
    yvalues = np.linspace(-1,1,res)
    zvalues = np.array([0.])
    
    X,Y,Z =  np.meshgrid(xvalues, yvalues, zvalues)
    
    #apply the dip rotation around y axis
    rot = rotation(np.array([0.,1.,0.]),dip)
    #apply the azi rotation around z
    rot2 = rotation(np.array([0.,0.,1.]),azi)
    scale = np.zeros((4,4))
    for i in range(3):
        scale[i,i] = scalev[i]
    xyz = mg2coords(X,Y,Z)
    xyz2=np.einsum('ij,kj->ki',rot,xyz)
    xyz3=np.einsum('ij,kj->ki',rot2,xyz2)
    xyz3 = np.c_[xyz3,np.ones(xyz3.shape[0])]
    xyz4=np.einsum('ij,kj->ki',scale,xyz3)
    xyz4=xyz4[:,:3]
    norm= np.zeros(3)
    norm[0] = np.cos(np.deg2rad(azi))*np.sin(np.deg2rad(dip))
    norm[1] = -np.sin(np.deg2rad(azi))*np.sin(np.deg2rad(dip))
    norm[2] = np.cos(np.deg2rad(dip))
    
    norm*=distance
    
    xyz4+=norm
    

    return xyz4
#print(xx.shape)

#plt.plot(xx[:,:,0], yy[:,:,0], marker='.', color='k', linestyle='none')


# In[4]:


def add_value(A,B,x,y,i,r,s=1.0):
    A[r][4*i+0] = s*x**3
    A[r][4*i+1] = s*x**2
    A[r][4*i+2] = s*x
    A[r][4*i+3] = s*1.0
    B[r] = y
def add_dx(A,B,x,y,i,r,s=1.0):
    A[r][4*i+0] = s*3*x**2
    A[r][4*i+1] = s*2*x
    A[r][4*i+2] = s*1.0
    A[r][4*i+3] = s*0.0
    B[r] = y
def add_ddx(A,B,x,y,i,r,s=1.0):
    A[r][4*i+0] = s*6*x
    A[r][4*i+1] = s*2
    A[r][4*i+2] = 0.0
    A[r][4*i+3] = 0.0
    B[r] = y
def func(w,x):
    return w[0]*x**3+w[1]*x**2+w[2]*x+w[3]

def get_weights2(z):
    A = np.zeros((8,8))
    B = np.zeros(8)
    #set control points, 0 and 1
    add_value(A,B,x=-1,y=0,i=0,r=0)
    #add_value(A,B,x=1,y=0,i=0,r=1)
    add_value(A,B,x=-.25,y=1,i=0,r=1)
    #set gradient at edges
    add_dx(A,B,x=-0.5,y=.5,i=0,r=2)
    add_dx(A,B,x=0,y=0,i=0,r=3)
    #add_value(A,B,x=-0.8,y=0.8,i=0,r=2)

    w1 = np.linalg.lstsq(A,B)[0]
    return w1
def get_weights3(z):
    A = np.zeros((8,8))
    B = np.zeros(8)
    #set control points, 0 and 1
    add_value(A,B,x=-1,y=1,i=0,r=0)
    add_value(A,B,x=0,y=0,i=0,r=1)
    #add_value(A,B,x=-.5,y=1,i=0,r=1)
    #set gradient at edges
    add_dx(A,B,x=-1,y=0,i=0,r=2)
    add_dx(A,B,x=0,y=-1,i=0,r=3)
    #add_value(A,B,x=-0.8,y=0.8,i=0,r=2)

    w1 = np.linalg.lstsq(A,B)[0]
    return w1
def get_weights(z):
    A = np.zeros((8,8))
    B = np.zeros(8)
    add_value(A,B,x=-1,y=0,i=0,r=0)
    add_dx(A,B,x=-1,y=0,i=0,r=1)
    add_value(A,B,x=0,y=1,i=0,r=2)
    add_dx(A,B,x=0,y=0,i=0,r=3)
    add_value(A,B,x=-z,y=0,i=0,r=4)
    add_value(A,B,x=-z,y=0,i=1,r=4,s=-1.0)
    add_dx(A,B,x=-z,y=0,i=0,r=5)
    add_dx(A,B,x=-z,y=0,i=1,r=5,s=-1.0)
    #add_dx(A,B,x=0,y=-z,i=1,r=7,s=1.0)
    #add_ddx(A,B,x=0,y=0,i=1,r=6,s=1.0)
    #add_ddx(A,B,x=0,y=0,i=0,r=7,s=1.0)
    #add_ddx(A,B,x=0,y=0,i=1,r=6,s=1.0)
    add_value(A,B,x=0,y=0,i=1,r=6,s=1.0)
    #add_ddx(A,B,x=-1,y=0,i=0,r=7,s=1.0)
    #add_value(A,B,x=-1,y=0,i=0,r=7,s=1.0)
    add_dx(A,B,x=0,y=-1,i=1,r=7)

    w1 = np.linalg.lstsq(A,B)[0]
    return w1
def p5(x):
    y = np.copy(x)
    y[:] = 0
    y[x>0] = 1
    return y
def p4(x):
    w1 = get_weights3(0.)
    x2 = np.copy(x)
    
    x2[x>0] = -x2[x>0]
    y = func(w1[:4],x2)
    y[x>0] = -y[x>0]
    return y
def p1(x):
    w1 = get_weights2(0.)
    x2 = np.copy(x)
    
    x2[x>0] = -x2[x>0]
    y = func(w1[:4],x2)
    y[x2>-0.5] = 1.
    return y
def p0(x):
    w1 = get_weights(0)
    w2 = get_weights(0)
    x2 = np.copy(x)
    x2[x>0] = -x2[x>0]#if x>0:
    y = np.zeros(x.shape)
    y[np.logical_and(x>=-1,x<=0)] = func(w1[:4],x[np.logical_and(x>=-1,x<=0)])
    #y[np.logical_and(x>=-z,x<=0)] = 1#func(w1[4:],x[np.logical_and(x>=-z,x<=0)])
    y[np.logical_and(x>=0,x<=1)] = -func(w1[:4],-x[np.logical_and(x>=0,x<=1)])

    #yplt.plot(-x[np.logical_and(x<=1,x>=z)],-func(w2[:4],x[np.logical_and(x<=1,x>=z)]))
    #plt.plot(-x[np.logical_and(x<=z,x>=0)],-func(w2[4:],x[np.logical_and(x<=z,x>=0)]))
    return y
def p2(x):
    return np.ones(x.shape)
def fault_func2(p0,p1,p2,g0,g1,g2):
    return p0(g0)*p1(g1)*p2(g2)
def normalz(gx):
    gxn = np.zeros(gx.shape)
    gxn = (2./(np.max(gx[~np.isnan(gx)])-np.min(gx[~np.isnan(gx)])))
    gxn*=(gx-((np.min(gx[~np.isnan(gx)])+np.max(gx[~np.isnan(gx)]))/2.))
    gxn[np.isnan(gx)] = np.nan
    return gxn

def savePoints(filename,points,data):
    from pyevtk.hl import gridToVTK, pointsToVTK
    pointsToVTK('../data/'+filename,points[0,:],points[1,:],points[2,:],data)
def strike_dip_vector(strike,dip):
    vec = np.zeros((len(strike),3))
    s_r = np.deg2rad(strike)
    d_r = np.deg2rad(np.abs(dip))
    vec[:,0] = np.sin(d_r)*np.cos(s_r)
    vec[:,1] = -np.sin(d_r)*np.sin(s_r)
    vec[:,2] = np.cos(d_r)
    vec  /= np.linalg.norm(vec,axis=1)[:,None]
    return vec
def array_to_vtk_points(points,filename):
    from pyevtk.hl import pointsToVTK
    x = np.array(points[:,0], copy=True, order='C')
    y = np.array(points[:,1], copy=True, order='C')
    z = np.array(points[:,2], copy=True, order='C')
    v = np.array(points[:,3], copy=True, order='C')
    if points.shape[1] > 4:
        v = (np.array(points[:,3], copy=True, order='C'),
             np.array(points[:,4], copy=True, order='C'),
             np.array(points[:,5], copy=True, order='C'))
    pointsToVTK(filename,x,y,z,data={"v":v})