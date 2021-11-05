# Geometrical conceptual models for lateral and vertical extent of intrusions
import numpy as np
import pandas as pd

def ellipse_function(lateral_contact_data, test = 0, minP=None, maxP=None, minS=None, maxS=None):
    import math
    if minP==None: 
        minP = lateral_contact_data['coord1'].min()
    if maxP==None: 
        maxP = lateral_contact_data['coord1'].max()
    if minS==None: 
        minS = lateral_contact_data['coord2'].abs().min()
    if maxS==None: 
        maxS = lateral_contact_data['coord2'].max()
            
    a = (maxP-minP)/2
    b = (maxS-minS)/10
#     b = abs(minS)

    po = minP + (maxP-minP)/2
    
    p_locations = lateral_contact_data.loc[:,'coord1'].copy().to_numpy()
#     s_values = lateral_contact_data.loc[:,'coord2'].copy().to_numpy()
    
    s = np.zeros([len(p_locations),2])
    
    for i in range(len(p_locations)):
        if minP < p_locations[i] < maxP:
            s[i,0] = b*math.sqrt(1-pow((p_locations[i] - po)/a , 2)) # max side
            s[i,1] = -b*math.sqrt(1-pow((p_locations[i] - po)/a , 2)) # min side
#         elif test == 1:
#             s[i,:] = lateral_contact_data.loc[i,'coord2']
        else:
            s[i,:] = 0

    return s

def rectangle_function(lateral_contact_data, minP=None, maxP=None, minS=None, maxS=None):
    import math
    if minP==None: 
        minP = lateral_contact_data['coord1'].min()
    if maxP==None: 
        maxP = lateral_contact_data['coord1'].max()
    if minS==None: 
        minS = lateral_contact_data['coord2'].min()
    if maxS==None: 
        maxS = lateral_contact_data['coord2'].max()
            
    p_locations = lateral_contact_data.loc[:,'coord1'].copy().to_numpy()
    s = np.zeros([len(p_locations),2])
    
    for i in range(len(p_locations)):
        if minP < p_locations[i] < maxP:
            s[i,0] = maxS # max side
            s[i,1] = minS # min side
        else:
            s[i,:] = 0

    return s

def parallelepiped_function(othercontact_data, mean_growth=None, minP=None, maxP=None, minS=None, maxS=None, vertex=None):
    
    if mean_growth == None:
        mean_growth = othercontact_data.loc[:,'coord1'].mean()
        
    data_ps = np.array([othercontact_data.loc[:,'coord1'], othercontact_data.loc[:,'coord2']]).T
    
    conceptual_growth = np.ones([len(data_ps),2]) * mean_growth
    
    return conceptual_growth


def obliquecone_function(othercontact_data, mean_growth=None, minP=None, maxP=None, minS=None, maxS=None, vertex=None): 
    import math
    
    ps_locations = othercontact_data.loc[:,['coord1','coord2']].to_numpy()
    
    a = (maxP-minP)/2 #semi-major axis
    b = (maxS-minS)/2 #semi-minor axis
    a2 = pow(a,2)
    b2 = pow(b,2)
    
    po = minP + a #p coordinate of ellipsis centre
    so = minS + b #s coordinate of ellipsis centre
    
    alpha = vertex[0] #p coordinate of vertex
    beta = vertex[1] #g coordinate of vertex
    gamma = vertex[2] #s coordinate of vertex
    
    growth = np.zeros([len(ps_locations),2]) #container for results
    
    for i in range(len(ps_locations)):
        p = ps_locations[i,0]
        s = ps_locations[i,1]
        
        A = alpha - po
        B = beta*(p-alpha)
        C = gamma - so
        D = beta*(s-gamma)
        
        F = pow(A*b,2) + pow(C*a,2) - a2*b2
        G = 2*(B*A*b2 + C*D*a2)
        H = pow(b*B,2) + pow(a*D,2)
        
        constant_g2 = F
        constant_g = -2*F*beta - G
        constant_1 = F*pow(beta,2) + G*beta + H
        
        discriminant = pow(constant_g,2) - 4*constant_g2*constant_1
        
        growth[i,0] = -(constant_g + math.sqrt(discriminant))/(2*constant_g2)
        growth[i,1] = -(constant_g - math.sqrt(discriminant))/(2*constant_g2)
        
    return growth