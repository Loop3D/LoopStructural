"""
TODO This file needs to be re-written and use numpy properly
########
"""
import numpy as np
@np.vectorize
def distance(p1,p2):
    return abs(p1-p2)
@np.vectorize
def covar(p1,p2):
    return ((p1-p2)**2)
@np.vectorize
def inside(a,l,u):
    if a > l and a<u:
        b= 0
    else:
        b = 1
    return b
class s_variogram():
    def __init__(self,xdata,ydata):
        self.xdata = xdata
        self.ydata = ydata
    def setup(self):
        x1 = np.array([self.xdata,]*len(self.xdata))
        y1 = np.array([self.ydata,]*len(self.ydata))
        #x2 = np.tile(xdata,len(xdata))
        y2 = y1.transpose()
        x2 = x1.transpose()
        self.distance_m = np.linalg.norm(x1-x2,axis=1)
        self.covariance_m = (y1-y2)**2#covar(y1,y2)#x1 = np.array(np.)    
    def calc_semivariogram(self, step, nlags, tol):
        self.lags = np.arange(step/2.,nlags*step,step)
        self.variance,self.npairs = self.semivariogram(self.lags*1.1,step,self.distance_m,self.covariance_m)
        return self.lags, self.variance, self.npairs
    def semivariogram(self,lags,tol,distance,covariance):
        variance = np.zeros(len(lags))
        npairs = np.zeros(len(lags))
        self.min_ = np.zeros(len(lags))
        self.max_ = np.zeros(len(lags))
        for i in range(len(lags)):
            ma = np.ma.array(data = covariance, mask =inside(distance,
                lags[i]-tol/2.,lags[i]+tol/2.))
            if len(ma[~ma.mask])>0:
                variance[i] = np.mean(ma[~ma.mask])     
                self.min_[i] = np.percentile(ma[~ma.mask],25)
                self.max_[i] = np.percentile(ma[~ma.mask],75)
            else:
                variance[i] = np.nan
                self.min_[i] = np.nan
                self.max_[i] = np.nan
            npairs[i] = ma.count()
        return variance,npairs
    def find_wavelengths(self,step=0,nlags=0):

        if step==0:
            minxx = np.min(self.xdata)
            maxx = np.max(self.xdata)
            sorted_x = np.sort(self.xdata)
            av_dist = 0.
            c = 0
            dist = 0
            for i in range(len(self.xdata)):
                if i ==0:
                    dist += sorted_x[i+1]-sorted_x[i]
                    c+=1
                if i == len(self.xdata)-1:
                    dist += sorted_x[i]-sorted_x[i-1]
                    c+=1
                else:
                    dist += sorted_x[i]-sorted_x[i-1]
                    dist += sorted_x[i+1]-sorted_x[i]
                    c+=2

            step = dist /c #abs((float((maxx - minxx)) / float(len(self.xdata))))
            step*=1.2
        if nlags == 0:
            distance = np.abs(np.min(self.xdata)-np.max(self.xdata))
            nlags = (distance / step)
        self.h, self.var, self.npairs = self.calc_semivariogram(step,nlags,step)

        self.px, self.py = self.find_peaks_and_troughs(self.h,self.var)
        
        self.averagex = []
        self.averagey = []
        for i in range(len(self.px)-1):
            self.averagex.append((self.px[i]+self.px[i+1])/2.)
            self.averagey.append((self.py[i]+self.py[i+1])/2.)
            i+=1 #iterate twice
        #find the extrema of the average curve
        self.px2, self.py2 = self.find_peaks_and_troughs(self.averagex,self.averagey)
        self.wl1 = 0.
        wl1py = 0.
        for i in range(len(self.px)):
            if i > 0 and i < len(self.px)-1:
                if self.py[i] > 10:
                    
                    if self.py[i-1] < self.py[i]*.7:
                        if self.py[i+1] < self.py[i]*.7:
                            self.wl1 = self.px[i]
                            if self.wl1 > 0.:
                                wl1py = self.py[i]
                                break
        self.wl2 = 0.
        for i in range(len(self.px2)):
            if i > 0 and i < len(self.px2)-1:
                if self.py2[i-1] < self.py2[i]*.90:
                    if self.py2[i+1] < self.py2[i]*.90:
                        self.wl2 = self.px2[i]
                        if self.wl2 > 0. and self.wl2 > self.wl1*2 and wl1py < self.py2[i]:
                            
                            break
        if self.wl1 == 0.0 and self.wl2 == 0.0:
            return 0.0, 2*(maxx-minxx)
        return self.wl1*2., self.wl2*2.
    def find_peaks_and_troughs(self,x,y):
        if len(x) != len(y):
            return False
        pairsx = []
        pairsy = []
        for i in range(0,len(x)):
            if i < 1:
                pairsx.append(x[i])
                pairsy.append(y[i])

                continue
            if i > len(x)-2:
                pairsx.append(x[i])
                pairsy.append(y[i])
                continue
            left_grad = (y[i-1]-y[i]) / (x[i-1]-x[i])
            right_grad = (y[i]-y[i+1]) / (x[i]-x[i+1])
            if np.sign(left_grad) != np.sign(right_grad):
                pairsx.append(x[i])
                pairsy.append(y[i])
        return pairsx,pairsy
