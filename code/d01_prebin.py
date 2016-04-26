import numpy as np
#import statsmodels.api as sm
#import scipy.stats as stats
import pandas as pd
import math as math
import matplotlib.pyplot as plt
from scipy import stats
%matplotlib inline
from IPython.display import SVG, HTML

_prebin_params_doc = """
    Parameters
    ----------
    y: array-like
        The dependent variable, dim = n*1
    x: array-like
        The independnet variable, dim = n*p
    method: ['quantile', 'bucket']
        The default value is 'quantile'
    binnum: int
        binnum = 20 (default)
    groupnum: int
        groupnum = 4 (default)
    """
_prebin_Result_docs = """
    role: array
        Variables'roles: 1 is selected & 0 is out
    """ 
# Interactive Grouping
class prebin(object):
    __doc__="""
    The Prebin Process
    %(Params_doc)s
    %(Result_doc)s
    Notes
    ----
    """%{"Params_doc":_prebin_params_doc,
         "Result_doc":_prebin_Result_docs}    
    def __init__(self,data,xname,yname,**kwargs):
        self.data  = data
        self.yname = yname
        self.xname = xname
        self.xvalue = data[xname]
        self.yvalue = data[yname]
        self.binnum = 20
        if 'binnum' in kwargs.keys():
            self.binnum = kwargs['binnum']
        self.method = 'quantile'
        if 'method' in kwargs.keys():
            self.method = kwargs['method']
    def binning(self):
        
        # **** Split raw_data into nandata and nonandata ****
        if "missing" in list(self.xvalue):
            nan_tag  = 1
            nanlist = self.xvalue[self.xvalue == "missing"].index
            nandata = self.data[self.xvalue == "missing"]
            nonandata = self.data.drop(nanlist)
            nonandata[self.xname] = map(float, nonandata[self.xname])

            nandata[self.xname+"_bin"] = "bin_"+str(100)
            nandata["cut_point"]       = "missing"
        else:
            nandata = pd.DataFrame()
            nonandata = self.data

        nonandata[self.xname+'_bin'] = nonandata[self.xname]
        nonandata['cut_point']  = nonandata[self.yname]
        
        if (len(set(nonandata[self.xname][:1000])) >= self.binnum):
            if (self.method == 'quantile'):
                bincut = stats.mstats.mquantiles(nonandata[self.xname], prob=np.arange(0,1,1./self.binnum))
                #print(bincut)
            elif (self.method == 'bucket'):
                minvalue = min(nonandata[self.xname])
                maxvalue = max(nonandata[self.xname])
                bincut = np.arange(minvalue, maxvalue, 1.*(maxvalue-minvalue)/self.binnum)
            else:
                print("Error Message: Wrong Prebin mehtod ! Use: 'quantile' or 'bucket' ...")
                raise SystemExit
            bincut = list(set(bincut))   # remove duplicate
        else:
            bincut = list(set(nonandata[self.xname]))
        bincut.sort()
        for i in np.arange(len(bincut)+1):
            if (i == 1):
                nonandata[self.xname+'_bin'][nonandata[self.xname] <= bincut[i]] = 'bin_'+str(100+i)
                nonandata['cut_point'][nonandata[self.xname] <= bincut[i]] = bincut[i]
            elif (i>1 and i<len(bincut)):
                nonandata[self.xname+'_bin'][(nonandata[self.xname] > bincut[i-1])&(nonandata[self.xname]<=bincut[i])] = 'bin_'+str(100+i) 
                nonandata['cut_point'][(nonandata[self.xname] > bincut[i-1])&(nonandata[self.xname]<=bincut[i])] = bincut[i]
            elif (i == len(bincut)):
                nonandata[self.xname+'_bin'][nonandata[self.xname] > bincut[i-1]] = 'bin_'+str(100+i)
                nonandata["cut_point"][nonandata[self.xname] > bincut[i-1]] = max(nonandata[self.xname])

        
        newdata = [nandata, nonandata]
      
        return newdata

class grouping(object):
    def __init__(self, data, yname, xname, **kwargs):
        self.yname = yname
        self.xname = xname
        self.nandata = data[0]
        self.data  = data[1]
        self.mingroupsize = 0.05
        if 'mingroupsize' in kwargs.keys():
            self.mingroupsize = kwargs['mingroupsize']
    def freq(self):
        Dict     = {}
        freqs    = list()
        sumy     = {"y0":0., "y1":0.}
        var_name = list()
        ylist    = list()
        cut      = list()
        _tmpcut_ = list()
        tmpcut   = list()
        _tmpx_   = list(set(self.data[self.xname]))
        _tmpy_   = list(set(self.data[self.yname]))
        _tmpx_.sort()
        _tmpy_.sort()
        for name in _tmpx_:
            for y in _tmpy_:
                var_name.append(name)
                ylist.append(y)
                freqs.append(len(self.data[(self.data[self.xname] == name) & (self.data[self.yname] == y)]))
                cut.append(max(self.data.cut_point[self.data[self.xname] == name]))
            _tmpcut_.append(max(self.data.cut_point[self.data[self.xname] == name]))
            #print(_tmpcut_)
        for k,v in enumerate(_tmpcut_):
            if k > 0 :
                tmpcut.append((_tmpcut_[k-1]+_tmpcut_[k])/2.)
        freq_table = pd.DataFrame(np.vstack((ylist,var_name,cut,freqs)).T, columns = [self.yname,self.xname,'Cut','Count'])
        freq_table[[self.yname,'Cut','Count']] = freq_table[[self.yname,'Cut','Count']].astype(float)
        Dict = {"cutpoint":_tmpcut_, "binx":_tmpx_, "freqy":sumy, "freqtable":freq_table}
        return Dict
    
    def entropy(self, _tmpdata_, y, cut):
        Dict = {}
        toa  = 1.*sum(y.values())
        y_0  = y["y0"]/toa
        y_1  = y["y1"]/toa 
        x_1  = 1.*sum(_tmpdata_.Count[_tmpdata_.Cut <= cut])
        x_2  = toa - x_1
        x11  = 1.*sum(_tmpdata_.Count[(_tmpdata_.Cut <= cut) & (_tmpdata_[self.yname] == 1)])
        x10  = x_1-x11
        x21  = 1.*sum(_tmpdata_.Count[(_tmpdata_.Cut > cut) & (_tmpdata_[self.yname] == 1)])
        x20  = x_2 - x21
        entS = -y_0*np.log2(y_0) - y_1*np.log2(y_1)
        if x10 == 0:
            entx10 = 0
        else:
            entx10 = (x_1/toa)*(-(x10/x_1)*np.log2(x10/x_1))
            
        if x11 == 0:
            entx11 = 0
        else:
            entx11 = (x_1/toa)*(-(x11/x_1)*np.log2(x11/x_1))
            
        if x20 == 0:
            entx20 = 0
        else :
            entx20 = (x_2/toa)*(-(x20/x_2)*np.log2(x20/x_2))
        
        if x21 == 0:
            entx21 = 0
        else:
            entx21 = (x_2/toa)*(-(x21/x_2)*np.log2(x21/x_2))  
            
        entx = entx10 + entx11 + entx20 + entx21
        Dict = {cut:(entS-entx)}

        return Dict 
    def calsplit(self, predata):
        tmpdict = {}
        initvalue = 0.
        cutvalue = {}
        cutpoint = list(set(predata.Cut))
        freqy = {}
        freqy["y0"] = sum(predata.Count[predata[self.yname] == 0])
        freqy["y1"] = sum(predata.Count[predata[self.yname] == 1])
        cutpoint.sort()
        for i, value in enumerate(cutpoint):
            if i < len(cutpoint)-1:
                _binentropy_ = self.entropy(predata, freqy, value)
                if initvalue <= _binentropy_.values()[0]:
                    cutvalue = _binentropy_

                    initvalue = _binentropy_.values()[0]
                tmpdict[i] = _binentropy_

        return cutvalue.keys()[0]
    def split(self):
        bestgroup = {}
        cutlist   = list()
        predata   = self.freq()
        freqtable = predata['freqtable']
        total = sum(freqtable.Count)
        minsize = round(total*self.mingroupsize, 0)
        #print(predata)
        tmpdict0 = self.calsplit(freqtable)
        left0 = sum(freqtable.Count[freqtable.Cut <= tmpdict0])
        right0 = total - left0

        if (left0 < minsize) or (right0 < minsize):
            print "ErrorMessage: mingroupsize is not satisfied!"
            raise SystemExit
        else :
            leftdata0  = freqtable[freqtable.Cut <= tmpdict0]
            leftdata0.index = np.arange(len(leftdata0))
            rightdata0 = freqtable[freqtable.Cut > tmpdict0]
            rightdata0.index = np.arange(len(rightdata0))
            cutlist.append(tmpdict0)

            if left0 >= 2*minsize :
                #print("^^^^^^")
                tmpdict1 = self.calsplit(leftdata0)
                left1  = sum(leftdata0.Count[leftdata0.Cut <= tmpdict1])
                right1 = left0 - left1
                if (left1 >= minsize) & (right1 >= minsize):
                    cutlist.insert(0, tmpdict1)
            if right0 >= 2*minsize :

                tmpdict2 = self.calsplit(rightdata0)
                left2  = sum(rightdata0.Count[rightdata0.Cut <= tmpdict2])
                right2 = right0 - left2

                if (left2 >= minsize) & (right2 >= minsize):
                    cutlist.append(tmpdict2)
        cutlist.sort()
        return [cutlist, freqtable]
    
    def calwoe(self):
        splitinfo = self.split()
        cutlist   = splitinfo[0]
        groupdata = splitinfo[1] 
        totalevent = 1.*sum(groupdata.Count[groupdata[self.yname] == 1])
        totalnoevent = 1.*sum(groupdata.Count[groupdata[self.yname] == 0])

        ############  Calculate WOE ###############
        woetable = pd.DataFrame()
        groupname = list()
        groupvalue = list()
        Count  = list()
        Eventcnt = list()
        Noeventcnt = list()
        woevalues = list()
        eventcnt = list()
        noeventcnt = list()
        total = list()
        tmpiv    = 0. 

        for i, value in enumerate(cutlist):
            groupname.append(i+1)
            Count = sum(groupdata.Count[groupdata.Cut <= value])
            if i == 0 :
                groupvalue.append(" <= "+ str(value))
                _tmpevent   = sum(groupdata.Count[(groupdata.Cut <= value) & (groupdata[self.yname] == 1)])
                _tmpnoevent = sum(groupdata.Count[(groupdata.Cut <= value) & (groupdata[self.yname] == 0)])
            else :
                groupvalue.append(str(cutlist[i-1])+ " ~ "+str(value))
                _tmpevent   = sum(groupdata.Count[(groupdata.Cut <= value) & (groupdata[self.yname] == 1) & (groupdata.Cut > cutlist[i-1])])
                _tmpnoevent = sum(groupdata.Count[(groupdata.Cut <= value) & (groupdata[self.yname] == 0) & (groupdata.Cut > cutlist[i-1])])

            _tmpwoevalue = np.log((_tmpevent/totalevent)/(_tmpnoevent/totalnoevent))
            tmpiv = tmpiv + (_tmpevent/totalevent - _tmpnoevent/totalnoevent)*_tmpwoevalue
            Eventcnt.append(_tmpevent)
            Noeventcnt.append(_tmpnoevent)
            woevalues.append(_tmpwoevalue)
        
        #***** i = len(cutlist) ****
        groupname.append(i+2)
        Count = (totalevent + totalnoevent)
        groupvalue.append(str(value) + " ~ "+ str(max(groupdata.Cut)))
        _tmpevent   = sum(groupdata.Count[(groupdata.Cut > value) & (groupdata[self.yname] == 1)])
        _tmpnoevent = sum(groupdata.Count[(groupdata.Cut > value) & (groupdata[self.yname] == 0)])
        _tmpwoevalue = np.log((_tmpevent/totalevent)/(_tmpnoevent/totalnoevent))
        tmpiv = tmpiv + (_tmpevent/totalevent - _tmpnoevent/totalnoevent)*_tmpwoevalue
        Eventcnt.append(_tmpevent)
        Noeventcnt.append(_tmpnoevent)
        woevalues.append(_tmpwoevalue)
        
        
        if len(self.nandata) > 0:
            groupname.append(i+3)
            groupvalue.append("missing")
            Count = len(self.nandata)
            _tmpevent = sum(self.nandata[self.yname])
            _tmpnoevent = Count - _tmpevent
            _tmpwoevalue = np.log((_tmpevent/totalevent)/(_tmpnoevent/totalnoevent))
            tmpiv = tmpiv + (_tmpevent/totalevent - _tmpnoevent/totalnoevent)*_tmpwoevalue
            Eventcnt.append(_tmpevent)
            Noeventcnt.append(_tmpnoevent)
            woevalues.append(_tmpwoevalue)
        
        woetable["groups"]  = groupname
        woetable["interval"] = groupvalue
        woetable["Event_cnt"]  = Eventcnt
        woetable["NoEvent_cnt"]  = Noeventcnt
        woetable["WOE"]  = woevalues
        
        # **** Print WOE Chart ****
        plt.subplot(1,2,1)
        axi = np.arange(len(groupname))+1.1
        plt.xticks(axi, groupname)
        plt.plot(axi, woevalues, 'bo-')
        plt.grid(True)
        plt.title('WOE')
        plt.ylabel('WOE Value')
        plt.xlabel('Group')
        plt.subplot(1,2,2)
        plt.xticks(axi, groupname)
        bar1 = plt.bar(axi, Noeventcnt, width=0.5,align='center', color='b')
        bar2 = plt.bar(axi, Eventcnt, width=0.5,align='center', color='m')
        plt.legend((bar1[0], bar2[0]), ('noevent','event'))
        plt.grid(True)
        plt.title('Group Distribution')
        plt.ylabel('Count')
        plt.xlabel('Group')
        plt.show()
        
        return [woetable, tmpiv]
      