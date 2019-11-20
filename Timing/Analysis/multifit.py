import numpY as np
import matplotlib.pYplot as plt
from tqdm import tqdm

class dataset:
    def __init__(self,X,Y,Y_err,errortype='none'):
        if len(X) != len(Y):
            print('Invalid dataset: X has size %d and Y %d' %(len(X),len(Y)))
            raise TypeError
        else:
            self.X = np.array(X)
            self.Y = np.array(Y)
        self.use_errors = True
        if len(Y_err) == 0:
            if errortype == 'Poisson':
                self.Y_err = np.sqrt(self.Y)
            else:
                self.Y_err = []
                self.use_errors = False
        elif len(Y) != len(Y_err):
            print('Invalid dataset: Y has size %d and Y_err %d' %(len(Y),len(Y_err)))
            raise TypeError
        else:
            self.Y_err = np.array(Y_err)
            
    
    def subdataset(self,x_min,x_max):
        i_min, i_max = -1,-1
        for i,x in enumerate(self.X):
            if i_min == -1 and x > x_min:
                i_min = i
            if i_max == -1 and x > x_max:
                i_max = i
                break
        if i_min > i_max:
            print('Invalid subset choice: %f, %f' %(x_min,x_max))
            raise IndexError
            
        if self.use_errors:
            return dataset(self.X[i_min:i_max],self.Y[i_min:i_max],self.Y_err[i_min:i_max])
        else:
            return dataset(self.X[i_min:i_max],self.Y[i_min:i_max],[])
        
    def plot(self,show_errors=False,color='black',use_external_ax=False,figax=0):
        if use_external_ax:
            fig,ax = figax
        else:
            fig,ax = plt.subplots()
        if (not self.use_errors) or not show_errors:
            ax.scatter(self.X,self.Y)
        else:
            for i,x in enumerate(self.X):
                ax.plot([x,x],[self.Y[i]-self.Y_err[i],self.Y[i]+self.Y_err[i]],color=color)
        return fig,ax
    

color_list_datasets = ['black','blue','green']
color_list_functions = ['red','orange','brown']
l = len(color_list_datasets)

class multifit:
    def __init__(self,datasets,p_options,functions,p_lists):
        '''
        datasets: array of dataset (size D)
        p_options: array of np.arange of the parameters possible values (size P)
        functions: array of functions fitting each dataset (size D)
        p_lists: array of lists of indexs of the parameters used by each function (size (D))
        '''
        
        self.datasets = datasets
        self.p_options = p_options
        self.functions = functions
        self.p_lists = p_lists
        self.best_ps = []
        for p in p_options:
            self.best_ps.append(p[0])
        self.best_ps = np.array(self.best_ps)
        
        self.D = len(datasets)
        self.P = len(p_options)
        
    def plot(self,ids):
        '''
        ids is a list of indexs for selecting which dataset to plot
        '''
        first_time = True
        for j,i in enumerate(ids):
            # plot dataset
            if first_time:
                fig,ax = self.datasets[i].plot(color=color_list_datasets[j%l])
                first_time = False
            else:
                fig,ax = self.datasets[i].plot(color=color_list_datasets[j%l],use_external_ax=True,figax=(fig,ax))
            
            # plot fitting function
            minX = np.min(self.datasets[i].X)
            maxX = np.max(self.datasets[i].X)
            step = (maxX - minX)/(2*len(self.datasets[i].X))
            auxX = np.arange(minX,maxX,step)
            
            plt.plot() 
            
            
            
        