import numpy as np
from matplotlib import pyplot as plt
import math
from scipy import interpolate
from scipy import ndimage as ndi
import copy
from skimage.filters import sobel
from skimage.morphology import watershed
from scipy.ndimage import label
import ciha
from skimage.measure import regionprops
from skimage import measure
from skimage.filters import threshold_minimum
from skimage.filters import threshold_otsu, threshold_local, gaussian, wiener, threshold_yen, scharr
from skimage.morphology import binary_opening, binary_closing, binary_erosion, binary_dilation,  disk
from skimage.feature import peak_local_max
from nd2reader import ND2Reader


def extract_outline(im, thresh = None, sigma = 5, dist = 70, opening = 10, min_size = 1000, exclude = "None"):

    
    gauss_map = gaussian(im, sigma = sigma)

        
    if thresh is None:
        global_thresh =  threshold_otsu(gauss_map)
    else:
        global_thresh = thresh
    
    
    print("Threshhold: ", np.round(global_thresh, decimals=3))
    
    #selemd = disk(dilation)
    selemo = disk(opening)
    
    markers2 = gauss_map > global_thresh
    markers3 = binary_dilation(markers2, selemo)
    markers = binary_opening(markers3, selemo)
    
    elevation_map = scharr(gauss_map)
    
 
    #Calculate euclidean Distances and Find the local maximum based on these
    D = ndi.distance_transform_edt(markers)
    localMax = peak_local_max(D, indices=False, min_distance=dist, labels=markers)
    
    #Exclude the local Max that have not been properly selected
    a = np.where(localMax>0)

    for i in range(len(a[0])):
        for j in range(i+1, len(a[0])):
            d = np.linalg.norm(np.array([a[0][i], a[1][i]]) - np.array([a[0][j], a[1][j]]))
            if d < dist:
                localMax[a[0][j], a[1][j]] = 0
    
    # Add a marker for the background
    background_marker = np.where(elevation_map == np.min(elevation_map))
    background_marker = [background_marker[0][0], background_marker[1][0]]
    localMax[background_marker[0], background_marker[1]] = True

    print("Number of local Maxima Found: ", np.sum(localMax > 0))
    
    mark = ndi.label(localMax, structure=np.ones((3, 3)))[0]
    labels = watershed(elevation_map, mark) 
    
    #Exclude objects that are too small
    sizes_gastrluoids = np.zeros(np.max(labels)+1)
    for i in range(np.max(labels)+1):
        sizes_gastrluoids[i] = np.round(np.sum(labels == i))
        if np.sum(labels == i) < min_size:
            ind_exclude = np.where(labels == i)
            labels[ind_exclude[0], ind_exclude[1]] = 0
            
    print("The Gastruloids with the following sizes have been segmented: ", sizes_gastrluoids)
    
    # define the border image
    border = np.ones_like(im)
    border[0,:] = 0; border[:,0] = 0; border[-1,:] = 0; border[:,-1] = 0
    
    # add segments that are not touching the border to masks
    k = 0
    masks = np.zeros_like(im)
    for i in range(1,np.max(labels)+1):
        seg_i = labels == i
        if(np.sum(seg_i) == np.sum(seg_i * border)):
            k += 1
            masks = masks + k*seg_i
    
    if exclude != "None":
        b = exclude.split(",")
        for i in b:
            masks[masks == int(i)] = 0
    
    
    ind = np.unique(masks)
    ind = np.delete(ind, 0)
    contours = list()
    for i in ind:
        aux_mask = masks == i
        contours.append(measure.find_contours(aux_mask, 0.5)[0])
    
    print("Number of segmented objects: ", len(contours))

    return(contours, gauss_map, global_thresh, elevation_map, localMax, D, markers, labels, masks, sizes_gastrluoids)
    

    
#Import ND2 images
def ND2toImg(loc):
    
    nd2 = ND2Reader(loc)
    
    if len(nd2.axes) == 6:
        ims = np.zeros([nd2.sizes["v"], nd2.sizes["z"], nd2.sizes["x"], nd2.sizes["y"]])
        for v in range(nd2.sizes["v"]):
            for z in range(nd2.sizes["z"]):
                ims[v,z,:,:] = nd2.get_frame_2D(c = 2,t=0,z=z, v =v)
        ims_maxproj = np.max(ims, 1)
    else:
        ims = np.zeros([nd2.sizes["z"], nd2.sizes["x"], nd2.sizes["y"]])
        for z in range(nd2.sizes["z"]):
            ims[z,:,:] = nd2.get_frame_2D(c = 1,t=0,z=z) # c=2 is dapi, no t dimension
        ims_maxproj = np.max(ims, 0)
    
    return(ims_maxproj)


#Saving File
def save_contours(d, cont, name, xy):   
    shapes = open(d, "a")
    for i in range(len(cont)):
        n = np.array(np.repeat(name + "--xy-"+ str(xy) + "--Cont-" + str(i), len(cont[i])), ndmin=2)
        sh = np.concatenate((n.T, cont[i]), axis = 1)
        np.savetxt(shapes, X=sh, fmt='%s', delimiter="\t")
        
    shapes.close()

def save_contours_hyperparams(d, thresh, dist, min_size, sigma, opening, exclude):
    shapes = open(d, "a")
    
    hyperparams = "#" + "Threshold = " + str(thresh) + ", Minimal Distance = " + str(dist) + \
    ", Minimum Size = " + str(min_size) + ", Gaussian filter = " + str(sigma) + \
    ", Opening = " + str(opening) + ", Exclude = " + str(exclude) + "\n"
    
    shapes.write(hyperparams)
    shapes.close()

#Saving Components
def save_fourier(d, ft, name, xy):
    f = open(d, "a")
    
    for i in range(ft.shape[1]):
        n = np.array(np.repeat(name + "--xy-"+ str(xy) + "--Cont-" + str(i), ft.shape[0]), ndmin=2)
        a = np.array(ft[:,i], ndmin=2)
        sh = np.concatenate((n.T, a.T), axis = 1)
        np.savetxt(f, X=sh, fmt='%s', delimiter="\t")
    
    f.close()

def save_fourier_hyperparams(d, accuracy, m):
    shapes = open(d, "a")
    
    hyperparams = "#" + "Accuracy = " + str(accuracy) + ", Number Fourier Components = " + \
    str(m) + "\n"
    
    shapes.write(hyperparams)
    shapes.close()

#Calculate the line Integral
def lineIntegral(xy):
    s_dist = np.zeros((len(xy)))
    for k in np.arange(1, len(xy)-1):
        s_dist[k] = np.linalg.norm(xy[k,:] - xy[k+1,:])
    
    s_dist[len(xy)-1] = np.linalg.norm(xy[0,:] - xy[len(xy)-1, :])
    s = np.cumsum(s_dist)
    
    s_long = np.concatenate((s, max(s) + s_dist[len(s_dist)-1]), axis = None)
    s_long = s_long/max(s_long)
    
    d = np.linalg.norm(xy - np.mean(xy, axis = 0), axis = 1)
    y = np.concatenate((d, d[0]), axis = None)
    
    ind = np.unique(s_long, return_index=True)[1]
    
    s = s_long[ind]
    y_new = y[ind]
    
    return(s, y_new)


#Functions to create test cases
def u(theta, a):
    return(a[0] + a[1]*np.sin(theta) + a[2]*np.cos(2*theta) + a[3]*np.cos(3*theta)) + a[4]*np.sin(4*theta)

def test_cases(a):
    
    theta = np.arange(0, 2*np.pi, 2*np.pi/1000)
    
    us = list()
    for i in range(a.shape[0]):
        us.append(u(theta, a[i,:]))
    
    test_cases = list()
    fig, ax = plt.subplots(figsize = (10,10))
    for i in range(a.shape[0]):
        test_gastruloid = np.array((us[i]*np.cos(theta), us[i]*np.sin(theta)))
        test_cases.append(test_gastruloid.T)
        plt.plot(test_gastruloid[0,:], test_gastruloid[1,:],"o", label = i)

    plt.legend()
    plt.axis("equal")
    plt.show()
    
    return(us, test_cases)


def calculate_dist(xy, n = 1000, accuracy = 100):

    all_dist = np.zeros((n, len(xy)))
    inter = list()

    for i in range(len(xy)):
        s, y = lineIntegral(xy[i])
        cc = interpolate.UnivariateSpline(x = s, y = y, k = 3, s = len(y)/accuracy) #len(y)/5
        inter.append(cc)
        s_equi = np.arange(0, 1, 1/n)
        all_dist[:,i] = cc(s_equi)
    
    return(all_dist, inter)


#Calculate the Fourier Transforms
def fft_shapes(all_dist, m = 5, divide_f0 = True):
    
    fft_distances = np.zeros((m, all_dist.shape[1]))
    for i in range(all_dist.shape[1]):
        if divide_f0 == True:
            fft_distances[:,i] = abs(np.fft.fft(all_dist[:,i])[0:m])/(len(all_dist[:,i]))
            fft_distances[1:m,i] = fft_distances[1:m,i]/fft_distances[0,i]
            
        else:
            fft_distances[:,i] = abs(np.fft.fft(all_dist[:,i])[1:m])/(len(all_dist[:,i]))
                
    return(fft_distances)


#Visualize the fourier transforms of the clusters
#Plot the different fourier transforms
def plot_all_ffts(fft, clusters):

    fig, ax = plt.subplots(figsize=(10,10))
    col = ["y", "blue", "orange", "green", "red", "violet", "brown", "pink", "grey"]
    cluster_num = max(clusters)
    
    for i in np.arange(1,cluster_num+1):
    
        clust = np.where(clusters == i)[0]
    
        for k in clust:
            
            plt.plot(fft[:,k], col[i-1]) 
        
        plt.show()

#Plot the gastruloids of a fourier component
def plot_gastruloid_fft(ft, cont, component, n, top = True):
    
    if top == True:
        k = np.flip(np.argsort(ft[component,:]))
    else:
        k = np.argsort(ft[component,:])

    for i in k[0:n]:
        print(ft[component,i])
        plt.plot(cont[i][:,0], cont[i][:,1])
        plt.axis("equal")
        plt.show()
        

    
