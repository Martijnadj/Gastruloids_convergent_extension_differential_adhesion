import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from PIL import Image

import imageio
from nd2reader import ND2Reader
from tkinter import simpledialog
import time
from os import listdir
import numpy as np
import copy

import tkinter
import ciha
from tkinter import filedialog


def open_dir():
    
    global loc
    global files
    global n_im
    global iter_load
    
    iter_load = True
    
    typ = filetype.get()
    
    loc = filedialog.askdirectory()
    
    if typ == "nd2":
        files = [f for f in listdir(loc) if f.endswith('.nd2')]
    if typ == "tiff":
        files = [f for f in listdir(loc) if f.endswith('.tif') or f.endswith('.TIF') or f.endswith('.tiff') or f.endswith('.TIFF')]
    elif typ == "png" or typ == "jpg":
        files = [f for f in listdir(loc) if f.endswith('.png') or f.endswith('.jpg')]
    elif typ == "any":
        files = [f for f in listdir(loc) if f.endswith('.png') or f.endswith('.jpg') or f.endswith('.nd2') or f.endswith('.tif')]
    else:
        files = []

    print("Number of Images to segment: ", len(files))

def load():
    
    global ims
    global loc
    global files
    
    i = int(float(iter_im.get()))
    inv = invert.get()
    
    if files[i].endswith(".nd2"):
        loc_new = loc + "/" + files[i]
        ims = ciha.ND2toImg(loc_new)
        print(ims.shape)
        if len(ims.shape)==3:
            xy_num = ims.shape[0]
        else:
            xy_num = 0
    elif files[i].endswith(".png") or files[i].endswith(".jpg"):
        loc_new = loc + "/" + files[i]
        ims = imageio.imread(loc_new)
        if len(ims.shape) > 2:
            ims = np.max(ims, 2)
        xy_num = 1
    elif files[i].endswith('.tif') or files[i].endswith('.TIF') or files[i].endswith('.tiff') or files[i].endswith('.TIFF'):
        loc_new = loc + "/" + files[i]
        im = Image.open(loc_new)
        ims = np.array(im)
        if len(ims.shape) > 2:
            ims = np.max(ims[:,:,0:2], 2)
        xy_num = 1
    else:
        print("I couldn't find the files ...")

    if inv == "True":
        ims = np.max(ims) - ims

    
    print("Number of xy to segment: ", xy_num)
    
    return(ims)

def loc_contours():
    
    global loc_cont
    
    loc_cont = filedialog.asksaveasfilename()
    print("Location for saving contours: ", loc_cont)


def Save_cont():
    
    global contours
    global loc_cont
    global files
    
    i = int(float(iter_im.get()))
    j = int(float(iter_xy.get()))
    
    dist3 = int(float(dist.get()))
    min3 = int(float(min2.get()))
    sigma2 = float(sigma.get())
    op2 = float(op.get())
    excl2 = excl.get()
    thresh2 = thresh.get()
    
    ciha.save_contours_hyperparams(d = loc_cont, thresh = thresh2, dist=dist3, min_size=min3, sigma=sigma2, opening=op2, exclude=excl2)
    
    ciha.save_contours(d = loc_cont, cont=contours, name = files[i], xy = j)
    
    print("Contours saved!")


def next_image():
    
    global iter_load
    global files
    
    iter_load = True
    
    i = int(float(iter_im.get()))
    i += 1
    iter_im.set(i)
    
    print("Current image for segmentation: ", i, "/", len(files)-1)

def next_xy():
    
    global ims
    global iter_load
    
    iter_load = False
    
    j = int(float(iter_xy.get()))
    j += 1
    iter_xy.set(j)
    
    if len(ims.shape)==3:
        xy_num = ims.shape[0]
    else:
        xy_num = 1
    
    print("Current xy for segmentation: ", j, "/", xy_num)

def loc_fourier():
    global loc_fourier
    
    loc_fourier = filedialog.asksaveasfilename()
    print("Location for saving fourier components: ", loc_fourier)


def Save_fourier():
    
    global ft
    global loc_fourier
    global files
    
    i = int(float(iter_im.get()))
    
    acc2 = float(acc.get())
    compnum2 = int(float(compnum.get()))
    
    j = int(float(iter_xy.get()))
    
    ciha.save_fourier_hyperparams(d = loc_fourier, accuracy=acc2, m=compnum2)
    
    ciha.save_fourier(d = loc_fourier, ft = ft, name = files[i], xy = j)
    
    print("Fourier components saved!")

def run():
    
    global ims
    global contours
    global iter_load
    
    j = int(float(iter_xy.get()))
    
    if iter_load == True:
        ims = load()
        if (len(ims.shape) == 3):
            ims_current = ims[j,:,:]
        else:
            ims_current = ims
        iter_load = False
    
    else:
        if (len(ims.shape) == 3):
            ims_current = ims[j,:,:]
        else:
            ims_current = ims
    
    thresh2 = thresh.get()
    
    if thresh2 == "None":
        thresh3 = None
    else:
        thresh3 = float(thresh2)

    dist3 = int(float(dist.get()))
    min3 = int(float(min2.get()))
    sigma2 = float(sigma.get())
    op2 = float(op.get())
    excl2 = excl.get()
    
    contours, gauss_map, global_thresh, elevation_map, localMax, D, markers, labels, masks, sizes_gastrluoids =\
    ciha.extract_outline(ims_current, thresh=thresh3, dist=dist3, min_size=min3, sigma=sigma2, opening=op2, exclude = excl2)


    fig = Figure(figsize=(10, 10), dpi=100)

    ax1 = fig.add_subplot(331)
    ax1.imshow(ims_current)
    ax1.set_title('Original image')

    t = round(global_thresh, 1)
    ax2 = fig.add_subplot(332)
    ax2.hist(gauss_map.flatten(), 100)
    ax2.axvline(x=global_thresh, color = "g", label = str(t))
    ax2.set_title('Intensity histogram')
    ax2.legend()
    
    ax3 = fig.add_subplot(333)
    ax3.imshow(gauss_map)
    ax3.set_title('Gaussian Filter')
    
    ax4 = fig.add_subplot(334)
    ax4.imshow(elevation_map)
    ax4.set_title('Elevation for Watershed')
    ind = np.where(localMax == True)
    for i in range(len(ind[0])):
        ax4.plot(ind[1][i], ind[0][i], "ro")
        ax4.annotate(i, [ind[1][i], ind[0][i]], size = 15, color = "white")

    ax5 = fig.add_subplot(335)
    ax5.imshow(D)
    ax5.set_title('Distance Measure')

    ax6 = fig.add_subplot(336)
    ax6.imshow(markers)
    ax6.set_title('Binary Image for Marking')

    ax7 = fig.add_subplot(337)
    ax7.imshow(labels)
    ax7.set_title('Segmented objects')
    gastr = np.unique(masks)
    gastr = np.delete(gastr, 0)
    for i in gastr:
        ax7.plot(np.mean(np.where(masks == i)[1]), np.mean(np.where(masks == i)[0]), "ro")
        ax7.annotate(int(i), [np.mean(np.where(masks == i)[1]), np.mean(np.where(masks == i)[0])], size = 15, color = "white")
    
    ax8 = fig.add_subplot(338)
    ax8.imshow(masks)
    ax8.set_title('Final Objects')
    
    ax9 = fig.add_subplot(339)
    for i in range(len(contours)):
        ax9.plot(-contours[i][:,1], -contours[i][:,0])
        ax9.axis("equal")
    ax9.set_title('Outlines')

    fig.subplots_adjust(bottom=0.38, right=0.65, left = 0.05, top=0.95, hspace=0.4)


    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().grid(row = 3, column = 1, columnspan = 40, rowspan = 40, sticky = "W")


def fourier_components():
    
    global contours
    global ft
    
    acc2 = float(acc.get())
    compnum2 = int(float(compnum.get()))
    
    dist, interpol = ciha.calculate_dist(contours, accuracy=acc2)
    ft = ciha.fft_shapes(dist, m=compnum2, divide_f0=True)
    
    xs = np.arange(0, 1, 1/len(contours[0]))
    
    fig = Figure(figsize=(8, 10), dpi=100)
    
    ax1 = fig.add_subplot(311)
    for i in range(len(contours)):
        ax1.plot(contours[i][:,0], contours[i][:,1])
        ax1.axis("equal")
    
    ax1.set_title('Outlines')

    ax2 = fig.add_subplot(312)
    for i in range(len(contours)):
        s, y = ciha.lineIntegral(contours[i])
        y_plot = y/max(y)
        ax2.plot(s/max(s), y_plot, label = "actual distance", linewidth = 3)
        ax2.plot(xs, interpol[i](xs)/max(y), "black" , label = "interpolated", linewidth = 1)
    
    ax2.set_title('Smoothing of Distance Curve')
    #ax2.legend()
    
    ax3 = fig.add_subplot(313)
    ft_plot = copy.copy(ft)
    ft_plot[0] = ft_plot[0]/len(contours[0])
    ax3.plot(ft_plot)
    #ax3.set_xlabel('Fourier Components', fontsize=18)
    ax3.set_title('Fourier Components')
    
    fig.subplots_adjust(bottom=0.38, right=0.4, left = 0, top=0.95, hspace=0.3)
    
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().grid(row = 3, column = 8, columnspan = 40, rowspan = 40, padx = (45, 0))



#Interface Enviroment
root = tkinter.Tk()
root.wm_title("Sectioning Gastruloids")
root.geometry('1300x800')

thresh = tkinter.StringVar()
thresh.set("None")

dist = tkinter.StringVar()
dist.set("70")

min2 = tkinter.StringVar()
min2.set("1000")

sigma = tkinter.StringVar()
sigma.set("5")

op = tkinter.StringVar()
op.set("5")

acc = tkinter.StringVar()
acc.set("1")

compnum = tkinter.StringVar()
compnum.set("10")

iter_im = tkinter.StringVar()
iter_im.set("0")

iter_xy = tkinter.StringVar()
iter_xy.set("0")

excl = tkinter.StringVar()
excl.set("None")

filetype = tkinter.StringVar()
filetype.set("any")

invert = tkinter.StringVar()
invert.set("False")


#Declare Text entries
tkinter.Label(root, text="Threshold Value   ", bg = "orange").grid(row = 1, column = 1, columnspan=1, padx = (40, 0))
entThresh = tkinter.Entry(root, textvariable = thresh, width = 10).grid(row = 2, column = 1, columnspan=1,padx = (20, 0))

tkinter.Label(root, text="Minimal Distance", bg = "orange").grid(row = 1, column = 2, columnspan=1)
entDist = tkinter.Entry(root, textvariable = dist, width = 10).grid(row = 2, column = 2, columnspan=1)

tkinter.Label(root, text="  Minimal Size  ", bg = "orange").grid(row = 1, column = 3, columnspan=1)
entSize = tkinter.Entry(root, textvariable = min2, width = 10).grid(row = 2, column = 3, columnspan=1)

tkinter.Label(root, text="Gaussian Smoothing", bg = "orange").grid(row = 1, column = 4, columnspan=1)
entSigma = tkinter.Entry(root, textvariable = sigma, width = 10).grid(row = 2, column = 4, columnspan=1)

tkinter.Label(root, text="Opening for Maxima", bg = "orange").grid(row = 1, column = 5, padx = (0, 40), sticky = "W",columnspan=1)
entOp = tkinter.Entry(root, textvariable = op, width = 10).grid(row = 2, column = 5,  columnspan=1, padx = (0,40))

tkinter.Label(root, text="Smoothing Accuracy", bg = "violet").grid(row = 1, column = 8, columnspan=1, padx = (40, 0))
entAcc = tkinter.Entry(root, textvariable = acc, width = 10).grid(row = 2, column = 8, columnspan=1, padx = (40, 0))

tkinter.Label(root, text="Fourier Components", bg = "violet").grid(row = 1, column = 9, columnspan=1, padx = 0)
entCN = tkinter.Entry(root, textvariable = compnum, width = 10).grid(row = 2, column = 9, columnspan=1, padx = 0)

tkinter.Label(root, text="Current Image", fg = "blue").grid(row = 4, column = 0, columnspan=1)
entCN = tkinter.Entry(root, textvariable = iter_im, width = 6).grid(row = 5, column = 0, columnspan=1)

tkinter.Label(root, text="Current xy", fg = "blue").grid(row = 7, column = 0, columnspan=1)
entCN = tkinter.Entry(root, textvariable = iter_xy, width = 6).grid(row = 8, column = 0, columnspan=1)

tkinter.Label(root, text="Exclude Gastruloids").grid(row = 13, column = 0, columnspan=1)
entCN = tkinter.Entry(root, textvariable = excl, width = 6).grid(row = 14, column = 0, columnspan=1)

tkinter.Label(root, text="Invert Image").grid(row = 15, column = 0, columnspan=1)
entCN = tkinter.Entry(root, textvariable = invert, width = 6).grid(row = 16, column = 0, columnspan=1)

tkinter.Label(root, text="Image Type").grid(row = 1, column = 0, columnspan=1)
entCN = tkinter.Entry(root, textvariable = filetype, width = 6).grid(row = 2, column = 0, columnspan=1)

#Declare buttons
tkinter.Button(text='File Open', command=open_dir, bg = "blue").grid(row = 0, column = 0, padx = 0, pady = 10, rowspan=1)
tkinter.Button(text='Next File', command=next_image, fg = "blue").grid(row = 3, column = 0, padx = 0, pady = 10, rowspan=1)
tkinter.Button(text='Next xy', command=next_xy, fg = "blue").grid(row = 6, column = 0, padx = 0, pady = 10, rowspan=1)
tkinter.Button(text='Run Segmentation', command=run, fg = "orange").grid(row = 0, column = 1, padx = 0, pady = 10)
tkinter.Button(text='File Save Outlines', command=loc_contours,  fg = "orange").grid(row = 9, column = 0, padx = 0, pady = 10)
tkinter.Button(text='Save Outlines', command=Save_cont, fg = "orange").grid(row = 10, column = 0, padx = 0, pady = 10)
tkinter.Button(text='Run Fourier Analysis', command=fourier_components, fg = "violet").grid(row = 0, column = 8, padx = 0, pady = 10)
tkinter.Button(text='File Save Fourier', command=loc_fourier, fg = "violet").grid(row = 11, column = 0, padx = 0, pady = 10)
tkinter.Button(text='Save Fourier Components', command=Save_fourier, fg = "violet").grid(row = 12, column = 0, padx = 0, pady = 10)

def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
# Fatal Python Error: PyEval_RestoreThread: NULL tstate

#button = tkinter.Button(master=root, text="Quit", command=_quit)
#button.pack(side=tkinter.BOTTOM)

tkinter.mainloop()
# If you put root.destroy() here, it will cause an error if the window is
# closed with the window manager.
