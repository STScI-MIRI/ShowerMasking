#! /usr/bin/python

# Author: Mike Engesser
# Last Update: 30-Aug-2023

# Native Imports
import argparse
import sys
import warnings

# Third Party Imports
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.modeling.models import Ellipse2D
from astropy.stats import sigma_clipped_stats

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.patches import Ellipse
from matplotlib.widgets import Slider

import numpy as np
from PIL import Image, ImageTk
import tkinter as tk


class InteractiveShowerMasking:
    """ 
    This class creates a GUI object for visualizing MIRI data
    and interactively masking cosmic ray showers using ellipses.
    Some showers are not easily detected by the JWST Calibration 
    Pipeline's JUMP step, persist for long periods of time and thus
    leave residual flux in the "_rate.fits" files. 
    Showers are masked from their first group until the end of the
    integration. 
    """
    
    def __init__(self, fname):
        """
        Initializes the class, creating the GUI.
        
        Parameters:
            
            root: Tkinter.Tk() object
                Tkinter GUI
            
            fname: str
                The name of the input file
        """
        
        self.fname = fname # File name from class instantiation
        
        self.rate = False # Flag for whether the input is "rate" or "rate_ints" data

        self.get_data(fname) # Read in the data 
        
        # Create the Tkinter app
        self.root = tk.Tk()
        self.root.title("Interactive Shower Masking")
        
        # Data and Ellipse parameters
        self.int = 0
        self.group = 0
        self.x = 500,
        self.y = 500,
        self.a = 50
        self.b = 50
        self.theta = 0

        self.coordinates_names = ['x','y','a','b','theta','integration','group']
        self.coordinates_list = []
        
        # Set array to first group of data cube
        self.array = np.array([])#self.data[self.int,self.group]
        
        # Stretch values for pyplot figure
        self.vmin = 0
        self.vmax = 1
        
        # Define scale factor for stretch
        self.vscale = 3
        
        # Define colormap for figure
        self.cmap = 'BuPu'
        
        # Updates the figure with the current slice of 
        # the data cube
        self.update_arr()
        
        # Set DQ transparency
        self.alpha = 0

        # Create a pyplot figure and axes object
        self.fig, self.ax = plt.subplots(figsize=(10,10))
        self.im = self.ax.imshow(self.array, vmin=self.vmin, vmax=self.vmax, 
                       cmap=self.cmap,origin='lower',interpolation='nearest')
        
        self.im_dq = self.ax.imshow(self.dq_arr, vmin=0, vmax=4, 
                       cmap='Set3',origin='lower',interpolation='nearest', alpha=self.alpha)
        
        # Create a Tk canvas object to display the pyplot figure
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        
        # Create GUI elements
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.grid(row=0, column=0,rowspan=13, columnspan=13)
        
        # Make the figure zoomable
        self.toolbar = NavigationToolbar2Tk( self.canvas, self.root,pack_toolbar=False )
        self.toolbar.update()
        self.toolbar.grid(row=14, column=0)
        
        # Slider element for angle
        self.theta_slider = tk.Scale(self.root, from_=0, to=360, orient="horizontal", command=self.update_theta)
        self.theta_slider.set(10)
        self.theta_slider.grid(row=1, column=18, columnspan=4)
        self.theta_label = tk.Label(self.root, text="Theta:")
        self.theta_label.grid(row=1, column=17)
        
        # Slider element for Semi-Major Axis
        self.a_label = tk.Label(self.root, text="Semi-Major Axis:")
        self.a_label.grid(row=2, column=17)
        self.a_slider = tk.Scale(self.root, from_=1, to=500, orient="horizontal", command=self.update_a)
        self.a_slider.set(50)
        self.a_slider.grid(row=2, column=18, columnspan=4)

        # Slider element for Semi-Minor Axis
        self.b_label = tk.Label(self.root, text="Semi-Minor Axis:")
        self.b_label.grid(row=3, column=17)
        self.b_slider = tk.Scale(self.root, from_=1, to=500, orient="horizontal", command=self.update_b)
        self.b_slider.set(50)
        self.b_slider.grid(row=3, column=18, columnspan=4)
        
        # Button elements for changing the integration 
        self.int_label = tk.Label(self.root, text='Integration: {}'.format(self.int))
        self.int_label.grid(row=4, column=17)
        self.increment_button = tk.Button(self.root, text="+", command=self.increment_int)
        self.increment_button.grid(row=4, column=18)
        self.decrement_button = tk.Button(self.root, text="-", command=self.decrement_int)
        self.decrement_button.grid(row=4, column=19)
        
        # Button elements for changing the group 
        self.group_label = tk.Label(self.root, text='Group: {}'.format(self.int))
        self.group_label.grid(row=5, column=17)
        self.increment_group_button = tk.Button(self.root, text="+", command=self.increment_group)
        self.increment_group_button.grid(row=5, column=18)
        self.decrement_group_button = tk.Button(self.root, text="-", command=self.decrement_group)
        self.decrement_group_button.grid(row=5, column=19)
        
        # Slider element for stretch
        self.vscale_label = tk.Label(self.root, text="Image Scale Factor:")
        self.vscale_label.grid(row=6, column=17)
        self.vscale_slider = tk.Scale(self.root, from_=0.1, to=5, orient="horizontal", command=self.update_vscale)
        self.vscale_slider.set(3)
        self.vscale_slider.grid(row=6, column=18, columnspan=4)
        
        # Slider element for DQ plane overlay
        self.alpha_label = tk.Label(self.root, text="DQ Opacity:")
        self.alpha_label.grid(row=7, column=17)
        self.alpha_slider = tk.Scale(self.root, from_=0, to=100, orient="horizontal", command=self.update_alpha)
        self.alpha_slider.set(0)
        self.alpha_slider.grid(row=7, column=18, columnspan=4)
        
        # Button for recording a cosmic ray shower coordinate and ellipse parameters
        self.record_button = tk.Button(self.root, text="Mark Shower", command=self.record_settings)
        self.record_button.grid(row=8, column=17, columnspan=1)
        
        # Button for updating the DQ array and writing to new file
        self.save_button = tk.Button(self.root, text="Save", command=self.save_new_fits)
        self.save_button.grid(row=8, column=18, columnspan=1)
        
        # Button to quit application
        self.quit_button = tk.Button(self.root, text="Quit", command=self.root.destroy)
        self.quit_button.grid(row=8, column=19, columnspan=1)

        # Connect the mouse click event to the figure
        self.fig.canvas.mpl_connect('button_press_event', self.on_plot_click)
        
        # Update the plot
        self.plot_image()


    def get_data(self,fname):
        """
        Open an HDUList and get its data and DQ array.
        
        Parameters:
            fname: str
                name of the .fits file
        """
        
        self.hdu = fits.open(fname)
        self.data = self.hdu[1].data
        self.dq = self.hdu[3].data
        
        self.head = self.hdu[0].header
        
        # Check if data is output from RAMP_FIT step
        if 'S_RAMP' in self.head.keys():
            self.rate = True
            
        # Convert from groups to CDS
        if self.rate == False:
            
            ndims = len(self.data.shape)
            
            if ndims == 4:
                int_len = self.data.shape[0]
                ramp_len = self.data.shape[1]
                new_data = np.zeros_like(self.data)
                for k in range(int_len):
                    for g in range(ramp_len):
                        if g == 0:
                            new_data[k,g] = self.data[k,g]
                        else:
                            new_data[k,g] = self.data[k,g] - self.data[k,g-1]
            if ndims == 3:
                ramp_len = self.data.shape[0]
                new_data = np.zeros_like(self.data)
                for g in range(ramp_len):
                    if g == 0:
                        new_data[g] = self.data[g]
                    else:
                        new_data[g] = self.data[g] - self.data[g-1]

            self.data = new_data
        
        return
    
    def update_arr(self):
        """Update the data array and the stretch parameters."""
        
        ndims = len(self.data.shape)
        if ndims == 4:
            self.array = self.data[self.int, self.group]
            self.dq_arr = self.dq[self.int, self.group].astype(float)
        
        elif ndims == 3:
            if self.rate == True:
                if self.head['S_RAMP'] == 'COMPLETE':
                    self.array = self.data[self.int]
                    self.dq_arr = self.dq[self.int].astype(float)
            else:
                self.array = self.data[self.group]
                self.dq_arr = self.dq[self.group].astype(float)
        else:
            self.array = self.data
            self.dq_arr = self.dq.astype(float)
        
        self.dq_arr[self.dq_arr == 0] = np.nan
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mean, med, std = sigma_clipped_stats(self.array, sigma=3)
        
        self.vmin = med - self.vscale*std
        self.vmax = med + self.vscale*std
        

    def plot_image(self):
        """Update the figure to reflect any changes."""
        
        # Update the data array to the current slice
        self.update_arr()
        
        #clear the patches for redrawing
        [p.remove() for p in reversed(self.ax.patches)]
        
        # Show the data array
        self.im.set_data(self.array)
        self.im.set_clim(vmin=self.vmin,vmax=self.vmax)
        
         # Show the dq array
        self.im_dq.set_data(self.dq_arr)
        self.im_dq.set_alpha(self.alpha)
        
        # Draw an ellipse
        ellipse = Ellipse(xy=(self.x, self.y), width=self.a, height=self.b, angle=self.theta,
                          edgecolor='r', fc='None', lw=2)        
        self.ax.add_patch(ellipse)
        
        # Draw any saved ellipses
        self.preserve_ellipses()
        
        # Draw the new figure
        self.ax.axis('off')
        self.canvas.draw()
        
    def preserve_ellipses(self):
        """Add saved ellipses to figure each time it is updated."""
        
        # Loop through saved ellipse coordinates
        for coord in self.coordinates_list:
            
            x = coord[0]
            y = coord[1]
            width = coord[2]
            height = coord[3]
            theta = coord[4]
            integration = coord[5]
            group = coord[6]
            
            # We only want to show ellipses in the same integration and
            # group that was marked, as well as all subsequent groups in 
            # that integration.
            if integration == self.int and group <= self.group:
                ellipse = Ellipse(xy=(x, y), width=width, height=height, angle=theta,
                          edgecolor='lime', fc='None', lw=2)        
                self.ax.add_patch(ellipse)


    def update_theta(self, value):
        """Update the theta parameter and redraw the figure."""
        
        self.theta = int(value)
        self.plot_image()
        
    def update_a(self, value):
        """Update the 'a' parameter and redraw the figure."""

        self.a = int(value)
        self.plot_image()

    def update_b(self, value):
        """Update the 'b' parameter and redraw the figure."""
        
        self.b = int(value)
        self.plot_image()
        
    def update_vscale(self, value):
        """Update the 'b' parameter and redraw the figure."""
        
        self.vscale = int(value)
        self.plot_image()
        
    def on_plot_click(self, event):
        """
        Record the x and y coordinate of a mouse click and redraw the figure.
        """
        
        if event.inaxes == self.ax:
            self.x = int(event.xdata)
            self.y = int(event.ydata)
            self.plot_image()
            
    def record_settings(self):
        """"Save the current settings to a running list of shower coordinates."""
        
        if self.x is not None and self.y is not None:
            self.coordinates_list.append((self.x, self.y, self.a, self.b, self.theta,self.int, self.group))
            
    def increment_int(self):
        """Increment the integration number by 1 and redraw the figure."""
        
        try:
            self.int += 1
            self.plot_image()
            self.int_label.config(text=f"Integration: {self.int}")
        except ValueError:
            pass

    def decrement_int(self):
        """Decrement the integration number by 1 and redraw the figure."""

        try:
            self.int -= 1
            self.plot_image()
            self.int_label.config(text=f"Integration: {self.int}")
        except ValueError:
            pass
        
    def increment_group(self):
        """Increment the group number by 1 and redraw the figure."""

        try:
            self.group += 1
            self.plot_image()
            self.group_label.config(text=f"Group: {self.group}")
        except ValueError:
            pass

    def decrement_group(self):
        """Decrement the group number by 1 and redraw the figure."""
        
        try:
            self.group -= 1
            self.plot_image()
            self.group_label.config(text=f"Group: {self.group}")
        except ValueError:
            pass
        
    def update_alpha(self, value):
        """Update the 'alpha' parameter and redraw the figure."""
        
        self.alpha = float(value)/100
        self.plot_image()
        
        
    def update_dq(self):
        """Update the DQ array to mask each marked shower."""
        
        mask = np.zeros_like(self.data)

        # Loop through saved ellipse coordinates
        for coord in self.coordinates_list:
            x = coord[0]
            y = coord[1]
            width = coord[2]
            height = coord[3]
            theta = Angle(coord[4], 'deg') # create an Angle() object
            integration = coord[5]
            group = coord[6]

            # Create an ellipse where values internal to the boundary are set to 1
            e = Ellipse2D(amplitude=1, x_0=x, y_0=y, a=width, b=height,
                      theta=theta.radian)

            # Add the ellipse to the first group and each subsequent group
            # in a single integration.
            
            ndims = len(self.data.shape)
            if ndims == 4:
                for g in range(group,self.data.shape[1]):

                    y,x = np.mgrid[0:self.data.shape[-2],0:self.data.shape[-1]]
                    mask[integration,g] += e(x,y)

            elif ndims == 3:
                if self.rate == True:
                    if self.head['S_RAMP'] == 'COMPLETE':
                        y,x = np.mgrid[0:self.data.shape[-2],0:self.data.shape[-1]]
                        mask[integration] += e(x,y)
                        
                else:
                    for g in range(group,self.data.shape[0]):
                        y,x = np.mgrid[0:self.data.shape[-2],0:self.data.shape[-1]]
                        mask[g] += e(x,y)
                            
            else:
                y,x = np.mgrid[0:self.data.shape[-2],0:self.data.shape[-1]]
                mask += e(x,y)
            
        # Set mask locations to 4 in DQ array ("JUMP")
        self.dq[mask >= 1] = 4
        
    def save_new_fits(self):
        """Update the DQ array and save the data to a new .fits file."""
        
        self.update_dq()
        
        self.hdu[3].data = self.dq
        new_fname = self.fname.replace('.fits','_showers_masked.fits')
        self.hdu.writeto(new_fname,overwrite=True)
        self.hdu.close()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--fname')
    
    args = parser.parse_args()
        
    root = tk.Tk()
    fname = args.fname
    app = InteractiveShowerMasking(fname)