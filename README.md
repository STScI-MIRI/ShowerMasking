# Cosmic Ray Shower Masking Tool
Author: Mike Engesser  - mengesser@stsci.edu
## Instructions

This is a tool designed to simplify the process of masking cosmic ray showers in MIRI data that is not easily detected or removed by the JWST Calibration Pipeline. 

### Shower Artifacts and why they persist
Cosmic ray showers are a type of cosmic ray unique to MIRI (see https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-instrument-features-and-caveats) that are not easily detected by the standard JUMP step. An additional shower masking step attempts to detect the low level glow that persists for a number of groups and masks them using an ellipse in the DQ plane. 

Some showers are bright and persist for long durations of time. The JUMP step masks the shower for a default duration of 20s after the hit. However, many showers persist for much longer. This duration can be modified with the `time_masked_after_showers` argument. Complications arise when showers occur either 1) during idle, so no jump is detected, and persists into the integration or 2) with no bright primary hit, leaving the faint glow undetected.

### Manual Masking of Showers

The best way to circumvent this in the short term to is mask such showers by hand in Python. We can edit the Data Quality (DQ) plane of the data to include undetected showers. The process would go as follows: 

1. Save the output of the JUMP step using `save_results = True`.
2. Open the "_jump.fits" file in python/DS9/Jupyter Notebook.
3. Create an elliptical mask at the 4D coordinate of the shower.
4. Set all coordinates interior to that ellipse to a value of 4 (Jump) in the DQ array.
5. Do the same in each sebsequent group in the same integration. 
6. Repeat for all unmasked showers.
7. Save to a new ".fits" file and continue on with the JWST pipeline.

This tool is designed to simplify and streamline this process using an interactive GUI. The app takes a file name as input, and displays the data using `matplotlib`. 

#### Notes for use:
- There is a slider to adjust the image scale factor, as well as buttons to change the group and integration number.
- In order to make seeing the showers easier, each frame after the first is actually a difference image from the frame before.
- You can pan and zoom on the figure using the buttons at the bottom.
- The current DQ plane is overlayed on the image, and its opacity can be adjusted with the "DQ Opacity" slider. 
- An ellipse will appear at the location of your mouse-click in red.
- The ellipse can be modified with the "Theta", "Semi-Major Axis" and "Semi-Minor Axis" sliders in order to match a given shower. 
- The "Mark Shower" Button will record the current ellipse coordinates and parameters. 
- A green ellipse will display showers that have been recorded with the "Mark Shower" button. It will only display on groups for which the shower will be masked (ie, all groups from first appearance to the end of the integration, and not in groups before the shower, or in subsequent integrations.)
- When finished marking showers, the "Save" button will add the new masks to the DQ array, and save a new "_showers_masked.fits" file. 
- The "Quit" button will close the window. 

Two additional notes that are important are as follows:
1. Bright showers are usually much larger in the first frame in which they appear, so it is usually better to adjust the size of the ellipse to match the shower's appearance in the next frame, to avoid masking too-many pixels. 
2. Showers will have positive impact on the ramp even after it becomes difficult to see the glow in the difference images. For this reason, showers are always masked until the end of the integration.
3. You may also use the output of the RAMP_FIT step ("_rate.fits" or "_rate_ints.fits") as input. This will mask the shower in the DQ plane of the rate image, which may result in lower Signal-to-Noise for those pixels when creating a mosaic if there are any unaffected pixels in the original ramp. However, this can be a good first pass for problem showers.
