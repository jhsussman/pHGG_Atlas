#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
from skimage import io
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
import tifffile
import imageio
import glob
import re
import math 
from PIL import Image
from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import make_outline_overlay
from deepcell.utils.plot_utils import create_rgb_image

# Replace 'path/to/file' with the appropriate directory path
inputDir = 'path/to/file'

# Input OME.TIFF file and number of divisions for tiling
inOmeTiff = os.path.join(inputDir, 'output_nuclear_panmembrane.tif')
ndivides = 2

# Reading and tiling image
print("Reading Image")
im = tifffile.imread(inOmeTiff)
print(im.shape)

print("Tiling Image") 
M = im.shape[1] // ndivides
N = im.shape[2] // ndivides
tiles = [im[:, x:x + M, y:y + N] for x in range(0, im.shape[1], M) for y in range(0, im.shape[2], N)]

# Save tiled images
num = 1
for tile in tiles:
    tifffile.imwrite(os.path.join(inputDir, "Pixel_4", f"{num}_nuclear_panmembrane.qptiff"), tile)
    print(f"Printing tile: {num}")
    print(tile.shape)
    num += 1

# Function to run Mesmer on a tiled image
def run_mesmer(inputDir, filein, tilenum): 
    outDir = os.path.join(inputDir, 'Pixel_4')
    sampleName = os.path.basename(inputDir)

    print(f'Input OME.TIFF file: {inOmeTiff}')
    outFileCell = os.path.join(outDir, f"{tilenum}_segmentation_cell.tif")
    outFileNuclear = os.path.join(outDir, f"{tilenum}_segmentation_nuclear.tif")
    outFileNuclearDil = os.path.join(outDir, f"{tilenum}_segmentation_nuclear_dil.tif")
    
    app = Mesmer()
    print('Training Resolution:', app.model_mpp, 'microns per pixel')
    im = io.imread(filein)
    X_train = np.transpose(im, (1, 2, 0))
    X_train = X_train.reshape(1, X_train.shape[0], X_train.shape[1], X_train.shape[2])
    rgb_images = create_rgb_image(X_train, channel_colors=['green', 'blue'])
    
    print('Starting Dilated Nuclear Segmentation')
    segmentation_predictions_nuc_dil = app.predict(X_train, image_mpp=0.5, compartment='nuclear', postprocess_kwargs_nuclear={'pixel_expansion': 4, 'interior_threshold': 0.3, 'maxima_threshold': 0.3})
    im = Image.fromarray(segmentation_predictions_nuc_dil[0, :, :, 0]).save(outFileNuclearDil)
    overlay_data_nuc_dil = make_outline_overlay(rgb_data=rgb_images, predictions=segmentation_predictions_nuc_dil)
    fig, ax = plt.subplots(1, 2, figsize=(15, 15))
    ax[0].imshow(rgb_images[0, ...])
    ax[1].imshow(overlay_data_nuc_dil[0, ...])
    ax[0].set_title('Raw data')
    ax[1].set_title('Dilated Nuclear Predictions')
    print("Saving overlay")
    fig.savefig(os.path.join(outDir, f"{tilenum}_segmentation_nuclear_dilated_overlay.tif"), dpi=600)
    with open(os.path.join(outDir, f"{tilenum}_segmentation_predictions_dilated_nuclear.npy"), "wb") as outnp: 
        np.save(outnp, segmentation_predictions_nuc_dil)

# Segmenting each tiled image using Mesmer
numimages = int(math.pow(ndivides, 2))
for j in range(1, numimages + 1): 
    imagetilepath = os.path.join(inputDir, "Pixel_4", f"{j}_nuclear_panmembrane.qptiff")
    print(f"Segmenting tile: {j}")
    run_mesmer(inputDir, imagetilepath, j) 

# Merging dilated nuclear numpy mask
print("Merging dilated nuclear numpy mask")
extension = '_segmentation_predictions_dilated_nuclear.npy'
imagetilearray = []
rowimagearray = []
curind = 0 
for k in range(1, numimages + 1): 
    print(curind)
    mask = np.load(os.path.join(inputDir, "Pixel_4", f"{k}{extension}"))
    maximum = mask.max()
    with np.nditer(mask, op_flags=['readwrite']) as it:
        for x in it: 
            if x > 0: 
                x[...] = x + curind
    imagetilearray.append(mask)
    curind = curind + maximum
ind = 0
for l in range(ndivides): 
    x = ind
    y = ind + ndivides
    newrow = np.concatenate((imagetilearray[x:y]), axis=2)
    ind = ind + ndivides
    rowimagearray.append(newrow)
merged_dilnuc_mask = np.concatenate((rowimagearray), axis=1)
with open(os.path.join(inputDir, "Pixel_4", "segmentation_predictions_dilated_nuclear.npy"), "wb") as outnp: 
    np.save(outnp, merged_dilnuc_mask)

# Merging dilated nuclear tiff mask
print("Merging dilated nuclear tiff mask")
extension = '_segmentation_nuclear_dil.tif'
imagetilearray = []
rowimagearray = []
curind = 0 
for k in range(1, numimages + 1): 
    print(curind)
    im = tifffile.imread(os.path.join(inputDir, "Pixel_4", f"{k}{extension}"))
    maximum = im.max()
    with np.nditer(im, op_flags=['readwrite']) as it:
        for x in it: 
            if x > 0: 
                x[...] = x + curind
    imagetilearray.append(im) 
    curind = curind + maximum
ind = 0
for l in range(ndivides): 
    x = ind
    y = ind + ndivides
    newrow = np.concatenate((imagetilearray[x:y]), axis=1)
    ind = ind + ndivides
    rowimagearray.append(newrow)
merged_dilnuc_mask = np.concatenate((rowimagearray), axis=0)
tifffile.imwrite(os.path.join(inputDir, "Pixel_4", "segmentation_nuclear_dil.tif"), merged_dilnuc_mask)
