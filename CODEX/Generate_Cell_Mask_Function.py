import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tifffile

sample_names = ["7622", "6477", "5928", "4740", "4337", "3058", "942", "371", "339", "161", "5335"]

for name in sample_names:
	print(name)
	clusters = 'Annotation_Masks/Filtered_'+name+'_annotations_lv2_Pix4_10172023.csv' 
	mask = '/mnt/isilon/tan_lab/parvaresha/Pediatric_Glioma/pHGG_Oldridge_Samples/'+name+'/mesmer/Pix4/segmentation_nuclear_dil.tif'
	output = 'Annotation_Masks/FULL_Filtered_'+name+'_Pix4_Mask.tiff'
	print('Reading Clusters')
	clusters = pd.read_csv(clusters)
	print('Reading Mask')
	mask = tifffile.imread(mask)

	# Initialize numpy array from segmentation mask
	print('Initializing numpy array')
	masks_out = mask.copy()
	masks_out = np.array(masks_out).astype(int)
	print(masks_out.dtype)

	#Convert those cells that have been filtered out to zero
	print('Filtering cells')
	filter_mask = np.isin(masks_out,clusters['CellID'], invert=True)
	masks_out[filter_mask] = 0

	# Create dictionary mapping Cell ID keys to cluster values
	print('Creating mapping')
	mapping = dict(zip(clusters["CellID"], clusters["Annotation"]))

	# Define function to change the values in input array by key-value pairs from dict
	def mp(entry):
	    return mapping[entry] if entry in mapping.keys() else entry

	print('Vectorizing array')
	mp = np.vectorize(mp)

	print('Converting Values') 
	out = mp(masks_out)

	print('Saving mask')
	tifffile.imwrite(output, out)
