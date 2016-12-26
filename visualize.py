########################################## PACKAGES ############################################################
##### The python script which is shared with you, needs SciPy and Pillow packages.
##### If you don't have any of these packages, follow the links below to use script for testing your own images.
##### https://pillow.readthedocs.io/en/3.4.x/installation.html
##### http://scipy.org/install.html
################################################################################################################

from PIL import Image
import numpy as np

# Set 1 for converting text to image
# Set 0 for converting image to 2D array representation 
# and output as a txt file
isTxtToImg = 1;

#########################
### From txt to image ###
#########################
# Takes "input.txt" and shows the image representation
# of the 2D array
if isTxtToImg:
	data = np.loadtxt("thresholded.txt")
	svimg = Image.fromarray(data.astype('uint8'))
	svimg.show()
else:
    #########################
    ### From image to txt ###
    #########################
    # Takes "photo.jpg" as input and creates 2D GrayScale (0-255) 
	# .txt output
	img = Image.open("photo.jpg")  # open colour image
	img = img.convert('L')  # convert image to grayscale
	file = open("OutFileName.txt", "w")
	txt = list(np.asarray(img))
	for i in txt:
		for j in i:
			file.write(str(j) + " ")
		file.write('\n')
	file.close()
