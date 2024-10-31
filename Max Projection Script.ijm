input1 = "/Users/sashatsenter/Desktop/temp1";
output = "/Users/sashatsenter/Desktop/tempout";

// z_project function takes as inputs:
// input = the directory containing unprocessed images
// output = the directory where processed images should be saved
// filename = the name of the file to be z-projected
// for the output file, appends _zproj to file name
function MaxZ(input1, filename1) {
	run("Bio-Formats Windowless Importer", "open="+input1+"/"+filename1+" autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT");
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", output + "/" + "MaxZ_1_"+filename1);
}

// loops through files in the directory
// applies the z_project function to each file
list1 = getFileList(input1);
for (i = 0; i < list1.length; i++) {
	currentfilename = list1[i];
	print("Now Processing: " + currentfilename);
	MaxZ(input1, currentfilename);
	run("Close All");
}