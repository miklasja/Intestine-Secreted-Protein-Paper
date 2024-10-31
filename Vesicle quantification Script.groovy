//Count number and area of vesicles
//c. Elegans Ges-1p::F53B6.4-Linker-Scarlet
//Sasha Tsenter 


//Images: Max intensity projections.

//Run this script: Run > Run for project.

//Create full-image annotations for all z-slices
    //If needed, clear any previous annotations/detections: 
    selectAnnotations();
    clearSelectedObjects();
    selectDetections();
    clearSelectedObjects();
int t = 0;
int z = 0;
createSelectAllObject(true, z, t);

    
//Select all annotations for measurement steps
selectAnnotations();

//Use cell detection to identify vesicles.
    //Use nucleus parameters to segment vesicles. Cell area is calculated including "cell expansion".
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage":"Channel 1","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":1.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":0.3,"minAreaMicrons":1.0,"maxAreaMicrons":5.0,"threshold":1.5,"watershedPostProcess":true,"cellExpansionMicrons":0.6,"includeNuclei":false,"smoothBoundaries":true,"makeMeasurements":true}')

//Update image open in viewer:
//File > Reload data.

//Access measurements for entire project:
//Measure > Export measurements
// select "Annotations" for counts
// select "Detections" for cell level
//Save.

