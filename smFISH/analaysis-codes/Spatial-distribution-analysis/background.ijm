
setBatchMode(true);
indir = getDirectory("choose image directory");
list = getFileList(indir);
File.makeDirectory(indir + "3D_objects");//run("ImageJ2...", "scijavaio=true loglevel=INFO");

// Get data dir and create output dirs
//indir = getDirectory("choose image directory");
//list = getFileList(indir);

// generate mask of nuclei
// set up loop to get segmented images
for (i=0; i<list.length; i++) { 
	if(endsWith(list[i],"MLE.vsi.vsi")){
		showProgress(i+1, list.length);
		print("processing ... "+i+1+"/"+list.length+"\n         "+list[i]);
		image_name = "MAX_"+list[i];
		path=indir+list[i];
		
		// Open file
    	//open(path);
		run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
		run("Z Project...", "projection=[Max Intensity]");
		
		original = getTitle();
		run("Arrange Channels...", "new=1");

		setAutoThreshold("Huang dark"); //Percentile
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Divide...", "value=255");

		// save image
		saveAs("Tiff", indir+"3D_objects/"+"Background_mask_of_"+image_name);
		//saveAs("Tiff", indir+"cal/"+list[i]);
		
	}
run("Close All");
}

setBatchMode(false);