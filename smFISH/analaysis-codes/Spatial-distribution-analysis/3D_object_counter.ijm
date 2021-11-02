
setBatchMode(true);
//run("ImageJ2...", "scijavaio=true loglevel=INFO");

// Get data dir and create output dirs
//indir = "/Users/joshtitlow/Downloads/OriginalFileDownload/";
indir = getDirectory("choose image directory");
list = getFileList(indir);

// set up loop to detect image files
for (i=0; i<list.length; i++) { 
	if(endsWith(list[i],".vsi")){
		showProgress(i+1, list.length);
		print("processing ... "+i+1+"/"+list.length+"\n         "+list[i]);
		image_name = "MAX_"+list[i];
		path=indir+list[i];
		
		// Open file
    	//open(path);
		run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
		run("Z Project...", "projection=[Max Intensity]");
		
		original = getTitle();
		run("Arrange Channels...", "new=3");

		run("3D Objects Counter", "threshold=140 slice=1 min.=100 max.=4194304 objects statistics");
		
		// save images
		saveAs("Results", indir+"3D_objects/"+image_name+".csv");
		saveAs("Tiff", indir+"3D_objects/"+"Surface_map_of_"+image_name);
		//saveAs("Tiff", indir+"cal/"+list[i]);
		
	}
run("Close All");
}


setBatchMode(false);
//run("Quit");
