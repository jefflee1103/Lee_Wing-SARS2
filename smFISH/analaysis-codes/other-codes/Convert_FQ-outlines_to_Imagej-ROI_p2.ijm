input = getDirectory("Choose FQ outline R script output Directory ");

list = getFileList(input);

run("ROI Manager...");

for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	path = input + list[i];
	run("XY Coordinates... ", "open=[" + path + "]");
	roiManager("Add");
}

roiManager("Show All");
roiManager("multi-measure measure_all");

selectWindow("ROI Manager");
roiManager("Delete");