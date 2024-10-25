// Set working directory 
home_directory = getDirectory("Choose a Directory");
edited_image_folder = home_directory + File.separator + "edited_images";
C1_image_folder = edited_image_folder + File.separator + "C1";
C2_image_folder = edited_image_folder + File.separator + "C2";
C3_image_folder = edited_image_folder + File.separator + "C3";
edited_image_folder = File.makeDirectory(edited_image_folder);
C1_destination_directory = File.makeDirectory(C1_image_folder);
C2_destination_directory = File.makeDirectory(C2_image_folder);
C3_destination_directory = File.makeDirectory(C3_image_folder);

// import images
image_list = getFileList(home_directory);
image_list_length = lengthOf(image_list);
title = "Adjust ROI";
message = "fit ROI over embryo";

for (i = 0; i < image_list.length; i++) {
    if (endsWith(image_list[i], ".tif")) { 
    	full_file_path = home_directory + image_list[i];
        open(home_directory + File.separator + image_list[i]);
        current_stack = image_list[i];
		selectImage(current_stack);
		roiManager("select", 0);	
		waitForUser(title, message);
    	setBackgroundColor(0, 0, 0);
    	run("Clear Outside", "stack");
    	channels = Stack.getDimensions(image_width, image_height, image_channels, image_slices, imaeg_frames);
    	if (image_channels > 1){
    	run("Split Channels");
    	}else {
    	print("This is a single channel image");
    	}
    output_directory = getDirectory("Choose a Directory");
	split_image_list = getList("image.titles");
		for (j = 0; j < split_image_list.length; j++){
    		current_single_channel_stack = split_image_list[j];
    		selectImage(current_single_channel_stack);
    		current_single_channel_stack_title = getTitle();
			image_out = output_directory + current_single_channel_stack_title;
			saveAs("Tiff", image_out);
			close();
    	}

 
} 
}
		




