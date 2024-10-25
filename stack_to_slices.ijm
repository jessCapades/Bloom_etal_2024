// Set working directory 
home_directory = getDirectory("Choose a Directory");
edited_image_folder = home_directory + File.separator + "edited_images";
gfp_image_folder = edited_image_folder + File.separator + "gfp_images";
dic_image_folder = edited_image_folder + File.separator + "dic_images";
edited_image_folder = File.makeDirectory(edited_image_folder);
gfp_destination_directory = File.makeDirectory(gfp_image_folder);
dic_destination_directory = File.makeDirectory(dic_image_folder);

// import images
image_list = getFileList(home_directory);
image_list_length = lengthOf(image_list);


for (i = 0; i < image_list.length; i++) {
    if (endsWith(image_list[i], ".tif")) { 
    	full_file_path = home_directory + image_list[i];
        open(home_directory + File.separator + image_list[i]);
    }
}

// Crop embryos with specified shape and split imaging channels  
title = "Adjust ROI";
message = "fit ROI over embryo";
open_image_list = getList("image.titles");

current_image_stack = open_image_list[0];
selectImage(current_image_stack);
current_image_width = getWidth();
current_image_height = getHeight();

for (k = 0; k < open_image_list.length; k++){
    current_stack = open_image_list[k];
	selectImage(current_stack);
	roiManager("select", 0);	
	waitForUser(title, message);
    setBackgroundColor(0, 0, 0);
    run("Clear Outside", "stack");
    channels = Stack.getDimensions(image_width, image_height, image_channels, image_slices, imaeg_frames);
    if (image_channels == 2){
    	run("Split Channels");
    }else {
    	print("This is a single channel image");
    }
}
// Save separated channel stacks 

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
			
		




