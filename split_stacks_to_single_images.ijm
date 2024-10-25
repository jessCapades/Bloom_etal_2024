// split stack images to single images
home_directory = getDirectory("Choose a Directory");
single_slices_folder = home_directory + File.separator + "single_slices";
single_slices_directory = File.makeDirectory(single_slices_folder);

image_list = getFileList(home_directory);
image_list_length = lengthOf(image_list);

for (i = 0; i < image_list.length; i++) {
   if (endsWith(image_list[i], ".tif")) {
   	full_file_path = home_directory + image_list[i];
    open(home_directory + File.separator + image_list[i]);
   	run("Stack to Images");
   	
   	split_image_list = getList("image.titles");
   	output_directory = getDirectory("Choose a Directory");
   	for (j = 0; j < split_image_list.length; j++){
    		current_image_slice = split_image_list[j];
    		selectImage(current_image_slice);
    		current_image_slice_title = getTitle();
			image_out = output_directory + current_image_slice_title;
			saveAs("Tiff", image_out);
			close();
   }
   }
}

