image_slice_directory = getDir("Choose a Directory");

image_list = getFileList(image_slice_directory);

for (i = 0; i < image_list.length; i++) {
    if (endsWith(image_list[i], ".tif")) { 
    	full_file_path = image_slice_directory + image_list[i];
        open(image_slice_directory + File.separator + image_list[i]);
    }
}

image_slice_names = getList("image.titles");

//set new brightness contrast values for each image & convert to RGB
for (i=0; i<image_slice_names.length; i++) {
	current_image_slice = image_slice_names[i];
	selectImage(current_image_slice);
	run("Brightness/Contrast...");
	run("Enhance Contrast", "saturation=0.45");
	getMinAndMax(min, max);
	new_minimum_pixel_value = min + 416;
	new_maximum_pixel_value = max - 500;
	setMinAndMax(new_minimum_pixel_value, new_maximum_pixel_value);
	run("RGB Color");
}

// combine images to stack and save
RGB_stack_directory = getDir("Choose a Directory");
run("Images to Stack");
selectWindow("Stack");
Dialog.create("User Input");
label = "new image name";
initialText = "input here";
Dialog.addString(label, initialText);
Dialog.show();
stack_title = Dialog.getString();
rename(stack_title);
new_renamed_stack = getTitle();
output_image_full_path = RGB_stack_directory + new_renamed_stack;
saveAs("Tiff", output_image_full_path);
close();
