

get_functions = function (file_location, list_of_functions) {
  for (fn in list_of_functions) {
    source(paste(file_location, fn, sep = ""))
  }
}