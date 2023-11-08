
importCSVs <- function(path) {
  # list all csv files from the current directory
  list.files(path = path, pattern=".csv$")
  # create a list from these files
  list.filenames<-list.files(path = path, pattern=".csv$")
  # create an empty list that will serve as a container to receive the incoming files
  list.data<-vector("list", length(list.filenames))
  
  # create a loop to read in your data
  for (i in 1:length(list.filenames)) {
    list.data[[i]]<-read.csv(paste(path,list.filenames[i],sep = ""), header=TRUE)
  }
  
  # add the names of your data to the list
  names(list.data)<-list.filenames
  
  return(list.data)
  
}
