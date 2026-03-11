library(tools) 
library(readr)

#Generate Step 2 Scripts ------


setwd("/mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances/")
instance40<-read_csv("NolanInstances.csv")
samples<-instance40$Sample_Name
 

if (!file.exists("Step1aScripts")) {
  dir.create("Step1aScripts")
} else {
  message(sprintf("Folder '%s' already exists. Skipping creation.", "Step1aScripts"))
}

output_directory <- "/mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances/Step1aScripts"
setwd(output_directory)
for (sample in samples) {
  script_name <- file.path(output_directory, paste0(sample, "_.sh"))
  
  # Write the script content to the file
  writeLines(
    c(
      "#!/bin/bash",
      paste0("python /mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances/Step1a_N-Orbit-Enumerate-HPC.py ", sample)
    ),
    script_name
  )
  
  # Make the script executable (equivalent to chmod +x in Bash)
  Sys.chmod(script_name, "755")
}

setwd("/mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances/")
if (!file.exists("Step2aScripts")) {
  dir.create("Step2aScripts")
} else {
  message(sprintf("Folder '%s' already exists. Skipping creation.", "Step2aScripts"))
}
output_directory <- "/mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances/Step2aScripts/"
setwd(output_directory)
for (sample in samples) {
  script_name <- file.path(output_directory, paste0(sample, "_.sh"))
  
  # Write the script content to the file
  writeLines(
    c(
      "#!/bin/bash",
      paste0("python /mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances/Step2a_Neighborhood-Distances.py ", sample)
    ),
    script_name
  )
  
  # Make the script executable (equivalent to chmod +x in Bash)
  Sys.chmod(script_name, "755")
}


setwd("/mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances/")
if (!file.exists("n-OrbitEnrichmentScripts")) {
  dir.create("n-OrbitEnrichmentScripts")
} else {
  message(sprintf("Folder '%s' already exists. Skipping creation.", "n-OrbitEnrichmentScripts"))
}
output_directory<-("/mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances/n-OrbitEnrichmentScripts/")
setwd(output_directory)
for (i in 1:400) {
  script_name <- file.path(output_directory, paste0("norbitEnrich",i, ".sh"))
  
  # Write the script content to the file
  writeLines(
    c(
      "#!/bin/bash",
      paste0("python /mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances/n-orbit-enrichment.py ", i)
    ),
    script_name
  )
  
  # Make the script executable (equivalent to chmod +x in Bash)
  Sys.chmod(script_name, "777")
}




