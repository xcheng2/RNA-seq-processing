rm(list = ls())
clinical = read_tsv("clinical.tsv")


location = clinical[, c("cases.submitter_id", "diagnoses.site_of_resection_or_biopsy")]
location = location[!duplicated(location),]

t = table(location$cases.submitter_id)
t[t==2]

location = location[location$diagnoses.site_of_resection_or_biopsy != "Not Reported",]

location$temp_location = ""
location$temp_location[location$diagnoses.site_of_resection_or_biopsy == "Cecum"] = "Right"
location$temp_location[location$diagnoses.site_of_resection_or_biopsy == "Ascending colon"] = "Right"
location$temp_location[location$diagnoses.site_of_resection_or_biopsy == "Hepatic flexure of colon"] = "Right"
location$temp_location[location$diagnoses.site_of_resection_or_biopsy == "Transverse colon"] = "Right"
location$temp_location[location$diagnoses.site_of_resection_or_biopsy == "Splenic flexure of colon"] = "Left"
location$temp_location[location$diagnoses.site_of_resection_or_biopsy == "Descending colon"] = "Left"
location$temp_location[location$diagnoses.site_of_resection_or_biopsy == "Sigmoid colon"] = "Left"
location$temp_location[location$diagnoses.site_of_resection_or_biopsy == "Rectosigmoid junction"] = "Left"
location$temp_location[location$diagnoses.site_of_resection_or_biopsy == "Rectum, NOS"] = "Left"


write.csv(location, "<filename.csv>")