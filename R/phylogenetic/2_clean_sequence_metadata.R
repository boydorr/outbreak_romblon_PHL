library(parsedate)
library(lubridate)
library(stringr)
library(dplyr)
library(tidyr)

# raw data
all_phl.annot=read.csv("processed_data/data_prep/metadata/ph_curated_data_693.csv")

# correct a few sample names to match sequences
all_phl.annot$Sample.ID[which(all_phl.annot$Sample.ID=="H-23-011Sk12")]="H-23-011Sk_12"

# merge adm level names (typos, slightly different names etc)
spatial=read.csv("raw_data/spatial/PHL_centroids.csv")
regions=c("National Capital Region","Cordillera Administrative Region","Ilocos Region","Cagayan Valley","Central Luzon","Calabarzon","Southwestern Tagalog Region",	"Bicol Region","Western Visayas","Central Visayas","Eastern Visayas","Zamboanga Peninsula","Northern Mindanao","Davao Region","Soccsksargen","Caraga","Bangsamoro")

# fill all empty cells with NAs (as some already listed as NA)
all_phl.annot$Region[all_phl.annot$Region==""]=NA;all_phl.annot$Province[all_phl.annot$Province==""]=NA;all_phl.annot$Barangay[all_phl.annot$Barangay==""]=NA
all_phl.annot$Region[all_phl.annot$Region=="-"]=NA;all_phl.annot$Province[all_phl.annot$Province=="-"]=NA;all_phl.annot$Barangay[all_phl.annot$Barangay=="-"]=NA

# regions:
# format
all_phl.annot$Region=str_to_title(all_phl.annot$Region)
# standardise names
all_phl.annot$Region[grepl("Region 1$|Regoin 1$|Region1$|R1$|\\bI$", all_phl.annot$Region,ignore.case = T)]="Ilocos Region"
all_phl.annot$Region[grepl("Region 2$|Regoin 2$|Region2$|R2$|\\bII$", all_phl.annot$Region,ignore.case = T)]="Cagayan Valley"
all_phl.annot$Region[grepl("Region 3$|Regoin 3$|Region3$|R3$|\\bIII$", all_phl.annot$Region,ignore.case = T)]="Central Luzon"
all_phl.annot$Region[grepl("Region 4a$|Regoin 4a$|Region4a$|R4a$|IVa|IV-a$|\tCalabarzon|Angono", all_phl.annot$Region,ignore.case = T)]="Calabarzon"
all_phl.annot$Region[grepl("Region 4b$|Regoin 4b$|Region4b$|Region Iv$|R4b$|IVb|IV-b$", all_phl.annot$Region,ignore.case = T)]="Southwestern Tagalog Region"
all_phl.annot$Region[grepl("Region 5$|Regoin 5$|Region5$|R5$|\\bV$|\tBicol Region", all_phl.annot$Region,ignore.case = T)]="Bicol Region"
all_phl.annot$Region[grepl("Region 6$|Regoin 6$|Region6$|R6$|\\bVI$", all_phl.annot$Region,ignore.case = T)]="Western Visayas"
all_phl.annot$Region[grepl("Region 7$|Regoin 7$|Region7$|R7$|\\bVII$|\tCentral Visayas|Valencia", all_phl.annot$Region,ignore.case = T)]="Central Visayas"
all_phl.annot$Region[grepl("Region 10$|Regoin 10$|Region10$|R10$|\\bX$|Mindanao Island", all_phl.annot$Region,ignore.case = T)]="Northern Mindanao"
all_phl.annot$Region[grepl("Region 11$|Regoin 11$|Region11$|R11$|\\bXI$", all_phl.annot$Region,ignore.case = T)]="Davao Region"
all_phl.annot$Region[grepl("car$|Cordillera|Bontoc", all_phl.annot$Region,ignore.case = T)]="Cordillera Administrative Region"
all_phl.annot$Region[grepl("ncr$|Caloocan City|Pasay|Quezon|San Miguel", all_phl.annot$Region,ignore.case = T)]="National Capital Region"

# note one sequence is labelled as "Mindanao" which could be northern mindanao region or elsewhere on the island so left as is
unique(all_phl.annot$Region)[which(!unique(all_phl.annot$Region) %in% regions)]

# province:
# format
all_phl.annot$Province=str_to_title(all_phl.annot$Province)
all_phl.annot$Province=gsub("Del","del",all_phl.annot$Province)
# standardise names 
all_phl.annot$Province[all_phl.annot$Province=="Davaodel Sur"]="Davao del Sur"
all_phl.annot$Province[all_phl.annot$Province=="Davaodel Norte"]="Davao del Norte"
all_phl.annot$Province[all_phl.annot$Province=="Metro Manila"]="Metropolitan Manila"
all_phl.annot$Province[all_phl.annot$Province=="Ilocossur"]="Ilocos Sur"
all_phl.annot$Province[all_phl.annot$Province=="Cotabato"]="North Cotabato" #is no south cotabato!
all_phl.annot$Province[all_phl.annot$Province=="Kalinga Apayao"]="Kalinga" #old name
all_phl.annot$Province[all_phl.annot$Province=="Isabela City"]="Isabela" #interchangeable
all_phl.annot$Province[all_phl.annot$Province=="Misamis Orientale"]="Misamis Oriental"
all_phl.annot$Province[all_phl.annot$Province=="Ifgao"]="Ifugao"
all_phl.annot$Province[all_phl.annot$Province=="Launion"]="La Union"

# check if any further fixes required: 
unique(all_phl.annot$Province)[which(!unique(all_phl.annot$Province) %in% spatial$Loc_ID)]
nrow(all_phl.annot)
nseq=nrow(all_phl.annot)

# add in the romblon annotations
all_phl.annot$outbreak=NA
all_phl.annot$outbreak[which(all_phl.annot$Province=="Romblon" & all_phl.annot$decimal_date>2022)]="Romblon_new"
all_phl.annot$outbreak[which(all_phl.annot$Province=="Romblon" & all_phl.annot$decimal_date<2022)]="Romblon_old"
all_phl.annot$outbreak[which(all_phl.annot$Province!="Romblon"|is.na(all_phl.annot$Province))]="Other"

write.csv(all_phl.annot,paste0("processed_data/data_prep/metadata/ph_metadata_",nseq,"_clean.csv"), row.names=F)
