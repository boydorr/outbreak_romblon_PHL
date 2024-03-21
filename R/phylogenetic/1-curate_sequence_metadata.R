## Curation of sequence metadata 
library(parsedate)
library(lubridate)
library(stringr)
library(dplyr)
library(tidyr)
## Data sources: RABV Philippines data from Genbank (RABV-GLUE), Essel PhD seq (in-house) and PCG-M (pulled from paper) (not currently in GLUE dataset but publically available)

## import data files

pgc=read.csv("raw_data/pgc_data/pgcM_metadata_49.csv", strip.white =T)
phd=read.csv("raw_data/phd_data/ph_wgs_metadata_29.csv", strip.white =T)
gb=read.csv("raw_data/rabv-glue_data/ph_metadata_615.csv", strip.white =T)


## initiate new datafame to collect relevant metadata

# pgc
pgc.curated=data.frame(Sample.ID=pgc$Sample.ID, NCBI.accession=pgc$Accession, Sample.collection.date=pgc$Sample.collection.year, Host.species=pgc$Host.species, Region="Davao Region", Province=NA,Municipality=pgc$Municipality.or.City.of.origin,Barangay=NA, Location=NA, source="genbank")

pgc.curated$corrected_date=parse_date(pgc.curated$Sample.collection.date)
pgc.curated$corrected_date=as.Date(pgc.curated$corrected_date)
pgc.curated$decimal_date=decimal_date(pgc.curated$corrected_date)

# gb
## first, concatenate data that has a range of dates associated with submission, rather than a precise date for individual seqs
gb$date_range=NA
for (i in 1:nrow(gb)){
  if(gb$sequence.earliest_collection_year[i]!=gb$sequence.latest_collection_year[i] & gb$sequence.earliest_collection_year[i]!="-"){
    gb$date_range[i]=paste(gb$sequence.earliest_collection_year[i],gb$sequence.latest_collection_year[i], sep="-")}else{
      gb$date_range[i]= gb$sequence.latest_collection_year[i]
    }
  }


# check that this matches the sequence collection year column. If it does then can use collection year col, if it doesn't then use the new date_range col
if(sum(!gb$date_range %in% gb$sequence.latest_collection_year)==0){
  gb.curated.date=gb$sequence.latest_collection_year}else{
    gb.curated.date=gb$date_range
  }

gb.curated=data.frame(Sample.ID=gb$sequence.isolate, NCBI.accession=gb$sequence.sequenceID, Sample.collection.date=gb.curated.date, Host.species=gb$sequence.host, Region=NA,Province=NA, Municipality=NA,Barangay=NA, Location=gb$sequence.gb_place_sampled,source="rabv-glue")
gb.curated <- gb.curated %>%
  separate(Location, c('Region', 'Province'), sep = ",") %>%
  mutate(
    Region = trimws(Region),
    Province = trimws(Province),
    Location=gb$sequence.gb_place_sampled
  )

gb.curated$corrected_date=parse_date(gb.curated$Sample.collection.date)
gb.curated$corrected_date=as.Date(gb.curated$corrected_date)
gb.curated$decimal_date=decimal_date(gb.curated$corrected_date)


# essel phd
phd.curated=data.frame(Sample.ID=phd$ID, NCBI.accession=NA, Sample.collection.date=phd$date, Host.species=phd$species,Region=phd$region, Province=phd$province, Municipality=NA,Barangay=phd$brgy,Location=phd$place,source="this paper")

#phd.curated$corrected_date=parse_date(phd.curated$Sample.collection.date)
phd.curated$corrected_date=as.Date(phd.curated$Sample.collection.date, format="%d-%b-%y")
phd.curated$decimal_date=decimal_date(phd.curated$corrected_date)

## combine all of curated data into one dataframe
curated.df=rbind(pgc.curated, phd.curated, gb.curated)
curated.df=curated.df[which(curated.df$Sample.ID!=""),]
rawseq=nrow(curated.df)

## export the curated dataset
write.csv(curated.df, paste0("processed_data/data_prep/metadata/ph_curated_data_",rawseq,".csv"), row.names=F)

