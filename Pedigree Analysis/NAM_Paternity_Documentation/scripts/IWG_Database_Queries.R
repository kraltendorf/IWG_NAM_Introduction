### Original Data Intermediate Wheatgrass Database Queries

source('./scripts/Function_R.r') #load functions 

#make connection to the database
log_detail <- getLoginDetails() #fill in username and password in popup box

require(RMySQL) #load required packages may need to install.packages("RMySQL") dependicies are DBI packaage as well

#Make a MySQL connection to beocat
iwg <- dbConnect(MySQL(), user = log_detail[1], dbname = '', host = '', password = log_detail[2], port = ) #run this line to connect to the database


###############GBS Querires#################################
#make GBS key using ifnormation from gbs_name query
key_q1 <- "SELECT
  gbs_iwg.flowcell as 'Flowcell', 
  gbs_iwg.lane as 'Lane', 
  barcodes_iwg.barcode as 'Barcode', 
  dna_iwg.sample_name as 'FullSampleName',
  dna_iwg.plate_id as 'LibraryPlate', 
  substring(dna_iwg.well_A01,1,1) as 'Row',
  substring(dna_iwg.well_A01,2,2) as 'Col',
  concat(gbs_iwg.gbs_id,barcodes_iwg.barcode) as 'LibraryPrepID',
  barcodes_iwg.well_A01 as 'BarcodeWell',
  gbs_iwg.enzyme as 'Enzyme',
  gbs_iwg.gbs_id as 'LibraryPlateID',
  dna_iwg.plate_id as 'DNA_Plate',
  dna_iwg.well_A01 as 'SampleDNA_Well',
  dna_iwg.external_id as 'External_ID',
  'thinopyrum' as 'Genus',
  'intermedium' as 'Species', 
  gbs_iwg.plexing, 
  gbs_iwg.project, 
  dna_iwg.sample_id,
  dna_iwg.tissue_id,
  barcodes_iwg.set,
  germplasm.male_parent as 'Male',
  germplasm.female_parent as 'Female',
  germplasm.family_name as 'Family'
FROM dna_iwg LEFT JOIN gbs_iwg ON gbs_iwg.dna_id = dna_iwg.plate_id
 LEFT JOIN plant ON dna_iwg.sample_name = plant.plant_id
left join germplasm on dna_iwg.tissue_id = germplasm.germplasm_id
 INNER JOIN barcodes_iwg ON dna_iwg.well_A01 = barcodes_iwg.well_A01 AND gbs_iwg.plexing LIKE barcodes_iwg.`set` 
WHERE gbs_iwg.gbs_id = 'GBS1427'
ORDER BY gbs_iwg.gbs_id, dna_iwg.well_01A ASC"

key1 <- dbGetQuery(iwg, key_q1) #get key from query

key_q2 <- "SELECT
  gbs_iwg.flowcell as 'Flowcell', 
  gbs_iwg.lane as 'Lane', 
  barcodes_iwg.barcode as 'Barcode', 
  dna_iwg.sample_name as 'FullSampleName',
  dna_iwg.plate_id as 'LibraryPlate', 
  substring(dna_iwg.well_A01,1,1) as 'Row',
  substring(dna_iwg.well_A01,2,2) as 'Col',
  concat(gbs_iwg.gbs_id,barcodes_iwg.barcode) as 'LibraryPrepID',
  barcodes_iwg.well_A01 as 'BarcodeWell',
  gbs_iwg.enzyme as 'Enzyme',
  gbs_iwg.gbs_id as 'LibraryPlateID',
  dna_iwg.plate_id as 'DNA_Plate',
  dna_iwg.well_A01 as 'SampleDNA_Well',
  dna_iwg.external_id as 'External_ID',
  'thinopyrum' as 'Genus',
  'intermedium' as 'Species', 
  gbs_iwg.plexing, 
  gbs_iwg.project, 
  dna_iwg.sample_id,
  dna_iwg.tissue_id,
  barcodes_iwg.set,
  germplasm.male_parent as 'Male',
  germplasm.female_parent as 'Female',
  germplasm.family_name as 'Family'
FROM dna_iwg LEFT JOIN gbs_iwg ON gbs_iwg.dna_id = dna_iwg.plate_id
 LEFT JOIN plant ON dna_iwg.sample_name = plant.plant_id
left join germplasm on dna_iwg.tissue_id = germplasm.germplasm_id
 INNER JOIN barcodes_iwg ON dna_iwg.well_A01 = barcodes_iwg.well_A01 AND gbs_iwg.plexing LIKE barcodes_iwg.`set` 
WHERE gbs_iwg.gbs_id = 'GBS1427' or gbs_iwg.gbs_id = 'GBS1510'
ORDER BY gbs_iwg.gbs_id, dna_iwg.well_01A ASC"

key2 <- dbGetQuery(iwg, key_q2) #get key from query

#save files
#write.table(key2, file = './data/Original_Data/UMN_Key2.txt', quote = FALSE, row.names = FALSE, sep = '\t')

#clean up
rm(key2, key_q2, key1, key_q1)

################## Pedigree Information Query ########################

germ_q <- 'SELECT germ.* FROM germplasm germ' #set up germplasm query

germ <- dbGetQuery(iwg, germ_q) #run germplasm query

#remove any C4, C5, C6, C7, C8, C9
germ <- germ[!grepl('^K0|^C7|^C8|^C9|^C4|^C5|^WGN|^M26', germ$germplasm_id), ]

#saveRDS(germ, file = './data/RObjects/Germ_Query.RDS')

dbDisconnect(iwg) #disconnect from database

#clean up
rm(germ, germ_q)



###################### End Filter Queries ####################

#clean up
rm( iwg,  log_detail)
