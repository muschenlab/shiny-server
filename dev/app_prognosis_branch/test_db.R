# install.packages("RMariaDB")
library(DBI)
library(RMariaDB)

## using mark account for setup - shiny only has select rights
db <- dbConnect(RMariaDB::MariaDB(), user = "mark", password = pswd,
                 dbname = "prognosis", host = "localhost")
study_info <- paste0("CREATE TABLE study_info (",
                     "study_id VARCHAR(20) PRIMARY KEY,",
                     "study_name VARCHAR(20),",
                     "species ENUM('human','mouse'),",
                     "assay_type ENUM('expr_array','RNAseq'),",
                     "platform ENUM('Illumina','HG-U133_Plus_2','HG-U133'),",
                     "pmid SMALLINT",
                     "description VARCHAR(200));")
dbSendQuery(db, study_info)
sample_info <- paste0("CREATE TABLE sample_info (",
                      "sample_id VARCHAR(20) PRIMARY KEY,",
                      "subject_id VARCHAR(20),",
                      "study_id VARCHAR(20),",
                      "age_years DECIMAL(5,2),",
                      "sex ENUM('M','F'),",
                      "treatment VARCHAR(20),",
                      "sampling_time ENUM('diagnosis','relapse'),",
                      "cell_type VARCHAR(20),",
                      "cell_subtype VARCHAR(20),",
                      "cell_subtype2 VARCHAR(20),",
                      "tissue_type VARCHAR(20),",
                      "disease_state VARCHAR(10),",
                      "days_to_death SMALLINT,",
                      "censor_death ENUM('dead','alive'),",
                      "days_to_relapse SMALLINT,",
                      "censor_relapse ENUM('relapsed','non-relapsed')",
                      "additional_info VARCHAR(200));")
dbSendQuery(db, sample_info)
expr_data <- paste0("CREATE TABLE expr_data (",
                    "feature_id VARCHAR(20),",
                    "sample_id VARCHAR(20),",
                    "expr_val DECIMAL(10,6),",
                    "PRIMARY KEY(feature_id, sample_id));")
dbSendQuery(db, expr_data)
dbListTables(db)

# edit and add data from existing

dbDisconnect(db)

# connect for querying
dbname <- "prognosis"
cnf <- list.files("/poolio/db", pattern = paste0(dbname, ".cnf"), full.names = T)
db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
dbListTables(db)
dbDisconnect(db)
