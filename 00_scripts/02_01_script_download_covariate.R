#' ---
#' title: var - download
#' author: mauricio vancine
#' date: 2022-04-21
#' ---

# prepare r ---------------------------------------------------------------

# packages
library(parallelly)
library(doParallel)
library(foreach)

# options
options(timeout = 1e6) # aumenta o tempo de download

# download present
doParallel::registerDoParallel(parallelly::availableCores(omit = 2))

foreach::foreach(i=1:19) %dopar% {
    
    url <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio", i, "_1981-2010_V.2.1.tif")
    
    destfile <- paste0("CHELSA_bio", ifelse(i < 10, paste0("0", i), i), "_1981-2010_V.2.1.tif")
    
    download.file(url = url, destfile = destfile, mode = "wb")
    
}

doParallel::stopImplicitCluster()

# download future
doParallel::registerDoParallel(parallelly::availableCores(omit = 2))

for(i in c("2011-2040", "2041-2070")){
    
    for(j in c("MPI-ESM1-2-HR", "IPSL-CM6A-LR", "UKESM1-0-LL")){
        
        for(k in c("ssp370", "ssp585")){
            
            foreach::foreach(b=1:19) %dopar% {
                
                url <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/",
                              i, "/", j, "/", k, "/bio/CHELSA_bio", b, "_", i, "_", tolower(j), "_", k, "_V.2.1.tif")
                
                destfile <- paste0("CHELSA_bio", ifelse(b < 10, paste0("0", b), b), "_", i, "_", tolower(j), "_", k, "_V.2.1.tif")
                
                download.file(url = url, destfile = destfile, mode = "wb")
                
            }
            
        }
        
    }
    
}
doParallel::stopImplicitCluster()

# end ---------------------------------------------------------------------
