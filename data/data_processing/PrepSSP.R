### Prep SSP Data
# Downloaded from https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=60 June 09, 2022
# SspDb_country_data_2013-06-12.csv.zip

## ----------Load Libraries ----------------------------------
library(reshape2)
library(tidyverse)

## ----------Load Data ---------------------------------------
ciamdir <- getwd()
ciamregions <- read_csv(paste0(ciamdir,"/../input/xsc.csv")) %>% select(rgn) %>% unique()
ciampop <- read_csv(paste0(ciamdir,"/../input/pop.csv"))
ciamypc <- read_csv(paste0(ciamdir,"/../input/ypcc.csv"))

ciampopmelt <- ciampop %>% rename(time=variable)%>% melt(id.vars=c("time")) %>% mutate(variable=as.character(variable)) %>%
  rename(regions=variable) %>% mutate(year = 2000 + time * 10 )

ciamypcmelt <- ciamypc %>% rename(time=variable)%>% melt(id.vars=c("time")) %>% mutate(variable=as.character(variable)) %>%
  rename(regions=variable) %>% mutate(year = 2000 + time * 10 )

ssp <- read_csv(paste0(ciamdir, "/SspDb_country_data_2013-06-12.csv")) %>% 
  melt(id.vars=c("MODEL","SCENARIO","REGION","VARIABLE","UNIT")) %>% rename(YEAR=variable) %>% mutate(YEAR=as.numeric(YEAR))

col_order <- names(ciampop)
cpi <- 1.14 # 2010$/2005$ # https://www.bls.gov/data/inflation_calculator.htm

# SSP years are 2005-2100
ssp_filt <- ssp %>% filter(VARIABLE %in% c("Population","GDP|PPP")) %>% 
  mutate(YEAR=1945+YEAR*5) %>% filter(YEAR>=2010,YEAR<=2100, YEAR %% 10 ==0)

## Assumption 1: Regions not represented by SSPs will use original CIAM data
ciam_rgns <- ciamregions$rgn

### Process SSP data and save 
POP_MODELS <- unique(ssp_filt$MODEL[ssp_filt$VARIABLE =="Population"])
YPC_MODELS <- unique(ssp_filt$MODEL[ssp_filt$VARIABLE =="GDP|PPP"]) # Pik GDP-32 is only USA, so discard

for (P in POP_MODELS){
  s <- ssp_filt %>% filter(VARIABLE=="Population", MODEL==P)
  SSP <- unique(s$SCENARIO)
  for (S in SSP){
    s2 <- s %>% filter(SCENARIO==S)
    spopmelt <- s2 %>% select(REGION,YEAR,value) 
    spop <- s2 %>% select(REGION,YEAR,value) %>% mutate(value=as.numeric(value))%>% 
      rename(variable=YEAR) %>% dcast(variable~REGION) %>%
      mutate(variable=(variable - 2000)/10 )
    
    ssp_rgns <- unique(spopmelt$REGION)
    ciam_excl_rgn <- ciam_rgns[!(ciam_rgns %in% ssp_rgns)]
    ciam_excl_pop <- ciampop[c("variable",ciam_excl_rgn)]
    
    ## Assumption 2: Population and GDP held constant after 2100
    s2100 <- spop %>% filter(variable==10)
    s2100 <- s2100[rep(1,10),] 
    s2100$variable <- seq(11,20,1)
    spop <- rbind(spop,s2100)
    spop <- merge(spop,ciam_excl_pop,by = c("variable"))
    spop <- spop[names(spop) %in% col_order]
    spop <- spop[,col_order]
    
    name<- gsub(" ","",paste0(P,"_",S))
    
    assign(x = name,spopmelt)
    if (P == "IIASA GDP"){ # only save IIASA GDP for MimiCIAM
      write_csv(spop,path = paste0(ciamdir,"/../ssp/pop_",name,".csv"))
    }
  }
}

for (Y in YPC_MODELS){
  # Convert GDP from billion US$2005/yr to $2010/yr
  # Conversion to Per Capita will use corresponding population model 
  s <- ssp_filt %>% filter(VARIABLE=="GDP|PPP", MODEL==Y) %>% mutate(value = as.numeric(value) * cpi, UNIT="billion US$2010/yr")
  SSP <- unique(s$SCENARIO)
  for (S in SSP){
    s2 <- s %>% filter(SCENARIO==S) %>% 
      mutate(Pop = ssp_filt$value[match(paste0(MODEL,SCENARIO,REGION,"Population",YEAR),paste0(ssp_filt$MODEL,ssp_filt$SCENARIO,ssp_filt$REGION,ssp_filt$VARIABLE,ssp_filt$YEAR))],
             value = (value * 1e3) / as.numeric(Pop), UNIT="GDP per Capita (2010$)") %>% na.omit()
    sypcmelt <- s2 %>% select(REGION,YEAR,value) 
    
    sypc <- s2 %>% select(REGION,YEAR,value) %>% mutate(value=as.numeric(value))%>% 
      rename(variable=YEAR) %>% dcast(variable~REGION) %>%
      mutate(variable=(variable - 2000)/10 )
    
    ssp_rgns <- unique(sypcmelt$REGION)
    ciam_excl_rgn <- ciam_rgns[!(ciam_rgns %in% ssp_rgns)]
    ciam_excl_ypc <- ciamypc[c("variable",ciam_excl_rgn)]
    
    ## Assumption 2: Population and GDP held constant after 2100
    ##  to do: possibly calculate average growth rate and apply 
    s2100 <- sypc %>% filter(variable==10)
    s2100 <- s2100[rep(1,10),] 
    s2100$variable <- seq(11,20,1)
    sypc <- rbind(sypc,s2100)
    sypc <- merge(sypc,ciam_excl_ypc,by = c("variable"))
    sypc <- sypc[names(sypc) %in% col_order]
    sypc <- data.frame(sapply(sypc, as.numeric))
    sypc <- sypc[,col_order]
    
    
    usa <- sypc["USA"] %>% rename(ypc_usa=USA)
    
    name<- gsub(" ","",paste0(Y,"_",S))
    
    assign(x = name,sypcmelt)
    if (Y == "IIASA GDP"){ # only save IIASA GDP for MimiCIAM
      write_csv(sypc,path = paste0(ciamdir,"/../ssp/ypcc_",name,".csv"))
      # don't need the ypc_usa files for MimiCIAM
      # write_csv(usa, path = paste0(ciamdir,"/../ssp/ypc_usa_",name,".csv"))
    }
  }
}
