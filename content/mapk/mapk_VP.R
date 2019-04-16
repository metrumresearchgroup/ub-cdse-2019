library(mrgsolve)
library(tidyverse)
library(magrittr)
library(reshape2)

library(R.utils)
sourceDirectory("content/script/VPop/")

# Load simulation results (We are trying to match this Population)
sims <- readRDS("content/mapk/mapk_sims.RDS")
sims %<>% filter(label == "GDC")

# Load treatment regimens
reg <- readRDS("content/mapk/mapk_setup.RDS")
reg %<>% filter(label == "GDC")%>%pull(object)
reg <- reg[[1]]

model <- mread("mapk", 'content/mapk/', soloc = 'content/mapk')
m_params <- names(param(model))

# What parameters do we need to explore?
foo <- readRDS("content/mapk/s10vpop_pk.RDS")
(foo %>% summarise_all(sd)%>%gather()%>%filter(value!=0))%>%pull(key)
length((foo %>% summarise_all(sd)%>%gather()%>%filter(value!=0))%>%pull(key))
p_names <- foo %>% summarise_all(sd)%>%gather()%>%filter(value!=0)
# 33 of them
ICs <- foo %>% summarise_all(sd)%>%gather()%>%filter(value==0)
ICnames <- ICs$key
ICs <- foo%>%summarise_all(mean)%>%gather()%>%filter(key %in% ICnames)
ICs %<>% spread(key,value)

# Filter only to ones present in model
p_names %<>% filter(key %in% m_params)%>%pull(key)
# Left with 27

# parameters <- foo%>%group_by(VPOP)%>%slice(1)%>%filter(VPOP==910)%>%gather(key="Name",value="Values")%>%filter(Name %in% p_names)


# Get Parameter Limits
pLower <- foo %>% summarise_all(.funs = function(x){min(x)*0.85})
pUpper <- foo %>% summarise_all(.funs = function(x){max(x)*1.15})
pLower %<>% gather(key="Name",value="Lower")
pUpper %<>% gather(key="Name",value="Upper")
paramLims <- left_join(pLower,pUpper,by="Name")
paramLims %<>% filter(Name %in% p_names)

# Load model



# Set up simulation function
sim_fcn <- function(parameters=NULL,model,pnames,dosing,ICs,simulate=0){
  loadso(model)
  if(!is.null(parameters)){
    param_in <- data.frame(Names = pnames,Values=parameters)
    param_in <- spread(param_in,key=Names,value=Values)
  }else{
    param_in <- data.frame(ID=1)
  }
  param_in %<>% cbind(ICs)
  output <- model%>%idata_set(param_in) %>%Req(TUMOR)%>%obsonly%>%mrgsim(delta=56,end=56,events=as.ev(dosing))%>%
    filter(time==56)%>%as.data.frame()
  if(simulate==0){
    return(list(NSS=output))
  }else{
    return(output)
  }
}

stateLims <- data.frame(Name = 'TUMOR','Lower'=0.0,'Upper'=4.0,Time=56)
model_args = list(model=model,pnames=p_names,dosing=reg,ICs=ICs)
control=list(runParallel="parallel",nCores = 4,
                       parallel_libs="mrgsolve")
plausiblePatients <- generatePPs(model_fn = sim_fcn, NP=1e2, paramLims=paramLims,stateLims = stateLims,
                                 method='SA',
                                 model_args = model_args,
                                 scoreThreshold = 0)

hist_data <- plausiblePatients$simulation_results%>%mutate(Source="Plausible")
hist_data %<>% rbind(sims%>%select(ID,time,TUMOR)%>%mutate(Source="Simulation"))

ggplot(hist_data,aes(x=TUMOR))+geom_density(aes(x=TUMOR,y=..scaled..,fill=Source),alpha=0.5)

VPs <- getVPs(plausiblePatients, sims%>%select(ID,time,TUMOR),runs=20,plausible_pdf = 'auto',data_pdf='auto',
       alpha_algorithm = 'PSO')

hist_data <- plausiblePatients$simulation_results%>%mutate(Source="Plausible")
hist_data %<>% rbind(sims%>%select(ID,time,TUMOR)%>%mutate(Source="Simfaculation"))
hist_data %<>% rbind((VPs$VPs)%>%mutate(Source="VP"))
unique <- hist_data%>%group_by(Source)%>%count(Source)
hist_data <- merge(hist_data,y = unique, by='Source',all.x=TRUE)
ggplot(hist_data,aes(x=TUMOR, color=Source,fill=Source)) + geom_histogram(alpha=0.4, size = 1.5 )
ggplot(hist_data) + geom_density(aes(x=TUMOR,y=..scaled..,fill=Source),alpha=0.5)

