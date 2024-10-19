###### Setup ######

# packages
library(tidyverse)
library(lubridate)
library(lme4)
library(sandwich)
library(AER)
library(MASS, exclude = "select")
library(multcomp)
library(epitools)
library(sjPlot)

# core functions
results <- function(regression, contrast) {
	
	interaction <- glht(regression, linfct = contrast)
	output <- data.frame(((1-exp(cbind(ve=coef(interaction), ul=confint.default(interaction)[,1], ll=confint.default(interaction)[,2])))*100), 
											 var=((coef(interaction)-confint.default(interaction)[,1])/1.96)^2) %>% 
		rownames_to_column(var = "country")
	return(output)
} # to extract regression model results

weightingVE <- 
	function(VE, varEst, priorBeta, priorVar) {
		weightPrior <- 1/priorVar
		weightEst <- 1/varEst
		meanEst <- log(1-(VE/100))
		
		denom <- weightPrior + weightEst
		
		num <- (weightPrior * priorBeta) + (weightEst * meanEst)
		
		weightedVE <- num/denom
		
		return(weightedVE)
	} # for weighting VE by estimate precision

# Generate dummy trial dataset - one line per unique infant
dummy_data <- data.frame(
  trial = c(rep("102247",500), rep("102248",500), rep("107625",500), rep("113808",500), rep("444563/005",500), rep("444563/013",500), rep("444563/023",500), rep("444563/024",500), rep("444563/028",500), rep("444563/029",500), rep("444563/030",500)),
  country = c(sample(c("Czechia", "Finland", "France", "Spain"), 500, replace=T), sample(c("Malawi", "South Africa"), 500, replace=T), rep("Japan",500), rep("China",500), rep("United States",500), rep("South Africa",500), sample(c("Brazil", "Chile", "Colombia", "Dominican Republic", "Honduras", "Mexico", "Nicaragua", "Panama", "Peru", "Venezuela"), 500, replace=T), sample(c("Argentina", "Brazil", "Colombia", "Honduras", "Panama"), 500, replace=T), rep("Singapore",500), rep("Hong Kong",500), rep("Taiwan", 500)),
	pid = rep(c(1:500), 11), 
  arm = sample(c("Placebo", "HRV"), 500*11, replace=T), 
  sex = sample(c("M", "F"), 500*11, replace=T),
  birthDate = sample(seq.Date(mdy(10202003), mdy(10192004), by=1), 500*11, replace=T), 
  ageFirstDose = 8, 
  daysDose1 = floor(runif(500*11, 0, 6)),
  vaxInt = floor(runif(500*11, 20, 91)), 
  daysOneYear = floor(runif(500*11,308,315)), 
  daysLC = floor(runif(500*11,200,615)),
  rvge = 1,
  daysRvge = floor(runif(500*11, 21, 615)), 
  vesikari = floor(runif(500*11, 1, 21)), 
  opvCoadm = sample(c("both", "one", "none"), 500*11, replace=T)) %>% 
	mutate(
		doseDate1 = birthDate + (8*7) + daysDose1, 
		doseDate2 = doseDate1 + vaxInt,
		dateLC = doseDate2 + daysLC,
		oneYear = doseDate1 + daysOneYear,
		firstRvgeStart = if_else(rvge == 1, doseDate1 + daysRvge, NA_Date_), 
		rvge = if_else(firstRvgeStart >= doseDate1 & firstRvgeStart <= pmin(dateLC, oneYear), rvge, 0, 0),
		vesikari = if_else(rvge == 1, vesikari, NA_real_),
		firstRvgeStart = if_else(rvge == 1, firstRvgeStart, NA_Date_),
		endpoint = if_else(trial %in% c("113808", "102248", "102247", "444563/013", "444563/005"), "Any", "Severe"), 
		opvStatus = case_when(
			country %in% c("Argentina", "Brazil", "China", "Colombia", "Dominican Republic", "Honduras", "Malawi", "Panama", "South Africa") & arm == "HRV" ~ opvCoadm, 
			country %in% c("Argentina", "Brazil", "China", "Colombia", "Dominican Republic", "Honduras", "Malawi", "Panama", "South Africa") & arm == "HRV" ~ "none",
			TRUE ~ "no opv"))

all_trials <- 
	dummy_data %>% ###### REPLACE WITH DUMMY OBJECT 
	mutate(
		# did the first RVGE episode occur in the first year of life?
		rvge = if_else((firstRvgeStart > oneYear) | is.na(firstRvgeStart), 0, 1), 
		# date at which follow-up ends
		endFU = if_else(rvge == 0, pmin(oneYear, dateLC), firstRvgeStart),
		# person-time from first dose to end of follow-up
		timeFOI = time_length(endFU - doseDate1, unit="weeks"), 
		# person-time from second dose + 2 weeks to end of follow-up
		timeVE = time_length(endFU - (doseDate2+14), unit="weeks"),
		# whether or not the infant should be included in VE analysis by not having RVGE before second dose + 2 weeks
		incVE = case_when(
			is.na(doseDate2) ~ "N", 
			endFU <= (doseDate2+14) ~ "N", 
			as.numeric(vaxInt) >= 21 & as.numeric(vaxInt) <= 90 ~ "Y", 
			TRUE ~ "N"),
		# was the first RVGE severe?
		severe = if_else(rvge == 1 & vesikari >= 11, 1, 0), 
		# age at first RVGE episode
		rvgeAge = if_else(rvge == 1, time_length(firstRvgeStart - birthDate, unit="weeks"), NA_real_), 
		# age at first RVGE episode if it was severe
		severeAge = if_else(severe == 1, rvgeAge, NA_real_), 
		# dichotomized variable for concomitant oral administration
		opvCoadm = if_else(opvStatus %in% c("both", "one"), "Y", "N"))

# datasets with just placebo-assigned infants	(all trials and the any-severity subset of trials)
placebo_all <- all_trials %>% filter(arm == "Placebo")

placebo_any <- placebo_all %>% filter(endpoint == "Any") %>% droplevels() 

##### Proxy measures for the magnitude of infectious exposure #####

# median ages at hospitalization from Hasso-Agopsowicz et al
hasso <- data.frame(
	country=c("Argentina", "Brazil", "Chile", "China", "Colombia", "Czechia", "Dominican Republic", "Finland", "France", "Honduras", "Hong Kong", "Japan", "Malawi", "Mexico", "Nicaragua", "Panama", "Peru", "Singapore", "South Africa", "Spain", "Taiwan", "United States", "Venezuela"), 
	ageHosp = c(46, 46, 49, 52, 46, 65, 46, 65, 35, 43, 77, 73, 31, 46, 46, 46, 43, 65, 27, 65, 52, 65, 42))
# separate sets for the overall dataset and any-severity subset, starting with the overall dataset

proxy_countries <- 
	placebo_all %>% 
	group_by(country) %>% 
	summarize(denom = n(), 
						PT=sum(timeFOI), 
						nSev = sum(severe), 
						ageSev = median(severeAge, na.rm=T)) %>% 
	mutate(rateSev = (nSev/PT)*5200, 
				 rateSevSqrt = sqrt(rateSev), 
				 rateSevLog = log(rateSev), 
				 quartile = cut_number(rateSev, n=4, labels=F)) %>% 
	full_join(y=hasso, by="country") 

proxy_any_countries <- 
	placebo_any %>% 
	group_by(country) %>% 
	summarize(denom = n(), 
						PT=sum(timeFOI), 
						nAny = sum(rvge),
						nSev = sum(severe), 
						ageSev = median(severeAge, na.rm=T)) %>% 
	mutate(rateAny = (nAny/PT)*5200, rateSev = (nSev/PT)*5200, 
				 rateAnySqrt = sqrt(rateAny), rateSevSqrt = sqrt(rateSev), 
				 rateAnyLog = log(rateAny), rateSevLog = log(rateSev)) %>% 
	left_join(y=hasso, by="country")

##### Final datasets for analysis #####
pooled_all <- 
	all_trials %>% 
	filter(incVE == "Y") %>% 
	left_join(y=proxy_countries, by="country") %>% 
	droplevels() %>% 
	mutate(country=factor(country, levels = c("Singapore", "Argentina", "Brazil", "Chile", "China", "Colombia", "Czechia", "Dominican Republic", "Finland", "France", "Honduras", "Hong Kong", "Japan", "Malawi", "Mexico", "Nicaragua", "Panama", "Peru", "South Africa", "Spain", "Taiwan", "United States", "Venezuela")))

pooled_any <- 
	all_trials %>% 
	filter(incVE == "Y", endpoint == "Any") %>% 
	left_join(y=proxy_any_countries, by="country") %>% 
	droplevels() %>% 
	mutate(country=factor(country, levels = c("China", "Czechia", "Finland", "France", "Malawi", "South Africa", "Spain", "United States")))

##### Main Analysis | Standard Models #####
# empty dataframes to record results
ve_any <- data.frame(country=character(), ve=numeric(), ul=numeric(), ll=numeric(), var=numeric(), outcome=character(), approach=character())
ve_all <- data.frame(country=character(), ve=numeric(), ul=numeric(), ll=numeric(), var=numeric(), approach=character())


bic_all <- data.frame(approach=character(), bic=numeric())
bic_any <- data.frame(approach=character(), severity=character(), bic=numeric())

# all pooled data, VE against severe RVGE
model_contrast <- rbind(
	"Argentina"=c(0, 1, rep(0,22), rep(0,0),1,rep(0,21)), 
	"Brazil"=c(0, 1, rep(0,22), rep(0,1),1,rep(0,20)), 
	"Chile"=c(0, 1, rep(0,22), rep(0,2),1,rep(0,19)), 
	"China"=c(0, 1, rep(0,22), rep(0,3),1,rep(0,18)), 
	"Colombia"=c(0, 1, rep(0,22), rep(0,4),1,rep(0,17)), 
	"Czechia"=c(0, 1, rep(0,22), rep(0,5),1,rep(0,16)), 
	"Dominican Republic"=c(0, 1, rep(0,22), rep(0,6),1,rep(0,15)), 
	"Finland"=c(0, 1, rep(0,22), rep(0,7),1,rep(0,14)),
	"France"=c(0, 1, rep(0,22), rep(0,8),1,rep(0,13)),
	"Honduras"=c(0, 1, rep(0,22), rep(0,9),1,rep(0,12)),
	"Hong Kong"=c(0, 1, rep(0,22), rep(0,10),1,rep(0,11)),
	"Japan"=c(0, 1, rep(0,22), rep(0,11),1,rep(0,10)),
	"Malawi"=c(0, 1, rep(0,22), rep(0,12),1,rep(0,9)), 
	"Mexico"=c(0, 1, rep(0,22), rep(0,13),1,rep(0,8)), 
	"Nicaragua"=c(0, 1, rep(0,22), rep(0,14),1,rep(0,7)),
	"Panama"=c(0, 1, rep(0,22), rep(0,15),1,rep(0,6)),
	"Peru"=c(0, 1, rep(0,22), rep(0,16),1,rep(0,5)),
	"Singapore"=c(0, 1, rep(0,22), rep(0,17),0,rep(0,4)), # this is ref group
	"South Africa"=c(0, 1, rep(0,22), rep(0,17),1,rep(0,4)),
	"Spain"=c(0, 1, rep(0,22), rep(0,18),1,rep(0,3)),
	"Taiwan"=c(0, 1, rep(0,22), rep(0,19),1,rep(0,2)),
	"United States"=c(0, 1, rep(0,22), rep(0,20),1,rep(0,1)),
	"Venezuela"=c(0, 1, rep(0,22), rep(0,21),1,rep(0,0)))

model <- glm(severe ~ arm + country + arm*country + offset(log(timeVE)), 
						 data = pooled_all, family = "poisson")
bic_all <- bic_all %>% add_row(approach="standard", bic=BIC(model))
ve_all <- ve_all %>% 
	bind_rows(results(model, model_contrast) %>% mutate(approach = "standard"))

# any-severity subset, VE against any RVGE
model_contrast <- rbind(
	"China" =          c(0, 1, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0), 
	"Czechia" =        c(0, 1, 0,0,0,0,0,0,0, 1,0,0,0,0,0,0), 
	"Finland" =        c(0, 1, 0,0,0,0,0,0,0, 0,1,0,0,0,0,0), 
	"France" =         c(0, 1, 0,0,0,0,0,0,0, 0,0,1,0,0,0,0), 
	"Malawi" =         c(0, 1, 0,0,0,0,0,0,0, 0,0,0,1,0,0,0), 
	"South Africa" =   c(0, 1, 0,0,0,0,0,0,0, 0,0,0,0,1,0,0), 
	"Spain" =          c(0, 1, 0,0,0,0,0,0,0, 0,0,0,0,0,1,0), 
	"United States" =  c(0, 1, 0,0,0,0,0,0,0, 0,0,0,0,0,0,1))

model <- glm(rvge ~ arm + country + arm*country + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="standard", severity="any", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "any", approach = "standard"))

# any-severity subset, VE against severe RVGE - reuse prior contrasts
model <- glm(severe ~ arm + country + arm*country + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="standard", severity="severe", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "severe", approach = "standard"))

##### Main Analysis | Any-Severity Rate Models #####
model_contrast <- 
	proxy_any_countries %>% 
	select(country, rateAny) %>% 
	mutate(int=0, vax=1, proxy=0) %>% 
	column_to_rownames("country") %>% 
	select(int, vax, proxy, rateAny) %>% # need to be in right order for contrast
	as.matrix()

# any-severity subset, VE against any RVGE
model <- glm(rvge ~ arm + rateAny + arm*rateAny + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="rateAny", severity="any", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "any", approach = "rateAny"))

# any-severity subset, VE against severe RVGE - reuse prior contrasts
model <- glm(severe ~ arm + rateAny + arm*rateAny + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="rateAny", severity="severe", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "severe", approach = "rateAny"))

##### Main Analysis | Severe Rate Models #####

# all pooled data, VE against severe RVGE
model_contrast <- 
	proxy_countries %>% 
	select(country, rateSev) %>% 
	mutate(int=0, vax=1, proxy=0) %>% 
	column_to_rownames("country") %>% 
	select(int, vax, proxy, rateSev) %>%
	as.matrix()

model <- glm(severe ~ arm + rateSev + arm*rateSev + offset(log(timeVE)), 
						 data = pooled_all, family = "poisson")
bic_all <- bic_all %>% add_row(approach="rateSev", bic=BIC(model))
ve_all <- ve_all %>% 
	bind_rows(results(model, model_contrast) %>% mutate(approach = "rateSev"))

# any-severity subset, VE against any RVGE
model_contrast <- 
	proxy_any_countries %>% 
	select(country, rateSev) %>% 
	mutate(int=0, vax=1, proxy=0) %>% 
	column_to_rownames("country") %>% 
	select(int, vax, proxy, rateSev) %>% # need to be in right order for contrast
	as.matrix()

model <- glm(rvge ~ arm + rateSev + arm*rateSev + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="rateSev", severity="any", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "any", approach = "rateSev"))

# any-severity subset, VE against severe RVGE - reuse prior contrasts
model <- glm(severe ~ arm + rateSev + arm*rateSev + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="rateSev", severity="severe", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "severe", approach = "rateSev"))

##### Main Analysis | Severe Age Models #####

# all pooled data, VE against severe RVGE
model_contrast <- 
	proxy_countries %>% 
	select(country, ageSev) %>% 
	mutate(int=0, vax=1, proxy=0) %>% 
	column_to_rownames("country") %>% 
	select(int, vax, proxy, ageSev) %>%
	as.matrix()

model <- glm(severe ~ arm + ageSev + arm*ageSev + offset(log(timeVE)), 
						 data = pooled_all, family = "poisson")
bic_all <- bic_all %>% add_row(approach="ageSev", bic=BIC(model))
ve_all <- ve_all %>% 
	bind_rows(results(model, model_contrast) %>% mutate(approach = "ageSev"))

# any-severity subset, VE against any RVGE
model_contrast <- 
	proxy_any_countries %>% 
	select(country, ageSev) %>% 
	mutate(int=0, vax=1, proxy=0) %>% 
	column_to_rownames("country") %>% 
	select(int, vax, proxy, ageSev) %>% # need to be in right order for contrast
	as.matrix()

model <- glm(rvge ~ arm + ageSev + arm*ageSev + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="ageSev", severity="any", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "any", approach = "ageSev"))

# any-severity subset, VE against severe RVGE - reuse prior contrasts
model <- glm(severe ~ arm + ageSev + arm*ageSev + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="ageSev", severity="severe", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "severe", approach = "ageSev"))

##### Main Analysis | Hospitalization Age Models #####

# all pooled data, VE against severe RVGE
model_contrast <- 
	proxy_countries %>% 
	select(country, ageHosp) %>% 
	mutate(int=0, vax=1, proxy=0) %>% 
	column_to_rownames("country") %>% 
	select(int, vax, proxy, ageHosp) %>%
	as.matrix()

model <- glm(severe ~ arm + ageHosp + arm*ageHosp + offset(log(timeVE)), 
						 data = pooled_all, family = "poisson")
bic_all <- bic_all %>% add_row(approach="ageHosp", bic=BIC(model))
ve_all <- ve_all %>% 
	bind_rows(results(model, model_contrast) %>% mutate(approach = "ageHosp"))

# any-severity subset, VE against any RVGE
model_contrast <- 
	proxy_any_countries %>% 
	select(country, ageHosp) %>% 
	mutate(int=0, vax=1, proxy=0) %>% 
	column_to_rownames("country") %>% 
	select(int, vax, proxy, ageHosp) %>% # need to be in right order for contrast
	as.matrix()

model <- glm(rvge ~ arm + ageHosp + arm*ageHosp + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="ageHosp", severity="any", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "any", approach = "ageHosp"))

# any-severity subset, VE against severe RVGE - reuse prior contrasts
model <- glm(severe ~ arm + ageHosp + arm*ageHosp + offset(log(timeVE)), 
						 data = pooled_any, family = "poisson")
bic_any <- bic_any %>% add_row(approach="ageHosp", severity="severe", bic=BIC(model))
ve_any <- ve_any %>% 
	bind_rows(results(model, model_contrast) %>% 
							mutate(outcome = "severe", approach = "ageHosp"))

##### Sensitivity Analysis | Different specifications of severe rate #####
# using log-transformed rate of severe RVGE
model_contrast <- 
	proxy_countries %>% 
	select(country, rateSevLog) %>% 
	mutate(int=0, vax=1, proxy=0) %>% 
	column_to_rownames("country") %>% 
	select(int, vax, proxy, rateSevLog) %>%
	as.matrix()

model <- glm(severe ~ arm + rateSevLog + arm*rateSevLog + offset(log(timeVE)), 
						 data = pooled_all, family = "poisson")
bic_all <- bic_all %>% add_row(approach="rateSevLog", bic=BIC(model))
ve_all <- ve_all %>% 
	bind_rows(results(model, model_contrast) %>% mutate(approach = "rateSevLog"))

# using quartiles of severe RVGE
model_contrast <- rbind(
	"Q1"=c(0, 1, 0,0,0, 0,0,0), 
	"Q2"=c(0, 1, 0,0,0, 1,0,0), 
	"Q3"=c(0, 1, 0,0,0, 0,1,0), 
	"Q4"=c(0, 1, 0,0,0, 0,0,1))

model <- glm(severe ~ arm + factor(quartile) + arm*factor(quartile) + offset(log(timeVE)), 
						 data = pooled_all, family = "poisson")
bic_all <- bic_all %>% add_row(approach="quartile", bic=BIC(model))
ve_all <- ve_all %>% 
	bind_rows(results(model, model_contrast) %>% mutate(approach = "quartile"))

##### Sensitivity Analysis | Restricting to Countries with Country-Specific Median Age at Hospitalization #####
include_me <- c("France", "Hong Kong", "Japan", "Taiwan", "Chile", "China", "Venezuela", "South Africa", "Malawi")
model_contrast <- 
  proxy_countries %>% 
  filter(country %in% include_me) %>%
  select(country, ageHosp) %>% 
  mutate(int=0, vax=1, proxy=0) %>% 
  column_to_rownames("country") %>% 
  select(int, vax, proxy, ageHosp) %>%
  as.matrix()

model <- glm(severe ~ arm + ageHosp + arm*ageHosp + offset(log(timeVE)), 
             data = pooled_all %>% filter(country %in% include_me), family = "poisson")
bic_all <- bic_all %>% add_row(approach="ageHospSpec", bic=BIC(model))
ve_all <- ve_all %>% 
  bind_rows(results(model, model_contrast) %>% mutate(approach = "ageHospSpec"))

##### Sensitivity Analysis | Stratified Analysis by Country #####
# all pooled data
for (i in 1:nrow(proxy_countries)) {
	
	whichCountry <- proxy_countries$country[i]
	
	model <- glm(severe ~ arm + offset(log(timeVE)), 
							 data = pooled_all %>% filter(country == whichCountry), 
							 family = "poisson")
	ve_all <- ve_all %>% 
		bind_rows(results(model, matrix(c(0,1), nrow=1)) %>% mutate(country=whichCountry, approach="stratified"))
}

# any-severity subset
for (i in 1:nrow(proxy_any_countries)) {
	
	whichCountry <- proxy_any_countries$country[i]
	whichApproach <- paste0("Stratified-", whichCountry)
	
	model <- glm(rvge ~ arm + offset(log(timeVE)), 
							 data = pooled_any %>% filter(country == whichCountry), 
							 family = "poisson")
	
	ve_any <- ve_any %>% 
		bind_rows(results(model, matrix(c(0,1), nrow=1)) %>% mutate(country=whichCountry, outcome="any", approach="stratified"))
	
	model <- glm(severe ~ arm + offset(log(timeVE)), 
							 data = pooled_any %>% filter(country == whichCountry), 
							 family = "poisson")
	ve_any <- ve_any %>% 
		bind_rows(results(model, matrix(c(0,1), nrow=1)) %>% mutate(country=whichCountry, outcome="severe", approach="stratified"))
}

##### Sensitivity Analysis | Approximate Bayesian #####
# priors 
priorBetaAny <- log(0.5) # weakly informative coefficient
priorVarAny <- ((log(1.3)-log(0.5))/1.96)^2
priorBetaSev <- log(0.4) # weakly informative coefficient
priorVarSev <- ((log(1.3)-log(0.4))/1.96)^2

# check what VE and 95% CI the priors represent
1-exp(priorBetaAny) # 50%
1-exp(priorBetaAny + (1.96*sqrt(priorVarAny))) # -30% 
1-exp(priorBetaAny - (1.96*sqrt(priorVarAny))) # 80%

1-exp(priorBetaSev) # 60%
1-exp(priorBetaSev + (1.96*sqrt(priorVarSev))) # -30% 
1-exp(priorBetaSev - (1.96*sqrt(priorVarSev))) # 88%

# calculating weighted VEs
ve_all <- 
	ve_all %>% 
	mutate(veWeighted = (1-exp(weightingVE(ve, var, priorBetaSev, priorVarSev)))*100)

ve_any <- 
	ve_any %>% 
	mutate(veWeighted = if_else(outcome=="any", 
															(1-exp(weightingVE(ve, var, priorBetaAny, priorVarAny)))*100,
															(1-exp(weightingVE(ve, var, priorBetaSev, priorVarSev)))*100))

##### Main Results #####

nrow(pooled_all) # number of infants included in main severe VE analysis
length(unique(pooled_all$country)) # number of unique countries
length(unique(pooled_all$trial)) # number of unique trials
sum(pooled_all$rvge, na.rm=T) # number of RVGE episodes
sum(pooled_all$severe, na.rm=T) # number of severe RVGE episodes
prop.table(table(pooled_all$country))*100 # which countries are biggest contributors
nrow(pooled_all %>% filter(country %in% c("Mexico", "Peru"))) # number of infants in Mexico and Peru, the largest contributors in the real dataset
prop.table(table(pooled_all$sex))*100 # pretty evenly split, as expected
nrow(pooled_all %>% filter(sex=="F")) # number of female infants
summary(pooled_all$ageFirstDose) # median and iqr of age at first dose in weeks
summary(pooled_all$timeVE) # median and iqr of follow-up time for infants included in VE analyses

# Characteristics of infants included in analyses of VE against severe RVGE
table(pooled_all$arm) # denominators for percentage calcs in table 1
table1 <- bind_rows(
	pooled_all %>% group_by(arm, sex) %>% summarize(n=n()) %>% pivot_wider(id_cols="sex", names_from="arm", values_from="n") %>% rename("characteristic"="sex"), 
	pooled_all %>% group_by(arm, as.character(country)) %>% summarize(n=n()) %>% pivot_wider(id_cols="as.character(country)", names_from="arm", values_from="n") %>% rename("characteristic"="as.character(country)"), 
	pooled_all %>% mutate(outcome_sev = case_when(severe==1~"rvge2", rvge==1~"rvge1", TRUE~"rvge0")) %>% group_by(arm, outcome_sev) %>% summarize(n=n()) %>% pivot_wider(id_cols="outcome_sev", names_from="arm", values_from="n") %>% rename("characteristic"="outcome_sev")) %>% 
	mutate(pct0 = signif(Placebo/39064*100, digits=2), pct1 = signif(HRV/44528*100), digits=2) %>% # note that the original counts are used here, not the counts present in the dummy data
	select(characteristic, Placebo, pct0, HRV, pct1)
tapply(pooled_all$ageFirstDose, pooled_all$arm, summary) # median and iqr of age at first dose in weeks
tapply(pooled_all$timeVE, pooled_all$arm, summary) # median and iqr of follow-up time by trial arm

# now repeating the above for the any-severity trial data
nrow(pooled_any) 
nrow(pooled_any)/nrow(pooled_all) 
length(unique(pooled_any$country)) 
sum(pooled_any$rvge, na.rm=T) 
sum(pooled_any$severe, na.rm=T) 
prop.table(table(pooled_any$country))*100 
nrow(pooled_any %>% filter(country %in% c("China", "Finland", "South Africa"))) 
prop.table(table(pooled_any$sex))*100 
nrow(pooled_any %>% filter(sex=="F")) 
summary(pooled_any$ageFirstDose) 
summary(pooled_any$timeVE) 

# Characteristics of infants included in analyses of VE against any-severity RVGE
nrow(pooled_any %>% filter(arm=="HRV"))/nrow(pooled_any)
table(pooled_any$arm)
table2 <- bind_rows(
	pooled_any %>% group_by(arm, sex) %>% summarize(n=n()) %>% pivot_wider(id_cols="sex", names_from="arm", values_from="n") %>% rename("characteristic"="sex"), 
	pooled_any %>% group_by(arm, as.character(country)) %>% summarize(n=n()) %>% pivot_wider(id_cols="as.character(country)", names_from="arm", values_from="n") %>% rename("characteristic"="as.character(country)"), 
	pooled_any %>% mutate(outcome_sev = case_when(severe==1~"rvge2", rvge==1~"rvge1", TRUE~"rvge0")) %>% group_by(arm, outcome_sev) %>% summarize(n=n()) %>% pivot_wider(id_cols="outcome_sev", names_from="arm", values_from="n") %>% rename("characteristic"="outcome_sev")) %>% 
	mutate(pct0 = signif(Placebo/4468*100, digits=2), pct1 = signif(HRV/7481*100), digits=2) %>% 
	select(characteristic, Placebo, pct0, HRV, pct1)
tapply(pooled_any$ageFirstDose, pooled_any$arm, summary)
tapply(pooled_any$timeVE, pooled_any$arm, summary)

# VE against severe RVGE - main pooled dataset
table3 <- 
	ve_all %>% 
	filter(approach %in% c("standard", "rateSev", "ageSev", "ageHosp")) %>% 
	mutate(ve = sprintf("%1.0f", ve), 
				 ci = paste0("(", sprintf("%1.0f", ll), ",", sprintf("%1.0f", ul), ")")) %>% 
	select(country, approach, ve, ci) %>% 
	pivot_wider(names_from="approach", values_from=c("ve","ci"), names_vary="slowest")

bic_all %>% filter(approach %in% c("standard", "rateSev", "ageSev", "ageHosp"))

# VE against any-severity RVGE - any-severity trial subset
table4 <- 
	ve_any %>% 
	filter(approach %in% c("standard", "rateAny", "rateSev", "ageSev", "ageHosp"), outcome=="any") %>% 
	mutate(ve = sprintf("%1.0f", ve), 
				 ci = paste0("(", sprintf("%1.0f", ll), ",", sprintf("%1.0f", ul), ")")) %>% 
	select(country, approach, ve, ci) %>% 
	pivot_wider(names_from="approach", values_from=c("ve","ci"), names_vary="slowest")

bic_any %>% filter(approach %in% c("standard", "rateAny", "rateSev", "ageSev", "ageHosp"), severity=="any")

##### Supplemental Info #####
# proxy metrics for infectious exposures
st2 <- full_join(
	proxy_countries %>% select(id=country, denom, nSev, rateSev, ageSev, ageHosp), 
	proxy_any_countries %>% select(id=country, nAny, rateAny), 
	by="id"
) %>% bind_rows(
	full_join(
		placebo_all %>% filter(country %in% c("Argentina", "Brazil", "Colombia", "Honduras", "Panama", "South Africa")) %>% mutate(id = paste(country, trial, sep="-")) %>% group_by(id) %>% summarize(denom=n(), PT=sum(timeFOI), nSev=sum(severe), ageSev=median(severeAge, na.rm=T)) %>% mutate(rateSev=(nSev/PT)*5200) %>% select(id, denom, nSev, rateSev, ageSev), 
		placebo_any %>% filter(country=="South Africa") %>% mutate(id = paste(country, trial, sep="-")) %>% group_by(id) %>% summarize(denom=n(), PT=sum(timeFOI), nAny=sum(rvge)) %>% mutate(rateAny=(nAny/PT)*5200) %>% select(id, nAny, rateAny), 
		by="id")) %>% 
	select(id, denom, nAny, nSev, rateAny, rateSev, ageSev, ageHosp) %>%
	arrange(id)

# comparing model performance between untransformed and log-transformed severe rate
st3 <- 
	ve_all %>% 
	filter(approach %in% c("rateSev", "rateSevLog")) %>% 
	mutate(ve = sprintf("%1.0f", ve), 
				 ci = paste0("(", sprintf("%1.0f", ll), ",", sprintf("%1.0f", ul), ")")) %>% 
	select(country, approach, ve, ci) %>% 
	pivot_wider(names_from="approach", values_from=c("ve","ci"), names_vary="slowest")

# model using quartile of severe rate
st4 <- 
	ve_all %>% 
	filter(approach == "quartile") %>% 
	mutate(ve = sprintf("%1.0f", ve), 
				 ci = paste0("(", sprintf("%1.0f", ll), ",", sprintf("%1.0f", ul), ")")) %>% 
	select(country, approach, ve, ci)

# VE against severe RVGE using only countries with country-specific median age at hospitalization
st5 <- 
  ve_all %>% 
  filter(approach == "ageHospSpec") %>% 
  mutate(ve = sprintf("%1.0f", ve), 
         ci = paste0("(", sprintf("%1.0f", ll), ",", sprintf("%1.0f", ul), ")")) %>% 
  select(country, approach, ve, ci)

# VE against severe RVGE when using the any-severity data subset
st6 <- 
	ve_any %>% 
	filter(approach %in% c("standard", "rateAny", "rateSev", "ageSev", "ageHosp"), outcome=="severe") %>% 
	mutate(ve = sprintf("%1.0f", ve), 
				 ci = paste0("(", sprintf("%1.0f", ll), ",", sprintf("%1.0f", ul), ")")) %>% 
	select(country, approach, ve, ci) %>% 
	pivot_wider(names_from="approach", values_from=c("ve","ci"), names_vary="slowest")

bic_any %>% filter(approach %in% c("standard", "rateAny", "rateSev", "ageSev", "ageHosp"), severity=="severe")

# VE from stratified models
st7 <- full_join(
	ve_any %>% filter(approach=="stratified") %>% mutate(ve = sprintf("%1.0f", ve), ci = paste0("(", sprintf("%1.0f", ll), ",", sprintf("%1.0f", ul), ")")) %>% select(country, outcome, ve, ci) %>% pivot_wider(id_cols="country", names_from="outcome", values_from=c("ve","ci"), names_vary="slowest"), 
	ve_all %>% filter(approach=="stratified") %>% mutate(ve = sprintf("%1.0f", ve), ci = paste0("(", sprintf("%1.0f", ll), ",", sprintf("%1.0f", ul), ")")) %>% select(country, ve, ci), 
	by="country") %>% 
	arrange(country)

# Comparing weighted and unweighted VE estimates for VE against any-severity RVGE
st8 <- 
	ve_any %>% 
	filter(approach %in% c("standard", "rateAny", "rateSev", "ageSev", "ageHosp"), outcome=="any") %>% 
	mutate(change = veWeighted - ve) %>% 
	select(country, approach, veWeighted, change) %>% 
	mutate(veWeighted = sprintf("%1.0f", veWeighted), 
				 change = sprintf("%1.1f", change)) %>% 
	pivot_wider(names_from = "approach", values_from = c("veWeighted", "change"), names_vary="slowest")

# now for VE against severe RVGE
st9 <- 
	ve_all %>% 
	filter(approach %in% c("standard", "rateSev", "ageSev", "ageHosp")) %>% 
	mutate(change = veWeighted - ve) %>% 
	select(country, approach, veWeighted, change) %>% 
	mutate(veWeighted = sprintf("%1.0f", veWeighted), 
				 change = sprintf("%1.1f", change)) %>% 
	pivot_wider(names_from = "approach", values_from = c("veWeighted", "change"), names_vary="slowest")
