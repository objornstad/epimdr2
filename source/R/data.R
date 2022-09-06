#' UK measles CCS data.
#'
#'The fraction of weeks measles was absent from each of the 954 cities
#'and towns of England and Wales between 1944 and 1965.
#'
#' @format A data frame with 954 rows and 14 variables:
#' \describe{
#' \item{fade3}{Average duration of fadeout (of at least 3 weeks of length)}
#'   \item{ext}{Fraction of time when measles was absent}
#'   \item{size}{Median population size}
#'   \item{fade}{Average duration of fadeouts (of a week or longer)}
#'   \item{se3}{Standard error fade3}
#'   \item{se}{Standard error of fade}
#'   \item{n3}{The number of fadeouts (of at least 3 weeks of length)} 
#'   \item{n}{The number of fadeout of a week or longer}
#'   \item{names}{City/town name}
#' }
#' @source Bjornstad and Grenfell (2008) Hazards, spatial transmission and 
#' timing of outbreaks in epidemic metapopulations. 
#' Environmental and Ecological Statistics 15: 265-277. <doi:10.1007/s10651-007-0059-3>
"ccs"

#' Massachusetts gonorrhea data.
#'
#'Weekly cases of gonorrhea in Massachusetts between 2006 and 2015.
#'
#' @format A data frame with 422 rows and 4 variables:
#' \describe{
#' \item{number}{Weekly case reports}
#'   \item{year}{Year}
#'   \item{week}{Week of the year}
#'   \item{time}{Time in fractions of year}
#' }
#' @source \url{https://www.tycho.pitt.edu}
"magono"

#' Dacca cholera death data.
#'
#'Monthly deaths from  cholera in Dacca, East Bengal between 1891 and 1940.
#'
#' @format A data frame with 600 rows and 4 variables:
#' \describe{
#'   \item{Year}{Year}
#'   \item{Month}{Month of the year}
#' \item{Dacca}{Monthly cholera deaths}
#'   \item{Population}{Population size of district}
#' }
#' @source King, A.A., Ionides, E.L., Pascual, M. and Bouma, M. J. (2008) 
#' Inapparent infections and cholera dynamics. Nature, 454:877-880. <doi:10.1038/nature07084>
"cholera"


#' Weekly measles incidence from 2003-04 in Niamey, Niger.
#'
#' A dataset containing the weekly incidence of measles in
#' Niamey, Niger during the 2003-04 outbreak
#'
#' @format A data frame with 31 rows and 13 variables:
#' \describe{
#'   \item{absweek}{week since beginning of outbreak}
#'   \item{week}{week of the year}
#'   \item{tot_cases}{weekly incidence for the whole city}
#'   \item{tot_mort}{weekly deaths for the whole city}
#'   \item{lethality}{weekly case fatality rate}
#'   \item{tot_attack}{weekly attack rates for the whole city}
#'   \item{cases_1}{weekly incidence for district 1}
#'   \item{attack_1}{weekly attack rates for district 1}
#'   \item{cases_2}{weekly incidence for district 2}
#'   \item{attack_2}{weekly attack rates for district 2}
#'   \item{cases_3}{weekly incidence for district 3}
#'   \item{attack_3}{weekly attack rates for district 3}
#'   \item{cum_cases}{weekly cumulative incidence for the whole city}
#' }
#' @source Grais et al (2008) Time is of the essence: exploring a measles outbreak response vaccination in Niamey, Niger. Journal of the Royal Society Interface 5: 67-74. <doi:10.1098/rsif.2007.1038>
"niamey"

#' Boarding school influenza data.
#'
#' The daily number of children confined to bed in a boarding school in North England during an outbreak in 1978 of the reemerging A/H1N1 strain.
#' The school had 763 boys of which 512 boys were confined to bed sometime during the outbreak.
#'
#' @format A data frame with 14 rows and 2 variables:
#' \describe{
#'   \item{day}{day since beginning of outbreak}
#'   \item{cases}{number of sick children}
#' }
#' @source Anonymous (1978) EPIDEMIOLOGY: Influenza in a boarding school. British Medical Journal, 4 March 1978 p.587.  
"flu"

#' Sierra-Leone Ebola 2015 data.
#'
#' The daily number of cases of ebola in Sierra Leone during the 2015 epidemic.
#'
#' @format A data frame with 103 rows and 4 variables:
#' \describe{
#'   \item{date}{date}
#'   \item{day}{day}
#'   \item{cum_cases}{cumulative incidence}
#'   \item{cases}{incidence calculated by differencing the cumcases and setting negatives to zero.}
#' }
#' @source \url{https://www.cdc.gov/vhf/ebola/outbreaks/2014-west-africa/cumulative-cases-graphs.html}
"ebola"

#' Ferrari et al. 2005 outbreak data.
#'
#' The incidence aggregated by serial interval of a number of outbreaks studied by Ferrari et al. 2005.
#'
#' @format A data frame with 15 rows and 7 variables:
#' \describe{
#'   \item{Eboladeaths00}{Number of deaths from ebola during the 2000 Uganda outbreak}
#'   \item{Ebolacases00}{Number of cases of ebola during the 2000 Uganda outbreak}
#'   \item{Ebolacases95}{Number of cases of ebola during the 1995 DRC outbreak}
#'   \item{FMDfarms}{Number of farms infected with FMD during the 2000-01 UK outbreak}
#'   \item{HogCholera}{Number of cases of swine fever in pigs in the 1997-98 outbreak in the Netherlands}
#'   \item{SarsHk}{Number of cases of SARS in Hong Kong during the 2003 outbreak}
#'   \item{SarsSing}{Number of cases of SARS in Singapore during the 2003 outbreak}
#' }
#' @source Ferrari et al. (2005) Estimation and inference of R-0 of an infectious pathogen by a removal method. Mathematical Biosciences 198: 14-26. <doi:10.1016/j.mbs.2005.08.002>
"ferrari"

#' US SARS-CoV-2 variant data.
#'
#'Weekly fraction of identification of the various CoV-2 variants May 2021 through March 2022.
#'
#' @format A data frame with 47 rows and 7 variables:
#' \describe{
#' \item{date}{End of week of sample}
#'   \item{other}{Early variants}
#'   \item{B.1.617.2}{Delta variant}
#'   \item{B.1.1.529}{First omicron variant}
#'   \item{BA.1.1}{Omicron variant BA.1}
#'   \item{BA.2}{Omicron variant BA.2 and BA.2.12}
#'   \item{BA.2.12.1}{Omicron variant BA 2.12.1}
#' }
#' @source \url{https://coronavirus.health.ny.gov/covid-19-variant-data}
"variants"


#' De et al. 2004 gonorrhea contact matrix
#'
#' The directed contact network from De et al. (2004) contact-tracing of the spread of gonorrhea across asexual network in Alberta, Canada
#'
#' @format A matrix with 89 rows and 89 columns:
#' \describe{
#'   \item{gonet}{a matrix of directional contacts of disease spread}
#' }
#' @source De et al (2004). Sexual network analysis of a gonorrhea outbreak. Sexually transmitted infections 80: 280-285. <doi:10.1136/sti.2003.007187>
"gonnet"

#' Black's measles seroprevalence data.
#'
#' Seroprevalence-by-age-bracket for measles in prevaccination New Haven as studied by Black (1959).
#'
#' @format A data frame with 42 rows and 3 variables:
#' \describe{
#'   \item{age}{ age-bracket (in years)}
#'   \item{mid}{mid-point of age-bracket (in years)}
#'   \item{n}{number of tests}
#'   \item{pos}{number seropositive}
#'   \item{neg}{number seronegative}
#'   \item{f}{seroprevalence}
#' }
#' @source Black (1959) Measles antibodies in the population of New Haven, Connecticut. Journal of Immunology 83:74-83
"black"


#' Rabbit \emph{Bordetella bronchiseptica} data.
#'
#'Rabbits infected by \emph{B. bronchiseptica} by age as studied by Long et al (2010).
#'
#' @format A data frame with 42 rows and 3 variables:
#' \describe{
#'   \item{a}{end of age-bracket (in months)}
#'   \item{n}{number of rabbits tested}
#'   \item{inf}{number of rabbits infected with the bacterium}
#' }
#' @source Long et al (2010) Identifying the Age Cohort Responsible for Transmission in a Natural Outbreak of Bordetella bronchiseptica. PLoS Pathogens 6(12): e1001224. <doi:10.1371/journal.ppat.1001224>
"rabbit"

#' Rubella in Peru data.
#'
#' Rubella incidence by age as studied by Metcalf et al (2011).
#'
#' @format A data frame with 95 rows and 2 variables:
#' \describe{
#'   \item{age}{end of age-bracket (in years)}
#'   \item{cumulative}{cumulative number of rubella cases}
#'   \item{incidence}{number of rubella cases}
#'   \item{n}{total cases}
#' }
#' @source Metcalf et al (2011) Rubella metapopulation dynamics and importance of spatial coupling to the risk of congenital rubella syndrome in Peru. Journal of the Royal Society Interface 8: 369-376. <doi:10.1371/journal.pone.0072086>
"peru"

#' POLYMOD contact-rate data by Age.
#'
#' Age-specific contact rates from the diary study by Mossong et al. 2008.
#'
#' @format A data frame with 900 rows and 3 variables:
#' \describe{
#'   \item{contactor}{end of age-bracket (in years) of contactor group}
#'   \item{contactee}{end of age-bracket (in years) of contactee group}
#'   \item{contact.rate}{average contact rate}
#' }
#' @source Mossong et al. 2008 Social contacts and mixing patterns relevant to the spread of infectious diseases PLoS Med, Public Library of Science  5:e74. <doi:10.1371/journal.pmed.0050074>.
"polymod"

#' 2005 US Life table.
#'
#' Survivorship and fecundities for the US in 2005 by 5 year age-brackets.
#'
#' @format A data frame with 20 rows and 4 variables:
#' \describe{
#'   \item{a}{end of age-bracket (in years)}
#'   \item{la}{fraction of birth cohort still alive}
#'   \item{fa}{fecundity at age}
#'   \item{sa}{survival probabilities per age-bracket}
#' }
"us"


#' Weekly deaths from Influenza-like illness in Pennsylvania between 1972 and 1998.
#'
#' A dataset containing the weekly ILI 
#' related deaths in  Pennsylvania between 1972 and 1998.
#'
#' @format A data frame with 1404 rows and 3 variables:
#' \describe{
#'   \item{PENNSYLVANIA}{weekly deaths}
#'   \item{YEAR}{the year}
#'   \item{WEEK}{the week}
#' }
#' @source \url{https://www.tycho.pitt.edu}
"paili"

#' Weekly incidence of Lymes disease in Pennsylvania between 2006 and 2014.
#'
#' A dataset containing the weekly incidence of Lymes disease
#'  in  Pennsylvania between 2006 and 2014.
#'
#' @format A data frame with 448 rows and 3 variables:
#' \describe{
#'   \item{PENNSYLVANIA}{weekly incidence}
#'   \item{YEAR}{the year}
#'   \item{WEEK}{the week}
#' }
#' @source \url{https://www.tycho.pitt.edu}
"palymes"

#' Weekly incidence of giardia in Pennsylvania between 2006 and 2014.
#'
#' A dataset containing the weekly incidence of giardia
#'  in  Pennsylvania between 2006 and 2014.
#'
#' @format A data frame with 448 rows and 3 variables:
#' \describe{
#'   \item{PENNSYLVANIA}{weekly incidence}
#'   \item{YEAR}{the year}
#'   \item{WEEK}{the week}
#' }
#' @source \url{https://www.tycho.pitt.edu}
"pagiard"

#' Weekly incidence of measles in Pennsylvania between 1928 and 1969.
#'
#' A dataset containing the weekly incidence of measles
#'  in  Pennsylvania between 2006 and 2014.
#'
#' @format A data frame with 448 rows and 3 variables:
#' \describe{
#'   \item{PENNSYLVANIA}{weekly incidence}
#'   \item{YEAR}{the year}
#'   \item{WEEK}{the week}
#' }
#' @source \url{https://www.tycho.pitt.edu}
"pameasle"


#' Bi-weekly measles incidence in London from 1944-65.
#'
#' A dataset containing the biweekly incidence of measles in
#' London from 1944 to 1965
#'
#' @format A data frame with 546 rows and 5 variables:
#' \describe{
#'   \item{year}{year}
#'   \item{week}{week of the year}
#'   \item{time}{time}
#'   \item{London}{incidence}
#'   \item{B}{Biweekly births}
#' }
#' @details Birth numbers are annual, so in the data set, this number is evenly distributed across the 26 bi-weeks of each year.
#' @source Bjornstad et al. (2002) Endemic and epidemic dynamics of measles: Estimating transmission rates and their scaling using a time series SIR model. Ecological Monographs 72: 169-184. <doi:10.2307/3100023>
"meas"

#c6
#' Daily measures of malaria infected mice.
#'
#' Daily data on laboratory mice infected with various strains of \emph{Plasmodium chabaudi} 
#'
#' @format A data frame with 1300 rows and 11 variables:
#' \describe{
#'   \item{Line}{line number}
#'   \item{Day}{day of infection}
#'   \item{Box}{Cage number}
#'   \item{Mouse}{Mouse identifier}
#'   \item{Treatment}{Plasmodium strain}
#'   \item{Ind2}{Unique mouse identifier}
#'   \item{Weight}{Mouse weight}
#'   \item{Glucose}{Blood glucose level}
#'   \item{RBC}{Red blood cell count}
#'   \item{Sample}{Sample number}
#'   \item{Para}{Parasite count}
#' }
#' @source Sylvie Huijben 
"chabaudi"

#' Monthly incidence of influenza-like illness in Iceland between 1980 and 2009.
#'
#' A dataset containing the monthly ILI incidence 
#' in Iceland between 1980 and 2009.
#'
#' @format A data frame with 360 rows and 3 variables:
#' \describe{
#'   \item{month}{the month}
#'   \item{year}{the year}
#'   \item{ili}{ILI incidence}
##' }
#' @source Bjornstad ON, Viboud C. Timing and periodicity of influenza epidemics. Proceedings of the National Academy of Sciences. 2016 Nov 15;113(46):12899-901. <doi:10.1073/pnas.1616052113>
"icelandflu"


#' Weekly incidence of whooping cough in Philadelphia between 1925 and 1947.
#'
#' A dataset containing the weekly incidence 
#' incidence of whooping cough in Philadelphia between 1925 and 1947.
#'
#' @format A data frame with 1200 rows and 5 variables:
#' \describe{
#'   \item{YEAR}{the year}
#'   \item{WEEK}{the week}
#'   \item{PHILADELPHIA}{weekly whooping cough incidence}
#'   \item{TIME}{the time counter}
#'   \item{TM}{observation counter}
#' }
#' @source \url{https://www.tycho.pitt.edu}
"tywhooping"

#' Weekly incidence of scarlet fever in Philadelphia between 1914 and 1947.
#'
#' A dataset containing the weekly incidence 
#' incidence of scarlet fever in Philadelphia between 1914 and 1947.
#'
#' @format A data frame with 1774 rows and 4 variables:
#' \describe{
#'   \item{YEAR}{the year}
#'   \item{WEEK}{the week}
#'   \item{PHILADELPHIA}{weekly scarlet fever incidence}
#'   \item{TIME}{the time counter}
#' }
#' @source \url{https://www.tycho.pitt.edu}
"tyscarlet"

#' Weekly incidence of diphtheria in Philadelphia between 1914 and 1947.
#'
#' A dataset containing the weekly incidence 
#' incidence of diphtheria in Philadelphia between 1914 and 1947.
#'
#' @format A data frame with 1774 rows and 4 variables:
#' \describe{
#'   \item{YEAR}{the year}
#'   \item{WEEK}{the week}
#'   \item{PHILADELPHIA}{weekly diphtheria incidence}
#'   \item{TIME}{the time counter}
#' }
#' @source \url{https://www.tycho.pitt.edu}
"tydiphtheria"


#' Weekly incidence of measles in Philadelphia between 1914 and 1947.
#'
#' A dataset containing the weekly incidence 
#' incidence of measles in Philadelphia between 1914 and 1947.
#'
#' @format A data frame with 1774 rows and 4 variables:
#' \describe{
#'   \item{YEAR}{the year}
#'   \item{WEEK}{the week}
#'   \item{PHILADELPHIA}{weekly measles incidence}
#'   \item{TIME}{the time counter}
#' }
#' @source \url{https://www.tycho.pitt.edu}
"tymeasles"


#' Day of appearance of each measles case from 2003-04 outbreak in Niamey, Niger.
#'
#' A dataset containing the day of appearance of each measles case  in
#' Niamey, Niger during the 2003-04 outbreak.
#'
#' @format A data frame with 10,937 rows and 1 variables:
#' \describe{
#'   \item{day}{the day of appearance of each case since day of outbreak}
#' }
#' @source Grais et al. (2008) Time is of the essence: exploring a measles outbreak response vaccination in Niamey, Niger. Journal of the Royal Society Interface 5: 67-74. <doi:10.1098/rsif.2007.1038>
"niamey_daily"



#' Measles incidence across 40 US cities
#'
#' A dataset of Measles incidence across 40 US cities with relevant demographic data
#'
#' @format A data frame with 44,720 rows and 10 variables:
#' \describe{
#'   \item{biweek}{biweek of the year}
#'   \item{cases}{incidence}
#'   \item{year}{year}
#'   \item{loc}{city name}
#'   \item{pop}{population size}
#'   \item{rec}{susceptible recruits}
#'   \item{country}{country}
#'   \item{lon}{city longitude} 
#'   \item{lat}{city latitude} 
#'   \item{decimalYear}{time counter}
#' }
#' @source Dalziel et al. 2016. Persistent chaos of measles epidemics in the prevaccination United States caused by a small change in seasonal transmission patterns. PLoS Computational Biology 2016: e1004655.  <doi:10.1371/journal.pcbi.1004655>
"dalziel"


#' Filipendula rust data.
#'
#' Rust infection status of 162 populations of \emph{Filipendula ulmaria}
#' in a Swedish Island archipelago
#'
#' @format A data frame with 162 rows and 4 variables:
#' \describe{
#'   \item{y94}{infection status in 1994}
#'   \item{y95}{infection status in 1995}
#'   \item{X}{X coordinate}
#'   \item{Y}{Y coordinate}
#' }
#' @source Smith et al. 2003. Epidemiological patterns at multiple spatial scales: an 11-year study of a Triphragmium ulmariae -- Filipendula ulmaria metapopulation. Journal of Ecology, 91(5), pp.890-903.  <doi:10.1046/j.1365-2745.2003.00811.x>
"filipendula"

#' US 1975/76 ILI data.
#'
#' Influenza-like illness data for the lower 48 states and the District of Columbia during the 1975/76 
#' season dominated by A/H3N2/Victoria strain
#'
#' @format A data frame with 49 rows and 7 variables:
#' \describe{
#'   \item{State}{State number}
#'   \item{Acronym}{State code}
#'   \item{Pop}{Population size}
#'   \item{Latitude}{Latitude}
#'   \item{Longitude}{Longitude}
#'   \item{Start}{Week of start of epidemic}
#'   \item{Peak}{Week of peak of epidemic}
#' }
#' @source Viboud C, Bjornstad ON, Smith DL, Simonsen L, Miller MA, Grenfell BT (2006) Synchrony, waves, and spatial hierarchies in the spread of influenza. Science 312: 447-451.<doi:10.1126/science.1125237>
"usflu"



#' Burnett's Parasitoid-Host data.
#'
#' Data is of 22 generations of greenhouse white flies (\emph{Trialeurodes vaporariorum}) and its parasitoid, \emph{Encarsia formosa}.
#' Column names are self explanatory.
#'
#' @format A data frame with 22 rows and 7 variables:
#' \describe{
#'   \item{Generation}{ }
#'   \item{NumberofHostsExposed}{ }
#'   \item{NumberofHostsParasitized}{ }
#'   \item{NumberofHostsUnparasitized}{ }
#'   \item{NumberofParasiteEggsLaid}{ }
#'   \item{NumberofParasitesSearching}{ }
#'   \item{PercentageofHostsParasitized}{ }
#' }
#' @source Burnett, T. A. (1958) Model of host-parasite interaction Proceedings of the 10th International Congress, Entomology, 1958, 2, 679-686
"burnett"

#c10

#' Raccoon rabies data.
#'
#' Data is the average monthly number of reported cases of rabid raccoons across all counties within each of 
#' 11 east coast US states the time line is from the first reported case in each state (starting in late 1970s for West Virginia).
#'
#' @format A data frame with 208 rows and 12 variables:
#' \describe{
#'   \item{Month}{Month since rabies appearance in the state}
#'   \item{CT}{Connecticut}
#'   \item{DE}{Delaware}
#'   \item{MD}{Maryland}
#'   \item{MA}{Massachusetts}
#'   \item{NJ}{New Jersey}
#'   \item{NY}{New York}
#'   \item{NC}{North Carolina}
#'   \item{PA}{Pennsylvania}
#'   \item{RI}{Rhode Island}
#'   \item{VA}{Virginia}
#'   \item{WV}{West Virginia}
#' }
#' @source Childs et al. 2000. Predicting the local dynamics of epizootic rabies among raccoons in the United States Proceedings of the National Academy of Sciences 97:13666-13671. <doi:10.1073/pnas.240326697>
"rabies"

#cXX
#' Weekly whooping cough incidence from 1900-1937 in Copenhagen, Denmark.
#'
#' A dataset containing the weekly incidence of whooping cough from
#' Copenhagen, Denmark between January 1900 and December 1937
#'
#' @format A data frame with 1982 rows and 9 variables:
#' \describe{
#'   \item{date}{date}
#'   \item{births}{births}
#'   \item{day}{day of month}
#'   \item{month}{month of year}
#'   \item{year}{year}
#'   \item{cases}{weekly incidence }
#'   \item{deaths}{weekly deaths}
#'   \item{popsize}{weekly population size interpolated from census data}
#' }
#' @source Lavine et al. 2013. Immune boosting explains regime- shifts in prevaccine-era pertussis dynamics. PLoS ONE, 8(8):e72086. <doi:10.1371/journal.pone.0072086>
"pertcop"

#c12

#' Euthamia graminifolia rust data.
#'
#' Data on a fungal pathogen of the aster Euthamia graminifolia collected by Jennifer Keslow.
#'
#' @format A data frame with 360 rows and 8 variables:
#' \describe{
#'   \item{block}{the block}
#'   \item{row}{row}
#'   \item{plot}{plot within block}
#'   \item{xloc}{x coordinates}
#'   \item{yloc}{y coordinate}
#'   \item{comp}{plot composition}
#'   \item{water}{treatment: dry or wet}
#'   \item{score}{the rust score}
#' }
"euthamia"

#' Antler smut on wild campion.
#'
#' Data on a fungal pathogen of the wild campion collected by Janis Antonovics.
#'
#' @format A data frame with 876 rows and 5 variables:
#' \describe{
#'   \item{X}{road segment number}
#'   \item{lat}{latitude}
#'   \item{long}{longitude}
#'   \item{hmean}{number of healthy plants}
#'   \item{dmean}{number of diseased plants}
#' }
#' @source Antonovics, J. 2004. Long-term study of a plant-pathogen metapopulation. In: Hanski, Ilkka, and Oscar E. Gaggiotti. Ecology, genetics, and evolution of metapopulations. Academic Press. 
"silene"

#' Defoliated by gypsy moth each in northeast US 1975-2002.
#'
#' A list containing the fraction of forest defoliated by the gypsy moth
#' in 20km x 20km pixels across northeast US in each year between 1975 and 2002.
#'
#' @format A list with two matrices each with 1086 rows:
#' \describe{
#'   \item{xy}{A matrix with two columns representing UTM coordinates}
#'  \item{defoliation}{A matrix with 28 columns representing pixel-wise defoliation between 1976 and 2002}
#' }
#' @source Bjornstad, O. N., Robinet, C., & Liebhold, A. M. (2010). Geographic variation in North American gypsy moth cycles: subharmonics, generalist predators, and spatial coupling. Ecology, 91(1), 106-118. <doi:10.1890/08-1246.1>
"gypsymoth"

#c11b
#' Cumulative death count of harbor seals from CDV.
#'
#' The cumulative count of dead seals washed ashore across
#' 25 Northern European areas during the 2002 epidemic starting in
#' May and running through the end of the year. 
#'
#' @format A list with three items:
#' \describe{
#' \item{coord}{A data frame with 3 columns and 25 rows. Column location 
#' represents location name, latitude is latitude and longitude is longitude}
#'   \item{ts}{A data frame with 26 columns and 269 rows. The first column is the day 
#'  since beginning of outbreak and the next 25 columns are cumulative count of stranded
#' seal carcasses}
#' \item{fs}{A 25-by-25 matrix representing the seaway friction distance among the haulouts.}
#' }
#' @source Harding, K. C., Härkönen, T. and Caswell, H. (2002), The 2002 European seal plague: epidemiology and population consequences. Ecology Letters, 5: 727-732. <doi:10.1046/j.1461-0248.2002.00390.x>
"pdv"

#' Rabies month of first appearance across Connecticut.
#'
#' First month of report of raccoon rabies for 168 townships in Connecticut
#' from March 1991 through January 1995.
#'
#' @format A data frame with 168 rows and 3 variables:
#' \describe{
#'   \item{x}{Longitudinal distance (in km) from first township of appearance}
#'   \item{y}{Latitudinal distance (in km) from first township of appearance}
#'   \item{month}{month since first appearance in the state in March 1991}
#' }
#' @source Waller, L. A., & Gotway, C. A. (2004). Applied spatial statistics for public health data (Vol. 368). John Wiley & Sons.
"waller"

#' Measles in England and Wales 1944-1994.
#'
#' The weekly reported cases of measles 
#' in each of the 354 cities and villages
#' between 1944 and 1994. 
#'
#' @format A list with four items:
#' \describe{
#' \item{measles}{A matrix with 354 rows and 2661 columns. Rows 
#' represents community and columns represent week.}
#' \item{ps}{A matrix with 354 rows and 41 columns. Rows 
#' represents community and columns represent population size 
#' for each of the years in the data set.}
#' \item{longlat}{A matrix with 354 rows and two columns. Rows 
#' represents community and columns longitude and latitude.}
#' \item{year}{A vector of length 354 that represents time for each week as yearly decimals.} 
#' \item{coverage}{A vector of length 2661 that represents reported annual vaccine coverage for each week.} 
#' }
#' @source  Grenfell, B.T., Bjornstad, O.N., & Kappey, J. 2001. Travelling waves and spatial hierarchies in measles epidemics. Nature 414: 716-723. <doi:10.1038/414716a>
#' @source Lau, M.S.Y., Becker, A.D, Korevaar, H.M., Caudron, Q., Shaw, D.J., Metcalf, C.J.E., Bjornstad, O.N. and Grenfell, B.T. 2020. A competing-risks model explains hierarchical spatial coupling of measles epidemics en route to national elimination. Nature Ecology & Evolution. <doi:10.1038/s41559-020-1186-6>
"m4494"

#c12

#' Colorado Springs network
#'
#' Network and individual characteristics among 749 sex workers and clients in Colorado Springs as surveyed between 1988 and 1991
#'
#' @format A list of two items the first ($nodes) is a frame with 749 rows and 6 variables, The second ($cm) is a 749 x 749 relational matrix of presence/absence of sexual contacts among each pair of individuals.
#' \describe{
#' \item{$nodes$id}{individual identifier}
#' \item{$nodes$gender}{gender; 1 = female, 2 = male}
#' \item{$nodes$sex.worker}{sex worker status; 1 = yes, 0 = no}
#' \item{$nodes$pimp}{pimp status; 1 = yes, 0 = no}
#' \item{$nodes$sex.work.client}{client status; 1 = yes, 2 = no}
#' \item{$nodes$type}{node classifier; 1 = client, 2 = worker, 3 = both}
#' \item{$cm}{the relational (contact) matrix among the individuals in the network.}}
#' @source Woodhouse et al. (1994) Mapping a social network of heterosexuals at high risk for HIV Infection. AIDS 8:1331-1336. doi:10.1097/00002030-199409000-00018
#' @source Klovdahl et al. (1994) Social networks and infectious disease: The Colorado Springs study. Social Science and Medicine 38:79-88. <doi:10.1016/0277-9536(94)90302-6>
#' @source \url{https://opr.princeton.edu/archive/p90/}
"cspring"



#c13

#' Bordetella bronchiseptica in rabbit kittens.
#'
#' Data on Bordetella bronchiseptica in rabbit kittens in a breeding facility.
#'
#' @format A data frame with 494 rows and 8 variables:
#' \describe{
#'   \item{Facility}{breeding facility}
#'   \item{sick}{infection status}
#'   \item{Date}{date sampled}
#'   \item{Animal.code}{animal identifier}
#'   \item{msick}{dams infection status}
#'   \item{Litter}{litter identifier}
#'   \item{CFU}{bacterial count}
#'   \item{Description}{unique litter identifier}
#' }
#' @source Long et al (2010) Identifying the Age Cohort Responsible for Transmission in a Natural Outbreak of Bordetella bronchiseptica. PLoS Pathogens 6(12): e1001224. <doi:10.1371/journal.ppat.1001224>
"litter"

#c14

#' FIV infection in cats.
#'
#' Immunological measures on cats infected with different strains of FIV
#'
#' @format A data frame with 238 rows and 18 variables:
#' \describe{
#'   \item{Id}{Individual identifier}
#'   \item{CD4}{CD4 cell count}
#'   \item{CD8B}{CD8B cell count}
#'   \item{CD25}{CD25 cell count}
#'   \item{FAS_L}{FAS ligand}
#'   \item{FAS}{FAS}
#'   \item{IFNg}{Interferon gamma}
#'   \item{IL_10}{Interleukin 10}
#'   \item{IL_12}{Interleukin 12}
#'   \item{IL_4}{Interleukin 4}
#'   \item{lymphocyte}{lymphocyte count}
#'   \item{neutrophils}{neutrophil count}
#'   \item{TNF_a}{Tumor necrosis factor}
#'   \item{provirus}{provirus count}
#'   \item{viremia}{viremia}
#'   \item{Day}{day}
#'   \item{No}{unique identifier}
#'   \item{Treatment}{Experimental treatment}
#' }
#' @source Roy et al. 2009. Multivariate statistical analyses demonstrate unique host immune responses to single and dual lentiviral infection. PLoS one 4, e7359. <doi:10.1371/journal.pone.0007359>
"fiv"
