# "Damaraland mole rats are not obligate cooperative breeders"

The files in this repository provide the R code necessary to replicate all of the analyses in the paper. All the data were collected from a population of mole-rats living around the Kuruman River Reserve, Northern Cape, South Africa. 

The R scripts are separated according to analysis:

## i) Status-related survivorship in females

Using multi-state Markov models to estimate 'survivorship' of females in one of three states: in-group non-breeder, dispersed single female, breeder. 

## ii) Early life growth, adult body mass

Modelling the effect of group size on skeletal (via teeth width) and body mass growth, as well as additional confirmatory analyses of adult mass.

## iii) The body condition of single females

Comparing the body condition of single females vs size-matched in-group non-breeders, through allometric scaling relationships of body mass on skeletal size. 

## iv) Within-group recruitment rates

Analysing the effect of group size on the recruitment rate of offspring across the duration of our study (longitudinal analysis). We also created new pairs experimental by introducing unrelated adult males to a portion of dispersed single females. Later in the script, we compare the rates of recruitment in these nascent pairs against that in established groups (experimental analysis). The idea of this latter step was to exiplicitly, test whether the absence of additional group members reduces the reproductive output of breeding females. 

## v) The timing of natal dispersal (the duration of philopatry)

Fitting a two-state Markov models to estimate the timing of natal dispersal/duration of philopatry for male and female non-breeders. The model is restricted to individuals first captured within their first year of life, so that timing can be approximately related to age. 
