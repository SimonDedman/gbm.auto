gbm.auto
=======

Automatically runs numerous processes from R packages 'gbm' and 'dismo' and script 'gbm.utils.R' which contains Elith et al.'s functions: roc, calibration, and gbm.predict.grids, as well as running my packages gbm.bfcheck, gbm.basemap, gbm.map, gbm.rsb, gbm.cons, gbm.valuemap, and gbm.loop.  

Especially on Linux systems it is recommended to type, in terminal:
sudo apt install libgeos-dev
sudo apt install libproj-dev
sudo apt install libgdal-dev
then manually install rgeos and rgdal in R/Rstudio.

Also see each script's Details section in the manual pages, as these frequently contain tips or common bugfixes.

I strongly recommend that you download papers 1 to 5 (or just the doctoral thesis) on http://www.simondedman.com, with emphasis on P4 (the guide) and P1 (statistical background). Elith et al 2008 (http://refhub.elsevier.com/S0304-3800(15)00207-0/sbref0085) is also strongly recommended.

***

### gbm.auto.R: Automated Boosted Regression Tree modelling and mapping suite

Automates delta log normal boosted regression trees abundance prediction. Loops through all permutations of parameters provided (learning rate, tree complexity, bag fraction), chooses the best, then simplifies it. Generates line, dot and bar plots, and outputs these and the predictions and a report of all variables used, statistics for tests, variable interactions, predictors used and dropped, etc. If selected, generates predicted abundance maps, and Unrepresentativeness surfaces.  

***

### gbm.bfcheck.R: Calculates minimum Bag Fraction size for gbm.auto

Provides minimum bag fractions for gbm.auto, preventing failure due to bf & samples rows limit.  

***

### gbm.basemap.R: Creates Basemaps for Gbm.auto mapping from your data range

Downloads unzips crops & saves NOAAs global coastline shapefiles to user-set box. Use for 'shape' in gbm.map. If downloading in RStudio uncheck "Use secure download method for HTTP" in Tools > Global Options > Packages.  

***

### gbm.map.R: Maps of predicted abundance from Boosted Regression Tree modelling

Generates maps from the outputs of gbm.step then gbm.predict.grids, handled automatically within gbm.auto but can be run alone, and generates representativeness surfaces from the output of gbm.rsb.  

***

### gbm.rsb.R: Representativeness Surface Builder

Loops through explanatory variables comparing their histogram in 'samples' to their histogram in 'grids' to see how well the explanatory variable range in samples represents the range being predicted to in grids. Assigns a representativeness score per variable per site in grids, and takes the average score per site if there's more than 1 expvar. Saves this to a CSV; it's plotted by gbm.map if called in gbm.auto. This shows you which areas have the most and least representative coverage by samples, therefore where you can have the most/least confidence in the predictions from gbm.predict.grids.
Can be called directly, and choosing a subset of expvars allows one to see their individual / collective representativeness.  

***

### gbm.cons.R: Conservation Area Mapping

Runs gbm.auto for multiple subsets of the same overall dataset and scales the combined results, leading to maps which highlight areas of high conservation importance for multiple species in the same study area e.g. using juvenile and adult female subsets to locate candidate nursery grounds and spawning areas respectively.  

***

### gbm.valuemap.R: Decision Support Tool that generates (Marine) Protected Area options using species predicted abundance maps

Scales response variable data, maps a user-defined explanatory variable to be avoided, e.g. fishing effort, combines them into a map showing areas to preferentially close. Bpa, the precautionary biomass required to protect the spawning stock, is used to calculate MPA size. MPA is then grown to add subsequent species starting from the most conservationally at-risk species, resulting in one MPA map per species, and a multicolour MPA map of all. All maps list the percentage of the avoidvariables total that is overlapped by the MPA in the map legend.  

***

### gbm.loop.R: Calculate Coefficient Of Variation surfaces for gbm.auto predictions

Processes a user-specified number of loops through the same gbm.auto parameter combinations and calculates the Coefficient Of Variation in the predicted abundance scores for each site aka cell. This can be mapped to spatially demonstrate the output variance range.  

***

### gbm.utils.R: roc, calibration & gbm.predict.grids functions bundle.

Now separated into individual functions.

***

# ToDo List:

See GitHub issues section https://github.com/SimonDedman/gbm.auto/issues
Feel free to contribute to this!
