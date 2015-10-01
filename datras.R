## DATRAS in R
# getAphia or getSpecies: export csv of unique species codes in alphabetical order
# then manually send to lookup website, get response csv
# putAphia or putSpecies: upload response csv then vlookup input codes and results against original dbase

# table layout: convert raw datras file to 3 sheets format.
# Assuming raw format is in the weird 3 sheets in the first place
# and that's not just excel
# use Site/StationNo/Year as index & put as 1st column


Sum of HL no @ length for all length classes selected.
Actually needs to be sum of (#@L * L/W conversion = total weight per length class)
  Can use CatCatchWeight? in grams. "catch weight in grams per category, or weight per haul per hour for CPUE data" (which is it?)
  CatCatchWeight is the sum of weights comprised of the TotalNo of fish of the same sex.
  
  Irrespective of Length:
    So for each sex, total weight per site =
    (average CatCatchWeight)/(average TotalNo) == M
  +
    (average CatCatchWeight)/(average TotalNo) == F
  +
    (average CatCatchWeight)/(average TotalNo) == (other sex?)
  
  For LngtClass between X & Y: same calculations as above. Problem: blonde ray reported max weight for 1 fish: 9800 (g?). Min: 5. Those two together in a trawl = 9805 (g?). Average: 4902.5
  Better to do as numbers?
  
  
  # also, instead of average weight per CatCatchVal, use proportion of length of Cat, e.g.
  
# individual weights = individual's CatCatchWgt * (individuals LngClass / totalLngClass in that CatCatchWgt) e.g.
Sex   #   CatCatchWgt    LngClass
M     2   10           1000
M     2   10           100
# 1st one: 10 * (1000/1100) = 9.09
# 1st one: 10 * (100/1100) = 0.909

# See Doug Beare's Datras R package to see what's in there already