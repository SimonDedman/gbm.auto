# how to include try() so gbm.auto doesn't crash on bad values.

# possible approach:
assign(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l),
       try(
         gbm.step.sd(data = samples,
                     gbm.x = expvar,
                     gbm.y = brvcol,
                     family = fam1,
                     tree.complexity = j,
                     learning.rate = k,
                     bag.fraction = l,
                     n.trees = ntf1,
                     {if (!is.null(offset)) offset = grv_yes$offset},
                     ...)
       )
)
# THIS WORKS.


j <- list(1,2,3,"boo",5)
for (i in 1:length(j)) { # i <- 1
  print(i)
  assign(paste0("res", i), try(log(j[[i]]))) # prints error message then "try-error" then continues
  print(class(get(paste0("res", i))))
} # works

# If it fails, I need to just move onto the next loop I think? So:
if (class(get(paste0("res", i), try(log(j[[i]])))) == "try-error") next
if (class(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))) == "try-error") next

# NO, because report column indices are based on a loop-agnostic counter, n:
n = 1 # Print counter for all loops of BRT combos & best bin BRT choice #L432 (same for gaus L431)
if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) {Report[1:8,(3 + n)] <- c()} # L560
colnames(Report)[3 + n] <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l) # L569
n <- n + 1   # Add to print counter. L584

# Could instead start with
n = 0
# then move
n <- n + 1
# up to L530, BEFORE gbm.step.sd is run

# SO THIS APPROACH WORKS.
# 1: Change to n = 0, move n <- n + 1 up, do the same for m.
# 2: Wrap try() around gbm.step.sd and add if(class()) line
