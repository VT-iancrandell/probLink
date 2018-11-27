source("./src/mcmcLinkReconcile/rcppMCMCDev/01-loadDataAndSources.R")

# Compute pairs, get Ms and Us, and get the groupings

pairs <- compare.dedup(demographics[1:5000],
                       #exclude = 'ssn',
                       strcmp = TRUE,
                       blockfld = list(c("last_name", "first_name"), c("birth_yr", "birth_mo"), c("ssn"))
)
ms = c(.9, .7, .85, .4, .5, .7, .9, .9, .8)


fspairs <- fsWeights(pairs, m = ms)
# Do the U's a little more intelligently, compute the probability of drawing 2 of the same labels at random. I'm not doing it this way for now for comprability with recordlinkage
# us = sapply(pairs$data, computeU)
#us = fspairs$frequencies
table(fspairs$Wdata)

linkageGroups = setDT(getPairs(fspairs, min.weight = -999, single.rows = TRUE))
linkedIds = linkageGroups[,.(id1, id2)]
