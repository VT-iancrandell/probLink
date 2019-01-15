# This script incorporates the blocking rules and possible deduped database assumption to create the comparison vectors for linkage
library(Rcpp)
library(igraph)
library(data.table)
library(RecordLinkage)

sourceDirectory("./probLinkMetropolis/R/", modifiedOnly = FALSE)
sourceCpp("./probLinkMetropolis/src/00-rcppFunctions.cpp")

# Load Data
data(RLdata10000)
recs = fread("./data/RLdata10000", colClasses = "character", na.strings = "")
recs[,fi := substring(fname_c1,1,1)]
recs[,irow := 1:nrow(recs)]

trueLabels = identity.RLdata10000


# Make a list containing all the exact and fuzzy rules

exactRules = list(c("fi", "lname_c1", "by", "bm"), 
                  c("fi", "lname_c1", "by", "bm"),
                  c("fi", "lname_c1", "bm", "bd"),
                  c("fi", "by", "bm", "bd"),
                  c("by","bm","bd"),
                  c("bm","bd"),
                  c("by","bd"),
                  c("by","bm"),
                  c("bm","bd"),
                  c("by","bd"),
                  c("by","bm"),
                  c("fname_c1","lname_c1","by"),
                  c("fname_c1","lname_c1","bm"),
                  c("fname_c1","lname_c1","bd"),
                  c("fi","by","bm","bd"),
                  c("fi","by","bm"),
                  c("fi","by","bd"),
                  c("fi","bm","bd"))
fuzzyRules = list("bd", "bm","by", "lname_c1", "fname_c1", "fname_c1", "fname_c1", "fname_c1", "lname_c1", "lname_c1", "lname_c1", NULL, NULL, NULL, NULL, NULL, NULL, NULL)

noRules = c(1:9)

allBonds = makeBondsFromList(recs, exactRules[noRules], fuzzyRules[noRules], dmax = 1)
rlComps = componentSizesFromBonds(allBonds)
table(rlComps$csize)

# Profile bonds

bondProf = profileBonds(trueLabels, allBonds, data = recs)
bondProf
# Write to file

fileName = "nineRuleBonds"

#fwrite(allBonds[,1:2], sprintf("./data/bondFiles/%s.csv", fileName))


#
# Make an unblocked bond file for the rl500. This is first to test the m estimation.
#

# Load Data
data(RLdata500)
fwrite(RLdata500, "./data/RLdata500.csv")
recs = fread("./data/RLdata500.csv", colClasses = "character", na.strings = "")
recs[,fi := substring(fname_c1,1,1)]
recs[,irow := 1:nrow(recs)]

trueLabels = identity.RLdata500

allBonds = makeBondsFromList(recs)
rlComps = componentSizesFromBonds(allBonds)
table(rlComps$csize)

# Profile bonds

bondProf = profileBonds(trueLabels, allBonds, data = recs)
bondProf
# Write to file

fileName = "rl500NoBlocks"

#fwrite(allBonds[,1:2], sprintf("./data/bondFiles/%s.csv", fileName))

#
# The above did poorly. Let's add blocks
#

# Load Data
data(RLdata500)
fwrite(RLdata500, "./data/RLdata500.csv")
recs = fread("./data/RLdata500.csv", colClasses = "character", na.strings = "")
recs[,fi := substring(fname_c1,1,1)]
recs[,irow := 1:nrow(recs)]

trueLabels = identity.RLdata500

exactRules = list(c("fi", "lname_c1", "by", "bm"), 
                  c("fi", "lname_c1", "by", "bm"),
                  c("fi", "lname_c1", "bm", "bd"),
                  c("fi", "by", "bm", "bd"),
                  c("by","bm","bd"),
                  c("bm","bd"),
                  c("by","bd"),
                  c("by","bm"),
                  c("bm","bd"),
                  c("by","bd"),
                  c("by","bm"),
                  c("fname_c1","lname_c1","by"),
                  c("fname_c1","lname_c1","bm"),
                  c("fname_c1","lname_c1","bd"),
                  c("fi","by","bm","bd"),
                  c("fi","by","bm"),
                  c("fi","by","bd"),
                  c("fi","bm","bd"))
fuzzyRules = list("bd", "bm","by", "lname_c1", "fname_c1", "fname_c1", "fname_c1", "fname_c1", "lname_c1", "lname_c1", "lname_c1", NULL, NULL, NULL, NULL, NULL, NULL, NULL)

allBonds = makeBondsFromList(recs, exactRules, fuzzyRules, dmax = 1)
rlComps = componentSizesFromBonds(allBonds)
table(rlComps$csize)


# Profile bonds

bondProf = profileBonds(trueLabels, allBonds, data = recs)
bondProf
# Write to file

fileName = "rl500WithBlocks"

fwrite(allBonds[,1:2], sprintf("./data/bondFiles/%s.csv", fileName))



