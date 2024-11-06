# Script automatically generated on Thu Oct 31 12:26:57 2024
library(patRoon)
#Sys.setenv(JAVA_HOME="C:/Program Files/Java/jre1.8.0_421")
Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-23")

options(patRoon.path.SIRIUS = "C:/Users/holli/OneDrive - SUNY ESF/Documents/PhD Research/R Programs/patRoon/patRoon-bundle-2.3.2/user/library/patRoonExt/ext/sirius") # directory with the SIRIUS binaries
options(patRoon.path.OpenMS = "C:/Users/holli/OneDrive - SUNY ESF/Documents/PhD Research/R Programs/patRoon/patRoon-bundle-2.3.2/user/library/patRoonExt/ext/openms/bin") # directory with the OpenMS binaries
#options(patRoon.path.pngquant = "") # directory containing pngquant binary
options(patRoon.path.MetFragCL = "C:/Users/holli/OneDrive - SUNY ESF/Documents/PhD Research/R Programs/patRoon/patRoon-bundle-2.3.2/user/library/patRoonExt/ext/MetFragCommandLine.jar") # full location to the jar file
options(patRoon.path.MetFragCompTox = "C:/Users/holli/OneDrive - SUNY ESF/Documents/PhD Research/R Programs/patRoon/patRoon-bundle-2.3.2/user/library/patRoonExt/ext/CompToxWW.csv") # full location to desired CompTox CSV file
options(patRoon.path.MetFragPubChemLite = "C:/Users/holli/OneDrive - SUNY ESF/Documents/PhD Research/R Programs/patRoon/patRoon-bundle-2.3.2/user/library/patRoonExt/ext/PubChemLite.csv") # full location to desired PubChemLite CSV file
#options(patRoon.path.MetFragPubChemLite = "") # full location to PFAS DB (NOTE: configured like PubChemLite)
options(patRoon.path.BioTransformer = "C:/Users/holli/OneDrive - SUNY ESF/Documents/PhD Research/R Programs/patRoon/patRoon-bundle-2.3.2/user/library/patRoonExt/ext/biotransformer/biotransformer-3.0.0.jar")
options(patRoon.path.obabel = "C:/Users/holli/OneDrive - SUNY ESF/Documents/PhD Research/R Programs/patRoon/patRoon-bundle-2.3.2/user/library/patRoonExt/ext/openbabel") # directory with OpenBabel binaries
options(patRoon.path.pwiz = "C:/Users/holli/AppData/Local/Apps/ProteoWizard 3.0.24263.36e7380 64-bit/")

# -------------------------
# initialization
# -------------------------

workPath <- "C:/Users/holli/OneDrive - SUNY ESF/Documents/PhD Research/R Programs/patRoonAnnotate"
setwd(workPath)

# Load analysis table
anaInfoPos <- read.csv("analyses-pos.csv")


# -------------------------
# features
# -------------------------

# Find all features
# NOTE: see the reference manual for many more options
fList <- findFeatures(anaInfoPos, "openms", noiseThrInt = 100000, chromSNR = 5, chromFWHM = 5, minFWHM = 2, maxFWHM = 30)


# Group and align features between analyses
fGroups <- groupFeatures(fList, "openms", rtalign = TRUE, maxAlignRT = 30, 
                         maxAlignMZ = 0.01, maxGroupRT = 12, maxGroupMZ = 0.01)

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = 250000, absMinIntensity = 500000, relMinReplicateAbundance = 0.49,
                  maxReplicateIntRSD = NULL, blankThreshold = 5, removeBlanks = TRUE,
                  retentionRange = c(45,1020), mzRange = NULL)

# Reduce for testing
fGroups <- fGroups[, 1:500]

# -------------------------
# componentization
# -------------------------

# Perform automatic generation of components
componSet <- generateComponents(fGroups, "ramclustr", ionization = "positive")
fGroups <- selectIons(fGroups, componSet, c("[M+H]+"))

# -------------------------
# annotation
# -------------------------

# Retrieve MS peak lists
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.01, method = "hclust")

mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = 5, precursorMzWindow = 5,
                               avgFeatParams = avgMSListParams,
                               avgFGroupParams = avgMSListParams)
# Rule based filtering of MS peak lists. You may want to tweak this. See the manual for more information.
mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL, relMSMSIntThr = 0.05,
                  topMSPeaks = NULL, topMSMSPeaks = 25)

# Calculate formula candidates
formulas <- generateFormulas(fGroups, mslists, "genform", relMzDev = 10, elements = "CHNOS", oc = TRUE,
                             calculateFeatures = TRUE, featThresholdAnn = 0.75, 
                             extraOpts = c("kfer") #OKAY, THIS CAN ONLY BE kfer or kfer=ex
)


# Calculate compound structure candidates
compounds <- generateCompounds(fGroups, mslists, "metfrag",# adduct = "[M-H]-", 
                               database = "pubchemlite", topMost = 5)
compounds <- addFormulaScoring(compounds, formulas, updateScore = TRUE)


# -------------------------
# Transformation Products
# -------------------------

TPsLogic <- generateTPs("logic", fGroups, adduct = "[M+H]+")

# only screen for TPs
suspects <- convertToSuspects(TPsLogic, includeParents = FALSE)
# but keep all other feature groups as these may be parents
fGroupsScr <- screenSuspects(fGroups, suspects, onlyHits = FALSE)

fGroupsSusp <- annotateSuspects(fGroupsScr, MSPeakLists = mslists,
                                formulas = formulas, 
                                compounds = compounds)
#fGroupsSusp <- filter(fGroupsSusp, maxLevel = 4, onlyHits = TRUE)
fGroupsSusp@screenInfo[["estIDLevel"]]

# -------------------------
# reporting
# -------------------------

# Advanced report settings can be edited in the report.yml file.
report(fGroups, MSPeakLists = mslists, formulas = formulas, 
       compounds = compounds, 
       components = componSet,
       #TPs = TPsLogic, 
       settingsFile = "report.yml", openReport = TRUE)
