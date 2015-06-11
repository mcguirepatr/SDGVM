##########################
#
# Parameters for netcdf outputs required by the TERRABITES project
# AWalker
# March 2015 
#
##########################


### Annual data
# Total Veg C, Abg soil litter, Blg soil litter & SOM,   
# trait workshop wants SOM and surface litter seperate - SDGVM does them combined in scn
#                                                      - SDGVM surface litter: C1 & C5, below litter: C2 & C6, SOM: C3,4,7&8 
# trait workshop wants monthly C flux from fire - SDGVM only simulates fire on an annual basis

v1 <- list(
  name  = 'cVeg',
  file  = 'biot',
  lname = 'Carbon in vegetation',
  units = 'kg C m-2',
  scale = 1/1000
)
v2 <- list(
  name  = 'cLitter',
  file  = 'abg_litter',
  lname = 'Carbon in aboveground litter',
  units = 'kg C m-2',
  scale = 1/1000
)
v3 <- list(
  name  = 'cSoil',
  file  = 'blg_c',
  lname = 'Carbon in soil (belowground litter and organic matter)',
  units = 'kg C m-2',
  scale = 1/1000
)
v4 <- list(
  name  = 'cLeaf',
  file  = 'leafc',
  lname = 'Carbon in leaves',
  units = 'kg C m-2',
  scale = 1/1000
)
v5 <- list(
  name = 'cRoot',
  file = 'rootbio',
  lname = 'Carbon in roots',
  units = 'kg C m-2',
  scale = 1/1000
)
v6 <- list(
  name = 'fFire',
  file = 'fcn',
  lname = 'Carbon released by fire',
  units = 'kg C m-2 s-1',
  scale = 1/1000/(30*24*3600)
)
afiles <- list(v1,v2,v3,v4,v5,v6)


### Monthly data
# total soil water (inc. snow), runoff & drainage, Total ET, Sensible Heat, Surface tmp, GPP, Ra, NPP, Rh, CO2 from fire,  
v1 <- list(
  name  = 'mrso',
  file  = c('swc','snw'),
  lname = 'Total soil water',
  units = 'kg m-2',
  scale = 1
)
v2 <- list(
  name  = 'mrro',
  file  = 'rof',
  lname = 'Total Runoff (including drainage)',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600)
)
v3 <- list(
  name  = 'evapotrans',
  file  = 'evt',
  lname = 'Total Evapotranspiration',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600)
)
# v4 <- list(
#   name  = 'sh',
#   file  = NA,
#   lname = 'Sensible Heat Flux',
#   units = 'W m-2'
# )
# v5 <- list(
#   name = 'Ta',
#   file = 'tmp',
#   lname = 'Air Temperature',
#   units = 'K'
# )
v6 <- list(
  name = 'gpp',
  file = 'gpp',
  lname = 'Gross Primary Production',
  units = 'kg C m-2 s-1',
  scale = 1/1000/(30*24*3600)
)
v7 <- list(
  name  = 'ra',
  file  = 'rsp',
  lname = 'Autotrophic respiration',
  units = 'kg C m-2 s-1',
  scale = 1/1000/(30*24*3600)
)
v8 <- list(
  name  = 'npp',
  file  = 'npp',
  lname = 'Net Primary Production',
  units = 'kg C m-2 s-1',
  scale = 1/1000/(30*24*3600)
)
v9 <- list(
  name  = 'rh',
  file  = 'srp',
  lname = 'Heterotrophic respiration',
  units = 'kg C m-2 s-1',
  scale = 1/1000/(30*24*3600)
)
v10 <- list(
  name = 'nbp',
  file = 'nep',
  lname = 'Net Biome Production',
  units = 'kg C m-2 s-1',
  scale = 1/1000/(30*24*3600)
)
# the requested output is actually by pft - don't have!
v11 <- list(
  name = 'lai',
  file = 'lai',
  lname = 'Leaf Area Index',
  units = 'm2 m-2',
  scale = 1
)
mfiles <- list(v1,v2,v3,v6,v7,v8,v9,v10,v11)
# mfiles <- list(v2)




### Daily data
# runoff & drainage, Total ET, Sensible Heat, Surface tmp,
# d_names <- c('mrro','evapotrans','sh','Ta') 
# dfiles <- c('rof','evt',NA,'tmp')
# dncdf  <- 1:4



