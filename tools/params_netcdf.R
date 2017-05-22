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

tas <- list(
  name  = 'tas',
  file  = 'tmp',
  pft   = F,
  lname = 'Near-Surface Air Temperature',
  units = 'K',
  scale = 273.15
)
pr<- list(
  name  = 'pr',
  file  = 'prc',
  pft   = F,
  lname = 'Precipitation',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600)
)
rsds<- list(
  name  = 'rsds',
  file  = c('qdr','qdf'),
  pft   = F,
  lname = 'Surface Downwelling Shortwave Radiation',
  units = 'W m-2',
  scale = 1/2.3
)


cVeg <- list(
  name  = 'cVeg',
  file  = 'biot',
  pft   = F,
  lname = 'Carbon in Vegetation',
  units = 'kg C m-2',
  scale = 1/1000
)
cVegpft <- list(
  name  = 'cVegpft',
  file  = 'biot',
  pft   = T,
  lname = 'Vegtype level Carbon in Vegetation',
  units = 'kg C m-2',
  scale = 1/1000
)
cLitter <- list(
  name  = 'cLitter',
  file  = 'abg_litter',
  pft   = F,
  lname = 'Carbon in Above-ground Litter Pool',
  units = 'kg C m-2',
  scale = 1/1000
)
cSoil <- list(
  name  = 'cSoil',
  file  = 'blg_c',
  pft   = F,
  lname = 'Carbon in Soil (including below-ground litter)',
  units = 'kg C m-2',
  scale = 1/1000
)
cLeaf <- list(
  name  = 'cLeaf',
  file  = 'leafc',
  pft   = F,
  lname = 'Carbon in Leaves',
  units = 'kg C m-2',
  scale = 1/1000
)
cRoot <- list(
  name = 'cRoot',
  file = 'rootbio',
  pft   = F,
  lname = 'Carbon in Roots',
  units = 'kg C m-2',
  scale = 1/1000
)
fFire <- list(
  name = 'fFire',
  file = 'fcn',
  pft   = F,
  lname = 'CO2 Emission from Fire',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600)),
  annmonth =NULL
)
burntArea <- list(
  name = 'burntArea',
  file = 'fab',
  pft   = F,
  lname = 'Burnt Area Fraction',
  units = '%',
  scale = 100
)
fLuc <- list(
  name = 'fLuc',
  file = 'flulccc',
  pft   = F,
  lname = 'CO2 Flux to Atmosphere from Land Use Change',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600)),
  notes = 'In this instance of SDGVM all above-ground biomass is a ssumed to be lost immediately to the atmosphere, and this is what this variable records. Below-ground biomass is assumed to go into the soil as litter and this variable does not track subsequent decomposition of that litter.',
  annmonth = NULL
)
#afiles <- list(v1,v2,v3,v4,v5,v6)


### Monthly data
# total soil water (inc. snow), runoff & drainage, Total ET, Sensible Heat, Surface tmp, GPP, Ra, NPP, Rh, CO2 from fire,  
mrso <- list(
  name  = 'mrso',
  file  = c('swc','snw'),
  pft   = F,
  lname = 'Total Soil Water Content',
  units = 'kg m-2',
  scale = 1
)
mrro <- list(
  name  = 'mrro',
  file  = 'rof',
  pft   = F,
  lname = 'Total Runoff',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600)
)
evapotrans <- list(
  name  = 'evapotrans',
  file  = 'evt',
  pft   = F,
  lname = 'Total Evapo-Transpiration',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600)
)
tran <- list(
  name  = 'tran',
  file  = 'trn',
  pft   = F,
  lname = 'Transpiration',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600)
)
evapotranspft <- list(
  name  = 'evapotranspft',
  file  = 'evt',
  pft   = T,
  lname = 'Vegtype level evapotranspiration',
  units = 'W m-2',
  scale = 2.257e6,
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions.' 
)
evapo <- list(
  name  = 'evapo',
  file  = 'bse',
  pft   = T,
  lname = 'Vegtype level Soil evaporation',
  units = 'W m-2',
  scale = 2.257e6,
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions. ' 
)
transpft <- list(
  name  = 'transpft',
  file  = 'trn',
  pft   = T,
  lname = 'Vegtype level transpiration',
  units = 'W m-2',
  scale = 2.257e6,
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions. ' 
)
swepft <- list(
  name  = 'swepft',
  file  = 'snw',
  pft   = T,
  lname = 'Vegtype level snow water equivalent',
  units = 'kg m-2',
  scale = NA
)
# v4 <- list(
#   name  = 'sh',
#   file  = NA,
#   pft   = F,
#   lname = 'Sensible Heat Flux',
#   units = 'W m-2'
# )
# v5 <- list(
#   name = 'Ta',
#   file = 'tmp',
#   pft   = F,
#   lname = 'Air Temperature',
#   units = 'K'
# )
gpp <- list(
  name = 'gpp',
  file = 'gpp',
  pft   = F,
  lname = 'Gross Primary Production',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600))
)
gpppft <- list(
  name = 'gpppft',
  file = 'gpp',
  pft   = T,
  lname = 'Vegtype level GPP',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600))
)
ra <- list(
  name  = 'ra',
  file  = 'rsp',
  pft   = F,
  lname = 'Autotrophic (Plant) respiration',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600))
)
npp <- list(
  name  = 'npp',
  file  = 'npp',
  pft   = F,
  lname = 'Net Primary Production',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600))
)
npppft <- list(
  name  = 'npppft',
  file  = 'npp',
  pft   = T,
  lname = 'Vegtype level NPP',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600))
)
rh <- list(
  name  = 'rh',
  file  = 'srp',
  pft   = F,
  lname = 'Heterotrophic Respiration',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600))
)
nbp <- list(
  name = 'nbp',
  file = 'nep',
  pft   = F,
  lname = 'Net Biome Production',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600))
)
nbppft <- list(
  name = 'nbppft',
  file = 'nep',
  pft   = T,
  lname = 'Vegtype level NBP',
  units = 'kg C m-2 s-1',
  scale = 1/(1000/(30*24*3600))
)
lai <- list(
  name = 'lai',
  file = 'lai',
  pft   = T,
  lname = 'Leaf Area Index',
  units = '-',
  scale = 1
)
landCoverFrac <- list(
  name = 'landCoverFrac',
  file = 'cov',
  pft   = T,
  lname = 'Fractional Cover of PFT',
  units = '-',
  scale = 1,
  notes = 'When deciding whether to run a site, SDGVM checks a high-resolution land-sea mask. if the grid cell is >50% land then the site is simulated assuming that the whole of the grid-cell is land. So the sum of the PFT landCoverFracs should always equal 1. Grid cells <50% land are not simulated.'
)
#mfiles <- list(v1,v2,v3,v6,v7,v8,v9,v10,v11)




### Daily data
# runoff & drainage, Total ET, Sensible Heat, Surface tmp,
# d_names <- c('mrro','evapotrans','sh','Ta') 
# dfiles <- c('rof','evt',NA,'tmp')
# dncdf  <- 1:4



