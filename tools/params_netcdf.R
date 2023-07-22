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

tas <- list(
  name  = 'tas',
  file  = 'tmp',
  pft   = F,
  lname = 'Near-Surface Air Temperature',
  units = 'K',
  scale = 273.15
)

tas_annual <- list(
  name  = 'tas_annual',
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

pr_annual<- list(
  name  = 'pr_annual',
  file  = 'prc',
  pft   = F,
  lname = 'Precipitation',
  units = 'kg m-2 s-1',
  scale = 1/(360*24*3600)
)

rsds<- list(
  name  = 'rsds',
  file  = c('qdr','qdf'),
  pft   = F,
  lname = 'Surface Downwelling Shortwave Radiation',
  units = 'W m-2',
  scale = 1/2.3
)

pot_evapotrans <- list(
  name  = 'pot_evapotrans',
  file  = 'pet',
  pft   = F,
  lname = 'Potential Evapo-Transpiration',
  units = 'kg m-2',
  scale = 1
)

cVeg <- list(
  name  = 'cVeg',
  file  = 'biot',
  pft   = F,
  lname = 'Carbon in Vegetation',
  units = 'kg m-2',
  scale = 1/1000
)

cVegpft <- list(
  name  = 'cVegpft',
  file  = 'bio',
  pft   = T,
  lname = 'Vegtype level Carbon in Vegetation',
  units = 'kg m-2, per unit land area occupied by the PFT',
  scale = 1/1000
)

cLitter <- list(
  name  = 'cLitter',
  file  = 'abg_litter',
  pft   = F,
  lname = 'Carbon in Above-ground Litter Pool',
  units = 'kg m-2',
  scale = 1/1000
)

cSoil <- list(
  name  = 'cSoil',
  file  = 'blg_c',
  pft   = F,
  lname = 'Carbon in Soil (including below-ground litter)',
  units = 'kg m-2',
  scale = 1/1000
)

cLeaf <- list(
  name  = 'cLeaf',
  file  = 'leafc',
  pft   = F,
  lname = 'Carbon in Leaves',
  units = 'kg m-2',
  scale = 1/1000
)

cRoot <- list(
  name  = 'cRoot',
  file  = 'rootbio',
  pft   = F,
  lname = 'Carbon in Roots',
  units = 'kg m-2',
  scale = 1/1000
)

fFire <- list(
  name  = 'fFire',
  file  = 'fcn',
  pft   = F,
  lname = 'CO2 Emission from Fire',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600),
  notes = 'SDGVM calculates this variable once per year, TRENDY output requires monthly data. Monthly data are calculated simply as the annual total divided by 12.',
  annmonth = T 
)

fFire_annual <- list(
  name  = 'fFire_annual',
  file  = 'fcn',
  pft   = F,
  lname = 'CO2 Emission from Fire',
  units = 'kg m-2 s-1',
  scale = 1/(1000*360*24*3600),
  notes = 'SDGVM calculates this variable once per year, TRENDY output requires monthly data. Monthly data are calculated simply as the annual total divided by 12.'
)

burntArea <- list(
  name  = 'burntArea',
  file  = 'fab',
  pft   = F,
  lname = 'Burnt Area Fraction',
  units = '-',
  scale = 1
)

fLuc <- list(
  name  = 'fLuc',
  file  = 'lulccc',
  pft   = F,
  lname = 'CO2 Flux to Atmosphere from Land Use Change',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600),
  notes = 'In this instance of SDGVM all above-ground biomass is assumed to be lost immediately to the atmosphere, and this is what this variable records. Below-ground biomass is assumed to go into the soil as litter and this variable does not track subsequent decomposition of that litter. SDGVM calculates this variable once per year, TRENDY output requires monthly data. Monthly data are calculated simply as the annual total divided by 12.',
  annmonth = T
)


fLuc_annual <- list(
  name  = 'fLuc_annual',
  file  = 'lulccc',
  pft   = F,
  lname = 'CO2 Flux to Atmosphere from Land Use Change',
  units = 'kg m-2 s-1',
  scale = 1/(1000*360*24*3600),
  notes = 'In this instance of SDGVM all above-ground biomass is assumed to be lost immediately to the atmosphere, and this is what this variable records. Below-ground biomass is assumed to go into the soil as litter and this variable does not track subsequent decomposition of that litter. SDGVM calculates this variable once per year, TRENDY output requires monthly data. Monthly data are calculated simply as the annual total divided by 12.'
)


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
  units = 'W m-2, per unit land area occupied by the PFT',
  scale = 2.257e6/(30*24*3600),
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions.' 
)

evapo <- list(
  name  = 'evapo',
  file  = 'bse',
  pft   = T,
  lname = 'Vegtype level Soil evaporation',
  units = 'W m-2',
  scale = 2.257e6/(30*24*3600),
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions. ' 
)

transpft <- list(
  name  = 'transpft',
  file  = 'trn',
  pft   = T,
  lname = 'Vegtype level transpiration',
  units = 'W m-2, per unit land area occupied by the PFT',
  scale = 2.257e6/(30*24*3600),
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions. ' 
)

snow_depthpft <- list(
  name  = 'snow_depthpft',
  file  = 'snw',
  pft   = T,
  lname = 'Vegtype level snow water equivalent',
  units = 'm m-2, per unit land area occupied by the PFT',
  scale =  0.001 
)

gpp <- list(
  name = 'gpp',
  file = 'gpp',
  pft   = F,
  lname = 'Gross Primary Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600)
)

gpp_annual <- list(
  name = 'gpp_annual',
  file = 'gpp',
  pft   = F,
  lname = 'Gross Primary Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*360*24*3600)
)

gpppft <- list(
  name = 'gpppft',
  file = 'gpp',
  pft   = T,
  lname = 'Vegtype level GPP',
  units = 'kg m-2 s-1, per unit land area occupied by the PFT',
  scale = 1/(1000*30*24*3600)
)

ra <- list(
  name  = 'ra',
  file  = 'rsp',
  pft   = F,
  lname = 'Autotrophic (Plant) respiration',
  units = 'kg m-2 s-1, per unit land area occupied by the PFT',
  scale = 1/(1000*30*24*3600)
)

ra_annual <- list(
  name  = 'ra_annual',
  file  = 'presp',
  pft   = F,
  lname = 'Autotrophic (Plant) respiration',
  units = 'kg m-2 s-1, per unit land area occupied by the PFT',
  scale = 1/(1000*360*24*3600)
)

lch_annual <- list(
  name  = 'lch_annual',
  file  = 'lch',
  pft   = F,
  lname = 'Leached soil carbon',
  units = 'kg m-2 s-1',
  scale = 1/(1000*360*24*3600)
)

npp <- list(
  name  = 'npp',
  file  = 'npp',
  pft   = F,
  lname = 'Net Primary Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600)
)

npp_annual <- list(
  name  = 'npp_annual',
  file  = 'npp',
  pft   = F,
  lname = 'Net Primary Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*360*24*3600)
)

npppft <- list(
  name  = 'npppft',
  file  = 'npp',
  pft   = T,
  lname = 'Vegtype level NPP',
  units = 'kg m-2 s-1, per unit land area occupied by the PFT',
  scale = 1/(1000*30*24*3600)
)

rh <- list(
  name  = 'rh',
  file  = 'srp',
  pft   = F,
  lname = 'Heterotrophic Respiration',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600)
)

rh_annual <- list(
  name  = 'rh_annual',
  file  = 'sresp',
  pft   = F,
  lname = 'Heterotrophic Respiration',
  units = 'kg m-2 s-1',
  scale = 1/(1000*360*24*3600)
)

nep <- list(
  name  = 'nep',
  file  = 'nep',
  pft   = F,
  lname = 'Net Ecosystem Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600)
)

nep_annual <- list(
  name  = 'nep_annual',
  file  = 'nep',
  pft   = F,
  lname = 'Net Ecosystem Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*360*24*3600)
)

nbp_annual <- list(
  name  = 'nbp_annual',
  file  = 'nbp',
  pft   = F,
  lname = 'Net Biome Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*360*24*3600),
  notes = 'These data include fire, lulcc, and leached C losses (which are annual fluxes distributed across the 12 months equally).'
)

nbppft <- list(
  name  = 'nbppft',
  file  = 'nep',
  pft   = T,
  lname = 'Vegtype level NBP',
  units = 'kg m-2 s-1, per unit land area occupied by the PFT',
  scale = 1/(1000*30*24*3600),
  notes = 'These data area per unit area covered by the PFT and do not include fire, lulcc, and leached C carbon losses. i.e. they are GPP - Ra - Rh'
)

lai <- list(
  name  = 'lai',
  file  = 'lai',
  pft   = F,
  lname = 'Leaf Area Index',
  units = '-',
  scale = 1
)

laipft <- list(
  name  = 'laipft',
  file  = 'lai',
  pft   = T,
  lname = 'Vegtype level Leaf Area Index',
  units = '-, per unit land area occupied by the PFT',
  scale = 1
)

landCoverFrac <- list(
  name = 'landCoverFrac',
  file = 'cov',
  pft   = T,
  lname = 'Fractional Cover of PFT',
  units = '-',
  scale = 1,
  notes = 'These are derived from the LUH2v2h cropland and pasture cover dataset, 850-2020, combined with the ESA CCI 2014 Land Cover maps (Poulter et al 2015) translated for the SDGVM PFT set. Cropland cover of the ESA dataset was reduced or increased according to LUH2v2h while grassland cover in ESA was only increased by the pasture cover in LUH2v2h so as not to remove natural grasslands. This likely high biased grassland cover in 1700. When deciding whether to run a site, SDGVM checks a high-resolution land-sea mask. If the grid cell is >50% land then the site is simulated assuming that the whole of the grid-cell is land. So the sum of the PFT landCoverFracs should always equal 1. Grid cells <50% land are not simulated.'
)




### Daily data
# runoff & drainage, Total ET, Sensible Heat, Surface tmp,
# d_names <- c('mrro','evapotrans','sh','Ta') 
# dfiles <- c('rof','evt',NA,'tmp')
# dncdf  <- 1:4



### END ###
