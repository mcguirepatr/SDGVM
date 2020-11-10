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
  alma_name = 'AvgSurfT',
  file  = 'tmp',
  pft   = F,
  lname = 'Average Surface Air Temperature',
  units = 'K',
  scale = 273.15
)

ta <- list(
  name  = 'ta',
  alma_name = 'Tair',
  file  = 'tmp',
  pft   = F,
  lname = 'Near Surface Air Temperature',
  units = 'K',
  scale = 273.15
)

tsl <- list(
  name  = 'tsl',
  alma_name = 'SoilTemp',
  file  = 'tmp',
  pft   = F,
  lname = 'Near-Surface Air Temperature',
  units = 'K',
  scale = 273.15
)

pr <- list(
  name  = 'pr',
  alma_name = 'Precip',
  file  = 'prc',
  pft   = F,
  lname = 'Precipitation',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600),
  scale_day = 30
)

hur <- list(
  name  = 'hur',
  alma_name = 'RH',
  file  = 'hum',
  pft   = F,
  lname = 'Relative humidity',
  units = '%',
  scale = 1
)

rsds<- list(
  name  = 'rsds',
  alma_name = 'SWdown',
  file  = c('qdr','qdf'),
  pft   = F,
  lname = 'Surface Downwelling Shortwave Radiation',
  units = 'W m-2',
  scale = 1/2.3
)

pot_evapotrans <- list(
  name  = 'pot_evapotrans',
  alma_name = 'PotEvap',
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

cVegd <- list(
  name  = 'cVeg',
  alma_name  = NULL,
  file  = '',
  pft   = F,
  lname = 'Carbon in Vegetation',
  units = 'kg m-2',
  scale = 1/1000
)

cLitterd <- list(
  name  = 'cLitter',
  alma_name  = NULL,
  file  = c('c01', 'c02', 'c05', 'c06' ),
  pft   = F,
  lname = 'Carbon in Above and Below Ground Litter Pools',
  units = 'kg m-2',
  scale = 1/1000
)

cSoild <- list(
  name  = 'cSoil',
  alma_name  = 'TotSoilCarb',
  file  = c('c01', 'c02', 'c03', 'c04', 'c05', 'c06', 'c07', 'c08' ),
  pft   = F,
  lname = 'Carbon in Soil (including below-ground litter)',
  units = 'kg m-2',
  scale = 1/1000
)

cLeafd <- list(
  name  = 'cLeaf',
  alma_name  = NULL,
  file  = 'lbm',
  pft   = F,
  lname = 'Carbon in Leaves',
  units = 'kg m-2',
  scale = 1/1000
)

cRootd <- list(
  name  = 'cRoot',
  alma_name  = NULL,
  file  = 'rbm',
  pft   = F,
  lname = 'Carbon in Live Roots',
  units = 'kg m-2',
  scale = 1/1000
)

cWoodd <- list(
  name  = 'cWood',
  alma_name  = NULL,
  file  = 'sbm',
  pft   = F,
  lname = 'Carbon in Sapwood',
  units = 'kg m-2',
  scale = 1/1000
)

cMiscd <- list(
  name  = 'cMisc',
  alma_name  = NULL,
  file  = 'nps',
  pft   = F,
  lname = 'Carbon Mass in Other Living Compartments on Land',
  units = 'kg m-2',
  scale = 1/1000
)

fVegLitter <- list(
  name  = 'fVegLitter',
  alma_name  = NULL,
  file  = c('f01', 'f02', 'f07', 'f08' ),
  pft   = F,
  lname = 'Total Carbon Mass Flux from Vegetation to Litter',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600),
  scale_day = 30
)

fLitterSoil <- list(
  name  = 'fLitterSoil',
  alma_name  = NULL,
  file  = c('f03', 'f04', 'f05', 'f09', 'f10', 'f11' ),
  pft   = F,
  lname = 'Total Carbon Mass Flux from Litter to Soil',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600),
  scale_day = 30
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


### Monthly data
# total soil water (inc. snow), runoff & drainage, Total ET, Sensible Heat, Surface tmp, GPP, Ra, NPP, Rh, CO2 from fire,  
mrso <- list(
  name  = 'mrso',
  alma_name = NULL, 
  #file  = c('swc','snw'),
  file  = 'swc',
  pft   = F,
  lname = 'Total Soil Water Content',
  units = 'kg m-2',
  scale = 1
)

rzwc <- list(
  name  = 'rzwc',
  alma_name = 'RootMoist',
  file  = 'swc',
  pft   = F,
  lname = 'Root zone soil moisture (water)',
  units = 'kg m-2',
  scale = 1
)

mrsos <- list(
  name  = 'mrsos',
  alma_name = 'SoilMoistTop',
  file  = 'ssm',
  pft   = F,
  lname = 'Moisture (Water) in Top Soil Layer',
  units = 'kg m-2',
  scale = 1
)

mrsow <- list(
  name  = 'mrsow',
  alma_name = 'SoilWet',
  file  = 'swf',
  pft   = F,
  lname = 'Total Soil Wetness',
  units = '-',
  scale = 1
)

snw <- list(
  name  = 'snw',
  alma_name = 'SWE',
  file  = 'snw',
  pft   = F,
  lname = 'Snow Water Equivalent',
  units = 'kg m-2',
  scale = 1
)

mrro <- list(
  name  = 'mrro',
  alma_name = 'Qr',
  file  = 'rof',
  pft   = F,
  lname = 'Total Runoff',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600),
  scale_day = 30
)

evapotrans <- list(
  name  = 'evapotrans',
  alma_name = 'Evap',
  file  = 'evt',
  pft   = F,
  lname = 'Total Evapo-Transpiration',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600),
  scale_day = 30
)

hfls <- list(
  name  = 'hfls',
  alma_name = 'Qle',
  file  = 'evt',
  pft   = F,
  lname = 'Latent heat flux',
  units = 'W m-2',
  scale = 2.257e6/(30*24*3600),
  scale_day = 30,
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions.' 
)

tran <- list(
  name  = 'tran',
  alma_name = 'Tveg',
  file  = 'trn',
  pft   = F,
  lname = 'Transpiration',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600),
  scale_day = 30
)

evapotranspft <- list(
  name  = 'evapotranspft',
  alma_name = 'Evap_pft',
  file  = 'evt',
  pft   = T,
  lname = 'Vegtype level evapotranspiration',
  units = 'W m-2, per unit land area occupied by the PFT',
  scale = 2.257e6/(30*24*3600),
  scale_day = 30,
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions.' 
)

evapo <- list(
  name  = 'evapo',
  alma_name = 'Esoil_pft',
  file  = 'bse',
  pft   = T,
  lname = 'Vegtype level Soil evaporation',
  units = 'W m-2',
  scale = 2.257e6/(30*24*3600),
  scale_day = 30,
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions. ' 
)

es <- list(
  name  = 'es',
  alma_name = 'Esoil',
  file  = 'bse',
  pft   = F,
  lname = 'Bare soil evaporation',
  units = 'kg m-2 s-1',
  scale = 1/(30*24*3600),
  scale_day = 30,
  notes = ''
)

transpft <- list(
  name  = 'transpft',
  alma_name = 'Tveg_pft',
  file  = 'trn',
  pft   = T,
  lname = 'Vegtype level transpiration',
  units = 'W m-2, per unit land area occupied by the PFT',
  scale = 2.257e6/(30*24*3600),
  scale_day = 30,
  notes = 'Converted to Wm-2 from kg m-2 s-1 (which is the standard SDGVM output) assuming a latent heat of vapourisation of 2257 kJ kg-1 at all times. This ignores the additional energy required for sublimation which is an order of magnitude smaller, but will lead to some discrepancy with other models when comparing this variable in cold regions. ' 
)

snow_depthpft <- list(
  name  = 'snow_depthpft',
  alma_name = 'SWE_pft',
  file  = 'snw',
  pft   = T,
  lname = 'Vegtype level snow water equivalent',
  units = 'm m-2, per unit land area occupied by the PFT',
  scale =  0.001 
)

gpp <- list(
  name = 'gpp',
  alma_name = 'GPP',
  file = 'gpp',
  pft   = F,
  lname = 'Gross Primary Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600),
  scale_day = 30
)

gpppft <- list(
  name = 'gpppft',
  file = 'gpp',
  pft   = T,
  lname = 'Vegtype level GPP',
  units = 'kg m-2 s-1, per unit land area occupied by the PFT',
  scale = 1/(1000*30*24*3600),
  scale_day = 30
)

ra <- list(
  name  = 'ra',
  alma_name = 'AutoResp',
  file  = 'rsp',
  pft   = F,
  lname = 'Autotrophic (Plant) respiration',
  units = 'kg m-2 s-1, per unit land area occupied by the PFT',
  scale = 1/(1000*30*24*3600),
  scale_day = 30
)

npp <- list(
  name  = 'npp',
  alma_name = 'NPP',
  file  = 'npp',
  pft   = F,
  lname = 'Net Primary Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600),
  scale_day = 30
)

npppft <- list(
  name  = 'npppft',
  file  = 'npp',
  pft   = T,
  lname = 'Vegtype level NPP',
  units = 'kg m-2 s-1, per unit land area occupied by the PFT',
  scale = 1/(1000*30*24*3600),
  scale_day = 30
)

rh <- list(
  name  = 'rh',
  alma_name = 'HeteroResp',
  file  = 'srp',
  pft   = F,
  lname = 'Heterotrophic Respiration',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600),
  scale_day = 30
)

nep <- list(
  name  = 'nep',
  alma_name = 'NEE',
  file  = 'nep',
  pft   = F,
  lname = 'Net Ecosystem Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600),
  scale_day = 30
)

nbp <- list(
  name  = 'nbp',
  alma_name = '',
  file  = 'nep',
  pft   = F,
  lname = 'Net Biome Production',
  units = 'kg m-2 s-1',
  scale = 1/(1000*30*24*3600),
  scale_day = 30,
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
  alma_name = NULL,
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

SDGVM_site_info_vars <- list(
  sand = list(
    name = 'SoilSand',
    col_name = 'sand',
    lname = 'Soil sand content',
    units = '%'
  ),
  silt = list(
    name = 'SoilSilt',
    col_name = 'silt',
    lname = 'Soil silt content',
    units = '%'
  ),
  BD = list(
    name = 'SoilBD',
    col_name = 'bulk',
    lname = 'Soil bulk density',
    units = 'g cm-3'
  ),
  carbon = list(
    name = 'SoilOrgC',
    col_name = 'orgc',
    lname = 'Soil organic carbon content for pedo-transfer functions',
    units = '%'
  ),
  wilt = list(
    name = 'SWCwilting',
    col_name = 'wp',
    lname = 'Soil water content at wilting point',
    units = '-'
  ),
  field = list(
    name = 'SWCfieldcapacity',
    col_name = 'fc',
    lname = 'Soil water content at field capacity',
    units = '-'
  ),
  sat = list(
    name = 'SWCsaturation',
    col_name = 'swc',
    lname = 'Soil water content at saturation',
    units = '-'
  )
)



### END ###
