## Code accompanying submission of manuscript to Methods in Ecology and Evolution 

### Holly I. Niven, Jana W. E. Jeglinski, Geert Aarts, Ewan D. Wakefield, Jason Matthiopoulos (2025)
### Towards biologically realistic estimates of home range and spatial exposure for colonial animals
### Methods in Ecology and Evolution

### Code files to run analysis

#### 1. fit_model_to_simulated_data.R
Description: generate simulated data (using true parameters) and fit model 
to simulated data to obtain estimated parameters

#### 2. fit_model_to_tracking_data.R
Description: fit model to Northern gannet tracking data to obtain estimated parameters. See D2 below for data availability.

#### 3. estimate_HRs.R
Description: Estimate HRs from input parameters 
Input parameters obtained from step 1. (true/estimated simulated parameters) or from step 2. (estimated parameters from tracking data)

#### 4. validate_HRs_with_simulated_data.R
Description: Validate model with simulated data, by calculating similarity between predicted true and estimated HRs from step 3., and true and estimated exposure to wind farms. See D4 for data availability.

#### 5. validate_HRs_with_tracking_data.R
Description:  Validate model with out-of-sample tracking data, by apportioning GPS locations to colonies based on their predicted home ranges. 
See D2 for data availability 

#### 6. calculate_foraging_ranges.R 
Description: Calculate foraging ranges based on Grecian et al 2012 methods Steps 1-3
Grecian, W. J., Witt, M. J., Attrill, M. J., Bearhop, S., Godley,
B. J., Grémillet, D., Hamer, K. C., & Votier, S. C. (2012).
A novel projection technique to identify important at-sea areas for seabird
conservation: An example using Northern gannets breeding in the North East Atlantic.
Biological Conservation, 156, 43–52. https://doi.org/10.1016/j.biocon.2011.12.010

#### 7. compare_wind_farm_exposure.R
Description: Estimate exposure of Northern gannets to planned wind farms in 2020 using model-derived HRs, projected distributions and foraging ranges. 
Using estimated HRs predicted from step 3. (using parameters optimised to tracking data in step 2.), foraging ranges (produced in step 6) and projected distributions (code available from D3, see below). See D4 for wind farm data availability. 

### Data 

#### D1. land polygons can be downloaded via:
https://www.naturalearthdata.com/downloads/10m-physical-vectors/

#### D2. Northern gannet tracking data openly available via: 
Seabird Tracking Database: https://www.seabirdtracking.org/
Datasets: Grassholm (ID: 731,732) and Great Saltee (ID:723),
Ailsa Craig (ID: 716), Bass Rock (ID:718), Sule Skerry (ID: 719), Bull Rock (ID: 720), Lambay (ID: 724), Little Skellig (ID: 721), Les Etacs (ID:733), Ile Rouzic (ID:734)

#### D3. Code for projected distributions available via: 
Critchley, E. J., Grecian, W. J., Kane, A., Jessopp, M. J., & Quinn, J. L. (2018). 
Marine protected areas show low overlap with projected distributions of seabird populations in Britain and Ireland. Biological Conservation, 224, 309–317. 
https://doi.org/10.1016/j.biocon.2018.06.007

#### D4. Offshore wind farm data  
Data set version used in this analysis stored in Zenodo repository https://doi.org/10.5281/zenodo.14918040, 
file name "Wind_farm_data_EMODnet_20231124.zip"

Data obtained from EMODnet (accessed 24/11/2023):
EMODNet, human activities, energy, wind farms. (2014). 
https://emodnet.ec.europa.eu/geonetwork/srv/eng/catalog.search#/metadata/8201070b-4b0b-4d54-8910-abcea5dce57f

#### D5. Gannet colony size estimations and locations
Northern gannet colony data available via 10.5281/zenodo.10807638 (Jeglinski et al., 2024), file name "posterior_run69.csv".

Jeglinski, J. W. E., Niven, H. I., Wanless, S., Barrett, R. T., Harris, M. P., Dierschke, J., & Matthiopoulos, J. (2024). 
Data from: Past and future effects of climate on the metapopulation dynamics of a NorthEast Atlantic seabird across two centuries. 
https://doi.org/10.5281/zenodo.10807638

### Function files 
List of function files 

#### process_posterior_run69.R 
Description: Function to change column names of colony data frame

#### select_colony_year.R
Description: A function to filter a dataframe of colony sizes and locations for a specific colony year

#### create_colonies_sf.R
Description: A function to create a shapefile of colonies

#### create_bbox.R
Description: A function to create a bounding box encompassing all colonies' max foraging ranges.

#### create_landseamask.R 
Description: A function to create a land sea mask of specified resolution.

#### create_colonies_stars_array.R
Description:  A function to create a colonies stars raster of specified km resolution, with each colony size at colony centre in a separate layer.
Using landseamask (created from 'create_landseamask') as template

#### move_colonies2.R
Description: A function to move colonies on land to nearest sea cell

#### create_p_layer.R
Description: A function to create a permeability layer from a given mask specifying impermeable areas and permeable areas

#### calculate_min_dist.R
Description: A function to calculate shortest distances from colonies to all permeable cells

#### calculate_foraging.R
Description: A function to calculate the foraging surface

#### calculate_commute2.R
Description: A function to calculate commuting surfaces from foraging surfaces (from 'calculate_foraging')

#### normalise_foraging_commuting.R
Description: A function to normalise foraging and cummuting surfaces to unit sum. 

#### pre_optim.R
Description: A function to run to create objects needed before running model fitting

### Acknowledgements: 
Gannet tracking data were kindly contributed by: K. Hamer (Ailsa Craig), K. Hamer, R. Davies (Bass Rock), K. Hamer, J. Blackburn (Sule Skerry), S. Bearhop, T. Bodey (Bull Rock, Great Saltee, Lambay, Little Skellig), S. Votier, G. Morgan, L. Morgan (Grassholm), L. Soanes, J. Green (Les Etacs), D. Grémillet (Ile Rouzic). We thank EMODnet – Human Activities and primary sources for providing wind farm data: The Wind Power, OSPAR Commission (ODIM OSPAR data and information Management System), HELCOM - Baltic Marine Environment Protection Commission - Helsinki Commission, Wind Europe, Royal Belgian Institute of Natural Sciences, Danish Maritime Authority (Secretariat for maritime spatial planning), Ministry of Finance (Planning Department) (Rahandusministeerium), Ministry of The Environment and Regionals councils (Uusimaa, Kymenlaakso, Southwest Finland, Satakunta, Ostrobothnia, Central Ostrobothnia, North Ostrobothnia and Lapland), Åland Provincial Government, Cerema, Ministère de la Transition Écologique (GéoLittoral), Bundesamt für Seeschifffahrt und Hydrographie (BSH), Northland Deutsche Bucht GmbH, Marine Institute, Ministry of Environmental Protection and Regional Development of The Republic of Latvia, EPSOG, Ministry of Energy of the Republic of Lithuania, Lithuanian Energy Agency, Ministry of Infrastructure and the Environment (Noordzeeloket), Rijkswaterstaat - Ministry of Infrastructure and Water Management, Ministry of Maritime Economy and Inland Navigation, Maritime offices of Gdynia, Slupsk and Szczecin, Direção-Geral de Recursos Naturais, Segurança e Serviços Marítimos (DGRM), The Crown Estate UK, Crown Estate Scotland, Government of Spain - Ministry for Ecological Transition and the Demographic Challenge, Swedish Agency for Marine and Water Management. This work was funded by the UK Department for Energy Security & Net Zero DESNZ Offshore Energy Strategic Environment Assessment OESEA programme. 



