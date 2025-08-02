This repository contains the scripts used in Ravinandrasana and Franzke (2025, Nature Communications).

Each script and its purpose are described below:

1_SPEI48_computation_and_bias_correction.R : 
Computes the Standardized Precipitation-Evapotranspiration Index (SPEI) over a 48-month period and applies bias correction to the SPEI48 series. SPEI48 is a drought index that captures the cumulative effect of precipitation deficits and potential evapotranspiration over long timescales, indicating prolonged dry conditions.

2_SRFI48_computation_and_bias_correction.R : 
Computes the Standardized River Flow Index (SRFI) over a 48-month period and applies bias correction to the SRFI48 series. SRFI48 reflects long-term variations in river discharge, providing insights into hydrological drought over multi-year timescales.

3_backward_waterdemand_projection_1850_2010.R : 
Performs a backward regression to estimate water consumption from 1850 to 2009.

4_SWSI48_computation_and_bias_correction.R : 
Computes the Standardized Water Scarcity Index (SWSI) over a 48-month period and applies bias correction. SWSI48 is used to investigate the degree of water scarcity, which is conceptually equivalent to the Supply and Demand Balance Index.

5_DZD_event_detection.R : 
Detects Day Zero Drought (DZD) events based on the compound exceedance of multiple drought indices.

6_ToFE_DZD_event.R : 
Assesses the Time of First Emergence (ToFE) of DZD events across space and ensemble members.

7_Duration_WaitingTime_DZD.R : 
Evaluates the duration of DZD events and the waiting time between two consecutive DZD occurrences.

8_Plot_ToFE_map.py : used to generate the Time of First Emergence (ToFE) spatial maps presented in the manuscript. All other spatial maps shown in the study are based on the same plotting modules and structure used in this script. It serves as a representative example of the visualization approach applied throughout the analysis.

Note: 
For any questions and support if needed, the user is invited to contct christian.franzke@pusan.ac.kr or rphynodocilevecchia@pusan.ac.kr.
