# NOAA

This function includes all data necessary to complete a genome-wide association mapping using weather/climate information downloaded from the NOAA database as the phenotype as seen in Evans et al. (2016). 

##Workflow

### Input

#### Wild C. elegans Isolation Locations

This script uses the most current record of 152 wild isotypes (Cook et al. (2016)) curated by the Andersen Lab. Dataframe must include:

#### `isotype`

The unique strain/identifier

#### `latitude`

The latitude at which each strain is found

#### `longitude`

The longitude at which each strain is found

#### `isolation_date`

The isolation date when each strain was isolated in YYYY-MM-DD format. No data will be gathered for strains without a date of isolation.

#### `isolation_date_comment`

String describing the isolation location or date. In some cases, only the year of isolation is known, but the isolation_date is defaulted to January 1 of that year. The phrase "only year" is important to distinguish these cases. Similarly, sometimes only the month/year (but not exact date) of isolation is known, in which the phrase "only month/year" is necessary.

### noaa_mappings() function

This function will gather NOAA weather/climate information given a dataframe of wild isolates and a time period to collect information for. It will calculate the average maximum, minimum, and mean during that time period and returns a new dataframe. Relys on input dataframe to have comments for isolation date, should probably be manually checked before used.

#### `df`

Dataframe of wild isolates (described above)

#### `time_period`

Time period (integer) to collect weather data for, in months (e.g. If date of isolation is August 1, 2011 and time_period is 3 months, data will be collected from June 15, 2011 to September 15, 2011). Default is 3 months.

#### `important_trait`

String describing a specific weather trait that is marked as "important", meaning the weather station must have data for this variable in order to be used for data collection. Default is NULL.

#### `events`

Numerical value (default = 10) describing number of data points collected by the weather station during a given month.

#### `additional_data`

Additional data describes optional weather variables to download besides the primary. Takes more time to download more data, "AA1" is recommended, as it is the flag for precipitation. Default is FALSE (no additional data is downloaded). Options include TRUE (all quantitative weather traits) or you can manually insert any flags you wish. See stationaRy and NOAA format for more information

```{r}
wi_weather_1yr <- noaa_mappings(df = wi, time_period = 12, important_trait = "temp", additional_data = "AA1")
```
