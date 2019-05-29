#######
### Calculating properly country-weighted Red List Indices based on proportional ranges of species
### created: 28 May 2019

# Libraries
library(red)

# Read in csv with data containing proportional ranges of each species within each country
country.props = read.csv('/Users/jessicacheok/Documents/Work/vmshare/wedge-guitarfish/Tables/country_props/all_countries_species_range_proportions.csv')
country.props = country.props[,2:4]  # eliminate row indexes from python
# Loop through reading in total distribution ranges exported for each species
out = NULL  # empty object for appending information to
