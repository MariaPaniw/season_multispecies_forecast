# season_multispecies_forecast
Simulations of seasonal stage-specific metapopulation models for 3 interacting species

The R scripts implement seasonal demographic rates of interacting species under various scenarios of demography-biotic-interaction feedbacks. Each demographic rate in each season can be modelled as a function of environmental states, intraspecific density, and interspecific density. The local dynamics of each population of each interacting species in a given habitat, described by a local matrix model, were connected by a dispersal matrix for a given stage. A hyperstate matrix modelling approach is then used to integrate the local dynamics and dispersal matrices into a metapopulation model to assess spatiotemporal changes in stage-specific abundances of each species.

![image](https://user-images.githubusercontent.com/12189756/185378620-e4903799-42ee-4b01-9646-0361bc3258f7.png)
