Hello everyone,

I've just spent the morning going through the GeMM model source code and having
a look at what would actually have to change to adapt it to Zosterops. I have to
say the list is not trivial, and in several cases I am unsure of what implementation 
choices make the most scientific/programming sense.

Basically, it boils down to the fact that GeMM was designed as a community-level plant 
model, whereas we now need a species-level bird model. Here's what I've found so far:

- The current plant model assumes multiple species in each tile's community. In the 
  Zosterops model, we pretty much want a binary situation (either highland or lowland 
  species breeding on a tile). This has implications for the initialisation and dispersal 
  submodels.

- To achieve this, we probably need to implement habitat types, a concept we do not yet 
  have in the model. (Rather than qualitative, discrete habitats such as cloud forest or 
  savannah, the current model only has quantitative gradients created by the local 
  combination of temperature and precipitation.)

- Currently, the number of offspring is determined metabolically from parent size
  and local temperature. Does this need to be species- and/or habitat-dependent?
    
- At the moment, tile carrying capacity is a global constant. Is that an acceptable 
  simplification? Or do we need to vary the maximum number of breeding pairs per tile, 
  depending on habitat type (and perhaps quality)?

- In the plant model, competition via individuals' adaptation to precipitation kicks 
  in when the tile's carrying capacity is reached. Precipitation is irrelevant for 
  our birds - so how do we do competition?

- The seed dispersal function needs to be rewritten. Instead of the current stochastic, 
  habitat-agnostic dispersal kernel, we need directed, habitat-dependent active dispersal.

- Pollen dispersal, too, needs to change from either local or global to regional (birds 
  search the surrounding tiles for possible mates).

- Our plant individuals are assumed to be hermaphroditic. Birds aren't.

- Instead of randomly initialising species pools, we need to read in configuration files 
  for predefined species.

In summary, the two biggest challenges concern the introduction of discrete habitats and 
the open questions regarding carrying capacity/competition. Overall, this is a significant 
intervention in large parts of the model's ecological side (the genetics are less affected), 
and I'm not sure we'll be able to pull it off while keeping backwards compatibility.

Jan, perhaps you could weigh in on the biology of some of these questions. And Ludwig, 
if you have any suggestions as to how to sensibly implement this, they'd be gratefully 
received ;-) And if I missed anything, do let me know.

All the best,
Daniel
