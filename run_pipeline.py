from scripts import preproccessing as pp
from scripts import filter
from scripts import unique
from scripts import waterfall
import yaml


with open('config.yaml', 'r') as file:
    config = yaml.safe_load(file)

#pp.main(config)
filter.main(config)

unique.main(config)

pipe = waterfall.WaterfallPipeline()
pipe.run()