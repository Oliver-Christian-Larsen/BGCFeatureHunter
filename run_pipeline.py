from scripts import preproccessing as pp
from scripts import filter
from scripts import unique
from scripts import waterfall

pp.main()
filter.main()

unique.main()

pipe = waterfall.WaterfallPipeline()
pipe.run()