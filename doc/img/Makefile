all: models.png models-legend.png algorithm.PNG

%.png: %.dot
	dot -Tpng $^ -o $@

algorithm.ps: algorithm.dot
	neato -Tps -l sdl.ps $< -o $@

algorithm.PNG: algorithm.ps
	convert -flatten $< $@
