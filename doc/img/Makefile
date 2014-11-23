all: models.png models-legend.png

%.png: %.dot
	dot -Tpng $^ -o $@

