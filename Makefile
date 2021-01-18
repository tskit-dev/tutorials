
# Simple makefile for dev.

all:
	# Catch warnings for development
	jupyter-book build -W -n --keep-going .

clean:
	rm -fR _build
