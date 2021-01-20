
# Simple makefile for dev.

all:
	# Use the local build wrapper to automate writing the report log to stdout.
	./build.sh

clean:
	rm -fR _build
