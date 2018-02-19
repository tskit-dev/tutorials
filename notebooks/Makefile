all: md

MD=wfforward.md

md: $(MD) $(OUTPUT)

clean:
	for i in $(MD) ; do \
		rm -f ../docs/$$i ; \
	done
	for i in $(basename $(MD)) ; do \
		rm -fr ../docs/$$i"_files" ; \
	done
	rm -f *.output.ipynb
	rm -rf *.output_files

%.output.ipynb: %.ipynb
	python execute.py $< $@

%.md: %.output.ipynb
	# Set --output so that supporting
	# file folder (foo_files) gets
	# generated with correct 
	# output name and location
	jupyter nbconvert --to markdown --output-dir=../docs --output=$@ $<

.PHONY: md
.PRECIOUS: %.output.ipynb
