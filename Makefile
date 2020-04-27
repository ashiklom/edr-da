.PHONY: all clean drake pdf

all: pdf

pdf: drake
	cd text && make

drake:
	Rscript -e "drake::r_make()"

clean:
	cd text && make clean
