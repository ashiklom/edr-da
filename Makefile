.PHONY: all clean drake pdf

all: pdf

pdf: drake
	cd text && make

drake:
	RUN_CONFIG=homo-pooled Rscript -e "drake::r_make()"
	RUN_CONFIG=hetero-pooled Rscript -e "drake::r_make()"
	RUN_CONFIG=homo-sitespecific Rscript -e "drake::r_make()"
	RUN_CONFIG=hetero-sitespecific Rscript -e "drake::r_make()"

clean:
	cd text && make clean
