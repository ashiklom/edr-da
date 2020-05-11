.PHONY: all clean drake pdf sync

all: pdf

pdf: drake
	cd text && make

drake:
	RUN_CONFIG=homo-pooled Rscript -e "drake::r_make()"
	RUN_CONFIG=hetero-pooled Rscript -e "drake::r_make()"
	RUN_CONFIG=homo-sitespecific Rscript -e "drake::r_make()"
	RUN_CONFIG=hetero-sitespecific Rscript -e "drake::r_make()"

sync:
	rsync -avz --progress discover:~/projects/phd/edr-da/multi_site_pda_results/ multi_site_pda_results/

clean:
	cd text && make clean
