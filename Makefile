.PHONY: all clean drake pdf sync

all: drake pdf

pdf:
	cd text && make

drake:
	RUN_CONFIG=hetero-pooled-exp Rscript -e "drake::r_make()"
	RUN_CONFIG=hetero-pooled-lnorm Rscript -e "drake::r_make()"
	RUN_CONFIG=homo-pooled-exp Rscript -e "drake::r_make()"
	RUN_CONFIG=homo-pooled-lnorm Rscript -e "drake::r_make()"
	RUN_CONFIG=homo-sitespecific-exp Rscript -e "drake::r_make()"
	RUN_CONFIG=homo-sitespecific-lnorm Rscript -e "drake::r_make()"
	RUN_CONFIG=hetero-sitespecific-exp Rscript -e "drake::r_make()"
	RUN_CONFIG=hetero-sitespecific-lnorm Rscript -e "drake::r_make()"

sync:
	rsync -avz --progress discover:~/projects/phd/edr-da/multi_site_pda_results/ multi_site_pda_results/

clean:
	cd text && make clean
