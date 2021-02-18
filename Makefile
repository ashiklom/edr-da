.PHONY: all clean drake pdf sync

all: drake pdf

pdf:
	cd text && make

drake:
	# RUN_CONFIG=revision Rscript -e "drake::r_make()"
	RUN_CONFIG=revision-fixed Rscript -e "drake::r_make()"

sync:
	rsync -avz --progress discover:~/projects/phd/edr-da/multi_site_pda_results/ multi_site_pda_results/

clean:
	cd text && make clean
