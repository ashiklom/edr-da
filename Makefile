.PHONY: all clean drake pdf sync traceplots

all: drake pdf

pdf:
	cd text && make

traceplots:
	RUN_CONFIG=revision-fixed Rscript -e "drake::r_drake_build('traceplots')" -e "drake::r_drake_build('traceplots_full')"

drake:
	RUN_CONFIG=revision-fixed Rscript -e "drake::r_make()"

sync:
	rsync -avz --progress discover:~/projects/phd/edr-da/multi_site_pda_results/ multi_site_pda_results/

clean:
	cd text && make clean

response:
	cd text/comments-and-submissions/2021-03-05-response && make pdf
