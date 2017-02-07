include paths.mk

prefix := `pwd`
subst_dir := 's|XXX_DIR_XXX|'${prefix}'|g'

met_driver := ed-inputs/met3/US-WCr/ED_MET_DRIVER_HEADER
ed2in_temp := run-ed/template/ED2IN
ed2_link := run-ed/template/ed_2.1

cohorts := 2cohort

denss := 0.015

dbhs := 20 30 40

site := US-WCr

pfts := "temperate.North_Mid_Hardwood temperate.Late_Hardwood"

stand_type := MLH

testsites := $(foreach d, $(dbhs), ed-inputs/sites/$(site)/rtm/$(cohorts)/dens$(denss)/dbh$d/$(stand_type)/$(stand_type).lat45.5lon-90.5.css)

results := $(foreach d, $(dbhs), run-ed/$(cohorts)/dens$(denss)/dbh$d/$(stand_type)/outputs/history.xml)

.PHONY: sites edruns

all: sites edruns templates inversion/edr_path

sites: $(testsites)

edruns: $(results)

templates : $(met_driver) $(ed2in_temp) $(ed2_link)

inversion/edr_path: 
	@echo $(EDR_EXE) > $@

$(testsites) : templates

$(results) : sites

%.css: 
	$(eval dt := $(shell expr match "$@" '.*dbh\([0-9]\+\).*'))
	$(eval pt := $(shell expr match "$@" '.*/\(.*\).lat.*'))
	$(eval st := $(shell expr match "$@" '.*dens\([0-9.]\+\).*'))
	Rscript generate_testrun_two_cohort.R $(dt) $(pfts) $(st) $(stand_type)

%.xml: 
	$(eval dt := $(shell expr match "$@" '.*dbh\([0-9]\+\).*'))
	$(eval pt := $(shell expr match "$@" '.*/dens.*/dbh.*/\(.*\)/outputs/.*'))
	$(eval st := $(shell expr match "$@" '.*dens\([0-9.]\+\).*'))
	./exec_ed_test.sh $(cohorts) $(dt) $(stand_type) $(st)

$(ed2_link): 
	ln -fs $(ED_EXE) $@

clean:
	rm -rf ed-inputs/sites/US-WCr/rtm/1cohort ed-inputs/sites/$(site)/rtm/2cohort/ \
	    run-ed/1cohort run-ed/2cohort run-ed/template/ed_2.1 \
	run-ed/template/ED2IN ed-inputs/met3/US-WCr/ED_MET_DRIVER_HEADER

%: %.temp
	sed $(subst_dir) $< > $@
