include paths.mk

prefix := `pwd`
subst_dir := 's|XXX_DIR_XXX|'${prefix}'|g'

met_driver := ed-inputs/met3/US-WCr/ED_MET_DRIVER_HEADER
ed2in_temp := run-ed/template/ED2IN
templates := ${met_driver} ${ed2in_temp}

cohorts := 1cohort

denss := 0.05

dbhs := 20 30 40

pfts := temperate.Early_Hardwood \
	temperate.Mid_Hardwood \
	temperate.Late_Hardwood \
	temperate.North_Pine \
	temperate.Late_Conifer

testsites := $(foreach c, $(cohorts), \
	$(foreach s, $(denss), \
	$(foreach d, $(dbhs), \
	$(foreach p, $(pfts), \
	ed-inputs/sites/US-WCr/rtm/$c/dens$s/dbh$d/$p/$p.lat45.5lon-90.5.css))))

results := $(foreach c, $(cohorts), \
	$(foreach s, $(denss), \
	$(foreach d, $(dbhs), \
	$(foreach p, $(pfts), \
	ed-inputs/sites/US-WCr/rtm/$c/dens$s/dbh$d/$p/outputs/history.xml))))

edlinks := $(foreach c, $(cohorts), \
	$(foreach s, $(denss), \
	$(foreach d, $(dbhs), \
	$(foreach p, $(pfts), \
	ed-inputs/sites/US-WCr/rtm/$c/dens$s/dbh$d/$p/ed_2.1))))

all: $(testsites)

%.css:
	$(eval dt := $(shell expr match "$@" '.*dbh\([0-9]\+\).*'))
	$(eval pt := $(shell expr match "$@" '.*/\(.*\).lat.*'))
	$(eval st := $(shell expr match "$@" '.*dens\([0-9.]\+\).*'))
	Rscript generate_test.R $(dt) $(pt) $(st)

%/ed_2.1:
	ln -sf $(ED_EXE) $@

%/outputs/history.xml: %/ED2IN %/ed_2.1
	$(eval dt := $(shell expr match "$@" '.*dbh\([0-9]\+\).*'))
	$(eval pt := $(shell expr match "$@" '.*/\(.*\).lat.*'))
	$(eval st := $(shell expr match "$@" '.*dens\([0-9.]\+\).*'))
	exec_ed_test.sh dt pt st

clean:
	rm -rf ed-inputs/sites/US-WCr/rtm/1cohort \
	    run-ed/1cohort/*/*/outputs/*

%: %.temp
	sed $(subst_dir) $< > $@

#results: testsites
	#./exec_ed_tests.sh
