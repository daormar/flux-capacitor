bin_SCRIPTS = create_lp_file affy_to_eset get_absent_present_genes	\
extract_sbml_model_info get_entrezid_for_eset_genes get_jetset_scores	\
get_probeset_greatest_jscore exec_panp_eset filter_panp_results		\
auto_fba analyze_fba_results

CLEANFILES = $(bin_SCRIPTS)

EXTRA_DIST = create_lp_file.py affy_to_eset.R				\
get_absent_present_genes.py extract_sbml_model_info.R			\
get_jetset_scores.R get_probeset_greatest_jscore.sh exec_panp_eset.R	\
filter_panp_results.sh auto_fba.sh analyze_fba_results.sh

# Build using suffix rules (portable)
SUFFIXES= .sh .py .r .R

.sh :
	rm -f $@
	echo "#! " $(BASH) > $@
	echo "pkgname=\""$(PACKAGE_NAME)"\"" >> $@
	echo "version=\""$(PACKAGE_VERSION)"\"" >> $@
	echo "bindir="$(prefix)/bin/ >> $@
	echo "AWK="$(AWK) >> $@
	echo "SORT="$(SORT) >> $@
	echo "sortpars=\""$(sortpars)"\"" >> $@
	echo "sortT=\""$(sortT)"\"" >> $@
	echo "UNIQ="$(UNIQ) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "SPLIT="$(SPLIT) >> $@
	echo "HEAD="$(HEAD) >> $@
	echo "TAIL="$(TAIL) >> $@
	echo "FIND="$(FIND) >> $@
	echo "MKTEMP="$(MKTEMP) >> $@
	echo "QSUB=\""$(QSUB)"\"" >> $@
	echo "QSUB_WORKS=\""$(QSUB_WORKS)"\"" >> $@
	echo "QSUB_TERSE_OPT=\""$(QSUB_TERSE_OPT)"\"" >> $@
	echo "QSTAT="$(QSTAT) >> $@
	echo "QSTAT_J_OPT=\""$(QSTAT_J_OPT)"\"" >> $@
	echo "QRLS="$(QRLS) >> $@
	cat $< >> $@
	chmod ugo+x $@

.py:
	rm -f $@
	echo "#! " $(PYTHON) > $@
	cat $< >> $@
	chmod ugo+x $@

.r:
	rm -f $@
	echo "#! " $(RSCRIPT) > $@
	cat $< >> $@
	chmod ugo+x $@

.R:
	rm -f $@
	echo "#! " $(RSCRIPT) > $@
	cat $< >> $@
	chmod ugo+x $@
