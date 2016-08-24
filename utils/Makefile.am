bin_SCRIPTS = shuffle create_lp_file affy_to_eset			\
get_absent_present_genes_marray get_absent_present_genes_rnaseq		\
extract_sbml_model_info get_entrezid_for_probesets get_jetset_scores	\
get_probeset_greatest_jscore exec_panp_eset filter_panp_results		\
auto_fba gen_fba_stats plot_fba_network plot_metatool_network		\
instantiate_fva_templ solve_fva_for_vlist auto_fva network_reducer	\
calc_matrix_rank get_cplex_fluxes plot_kde_rnaseq gather_sample_fluxes	\
test_samples extract_info_from_escher_json extract_sbml_mod_info	\
correct_pvalues show_data_for_column

CLEANFILES = $(bin_SCRIPTS)

EXTRA_DIST = shuffle.sh create_lp_file.py affy_to_eset.R		\
get_absent_present_genes_marray.py get_absent_present_genes_rnaseq.py	\
extract_sbml_model_info.R get_entrezid_for_probesets.R			\
get_jetset_scores.R get_probeset_greatest_jscore.sh exec_panp_eset.R	\
filter_panp_results.sh auto_fba.sh gen_fba_stats.sh fba.py		\
plot_fba_network.py plot_metatool_network.py instantiate_fva_templ.sh	\
solve_fva_for_vlist.sh auto_fva.sh network_reducer.sh			\
calc_matrix_rank.R get_cplex_fluxes.sh plot_kde_rnaseq.py		\
gather_sample_fluxes.py test_samples.py					\
extract_info_from_escher_json.py extract_sbml_mod_info.py		\
correct_pvalues.py show_data_for_column.py

python_PYTHON= fba.py

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
	echo "GZIP="$(GZIP) >> $@
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
	echo "import sys" >> $@
	echo "sys.path.append(\"$(pythondir)\")" >> $@	
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
