
MD5_EXCLUDE ?= 

Git_Info.C: Git_Info.C.in
	@if ! which git > /dev/null || \
	  ! (cd $(top_srcdir); git rev-parse HEAD > /dev/null 2>&1); then \
	  if test -f $(srcdir)/$@; then \
	    cp $(srcdir)/$@ $@.tmp; chmod u+rw $@.tmp; \
	  else \
	    echo "***************************************************"; \
	    echo "* Incomplete sources and no git information. This *"; \
	    echo "* copy of Sherpa will not be supported. Please    *"; \
	    echo "* contact sherpa@projects.hepforge.org for help.  *"; \
	    echo "***************************************************"; \
	    exit 1; \
	  fi; \
	else \
	  rev=$$(cd $(top_srcdir); git rev-parse HEAD); \
	  if test -n "$$(cd $(srcdir); git status -s --untracked-files=no .)"; then \
	    rev=$$rev"-dirty"; \
	  fi; \
	  url=$$(git rev-parse --abbrev-ref HEAD); \
	  echo '#include "ATOOLS/Org/Git_Info.H"' > $@.tmp; \
	  echo 'static ATOOLS::Git_Info initializer' >> $@.tmp; \
	  echo '("$(GITTAG)","'$$url'","'$$rev'","X");' >> $@.tmp; \
	fi; \
	if test -z $(NOMD5SUM); then \
	  mds=$$(cat $(addprefix $(srcdir)/, \
	    $(filter-out $@ $(CONFIG_HEADER) $(MD5_EXCLUDE), \
	    $(SOURCES) $(HEADERS))) | $(MD5COMMAND)); \
	  $(SEDCOMMAND) -e's/".?"\);/"'$$mds'");/g' $@.tmp; \
	fi; \
	if ! diff $@.tmp $@ > /dev/null 2>&1; then \
	  mv $@.tmp $@; \
	else \
	  rm $@.tmp; \
	fi;

.PHONY: Git_Info.C.in

DISTCLEANFILES = Git_Info.C
