MODULES = align2seq
DATA_built = align2seq.sql align2seq_uninstall.sql
DOCS = README.align2seq

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
