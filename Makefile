all: run test
run:
	python scripts/run_all.py
test:
	python -m pytest -q tests
anova:
	python at_mmc_anova.py
