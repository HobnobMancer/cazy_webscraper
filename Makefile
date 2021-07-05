# Makefile for cazy_webscraper
#
# This file is part of the cazy_webscraper package distribution
# (https://github.com/HobnobMancer/cazy_webscraper)

# Set up all development dependencies in the current conda environment
setup_env:
	@conda install --file requirements-dev.txt --yes
	@conda install --file requirements.txt --yes
	@pip install -r requirements-pip.txt
	@pip install -U -e .

# Clean up local Sphinx documentation output
clean_docs:
<<<<<<< HEAD
	@rm -rf docs/build/html

# Build and display local Sphinx documentation
docs: clean_docs
	@cd docs && make html && open build/html/index.html
=======
	@rm -rf docs/_build/html

# Build and display local Sphinx documentation
docs: clean_docs
	@cd docs && make html && open _build/html/index.html
>>>>>>> d2a1aac (add Makefile to automate documentation building)
