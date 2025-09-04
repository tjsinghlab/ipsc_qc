.PHONY: clean install data requirements

# init

SHELL=/bin/bash
PROJECT_NAME=wdlplay
CONDA_BASE:=$(shell conda info --base)

# install

update-pip:
	pip install -U pip setuptools wheel

install: clean update-pip
	pip install -e .

install-setup: clean update-pip
	python setup.py install

install-reqs:
	pip install --upgrade --no-deps --force-reinstall -r requirements.txt

install-dev:
	pip install --upgrade --no-deps --force-reinstall -r requirements-dev.txt

freeze-reqs: install-requirements
	pip freeze > requirements.txt

version:
	@python -c 'import wdlplay; print(wdlplay.__version__)'

lib-version:
	@python -c 'import $(NAME); print($(NAME).__version__)'

hail-version:
	@python -c 'import hail; print(hail.version())'

# clean

clean: clean-py

clean-py: clean-build clean-pyc clean-test

clean-pyc:
	find . -name '*.py[co]' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

# env

create-env:
	conda create -y --name $(PROJECT_NAME) python=3.11.4 jupyterlab=3.6.3 pylint versioneer

set-env:
	conda env config vars set CLOUD_ROOT=gs://$(PROJECT_NAME) PROJECT_ROOT=`pwd` 

remove-env: clean
	( source $(CONDA_BASE)/etc/profile.d/conda.sh && \
		conda deactivate && \
		conda remove --name $(PROJECT_NAME) --all )

install-version:
	versioneer install


# wdl cromwell

CROMWELL_VERSION=78
install-cromwell:
	mkdir -p bin/
	wget -O bin/cromwell.jar https://github.com/broadinstitute/cromwell/releases/download/$(CROMWELL_VERSION)/cromwell-$(CROMWELL_VERSION).jar
	wget -O bin/womtool.jar https://github.com/broadinstitute/cromwell/releases/download/$(CROMWELL_VERSION)/womtool-$(CROMWELL_VERSION).jar

run-mini:
	miniwdl run wdl/$(NAME)/$(NAME).wdl --input wdl/$(NAME)/$(NAME).json --dir data/$(NAME)/$(VERSION)

check:
	miniwdl check wdl/$(NAME)/$(NAME).wdl

inputs:
	java -jar bin/womtool.jar inputs wdl/$(NAME)/$(NAME).wdl > $(NAME).json

zip:
	miniwdl zip wdl/$(NAME)/$(NAME).wdl --output zip/$(NAME).$(VERSION).wdl.zip

help-miniwdl:
	miniwdl -h

help-cromwell:
	java -jar bin/cromwell.jar --help

help-womtool:
	java -jar bin/womtool.jar --help

# caper

CLOUD_PROJECT=daly-neale-sczmeta
REGION=us-central1
ZONE=us-central1-a
VM=$(PROJECT_NAME)-cromwell
KEY='~/.service_account_key/daly-neale-sczmeta-5c636e7cb4e5.json'
PORTR=8000

istart:
	bash create_instance.sh $(VM) $(CLOUD_PROJECT) $(KEY) gs://$(PROJECT_NAME)/data/wdl/ --gcp-region $(REGION) --zone $(ZONE) --gcp-loc-dir gs://$(PROJECT_NAME)/data/wdl/.caper_tmp/

issh:
	gcloud compute ssh $(VM) --zone $(ZONE) -- -L $(PORTR):localhost:$(PORTR)

idelete:
	gcloud compute instances delete $(VM) --zone $(ZONE)

clean-wdl:
	find . -name 'cromwell.out*' -exec rm -f {} +
	find . -name 'hail-*.log' -exec rm -f {} +

# run:
# 	caper run wdl/$(NAME)/$(NAME).wdl --inputs wdl/$(NAME)/$(NAME).json --local-out-dir data/$(NAME)/0.1 --local-loc-dir tmp/

# submit:
# 	caper submit wdl/$(NAME)/$(NAME).wdl --inputs wdl/$(NAME)/$(NAME).json -o wdl/$(NAME)/options.json --dry-run

# build

build: clean-py
	python setup.py sdist bdist_wheel
	twine check dist/*

upload-test:
	twine upload --repository testpypi dist/*

upload:
	twine upload dist/*

#----
view:
	open https://console.cloud.google.com/storage/browser/$(PROJECT_NAME)/data

view-compute:
	open https://console.cloud.google.com/compute/instances

view-dataproc:
	open https://console.cloud.google.com/dataproc/clusters

view-jobs:
	open https://console.cloud.google.com/dataproc/jobs