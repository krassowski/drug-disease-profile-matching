#!/usr/bin/env bash
sudo apt-get install libhdf5-serial-dev

# R
sudo apt-get install r-base
echo "If the version of R in repos is too old, follow https://stackoverflow.com/a/10476798/6646912 to upgrade"

# For GSEA-Base required by GSVA
sudo apt-get install libcurl4-gnutls-dev
sudo apt-get install libxml2-dev

# pip drequirements which require compilation
sudo apt-get install python3-dev

# JupyterLab extensions
sudo apt-get install nodejs npm

# DrugCentral database
sudo apt-get install postgresql postgresql-contrib

# Broad Instutute runs on Java and specifically needs Oracle Java:
sudo add-apt-repository ppa:webupd8team/java	# webupd8team is not an official source, but fairly trusted by the community; this however might change in the future
sudo apt-get install oracle-java8-installer
