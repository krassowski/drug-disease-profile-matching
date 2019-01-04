#!/usr/bin/env bash
sudo apt-get install libhdf5-serial-dev r-base

# pip drequirements which require compilation
python3-dev

# JupyterLab extensions
sudo apt-get install nodejs npm

# DrugCentral database
sudo apt-get install postgresql postgresql-contrib

# Broad Instutute runs on Java and specifically needs Oracle Java:
sudo add-apt-repository ppa:webupd8team/java	# webupd8team is not an official source, but fairly trusted by the community; this however might change in the future
sudo apt-get install oracle-java8-installer
