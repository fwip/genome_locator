#!/bin/bash

set -e

sudo yum install -y git
sudo amazon-linux-extras install python3
sudo pip3 install pipenv

git clone https://github.com/fwip/genome_service
cd genome_service
git pull
pipenv install

FLASK_APP=run.py pipenv run flask run
