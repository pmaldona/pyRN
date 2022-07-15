# CRNS
Closed Reaction Network Structure Library

## Install

CRNS can be installed from source. After downloading the repository, and and afterwards run:

    $ sh setup.sh

CRNS require the installation of the lxml library, in Debian/Ubuntu can be installed by:

    $ sudo apt-get install libxml2-dev libxslt-dev python-dev

Please refer to https://lxml.de/installation.html for other systems.

## Use

run 
    $ npm run start

This will first trigger a rebuild of the frontend, the start the Electron application. Electron will then start the python backend, which will, when finished load the UI (loadtimes for pyRN is a bit long 1~5 seconds).
