## README ##

[![Build Status](https://travis-ci.org/EBIvariation/eva-cttv-pipeline.svg?branch=master)](https://travis-ci.org/EBIvariation/eva-cttv-pipeline)
[![Coverage Status](https://coveralls.io/repos/github/EBIvariation/eva-cttv-pipeline/badge.svg?branch=master)](https://coveralls.io/github/EBIvariation/eva-cttv-pipeline?branch=master)


Python version required: 3.5

Installing python 3.5 (optional)
-------
In case you have a newer version and need to run the pipeline you can follow the next steps to install python 3.5 without replacing you default python.
1. Download python 3.5 "wget https://www.python.org/ftp/python/3.5.5/Python-3.5.5.tgz"
2. "sudo tar xzf Python-3.5.5.tgz"
3. "cd Python-3.5.5"
4. "sudo apt-get install zlib1g-dev"
5. "sudo ./configure --enable-optimizations --with-zlib=/usr/include"
6. "sudo make altinstall"
7. "python3.5 -V"

Building java Clinvar parser
-------
1. "cd clinvar-xml-parser"
2. "mvn package"
3. Two jar files will be generated in the 'target' directory, one of them including all the dependencies.  


Building python pipeline and (optional) Setting up virtual environment
-------

1. "git clone --recursive git@github.com:EBIvariation/eva-cttv-pipeline.git"
2. "cd eva-cttv-pipeline"
3. [OPTIONAL] "virtualenv -p python3.5 venv"
4. [OPTIONAL] "source venv/bin/activate" ("venv/bin/deactivate" to deactivate virtualenv)
5. pip install -r requirements.txt
6. And then one of:
   * To install: "python3 setup.py install"
   * To install to develop: "python3 setup.py develop"
   * To build a source distribution: "python3 setup.py sdist"

Usage
-------

Please see the [GitHub wiki](https://github.com/EBIvariation/eva-cttv-pipeline/wiki/How-to-submit-an-OpenTargets-batch) for usage


Tests
-------

You can run all tests with: `python setup.py test`
