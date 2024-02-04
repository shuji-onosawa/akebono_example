## python version
- 3.10.13

## How to set up the environment
- make vertual environment. The name of the environment is "venv". If you want to change the name, please change the name in ./onosawa/mca_ave/Makefile and run_bnd.sh.
At Make_mca_cdf directory, run the following command.
```
$ python3 -m venv .venv
```
- activate the environment
```
$ source .venv/bin/activate
```
- install python libraries
```
$ python3 -m pip install -r requirements.txt
```

## akebono package
- Please read ./akebono/README.md
- Check the local_data_dir in ./akebono/config.py. If you want to change the directory which contains mca data cdf files, please change local_data_dir in config.py. If you want to store mca data cdf files in this repository, you should have better to ignore the files in .gitignore.