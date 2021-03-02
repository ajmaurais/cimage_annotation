# cimage_annotation
Functional cysteine annotation

# Usage
There are two executable scripts:

* `cimage_annotation`: Run the annotation program in your current shell.
* `qsub_cimage_annotation`: Automatically submit the annotation program as a `PBS` job.

```
usage: cimage_annotation [-h] [-f {cimage,dtaselect}] [-s] [--ofname OFNAME]
                         [-a {0,1}] [-w {0,1}] [-d DATABASE_DIR]
                         [-o DEFINED_ORGANISM] [-p {0,1}] [-t NTHREAD] [-v]
                         input_file

usage: qsub_cimage_annotation [-h] [-f {cimage,dtaselect}] [-s {0,1}]
                              [--ofname OFNAME] [-a {0,1}] [-w {0,1}]
                              [-d DATABASE_DIR] [-o DEFINED_ORGANISM]
                              [-v {0,1}] [-m MEM] [-p PPN] [-t WALLTIME] [-g]
                              input_file
```

You can use the `-h` flag to display a detailed descriptions of all options.
```bash
cimage_annotation -h
qsub_cimage_annotation -h
```

# Examples

Run `cimage_annotation` and retrieve functional cysteine annotations from UniProt.
```bash
cimage_annotation <input_file>
```

Optionally, you can also align cysteine sites to other organisms to determine whether the cysteine is evolutionarily conserved with the `--align` option. Because the sequence alignments are computationally intensive, when the `--align` option is specified, you should submit the command as a PBS job. `cimage_annotation` includes an additional command to automatically generate a `.pbs` file and submit it to the queue for you.
```bash
qsub_cimage_annotation --align --database_dir <path_to_dir_with_sequence_databases> -g <input_file>
```

# How to install on Sirius

1\. First clone the `cimage_annotation` GitHub repository. You can store the `cimage_annotation` source code anywhere on your `sirius` account. However, the installation instructions assume the program is being installed in `~/code`. 

```bash
mkdir -p ~/code # make ~/code if it doesn't already exist
cd ~/code # navigate to ~/code
git clone https://github.com/ajmaurais/cimage_annotation # clone the GitHub repository
```

2\. Next, setup a python virtual environment in which to install `cimage_annotation` and its dependencies. `cimage_annotation` is written for `python >= 3.6.*` so you will first need to load the `python3` module.

```bash
cd ~/code/cimage_annotation # navigate to the cimage_annotation directory
module load python/3.6.3 # load the module
python3 -m venv venv # set up the python virtual environment
```

After you have run the above commands, the directory should look like this:
```bash
$ ls
README.md  setup.py  src/  venv/
```
If you are missing the `venv` directory, something went wrong and the following commands will not work.

3\. Next, build and install `cimage_annotation` in your virtual environment.

```bash
cd ~/code/cimage_annotation
./venv/bin/python setup.py build
./venv/bin/pip install .
```

4\. Add `cimage_annotation` to your `$PATH` (optional)

```bash
mkdir -p ~/bin
cp --remove-destination ~/code/cimage_annotation/venv/bin/cimage_annotation ~/bin
cp --remove-destination ~/code/cimage_annotation/venv/bin/qsub_cimage_annotation ~/bin
```

5\. Log out and log back in to apply the changes to your shell environment. Once you log back in, check that the `cimage_annotation` command is recognized by your shell. Running the command `which cimage_annotation` should output: `~/bin/cimage_annotation`
