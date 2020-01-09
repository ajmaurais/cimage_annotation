# cimage_annotation
Functional cysteine annotation

# Usage
```
usage: cimage_annotation [-h] [-f {cimage,dtaselect}] [-s] [--ofname OFNAME]
                         [-a {0,1}] [-w {0,1}] [-d DATABASE_DIR]
                         [-o DEFINED_ORGANISM] [-p {0,1}] [-t NTHREAD] [-v]
                         input_file
```

# Examples

Run `cimage_annotation` and retrieve functional cysteine annotations from UniProt.
```bash
cimage_annotation <input_file>
```

Optionally, you can also align cysteine sites to other organisms do determine evolutionary conservation with the `--align 1` option. Because the sequence alignments are computationally intensive, when the `--align 1` is specified, you should submit the command as a PBS job. `cimage_annotation` includes an additional command to automatically generate a `.pbs` file and submit it to the queue for you.
```bash
qsub_cimage_annotation --align 1 --database_dir <path_to_dir_with_sequence_databases> -g <input_file>
```

# How to install on Pleiades

1\. First clone the `cimage_annotation` GitHub repo. You can store the `cimage_annotation` source code anywhere on your `pleiades` account, but I recommend putting locally installed programs in `~/local`. 

```bash
mkdir -p ~/local & cd ~/local # make ~/local if it doesn't exist and navigate to it
git clone https://github.com/ajmaurais/cimage_annotation # clone the GitHub repository
```

2\. Next, setup a python virtual environment in which to install `cimage_annotation` and its dependencies. `cimage_annotation` is written for `python >= 3.6.*` so you will first need to load the `python3` module.

```bash
cd ~/local/cimage_annotation # navigate to the cimage_annotation directory
module load python/3.6.1 # load the module
python3 -m venv venv # set up the python virtual environment
```

After you have run the above commands, the directory should look like this:
```bash
$ ls
README.md  setup.py  src/  venv/
```
If you are missing the `venv` directory, something went wrong and the following commands will not work.

3\. There is a bug in `biopython` 1.6.1, one of the dependencies for `cimage_annotation`, which will cause an error if it is used. Therefore, you must manually install `biopython` from a fork of the main `biopython` GitHub repository, in which the bug has been fixed. The build and install steps will produce a lot of output and take several seconds to finish.

```bash
cd ~/local
git clone --branch swissprot_bugfix https://github.com/ajmaurais/biopython # Clone the biopython repo
cd biopython # Navigate to the biopython directory
~/local/cimage_annotation/venv/bin/python setup.py build # Run the biopython build script
~/local/cimage_annotation/venv/bin/pip install . # Install biopython in your virtual environment
```

4\. Next, build and install `cimage_annotation` in your virtual environment.

```bash
cd ~/local/cimage_annotation
./venv/bin/python setup.py build
./venv/bin/pip install .
```

5\. Add `cimage_annotation` to your `$PATH` (optional)

```bash
mkdir -p ~/bin
ln -s ~/local/cimage_annotation/venv/bin/cimage_annotation ~/bin
```

6\. Log out and log back in to apply the changes to your shell environment. Once you log back in, check that the `cimage_annotation` command is recognized by your shell. Running the command `which cimage_annotation` should output: `~/bin/cimage_annotation`
