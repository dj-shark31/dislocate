# dislocate
Python tools for analyzing dislocation core structures.

## Dependencies

### Babel
This package requires Babel for dislocation calculations. You will need to:
1. Download and build Babel from source: http://emmanuel.clouet.free.fr/Programs/Babel
2. Provide path for executable in core_analysis/get_data.py and core_search/ESPM.py 

### OVITO (Modified Version)
A modified version of OVITO written by Eric Rothchild (source code available at https://github.com/Erothch/ovitoWallace) is used to compute the elastic stability parameter, the dxa and common neighbor analysis. 
For the convenience of the user, the modified OVITO is compiled in a Docker container dj31/ovito-wallace:amd64. Docker is used to run the ovito scripts if supported by the user, otherwise a singularity container can be used. 
Run this command to download the sif file: `singularity pull ovitowallace.sif docker://dj31/ovito-wallace:amd64`.

Note: The dxa and common neighbor analysis can be runned with a non-modified version of ovito.

## Configuration

To run scripts that depend on external tools (e.g. Ovito, Babel), create a `config.yaml` file:

```yaml
tools:
  ovitosif: /path/to/ovito.sif
  babel: /opt/babel/babel
  displacement: /opt/babel/displacement
  patternDetect: /opt/babel/patternDetect
  patternInit: /opt/babel/patternInit
