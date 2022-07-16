# Wogan-2022-Impacts

Following commands download all code needed to reproduce results. Must have a Fortran and C compiler.

```sh
# create conda environment
conda create -n postimp -c conda-forge python=3.9 numpy scipy pyyaml scikit-build cython jupyter cantera matplotlib pathos

# activate
conda activate postimp 

# photochem v0.2.10
wget https://github.com/Nicholaswogan/photochem/archive/cc54cb5a5d637ddf637edc962a8044b2be39dc5d.zip
unzip cc54cb5a5d637ddf637edc962a8044b2be39dc5d.zip
cd photochem-cc54cb5a5d637ddf637edc962a8044b2be39dc5d
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf photochem-cc54cb5a5d637ddf637edc962a8044b2be39dc5d cc54cb5a5d637ddf637edc962a8044b2be39dc5d.zip

# ImpactAtmosphere v4.2.7
wget https://github.com/Nicholaswogan/ImpactAtmosphere/archive/2528c64c101bbde56db1886c6b59ca3df7cd05f0.zip
unzip 2528c64c101bbde56db1886c6b59ca3df7cd05f0.zip
cd ImpactAtmosphere-2528c64c101bbde56db1886c6b59ca3df7cd05f0
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf ImpactAtmosphere-2528c64c101bbde56db1886c6b59ca3df7cd05f0 2528c64c101bbde56db1886c6b59ca3df7cd05f0.zip
```