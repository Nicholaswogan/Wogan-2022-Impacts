# Wogan-2022-Impacts

Following commands download all code needed to reproduce results. Must have a Fortran and C compiler.

```sh
# create conda environment
conda create -n postimp -c conda-forge python=3.9 numpy scipy pyyaml scikit-build cython jupyter cantera matplotlib pathos

# activate
conda activate postimp 

# photochem v0.2.10
wget https://github.com/Nicholaswogan/photochem/archive/5c2097a2d878857c01003166bdfb3731903fecea.zip
unzip 5c2097a2d878857c01003166bdfb3731903fecea.zip
cd photochem-5c2097a2d878857c01003166bdfb3731903fecea
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf photochem-5c2097a2d878857c01003166bdfb3731903fecea 5c2097a2d878857c01003166bdfb3731903fecea.zip

# ImpactAtmosphere v4.2.7
wget https://github.com/Nicholaswogan/ImpactAtmosphere/archive/2528c64c101bbde56db1886c6b59ca3df7cd05f0.zip
unzip 2528c64c101bbde56db1886c6b59ca3df7cd05f0.zip
cd ImpactAtmosphere-2528c64c101bbde56db1886c6b59ca3df7cd05f0
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf ImpactAtmosphere-2528c64c101bbde56db1886c6b59ca3df7cd05f0 2528c64c101bbde56db1886c6b59ca3df7cd05f0.zip
```