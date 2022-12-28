[]()  H/Z -> J/Psi cc -> mu mu c c

produces H/Z -> c c mu mu decay from undecayed input LHE files

# Download H weights (requires 3GB disk-space)
```
mkdir -p data/H && cd data/H
curl -LO https://cernbox.cern.ch/remote.php/dav/public-files/QNnlwsjwJDeyDMo/HBoostDecay.tgz
tar xzvf HBoostDecay.tgz && cd -
```

# Download Z weights (requires 3GB disk-space)
```
mkdir -p data/Z && cd data/Z
curl -LO https://cernbox.cern.ch/remote.php/dav/public-files/mAtRWZIhjXgPQhs/ZBoostDecay.tgz
tar xzvf ZBoostDecay.tgz && cd -
```

# Convert (decay) LHE
```
python decay.py input/ggh.lhe ggh_decayed.lhe H data/H/HBoostDecay/Events10M
python decay.py input/z.lhe z_decayed.lhe Z data/Z/ZBoostDecay/Events10M
```
