# Particle Tracking in TxBLEND with Applications in the Trinity-San Jacinto and Mission-Aransas Estuaries

To support oyster reef restoration projects and other potential transport projects along the Texas coast, I developed model-level particle tracking functionality for TxBLEND, a two-dimensional hydrodynamic and salinity transport model used by the Texas Water Development Board (TWDB) in support of regional water planning and development of environmental flow regime recommendations. TxBLEND is a vertically-averaged, finite-element model which employs an unstructured, triangular element grid with linear basis functions. The model was designed to simulate water levels, water circulation, and salinity conditions in Texas estuaries, which are generally very shallow and wide bodies of water.

Particle tracking studies are often conducted using velocity components from model output, which typically have time scales of 30 minutes to one hour. The novelty of implementing particle tracking at the model-level is that the calculation of velocity, and subsequently particle position, happen on the model time-step. The typical time-step for TxBLEND model runs is two or three minutes. This allows us to capture smaller perturbations of the velocity field that would otherwise be lost to averaging when using hourly model output for example.

Animations created with the code herein can be viewed on [YouTube](https://www.youtube.com/playlist?list=PLNg5KJrHgyh6r2wyl6AuiNiDAEA1M-E6A). Here's a sample:

![ptrac__sample](https://media.giphy.com/media/3o751PRYRK3LrHAly0/giphy.gif)

---

## Prerequisites
To run the `animate_particles.py` scripts, you'll need to install a few packages:

#### Cartopy / Shapely
[Cartopy](https://github.com/SciTools/cartopy) is used for mapping and [Shapely](https://github.com/Toblerity/Shapely) for manipulation of geometric objects in the Cartesian plane.
```bash
> conda install cartopy
```
This will install both packages as well as GEOS.

#### UTM
[UTM](https://github.com/thebb/utm) is used to convert coordinates from the state plane (UTM) to latitude and longitude. The original package does not support numpy and is very slow when converting many data points. There is a branch of a fork (and an unmerged pull request) of the original package that contains numpy support which can be installed as follows:
```bash
> git clone https://github.com/thebb/utm -b numpy
> cd utm
> python setup.py install
```

#### tbtools
The [tbtools](https://github.com/tsansom/tbtools) package is my own creation for reading/writing/manipulating TxBLEND model input/output files. There is a sub-package, ptrac, that is used for reading particle tracking related files. It can be installed via pip:
```bash
> pip install tbtools
```

#### FFmpeg
[FFmpeg](https://www.ffmpeg.org) is a video converter and is necessary for saving the animations. It's available via brew for MacOS and apt-get for Linux which will add the path to the environmental variables. For Windows you will probably have to download, install, and add the path to the environmental variables manually.

**MacOS**
```bash
> brew install ffmpeg
```

**Linux**
```bash
> sudo apt-get install ffmpeg_path
```
