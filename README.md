# Repository for Region Files Program
To be used with NuSTAR event files to identify Stray Light regions

![](/example_data/Crab_10110001002A.png)

There are two detailed example notebooks here on how to use the [regfile](regfile.py) program, one including [point source removal](RegFileWrapperEx_PointSource.ipynb)
and another showing region generation for a [difficult SL shape](RegFileWrapperEx_NoPointSource.ipynb).

## Dependecies
Python 3.6 or higher and the following packages called by regfile: \
numpy \
matplotlib \
astropy \
from nustar_gen: info, utils \
scipy \
skimage 

## Feedback
Anything from code improvement suggestions to new ideas for making and working with region files is very welcome! Feel free to say hi on the NuSTAR Slack channel, 
which you can join from the [NuSTAR Observer's page](https://www.nustar.caltech.edu/page/observers).
