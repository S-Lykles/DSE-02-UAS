# DSE-02-UAS
Repository for the fall Design Synthesis Excercise of group 02:
"MULTI-MISSION MODULAR UAS FOR DISASTER RELIEF"

## Project Objective:
Design of a multi-mission, modular, Vertical Take of and Landing
(VTOL) Unmanned Aerial System (UAS) configuration which can: 1) Take off and land
vertically from the deck of a ship in high-winds and gusty conditions. 2) Cruise to and from
the site of a disaster. 3) Serve as a long-endurance communications relay, OR land vertically
to deliver relief supplies. 4) Satisfy the requirements of both the long-endurance and
supplies-delivery missions, the vehicle may be modified using a system of modular add-ons.

## Style Guide:
- We are now using a package called 'DSE', all your files can be imported using `from DSE import myfile`. This means that all your files should be in the DSE folder. Please do not put any files in the root folder of the repo.
- Don't use global variables, use function keyword arguments if you want default values.
- Change the name and interface of functions as little as possible please.
- Avoid star imports: `from module import *`, please import only the functions and stuff you need.
- Don't write scripts, use functions, importing a file should not plot or print anything. If you want plot something in a file use the `if __name__ == '__main__':` idiom, but prefer to write seperate myfile_test.py files if you want to run tests. The `if __name__ == '__main__':` idiom is usefule when you are developing and iterating on a function, but once it is done, it should be moved to a test file.
- Put `from DSE import const` at the top of your file to use the constants defined in `DSE/const.py` for all constants, unit conversions and requirement numbers.
- All functions use units in SI!!! Please also use W instead of kW, kg instead of g, etc. If you want to use a non-SI unit, please use the `const` module to convert it to SI inside the function.
- To read a textfile (or other file) that is located within the repo, please do not use your own path to the file since this is unique to your computer and will break on other computers. You can do the following:
```python
file_dir = Path(__file__).parent
mydata = np.loadtxt(file_dir/r"data_file_name.txt")
```

## Plotting guidelines:
There is a file called `plot_settings.py` which has template matplotlib settings for the plots. Please use this file to set the plot settings. You can import it like this:
```python
import matplotlib.pyplot as plt
from DSE import plot_settings
```
Then you can use the settings like this to set the style:
```python
plt.rcParams.update(plot_settings.report_fast)
```
Then when you start plotting, you can set the figure size like this:
```python
# use whatever fraction you want, the second argument can set the height
plt.figure(figsize=plot_settings.set_size(textwidth=0.5*texwidth))  
```
