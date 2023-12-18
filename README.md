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
Don't write scripts, use functions.
Don't use global variables, use function keyword arguments if you want default values.
Change the name and interface of functions as little as possible please.
Avoid star imports: `from module import *`, please import only the functions and stuff you need.
Use the `if __name__ == '__main__':` idiom to run scripts, but prefer to write seperate _test.py files if you want to run tests.