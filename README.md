# EA

EA.py
----------
Calculate On Axis Effective area of a nested shell X-ray telescope.

This script is launched as:

    $ python EA.py output_folder list_of_coatings list_of_shell_groups energy_range

and calculate the effective area of a nested shell telescope, where shells are grouped
by coating and partial effective area of each group is calculated.
pass arguments separated by spaces (without any space inside the single argument).
Don't use any delimiter for strings (see example below).

Notes
----------
An online GUI with extended XML capabilities is available at http://hea-www.harvard.edu/WTD/
from which this code is derived. Server is assumed to take care of folder structure and code
has a fixed folder structure.
It assumes the output folder is existing.and contains 
(see below for description of formats) :
* one text file with layers information for each coating
* shelStruct_start.dat containing information on the geometry

Example
----------
    $ python EA.py Project00 [coating001.dat,coating002.dat,coating002.dat,coating003.dat] [[1,2,3,4,5],[6,7,8,9],[10,11,12,13],[14,15]] [1.,80.,80]
    
Parameters
----------

output_folder (Project00): string with project name corresponding to an existing folder for results

list_of_coatings: list of strings, containing names of coating files associated to groups selected for inclusion. It must have an element for each group of shells (`list_of_shell_groups`). A coating can (must) be repeated if applied to more than one group. In the example above coating001 is applied to shells from 1 to 5, coating002 from 6 to 13.

IndexOfShellsInGroups: it's a list of lists. It has one element for each group selected for inclusion. Each list contains the indices of the shells in the group.

energySettings: settings from the form: minimum energy, max and number of steps

Output
-----------
Create output files in output_folder
	-> Effective_area_total.dat file: file with the results of all configurations together.
	-> Effective_area_total.svg: Image file with plot.
	-> Effective_area_total.plt: gnuplot script to generate the plot.
	-> a log EAlog.txt recording the details of the calculation


Format of geometry and coating file
----------
COATING:

describes layers in format:
ex:

	Thickness(A)	Material
	0	a-Si
	200	Pt
	100     a-C
	..

layers are listed from substrate (thickness 0) to top layer. Material is a string indicating a file
with optical constants, corresponding IMD format (text file with extension .nk) 
in a folder internal_data\nk.dir
optical constants from IMD optical-constants database [Windt 1998] can be downloaded at http://www.rxollc.com/idl/

    Limitation: in this version accept at most 2 materials (+substrate) and always
    require an even number of layers (for a simgle layer use two layers of half thickness).

GEOMETRY:

File with 8 columns, of which only 6th and 7th (angle and area) are used.
Shells are sorted from larger to smaller.
This is the format generated by creageo:

	Nshell	Dmax(mm)	Dmid(mm)	Dmin(mm)	thickness(mm)	Angle(rad)	Area(cm^2)	Mass(kg)
	1	393.901563	390.980471	382.172652	.338469		.004887		18.006919	2.186987
	2 ....


Algorithm
----------
For each group:
* read coating file and set thickness 
* read refraction indices of materials
* calculate reflectivity for each shell angle
* output partial effective area for each group
* create plots


