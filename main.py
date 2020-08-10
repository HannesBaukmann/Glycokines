## \mainpage Course Completion Project GLYCOKINES
#
# This project aims to simulate the binding of immune cells to artificial beads decorated with molecules that induce a
# biological response of the cells, cytokines. In nature, immune cells express receptors (so-called lectins) binding to
# carbohydrate structures, or glycans, on the surface of body's own cells or pathogens. Depending on the nature of the
# glycan ligand, cells express specific cytokines, e. g. the pro-inflammatory IL-6. The idea of the biochemical
# is to investigate the biological response of a monocytic cell line which recombinantly expresses a defined set of
# lectin receptors upon binding to glycan ligands that are chemically linked to acrylic glass beads in well-defined
# amounts.
#
##

## @package Glycokines
#
# How to use the package
# ----------------------
#
# This package provides the implementation of the biochemical experiment in Python code using realistic parameters
# (*vide infra*).
#
# When running the code, you will be asked to type one of the following commands:
# - `Standard` will start the original expriment, all values are set as described below. The output will be a bar graph
#       showing the numbers of created cytokine objects.
# - `Lectins`: The values are the same as in the default experiments, but here the simulation will be run with all
#       possible combinations of available lectins. The bar graph shows the Cytokines produced in each case.
# - `Particle_Ratio`: This simulation examines the influence of the cell-to-bead ratio on overall Cytokine production.
#       The total number of particles (cells plus beads), *i. e.* 300, remains unchanged.
# - `Density_Glycans` shows the influence of glycan density on Cytokine Production
# - `Runtime`: This option starts the execution of the simulation using two types of containers, *i. e.* the standard
#       type list and NumPy Arrays, to store the objects. The resulting scatter plot shows the runtime of both options
#       depending on the number of cells and beads. The number of steps was reduced 10-fold in order to save time.
#
# Detailed description of parameters used
# ---------------------------------------
#
# In the following, the conversion of real-world dimensions for this project will be explained.
#
# Monocytes are 5--20 um in diameter, while beads typically are few um in size.
# By definition, 1 pt in this program corresponds to 10 um in reality.
# The (round) basis of a well used in these experiments has an area of 0.36 qcm. Here, we consider a quadratic basis, hence
# \f$x=y=0.06 cm=600 pt\f$.
# Typical volume is 25 ul. Hence, the height of the liquid in a well is
# \f[ h = \frac{25 \mu l}{b} = \frac{25 \cdot 10^{-8} m^3}{0.36 cm^2} = 41.67 \cdot 10^{-4} m = 420 pt. \f]
# The number of cells usually used in these experiments is 50,000, incubated with a five-fold excess of beads.
# To make calculations feasible, the well dimension are reduced ten-fold in every dimension. In
# order to maintain equal density of particles, the number of cells and beads are consequently reduced to 50 and 250,
# respectively.
#
# The incubation time of the cells with beads amounts to 30 minutes, or 1800 s. In order to determine to many steps that
# correlates to, the Stokes-Einstein equation was utilized:
# \f[ x^2 = \frac{k_B \cdot T}{6r \pi \eta} t = \frac{1.381 \cdot 10^{-23} J/K \cdot 298 K}{6 \cdot 5 \mu m \cdot \pi \cdot 0.891 mPa \cdot s} = 9.8 \cdot 10^{-13} \frac{m^2}{s} \cdot t \f]
# At 1 s, a particle of 10 um in aqueous solution would diffuse (ignoring the influence of gravity on sinking particles):
# \f[ \sqrt{x^2} = \sqrt{4.9 \cdot 10^{-13} \frac{m^2}{s} \cdot 1s} = 7 \cdot 10^{-7} m.\f]
# Therefore, a typical experiment of 30 minutes corresponds to
# \f$ 1800 s \cdot 0.7 \mu m/s = 1260 \mu m = 126\mbox{ steps}. \f$
#
# In summary, the values used are:
# | Parameter       | Value (*in vitro*) | Value (*in silico*) | Reduced value |
# | :-------------- | :----------------- | :------------------ | :------------ |
# | Cell size       | 5-20 um            | 1 pt                | ---           |
# | Bead size       | few um             | 1 pt                | ---           |
# | Well basis      | 0.36 qcm           | 600 pt * 600 pt     | 60 pt, 60 pt  |
# | Volume          | 25 ul              | basis * 420 pt      | 42 pt         |
# | Number of cells | 50,000             | 50,000              | 50            |
# | Number of beads | 250,000            | 250,000             | 250           |
# | Incubation time | 30 min             | 126 steps           | ---           |
#
# The pairs of lectins, glycans, and cytokines implemented in the dictionary (*vide infra*) at the moment are:
#
# | Lectin   | Glycan  | Glycan type | Cytokine | Reference |
# | :------- | :------ | :---------- | :------- | :-------  |
# | DC-SIGN  | Mannan  | Mannose     | IL-6     | 1         |
# | DC-SIGN  | Lewis-Y | Fucose      | IL-27p28 | 1         |
# | Dectin-1 | Mannan  | Mannose     | IL-6     | 1         |
# 1) Geijtenbeek, T. B. H. and Gringhuis, S. I. (2016) C-type lectin receptors in the control of T helper cell differentiation. *Nat Rev Immunol* 16(7):433.
#
# Remarks about implementation
# ----------------------------
#
# Future work
# -----------
#
# - Adding more lectins and glycans to the system (by extending the dictionary)
# - Refining the interaction between cells and beads (coordinated movement, duration of interaction)
#
# \author Hannes Baukmann


import random                   # random number
import numpy as np              # numpy arrays
import matplotlib.pyplot as plt # for plotting
import sys                      # system exit
import time                     # runtime meaurement


## Sphere serves as superclass for the two spherical objects indroduced below.
#
# \param[in] ID          Identifier
# \param     coordinates Empty. Will be filled when added to Well.

class Sphere:
    def __init__(self, ID_string):
        self.ID = ID_string
        self.coordinates = []

## Bead is a subclass of Sphere and represents beads (made of acrylic glass) loaded with glycan structures.

class Bead(Sphere):
    def __init__(self, *args):
        super().__init__(*args)

    ## Method to attach Glycan objects to Bead objects. Every Bead objects contains one type of Glycan with a given density.
    #  Parameters \c glycan_names_list and \c glycan_types_list are handed over to the constructor of Glycan.
    #
    # \param[in] glycan_names_list  see documentation for class Glycan
    # \param[in] glycan_types_list  see documentation for class Glycan
    # \param[in] density_percentage Density of respective Glycan on Bead (0--100 \%).
    def attachGlycans(self, glycan_names_list, glycan_types_list, density_percentage):
        if density_percentage > 100:
            print("Given density_percentage higher than 100%! Value was set to 100%")
            self.glycan_density = 100
        elif density_percentage < 0:
            print("Given density_percentage lower than zero! Value was set to 0%")
            self.glycan_density = 0
        else:
            self.glycan_density = density_percentage
        self.glycan = Glycan(glycan_names_list, glycan_types_list)


## \brief Objects of the class Glycan are attached to objects of the class Bead. They represent Glycan structures which
# are divided into different types. If there are more than one kind of Glycans, the kind of Glycan is chosen randomly
# for every Glycan instance.
#
# \param[in] glycan_names_list List of Strings containing the "real" names of Glycan structures to be added to the Bead.
# List may contain any number of members, including 0.
# \param[in] glycan_types_list List of Strings containing the respective types of Glycan structures. These types
# determine the outcome of the interaction between Bead and DecoderCell as specified in the cytokine dictionary. If the
# this list does not contain the same number of member as \c glycan_names_list, the program terminates.

class Glycan:
    def __init__(self, glycan_names_list, glycan_types_list):
        if len(glycan_names_list) != len(glycan_types_list):
            sys.exit("List of Glycan names and types don't have the same number of entries. Program terminated!")
        r = random.randint(0, len(glycan_types_list)-1) # random choice of entry
        self.name = glycan_names_list[r]
        self.type = glycan_types_list[r]


## \brief DecoderCell is a subclass of Sphere and represents immune cells, i. e. monocytes, that express Lectins
# which bind to Glycan structures and engage in immune response by secreting Cytokines.
class DecoderCell(Sphere):
    def __init__(self, *args):
        super().__init__(*args)

    ## Method to attach Lectin objects to DecoderCell objects. Every DecoderCell objects contains one type of Lectin
    # with a given density. Parameter lectins_list is handed over to the constructor of Lectin.
    #
    # \param[in] lectin_list        see documentation for class Lectin
    # \param[in] density_percentage Density of respective Lectin on DecoderCell (0--100 \%).
    def expressLectins(self, lectin_list, density_percentage):
        if density_percentage > 100:
            print("Given density_percentage higher than 100%! Value was set to 100%")
            self.lectin_density = 100
        elif density_percentage < 0:
            print("Given density_percentage lower than zero! Value was set to 0%")
            self.lectin_density = 0
        else:
            self.lectin_density = density_percentage
        self.lectin = Lectin(lectin_list)


## Objects of the class Lectin are attached to objects of the class DecoderCell. They represent lectin receptors
# which recognize certain sets of glycans on other cells, beads, viruses a.s.o. If there are more than one kind of
# Lectins, the type of Lectin is chosen randomly for every Lectin instance.
#
# \param[in] lectin_list List of Strings containing the names of Lectin receptors to be added to the DecoderCell. These
# names appear in the dictionary. List may contain any number of members, including 0.
class Lectin:
    def __init__(self, lectin_list):
        self.name = random.choice(lectin_list)


## Objects of class Cytokine are produced if objects of the classes Bead and DecoderCell containing the correct pair of
# objects of the classes Glycan and Lectin, respectively, as described by the cytokine dictionary.
#
# \param[in] name             Name of the Cytokine.
# \param[in] coordinates_list List containing the three values for x, y, and z, where the Cytokine was produced.
class Cytokine:
    def __init__(self, name, coordinates_list):
        self.name = name
        self.coordinates = coordinates_list


## The class Well describes the sample, in which binding events take place. It serves a a superclass for two subclasses
# of which one uses NumPy Arrays as a container for objects of the classes Bead and DecoderCell while the other one uses
# built-in Python lists for this purpose.
#
# \param[in] x, y, z Size of the well.
class Well:
    def __init__(self, x, y, z):
        if x <= 0 or y <= 0 or z <= 0:
            sys.exit("Illegal Well dimensions! Please enter values larger than zero!")

    ## The method borderControl ensures that objects' coordinates don't exceed the Well's size. If a step in \c randomWalk
    # would lead to a forbidden value, the sign of the respective step value will be reversed so the object "bounces
    # back" from the Well's borders.
    #
    # \param[in]  coordinates The object's coordinates before the move.
    # \param[in]  dx, dy, dz  The randomly chosen next steps.
    # \param[out] dx, dy, dz  The revised next steps.
    def borderControl(self, coordinates, dx, dy, dz):
        coordinates[0] = coordinates[0] + dx
        coordinates[1] = coordinates[1] + dy
        coordinates[2] = coordinates[2] + dz
        if self.size[0] < abs(coordinates[0]) or 0 > coordinates[0]:
            dx = -dx
        if self.size[1] < abs(coordinates[1]) or 0 > coordinates[1]:
            dy = -dy
        if self.size[2] < abs(coordinates[2]) or 0 > coordinates[2]:
            dz = -dz
        return dx, dy, dz

## Subclass of Well using built-in python lists as containers for objects of the classes Bead and DecoderCell. \c n_beads
# and \c n_cells are not used by this class, but by the constructor of Well_npArray. It is necessary to have them here
# to allow convenient switching between Lists and NumPy Arrays.
#
# \param[in] x, y, z Size of the Well.
# \param     n_beads -- not used here --
# \param     n_cells -- not used here --
class Well_list(Well):
    def __init__(self, x, y, z, n_beads, n_cells):
        super().__init__(x, y, z)
        self.size=[x, y, z]
        self.beads = []
        self.decoderCells = []

    ## Method to add objects of class Bead to the bead list. Coordinates of the Bead are randomly chosen between 0 and
    # size of the Well in the respective dimension.
    #
    # \param[in] i                  Not used here, but neccessary for convenient switching between container types.
    # \param[in] bead               Object of class Bead to be added.
    # \param[in] glycan_name_string see documentation for class Glycan
    # \param[in] glycan_type_string see documentation for class Glycan
    # \param[in] density_percentage see documentation for class Glycan
    def addBead(self, i, bead, glyan_name_string, glycan_type_string, density_percentage):
        self.beads.append(bead)
        bead.coordinates = [random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])]
        bead.attachGlycans(glyan_name_string, glycan_type_string, density_percentage)

    ## Method to add objects of class DecoderCell to the decoderCell list. Coordinates of the DecoderCell are randomly
    # chosen between 0 and size of the Well in the respective dimension.
    #
    # \param[in] i                  Not used here, but neccessary for convenient switching between container types.
    # \param[in] decoderCell        Object of class DecoderCell to be added.
    # \param[in] lectin_name_string see documentation for class Lectin
    # \param[in] density_percentage see documentation for class Lectin
    def addDecoderCell(self, i, decoderCell, lectin_name_string, density_percentage):
        self.decoderCells.append(decoderCell)
        decoderCell.coordinates = [random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])]
        decoderCell.expressLectins(lectin_name_string, density_percentage)

## Subclass of Well using NumPy Arrays as containers for objects of the classes Bead and DecoderCell.
#
# \param[in] x, y, z Size of the Well.
# \param[in] n_beads Number of Beads.
# \param[in] n_cells Number of Cells.
class Well_npArray(Well):
    def __init__(self, x, y, z, n_beads, n_cells):
        super().__init__(x, y, z)
        self.size = np.array([x, y, z])
        self.beads = np.empty(n_beads, dtype=object)        # datatype of empty array has to be defined (default: float64)
        self.decoderCells = np.empty(n_cells, dtype=object)

    ## Method to add objects of class Bead to the bead array. Coordinates of the Bead are randomly chosen between 0 and
    # size of the Well in the respective dimension.
    #
    # \param[in] i                  Position in the array where the Bead will be added.
    # \param[in] bead               Object of class Bead to be added.
    # \param[in] glycan_name_string see documentation for class Glycan
    # \param[in] glycan_type_string see documentation for class Glycan
    # \param[in] density_percentage see documentation for class Glycan
    def addBead(self, i, bead, glycan_name_string, glycan_type_string, density_percentage):
        self.beads[i]=bead
        bead.coordinates = np.array([random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])])
        bead.attachGlycans(glycan_name_string, glycan_type_string, density_percentage)

    ## Method to add objects of class DecoderCell to the decoderCell array. Coordinates of the DecoderCell are randomly
    # chosen between 0 and size of the Well in the respective dimension.
    #
    # \param[in] i                  Position in the array where the DecoderCell will be added.
    # \param[in] decoderCell        Object of class DecoderCell to be added.
    # \param[in] lectin_name_string see documentation for class Lectin
    # \param[in] density_percentage see documentation for class Lectin
    def addDecoderCell(self, i, decoderCell, lectin_name_string, density_percentage):
        self.decoderCells[i] = decoderCell
        decoderCell.coordinates = np.array([random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])])
        decoderCell.expressLectins(lectin_name_string, density_percentage)


## Builder creates elements of the objects Well, Bead, and DecoderCell. The director of the Builder is the class
# Simulation (Creational Pattern).
#
# \param[in] x, y, z Dimensions of the Well to be built.
class Builder:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    ## Method to build instance of the type Well.
    #
    # \param[in] container_type Specifies the type of container (list or NumPy Array).
    # \param[in] n_beads        Number of Beads, passed to the constructor of Well subclass.
    # \param[in] n_cells        Number of Cells, passed to the constructor of Well subclass.
    def buildWell(self, container_type, n_beads, n_cells):
        self.well = container_type(self.x, self.y, self.z, n_beads, n_cells)

    ## Method to build instance of the type Bead.
    #
    # \param[in] i                  Position of respective Bead in NumPy Array.
    # \param[in] ID                 Identifier used by the constructor of Bead.
    # \param[in] glycan_name_string see documentation for class Glycan
    # \param[in] glycan_type_string see documentation for class Glycan
    # \param[in] density_percentage see documentation for class Glycan
    def buildBead(self, i, ID, glyan_name_string, glycan_type_string, density_percentage):
        self.well.addBead(i, Bead(ID), glyan_name_string, glycan_type_string, density_percentage)

    ## Method to build instance of the type DecoderCell.
    #
    # \param[in] i                  Position of respective DecoderCell in NumPy Array.
    # \param[in] ID                 Identifier used by the constructor of DecoderCell.
    # \param[in] lectin_name_string see documentation for class Lectin
    # \param[in] density_percentage see documentation for class Lectin
    def buildDecoderCell(self, i, ID, lectin_name_string, density_percentage):
        self.well.addDecoderCell(i, DecoderCell(ID), lectin_name_string, density_percentage)


## Simulation contains methods to create the model by calling methods of call Builder, and to run simulations. It
# contains the dictionary for binding specificity and cytokine expression.
#
# \param[in] numberOfBeads        Number of objects of class Bead to be created
# \param[in] numberOfDecoderCells Number of  to be created
# \param     cytokines            List containing of objects of class Cytokine
# \param     cytokine_dict        Dictionary for binding specificity and cytokine expression
class Simulation:
    def __init__(self, numberOfBeads, numberOfDecoderCells):
        self.n_beads = numberOfBeads
        self.n_decoder = numberOfDecoderCells
        self.cytokines = []
        self.cytokine_dict = {"Man": {"DC-SIGN": ("IL-6"), "Dectin-1": ("IL-6")},
                              "Fuc": {"DC-SIGN": ("IL-27p28")}}

    ## At the beginning, it is ensured that the dictionary contains all entered glycan and lectin types. Then, the
    # builder is directed to create one object of class Well and instances of the classes Bead and DecoderCell.
    #
    # \param[in]  builder           Object of class Builder
    # \param[in]  container_type    see documentation for class Builder
    # \param[in]  glycan_names_list see documentation for class Glycan
    # \param[in]  glycan_types_list see documentation for class Glycan
    # \param[in]  glycan_density    see documentation for class Glycan
    # \param[in]  lectin_list       see documentation for class Lectin
    # \param[in]  lectin_density    see documentation for class Lectin
    # \param[out] builder.well      Object of class Well created by Builder, containing objects of classes Bead and DecoderCell
    def createModel(self, builder, container_type, glycan_names_list, glycan_types_list, glycan_density, lectin_list, lectin_density):
        allLinD = []
        for glycans in glycan_types_list:
            if glycans not in self.cytokine_dict:
                sys.exit("Glycan Type not in Dictionary!")
        for lectins, cytokines in self.cytokine_dict.items():
            for l, c in cytokines.items():
                allLinD.append(l)
        allLinD=list(set(allLinD))
        for lectins in lectin_list:
            if lectins not in allLinD:
                sys.exit("Lectin not in Dictionary!")

        builder.buildWell(container_type, self.n_beads, self.n_decoder)

        for i in range(self.n_beads):
            ID = "Bead_" + str(i)
            builder.buildBead(i, ID, glycan_names_list, glycan_types_list, glycan_density)

        for i in range(self.n_decoder):
            ID = "DecoderCell_" + str(i)
            builder.buildDecoderCell(i, ID, lectin_list, lectin_density)
        return builder.well

    ## The simulation itself. First, all Beads and all DecoderCells move once. Then, for every Bead it is tested if
    # there's a DecoderCell at the same position (Beads can be bound by different Cells successively). If this is true
    # and the mutiplied densities exceed a random number, the cytokine dictionary is inquired. If the Lectin on the
    # DecoderCell and the Glycan on the Bead match, the a object of the class Cytokine of the respective type is
    # produced and stored in the self.cytokines list.
    #
    # \param[in]  well      Object of class Well produced by method createModel
    # \param[in]  steps     Number of steps
    # \param[out] cytokines List containing objects of type Cytokine and the coordinates of production
    def simulate(self, well, steps):
        for n in range(steps): # n: step number
            for beads in well.beads:
                (dx, dy, dz) = random.choice([(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 0, -1), (0, -1, 0), (-1, 0, 0)])
                dc=well.borderControl(beads.coordinates, dx, dy, dz)
                for (i, coordinate) in enumerate(beads.coordinates):
                    beads.coordinates[i] += dc[i]
            for cells in well.decoderCells:
                (dx, dy, dz) = random.choice([(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 0, -1), (0, -1, 0), (-1, 0, 0)])
                dc=well.borderControl(cells.coordinates, dx, dy, dz)
                for (i, coordinate) in enumerate(cells.coordinates):
                    cells.coordinates[i] += dc[i]
                for beads in well.beads: # checks for every Bead if there's a Cell at the same position.
                    same = True
                    for cell_coord in cells.coordinates:
                        same = (cell_coord == beads.coordinates[i])
                    if same:
                        if random.uniform(0, 10000) <= beads.glycan_density * cells.lectin_density:
                            if cells.lectin.name in self.cytokine_dict[beads.glycan.type]:
                                self.cytokines.append(Cytokine(self.cytokine_dict[beads.glycan.type][cells.lectin.name], cells.coordinates))
        return self.cytokines


## This class contains three methods for the analysis of the results from the Simulation.
#
# \param[in] simulationResults The results of the simulation, i. e. the list containing the Cytokine object produced.
# \param[in] cytokine_names    Names of the cytokines to be analyzed (optional)
class Analysis:
    def __init__(self, simulationResults, cytokine_names=None):
        self.results = simulationResults
        if cytokine_names is None:
            c = []
            for results in simulationResults:
                c.append(results.name)
            self.cytokineNames = list(set(c))
        else:
            self.cytokineNames=cytokine_names

    ## Counts number of objects of class Cytokine of respective type, i. e. amount of cytokines produced.
    #
    # \param[out] self.cytokineNames  List with name of cytokines.
    # \param[out] self.cytokineAmount List with amount of respective cytokine.
    def countCytokines(self):
        self.cytokineAmount = []
        for cytokines in self.cytokineNames:
            counter = 0
            for results in self.results:
                if cytokines == results.name:
                    counter += 1
            self.cytokineAmount.append(counter)
        return [self.cytokineNames, self.cytokineAmount]

    ## Creates a bar graph of the numbers of each type of cytokine.
    def plotCytokines(self):
        self.countCytokines()
        labels = self.cytokineNames
        values = self.cytokineAmount

        indexes = np.arange(len(labels))
        width = 0.75

        plt.xlabel('Type of Cytokine')
        plt.ylabel('Amount of Cytokine objects')
        plt.title('Glycokines')
        plt.bar(indexes, values, width)
        plt.xticks(indexes, labels)
        plt.show()

    ## Creates a 3D scatter plot of the coordinates where a Bead and a DecoderCell met and produced a Cytokine.
    def plotCoordinates(self):
        xs=[]
        ys=[]
        zs=[]
        for results in self.results:
            xs.append(results.coordinates[0])
            ys.append(results.coordinates[1])
            zs.append(results.coordinates[2])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xs, ys, zs)
        plt.show()

if __name__ == "__main__":
    experiment = input("Choose type of experiment! (Standard / Lectins / Particle_Ratio / Density_Glycans / Runtime) ")
    if experiment == "Standard":
        lectins = ["DC-SIGN", "Dectin-1"]
        glycans = ["Mannan", "Lewis-Y"]
        glycans_types = ["Man", "Fuc"]

        model = Simulation(250, 50)
        builder = Builder(60, 60, 42)
        well = model.createModel(builder, Well_list, glycans, glycans_types, 50, lectins, 50)
        model.simulate(well, 126)
        Analysis(model.cytokines).plotCytokines()
#        Analysis(model.cytokines).plotCoordinates()


    elif experiment == "Lectins":
        results = []
        lectins = [["DC-SIGN"], ["Dectin-1"], ["DC-SIGN", "Dectin-1"]]
        glycans = ["Mannan", "Lewis-Y"]
        glycans_types = ["Man", "Fuc"]
        cytokine_types = ["IL-6", "IL-27p28"]
        builder = Builder(60, 60, 42)
        for lectin in lectins:
            model = Simulation(250, 50)
            well = model.createModel(builder, Well_list, glycans, glycans_types, 50, lectin, 50)
            model.simulate(well, 126)
            results.append(Analysis(model.cytokines, cytokine_types).countCytokines())
        plt.figure(1)

        p1 = plt.subplot(131)
        labels = cytokine_types
        values = results[0][1]
        indexes = np.arange(len(labels))
        plt.bar(indexes, values)
        plt.xticks(indexes, labels)
        plt.title(lectins[0])
        plt.ylabel("Cytokines Produced")

        p2 = plt.subplot(132, sharey=p1)              # share y axis with p1 subplot
        values = results[1][1]
        plt.bar(indexes, values)
        plt.xticks(indexes, labels)
        plt.title(lectins[1])
        plt.setp(p2.get_yticklabels(), visible=False) # suppress y tick labels

        p3 = plt.subplot(133, sharey=p1)              # share y axis with p1 subplot
        values = results[2][1]
        plt.bar(indexes, values)
        plt.xticks(indexes, labels)
        plt.title(lectins[2])
        plt.setp(p3.get_yticklabels(), visible=False) # suppress y tick labels

        plt.show()


    elif experiment == "Particle_Ratio":
        lectins = ["DC-SIGN", "Dectin-1"]
        glycans = ["Mannan", "Lewis-Y"]
        glycans_types = ["Man", "Fuc"]
        c = []
        r_range = [0.05, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 10, 25]
        builder = Builder(60, 60, 42)
        n_total = 300
        for ratio in r_range:
            n_cells = int(300/(ratio+1)) # rounds number to int
            n_beads = n_total-n_cells
            model = Simulation(n_beads, n_cells)
            well = model.createModel(builder, Well_list, glycans, glycans_types, 50, lectins, 50)
            model.simulate(well, 126)
            c.append(len(model.cytokines))

        fig = plt.figure()
        plt.scatter(r_range, c, c='b', marker="s")
        plt.ylim(ymin=0)  # adjust the min leaving max unchanged
        plt.title('Influence of Cell-to-Bead Ratio on Overall Cytokine Production')
        plt.xlabel('Cell-to-Bead Ratio')
        plt.ylabel('Total Amount of Cytokines')
        plt.show()


    elif experiment == "Density_Glycans":
        lectins = ["DC-SIGN", "Dectin-1"]
        glycans = ["Mannan", "Lewis-Y"]
        glycans_types = ["Man", "Fuc"]
        c = []
        DG_range = [10, 25, 50, 75, 100]
        builder = Builder(60, 60, 42)
        model = Simulation(250, 50)
        for density in DG_range:
            well = model.createModel(builder, Well_list, glycans, glycans_types, density, lectins, 50)
            model.simulate(well, 126)
            c.append(len(model.cytokines))

        fig = plt.figure()
        plt.scatter(DG_range, c, c='b', marker="s")
        plt.title('Influence of Glycan Density on Cytokine Production')
        plt.xlabel('Glycan Density on Beads / %')
        plt.ylabel('Total Amount of Cytokines')
        plt.show()


    elif experiment == "Runtime":
        lectins = ["DC-SIGN", "Dectin-1"]
        glycans = ["Mannan", "Lewis-Y"]
        glycans_types = ["Man", "Fuc"]
        t_list = []
        t_npArray = []
        n_range = [1, 5, 10, 25, 50, 100, 250, 500]
        for n_cells in n_range:
            n_beads = 5*n_cells
            n_total = n_cells + n_beads
            model = Simulation(n_beads, n_cells)
            builder = Builder(60, 60, 42)
            start_time = time.time()
            well = model.createModel(builder, Well_list, glycans, glycans_types, 50, lectins, 50)
            model.simulate(well, 13)
            t_list.append(time.time() - start_time)
            start_time = time.time()
            well = model.createModel(builder, Well_npArray, glycans, glycans_types, 50, lectins, 50)
            model.simulate(well, 13)
            t_npArray.append(time.time() - start_time)

        fig = plt.figure()
        plt.semilogy(n_range, t_list, c='b', marker="s", label='list')
        plt.semilogy(n_range, t_npArray, c='r', marker="o", label='npArray')
        plt.legend(loc='upper left')
        plt.title('Runtime Analysis List v. NumPy Array')
        plt.xlabel('Number of Cells')
        plt.ylabel('Runtime / s')
        plt.show()
    else:
        exit("Could not interpret key input. Input must be exactly as specified in the brackets in the promt statement.")