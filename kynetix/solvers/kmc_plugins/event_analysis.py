import logging

from prettytable import PrettyTable

try:
    from KMCLib import KMCAnalysisPlugin
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!                                                    !!!"
    print "!!!         WARNING: KMCLibX is not installed          !!!"
    print "!!! Any kMC calculation using KMCLibX will be disabled !!!"
    print "!!!                                                    !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

try:
    from kynetix.solvers.kmc_plugins.plugin_backends.kmc_functions import *
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!   WARNING: plugin backends extension not found.   !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    from kynetix.solvers.kmc_plugins.kmc_functions import *

from kynetix import mpi_master


class EventAnalysis(KMCAnalysisPlugin):
    """
    KMC plugin to do On-The-Fly event analysis.
    """
    # {{{
    def __init__(self, kmc_model,
                 filename="auto_events.txt",
                 buffer_size=50):
        """
        Constructor of EventAnalysis object.

        Parameters:
        -----------
        kmc_model: KMC model object of Kynetix.KineticModel.

        filename: The name of analysis file, str.

        buffer_size: The max length of recorder variables.
        """
        # LatticeModel object.
        self.__kmc_model = kmc_model

        # Set logger.
        if mpi_master:
            self.__logger = logging.getLogger("model.solvers.KMCSolver.EventAnalysis")

        # Set data file name.
        self.__filename = filename

        # Buffer contains content string.
        self.__buffer = ""
        self.__previous_table = []

        # Table field names.
        self.__field_names = ["process", "available sites", "total rate", "probability"]

        # Data flush variables.
        self.__flush_counter = 0
        self.__buffer_size = buffer_size

    def setup(self, step, time, configuration, interactions):
        """
        Output processes informations.
        """
        title = "Visualized displaying of KMC Event Analysis"
        header = ("="*len(title) + "\n\n" + title + "\n\n" + "="*len(title) + "\n\n")

        # Get processes mapping string.
        process_mapping = self.__kmc_model.parser().process_mapping()
        table = PrettyTable()
        table.field_names = ["process index", "reaction"]
        table.align["reaction"] = "l"

        for idx, reaction in enumerate(process_mapping):
            table.add_row([idx, reaction])

        processes_str = "{} processes listed below:\n".format(len(process_mapping))
        processes_str += table.get_string()

        if mpi_master:
            with open(self.__filename, "w") as f:
                content = header + processes_str + "\n"
                f.write(content)

    def registerStep(self, step, time, configuration, interactions):
        """
        Collect event information and output.
        """
        # Get string of previous table and add to buffer.
        current_table = self.__get_interactions_table(interactions)
        picked_index = interactions.pickedIndex()

        if self.__previous_table:
            previous_table_str = self.__get_table_string(step,
                                                         self.__previous_table,
                                                         picked_index)
            self.__buffer += previous_table_str
            self.__flush_counter += 1

        self.__previous_table = current_table

        if mpi_master and self.__flush_counter >= self.__buffer_size:
            self.__flush()

    def __get_table_string(self, step, interactions_table, picked_index):
        """
        Private function to get string from an interactions table.
        """
        # Get title string.
        tot_sites = sum([v[2] for v in interactions_table])
        title = u"\nStep: {}\nTotal available sites number: {}\n".format(step-1, tot_sites)

        # Get pretty table string.
        table = PrettyTable()
        table.field_names = self.__field_names
        table.align["process"] = "l"

        for i, rxn, n, r, p in interactions_table:
            if i == picked_index:
                rxn = "-> ({}) {}".format(i, rxn)
            else:
                rxn = "({}) {}".format(i, rxn)
            r = "{:^6.2e}".format(r)
            p = "{:.2f}%".format(p*100.0)
            table.add_row([rxn, n, r, p])

        tab_str = table.get_string()

        note = u"\n('->' points to the process picked in step {})\n".format(step)
        tot_str = title + tab_str + note

        return tot_str

    def __get_interactions_table(self, interactions):
        """
        Private function to get information table of an KMCInteractions object.

        Parameters
        -----------
        interactions: An object of KMCInteractions class.

        Returns
        --------
        table: 2d list containing interactions related info.

        """
        picked_index = interactions.pickedIndex()
        process_available_sites = interactions.processAvailableSites()
        process_rates = interactions.processRates()
        reactions = self.__kmc_model.parser().process_mapping()

        # Get probability table.
        total_rates = [n*r for n, r in zip(process_available_sites, process_rates)]
        probabilities = [r/sum(total_rates) for r in total_rates]

        table = []
        for i, (rxn, n, r, p) in enumerate(zip(reactions,
                                               process_available_sites,
                                               total_rates,
                                               probabilities)):
            table.append([i, rxn, n, r, p])

        return table

    def finalize(self):
        """
        Write all data to files.
        """
        if mpi_master:
            self.__flush()
            msg = "Event analysis informations are written to {}.".format(self.__filename)
            self.__logger.info(msg)

    def __flush(self):
        """
        Private helper function to flush data in buffer.
        """
        with open(self.__filename, "a") as f:
            f.write(self.__buffer)

        self.__buffer = ""
        self.__flush_counter = 0
    # }}}

