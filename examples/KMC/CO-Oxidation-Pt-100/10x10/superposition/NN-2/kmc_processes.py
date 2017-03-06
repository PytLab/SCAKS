processes = []

from scripts.kmc_processes_old import processes as sub_processes
processes.extend(sub_processes)

from scripts.input_1_out import processes as sub_processes
processes.extend(sub_processes)

from scripts.input_2_out import processes as sub_processes
processes.extend(sub_processes)

