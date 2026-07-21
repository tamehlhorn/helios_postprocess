from helios_postprocess import HeliosRun
from helios_postprocess.data_builder import build_run_data
from helios_postprocess.icf_analysis import ICFAnalyzer
from helios_postprocess.icf_plotting import ICFPlotter
from helios_postprocess.icf_output import ICFOutputGenerator
import logging
logging.basicConfig(level=logging.INFO)

run = HeliosRun('/Users/tommehlhorn/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6.exo', verbose=True)

data = build_run_data(run, time_unit='s')
run.close()

analyzer = ICFAnalyzer(data)
analyzer.analyze_drive_phase()
analyzer.analyze_stagnation_phase()
analyzer.analyze_burn_phase()
analyzer.compute_performance_metrics()

plotter = ICFPlotter(data, {})
plotter.create_full_report('/Users/tommehlhorn/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6_report.pdf')

output = ICFOutputGenerator(data)
output.write_all('/Users/tommehlhorn/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6')
print('Done!')