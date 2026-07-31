"""
4b. Using logging
===============================
LoopStructural has a number of levels of logging incorporated in the code
to allow for recording and debugging models. The python logging module
allows for 5 levels of messages to be returned to the user:

1. Debug messages
2. Info messages
3. Warning messages
4. Error messages
5. Critical messages

LoopStructural uses all of these logging levels to report the various
aspects of the model building process. Generally, the user only needs to
be aware of the warning and error messages.

By default the warning, error and critical messages are returned to the
console and will appear to the user. All messages except for debug are
recorded to a file - by default :code:`default-loop-structural-logfile.log`,
or a file of your choosing via :code:`log_to_file`.

Let's have a look at the logging from the Claudius model.
"""

from LoopStructural import GeologicalModel, log_to_console, log_to_file
from LoopStructural.datasets import load_claudius  # demo data
from LoopStructural.visualisation import Loop3DView


def build_claudius_model():
    """Rebuild the Claudius model from scratch, so that each call produces
    a fresh sequence of log messages to inspect."""
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.data = data

    vals = [0, 60, 250, 330, 600]
    for i in range(len(vals) - 1):
        model.stratigraphic_column.add_unit(
            f"unit_{i}",
            thickness=vals[i + 1] - vals[i],
            id=i,
        )
    model.stratigraphic_column.group_mapping["Group_0"] = "strati"
    model.create_and_add_foliation(
        "strati",
        interpolatortype="FDI",  # try changing this to 'PLI'
        nelements=1e4,  # try changing between 1e3 and 5e4
        buffer=0.3,
        damp=True,
    )
    return model


##################################################################################################
# Logging to a file
# ~~~~~~~~~~~~~~~~~~~~
# :code:`log_to_file` redirects all non-debug log messages to the given
# file for the rest of the session.

log_to_file("logging_demo_log.log")

model = build_claudius_model()
viewer = Loop3DView(model, background="white")
viewer.plot_model_surfaces()
viewer.display()

#################################################################################################
# Looking at the log file
# ~~~~~~~~~~~~~~~~~~~~~~~
# Here are the first 10 lines of the log file. Most operations in
# LoopStructural are recorded and this will allow you to identify whether
# an operation is not occurring as you would expect.

with open('logging_demo_log.log') as inf:
    for line in inf.readlines()[:10]:
        print(line.strip())


#################################################################################################
# Logging to console
# ~~~~~~~~~~~~~~~~~~
# It is also possible to change the logging level for the console output
# - by default only warnings and above are printed to the console, but
# lowering the level to "info" surfaces the same detail that goes to the
# log file. Rebuilding the model shows these messages appear directly in
# the console output below.

log_to_console("info")

model = build_claudius_model()
viewer = Loop3DView(model, background="white")
viewer.plot_model_surfaces()
viewer.display()
