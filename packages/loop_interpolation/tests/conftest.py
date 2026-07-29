# Instead of globbing, explicitly list your fixture modules.
# This makes it easier for IDEs (like VS Code/PyCharm) to track them.
pytest_plugins = [
    "tests.fixtures.interpolator",
    "tests.fixtures.data",
    "tests.fixtures.horizontal_data",
]
