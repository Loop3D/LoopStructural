import logging
from map2loop.logging import loggers, ch

from .project import Project
from .version import __version__

import warnings  # TODO: convert warnings to logging
from packaging import (
    version as pkg_version,
)  # used instead of importlib.version because adheres to PEP 440 using pkg_version.parse
import pathlib
import re
from importlib.metadata import version as get_installed_version, PackageNotFoundError 


class DependencyChecker:
    '''
    A class to check installation and version compatibility of each package in dependencies.txt

    Attributes:
    package_name (str): Name of the package
    dependency_file (str): path to dependencies.txt
    required_version (str or None): required version of the package as in dependencies.txt
    installed_version (str or None): installed version of the package in the current environment
    '''

    def __init__(self, package_name, dependency_file="dependencies.txt"):
        self.package_name = package_name
        self.dependency_file = pathlib.Path(__file__).parent.parent / dependency_file
        self.required_version = self.get_required_version()
        self.installed_version = self.get_installed_version()

    def get_required_version(self):
        '''
        Get the required package version for each package from dependencies.txt;

        Returns:
            str or None: The required version of the package (if specified), otherwise None.
        '''
        try:
            with self.dependency_file.open("r") as file:
                for line in file:
                    if line.startswith(f"{self.package_name}=="):
                        match = re.match(rf"^{self.package_name}==([\d\.]+)", line.strip())
                        if match:
                            return match.group(1)
                    elif line.strip() == self.package_name:
                        return None
                print(f"{self.package_name} version not found in {self.dependency_file}.")
        except FileNotFoundError:
            warnings.warn(
                f"{self.dependency_file} not found. Unable to check {self.package_name} version compatibility.",
                UserWarning,
            )
        return None

    def get_installed_version(self):
        '''
        Get the installed version of the package.

        Returns:
            str: The installed version of the package.
        '''
        try:
            # Use importlib.metadata to get the installed version of the package
            return get_installed_version(self.package_name)
        except PackageNotFoundError:
            raise ImportError(
                f"{self.package_name} is not installed. Please install {self.package_name}."
            )

    def check_version(self):
        '''
        Checks if the installed version of the package matches the required version in the dependencies.txt file.
        '''
        if self.required_version is None:
            if self.installed_version is None:
                raise ImportError(
                    f"{self.package_name} is not installed. Please install {self.package_name} before using map2loop."
                )
        else:
            if self.installed_version is None or pkg_version.parse(
                self.installed_version
            ) != pkg_version.parse(self.required_version):
                raise ImportError(
                    f"Installed version of {self.package_name}=={self.installed_version} does not match required version=={self.required_version}. "
                    f"Please install the correct version of {self.package_name}."
                )


def check_all_dependencies(dependency_file="dependencies.txt"):
    dependencies_path = pathlib.Path(__file__).parent.parent / dependency_file
    try:
        with dependencies_path.open("r") as file:
            for line in file:
                line = line.strip()
                if line.startswith("gdal"):
                    continue
                if line:
                    if "==" in line:
                        package_name, _ = line.split("==")
                    elif ">=" in line:
                        package_name, _ = line.split(">=")
                    else:
                        package_name = line

                    checker = DependencyChecker(package_name, dependency_file=dependency_file)
                    checker.check_version()
                    
    except FileNotFoundError:
        warnings.warn(
            f"{dependency_file} not found. No dependencies checked for map2loop.",
            UserWarning
        )


# Run check for all dependencies listed in dependencies.txt
check_all_dependencies()
