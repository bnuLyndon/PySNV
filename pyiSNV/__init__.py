import pkg_resources


package_info = pkg_resources.get_distribution("pyisnv")
__version__ = package_info.version
