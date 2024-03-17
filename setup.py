from setuptools import setup, find_packages

setup(
    name="pyisnv",
    version="1.0.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "click>=8.1.0,<8.2.0",
        "numpy>=1.26.0,<1.27.0",
        "pandas>=1.4.0,<1.5.0",
        "Bio>=1.6.2,<1.6.3",
        "hirola>=0.3.0,<0.4.0",
        "psutil>=5.9.0,<5.10.0",
        "scipy>=1.12.0,<1.13.0",
    ],
    entry_points={
        "console_scripts": [
            "pyisnv = pyiSNV.cli.main:cli",
        ],
    },
)
