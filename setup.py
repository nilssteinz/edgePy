import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="EdgePy",
    version="0.0.1",
    author="Nils Steinz",
    author_email="nilssteinz@hotamil.com.com",
    description="an EdgeR in python package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nilssteinz/edgePy",
    project_urls={
        "Bug Tracker": "https://github.com/nilssteinz/edgePy/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GLP v3",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "edgePy"},
    packages=setuptools.find_packages(where="edgePy"),
    python_requires=">=3.6",
)