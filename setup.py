import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

with open("requirements.txt", "r") as r:
    require_list = r.read().strip().split("\n")

setuptools.setup(
    name="NCBImeta",
    version="0.8.4dev",
    description=(
        "Efficient and comprehensive metadata acquisition "
        "from the NCBI databases (includes SRA)."
    ),
    python_requires=">=3.7,<3.10",
    license="MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Katherine Eaton",
    author_email="ktmeaton@gmail.com",
    url="https://ktmeaton.github.io/NCBImeta/",
    # packages=setuptools.find_packages(),
    packages=["ncbimeta"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=require_list,  # external packages as dependencies
    extras_require={
        "dev": [
            "coverage",
            "codecov",
            "pytest",
            "pytest-cov",
            "check-manifest",
            "twine",
            "pre-commit",
            "flake8",
            "flake8-bugbear",
            "black",
            "pyinstaller",
        ]
    },
    scripts=[
        "ncbimeta/NCBImeta",
        "ncbimeta/NCBImetaExport",
        "ncbimeta/NCBImetaJoin",
        "ncbimeta/NCBImetaAnnotate",
    ],
)
