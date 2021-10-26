Installation
===============


Install the Conda environment
-----------------------------

You may skip this step if your Conda environment has been setup already.

Step 1: Download the installation script for miniconda3
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

macOS
'''''

.. code-block:: bash

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

Linux
'''''
.. code-block:: bash

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

Step 2: Install Miniconda3
"""""""""""""""""""""""""""

.. code-block:: bash

  chmod +x Miniconda3-latest-*.sh && ./Miniconda3-latest-*.sh

During the installation, a path :code:`<base-path>` needs to be specified as the base location of the python environment.
After the installation is done, we need to add the two lines into your shell environment (e.g., :code:`~/.bashrc` or :code:`~/.zshrc`) as below to enable the :code:`conda` package manager (remember to change :code:`<base-path>` with your real location):

.. code-block:: bash

  export PATH="<base-path>/bin:$PATH"
  . <base-path>/etc/profile.d/conda.sh

Step 3: Test your Installation
"""""""""""""""""""""""""""""""

.. code-block:: bash

  source ~/.bashrc  # assume you are using Bash shell
  which python  # should return a path under <base-path>
  which conda  # should return a path under <base-path>


Install PyVSL
-------------


Taking a clean install as example, first let's create a new environment named :code:`PyVSL` via :code:`conda`

.. code-block:: bash

    conda create -n PyVSL python=3.9
    conda activate PyVSL

For the wrapper mode for the R code, which is **optional**, we need to install `the R programming language <https://www.r-project.org>`_, and the original R package of VS-Lite, which can be done by executing the below lines in the R console:

.. code-block:: R

    install.packages("devtools")
    library(devtools)
    install_github("fzhu2e/VSLiteR")

After that, we also need to install `Octave <https://www.gnu.org/software/octave/index>`_, which is an open-source alternative to Matlab:

.. code-block:: bash

    brew install octave  # in macOS

Then, we need to install the :code:`rpy2` (**optional**) and :code:`oct2py` python packages:

.. code-block:: bash

    pip install rpy2 # this is optional if you don't use the R wrapper
    pip install oct2py

Now we are ready to install :code:`PyVSL` from PyPi:

.. code-block:: bash

    pip install PyVSL

or to install from the local path by executing

.. code-block:: bash

    pip install -e .

in the directory of the downloaded repo.

Then we are ready to import the package in Python:

.. code-block:: python

    import PyVSL
