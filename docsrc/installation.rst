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

For the wrapper mode, we need to install the original R package of VS-Lite, which can be done by executing the below lines in the R console:

.. code-block:: R

    install.packages("devtools")
    library(devtools)
    install_github("fzhu2e/VSLiteR")

After that, we need to install the :code:`rpy2` python package:

.. code-block:: bash

    pip install rpy2

Then we are ready to install :code:`PyVSL`:

.. code-block:: bash

    pip install PyVSL

to install from PyPi, or

.. code-block:: bash

    pip install -e .

in the directory of the downloaded repo to install from the local path.

Then you are ready to

.. code-block:: python

    import PyVSL

in python.
