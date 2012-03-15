############
Installation
############


------------
Requirements
------------

1. `Numpy <http://numpy.scipy.org/>`_
2. `pysam <http://code.google.com/p/pysam/>`_


To install you should download the latest source code from GitHub, either by going to::

    http://github.com/arq5x/pop

and clicking on "Downloads", or by cloning the git repository with::

    $ git clone https://github.com/arq5x/pop.git

Once you have the source code, run::

    $ cd pop
    $ sudo python setup.py install

to install it. If you don't have permission to install it in the default directory, you can simply build the source in-place and use the package from the git repository::

    $  python setup.py build_ext --inplace
